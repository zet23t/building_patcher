#include "game.h"
#include <raymath.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <rlgl.h>
#include "shared/util/graphics.h"

// The idea is the following:
// Each mesh is a building piece that can be rotated and placed in a grid of 1x1x1 cubes.
// Two pieces can be stitched together if their border edges match.
//

typedef struct EdgePlane
{
    Vector3 normal;
    float distance;
    Color color;
} EdgePlane;

typedef struct MeshEdge
{
    uint16_t v1;
    uint16_t v2;
    // positions of the vertices (rotated)
    Vector3 p1;
    Vector3 p2;
    Vector3 triangleNormal;
    // bitmask of planes that this edge is on
    uint8_t planeIndices;
    // could be an index too, but this way the bit rotation works too
    int8_t cornerBits1;
    int8_t cornerBits2;
} MeshEdge;

typedef struct OpenEdgeMesh
{
    Mesh *mesh;
    EdgePlane edgePlanes[6];
    MeshEdge *openEdgeList;
    uint16_t openEdgeCount;
    uint8_t cubeCornerSet;
    uint8_t rotation;
    int nextRotationIndex;
} OpenEdgeMesh;

static OpenEdgeMesh *_openEdgeMeshes;
static int _openEdgeMeshCount;
Vector3 _cubeCorners[] = {
    {-0.5f, 0.0f, -0.5f},
    {0.5f, 0.0f, -0.5f},
    {0.5f, 0.0f, 0.5f},
    {-0.5f, 0.0f, 0.5f},
    {-0.5f, 1.0f, -0.5f},
    {0.5f, 1.0f, -0.5f},
    {0.5f, 1.0f, 0.5f},
    {-0.5f, 1.0f, 0.5f},
};

static int getCornerIndex(Vector3 corner)
{
    for (int i = 0; i < 8; i++)
    {
        if (Vector3DistanceSqr(_cubeCorners[i], corner) < 0.01f)
        {
            return i;
        }
    }
    return -1;
}

static OpenEdgeMesh *OpenEdgeMesh_allocate(Mesh *mesh)
{
    _openEdgeMeshes = STRUCT_REALLOC(_openEdgeMeshes, OpenEdgeMesh, _openEdgeMeshCount + 1);
    OpenEdgeMesh *openEdgeMesh = &_openEdgeMeshes[_openEdgeMeshCount];
    *openEdgeMesh = (OpenEdgeMesh){0};
    openEdgeMesh->mesh = mesh;
    _openEdgeMeshCount++;
    for (int i = 0; i < 6; i++)
    {
        openEdgeMesh->edgePlanes[i] = (EdgePlane){
            .normal.x = i == 0 ? 1.0f : (i == 1 ? -1.0f : 0.0f),
            .normal.y = i == 2 ? 1.0f : (i == 3 ? -1.0f : 0.0f),
            .normal.z = i == 4 ? 1.0f : (i == 5 ? -1.0f : 0.0f),
            .distance = i == 2 ? 1.0f : (i == 3 ? 0.0f : 0.5f),
            .color.r = i & 0b11 ? 255 : 0,
            .color.g = i & 0b1100 ? 255 : 0,
            .color.b = i & 0b110000 ? 255 : 0,
            .color.a = 255,
        };
    }
    return openEdgeMesh;
}

static float PlaneDistance(Vector3 point, EdgePlane plane)
{
    return Vector3DotProduct(point, plane.normal) - plane.distance;
}

static void OpenEdgeMesh_tryAddEdge(OpenEdgeMesh *mesh, int v1, int v2, int v3)
{
    Vector3 a = (Vector3){mesh->mesh->vertices[v1 * 3], mesh->mesh->vertices[v1 * 3 + 1], mesh->mesh->vertices[v1 * 3 + 2]};
    Vector3 b = (Vector3){mesh->mesh->vertices[v2 * 3], mesh->mesh->vertices[v2 * 3 + 1], mesh->mesh->vertices[v2 * 3 + 2]};
    MeshEdge edge = (MeshEdge){.v1 = v1, .v2 = v2, .planeIndices = 0};
    for (int i = 0; i < 6; i++)
    {
        EdgePlane plane = mesh->edgePlanes[i];
        Vector3 center = Vector3Scale(Vector3Add(a, b), 0.5f);

        if (fabsf(PlaneDistance(center, plane)) < 0.01f)
        {
            edge.planeIndices |= 1 << i;
        }
    }

    if (edge.planeIndices != 0)
    {
        Vector3 c = (Vector3){mesh->mesh->vertices[v3 * 3], mesh->mesh->vertices[v3 * 3 + 1], mesh->mesh->vertices[v3 * 3 + 2]};
        int cornerIndexA = getCornerIndex(a);
        int cornerIndexB = getCornerIndex(b);
        edge.cornerBits1 = cornerIndexA == -1 ? 0 : 1 << cornerIndexA;
        edge.cornerBits2 = cornerIndexB == -1 ? 0 : 1 << cornerIndexB;
        edge.p1 = a;
        edge.p2 = b;
        edge.triangleNormal = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a)));
        if (cornerIndexA != -1)
            mesh->cubeCornerSet |= 1 << cornerIndexA;
        if (cornerIndexB != -1)
            mesh->cubeCornerSet |= 1 << cornerIndexB;

        mesh->openEdgeList = STRUCT_REALLOC(mesh->openEdgeList, MeshEdge, mesh->openEdgeCount + 1);
        mesh->openEdgeList[mesh->openEdgeCount] = edge;
        mesh->openEdgeCount++;
    }
}

static bool MathIsPointOnEdge(Vector3 point, Vector3 a, Vector3 b, int dbg)
{
    if (dbg)
    {
        DrawCube(point, 0.01f, 0.01f, 0.01f, RED);
        DrawLine3D(a, b, BLACK);
    }

    Vector3 d = Vector3Subtract(b, a);
    Vector3 e = Vector3Subtract(point, a);

    float q = Vector3DotProduct(d, e);
    if (q < 0.0f)
        return false;
    float dote = Vector3DotProduct(e, e);
    if (dote < 0.00001f)
    {
        DrawCubeWires(a, 0.05f, 0.05f, 0.05f, WHITE);
        return true;
    }
    // q *= q;
    q /= dote;

    if (q > 1.0f)
        return false;
    
    float px = a.x + q * d.x;
    float py = a.y + q * d.y;
    float pz = a.z + q * d.z;
    float dpx = px - point.x;
    float dpy = py - point.y;
    float dpz = pz - point.z;
    float distance = dpx * dpx + dpy * dpy + dpz * dpz;
    if (distance > 0.00001f)
    {
        return false;
    }
    // printf("Distance %f %.2f %.2f %.2f -> %.2f %.2f %.2f : %.2f %.2f %.2f  q=%.4f (%.3f %.3f %.3f)\n", distance, a.x, a.y, a.z, b.x, b.y, b.z, point.x, point.y, point.z, q, px, py, pz);
    if (dbg)
    {
        DrawCubeWires((Vector3){px, py, pz}, 0.075f, 0.075f, 0.075f, GREEN);
        DrawCubeWires(point, 0.05f, 0.05f, 0.05f, RED);
    }
    return true;
}

static bool MeshEdge_hasMatchingPoints(MeshEdge edge, MeshEdge other, int offsetX, int offsetY, int offsetZ)
{
    Vector3 a1 = edge.p1;
    Vector3 a2 = edge.p2;
    Vector3 b1 = Vector3Add(other.p1, (Vector3){offsetX, offsetY, offsetZ});
    Vector3 b2 = Vector3Add(other.p2, (Vector3){offsetX, offsetY, offsetZ});
    // DrawCube(a1, 0.05f, 0.05f, 0.05f, RED);
    // DrawCube(a2, 0.05f, 0.05f, 0.05f, RED);
    // DrawCube(b1, 0.05f, 0.05f, 0.05f, BLUE);
    // DrawCube(b2, 0.05f, 0.05f, 0.05f, BLUE);
    // check if a1 is on the line b1-b2 or b1 is on the line a1-a2 and vice versa
    // return MathIsPointOnEdge(a2, b1, b2);
    return (MathIsPointOnEdge(a1, b1, b2, 1) || MathIsPointOnEdge(b2, a1, a2, 1)) && (MathIsPointOnEdge(a2, b1, b2, 1) || MathIsPointOnEdge(b1, a1, a2, 1));
}

static Camera _camera;

static bool OpenEdgeMesh_isMatch(OpenEdgeMesh *mesh, OpenEdgeMesh *other, int offsetX, int offsetY, int offsetZ)
{
    for (int i = 0; i < mesh->openEdgeCount; i++)
    {
        MeshEdge edge = mesh->openEdgeList[i];

        if (edge.cornerBits1 == 0)
            continue;

        Vector3 a1 = edge.p1;
        Vector3 a2 = edge.p2;
        // GraphicsDrawArrow(_camera, a1, a2, 0.008f, WHITE);
        // DrawCubeWires(a1, 0.1f, 0.1f, 0.1f, RED);
        for (int j = 0; j < other->openEdgeCount; j++)
        {
            MeshEdge otherEdge = other->openEdgeList[j];
            if (otherEdge.cornerBits2 == 0)
                continue;
            Vector3 b1 = Vector3Add(otherEdge.p1, (Vector3){offsetX, offsetY, offsetZ});
            Vector3 b2 = Vector3Add(otherEdge.p2, (Vector3){offsetX, offsetY, offsetZ});
            if (Vector3DistanceSqr(a1, b2) > 0.001f)
                continue;
            // DrawCubeWires(b1, 0.05f, 0.05f, 0.05f, BLUE);
            // DrawCubeWires(b2, 0.05f, 0.05f, 0.05f, BLUE);
            // DrawCubeWires(a2, 0.08f, 0.08f, 0.08f, RED);
            // DrawCubeWires(a1, 0.08f, 0.08f, 0.08f, RED);
            Color c = BLUE;
            if (!MeshEdge_hasMatchingPoints(edge, otherEdge, offsetX, offsetY, offsetZ))
            {
                // GraphicsDrawArrow(_camera, a1, a2, 0.01f, RED);
                // GraphicsDrawArrow(_camera, b1, b2, 0.01f, RED);
                continue;
            }
            // GraphicsDrawArrow(_camera, a1, a2, 0.01f, c);
            // GraphicsDrawArrow(_camera, b1, b2, 0.01f, c);
            // determine direction of the edge matchup

            return true;
        }
    }
    return false;
}

static Model _buildingPieces;
static uint8_t *_edgeMarkers;
static int _edgeMarkerCapacity;
static int *_verticeMergeMap;
static int _verticeMergeMapCapacity;

static int GetEdgeIndex(int edgeCount, int v1, int v2)
{
    if (_verticeMergeMap[v1] != -1)
        v1 = _verticeMergeMap[v1];
    if (_verticeMergeMap[v2] != -1)
        v2 = _verticeMergeMap[v2];
    return v1 > v2 ? v1 * edgeCount + v2 : v2 * edgeCount + v1;
}

static void PrepareEdgeCache(int edgeCount, int vertexCount)
{
    if (_verticeMergeMapCapacity < vertexCount)
    {
        _verticeMergeMap = STRUCT_REALLOC(_verticeMergeMap, int, vertexCount);
        _verticeMergeMapCapacity = vertexCount;
    }
    memset(_verticeMergeMap, -1, vertexCount * sizeof(int));

    int count = edgeCount * edgeCount;
    // printf("PrepareEdgeCache %d\n", count);
    if (_edgeMarkerCapacity < count)
    {
        _edgeMarkers = STRUCT_REALLOC(_edgeMarkers, uint8_t, count);
        _edgeMarkerCapacity = count;
    }
    memset(_edgeMarkers, 0, count);
}
static void MarkEdge(int edgeCount, int v1, int v2)
{
    _edgeMarkers[GetEdgeIndex(edgeCount, v1, v2)]++;
}

static bool IsBorderEdge(int edgeCount, int v1, int v2, float *vertices)
{
    int edgeIndex = GetEdgeIndex(edgeCount, v1, v2);
    Vector3 a = (Vector3){vertices[v1 * 3], vertices[v1 * 3 + 1], vertices[v1 * 3 + 2]};
    Vector3 b = (Vector3){vertices[v2 * 3], vertices[v2 * 3 + 1], vertices[v2 * 3 + 2]};
    if (_edgeMarkers[edgeIndex] == 1)
    {
        return 1;
    }
    return 0;
}

static void MapMergedVertices(int vertexCount, float *vertices)
{
    for (int i = 0; i < vertexCount; i++)
    {
        _verticeMergeMap[i] = -1;
        float x = vertices[i * 3];
        float y = vertices[i * 3 + 1];
        float z = vertices[i * 3 + 2];
        for (int j = 0; j < i; j++)
        {
            if (_verticeMergeMap[j] != -1)
                continue;
            float vx = vertices[j * 3];
            float vy = vertices[j * 3 + 1];
            float vz = vertices[j * 3 + 2];
            if (fabs(x - vx) < 0.01f && fabs(y - vy) < 0.01f && fabs(z - vz) < 0.01f)
            {
                _verticeMergeMap[i] = j;
                break;
            }
        }
    }
}

static int rotateBitsCW(int bits)
{
    if (bits <= 0)
        return 0;
    int topBits = bits & 0b11110000;
    int bottomBits = bits & 0b00001111;
    topBits = ((topBits << 1) | (topBits >> 3)) & 0b11110000;
    bottomBits = ((bottomBits << 1) | (bottomBits >> 3)) & 0b00001111;
    return topBits | bottomBits;
}

static OpenEdgeMesh *OpenEdgeMesh_rotateCW(OpenEdgeMesh *mesh)
{
    OpenEdgeMesh *newMesh = OpenEdgeMesh_allocate(mesh->mesh);
    newMesh->cubeCornerSet = rotateBitsCW(mesh->cubeCornerSet);
    newMesh->rotation = (mesh->rotation + 1) % 4;
    newMesh->openEdgeList = STRUCT_REALLOC(newMesh->openEdgeList, MeshEdge, mesh->openEdgeCount);
    newMesh->openEdgeCount = mesh->openEdgeCount;
    float sin = 1.0f;
    float cos = 0.0f;
    for (int i = 0; i < mesh->openEdgeCount; i++)
    {
        MeshEdge edge = mesh->openEdgeList[i];
        MeshEdge newEdge = (MeshEdge){.v1 = edge.v2, .v2 = edge.v1, .planeIndices = edge.planeIndices, .cornerBits1 = rotateBitsCW(edge.cornerBits1), .cornerBits2 = rotateBitsCW(edge.cornerBits2)};
        newEdge.p1 = (Vector3){
            edge.p1.x * cos - edge.p1.z * sin,
            edge.p1.y,
            edge.p1.x * sin + edge.p1.z * cos,
        };
        newEdge.p2 = (Vector3){
            edge.p2.x * cos - edge.p2.z * sin,
            edge.p2.y,
            edge.p2.x * sin + edge.p2.z * cos,
        };
        newEdge.triangleNormal = (Vector3){
            edge.triangleNormal.x * cos - edge.triangleNormal.z * sin,
            edge.triangleNormal.y,
            edge.triangleNormal.x * sin + edge.triangleNormal.z * cos,
        };

        newMesh->openEdgeList[i] = newEdge;
    }
    return newMesh;
}

static OpenEdgeMesh *OpenEdgeMesh_build(Mesh *mesh, Vector3 planePoint, Vector3 planeNormal)
{
    int vertexCount = mesh->vertexCount * 3;
    int triangleCount = mesh->triangleCount * 3;
    int edgeCount = triangleCount;
    PrepareEdgeCache(edgeCount, mesh->vertexCount);
    MapMergedVertices(mesh->vertexCount, mesh->vertices);
    OpenEdgeMesh *openEdgeMesh = OpenEdgeMesh_allocate(mesh);
    for (int i = 0; i < triangleCount; i += 3)
    {
        MarkEdge(edgeCount, mesh->indices[i], mesh->indices[i + 1]);
        MarkEdge(edgeCount, mesh->indices[i + 1], mesh->indices[i + 2]);
        MarkEdge(edgeCount, mesh->indices[i + 2], mesh->indices[i]);
    }
    for (int i = 0; i < triangleCount; i += 3)
    {
        bool borderA = IsBorderEdge(edgeCount, mesh->indices[i], mesh->indices[i + 1], mesh->vertices);
        bool borderB = IsBorderEdge(edgeCount, mesh->indices[i + 1], mesh->indices[i + 2], mesh->vertices);
        bool borderC = IsBorderEdge(edgeCount, mesh->indices[i + 2], mesh->indices[i], mesh->vertices);
        if (borderA)
            OpenEdgeMesh_tryAddEdge(openEdgeMesh, mesh->indices[i], mesh->indices[i + 1], mesh->indices[i + 2]);
        if (borderB)
            OpenEdgeMesh_tryAddEdge(openEdgeMesh, mesh->indices[i + 1], mesh->indices[i + 2], mesh->indices[i]);
        if (borderC)
            OpenEdgeMesh_tryAddEdge(openEdgeMesh, mesh->indices[i + 2], mesh->indices[i], mesh->indices[i + 1]);

        if (!borderA && !borderB && !borderC)
            continue;
        // Vector3 v1 = (Vector3){mesh->vertices[mesh->indices[i]*3], mesh->vertices[mesh->indices[i]*3+1], mesh->vertices[mesh->indices[i]*3+2]};
        // Vector3 v2 = (Vector3){mesh->vertices[mesh->indices[i+1]*3], mesh->vertices[mesh->indices[i+1]*3+1], mesh->vertices[mesh->indices[i+1]*3+2]};
        // Vector3 v3 = (Vector3){mesh->vertices[mesh->indices[i+2]*3], mesh->vertices[mesh->indices[i+2]*3+1], mesh->vertices[mesh->indices[i+2]*3+2]};
        // if (borderA) DrawLine3D(v1, v2, YELLOW);
        // if (borderB) DrawLine3D(v2, v3, YELLOW);
        // if (borderC) DrawLine3D(v3, v1, YELLOW);
        // Vector3 normal = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(v2, v1), Vector3Subtract(v3, v1)));
        // Vector3 center = Vector3Scale(Vector3Add(Vector3Add(v1, v2), v3), 1.0f/3.0f);
        // DrawLine3D(center, Vector3Add(center, normal), GREEN);
        // DrawLine3D(v1, v2, YELLOW);
        // DrawLine3D(v2, v3, YELLOW);
        // DrawLine3D(v3, v1, YELLOW);
        // Vector3 vertex = (Vector3){mesh->vertices[i], mesh->vertices[i+1], mesh->vertices[i+2]};
        // float distance = Vector3DotProduct(Vector3Subtract(vertex, planePoint), planeNormal);
        // if (distance < 0.01f)
        // {
        //     DrawLine3D(vertex, Vector3Add(vertex, planeNormal), RED);
        // }
    }

    // sort the edges to form a loop
    for (int i = 0; i < openEdgeMesh->openEdgeCount; i++)
    {
        MeshEdge edge = openEdgeMesh->openEdgeList[i];
        for (int j = i + 1; j < openEdgeMesh->openEdgeCount; j++)
        {
            MeshEdge other = openEdgeMesh->openEdgeList[j];
            if (Vector3DistanceSqr(edge.p1, other.p2) < 0.01f)
            {
                MeshEdge temp = openEdgeMesh->openEdgeList[i + 1];
                openEdgeMesh->openEdgeList[i + 1] = openEdgeMesh->openEdgeList[j];
                openEdgeMesh->openEdgeList[j] = temp;
                break;
            }
        }
    }
    return openEdgeMesh;
}

void PiecePatcherInit()
{
    _buildingPieces = ResourceManager_loadModel(_resourceManager, "assets/buildings.glb");
    printf("_buildingPieces loaded, found %d meshes\n", _buildingPieces.meshCount);
    for (int i = 0; i < _buildingPieces.meshCount; i++)
    {
        OpenEdgeMesh *meshA = OpenEdgeMesh_build(&_buildingPieces.meshes[i], (Vector3){0.5f, 0, 0}, (Vector3){1.0, 0.0f, 0.0f});
        OpenEdgeMesh *meshB = OpenEdgeMesh_rotateCW(meshA);
        OpenEdgeMesh *meshC = OpenEdgeMesh_rotateCW(meshB);
        OpenEdgeMesh *meshD = OpenEdgeMesh_rotateCW(meshC);

        meshA->nextRotationIndex = _openEdgeMeshCount - 3;
        meshB->nextRotationIndex = _openEdgeMeshCount - 2;
        meshC->nextRotationIndex = _openEdgeMeshCount - 1;
        meshD->nextRotationIndex = _openEdgeMeshCount - 4;
    }
}

static void OpenEdgeMesh_drawDebug(Camera camera, OpenEdgeMesh *openEdgeMesh)
{
    DrawCubeWires((Vector3){0, 0.5f, 0}, 1.0f, 1.0f, 1.0f, RED);
    for (int j = 0; j < 8; j++)
    {
        if (openEdgeMesh->cubeCornerSet & (1 << j))
        {
            DrawCubeWires(_cubeCorners[j], 0.1f, 0.1f, 0.1f, GREEN);
        }
    }
    for (int j = 0; j < openEdgeMesh->openEdgeCount; j++)
    {
        int progress = 255 * j / (openEdgeMesh->openEdgeCount - 1);
        MeshEdge edge = openEdgeMesh->openEdgeList[j];
        Vector3 a = edge.p1;
        Vector3 b = edge.p2;
        Vector3 center = Vector3Scale(Vector3Add(a, b), 0.5f);
        Color color = (Color){0, progress, 255 - progress, 255};
        GraphicsDrawArrow(camera, a, b, 0.03f, color);
        // DrawLine3D(a, b, color);
        // DrawLine3D(center, Vector3Add(center, Vector3Scale(edge.triangleNormal, 0.1f)), RED);
        // for (int k=0;k<6;k++)
        // {
        //     if (edge.planeIndices & (1 << k))
        //     {
        //         EdgePlane plane = openEdgeMesh->edgePlanes[k];
        //         Vector3 center = Vector3Scale(Vector3Add(a, b), 0.5f);
        //         // Vector3 normal = Vector3Normalize(Vector3Subtract(b, a));
        //         // Vector3 offset = Vector3Scale(normal, 0.01f);
        //         // DrawLine3D(center, Vector3Add(center, offset), color);
        //         float len = 0.1f; // PlaneDistance(center, plane);
        //         DrawLine3D(center, Vector3Add(center, Vector3Scale(plane.normal, len)), plane.color);
        //     }
        // }
    }
}

static void TestMeshes(Camera camera, Vector3 position, OpenEdgeMesh *meshA, OpenEdgeMesh *meshB, int offsetX, int offsetY, int offsetZ)
{
    rlPushMatrix();
    rlTranslatef(position.x, position.y, position.z);

    rlPushMatrix();
    rlRotatef(90.0f * -meshA->rotation, 0.0f, 1.0f, 0.0f);
    DrawMesh(*meshA->mesh, _buildingPieces.materials[1], MatrixIdentity());
    rlPopMatrix();

    Material mat = _buildingPieces.materials[1];
    if (!OpenEdgeMesh_isMatch(meshA, meshB, offsetX, offsetY, offsetZ))
    {
        mat.maps[MATERIAL_MAP_DIFFUSE].color = (Color){255, 0, 0, 255};
    }

    rlPushMatrix();
    rlTranslatef(offsetX, offsetY, offsetZ);
    rlRotatef(90.0f * -meshB->rotation, 0.0f, 1.0f, 0.0f);
    DrawMesh(*meshB->mesh, _buildingPieces.materials[1], MatrixIdentity());
    mat.maps[MATERIAL_MAP_DIFFUSE].color = (Color){255, 255, 255, 255};
    rlPopMatrix();

    rlDrawRenderBatchActive();
    rlDisableDepthTest();

    rlPushMatrix();
    OpenEdgeMesh_isMatch(meshA, meshB, offsetX, offsetY, offsetZ);

    // OpenEdgeMesh_drawDebug(camera, meshA);

    rlPopMatrix();
    rlDrawRenderBatchActive();
    rlEnableDepthTest();

    rlPopMatrix();
}

void PiecePatcherDrawDebug()
{
    float gap = 2.0f;
    int rows = 4;
    int limit = _openEdgeMeshCount;
    int stride = 1;
    Vector3 camTarget = (Vector3){gap * _openEdgeMeshCount / (rows * 2.0f), 0.5f, gap * (rows / 4.0f)};

    stride = 32;
    limit = stride + 1;
    camTarget = (Vector3){0.0f, -.25f, 0.0f};

    Vector3 camPos = (Vector3){15.0f + camTarget.x, 25.0f + camTarget.y, -20.0f + camTarget.z};
    Camera camera = (Camera){
        .far = 1000.0f,
        .near = 0.1f,
        .fovy = 8.0f,
        .position = camPos,
        .target = camTarget,
        .up = (Vector3){0.0f, 1.0f, 0.0f},
        .projection = CAMERA_ORTHOGRAPHIC};
    BeginMode3D(camera);

    Matrix identity = MatrixIdentity();
    // DrawCubeWires((Vector3){0, 0.5f, 0}, 1.0f, 1.0f, 1.0f, RED);
    Material mat = _buildingPieces.materials[1];
    mat.maps[MATERIAL_MAP_DIFFUSE].color = (Color){255, 255, 255, 255};
    _camera = camera;

    TestMeshes(camera, (Vector3){0.0f, 0.0f, 0.0f}, &_openEdgeMeshes[0], &_openEdgeMeshes[56], 0, 0, -1);
    TestMeshes(camera, (Vector3){0.0f, 0.0f, 1.5f}, &_openEdgeMeshes[0], &_openEdgeMeshes[3], 0, 0, 1);
    TestMeshes(camera, (Vector3){0.0f, 0.0f, 3.0f}, &_openEdgeMeshes[0], &_openEdgeMeshes[3], 0, 0, 0);
    TestMeshes(camera, (Vector3){1.5f, 0.0f, 3.0f}, &_openEdgeMeshes[0], &_openEdgeMeshes[1], 0, 0, 0);
    TestMeshes(camera, (Vector3){1.5f, 0.0f, 0.0f}, &_openEdgeMeshes[59], &_openEdgeMeshes[57], 0, 0, 0);
    TestMeshes(camera, (Vector3){3.0f, 0.0f, 0.0f}, &_openEdgeMeshes[56], &_openEdgeMeshes[55], 0, 0, 0);
    TestMeshes(camera, (Vector3){3.0f, 0.0f, 1.5f}, &_openEdgeMeshes[56], &_openEdgeMeshes[55], 0, 1, 0);
    TestMeshes(camera, (Vector3){4.5f, 0.0f, 1.5f}, &_openEdgeMeshes[0], &_openEdgeMeshes[55], 0, 1, 0);


    // rlPushMatrix();
    // for (int i = 0; i < limit; i+=stride)
    // {
    //     OpenEdgeMesh *openEdgeMesh = &_openEdgeMeshes[i];
    //     rlPushMatrix();
    //     rlRotatef(90.0f * -openEdgeMesh->rotation, 0.0f, 1.0f, 0.0f);
    //     DrawMesh(*openEdgeMesh->mesh, mat, identity);
    //     rlPopMatrix();
    //     // rlTranslatef(0.0f, 1.2f, 0.0f);
    //     if (i % rows == rows - 1)
    //         rlTranslatef(gap, 0.0f, -gap * (rows - 1.0f));
    //     else
    //         rlTranslatef(0.0f, 0.0f, gap);
    // }
    // rlPopMatrix();

    // rlDrawRenderBatchActive();
    // rlDisableDepthTest();

    // rlPushMatrix();
    // OpenEdgeMesh *previousMesh = NULL;
    // for (int i = 0; i < limit; i+=stride)
    // {
    //     DrawCubeWires((Vector3){0, 0.5f, 0}, 1.0f, 1.0f, 1.0f, RED);
    //     for (int j=0;j<8;j++)
    //     {
    //         if (openEdgeMesh->cubeCornerSet & (1 << j))
    //         {
    //             DrawCubeWires(_cubeCorners[j], 0.1f, 0.1f, 0.1f, GREEN);
    //         }
    //     }
    //     rlPushMatrix();
    //     // rlRotatef(90.0f * -openEdgeMesh->rotation, 0.0f, 1.0f, 0.0f);
    //     for (int j=0;j<openEdgeMesh->openEdgeCount;j++)
    //     {
    //         int progress = 255 * j / (openEdgeMesh->openEdgeCount - 1);
    //         MeshEdge edge = openEdgeMesh->openEdgeList[j];
    //         Vector3 a = edge.p1;
    //         Vector3 b = edge.p2;
    //         Vector3 center = Vector3Scale(Vector3Add(a, b), 0.5f);
    //         Color color = (Color){0, progress, 255-progress, 255};
    //         GraphicsDrawArrow(camera, a, b, 0.03f, color);
    //         // DrawLine3D(a, b, color);
    //         // DrawLine3D(center, Vector3Add(center, Vector3Scale(edge.triangleNormal, 0.1f)), RED);
    //         // for (int k=0;k<6;k++)
    //         // {
    //         //     if (edge.planeIndices & (1 << k))
    //         //     {
    //         //         EdgePlane plane = openEdgeMesh->edgePlanes[k];
    //         //         Vector3 center = Vector3Scale(Vector3Add(a, b), 0.5f);
    //         //         // Vector3 normal = Vector3Normalize(Vector3Subtract(b, a));
    //         //         // Vector3 offset = Vector3Scale(normal, 0.01f);
    //         //         // DrawLine3D(center, Vector3Add(center, offset), color);
    //         //         float len = 0.1f; // PlaneDistance(center, plane);
    //         //         DrawLine3D(center, Vector3Add(center, Vector3Scale(plane.normal, len)), plane.color);
    //         //     }
    //         // }
    //     }
    //     rlPopMatrix();
    //     // rlTranslatef(0.0f, 1.2f, 0.0f);
    //     if (i % rows == rows - 1)
    //         rlTranslatef(gap, 0.0f, -gap * (rows - 1.0f));
    //     else
    //         rlTranslatef(0.0f, 0.0f, gap);
    // }
    // rlPopMatrix();
    // rlDrawRenderBatchActive();
    // rlEnableDepthTest();

    EndMode3D();
}