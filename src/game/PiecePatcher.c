#include "game.h"
#include <raymath.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <rlgl.h>
#include "shared/util/graphics.h"
#include "PiecePatcher.h"
// The idea is the following:
// Each mesh is a building piece that can be rotated and placed in a grid of 1x1x1 cubes.
// Two pieces can be stitched together if their border edges match.
//
static Camera _camera;

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
    float edgeLength;
    // bitmask of planes that this edge is on
    uint8_t planeIndices;
    // could be an index too, but this way the bit rotation works too
    uint8_t cornerBits1;
    uint8_t cornerBits2;
    uint8_t startPointMarker1;
    uint8_t startPointMarker2;
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
    int index;
    int marker;
    int compatible;
    int compatibleCount;
    int compatibleOrderId;
} OpenEdgeMesh;

static OpenEdgeMesh *_openEdgeMeshes;
static int _openEdgeMeshCount;
static int _openEdgeMeshCapacity;
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

static Model _buildingPieces;
static uint8_t *_edgeMarkers;
static int _edgeMarkerCapacity;
static int *_verticeMergeMap;
static int _verticeMergeMapCapacity;
static int *_verticeTriangleCounter;
static int _verticeTriangleCounterCapacity;

static int OpenEdgeMesh_isMatch(OpenEdgeMesh *mesh, OpenEdgeMesh *other, int offsetX, int offsetY, int offsetZ);

Mesh *PiecePatcher_getMesh(int index)
{
    return _openEdgeMeshes[index].mesh;
}

Material *PiecePatcher_getMaterial(int index)
{
    return &_buildingPieces.materials[1];
}

int PiecePatcher_getRotation(int index)
{
    return _openEdgeMeshes[index].rotation;
}

static int PiecePatcher_determineCompatibility(int i, PlacedPiece *pieces, int pieceCount, int validFlagMask)
{
    int marker = (int) ((void*)pieces) ^ pieceCount;
    if (_openEdgeMeshes[i].marker != marker)
    {
        int compatible = 0;
        int index = 0;
        int compatibleCount = 0;
        while (index < pieceCount)
        {
            PlacedPiece *piece = &pieces[index];
            int isCompatible = (OpenEdgeMesh_isMatch(&_openEdgeMeshes[i], &_openEdgeMeshes[piece->pieceIndex], piece->offsetX, piece->offsetY, piece->offsetZ) & validFlagMask);
            if (isCompatible)
            {
                compatibleCount++;
            }
            compatible = compatible | isCompatible;
            index++;
        }
        _openEdgeMeshes[i].compatible = compatible;
        _openEdgeMeshes[i].marker = marker;
        _openEdgeMeshes[i].compatibleCount = compatibleCount;
        _openEdgeMeshes[i].compatibleOrderId = -1;
    }
    return _openEdgeMeshes[i].compatible;
}

int PiecePatcher_getCompatibleMeshCountEx(PlacedPiece *pieces, int pieceCount, int validFlagMask)
{
    int result = 0;
    PlacedPiece *pieceList[256];

    for (int i = 0; i < _openEdgeMeshCount; i++)
    {
        if (PiecePatcher_determineCompatibility(i, pieces, pieceCount, validFlagMask))
        {
            result++;
        }
    }
    printf("Compatible count %d\n", result);
    // for (int i=0;i<pieceCount;i++)
    // {
    //     printf(" Piece %d %d %d %d\n", pieces[i].pieceIndex, pieces[i].offsetX, pieces[i].offsetY, pieces[i].offsetZ);
    // }
    return result;
}

int PiecePatcher_getCompatibleMeshByIndexEx(PlacedPiece *pieces, int pieceCount, int validFlagMask, int compatibleIndex)
{
    for (int i = 0; i < _openEdgeMeshCount; i++)
    {
        if (PiecePatcher_determineCompatibility(i, pieces, pieceCount, validFlagMask))
        {
            if (compatibleIndex == 0)
            {
                return i;
            }
            compatibleIndex--;
        }
    }

    return -1;
}

int PiecePatcher_getCompatibleMeshCount(int index, int offsetX, int offsetY, int offsetZ)
{
    int marker = index << 11 ^ (offsetX << 8) ^ (offsetY << 4) ^ offsetZ;
    int result = 0;
    for (int i = 0; i < _openEdgeMeshCount; i++)
    {
        if (_openEdgeMeshes[i].marker != marker)
        {
            _openEdgeMeshes[i].compatible = OpenEdgeMesh_isMatch(&_openEdgeMeshes[index], &_openEdgeMeshes[i], offsetX, offsetY, offsetZ);
            _openEdgeMeshes[i].marker = marker;
        }
        if (_openEdgeMeshes[i].compatible)
        {
            result++;
        }
    }
    return result;
}

int PiecePatcher_getCompatibleMeshByIndex(int index, int offsetX, int offsetY, int offsetZ, int compatibleIndex)
{
    int marker = index << 11 ^ (offsetX << 8) ^ (offsetY << 4) ^ offsetZ;
    for (int i = 0; i < _openEdgeMeshCount; i++)
    {
        if (_openEdgeMeshes[i].marker != marker)
        {
            _openEdgeMeshes[i].compatible = OpenEdgeMesh_isMatch(&_openEdgeMeshes[index], &_openEdgeMeshes[i], offsetX, offsetY, offsetZ);
            _openEdgeMeshes[i].marker = marker;
        }
        if (_openEdgeMeshes[i].compatible)
        {
            if (compatibleIndex == 0)
            {
                return i;
            }
            compatibleIndex--;
        }
    }
    return -1;
}


static int getCornerIndex(Vector3 corner)
{
    for (int i = 0; i < 8; i++)
    {
        if (Vector3DistanceSqr(_cubeCorners[i], corner) < 0.0001f)
        {
            return i;
        }
    }
    return -1;
}

static void OpenEdgeMesh_ensureCapacity(int capacity)
{
    if (_openEdgeMeshCapacity < capacity)
    {
        _openEdgeMeshCapacity = capacity;
        _openEdgeMeshes = STRUCT_REALLOC(_openEdgeMeshes, OpenEdgeMesh, _openEdgeMeshCapacity);
    }
}

static OpenEdgeMesh *OpenEdgeMesh_allocate(Mesh *mesh)
{
    OpenEdgeMesh_ensureCapacity(_openEdgeMeshCount + 1);
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

static void OpenEdgeMesh_tryAddEdge(OpenEdgeMesh *mesh, int v1, int v2, int v3, int hasMarkedCorners)
{
    Vector3 a = (Vector3){mesh->mesh->vertices[v1 * 3], mesh->mesh->vertices[v1 * 3 + 1], mesh->mesh->vertices[v1 * 3 + 2]};
    Vector3 b = (Vector3){mesh->mesh->vertices[v2 * 3], mesh->mesh->vertices[v2 * 3 + 1], mesh->mesh->vertices[v2 * 3 + 2]};
    MeshEdge edge = (MeshEdge){.v1 = v1, .v2 = v2, .planeIndices = 0};
    for (int i = 0; i < 6; i++)
    {
        EdgePlane plane = mesh->edgePlanes[i];
        Vector3 center = Vector3Scale(Vector3Add(a, b), 0.5f);

        // consider any edge that is near the cube boundaries or outside as a border edge
        if (PlaneDistance(center, plane) > -0.15f)
        {
            edge.planeIndices |= 1 << i;
        }
    }

    // if (mesh->mesh->texcoords2 && mesh->mesh->texcoords2[v1 * 2] > 0.5f)
    // {
    //     printf("Found edge with texcoord %f %f\n", mesh->mesh->texcoords2[v1 * 2], mesh->mesh->texcoords2[v1 * 2 + 1]);
    // }

    if (edge.planeIndices != 0)
    {
        Vector3 c = (Vector3){mesh->mesh->vertices[v3 * 3], mesh->mesh->vertices[v3 * 3 + 1], mesh->mesh->vertices[v3 * 3 + 2]};
        int cornerIndexA = getCornerIndex(a);
        int cornerIndexB = getCornerIndex(b);
        edge.cornerBits1 = cornerIndexA == -1 || hasMarkedCorners ? 0 : 1 << cornerIndexA;
        edge.cornerBits2 = cornerIndexB == -1 || hasMarkedCorners ? 0 : 1 << cornerIndexB;
        edge.startPointMarker1 = _verticeTriangleCounter[v1] < 0 || edge.cornerBits1 != 0;
        edge.startPointMarker2 = _verticeTriangleCounter[v2] < 0 || edge.cornerBits2 != 0;
        edge.p1 = a;
        edge.p2 = b;
        edge.edgeLength = Vector3Distance(a, b);
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

static bool MathIsPointOnLineSegment(Vector3 point, Vector3 a, Vector3 b, int dbg)
{
    if (dbg)
    {
        DrawCube(point, 0.01f, 0.01f, 0.01f, RED);
        DrawLine3D(a, b, BLACK);
    }

    Vector3 d = Vector3Subtract(b, a);
    Vector3 e = Vector3Subtract(point, a);

    float q = Vector3DotProduct(d, e);
    if (q < -0.0001f)
        return false;
    float dote = Vector3DotProduct(e, e);
    if (dote < 0.000001f)
    {
        if (dbg) DrawCubeWires(a, 0.15f, 0.15f, 0.15f, WHITE);
        return true;
    }
    // q *= q;
    float dotd = Vector3DotProduct(d, d); // Dot product of d with itself, equivalent to |d|^2
    q /= dotd;

    if (q > 1.0001f)
        return false;
    
    float px = a.x + q * d.x;
    float py = a.y + q * d.y;
    float pz = a.z + q * d.z;
    float dpx = px - point.x;
    float dpy = py - point.y;
    float dpz = pz - point.z;
    float distance = dpx * dpx + dpy * dpy + dpz * dpz;
    if (distance > 0.000001f)
    {
        
        if (dbg)
        {
            // printf("Distance %f %.2f | %.2f %.2f %.2f -> %.2f %.2f %.2f : %.2f %.2f %.2f  q=%.4f (%.3f %.3f %.3f)\n", distance, q, a.x, a.y, a.z, b.x, b.y, b.z, point.x, point.y, point.z, q, px, py, pz);
            DrawCube(point, 0.01f, 0.2f, 0.01f, RED);
            DrawCube((Vector3){px,py,pz}, 0.2f, 0.01f, 0.01f, RED);
            DrawLine3D(a, b, RED);
        }
        return false;
    }
    if (dbg)
    {
        DrawCubeWires((Vector3){px, py, pz}, 0.075f, 0.075f, 0.075f, GREEN);
        DrawCubeWires(point, 0.05f, 0.05f, 0.05f, RED);
    }
    return true;
}

static bool MeshEdge_hasMatchingPoints(MeshEdge edge, MeshEdge other, int offsetX, int offsetY, int offsetZ, float *overlap, int dbg)
{
    Vector3 a1 = edge.p1;
    Vector3 a2 = edge.p2;
    Vector3 b1 = Vector3Add(other.p1, (Vector3){offsetX, offsetY, offsetZ});
    Vector3 b2 = Vector3Add(other.p2, (Vector3){offsetX, offsetY, offsetZ});
    // if (dbg)
    // {
    //     DrawCubeWires(a1, 0.05f, 0.05f, 0.05f, RED);
    //     DrawCubeWires(a2, 0.05f, 0.05f, 0.05f, RED);
    //     DrawCubeWires(b1, 0.035f, 0.035f, 0.035f, BLUE);
    //     DrawCubeWires(b2, 0.035f, 0.035f, 0.035f, BLUE);
    // }
    // check if a1 is on the line b1-b2 or b1 is on the line a1-a2 and vice versa
    // return MathIsPointOnLineSegment(a2, b1, b2);

    bool isA1OnB = MathIsPointOnLineSegment(a1, b1, b2, dbg == 2);
    bool isA2OnB = MathIsPointOnLineSegment(a2, b1, b2, dbg == 2);
    bool isB1OnA = MathIsPointOnLineSegment(b1, a1, a2, dbg == 2);
    bool isB2OnA = MathIsPointOnLineSegment(b2, a1, a2, dbg == 2);
    if (overlap && (isA1OnB || isB2OnA) && (isA2OnB || isB1OnA))
    {
        if (isA1OnB && isA2OnB)
        {
            *overlap = edge.edgeLength;
            if (dbg) GraphicsDrawArrow(_camera, a1, a2, 0.01f, RED);
        }
        else if (isB1OnA && isB2OnA)
        {
            *overlap = other.edgeLength;
            if (dbg) GraphicsDrawArrow(_camera, b1, b2, 0.01f, BLUE);
        }
        else
        {
            Vector3 p1 = isA1OnB ? a1 : b2;
            Vector3 p2 = isA2OnB ? a2 : b1;
            *overlap = Vector3Distance(p1, p2);
            if (dbg)
            {
                DrawLine3D(a1, Vector3Add(a1, (Vector3){0, 0.2f, 0}), RED);
                DrawLine3D(a2, Vector3Add(a2, (Vector3){0, 0.2f, 0}), RED);
                DrawLine3D(b1, Vector3Add(b1, (Vector3){0, 0.2f, 0}), BLUE);
                DrawLine3D(b2, Vector3Add(b2, (Vector3){0, 0.2f, 0}), BLUE);
                GraphicsDrawArrow(_camera, p1, p2, 0.01f, GREEN);
            }
        }
    }
    // return (MathIsPointOnLineSegment(a1, b1, b2, dbg) || MathIsPointOnLineSegment(b2, a1, a2, dbg)) && (MathIsPointOnLineSegment(a2, b1, b2, dbg) || MathIsPointOnLineSegment(b1, a1, a2, dbg));

    return (isA1OnB || isB2OnA) && (isA2OnB || isB1OnA);
}


#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  ((byte) & 0x80 ? '1' : '0'), \
  ((byte) & 0x40 ? '1' : '0'), \
  ((byte) & 0x20 ? '1' : '0'), \
  ((byte) & 0x10 ? '1' : '0'), \
  ((byte) & 0x08 ? '1' : '0'), \
  ((byte) & 0x04 ? '1' : '0'), \
  ((byte) & 0x02 ? '1' : '0'), \
  ((byte) & 0x01 ? '1' : '0') 

void OpenEdgeMesh_drawCornerBitset(uint8_t cornerBits, Vector3 position, float size, Color color)
{
    for (int i = 0; i < 8; i++)
    {
        if (cornerBits & (1 << i))
        {
            DrawCubeWires(Vector3Add(position, _cubeCorners[i]), size, size, size, color);
        }
    }
}

static int OpenEdgeMesh_mapCornerBitSet(uint8_t cornerBits, int offsetX, int offsetY, int offsetZ)
{
    if (offsetX == 0 && offsetY == 0 && offsetZ == 0)
    {
        return cornerBits;
    }

    int result = 0;
    for (int i = 0; i < 8; i++)
    {
        if (cornerBits & (1 << i))
        {
            Vector3 corner = Vector3Add(_cubeCorners[i], (Vector3){offsetX, offsetY, offsetZ});
            int cornerIndex = getCornerIndex(corner);
            if (cornerIndex == -1)
            {
                // TraceLog(LOG_WARNING, "Corner not found %d %d %d", (int)corner.x, (int)corner.y, (int)corner.z);
            }
            else
            {
                result |= 1 << cornerIndex;
            }
        }
    }

    return result;
}

static int NormalToBitFlag(Vector3 dir)
{
    float ax = fabsf(dir.x);
    float ay = fabsf(dir.y);
    float az = fabsf(dir.z);
    if (ax > ay && ax > az)
    {
        return dir.x > 0 ? 1 : 2;
    }
    if (ay > ax && ay > az)
    {
        return dir.y > 0 ? 4 : 8;
    }
    return dir.z > 0 ? 16 : 32;
}

static int HighestBitOf(int bits)
{
    int result = 0;
    while (bits > 0)
    {
        bits >>= 1;
        result++;
    }
    return result;
}

static int OpenEdgeMesh_getStartingPointIndex(OpenEdgeMesh *mesh, Vector3 position)
{
    for (int i = 0; i < mesh->openEdgeCount; i++)
    {
        MeshEdge edge = mesh->openEdgeList[i];
        if (Vector3DistanceSqr(edge.p1, position) < 0.0001f)
        {
            return i;
        }
    }
    return -1;
}

static bool OpenMeshEdge_allEdgePointsOnLineSegment(OpenEdgeMesh *a, int startIndexA, int endIndexA, int directionA,
    OpenEdgeMesh *b, int startIndexB, int endIndexB, int directionB, Vector3 offset, int dbg)
{
    for (int i = startIndexA; ; )
    {
        Vector3 pointA = a->openEdgeList[i].p1;
        // if (dbg) DrawCube(pointA, 0.01f, 0.01f, 0.01f, RED);
        int j = startIndexB;
        while(1)
        {
            MeshEdge edgeB = b->openEdgeList[j];
            if (MathIsPointOnLineSegment(pointA, Vector3Add(edgeB.p1, offset), Vector3Add(edgeB.p2, offset), 0))
            {
                goto pointIsOnLine;
            }
            if (j == endIndexB) break;
            // if (dbg)
            // {
            //     DrawCubeWires(Vector3Add(edgeB.p1, offset), 0.01f, 0.01f, 0.01f, BLUE);
            //     DrawLine3D(Vector3Add(edgeB.p1, offset), Vector3Add(edgeB.p2, offset), BLUE);
            // }
            j = (j + b->openEdgeCount + directionB) % b->openEdgeCount;
        }
        if (dbg)
        {
            DrawCube(pointA, 0.01f, 0.01f, 0.01f, RED);
            MeshEdge edgeB = b->openEdgeList[j];
            DrawLine3D(Vector3Add(edgeB.p1, offset), Vector3Add(edgeB.p2, offset), BLUE);
        }
        return false;
        pointIsOnLine:
        if (i == endIndexA) break;
        i = (i + a->openEdgeCount + directionA) % a->openEdgeCount;
    }
    return true;
}

static int OpenEdgeMesh_isMatch(OpenEdgeMesh *mesh, OpenEdgeMesh *other, int offsetX, int offsetY, int offsetZ)
{
    int offsetDist = offsetX * offsetX + offsetY * offsetY + offsetZ * offsetZ;
    if (offsetDist > 1)
    {
        return MESH_TILE_MATCH_DISJUNCT;
    }
    int matchCount = 0;
    // int cornersToMatchMask = mesh->cubeCornerSet & OpenEdgeMesh_mapCornerBitSet(other->cubeCornerSet, offsetX, offsetY, offsetZ);
    
    // OpenEdgeMesh_drawCornerBitset(mesh->cubeCornerSet, (Vector3){0, 0, 0}, 0.1f, RED);
    // OpenEdgeMesh_drawCornerBitset(other->cubeCornerSet, (Vector3){offsetX, offsetY, offsetZ}, 0.075f, BLUE);
    // OpenEdgeMesh_drawCornerBitset(cornersToMatchMask, (Vector3){0, 0, 0}, 0.125f, YELLOW);

    // check if corner normals are opposing
    int cornerMatchCount = 0;
    for (int i=0;i<8;i++)
    {
        Vector3 corner = _cubeCorners[i];
        Vector3 normalA[8];
        Vector3 normalB[8];
        int normalACount = 0;
        int normalBCount = 0;
        // find normal of corner; if there are multiple normals, the edge is a valid match in any case
        for (int j = 0; j < mesh->openEdgeCount; j++)
        {
            MeshEdge edge = mesh->openEdgeList[j];

            if (Vector3DistanceSqr(edge.p1, corner) < 0.0001f || Vector3DistanceSqr(edge.p2, corner) < 0.0001f)
            {
                normalA[normalACount++] = edge.triangleNormal;
            }
        }

        for (int j = 0; j < other->openEdgeCount; j++)
        {
            MeshEdge edge = other->openEdgeList[j];

            if (Vector3DistanceSqr(Vector3Add(edge.p1, (Vector3){offsetX, offsetY, offsetZ}), corner) < 0.0001f || 
                Vector3DistanceSqr(Vector3Add(edge.p2, (Vector3){offsetX, offsetY, offsetZ}), corner) < 0.0001f)
            {
                normalB[normalBCount++] = edge.triangleNormal;
            }
        }

        for (int j = 0; j < normalACount; j++)
        {
            for (int k = 0; k < normalBCount; k++)
            {
                if (Vector3DotProduct(normalA[j], normalB[k]) < -0.5f)
                {
                    // printf("Corner %d mismatch %d %d\n", i, j, k);
                    return MESH_TILE_MATCH_MISMATCH;
                }
            }
        }
    }
    
    Vector3 offset = (Vector3){offsetX, offsetY, offsetZ};

    // loopedeloop, you're in for a ride!
    // foreach starting point in "mesh" ...
    for (int i = 0; i < mesh->openEdgeCount; i++)
    {
        MeshEdge edge = mesh->openEdgeList[i];
        if (edge.startPointMarker1 == 0)
            continue;
        Vector3 startA = edge.p1;
        // ... we look for a matching starting point in "other"
        for (int j = 0; j < other->openEdgeCount; j++)
        {
            MeshEdge otherEdge = other->openEdgeList[j];
            if (otherEdge.startPointMarker2 == 0 || Vector3DistanceSqr(edge.p1, Vector3Add(otherEdge.p2, offset)) > 0.0001f)
                continue;
            cornerMatchCount += 1;
            Vector3 startB = Vector3Add(otherEdge.p2, offset);

            DrawCube(startA, 0.01f, 0.01f, 0.01f, RED);
            
            // having found a pair, we can now start comparing:
            // What now follows is not insanity but a search for the end points of both edge loops
            for (int k = 0; k < mesh->openEdgeCount; k++)
            {
                // iterating here over "other" for a end point
                int indexK = (j + k) % other->openEdgeCount;
                MeshEdge otherNextEdge = other->openEdgeList[indexK];
                if (otherNextEdge.startPointMarker1 == 0)
                    continue;
                Vector3 endB = Vector3Add(otherNextEdge.p1, offset);
                // endB is found, now let's search for endA (end point of "mesh")
                for (int l = 1; l < mesh->openEdgeCount; l++)
                {
                    int indexL = (i - l + mesh->openEdgeCount) % mesh->openEdgeCount;
                    MeshEdge nextEdge2 = mesh->openEdgeList[indexL];

                    if (nextEdge2.startPointMarker1 == 0)
                        continue;
                    Vector3 endA = nextEdge2.p1;
                    // we've found two endpoints - do they match?
                    if (Vector3DistanceSqr(endA, endB) < 0.0001f)
                    {
                        // yes! So let's compare the loops individually if all their points are located on the other loop
                        if (!OpenMeshEdge_allEdgePointsOnLineSegment(mesh, i, indexL, -1, other, j, indexK, 1, offset, 0) ||
                            !OpenMeshEdge_allEdgePointsOnLineSegment(other, j, indexK, 1, mesh, i, indexL, -1, Vector3Negate(offset), 0))
                        {
                            return MESH_TILE_MATCH_MISMATCH;
                        }
                        matchCount++;
                    }
                    else if (
                        OpenEdgeMesh_getStartingPointIndex(mesh, endB) >= 0 || 
                        OpenEdgeMesh_getStartingPointIndex(other, Vector3Subtract(endA, offset)) >= 0)
                    {
                        
                        // GraphicsDrawArrow(_camera, startA, endA, 0.01f, RED);
                        // GraphicsDrawArrow(_camera, startB, endB, 0.01f, BLUE);
                        return MESH_TILE_MATCH_MISMATCH;
                    }
                    break;
                }
                break;
            }
        }
    }

    if (matchCount > 0)
    {
        return MESH_TILE_MATCH_COMPATIBLE;
    } 

    if (cornerMatchCount == 0)
    {
        return offsetDist == 0 ? MESH_TILE_MATCH_CONGRUENTLY_DISJUNCT : MESH_TILE_MATCH_DISJUNCT;
    }

    if (cornerMatchCount > 1)
    {
        return MESH_TILE_MATCH_MISMATCH;
    }

    // printf("Corner match count %d\n", cornerMatchCount);

    return offsetDist > 1 ? MESH_TILE_MATCH_PARTIALLY_CONGRUENTLY_DISJUNCT : MESH_TILE_MATCH_PARTIALLY_DISJUNCT;
}

static int GetMappedVertexIndex(int index)
{
    return _verticeMergeMap[index] == -1 ? index : _verticeMergeMap[index];
}

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
    if (_verticeTriangleCounterCapacity < vertexCount)
    {
        _verticeTriangleCounter = STRUCT_REALLOC(_verticeTriangleCounter, int, vertexCount);
        _verticeTriangleCounterCapacity = vertexCount;
    }
    memset(_verticeTriangleCounter, 0, vertexCount * sizeof(int));

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
    // Vector3 a = (Vector3){vertices[v1 * 3], vertices[v1 * 3 + 1], vertices[v1 * 3 + 2]};
    // Vector3 b = (Vector3){vertices[v2 * 3], vertices[v2 * 3 + 1], vertices[v2 * 3 + 2]};
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
        MeshEdge newEdge = (MeshEdge){
            .v1 = edge.v2, 
            .v2 = edge.v1, 
            .planeIndices = edge.planeIndices,
            .startPointMarker1 = edge.startPointMarker1,
            .startPointMarker2 = edge.startPointMarker2,
            .cornerBits1 = rotateBitsCW(edge.cornerBits1), 
            .cornerBits2 = rotateBitsCW(edge.cornerBits2)};
        if ((newEdge.cornerBits1 == 0) ^ (edge.cornerBits1 == 0))
        {
            TraceLog(LOG_WARNING, "Corner bits 1 zeroed of %s@%d %d %d", mesh->mesh->name, mesh->rotation, newEdge.cornerBits1, edge.cornerBits1);
        }

        if ((newEdge.cornerBits2 == 0) ^ (edge.cornerBits2 == 0))
        {
            TraceLog(LOG_WARNING, "Corner bits 2 zeroed of %s@%d %d %d", mesh->mesh->name, mesh->rotation, newEdge.cornerBits2, edge.cornerBits2);
        }

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
    for (int i = 0; i < newMesh->openEdgeCount; i++)
    {
        MeshEdge edge = newMesh->openEdgeList[i];
        if (edge.cornerBits1 < 0)
        {
            TraceLog(LOG_WARNING, "Corner bits 1 [%d] negative of %s@%d %d (was %d)", i, mesh->mesh->name, mesh->rotation, edge.cornerBits1,
                mesh->openEdgeList[i].cornerBits1);
        }
        if (edge.cornerBits2 < 0)
        {
            TraceLog(LOG_WARNING, "Corner bits 2 [%d] negative of %s@%d %d (was %d)", i, mesh->mesh->name, mesh->rotation, edge.cornerBits2,
                mesh->openEdgeList[i].cornerBits2);
        }
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
        int a = mesh->indices[i];
        int b = mesh->indices[i + 1];
        int c = mesh->indices[i + 2];
        MarkEdge(edgeCount, a, b);
        MarkEdge(edgeCount, b, c);
        MarkEdge(edgeCount, c, a);
        _verticeTriangleCounter[GetMappedVertexIndex(a)]++;
        _verticeTriangleCounter[GetMappedVertexIndex(b)]++;
        _verticeTriangleCounter[GetMappedVertexIndex(c)]++;
    }
    
    int markedCornerCount = 0;
    // flag vertex index of single triangle marked vertices with a negative counter so we can
    // flag them as corner points
    for (int i = 0; i < triangleCount; i += 3)
    {
        int a = GetMappedVertexIndex(mesh->indices[i]);
        int b = GetMappedVertexIndex(mesh->indices[i + 1]);
        int c = GetMappedVertexIndex(mesh->indices[i + 2]);
        int countA = _verticeTriangleCounter[a];
        int countB = _verticeTriangleCounter[b];
        int countC = _verticeTriangleCounter[c];
        if ((countA == 1 && countB == 1) || (countB == 1 && countC == 1) || (countC == 1 && countA == 1))
        {
            if (countA > 1) _verticeTriangleCounter[a] = -countA, markedCornerCount++;
            if (countB > 1) _verticeTriangleCounter[b] = -countB, markedCornerCount++;
            if (countC > 1) _verticeTriangleCounter[c] = -countC, markedCornerCount++;
        }

    }

    for (int i = 0; i < triangleCount; i += 3)
    {
        int a = GetMappedVertexIndex(mesh->indices[i]);
        int b = GetMappedVertexIndex(mesh->indices[i + 1]);
        int c = GetMappedVertexIndex(mesh->indices[i + 2]);
        int countA = _verticeTriangleCounter[a];
        int countB = _verticeTriangleCounter[b];
        int countC = _verticeTriangleCounter[c];
        if ((countA == 1 && countB == 1) || (countB == 1 && countC == 1) || (countC == 1 && countA == 1))
        {
            // skip triangle markers
            continue;
        }

        bool borderA = IsBorderEdge(edgeCount, a, b, mesh->vertices);
        bool borderB = IsBorderEdge(edgeCount, b, c, mesh->vertices);
        bool borderC = IsBorderEdge(edgeCount, c, a, mesh->vertices);
        if (borderA)
            OpenEdgeMesh_tryAddEdge(openEdgeMesh, a, b, c, markedCornerCount);
        if (borderB)
            OpenEdgeMesh_tryAddEdge(openEdgeMesh, b, c, a, markedCornerCount);
        if (borderC)
            OpenEdgeMesh_tryAddEdge(openEdgeMesh, c, a, b, markedCornerCount);

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
            if (Vector3DistanceSqr(edge.p1, other.p2) < 0.0001f)
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
        OpenEdgeMesh_ensureCapacity(_openEdgeMeshCount + 4);
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
    // for (int j = 0; j < 8; j++)
    // {
    //     if (openEdgeMesh->cubeCornerSet & (1 << j))
    //     {
    //         DrawCubeWires(_cubeCorners[j], 0.1f, 0.1f, 0.1f, GREEN);
    //     }
    // }
    for (int j = 0; j < openEdgeMesh->openEdgeCount; j++)
    {
        int progress = 255 * j / (openEdgeMesh->openEdgeCount - 1);
        MeshEdge edge = openEdgeMesh->openEdgeList[j];
        Vector3 a = edge.p1;
        Vector3 b = edge.p2;
        Vector3 center = Vector3Scale(Vector3Add(a, b), 0.5f);
        Color color = (Color){0, progress, 255 - progress, 255};
        GraphicsDrawArrow(camera, a, b, 0.005f, color);
        if (edge.startPointMarker1)
        {
            DrawCubeWires(a, 0.1f, 0.1f, 0.1f, WHITE);
        }
        if (edge.startPointMarker2)
        {
            DrawCubeWires(b, 0.075f, 0.075f, 0.075f, LIGHTGRAY);
        }
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

static int TestMeshes(Vector3 position, OpenEdgeMesh *meshA, OpenEdgeMesh *meshB, int offsetX, int offsetY, int offsetZ, int outlines, int expectedResult)
{
    int returnCode = 0;
    rlPushMatrix();
    rlTranslatef(position.x, position.y, position.z);

    rlPushMatrix();
    rlRotatef(90.0f * -meshA->rotation, 0.0f, 1.0f, 0.0f);
    DrawMesh(*meshA->mesh, _buildingPieces.materials[1], MatrixIdentity());
    rlPopMatrix();

    Material mat = _buildingPieces.materials[1];
    int result = OpenEdgeMesh_isMatch(meshA, meshB, offsetX, offsetY, offsetZ);
    if (result != expectedResult)
    {
        int fade = (int)((sinf(GetTime() * 4.0f) * .125f + .875f) * 255.0f);
        mat.maps[MATERIAL_MAP_DIFFUSE].color = (Color){255, fade, fade, 255};
        returnCode = 1;
    }
    else if (result == MESH_TILE_MATCH_MISMATCH)
    {
        mat.maps[MATERIAL_MAP_DIFFUSE].color = (Color){255, 255, 255, 64};
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
    
    if (outlines)
    {
        if (outlines & 1) OpenEdgeMesh_drawDebug(_camera, meshA);
        rlTranslatef(offsetX, offsetY, offsetZ);
        if (outlines & 2) OpenEdgeMesh_drawDebug(_camera, meshB);
    }
    
    rlPopMatrix();
    rlDrawRenderBatchActive();
    rlEnableDepthTest();

    rlPopMatrix();
    return returnCode;
}

void PlacedPieces_shift(PlacedPiece *pieces, int pieceCount, int dx, int dy, int dz)
{
    for (int i = 0; i < pieceCount; i++)
    {
        pieces[i].offsetX += dx;
        pieces[i].offsetY += dy;
        pieces[i].offsetZ += dz;
    }
}

void PiecePatcherDrawDebug()
{
    static float zoom = 1.0f;
    int debug = 0;

    float wheel = GetMouseWheelMove();
    if (wheel != 0)
    {
        zoom += wheel * 0.1f;
        zoom = fmaxf(0.1f, zoom);
        zoom = fminf(10.0f, zoom);
    }

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
        .fovy = 16.0f * zoom,
        .position = camPos,
        .target = camTarget,
        .up = (Vector3){0.0f, 1.0f, 0.0f},
        .projection = CAMERA_ORTHOGRAPHIC};
    BeginMode3D(camera);
    static float rotate = 0.0f;
    if (IsKeyDown(KEY_A))
    {
        rotate -= GetFrameTime() * 100.0f;
    }
    if (IsKeyDown(KEY_D))
    {
        rotate += GetFrameTime() * 100.0f;
    }
    rlPushMatrix();

    rlRotatef(rotate, 0, 1, 0);
    

    Matrix identity = MatrixIdentity();
    // DrawCubeWires((Vector3){0, 0.5f, 0}, 1.0f, 1.0f, 1.0f, RED);
    Material mat = _buildingPieces.materials[1];
    mat.maps[MATERIAL_MAP_DIFFUSE].color = (Color){255, 255, 255, 255};
    _camera = camera;

    static int mode = 4;
    static int outlines = 0;
    if (IsKeyPressed(KEY_O))
    {
        outlines = (outlines + 1) % 4;
    }
    if (IsKeyPressed(KEY_M))
    {
        debug = 1;
        mode = (mode + 1) % 5;
    }
    switch (mode)
    {
        case 0:
        {
        static int testIndex = 0;
        static int ox = 1, oy = 0, oz = 0;
        if (IsKeyPressed(KEY_ONE))
        {
            ox = 0;
            oy = 1;
            oz = 0;
            debug = 1;
        }
        if (IsKeyPressed(KEY_TWO))
        {
            ox = 0;
            oy = 0;
            oz = 1;
            debug = 1;
        }
        if (IsKeyPressed(KEY_THREE))
        {
            ox = 1;
            oy = 0;
            oz = 0;
            debug = 1;
        }
        if (IsKeyPressed(KEY_FOUR))
        {
            ox = 0;
            oy = -1;
            oz = 0;
            debug = 1;
        }
        if (IsKeyPressed(KEY_FIVE))
        {
            ox = 0;
            oy = 0;
            oz = 0;
            debug = 1;
        }

        if (IsKeyPressed(KEY_LEFT))
        {
            testIndex = (testIndex - 4 + _openEdgeMeshCount) % _openEdgeMeshCount;
            debug = 1;
        }
        if (IsKeyPressed(KEY_RIGHT))
        {
            testIndex = (testIndex + 4) % _openEdgeMeshCount;
            debug = 1;
        }


        OpenEdgeMesh *meshA = &_openEdgeMeshes[testIndex];
        if (debug) printf("Testing %s (%d) %d %d %d\n", meshA->mesh->name, testIndex, ox, oy, oz);
        int x = 0, z = 0;
        int matches = 0;
        for (int i = 0; i < _openEdgeMeshCount; i++)
        {
            OpenEdgeMesh *meshB = &_openEdgeMeshes[i];
            if (OpenEdgeMesh_isMatch(meshA, meshB, ox, oy, oz) == MESH_TILE_MATCH_COMPATIBLE)
            {
                matches++;
                if (debug) printf("  Match %s (%d)\n", meshB->mesh->name, i);
                TestMeshes((Vector3){x * 2.5f - 4.0f, 0.0f, z++ * 2.5f - 4.0f}, meshA, meshB, ox, oy, oz, outlines, MESH_TILE_MATCH_COMPATIBLE);
                if (z == 4)
                {
                    z = 0;
                    x++;
                }
            }
        }
        if (matches == 0)
        {
            DrawMesh(*meshA->mesh, mat, identity);
        }
    }
    break;
    case 1:
    {
        int testcases[] = {
            // 91, 87, -1, 0, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 84, 70, 0, -1, 0, MESH_TILE_MATCH_MISMATCH,
            95, 87, 0, 0, 1, MESH_TILE_MATCH_COMPATIBLE,
            159, 68, 0, -1, 0, MESH_TILE_MATCH_COMPATIBLE,
            68, 159, 0, 1, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 68, 70, 1, 0, 0, MESH_TILE_MATCH_MISMATCH,
            // 0, 0, 0, -1, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 44, 32, 0, 0, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 76, 60, 0, 1, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 0, 56, 0, 0, -1, MESH_TILE_MATCH_COMPATIBLE,
            // 0, 3, 0, 0, 1, MESH_TILE_MATCH_MISMATCH,
            // 0, 3, 0, 0, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 0, 1, 0, 0, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 59, 57, 0, 0, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 56, 55, 0, 0, 0, MESH_TILE_MATCH_MISMATCH,
            // 56, 55, 0, 1, 0, MESH_TILE_MATCH_MISMATCH,
            // 0, 55, 0, 1, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 24, 28, 0, 1, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 28, 60, 0, 0, 0, MESH_TILE_MATCH_MISMATCH,
            // 28, 50, 0, 0, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 71, 55, 0, 1, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 55, 55, 1, 0, 0, MESH_TILE_MATCH_MISMATCH, // <- broken, not sure if valid, ignore for now
            // 88, 96, 1, 0, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 84, 96, 1, 0, 0, MESH_TILE_MATCH_MISMATCH,
            // 152, 0, 0, 0, 1, MESH_TILE_MATCH_DISJUNCT,
            // 152, 132, 1, 0, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 136, 104, 1, 0, 0, MESH_TILE_MATCH_DISJUNCT,
            // 120, 119, 1, 0, 0, MESH_TILE_MATCH_COMPATIBLE,
            // 119, 115, 0, 0, 1, MESH_TILE_MATCH_COMPATIBLE,
            // 164, 87, 1, 0, 0, MESH_TILE_MATCH_DISJUNCT,
        };

        float width = fminf(sizeof(testcases) / sizeof(int) / 6 - 1, 5) * 2.5f;
        float height = 2.5f * (sizeof(testcases) / sizeof(int) / 6 - 1);
        rlDisableBackfaceCulling();
        for (int i=0;i<sizeof(testcases)/sizeof(int);i+=6)
        {
            int n = i / 6;
            int code = TestMeshes((Vector3){(n % 6) * 2.5f - width * .5f, 0.0f, (n / 6) * 2.5f - height * .125f}, &_openEdgeMeshes[testcases[i]], &_openEdgeMeshes[testcases[i+1]], testcases[i+2], testcases[i+3], testcases[i+4], outlines, testcases[i+5]);
            if (code && debug)
            {
                printf("Test %d failed\n", n);
            }
        }
        rlEnableBackfaceCulling();
    }
    break;
    case 2:
    {
            
        PlacedPiece placedPieces[4*4*4*4];
        int placedPieceCount = 0;
        
        #define ADD_PIECE(x, y, z, index) { \
            placedPieces[placedPieceCount++] = (PlacedPiece){ \
                .offsetX = x, \
                .offsetY = y, \
                .offsetZ = z, \
                .pieceIndex = index, \
            }; \
            DrawCubeWires((Vector3){x, y + .5f, z}, 1.0f, 1.0f, 1.0f, GREEN); \
        }
        #define SHIFT(dx, dy, dz) for (int _i = 0; _i < placedPieceCount; _i++) { \
            placedPieces[_i].offsetX -= dx; \
            placedPieces[_i].offsetY -= dy; \
            placedPieces[_i].offsetZ -= dz; \
        }
        #define ADD_RANDOM_PIECE(x, y, z) \
            { \
                SHIFT(x, y, z);\
                int _count = PiecePatcher_getCompatibleMeshCountEx(placedPieces, placedPieceCount, MESH_TILE_MATCH_BUILD_MASK); \
                if (_count > 0) { \
                    int _select = GetRandomValue(0, _count - 1); \
                    int _index = PiecePatcher_getCompatibleMeshByIndexEx(placedPieces, placedPieceCount, MESH_TILE_MATCH_BUILD_MASK, _select); \
                    SHIFT(-x, -y, -z);\
                    ADD_PIECE(x, y, z, _index); \
                } \
                else {\
                    SHIFT(-x, -y, -z);\
                    DrawCubeWires((Vector3){x, y + .5f, z}, 1.0f, 1.0f, 1.0f, RED); \
                }\
            }
        #define ADD_RANDOM_PIECE_AT(x,y,z) \
            SHIFT(x, y, z); ADD_RANDOM_PIECE(); SHIFT(-x, -y, -z);


        static int seed = 122;
        SetRandomSeed(seed);

        if (IsKeyPressed(KEY_R))
        {
            seed++;
        }
        if (IsKeyPressed(KEY_T))
        {
            seed--;
        }

        
        ADD_PIECE(0, 0, 0, 1);
        ADD_RANDOM_PIECE(0, 0, 0);
        ADD_RANDOM_PIECE(0, 0, 0);
        ADD_RANDOM_PIECE(0, 0, 0);
        ADD_RANDOM_PIECE(0, 1, 0);
        // ADD_RANDOM_PIECE(1, 0, 0);
        // ADD_RANDOM_PIECE(0, 1, 0);
        // ADD_RANDOM_PIECE(0, 0, 1);
        // ADD_RANDOM_PIECE(1, 1, 0);
        // ADD_RANDOM_PIECE(1, 0, 1);
        // ADD_RANDOM_PIECE(0, 0, 0);
        // ADD_RANDOM_PIECE(0, 1, 0);
        // ADD_RANDOM_PIECE(1, 0, 0);
        for (int i = 0; i < placedPieceCount; i++)
        {
            PlacedPiece piece = placedPieces[i];
            OpenEdgeMesh *mesh = &_openEdgeMeshes[piece.pieceIndex];
            Matrix m = MatrixIdentity();
            m = MatrixMultiply(m, MatrixRotateY(PI * 0.5f * -mesh->rotation));
            m = MatrixMultiply(m,MatrixTranslate(piece.offsetX, piece.offsetY, piece.offsetZ));
            DrawMesh(*mesh->mesh, mat, m);
        }

        #undef ADD_PIECE
        #undef SHIFT
        #undef ADD_RANDOM_PIECE
        #undef ADD_RANDOM_PIECE_AT
    }
    break;
    case 3:
    {
        int tiles[] = {
            0, 0, 0, 0,
            1, 0, 0, 0,
            2, 0, 0, 0,
            3, 0, 0, 0,
        };
        int tileCount = sizeof(tiles) / sizeof(int) / 4;
        for (int i = 0; i < tileCount; i++)
        {
            int index = tiles[i * 4];
            int x = tiles[i * 4 + 1];
            int y = tiles[i * 4 + 2];
            int z = tiles[i * 4 + 3];
            OpenEdgeMesh *mesh = &_openEdgeMeshes[index];
            Color color = WHITE;
            for (int j = 0; j < i; j++)
            {
                int ox = tiles[j * 4 + 1];
                int oy = tiles[j * 4 + 2];
                int oz = tiles[j * 4 + 3];
                if (OpenEdgeMesh_isMatch(mesh, &_openEdgeMeshes[tiles[j * 4]], ox - x, oy - y, oz - z) == MESH_TILE_MATCH_MISMATCH)
                {
                    color = RED;
                    break;
                }
            }

            Matrix m = MatrixIdentity();
            m = MatrixMultiply(m, MatrixRotateY(PI * 0.5f * -mesh->rotation));
            m = MatrixMultiply(m, MatrixTranslate(x, y, z));

            mat.maps[MATERIAL_MAP_DIFFUSE].color = color;
            DrawMesh(*mesh->mesh, mat, m);
        }
        mat.maps[MATERIAL_MAP_DIFFUSE].color = WHITE;
    }
    break;
    case 4:
    {
        static int cursor[3] = {0, 0, 0};
        static PlacedPiece placedPieces[512];
        static int placedPieceCount = 0;

        DrawCubeWires((Vector3){cursor[0], cursor[1] + .5f, cursor[2]}, 1.0f, 1.0f, 1.0f, GREEN);
        if (IsKeyPressed(KEY_UP))
        {
            cursor[IsKeyDown(KEY_LEFT_SHIFT) ? 1 : 2]++;
            printf("Cursor %d %d %d\n", cursor[0], cursor[1], cursor[2]);
        }
        if (IsKeyPressed(KEY_DOWN))
        {
            cursor[IsKeyDown(KEY_LEFT_SHIFT) ? 1 : 2]--;
            printf("Cursor %d %d %d\n", cursor[0], cursor[1], cursor[2]);

        }
        if (IsKeyPressed(KEY_LEFT))
        {
            cursor[0]++;
            printf("Cursor %d %d %d\n", cursor[0], cursor[1], cursor[2]);

        }
        if (IsKeyPressed(KEY_RIGHT))
        {
            cursor[0]--;
            printf("Cursor %d %d %d\n", cursor[0], cursor[1], cursor[2]);
        }
        if (IsKeyPressed(KEY_F1))
        {
            // save
            SaveFileData("placedPieces.bin", placedPieces, placedPieceCount * sizeof(PlacedPiece));
        }
        if (IsKeyPressed(KEY_F2))
        {
            // load
            int dataSize;
            const char *data = LoadFileData("placedPieces.bin", &dataSize);
            if (data)
            {
                placedPieceCount = dataSize / sizeof(PlacedPiece);
                memcpy(placedPieces, data, dataSize);
                UnloadFileData(data);
            }
        }

        if (IsKeyPressed(KEY_DELETE))
        {
            for (int i= placedPieceCount - 1; i>=0;i--)
            {
                if (placedPieces[i].offsetX == cursor[0] && placedPieces[i].offsetY == cursor[1] && placedPieces[i].offsetZ == cursor[2])
                {
                    placedPieces[i] = placedPieces[--placedPieceCount];
                    break;
                }
            }
        }

        if (IsKeyPressed(KEY_B) || IsKeyPressed(KEY_N))
        {
            int direction = IsKeyPressed(KEY_N) ? -1 : 1;
            for (int i= placedPieceCount - 1; i>=0;i--)
            {
                if (placedPieces[i].offsetX == cursor[0] && placedPieces[i].offsetY == cursor[1] && placedPieces[i].offsetZ == cursor[2])
                {
                    PlacedPiece piece = placedPieces[i];
                    int currentIndex = piece.pieceIndex;
                    piece.pieceIndex = (piece.pieceIndex + direction + _openEdgeMeshCount) % _openEdgeMeshCount;
                    placedPieces[i] = piece;
                    break;
                }
            }
        }

        if (IsKeyPressed(KEY_V))
        {
            for (int i= placedPieceCount - 1; i>=0;i--)
            {
                if (placedPieces[i].offsetX == cursor[0] && placedPieces[i].offsetY == cursor[1] && placedPieces[i].offsetZ == cursor[2])
                {
                    // print compatible pieces of neighbors
                    for (int j=0;j<placedPieceCount;j++)
                    {
                        if (j == i) continue;
                        PlacedPiece piece = placedPieces[j];
                        int dx = (piece.offsetX - cursor[0]);
                        int dy = (piece.offsetY - cursor[1]);
                        int dz = (piece.offsetZ - cursor[2]);
                        if (dx*dx+dy*dy+dz*dz > 1) continue;
                        int matchResult = OpenEdgeMesh_isMatch(&_openEdgeMeshes[placedPieces[i].pieceIndex], &_openEdgeMeshes[piece.pieceIndex], dx, dy, dz);
                        printf("Neighbor %d-%d (%d %d %d): %d\n", placedPieces[i].pieceIndex, piece.pieceIndex, dx, dy, dz, matchResult);
                    }
                    break;
                }
            }
        }

        if (IsKeyPressed(KEY_ENTER) || IsKeyPressed(KEY_BACKSPACE))
        {
            int direction = IsKeyPressed(KEY_BACKSPACE) ? -1 : 1;
            PlacedPieces_shift(placedPieces, placedPieceCount, -cursor[0], -cursor[1], -cursor[2]);
            int placedPieceIndex = -1;
            for (int i=placedPieceCount-1;i>=0 && !IsKeyDown(KEY_LEFT_SHIFT);i--)
            {
                if (placedPieces[i].offsetX == 0 && placedPieces[i].offsetY == 0 && placedPieces[i].offsetZ == 0)
                {
                    PlacedPiece piece = placedPieces[i];
                    placedPieceIndex = placedPieceCount;
                    placedPieces[i] = placedPieces[--placedPieceCount];
                    int currentIndex = placedPieces[i].pieceIndex;
                    int count = PiecePatcher_getCompatibleMeshCountEx(placedPieces, placedPieceCount, MESH_TILE_MATCH_BUILD_MASK);

                    for (int j=0;j<count;j++)
                    {
                        int index = PiecePatcher_getCompatibleMeshByIndexEx(placedPieces, placedPieceCount, MESH_TILE_MATCH_BUILD_MASK, j);

                        if (index == currentIndex)
                        {
                            placedPieces[i].pieceIndex = PiecePatcher_getCompatibleMeshByIndexEx(placedPieces, placedPieceCount, MESH_TILE_MATCH_BUILD_MASK, (j + direction + count) % count);
                            printf("Replacing piece %d: %d\n", i, placedPieces[i].pieceIndex);
                            placedPieceCount++;
                            break;
                        }
                    }

                    break;
                }
            }

            if (placedPieceIndex == -1)
            {
                int newIndex = PiecePatcher_getCompatibleMeshByIndexEx(placedPieces, placedPieceCount, MESH_TILE_MATCH_BUILD_MASK, 0);
                printf("Placing piece %d\n", newIndex);
                if (newIndex >= 0)
                {
                    placedPieces[placedPieceCount++] = (PlacedPiece){
                        .offsetX = 0,
                        .offsetY = 0,
                        .offsetZ = 0,
                        .pieceIndex = newIndex,
                    };
                }
            }
            
            PlacedPieces_shift(placedPieces, placedPieceCount, cursor[0], cursor[1], cursor[2]);
        }

        for (int i=0;i<placedPieceCount;i++)
        {
            PlacedPiece piece = placedPieces[i];
            OpenEdgeMesh *mesh = &_openEdgeMeshes[piece.pieceIndex];
            Matrix m = MatrixIdentity();
            m = MatrixMultiply(m, MatrixRotateY(PI * 0.5f * -mesh->rotation));
            m = MatrixMultiply(m, MatrixTranslate(piece.offsetX, piece.offsetY, piece.offsetZ));
            DrawMesh(*mesh->mesh, mat, m);
        }
    }
    break;
    }
    rlPopMatrix();

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