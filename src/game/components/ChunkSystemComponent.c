#include "game/game.h"
#include "game/PiecePatcher.h"


// add a random flappy bird level chunk geometry
SceneObjectId addChunk(SceneGraph *graph, int chunkX, int chunkLength, SceneComponentId cameraComponentId)
{
    SceneObjectId chunkId = SceneGraph_createObject(graph, "Chunk");
    SceneGraph_setLocalPosition(graph, chunkId, (Vector3){chunkX, 0.0f, 0.0f});
    // add ground

    // SceneObjectId groundId = SceneGraph_createObject(graph, "Ground");
    // SceneGraph_addComponent(graph, groundId, _componentIdMap.PrimitiveRendererComponentId,
    //     &(PrimitiveRendererComponent) {
    //         .primitiveType = PRIMITIVE_TYPE_CUBE,
    //         .color = (Color){0, 128, 0, 255},
    //         .size = (Vector3){chunkLength, 1.0f, 1.0f},
    //     });
    // SceneGraph_addComponent(graph, groundId, _componentIdMap.BoxColliderComponentId,
    //     &(BoxColliderComponent) {
    //         .size = (Vector3){chunkLength, 1.0f, 1.0f},
    //     });
    // SceneGraph_setLocalPosition(graph, groundId, (Vector3){chunkLength * 0.5f, -2.5f, 0.0f});
    // SceneGraph_setParent(graph, groundId, chunkId);

    // // add buildings
    SceneObjectId buildingId = SceneGraph_createObject(graph, "Building");
    // int buildingWidth = 1; //GetRandomValue(1, 3);
    // int buildingOffset = GetRandomValue(0, chunkLength - buildingWidth);
    // int buildingHeight = 1;// GetRandomValue(1, 4);

    // // SceneGraph_addComponent(graph, buildingId, _componentIdMap.PrimitiveRendererComponentId,
    // //     &(PrimitiveRendererComponent) {
    // //         .primitiveType = PRIMITIVE_TYPE_CUBE,
    // //         .color = (Color){128, 128, 128, 255},
    // //         .size = (Vector3){buildingWidth, buildingHeight, 1.0f},
    // //     });
    // SceneGraph_addComponent(graph, buildingId, _componentIdMap.BoxColliderComponentId,
    //     &(BoxColliderComponent) {
    //         .size = (Vector3){buildingWidth, buildingHeight, 1.0f},
    //         .offset = (Vector3){0.0f, buildingHeight * 0.5f, 0.0f},
    //     });
    // SceneGraph_setLocalPosition(graph, buildingId, (Vector3){buildingOffset, (buildingHeight * 0.5f) - 2.0f, 0.0f});
    SceneGraph_setParent(graph, buildingId, chunkId);

    PlacedPiece placedPieces[4*4*4*4];
    int placedPieceCount = 0;
    
    #define ADD_PIECE(x, y, z, index) { \
        placedPieces[placedPieceCount++] = (PlacedPiece){ \
            .offsetX = x, \
            .offsetY = y, \
            .offsetZ = z, \
            .pieceIndex = index, \
        }; \
    }
    #define SHIFT(dx, dy, dz) for (int _i = 0; _i < placedPieceCount; _i++) { \
        placedPieces[_i].offsetX += dx; \
        placedPieces[_i].offsetY += dy; \
        placedPieces[_i].offsetZ += dz; \
    }
    #define ADD_RANDOM_PIECE() \
        { int _count = PiecePatcher_getCompatibleMeshCountEx(placedPieces, placedPieceCount, MESH_TILE_MATCH_BUILD_MASK); \
        if (_count > 0) { \
            int _select = GetRandomValue(0, _count - 1); \
            int _index = PiecePatcher_getCompatibleMeshByIndexEx(placedPieces, placedPieceCount, MESH_TILE_MATCH_BUILD_MASK, _select); \
            ADD_PIECE(0, 0, 0, _index); \
        } \
        else TraceLog(LOG_WARNING, "No compatible pieces found %s:%d",__FILE__, __LINE__);\
        }
    #define ADD_RANDOM_PIECE_AT(x,y,z) \
        SHIFT(x, y, z); ADD_RANDOM_PIECE(); SHIFT(-x, -y, -z);

    ADD_PIECE(0, 0, 0, 1);
    ADD_RANDOM_PIECE_AT(1, 0, 0);
    ADD_RANDOM_PIECE_AT(1, 0, 0);
    ADD_RANDOM_PIECE_AT(0, 1, 0);
    ADD_RANDOM_PIECE_AT(1, 1, 0);
    

    #undef ADD_PIECE
    #undef SHIFT
    #undef ADD_RANDOM_PIECE
    #undef ADD_RANDOM_PIECE_AT

    for (int i=0; i < placedPieceCount; i++)
    {
        int index = placedPieces[i].pieceIndex;
        char name[64];
        Mesh *mesh = PiecePatcher_getMesh(index);
        snprintf(name, sizeof(name), "C:%s-%d", mesh->name, index);

        SceneObjectId element = SceneGraph_createObject(graph, name);
        SceneGraph_addComponent(graph, element, _componentIdMap.MeshRendererComponentId,
            &(MeshRendererComponent) {
                .material = PiecePatcher_getMaterial(index),
                .mesh = mesh,
            });
        SceneGraph_setParent(graph, element, buildingId);
        SceneGraph_setLocalPosition(graph, element, (Vector3){placedPieces[i].offsetX, placedPieces[i].offsetY, placedPieces[i].offsetZ});
        SceneGraph_setLocalRotation(graph, element, (Vector3){0.0f, -90.0f * PiecePatcher_getRotation(index), 0.0f});
    }



    // Model buildingMesh = ResourceManager_loadModel(_resourceManager, "assets/building.obj");
    // int index = 1;

    // SceneGraph_addComponent(graph, buildingId, _componentIdMap.MeshRendererComponentId,
    //     &(MeshRendererComponent) {
    //         .material = PiecePatcher_getMaterial(index),
    //         .mesh = PiecePatcher_getMesh(index),
    //     });
    // SceneGraph_setLocalRotation(graph, buildingId, (Vector3){0.0f, -90.0f * PiecePatcher_getRotation(index), 0.0f});
    // for (int i=0;i < 3;i++)
    // {
    //     int compatibleCount = PiecePatcher_getCompatibleMeshCount(index, 1, 0, 0);
    //     if (compatibleCount == 0) {
            
    //         return chunkId;
    //     }
    //     int indices[2];
    //     for (int z=0;z<2;z++)
    //     {
    //         int select = GetRandomValue(0, compatibleCount - 1);
    //         int compatibleIndex = PiecePatcher_getCompatibleMeshByIndex(index, 0, 0, z + 1, select);
    //         char name[64];
    //         Mesh *mesh = PiecePatcher_getMesh(compatibleIndex);
    //         snprintf(name, sizeof(name), "C:%s-%d", mesh->name, compatibleIndex);

    //         SceneObjectId element = SceneGraph_createObject(graph, name);
    //         SceneGraph_addComponent(graph, element, _componentIdMap.MeshRendererComponentId,
    //             &(MeshRendererComponent) {
    //                 .material = PiecePatcher_getMaterial(compatibleIndex),
    //                 .mesh = mesh,
    //             });
    //         SceneGraph_setParent(graph, element, chunkId);
    //         SceneGraph_setLocalPosition(graph, element, (Vector3){0.0f, 0.0f, 1.0f * (i + 1)});
    //         SceneGraph_setLocalRotation(graph, element, (Vector3){0.0f, -90.0f * PiecePatcher_getRotation(compatibleIndex), 0.0f});
    //         index = indices[z] = compatibleIndex;
    //     }


    //     index = indices[0];
    // }

    return chunkId;
}


void ChunkSystemComponent_update(SceneObject* sceneObject, SceneComponentId sceneComponent,
        float delta, void* componentData)
{
    SceneGraph *graph = sceneObject->graph;
    ChunkSystemComponent *chunk = (ChunkSystemComponent*)componentData;
    SceneComponent *cam = SceneGraph_getComponent(graph, chunk->cameraComponentId, NULL);
    if (cam == NULL) {
        return;
    }
    Vector3 camPos = SceneGraph_getWorldPosition(graph, cam->objectId);
    
    int minX = (int)(camPos.x / chunk->chunkLength) - 3;
    int maxX = (int)(camPos.x / chunk->chunkLength) + 3;

    for (int i=0; i < CHUNK_SYSTEM_COMPONENT_CHUNKCOUNT; i++)
    {
        if (chunk->chunkIds[i].version == 0) {
            continue;
        }

        Vector3 chunkPos = SceneGraph_getLocalPosition(graph, chunk->chunkIds[i]);
        if (chunkPos.x < minX * chunk->chunkLength || chunkPos.x > maxX * chunk->chunkLength)
        {
            SceneGraph_destroyObject(graph, chunk->chunkIds[i]);
            chunk->chunkIds[i] = (SceneObjectId){0};
        }
    }

    for (int x=minX; x <= maxX; x++)
    {
        int found = 0;
        for (int j=0; j < CHUNK_SYSTEM_COMPONENT_CHUNKCOUNT; j++)
        {
            if (chunk->chunkIds[j].version == 0)
            {
                continue;
            }
            Vector3 chunkPos = SceneGraph_getLocalPosition(graph, chunk->chunkIds[j]);
            if (chunkPos.x == x * chunk->chunkLength)
            {
                found = 1;
                break;
            }
        }

        if (found) continue;

        for (int j=0; j < CHUNK_SYSTEM_COMPONENT_CHUNKCOUNT; j++)
        {
            if (chunk->chunkIds[j].version == 0)
            {
                chunk->chunkIds[j] = addChunk(graph, x * chunk->chunkLength, chunk->chunkLength, chunk->cameraComponentId);
                break;
            }
        }
    }
}

BEGIN_COMPONENT_REGISTRATION(ChunkSystemComponent) {
    .updateTick = ChunkSystemComponent_update,
} END_COMPONENT_REGISTRATION