#include "game/game.h"


// add a random flappy bird level chunk geometry
SceneObjectId addChunk(SceneGraph *graph, int chunkX, int chunkLength, SceneComponentId cameraComponentId)
{
    SceneObjectId chunkId = SceneGraph_createObject(graph, "Chunk");
    SceneGraph_setLocalPosition(graph, chunkId, (Vector3){chunkX, 0.0f, 0.0f});
    // add ground

    SceneObjectId groundId = SceneGraph_createObject(graph, "Ground");
    SceneGraph_addComponent(graph, groundId, _componentIdMap.PrimitiveRendererComponentId,
        &(PrimitiveRendererComponent) {
            .primitiveType = PRIMITIVE_TYPE_CUBE,
            .color = (Color){0, 128, 0, 255},
            .size = (Vector3){chunkLength, 1.0f, 1.0f},
        });
    SceneGraph_addComponent(graph, groundId, _componentIdMap.BoxColliderComponentId,
        &(BoxColliderComponent) {
            .size = (Vector3){chunkLength, 1.0f, 1.0f},
        });
    SceneGraph_setLocalPosition(graph, groundId, (Vector3){chunkLength * 0.5f, -2.5f, 0.0f});
    SceneGraph_setParent(graph, groundId, chunkId);

    // add buildings
    SceneObjectId buildingId = SceneGraph_createObject(graph, "Building");
    int buildingWidth = 1; //GetRandomValue(1, 3);
    int buildingOffset = GetRandomValue(0, chunkLength - buildingWidth);
    int buildingHeight = 1;// GetRandomValue(1, 4);

    // SceneGraph_addComponent(graph, buildingId, _componentIdMap.PrimitiveRendererComponentId,
    //     &(PrimitiveRendererComponent) {
    //         .primitiveType = PRIMITIVE_TYPE_CUBE,
    //         .color = (Color){128, 128, 128, 255},
    //         .size = (Vector3){buildingWidth, buildingHeight, 1.0f},
    //     });
    SceneGraph_addComponent(graph, buildingId, _componentIdMap.BoxColliderComponentId,
        &(BoxColliderComponent) {
            .size = (Vector3){buildingWidth, buildingHeight, 1.0f},
            .offset = (Vector3){0.0f, buildingHeight * 0.5f, 0.0f},
        });
    SceneGraph_setLocalPosition(graph, buildingId, (Vector3){buildingOffset, (buildingHeight * 0.5f) - 2.0f, 0.0f});
    SceneGraph_setParent(graph, buildingId, chunkId);

    Model buildingMesh = ResourceManager_loadModel(_resourceManager, "assets/building.obj");
    SceneGraph_addComponent(graph, buildingId, _componentIdMap.MeshRendererComponentId,
        &(MeshRendererComponent) {
            .material = &buildingMesh.materials[buildingMesh.materialCount - 1],
            .mesh = &buildingMesh.meshes[0],
        });

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