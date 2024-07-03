#include "game/game.h"

#ifndef COMPONENT_IMPLEMENTATION
#include "BoxColliderComponent.h"
#endif

#define GRAVITY 5.0f

void GameplaySystemComponent_update(SceneObject* sceneObject, SceneComponentId sceneComponent,
        float delta, void* componentData)
{
    GameplaySystemComponent *gameplay = (GameplaySystemComponent*)componentData;
    SceneGraph *graph = sceneObject->graph;
    CameraComponent *camera;
    SceneComponent *cam = SceneGraph_getComponent(graph, gameplay->cameraComponentId, (void**)&camera);
    if (cam == NULL) {
        return;
    }
    
    gameplay->playerVelocity.y -= delta * GRAVITY;
    if (gameplay->playerVelocity.y < -2.0f) {
        gameplay->playerVelocity.y = -2.0f;
    }
    if (IsKeyPressed(KEY_SPACE)) {
        gameplay->playerVelocity.y += 4.0f;
        if (gameplay->playerVelocity.y > 3.0f) {
            gameplay->playerVelocity.y = 3.0f;
        }
    }
    Vector3 playerPos = SceneGraph_getLocalPosition(graph, gameplay->playerObjectId);
    playerPos.x += delta;
    playerPos.y += gameplay->playerVelocity.y * delta;

    SceneObjectId colliding = BoxColliderComponent_GetEntityIdByPoint(graph, playerPos, 0.5f);
    if (colliding.version != 0) {
        playerPos.y = 0;
        gameplay->playerVelocity.y = 0;
        playerPos.x = playerPos.x - 1.0f;
    }
    SceneGraph_setLocalPosition(graph, gameplay->playerObjectId, playerPos);

    camera->targetPoint = (Vector3){playerPos.x + 3.0f, 0.0f, 0.0f};


}

#undef GRAVITY

BEGIN_COMPONENT_REGISTRATION(GameplaySystemComponent) {
    .updateTick = GameplaySystemComponent_update,
} END_COMPONENT_REGISTRATION