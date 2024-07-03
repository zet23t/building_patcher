#include "game/game.h"
#include "raymath.h"

SceneObjectId BoxColliderComponent_GetEntityIdByPoint(SceneGraph *graph, Vector3 position, float radius)
{
    BoxColliderComponent *boxCollider;
    SceneComponent *boxComponent;
    int index = 0;
    while((boxComponent = SceneGraph_getComponentByType(graph, (SceneObjectId){0}, _componentIdMap.BoxColliderComponentId, (void**)&boxCollider, index++)))
    {
        SceneObject *boxObject = SceneGraph_getObject(graph, boxComponent->objectId);
        Matrix inv = SceneObject_getToLocalMatrix(boxObject);
        Vector3 local = Vector3Transform(position, inv);
        
        bool collides = CheckCollisionBoxSphere((BoundingBox){
            .max = (Vector3){boxCollider->size.x/2 + boxCollider->offset.x, boxCollider->size.y/2 + boxCollider->offset.y, boxCollider->size.z/2 + boxCollider->offset.z},
            .min = (Vector3){-boxCollider->size.x/2 + boxCollider->offset.x, -boxCollider->size.y/2 + boxCollider->offset.y, -boxCollider->size.z/2 + boxCollider->offset.z},
        }, (Vector3){local.x, local.y, local.z}, radius);
        if (collides)
        {
            return boxComponent->objectId;
        }
    }

    return (SceneObjectId){0};
}

BEGIN_COMPONENT_REGISTRATION(BoxColliderComponent) {
    0
} END_COMPONENT_REGISTRATION