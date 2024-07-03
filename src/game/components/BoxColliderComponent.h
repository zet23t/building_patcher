#ifdef COMPONENT

COMPONENT(BoxColliderComponent)

#elif defined(SERIALIZABLE_STRUCT_START)

SERIALIZABLE_STRUCT_START(BoxColliderComponent)
    SERIALIZABLE_FIELD(Vector3, size)
    SERIALIZABLE_FIELD(Vector3, offset)
SERIALIZABLE_STRUCT_END(BoxColliderComponent)


#elif defined(COMPONENT_IMPLEMENTATION)

#include "BoxColliderComponent.c"
#else

#ifndef BOX_COLLIDER_COMPONENT_H
#define BOX_COLLIDER_COMPONENT_H
SceneObjectId BoxColliderComponent_GetEntityIdByPoint(SceneGraph *graph, Vector3 position, float radius);
#endif

#endif