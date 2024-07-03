#ifdef COMPONENT

COMPONENT(GameplaySystemComponent)

#elif defined(SERIALIZABLE_STRUCT_START)

SERIALIZABLE_STRUCT_START(GameplaySystemComponent)
    SERIALIZABLE_FIELD(SceneObjectId, playerObjectId)
    SERIALIZABLE_FIELD(Vector3, playerVelocity)
    SERIALIZABLE_FIELD(SceneComponentId, cameraComponentId)
SERIALIZABLE_STRUCT_END(GameplaySystemComponent)

#elif defined(COMPONENT_IMPLEMENTATION)

#include "GameplaySystemComponent.c"

#endif