#ifndef CHUNK_SYSTEM_COMPONENT_CHUNKCOUNT
#define CHUNK_SYSTEM_COMPONENT_CHUNKCOUNT 16
#endif

#ifdef COMPONENT
COMPONENT(ChunkSystemComponent)

#elif defined(SERIALIZABLE_STRUCT_START)

SERIALIZABLE_STRUCT_START(ChunkSystemComponent)
    SERIALIZABLE_FIELD(int, chunkLength)
    SERIALIZABLE_FIELD(SceneComponentId, cameraComponentId)
    SERIALIZABLE_FIXED_ARRAY(SceneObjectId, chunkIds, CHUNK_SYSTEM_COMPONENT_CHUNKCOUNT)
SERIALIZABLE_STRUCT_END(ChunkSystemComponent)


#elif defined(COMPONENT_IMPLEMENTATION)

#include "ChunkSystemComponent.c"

#endif