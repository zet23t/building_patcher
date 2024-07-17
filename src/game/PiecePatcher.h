#ifndef __PIECEPATCHER_H__
#define __PIECEPATCHER_H__

#include <raylib.h>
#include <inttypes.h>

// The mesh is tile combination is not compatible at all
#define MESH_TILE_MATCH_MISMATCH 0
// The mesh is tile combination is compatible
#define MESH_TILE_MATCH_COMPATIBLE 1
// The mesh tiles are fully disjunct (either because being too far away or because sharing no common corners)
#define MESH_TILE_MATCH_DISJUNCT 2
// The mesh tiles are are not connected but share the same cube position (offsetXYZ=0)
#define MESH_TILE_MATCH_CONGRUENTLY_DISJUNCT 4
// The mesh tiles are partially disjunct because one corner is shared
#define MESH_TILE_MATCH_PARTIALLY_DISJUNCT 8
// Similar to MESH_TILE_MATCH_CONGRUENTLY_DISJUNCT but with one corner matching
#define MESH_TILE_MATCH_PARTIALLY_CONGRUENTLY_DISJUNCT 16

#define MESH_TILE_MATCH_BUILD_MASK (MESH_TILE_MATCH_COMPATIBLE)

typedef struct PlacedPiece {
    int8_t offsetX, offsetY, offsetZ;
    uint16_t pieceIndex;
} PlacedPiece;

Material *PiecePatcher_getMaterial(int index);
int PiecePatcher_getRotation(int index);
Mesh *PiecePatcher_getMesh(int index);

int PiecePatcher_clearCompatibleFlags();
int PiecePatcher_getCompatibleMeshCountEx(PlacedPiece *pieces, int pieceCount, int validFlagMask);
int PiecePatcher_getCompatibleMeshByIndexEx(PlacedPiece *pieces, int pieceCount, int validFlagMask, int compatibleIndex);
int PiecePatcher_getCompatibleMeshCount(int index, int offsetX, int offsetY, int offsetZ);
int PiecePatcher_getCompatibleMeshByIndex(int index, int offsetX, int offsetY, int offsetZ, int compatibleIndex);
void PiecePatcherInit();
void PiecePatcherDrawDebug();

#endif