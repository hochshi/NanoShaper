
//---------------------------------------------------------
/*    @file		tools.h
 *     @brief	tools.h Includes some common tools
 *							    							*/
//---------------------------------------------------------

#ifndef tools_h
#define tools_h

#include <globals.h>
#include <jama_eig.h>
#include <logging.h>
#include <sturm/solve.h>
#include <string>

namespace nanoshaper {
#ifdef DBGMEM_CRT
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW
#endif

using namespace TNT;
using namespace JAMA;

#define SMALL_SHIFT_MAP 7
#define SHIFT_MAP 27

/////////////////////////// MAPS ///////////////////////////////////

static const int shift_map[SHIFT_MAP][3] = {
    {0, 0, 0},    {0, 0, -1},  {0, 0, +1},   {0, +1, 0},   {0, +1, -1},
    {0, +1, +1},  {0, -1, 0},  {0, -1, -1},  {0, -1, +1},

    {+1, 0, 0},   {+1, 0, -1}, {+1, 0, +1},  {+1, +1, 0},  {+1, +1, -1},
    {+1, +1, +1}, {+1, -1, 0}, {+1, -1, -1}, {+1, -1, +1},

    {-1, 0, 0},   {-1, 0, -1}, {-1, 0, +1},  {-1, +1, 0},  {-1, +1, -1},
    {-1, +1, +1}, {-1, -1, 0}, {-1, -1, -1}, {-1, -1, +1}};

static const uint32_t MASK[32] = {
    0x00000001, 0x00000002, 0x00000004, 0x00000008, 0x00000010, 0x00000020,
    0x00000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800,
    0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000,
    0x00040000, 0x00080000, 0x00100000, 0x00200000, 0x00400000, 0x00800000,
    0x01000000, 0x02000000, 0x04000000, 0x08000000, 0x10000000, 0x20000000,
    0x40000000, 0x80000000,
};
static const int shift_map2[SMALL_SHIFT_MAP][3] = {
    {0, 0, 0},  {+1, 0, 0}, {0, +1, 0}, {0, 0, +1},
    {-1, 0, 0}, {0, -1, 0}, {0, 0, -1},
};

/// Marching cubes tables
static const int edgeTable[256] = {
    0x0,   0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f,
    0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190, 0x99,  0x393, 0x29a, 0x596, 0x49f,
    0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230,
    0x339, 0x33,  0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936,
    0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa,  0x7a6, 0x6af, 0x5a5,
    0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569,
    0x663, 0x76a, 0x66,  0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a,
    0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff,  0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453,
    0x55a, 0x256, 0x35f, 0x55,  0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53,
    0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,  0xfcc,
    0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9, 0xac3, 0xbca,
    0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc,  0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9,
    0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55,
    0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6,
    0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0xff,  0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f,
    0x66,  0x76a, 0x663, 0x569, 0x460, 0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af,
    0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa,  0x1a3, 0x2a9, 0x3a0, 0xd30,
    0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636,
    0x13a, 0x33,  0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895,
    0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99,  0x190, 0xf00, 0xe09,
    0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a,
    0x203, 0x109, 0x0};

// xg-hside = 0
// xg+hside = 1
// yg-hside = 2
// yg+hside = 3
// zg-hside = 4
// zg+hside = 5

static const int mulTable[8][3] = {{0, 2, 4}, {1, 2, 4}, {1, 3, 4}, {0, 2, 4},
                                   {0, 2, 5}, {1, 2, 5}, {1, 3, 5}, {0, 3, 5}};

static const int triTable[256][16] = {
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

///////////////////////////////////////////////////////////////////////////////////

#define MAX_POLY_DEGREE \
  10  // max polynomial degree supported by the bgp projection routine (6 is \
      // needed for skin, 4 by Connolly)

/////////////////////////////////SIMPLE CLASSES
//////////////////////////////////////////

/**@brief A point in R^3*/
class Point {
 public:
  double pos[3];
  Point() {
    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;
  }
  Point(const double x, const double y, const double z) {
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
  }
  Point(const Point& p) {
    pos[0] = p.pos[0];
    pos[1] = p.pos[1];
    pos[2] = p.pos[2];
  }

  double squaredDist(const Point* const x) {
    double acc = 0;
    for (int i = 0; i < 3; i++) {
      double t = x->pos[i] - pos[i];
      acc += t * t;
    }
    return acc;
  }

  double dot(const Point* const x) {
    double acc = 0;
    for (int i = 0; i < 3; i++)
      acc += x->pos[i] * pos[i];
    return acc;
  }

  double norm(const Point* const x) {
    double acc = 0;
    for (int i = 0; i < 3; i++)
      acc += (pos[i]) * (pos[i]);
    return sqrt(acc);
  }

  double dist(const Point* const x) {
    double acc = 0;
    for (int i = 0; i < 3; i++) {
      double t = x->pos[i] - pos[i];
      acc += t * t;
    }
    return sqrt(acc);
  }

  void sumInPlace(const Point* const x) {
    for (int i = 0; i < 3; i++) {
      pos[0] += x->pos[0];
      pos[1] += x->pos[1];
      pos[2] += x->pos[2];
    }
  }

  Point* sum(const Point* const x) {
    Point* p = new Point();
    for (int i = 0; i < 3; i++) {
      p->pos[0] = pos[0] + x->pos[0];
      p->pos[1] = pos[1] + x->pos[1];
      p->pos[2] = pos[2] + x->pos[2];
    }
    return p;
  }
};

class AtomInfo {
 protected:
  string name;
  int resNum;
  string resName;
  string chain;

 public:
  AtomInfo() {
    name = "X";
    resNum = -1;
    resName = "UNK";
    chain = "X";
  }

  AtomInfo(const string name, const int resNum, const string resName,
           const string chain) {
    this->name = name;
    this->resNum = resNum;
    this->resName = resName;
    this->chain = chain;
  }

  const string& getName() { return name; }

  void setName(const string name) { this->name = name; }

  const string& getResName() { return resName; }

  void setResName(const string res) { resName = res; }

  const int& getResNum() { return resNum; }

  void setResNum(const int res) { resNum = res; }

  void setChain(const string ch) { chain = ch; }

  const string& getChain() { return chain; }

  /** copy constructor*/
  AtomInfo(const AtomInfo& p) {
    name = p.name;
    resNum = p.resNum;
    resName = p.resName;
    chain = p.chain;
  }

  void operator=(const AtomInfo& p) {
    name = p.name;
    resNum = p.resNum;
    resName = p.resName;
    chain = p.chain;
  }
};

/**@brief Timer: timer class. If chrono is defined use high accuracy timer*/
class Timer {
 private:
#ifdef ENABLE_BOOST_CHRONO
  boost::chrono::high_resolution_clock::time_point start_t, end_t;
#else
  time_t start_t, end_t;
#endif

 public:
  Timer() {}

  void start() {
#ifdef ENABLE_BOOST_CHRONO
    start_t = boost::chrono::high_resolution_clock::now();
#else
    time(&start_t);
#endif
  }

  double stop() {
#ifdef ENABLE_BOOST_CHRONO
    end_t = boost::chrono::high_resolution_clock::now();
    boost::chrono::milliseconds ms =
        boost::chrono::duration_cast<boost::chrono::milliseconds>(end_t -
                                                                  start_t);
    return ms.count() / 1000.;
#else
    time(&end_t);
    return difftime(end_t, start_t);
#endif
  }
};

/**@brief An atom: an atom is a point in R^3 with radius, charge and dielectric
 * value*/
class Atom : public Point {
 public:
  double radius;
  double radius2;
  double charge;
  int dielectric;
  AtomInfo ai;

  Atom() : Point() {
    radius = 0;
    charge = 0;
    dielectric = 0;
  }

  Atom(const double x, const double y, const double z, const double r)
      : Point(x, y, z) {
    radius = r;
    radius2 = r * r;
    charge = 0;
    dielectric = 0;
  }

  Atom(const double x, const double y, const double z, const double r,
       const double c, const int d)
      : Point(x, y, z) {
    radius = r;
    radius2 = r * r;
    charge = c;
    dielectric = d;
  }

  Atom(const double x, const double y, const double z, const double r,
       const double c, const int d, const string name, const string resName,
       const int resNum, const string chain)
      : Point(x, y, z) {
    radius = r;
    radius2 = r * r;
    charge = c;
    dielectric = d;
    ai.setName(name);
    ai.setResName(resName);
    ai.setChain(chain);
    ai.setResNum(resNum);
  }

  // copy constructor
  Atom(Atom& p) : Point(p) {
    radius = p.radius;
    charge = p.charge;
    dielectric = p.dielectric;
    ai.setChain(p.ai.getChain());
    ai.setResNum(p.ai.getResNum());
    ai.setResName(p.ai.getResName());
    ai.setName(p.ai.getName());
  }

  void print() {
    logging::log<logging::level::info>("Atom Info:");
    logging::log<logging::level::info>("\tname -> {}", ai.getName());
    logging::log<logging::level::info>("\tresidue name -> {}", ai.getResName());
    logging::log<logging::level::info>("\tresidue number -> {}",
                                       ai.getResNum());
    logging::log<logging::level::info>("\tchain -> {}", ai.getChain());
    logging::log<logging::level::info>("\tposx -> {}", pos[0]);
    logging::log<logging::level::info>("\tposy -> {}", pos[1]);
    logging::log<logging::level::info>("\tposz -> {}", pos[2]);
    logging::log<logging::level::info>("\tcharge -> {}", charge);
    logging::log<logging::level::info>("\tradius -> {}", radius);
    logging::log<logging::level::info>("\tdielectric -> {}", dielectric);
  }
};

////////////////// Memory management routines //////////////////////

template <class T>
void cleanDelete(T*& t) {
  if (t != NULL) {
    delete t;
    t = NULL;
  } else {
    logging::log<logging::level::warn>(
        "Attempting to de-allocate a null object!");
  }
}

template <class T>
void cleanDeleteV(T*& t) {
  if (t != NULL) {
    delete[] t;
    t = NULL;
  } else {
    logging::log<logging::level::warn>(
        "Attempting to de-allocate a null vector object!");
  }
}

template <class T>
inline T read4DVector(const T* const var, const int i, const int j, const int k,
                      const int l, const int nx, const int ny, const int nz,
                      const int nl) {
#ifdef CHECK_BOUNDS
  if (i >= nx || j >= ny || k >= nz || l >= nl || i < 0 || j < 0 || k < 0 ||
      l < 0) {
    logging::log<logging::level::err>("Out of bound error in reading 4DVector");
    throw std::out_of_range("Out of bound error in reading 4DVector");
  }
#endif

  size_t index = i + nx * (j + ny * (k + nz * l));
  return var[index];
}

template <class T>
inline void write4DVector(T* const var, const T val, const int i, const int j,
                          const int k, const int l, const int nx, const int ny,
                          const int nz, const int nl) {
#ifdef CHECK_BOUNDS
  if (i >= nx || j >= ny || k >= nz || l >= nl || i < 0 || j < 0 || k < 0 ||
      l < 0) {
    logging::log<logging::level::err>("Out of bound error in writing 4DVector");
    throw std::out_of_range("Out of bound error in writing 4DVector");
  }
#endif

  size_t index = i + nx * (j + ny * (k + nz * l));
  var[index] = val;
}

template <class T>
inline T read3DVector(const T* const var, const int i, const int j, const int k,
                      const int nx, const int ny, const int nz) {
#ifdef CHECK_BOUNDS
  if (i >= nx || j >= ny || k >= nz || i < 0 || j < 0 || k < 0) {
    logging::log<logging::level::err>("Out of bound error in reading 3DVector");
    throw std::out_of_range("Out of bound error in reading 3DVector");
  }
#endif

  size_t index = i + nx * (j + ny * k);
  return var[index];
}

template <class T>
inline void write3DVector(T* const var, const T val, const int i, const int j,
                          const int k, const int nx, const int ny,
                          const int nz) {
#ifdef CHECK_BOUNDS
  if (i >= nx || j >= ny || k >= nz || i < 0 || j < 0 || k < 0) {
    logging::log<logging::level::err>("Out of bound error in writing 3DVector");
    throw std::out_of_range("Out of bound error in writing 3DVector");
  }
#endif

  size_t index = i + nx * (j + ny * k);
  var[index] = val;
}

template <class T>
inline T read2DVector(const T* const var, const int i, const int j,
                      const int nx, const int ny) {
#ifdef CHECK_BOUNDS
  if (i >= nx || j >= ny || i < 0 || j < 0) {
    logging::log<logging::level::err>("Out of bound error in reading 2DVector");
    throw std::out_of_range("Out of bound error in reading 2DVector");
  }
#endif

  size_t index = i + nx * j;
  return var[index];
}

template <class T>
inline void write2DVector(T* const var, const T val, const int i, const int j,
                          const int nx, const int ny) {
#ifdef CHECK_BOUNDS
  if (i >= nx || j >= ny || i < 0 || j < 0) {
    logging::log<logging::level::err>("Out of bound error in reading 2DVector");
    throw std::out_of_range("Out of bound error in reading 2DVector");
  }
#endif

  size_t index = i + nx * j;
  var[index] = val;
}

template <class T>
inline T readVector(const T* const var, const int i, const int nx) {
#ifdef CHECK_BOUNDS
  if (i >= nx || i < 0) {
    logging::log<logging::level::err>("Out of bound error in reading Vector");
    throw std::out_of_range("Out of bound error in reading Vector");
  }
#endif

  return var[i];
}

template <class T>
inline void writeVector(T* const var, const T val, const int i, const int nx) {
#ifdef CHECK_BOUNDS
  if (i >= nx || i < 0) {
    logging::log<logging::level::err>("Out of bound error in writing Vector");
    throw std::out_of_range("Out of bound error in writing Vector");
  }
#endif
  var[index] = val;
}

template <class T>
T* allocateVector(size_t n) {
  T* t = (T*)malloc(sizeof(T) * n);
  if (t == NULL) {
    logging::log<logging::level::err>("Not enough memory to allocate vector ");
    return NULL;
  }
  return t;
}

template <class T>
T** allocateMatrix2D(int nrows, int ncol) {
  T** t = (T**)malloc(sizeof(T*) * nrows);
  if (t == NULL) {
    logging::log<logging::level::err>(
        "Not enough memory to allocate 2D matrix ");
    return NULL;
  }
  for (int i = 0; i < nrows; i++) {
    t[i] = (T*)malloc(sizeof(T) * ncol);

    if (t[i] == NULL) {
      logging::log<logging::level::err>(
          "Not enough memory to allocate 2D matrix ");
      return NULL;
    }
  }
  return t;
}

template <class T>
T*** allocateMatrix3D(int nx, int ny, int nz) {
  T*** t = (T***)malloc(sizeof(T**) * nx);
  if (t == NULL) {
    logging::log<logging::level::err>(
        "Not enough memory to allocate 3D matrix ");
    return NULL;
  }
  for (int i = 0; i < nx; i++) {
    t[i] = (T**)malloc(sizeof(T*) * ny);
    if (t[i] == NULL) {
      logging::log<logging::level::err>(
          "Not enough memory to allocate 3D matrix ");
      return NULL;
    }

    for (int j = 0; j < ny; j++) {
      t[i][j] = (T*)malloc(sizeof(T) * nz);
      if (t[i][j] == NULL) {
        logging::log<logging::level::err>(
            "Not enough memory to allocate 3D matrix ");
        return NULL;
      }
    }
  }
  return t;
}

template <class T>
void deleteVector(T*& t) {
  if (t != NULL) {
    free(t);
    t = NULL;
  } else {
    logging::log<logging::level::warn>(
        "Attempting to de-allocate a null vector!");
  }
}

template <class T>
void deleteMatrix2D(int nrows, T**& t) {
  if (t != NULL) {
    for (int i = 0; i < nrows; i++)
      free(t[i]);
    free(t);
    t = NULL;
  } else {
    logging::log<logging::level::warn>(
        "Attempting to de-allocate a null vector!");
  }
}

template <class T>
void deleteMatrix3D(int nx, int ny, T***& t) {
  if (t != NULL) {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++)
        free(t[i][j]);
      free(t[i]);
    }
    free(t);
    t = NULL;
  } else {
    logging::log<logging::level::warn>(
        "Attempting to de-allocate a null vector!");
  }
}

///////////////////////////////////////////////////////////////////

/** a vector bool is a vector in which it is assumed that the information is
written/read in a bitwise way to minimize memory footprint. Thus the allocator
computes how many 32 bits words are needed to store n elements*/
inline uint32_t* allocateVectorBool(size_t n) {
  int nw = n / 32;
  if ((n % 32) != 0)
    nw++;
  uint32_t* ptr = (uint32_t*)malloc(nw * 4);
  if (ptr == NULL)
    logging::log<logging::level::err>(
        "Not enough memory to allocate bit vector");
  return ptr;
}

inline void deleteVectorBool(uint32_t*& t) {
  if (t != NULL) {
    free(t);
    t = NULL;
  } else {
    logging::log<logging::level::warn>(
        "Attempting to de-allocate a null bool vector!");
  }
}

inline bool read3DVectorBool(uint32_t* const var, const int i, const int j,
                             const int k, const int nx, const int ny,
                             const int nz) {
#ifdef CHECK_BOUNDS
  if (i >= nx || j >= ny || k >= nz || i < 0 || j < 0 || k < 0) {
    logging::log<logging::level::err>("Out of bound error in reading 3DVector");
    throw std::out_of_range("Out of bound error in reading 3DVector");
  }
#endif

  size_t index = i + nx * (j + ny * k);
  size_t word_index = index / 32;
  size_t position = index % 32;
  uint32_t mask = MASK[position];
  return ((var[word_index] & mask) == mask);
}

inline void write3DVectorBool(uint32_t* const var, const bool val, const int i,
                              const int j, const int k, const int nx,
                              const int ny, const int nz) {
#ifdef CHECK_BOUNDS
  if (i >= nx || j >= ny || k >= nz || i < 0 || j < 0 || k < 0) {
    logging::log<logging::level::err>("Out of bound error in reading 3DVector");
    throw std::out_of_range("Out of bound error in reading 3DVector");
  }
#endif

  size_t index = i + nx * (j + ny * k);
  size_t word_index = index / 32;
  size_t position = index % 32;
  uint32_t mask = MASK[position];
  var[word_index] = val ? (var[word_index] | mask) : (var[word_index] & ~mask);
}

/** portable aligned memory allocator.You will lose a bit of memory. For big
arrays that are heavily accessed this is the best allocation choice. For small
arrays very often allocated/deallocated this can introduce a signifcant memory
overhead. Free only using deleteAlignedVector. Typical bit alignement values are
32,64,128,256*/
template <class T>
T* allocateAlignedVector(size_t n, int bitAlignment) {
  T* mem = (T*)malloc(n * sizeof(T) + bitAlignment + sizeof(void*));
  T** ptr =
      (T**)((long)(mem + bitAlignment + sizeof(void*)) & ~(bitAlignment - 1));
  ptr[-1] = mem;

  if (mem == NULL) {
    logging::log<logging::level::err>(
        "Not enough memory to allocate aligned vector ");
    return NULL;
  }

  return (T*)ptr;
}

template <class T>
void deleteAlignedVector(T*& ptr) {
  if (ptr == NULL) {
    logging::log<logging::level::warn>(
        "Attempting to de-allocate a null aligned vector!");
    return;
  }

  free(((T**)ptr)[-1]);
  ptr = NULL;
}

///////////////////// STL comparators //////////////////////////////

/**@brief ascending on first double of pair<double,double*> comparator*/
bool compKeepIndex(pair<double, double*> a, pair<double, double*> b);

///////////////////////////////////////////////////////////////////

typedef std::pair<double, int> indexed_double;
bool index_double_comparator(const indexed_double& l, const indexed_double& r);

///////////////////// MISCELLANEAOUS ////////////////////////////////

/** @brief test if point proj is in cube whose center is point and side is side
 variable. If a toll is provided, the test is performed up to a tollerance
 value. By default the tollerance is 0*/
bool testInCube(const double proj0, const double proj1, const double proj2,
                const double point0, const double point1, const double point2,
                const double side, const double toll = 0);
double rintp(const double x);
string toLowerCase(string str);
void cleanLine();

static long int SEED = 1234;
// randnumber between 0 and 1
double randnum();

/** get the real roots by computing the companion matrix and then extracting the
eigenvalues. A root is real if its imaginary part in absolute value is less than
a given threshold. Usually this threshold is rather conservative such that
possibly unprecise real roots are not lost. Poly is modified in place to remove
possible zeros*/
void getRealRootsCompanion(double* const poly, const int degree,
                           double* const roots, int& numroots);
/** get real roots by using Sturm method. Directly search real roots. Much
faster, often less accurate than companion matrix*/
void getRealRootsSturm(const double* const polyy, const int degree,
                       double* const roots, int& numrootss);
/** plane by 3 points routine*/
void plane3points(const double p1[3], const double p2[3], const double p3[3],
                  double w[4], const bool normalize = true);
/** point to plane projection*/
void point2plane(const double p[3], double w[4], double* const dist,
                 double proj[3]);
/** in place inversion of 4x4 matrix*/
void inplace_invert4x4(double M[4][4]);
/** unrolled 4x4 matrix multiply*/
void Matrix4x4MultiplyBy4x4(const double src1[4][4], const double src2[4][4],
                            double dest[4][4]);
/** ray/sphere intersection. ray is o+t*dir */
bool raySphere(const double* const orig, const double* const dir,
               const double* const sphere_center, const double sphere_radius,
               double* const t1, double* const t2);
/** get the normal to a sphere*/
void getNormalToSphere(const double* const y, const double* const center,
                       const double radius, double* const normal);
/** project y to sphere and returns the normal*/
void projectToSphere(const double* const y, const double* const center,
                     const double radius, double* const proj, double& dist);

}  // namespace nanoshaper
#endif
