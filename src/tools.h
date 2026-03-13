
//---------------------------------------------------------
/*    @file		tools.h
*     @brief	tools.h Includes some common tools
*							    							*/
//---------------------------------------------------------

#ifndef tools_h
#define tools_h

#include "globals.h"

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif

#include "./sturm/solve.h"
#include <string>

#if defined(ENABLE_CGAL) && defined(CGAL_LINKED_WITH_TBB)
#include <tbb/global_control.h>
#endif

#include "jama_eig.h"

using namespace TNT;
using namespace JAMA;

#define SMALL_SHIFT_MAP 7
#define SHIFT_MAP 27

/////////////////////////// MAPS ///////////////////////////////////

static const int shift_map[SHIFT_MAP][3] =	{
						{0,0,0},
						{0,0,-1},
						{0,0,+1},
						{0,+1,0},
						{0,+1,-1},
						{0,+1,+1},
						{0,-1,0},
						{0,-1,-1},
						{0,-1,+1},

						{+1,0,0},
						{+1,0,-1},
						{+1,0,+1},
						{+1,+1,0},
						{+1,+1,-1},
						{+1,+1,+1},
						{+1,-1,0},
						{+1,-1,-1},
						{+1,-1,+1},

						{-1,0,0},
						{-1,0,-1},
						{-1,0,+1},
						{-1,+1,0},
						{-1,+1,-1},
						{-1,+1,+1},
						{-1,-1,0},
						{-1,-1,-1},
						{-1,-1,+1}	};

#ifdef ENABLE_BOOST_THREADS
static const boost::uint32_t MASK[32] = {
	0x00000001,
	0x00000002,
	0x00000004,
	0x00000008, 
	0x00000010, 
	0x00000020, 
	0x00000040, 
	0x00000080, 
	0x00000100, 
	0x00000200, 
	0x00000400, 
	0x00000800, 
	0x00001000, 
	0x00002000, 
	0x00004000, 
	0x00008000, 
	0x00010000, 
	0x00020000,
	0x00040000, 
	0x00080000,
	0x00100000,
	0x00200000,
	0x00400000,
	0x00800000,
	0x01000000,
	0x02000000,
	0x04000000,
	0x08000000,
	0x10000000,
	0x20000000,
	0x40000000,
	0x80000000,
};
#endif
static const int shift_map2[SMALL_SHIFT_MAP][3] = {
						{0,0,0},
						{+1,0,0},
						{0,+1,0},
						{0,0,+1},
						{-1,0,0},
						{0,-1,0},
						{0,0,-1}, };


/// Marching cubes tables
static const int edgeTable[256] = {
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

// xg-hside = 0
// xg+hside = 1
// yg-hside = 2
// yg+hside = 3
// zg-hside = 4
// zg+hside = 5

static const int mulTable[8][3] =
	{{0,2,4},
	{1,2,4},
	{1,3,4},
	{0,2,4},
	{0,2,5},
	{1,2,5},
	{1,3,5},
	{0,3,5}};

	// Marching cubes vertex indices and edges convention; vertInd is indicated in the cube
	//
	//		   v4_________e4_________v5
	//			/|  				/|
	//		e7 / |                 / |
	//		  /  |             e5 /  |
	//		 /   | e8            /	 | e9
	//	  v7/____|_____e6_______/v6	 |
	//		|	 |              |	 |
	//	 	|  v0|______e0______|____|v1
	//	e11 |	/        		|   /
	//		|  /			e10	|  / 
	//		| /	e3				| / e1
	//		|/					|/
	//	  v3/_________e2________/v2
	//            

static const int triTable[256][16] =
	{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
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

#define MAX_POLY_DEGREE 10 // max polynomial degree supported by the bgp projection routine (6 is needed for skin, 4 by Connolly)


/////////////////////////////////SIMPLE CLASSES ///////////////////////////////////////


/**@brief Timer: timer class. If chrono is defined use high accuracy timer*/
class Timer
{
	private:

	#ifdef ENABLE_BOOST_CHRONO
		boost::chrono::high_resolution_clock::time_point start_t,end_t;
	#else
		time_t start_t, end_t;
	#endif

	public:

		Timer()
		{
		}

		void start()
		{
			#ifdef ENABLE_BOOST_CHRONO
				start_t = boost::chrono::high_resolution_clock::now();
			#else
				time(&start_t);
			#endif
		}

		double stop()
		{
			#ifdef ENABLE_BOOST_CHRONO
				end_t = boost::chrono::high_resolution_clock::now();
				boost::chrono::milliseconds ms = boost::chrono::duration_cast<boost::chrono::milliseconds>(end_t-start_t);
				return ms.count()/1000.;
			#else
				time(&end_t);
				return difftime(end_t, start_t);
			#endif
		}
};


////////////////// Memory management routines //////////////////////

template<class T> void cleanDelete(T *&t)
{
	if (t != NULL)
	{
		delete t;
		t = NULL;
	}
	else
	{
		cout << endl << WARN << "Attempting to de-allocate a null object!";
	}
}

template<class T> void cleanDeleteV(T *&t)
{
	if (t != NULL)
	{
		delete[] t;
		t = NULL;
	}
	else
	{
		cout << endl << WARN << "Attempting to de-allocate a null vector object!";
	}
}

template<class T> inline T read4DVector(const T *const var,const int64_t i,const int64_t j,const int64_t k,const int64_t l,const int64_t nx,const int64_t ny,const int64_t nz,const int64_t nl)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || l>=nl || i<0 || j<0 || k<0 || l<0)
	{
		cout << endl << ERR << "Out of bound error in reading 4DVector";
		exit(-1);
	}
	#endif
	
	// The following is better than: index = i+nx*(j+ny*(k+nz*l)):
	// the following allows computations' interdependence/parallelism
	// whereas (ny*nx*nl), (nx*nl) and nl are deemed constant by the compiler (and computed once before looping...)
	int64_t index = k*ny*nx*nl + j*nx*nl + i*nl + l;
	return var[index];
}

template<class T> inline void write4DVector(T *const var,const T val,const int64_t i,const int64_t j,const int64_t k,const int64_t l,const int64_t nx,const int64_t ny,const int64_t nz,const int64_t nl)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || l>=nl || i<0 || j<0 || k<0 || l<0)
	{
		cout << endl << ERR << "Out of bound error in writing 4DVector";
		exit(-1);
	}
	#endif

	// The following is better than: index = i+nx*(j+ny*(k+nz*l)):
	// the following allows computations' interdependence/parallelism
	// while (ny*nx*nl), (nx*nl) and nl are deemed constant by the compiler (and computed once before looping...)
	int64_t index = k*ny*nx*nl + j*nx*nl + i*nl + l;
	var[index] = val;
}

template<class T> inline T read3DVector(const T *const var,const int64_t i,const int64_t j,const int64_t k,const int64_t nx,const int64_t ny,const int64_t nz)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in reading 3DVector";
		exit(-1);
	}
	#endif

	int64_t index = k*ny*nx + j*nx + i;
	return var[index];
}

template<class T> inline void write3DVector(T *const var,const T val,const int64_t i,const int64_t j,const int64_t k,const int64_t nx,const int64_t ny,const int64_t nz)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in writing 3DVector";
		exit(-1);
	}
	#endif

	int64_t index = k*ny*nx + j*nx + i;
	var[index] = val;
}

/** This function detemines the number of coarse macro-cells of a bilevel grid (along a direction) */
inline int64_t getCoarseN(const int64_t N)
{
	int64_t coarse_N = N >> 2;
	if ((coarse_N << 2) < N) ++coarse_N;

	return coarse_N;
}

/** This function detemines the coarse macro-cell ID of a bilevel grid from the global coordinate ID */
inline int getCoarseID(const int64_t coarse_N, const int64_t id)
{
	return id >> 2;
}

/** This function detemines the fine mini-cell ID of a bilevel grid from the global coordinate ID, id, and the previously computed macro-cell ID, coarse_id */
inline int getFineID(const int64_t coarse_id, const int64_t id)
{
	return id - (coarse_id << 2);
}

/** This function detemines the unrolled fine mini-cell ID of a bilevel grid from the 3 coordinates-wise IDs */
inline int getUnrolledFineID(const int64_t fine_x, const int64_t fine_y, const int64_t fine_z)
{
	return (fine_z<<4) | (fine_y<<2) | fine_x;
}

/** This function allocates data of bilevel grids in which each large coarse cell is comprised of
4^3 mini-cells. */
template<class T> T **allocateBilevelGridCells(const int64_t nx,const int64_t ny,const int64_t nz,const int64_t nl)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz * nl;

	T **var = (T **)malloc(sizeof(T *) * tot);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate bilevel grid";
		return NULL;
	}
	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		var[coarse_index] = NULL;
	}
	return var;
}

/** This function copies data of a bilevel grid into the small cells of another bilevel grid */
template<class T> T **copyBilevelGridCells(T **src,const int64_t nx,const int64_t ny,const int64_t nz,const int64_t nl)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz * nl;

	T **var = (T **)malloc(sizeof(T *) * tot);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate and copy bilevel grid";
		return NULL;
	}
	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		var[coarse_index] = NULL;

		if (src[coarse_index] != NULL)
		{
			var[coarse_index] = (T *)malloc(sizeof(T)*4*4*4);

			if (var[coarse_index] == NULL)
			{
				cout << endl << ERR << "Not enough memory to allocate and copy bilevel grid cells";
				return NULL;
			}
			for (int i = 0; i < 4*4*4; i++)
				var[coarse_index][i] = src[coarse_index][i];
		}
	}
	return var;
}


/** This function deallocates the data of a bilevel grid */
template<class T> void deleteBilevelGridCells(T **var,const int64_t nx,const int64_t ny,const int64_t nz,const int64_t nl)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz * nl;

	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		if (var[coarse_index] != NULL)
		{
			free(var[coarse_index]);
			var[coarse_index] = NULL;
		}
	}
}

/** This function allocates data of bilevel grids in which each large coarse cell is comprised of
4^3 mini-cells. */
template<class T> T **allocateBilevelGridCells(const int64_t nx,const int64_t ny,const int64_t nz)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz;

	T **var = (T **)malloc(sizeof(T *) * tot);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate bilevel grid";
		return NULL;
	}
	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		var[coarse_index] = NULL;
	}
	return var;
}

template<class T> T *allocateBilevelMinigridCells(const T ref_val)
{
	T *var = (T *)malloc(sizeof(T)*4*4*4);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate and write bilevel grid cells";
	}
	for (int l = 0; l < 4*4*4; l++)
		var[l] = ref_val;

	return var;
}

/** This function copies data of a bilevel grid into the small cells of another bilevel grid */
template<class T> T **copyBilevelGridCells(T **src,const int64_t nx,const int64_t ny,const int64_t nz)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz;

	T **var = (T **)malloc(sizeof(T *) * tot);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate and copy bilevel grid";
		return NULL;
	}
	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		var[coarse_index] = NULL;

		if (src[coarse_index] != NULL)
		{
			var[coarse_index] = (T *)malloc(sizeof(T)*4*4*4);

			if (var[coarse_index] == NULL)
			{
				cout << endl << ERR << "Not enough memory to allocate and copy bilevel grid cells";
				return NULL;
			}
			for (int64_t i = 0; i < 4*4*4; i++)
				var[coarse_index][i] = src[coarse_index][i];
		}
	}
	return var;
}

/** This function deallocates the data of a bilevel grid */
template<class T> void deleteBilevelGridCells(T **var,const int64_t nx,const int64_t ny,const int64_t nz)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz;

	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		if (var[coarse_index] != NULL)
		{
			free(var[coarse_index]);
			var[coarse_index] = NULL;
		}
	}
}

/** This function reads the datum of a small cell of a bilevel grid according to the
global cell coordinates */
template<class T> inline T readBilevelGrid(T **var,const T ref_val,
										   const int64_t i,const int64_t j,const int64_t k,
										   const int64_t nx,const int64_t ny,const int64_t nz)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in reading BilevelGrid";
		exit(-1);
	}
	#endif

	int64_t coarse_i = i >> 2;
	int64_t coarse_j = j >> 2;
	int64_t coarse_k = k >> 2;
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;

	int64_t coarse_index = coarse_k*coarse_ny*coarse_nx + coarse_j*coarse_nx + coarse_i;

	if (var[coarse_index] == NULL)
	{
		return ref_val;
	}
	int64_t fine_i = i - (coarse_i << 2);
	int64_t fine_j = j - (coarse_j << 2);
	int64_t fine_k = k - (coarse_k << 2);
	int64_t fine_index = (fine_k<<4) | (fine_j<<2) | fine_i;

	return var[coarse_index][fine_index];
}

/** This function writes the datum of a small cell within a large cell of a bilevel grid
according to the global cell coordinates. It is thread safe only if the algorithm avoids
the cases in whuch multiple threads are not allowed to write to the same large cell
which was not previously allocated */
template<class T> inline void writeBilevelGrid(T **var,const T ref_val,const T val,
											   const int64_t i,const int64_t j,const int64_t k,
											   const int64_t nx,const int64_t ny,const int64_t nz)
{
	if (val == ref_val)
	{
		// it is not necessary to write a background value if not previously allocated and modified
		return;
	}
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in writing BilevelGrid";
		exit(-1);
	}
	#endif

	int64_t coarse_i = i >> 2;
	int64_t coarse_j = j >> 2;
	int64_t coarse_k = k >> 2;
	int64_t coarse_nx = nx >> 2;
	int coarse_ny = ny >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;

	int64_t coarse_index = coarse_k*coarse_ny*coarse_nx + coarse_j*coarse_nx + coarse_i;

	if (var[coarse_index] == NULL)
	{
		var[coarse_index] = (T *)malloc(sizeof(T)*4*4*4);

		if (var[coarse_index] == NULL)
		{
			cout << endl << ERR << "Not enough memory to allocate and write bilevel grid cells";
		}
		for (int i = 0; i < 4*4*4; i++)
			var[coarse_index][i] = ref_val;
	}
	int64_t fine_i = i - (coarse_i<<2);
	int64_t fine_j = j - (coarse_j<<2);
	int64_t fine_k = k - (coarse_k<<2);
	int64_t fine_index = (fine_k<<4) | (fine_j<<2) | fine_i;

	var[coarse_index][fine_index] = val;
}

/** This function allocates data of bilevel grids in which each large coarse cell will be comprised of
4^3/8 mini-cells. */
inline unsigned int **allocate8xCompressedBilevelGridCells(const int64_t nx,const int64_t ny,const int64_t nz,const int64_t nl)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz * nl;

	unsigned int **var = (unsigned int **)malloc(sizeof(unsigned int *) * tot);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate 8x compressed bilevel grid";
		return NULL;
	}
	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		var[coarse_index] = NULL;
	}
	return var;
}

/** This function deallocates the data of a 8x compressed bilevel grid */
inline void delete8xCompressedBilevelGridCells(unsigned int **var,const int64_t nx,const int64_t ny,const int64_t nz,const int64_t nl)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz * nl;

	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		if (var[coarse_index] != NULL)
		{
			free(var[coarse_index]);
			var[coarse_index] = NULL;
		}
	}
}

/** This function reads the datum of a small cell of a compressed, bilevel grid according to the
global cell coordinates */
inline unsigned int read8xCompressedBilevelGrid(unsigned int **var,const unsigned int ref_val,
												const int64_t i,const int64_t j,const int64_t k,const int64_t l,
												const int64_t nx,const int64_t ny,const int64_t nz,const int64_t nl)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || l>=nl || i<0 || j<0 || k<0 || l<0)
	{
		cout << endl << ERR << "Out of bound error in reading 8xCompressedBilevelGrid";
		exit(-1);
	}
	#endif

	int64_t coarse_i = i >> 2;
	int64_t coarse_j = j >> 2;
	int64_t coarse_k = k >> 2;
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;

	int64_t coarse_index = l*coarse_nz*coarse_ny*coarse_nx + coarse_k*coarse_ny*coarse_nx + coarse_j*coarse_nx + coarse_i;

	if (var[coarse_index] == NULL) {
		return ref_val;
	}
	unsigned int fine_i = (unsigned int)(i - (coarse_i<<2));
	unsigned int fine_j = (unsigned int)(j - (coarse_j<<2));
	unsigned int fine_k = (unsigned int)(k - (coarse_k<<2));
	unsigned int fine_index = (fine_k<<4U) | (fine_j<<2U) | fine_i;

	// assumption that each var is between 0 and 15 so that there are 8 4-bit values compressed in one unsigned integer
	unsigned int cell_index = fine_index >> 3U;
	unsigned int shift_amount = fine_index - (cell_index<<3U);

	return (var[coarse_index][cell_index] >> (shift_amount<<2U)) & ((1U<<4U)-1U); // (1<<4)-1 = 15, which has unitary all the first 4 bits
}

/** This function writes the datum of a small cell within a large cell of a compressed bilevel grid
according to the global cell coordinates. It is thread safe only if the algorithm avoids
the cases in whuch multiple threads are not allowed to write to the same element of the
rightmost dimension */
inline void write8xCompressedBilevelGrid(unsigned int **var,const unsigned int ref_val,const unsigned int val,
										 const int64_t i,const int64_t j,const int64_t k,
										 const int64_t l,const int64_t nx,const int64_t ny,const int64_t nz,const int64_t nl)
{
	if (val == ref_val)
	{
		// it is not necessary to write a background value if not previously allocated and modified
		return;
	}
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || l>=nl || i<0 || j<0 || k<0 || l<0)
	{
		cout << endl << ERR << "Out of bound error in writing 8xCompressedBilevelGrid";
		exit(-1);
	}
	#endif

	int64_t coarse_i = i >> 2;
	int64_t coarse_j = j >> 2;
	int64_t coarse_k = k >> 2;
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;

	int64_t coarse_index = l*coarse_nz*coarse_ny*coarse_nx + coarse_k*coarse_ny*coarse_nx + coarse_j*coarse_nx + coarse_i;

	if (var[coarse_index] == NULL) {
		// assumption that val is between 0 and 15 so that there are 8 4-bit values compressed in one unsigned integer;
		// the number of subcells for each coarse cell is thus 4*4*4/8 = 8
		var[coarse_index] = (unsigned int *)malloc(sizeof(unsigned int)*8);

		// 8 reference values are stored in one unsigned integer starting at bits no. 0, 4, 8, ..., 28
		unsigned int compressed_ref_vals = ref_val;
		for (unsigned int i = 1U; i < 8U; i++)
			compressed_ref_vals |= ref_val << (i<<2U);

		for (int i = 0; i < 8; i++)
			var[coarse_index][i] = compressed_ref_vals;
	}
	if (var[coarse_index] == NULL)
		cout << endl << ERR << "Not enough memory to allocate cell in write8xCompressedBilevelGrid";

	unsigned int fine_i = (unsigned int)(i - (coarse_i<<2));
	unsigned int fine_j = (unsigned int)(j - (coarse_j<<2));
	unsigned int fine_k = (unsigned int)(k - (coarse_k<<2));
	unsigned int fine_index = (fine_k<<4U) + (fine_j<<2U) + fine_i;

	// assumption that val is between 0 and 15 so that there are 8 4-bit values compressed in one unsigned integer
	unsigned int cell_index = fine_index >> 3U;
	unsigned int shift_amount = fine_index - (cell_index<<3U);

	var[coarse_index][cell_index] |= val << (shift_amount<<2U);
}

/** This function allocates data of bilevel grids in which each large coarse cell will be comprised of
4^3/8 mini-cells. */
inline unsigned int **allocate8xCompressedBilevelGridCells(const int64_t nx,const int64_t ny,const int64_t nz)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz;

	unsigned int **var = (unsigned int **)malloc(sizeof(unsigned int *) * tot);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate 8x compressed bilevel grid";
		return NULL;
	}
	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		var[coarse_index] = NULL;
	}
	return var;
}

/** This function deallocates the data of a 8x compressed bilevel grid */
inline void delete8xCompressedBilevelGridCells(unsigned int **var,const int64_t nx,const int64_t ny,const int64_t nz)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz;

	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		if (var[coarse_index] != NULL)
		{
			free(var[coarse_index]);
			var[coarse_index] = NULL;
		}
	}
}

/** This function reads the datum of a small cell of a 8x compressed, bilevel grid according to the
global cell coordinates */
inline unsigned int read8xCompressedBilevelGrid(unsigned int **var,const unsigned int ref_val,
												const int64_t i,const int64_t j,const int64_t k,
												const int64_t nx,const int64_t ny,const int64_t nz)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in reading 8xCompressedBilevelGrid";
		exit(-1);
	}
	#endif

	int64_t coarse_i = i >> 2;
	int64_t coarse_j = j >> 2;
	int64_t coarse_k = k >> 2;
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;

	int64_t coarse_index = coarse_k*coarse_ny*coarse_nx + coarse_j*coarse_nx + coarse_i;

	if (var[coarse_index] == NULL) {
		return ref_val;
	}
	unsigned int fine_i = (unsigned int)(i - (coarse_i<<2));
	unsigned int fine_j = (unsigned int)(j - (coarse_j<<2));
	unsigned int fine_k = (unsigned int)(k - (coarse_k<<2));
	unsigned int fine_index = (fine_k<<4U) | (fine_j<<2U) | fine_i;

	// assumption that each var is between 0 and 15 so that there are 8 4-bit values compressed in one unsigned integer
	unsigned int cell_index = fine_index >> 3U;
	unsigned int shift_amount = fine_index - (cell_index<<3U);

	return (var[coarse_index][cell_index] >> (shift_amount<<2U)) & ((1U<<4U)-1U); // (1<<4)-1 = 15, which has unitary all the first 4 bits
}

/** This function writes the datum of a small cell of a 8x compressed, bilevel grid according to the
global cell coordinates. It is thread safe only if the algorithm avoids
the cases in whuch multiple threads are not allowed to write to the same large cell
which was not previously allocated or to the same element of the rightmost dimension */
inline void write8xCompressedBilevelGrid(unsigned int **var,const unsigned int ref_val,const unsigned int val,
										 const int64_t i,const int64_t j,const int64_t k,
										 const int64_t nx,const int64_t ny,const int64_t nz)
{
	if (val == ref_val)
	{
		// it is not necessary to write a background value if not previously allocated and modified
		return;
	}
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in writing 8xCompressedBilevelGrid";
		exit(-1);
	}
	#endif

	int64_t coarse_i = i >> 2;
	int64_t coarse_j = j >> 2;
	int64_t coarse_k = k >> 2;
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;

	int64_t coarse_index = coarse_k*coarse_ny*coarse_nx + coarse_j*coarse_nx + coarse_i;

	if (var[coarse_index] == NULL) {
		// assumption that val is between 0 and 15 so that there are 8 4-bit values compressed in one unsigned integer;
		// the number of subcells for each coarse cell is thus 4*4*4/8 = 8
		var[coarse_index] = (unsigned int *)malloc(sizeof(unsigned int)*8);

		// 8 reference values are stored in one unsigned integer starting at bits no. 0, 4, 8, ..., 28
		unsigned int compressed_ref_vals = ref_val;
		for (unsigned int i = 1U; i < 8U; i++)
			compressed_ref_vals |= ref_val << (i<<2U);

		for (int i = 0; i < 8; i++)
			var[coarse_index][i] = compressed_ref_vals;
	}
	if (var[coarse_index] == NULL)
		cout << endl << ERR << "Not enough memory to allocate cell in 8xCompressedBilevelGrid";

	unsigned int fine_i = (unsigned int)(i - (coarse_i<<2));
	unsigned int fine_j = (unsigned int)(j - (coarse_j<<2));
	unsigned int fine_k = (unsigned int)(k - (coarse_k<<2));
	unsigned int fine_index = (fine_k<<4U) | (fine_j<<2U) | fine_i;

	// assumption that val is between 0 and 15 so that there are 8 4-bit values compressed in one unsigned integer
	unsigned int cell_index = fine_index >> 3U;
	unsigned int shift_amount = fine_index - (cell_index<<3U);

	var[coarse_index][cell_index] |= val << (shift_amount<<2U);
}

/** This function allocates data of bilevel grids in which each large coarse cell will be comprised of
4^3/32 mini-cells. */
inline unsigned int **allocate32xCompressedBilevelGridCells(const int64_t nx,const int64_t ny,const int64_t nz)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz;

	unsigned int **var = (unsigned int **)malloc(sizeof(unsigned int *) * tot);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate 32x compressed bilevel grid";
		return NULL;
	}
	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		var[coarse_index] = NULL;
	}
	return var;
}

/** This function deallocates the data of a 32x compressed bilevel grid */
inline void delete32xCompressedBilevelGridCells(unsigned int **var,const int64_t nx,const int64_t ny,const int64_t nz)
{
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	int64_t coarse_nz = nz >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;
	if ((coarse_nz << 2) < nz) ++coarse_nz;
	int64_t tot = coarse_nx*coarse_ny*coarse_nz;

	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		if (var[coarse_index] != NULL)
		{
			free(var[coarse_index]);
			var[coarse_index] = NULL;
		}
	}
}

/** This function reads the datum of a small cell of a 32x compressed, bilevel grid according to the
global cell coordinates */
inline bool read32xCompressedBilevelGrid(unsigned int **var,const bool ref_val,
										 const int64_t i,const int64_t j,const int64_t k,
										 const int64_t nx,const int64_t ny,const int64_t nz)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in reading 32xCompressedBilevelGrid";
		exit(-1);
	}
	#endif

	int64_t coarse_i = i >> 2;
	int64_t coarse_j = j >> 2;
	int64_t coarse_k = k >> 2;
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;

	int64_t coarse_index = (int64_t)coarse_k*(int64_t)(coarse_ny*coarse_nx) + (int64_t)(coarse_j*coarse_nx + coarse_i);

	if (var[coarse_index] == NULL) {
		return ref_val;
	}
	unsigned int fine_i = (unsigned int)(i - (coarse_i<<2U));
	unsigned int fine_j = (unsigned int)(j - (coarse_j<<2U));
	unsigned int fine_k = (unsigned int)(k - (coarse_k<<2U));
	unsigned int fine_index = (fine_k<<4U) |  (fine_j<<2U) | fine_i;

	// assumption that each var is boolean so that there are 32 1-bit values compressed in one unsigned integer
	unsigned int cell_index = fine_index >> 5U;
	unsigned int shift_amount = fine_index - (cell_index<<5U);

	// extract and shift the pertinent bit
	return (bool)((var[coarse_index][cell_index] >> shift_amount) & 1U);
}

/** This function writes the datum of a small cell of a 8x compressed, bilevel grid according to the
global cell coordinates. It is thread safe only if the algorithm avoids                *
the cases in whuch multiple threads are not allowed to write to the same large cell
which was not previously allocated or to the same element of the rightmost dimension */
inline void write32xCompressedBilevelGrid(unsigned int **var,const bool ref_val,const bool val,
										  const int64_t i,const int64_t j,const int64_t k,
										  const int64_t nx,const int64_t ny,const int64_t nz)
{
	if (val == ref_val)
	{
		// it is not necessary to write a background value if not previously allocated and modified
		return;
	}
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in writing 32xCompressedBilevelGrid";
		exit(-1);
	}
	#endif

	int64_t coarse_i = i >> 2;
	int64_t coarse_j = j >> 2;
	int64_t coarse_k = k >> 2;
	int64_t coarse_nx = nx >> 2;
	int64_t coarse_ny = ny >> 2;
	if ((coarse_nx << 2) < nx) ++coarse_nx;
	if ((coarse_ny << 2) < ny) ++coarse_ny;

	int64_t coarse_index = coarse_k*(coarse_ny*coarse_nx) + coarse_j*coarse_nx + coarse_i;

	if (var[coarse_index] == NULL) {
		// val is boolean so that there are 32 1-bit values compressed in one unsigned integer;
		// the number of subcells for each coarse cell is thus 4*4*4/32 = 2
		var[coarse_index] = (unsigned int *)malloc(sizeof(unsigned int)*2);

		if (ref_val)
		{
			for (int i = 0; i < 2; i++)
				var[coarse_index][i] = ~0U;
		}
		else
		{
			for (int i = 0; i < 2; i++)
				var[coarse_index][i] = 0U;
		}
	}
	if (var[coarse_index] == NULL)
		cout << endl << ERR << "Not enough memory to allocate cell in 32xCompressedBilevelGrid";

	unsigned int fine_i = (unsigned int)(i - (coarse_i<<2));
	unsigned int fine_j = (unsigned int)(j - (coarse_j<<2));
	unsigned int fine_k = (unsigned int)(k - (coarse_k<<2));
	unsigned int fine_index = (fine_k<<4U) | (fine_j<<2U) | fine_i;

	// val is boolean so that there are 32 1-bit values compressed in one unsigned integer
	unsigned int cell_index = fine_index >> 5U;
	unsigned int shift_amount = fine_index - (cell_index<<5U);

	if (val)
		// switch ON the pertinent bit (if val is true)
		var[coarse_index][cell_index] |= (unsigned int)val << shift_amount;
	else
		// switch OFF the pertinent bit
		var[coarse_index][cell_index] &= ~(1U << shift_amount);
}

/** This function allocates the data of a 32x compressed grid */
inline unsigned int *allocate32xCompressedGrid(bool const ref_val,const int64_t nx,const int64_t ny,const int64_t nz)
{
	int64_t compressed_nx = nx >> 5L;
	if ((compressed_nx << 5L) < nx) ++compressed_nx;
	int64_t tot = compressed_nx*ny*nz;

	unsigned int *var = (unsigned int *)malloc(sizeof(unsigned int)*tot);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate 32x compressed grid ";
		return NULL;
	}
	if (ref_val)
	{
		for (int64_t i=0; i<tot; i++)
			var[i] = ~0U;
	}
	else
	{
		for (int64_t i=0; i<tot; i++)
			var[i] = 0U;
	}
	return var;
}

/** This function allocates the data of a 32x compressed bilevel grid */
inline void delete32xCompressedGrid(unsigned int **var,const int64_t nx,const int64_t ny,const int64_t nz)
{
	int64_t compressed_nx = nx >> 5L;
	if ((compressed_nx << 5L) < nx) ++compressed_nx;
	int64_t tot = compressed_nx*ny*nz;

	for (int64_t coarse_index = 0; coarse_index < tot; coarse_index++)
	{
		if (var[coarse_index] != NULL)
		{
			free(var[coarse_index]);
			var[coarse_index] = NULL;
		}
	}
}

/** This function reads the datum of a small cell of a compressed grid according to the
global cell coordinates */
inline bool read32xCompressedGrid(unsigned int *var,const int64_t i,const int64_t j,const int64_t k,
								  const int64_t nx,const int64_t ny,const int64_t nz)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in reading 32xCompressedGrid";
		exit(-1);
	}
	#endif

	int64_t fine_index = k*ny*nx + j*nx + i;
	int64_t cell_index = fine_index >> 5L;
	unsigned int shift_amount = (unsigned int)(fine_index - (cell_index<<5L));

	// extract and shift the pertinent bit
	return (bool)((var[cell_index] >> shift_amount) & 1U);
}

/** This function writes the datum of a cell of a compressed grid according to the
global cell coordinates. It is thread safe only if the algorithm avoids                *
the cases in whuch multiple threads are not allowed to the same element of the
rightmost dimension */
inline void write32xCompressedGrid(unsigned int *var,const bool val,
								   const int64_t i,const int64_t j,const int64_t k,
								   const int64_t nx,const int64_t ny,const int64_t nz)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in reading 32xCompressedGrid";
		exit(-1);
	}
	#endif

	int64_t fine_index = k*ny*nx + j*nx + i;
	int64_t cell_index = fine_index >> 5L;
	unsigned int shift_amount = (unsigned int)(fine_index - (cell_index<<5L));

	if (val)
		// switch ON the pertinent bit
		var[cell_index] |= (unsigned int)val << shift_amount;
	else
		// switch OFF the pertinent bit
		var[cell_index] &= ~(1U << shift_amount);
}

template<class T> inline T read2DVector(const T *const var,const int64_t i,const int64_t j,const int64_t nx,const int64_t ny)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || i<0 || j<0)
	{
		cout << endl << ERR << "Out of bound error in reading 2DVector";
		exit(-1);
	}
	#endif

	int64_t index = j*nx + i;
	return var[index];
}

template<class T> inline void write2DVector(T *const var,const T val,const int64_t i,const int64_t j,const int64_t nx,const int64_t ny)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || i<0 || j<0)
	{
		cout << endl << ERR << "Out of bound error in reading 2DVector";
		exit(-1);
	}
	#endif

	int64_t index = j*nx + i;
	var[index] = val;
}

template<class T> inline T readVector(const T *const var,const int64_t i,const int64_t nx)
{
	#ifdef CHECK_BOUNDS
	if (i >= nx || i < 0)
	{
		cout << endl << ERR << "Out of bound error in reading Vector";
		exit(-1);
	}
	#endif
	
	return var[i];
}

template<class T> inline void writeVector(T *const var,const T val,const int64_t i,const int64_t nx)
{
	#ifdef CHECK_BOUNDS
	if (i >= nx || i < 0)
	{
		cout << endl << ERR << "Out of bound error in writing Vector";
		exit(-1);
	}
	#endif
	var[index] = val;
}

template<class T> T *allocateVector(int64_t n)
{
	T *t = (T *)malloc(sizeof(T)*n);
	if (t == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate vector ";
		return NULL;
	}
	return t;
}

inline bool read32xCompressedVector(const unsigned int *const var,const int64_t i,const int64_t nx)
{
	#ifdef CHECK_BOUNDS
	if (i >= nx || i < 0)
	{
		cout << endl << ERR << "Out of bound error in reading 32x compressed vector";
		exit(-1);
	}
	#endif

	unsigned int shift_amount = (unsigned int)(i - ((i>>5L)<<5L));

	// extract and shift the pertinent bit
	return (bool)((var[ i>>5L ] >> shift_amount) & 1U);
}

inline void write32xCompressedVector(unsigned int *const var,const bool val,const int64_t i,const int64_t nx)
{
	#ifdef CHECK_BOUNDS
	if (i >= nx || i < 0)
	{
		cout << endl << ERR << "Out of bound error in writing 32x compressed vector";
		exit(-1);
	}
	#endif

	unsigned int shift_amount = (unsigned int)(i - ((i>>5L)<<5L));

	if (val)
		// switch ON the pertinent bit (if val is true)
		var[ i>>5L ] |= (unsigned int)val << shift_amount;
	else
		// switch OFF the pertinent bit
		var[ i>>5L ] &= ~(1U << shift_amount);
}

inline unsigned int *allocate32xCompressedVector(bool const ref_val,const int64_t nx)
{
	int64_t compressed_nx = nx >> 5L;
	if ((compressed_nx << 5L) < nx) ++compressed_nx;

	unsigned int *var = (unsigned int *)malloc(sizeof(unsigned int)*compressed_nx);

	if (var == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate 32x compressed vector ";
		return NULL;
	}
	if (ref_val)
	{
		for (int64_t i=0; i<compressed_nx; i++)
			var[i] = ~0U;
	}
	else
	{
		for (int64_t i=0; i<compressed_nx; i++)
			var[i] = 0U;
	}
	return var;
}

template<class T>  T **allocateMatrix2D(int64_t nrows,int64_t ncols)
{
	T **t = (T **)malloc(sizeof(T *)*nrows);
	if (t == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate 2D matrix ";
		return NULL;
	}
	for (int64_t i=0; i<nrows; i++)
	{
		t[i] = (T *)malloc(sizeof(T)*ncols);
		
		if (t[i] == NULL)
		{
			cout << endl << ERR << "Not enough memory to allocate 2D matrix ";
			return NULL;
		}	
	}
	return t;
}

template<class T>  T ***allocateMatrix3D(int64_t nx,int64_t ny,int64_t nz)
{
	T ***t = (T***)malloc(sizeof(T **)*nz);
	if (t == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate 3D matrix ";
		return NULL;
	}
	for (int64_t k=0; k<nz; k++)
	{
		t[k] = (T **)malloc(sizeof(T *)*ny);
		if (t[k] == NULL)
		{
			cout << endl << ERR << "Not enough memory to allocate 3D matrix ";
			return NULL;
		}

		for (int64_t j=0; j<ny; j++)
		{
			t[k][j] = (T *)malloc(sizeof(T)*nx);
			if (t[k][j] == NULL)
			{	
				cout << endl << ERR << "Not enough memory to allocate 3D matrix ";
				return NULL;
			}
		}
	}
	return t;
}

template<class T> void deleteVector(T *&t)
{
	if (t != NULL)
	{
		free(t);
		t = NULL;
	}
	else
	{
		cout << endl << WARN << "Attempting to de-allocate a null vector!";
	}
}

template<class T> void deleteMatrix2D(int64_t nrows,int64_t ncols,T **&t)
{
	if (t != NULL)
	{
		for (int64_t i=0; i<nrows; i++)
			free(t[i]);
		free(t);
		t = NULL;
	}
	else
	{
		cout << endl << WARN << "Attempting to de-allocate a null vector!";
	}
}

template<class T> void deleteMatrix3D(int64_t nx,int64_t ny,int64_t nz,T ***&t)
{
	if (t != NULL)
	{
		for (int64_t k=0; k<nz; k++)
		{
			for (int64_t j=0; j<ny; j++)
				free(t[k][j]);
			free(t[k]);
		}
		free(t);
		t = NULL;
	}
	else
	{
		cout << endl << WARN << "Attempting to de-allocate a null vector!";
	}
}

///////////////////////////////////////////////////////////////////

/** a vector bool is a vector in which it is assumed that the information is written/read in a bitwise way
to minimize memory footprint. Thus the allocator computes how many 32 bits words are needed to store n elements*/
#ifdef ENABLE_BOOST_THREADS
inline boost::uint32_t *allocateVectorBool(size_t n)
{
	size_t nw = n >> 5L;

	if ((nw<<5L) != n)
		nw++;
	boost::uint32_t *ptr = (boost::uint32_t*)malloc(nw*4);

	if (ptr == NULL)
		cout << endl << ERR << "Not enough memory to allocate bit vector";	
	return ptr;
}

inline void deleteVectorBool(boost::uint32_t *&t)
{
	if (t != NULL)
	{
		free(t);
		t = NULL;
	}
	else
	{
		cout << endl << WARN << "Attempting to de-allocate a null bool vector!";
	}
}

inline bool read3DVectorBool(boost::uint32_t *const var,const int i,const int64_t j,const int64_t k,const int64_t nx,const int64_t ny,const int64_t nz)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in reading 3DVector";
		exit(-1);
	}
	#endif

	int64_t index = k*ny*nx + j*nx + i;
	int64_t word_index = index >> 5L;
	int64_t position = index - (word_index<<5L);
	boost::uint32_t mask = MASK[position];
	return ((var[word_index] & mask) == mask);
}

inline void write3DVectorBool(boost::uint32_t *const var,const bool val, const int64_t i,const int64_t j,const int64_t k,const int64_t nx,const int64_t ny,const int64_t nz)
{
	#ifdef CHECK_BOUNDS
	if (i>=nx || j>=ny || k>=nz || i<0 || j<0 || k<0)
	{
		cout << endl << ERR << "Out of bound error in reading 3DVector";
		exit(-1);
	}
	#endif

	int64_t index = k*ny*nx + j*nx + i;
	int64_t word_index = index >> 5L;
	int64_t position = index - (word_index<<5L);
	boost::uint32_t mask = MASK[position];
	var[word_index] = val ? (var[word_index] | mask) : (var[word_index] & ~mask) ;
}
#endif // ENABLE_BOOST_THREADS

/** portable aligned memory allocator.You will lose a bit of memory. For big
arrays that are heavily accessed this is the best allocation choice. For small
arrays very often allocated/deallocated this can introduce a signifcant memory overhead.
Free only using deleteAlignedVector. Typical bit alignement values are 32,64,128,256*/
template<class T> T *allocateAlignedVector(size_t n,int bitAlignment)
{
    T *mem = (T *)malloc(n*sizeof(T) + bitAlignment + sizeof(void*));
    T **ptr = (T **)((int64_t)(mem + bitAlignment + sizeof(void*)) & ~(bitAlignment-1));
    ptr[-1] = mem;

	if (mem == NULL)
	{
		cout << endl << ERR << "Not enough memory to allocate aligned vector ";
		return NULL;
	}

    return (T *)ptr;
}

template<class T> void deleteAlignedVector(T *&ptr)
{
	if (ptr == NULL)
	{
		cout << endl << WARN << "Attempting to de-allocate a null aligned vector!";
		return;
	}

    free(((T **)ptr)[-1]);
	ptr = NULL;
}


inline double accurateSum3(double term1, double term2, double term3)
{
	double term[] = {term1, term2, term3};

	int id[3];

	id[0] = (fabs(term[   0 ]) < fabs(term[1])) ? 0     : 1;
	id[0] = (fabs(term[id[0]]) < fabs(term[2])) ? id[0] : 2;
	id[2] = (fabs(term[   0 ]) > fabs(term[1])) ? 0     : 1;
	id[2] = (fabs(term[id[2]]) > fabs(term[2])) ? id[2] : 2;

	id[1] = 3 - id[0] - id[2];

	return (term[id[0]] + term[id[1]]) + term[id[2]];
}


inline double accurateSum4(double term1, double term2, double term3, double term4)
{
	double term[] = {term1, term2, term3, term4};

	int id[4];

	id[0] = (fabs(term[   0 ]) < fabs(term[1])) ? 0     : 1;
	id[0] = (fabs(term[id[0]]) < fabs(term[2])) ? id[0] : 2;
	id[0] = (fabs(term[id[0]]) < fabs(term[3])) ? id[0] : 3;
	id[3] = (fabs(term[   0 ]) > fabs(term[1])) ? 0     : 1;
	id[3] = (fabs(term[id[3]]) > fabs(term[2])) ? id[3] : 2;
	id[3] = (fabs(term[id[3]]) > fabs(term[3])) ? id[3] : 3;

	if (id[0] + id[3] == 1)
	{
		id[1] = 2;
		id[2] = 3;
	}
	else if (id[0] + id[3] == 2)
	{
		id[1] = 1;
		id[2] = 3;
	}
	else if (id[0] * id[3] == 2)
	{
		id[1] = 0;
		id[2] = 3;
	}
	else if (id[0] + id[3] == 4)
	{
		id[1] = 0;
		id[2] = 2;
	}
	else if (id[0] + id[3] == 5)
	{
		id[1] = 0;
		id[2] = 1;
	}
	else
	{
		id[1] = 1;
		id[2] = 2;
	}
	if (fabs(term[id[1]]) > fabs(term[id[2]]))
	{
		int temp = id[2];
		id[2] = id[1];
		id[1] = temp;
	}
	return ((term[id[0]] + term[id[1]]) + term[id[2]]) + term[id[3]];
}


///////////////////// STL comparators //////////////////////////////

/**@brief ascending on first VERTEX_TYPE of pair<VERTEX_TYPE,VERTEX_TYPE*> comparator */
bool compKeepIndex(pair<VERTEX_TYPE,VERTEX_TYPE*> a, pair<VERTEX_TYPE,VERTEX_TYPE*> b);

///////////////////////////////////////////////////////////////////

typedef std::pair<double,int> indexed_double;
bool index_double_comparator( const indexed_double &l, const indexed_double &r);

///////////////////// MISCELLANEAOUS ////////////////////////////////

/** @brief test if point proj is in cube whose center is point and side is side variable.
 If a toll is provided, the test is performed up to a tollerance value.
 By default the tollerance is 0*/
bool testInCube(const double proj0,const double proj1,const double proj2,const double point0,const double point1,const double point2,const double side,const double toll=0);
double rintp(const double x);
string toLowerCase(string str);
void cleanLine();

static int64_t SEED = 1234;
// randnumber between 0 and 1
double randnum();

/** get the real roots by computing the companion matrix and then extracting the eigenvalues. A root is real
if its imaginary part in absolute value is less than a given threshold. Usually this threshold is rather conservative
such that possibly unprecise real roots are not lost. Poly is modified in place to remove possible zeros*/
void getRealRootsCompanion(double *const poly,const int degree,double *const roots,int &numroots);
/** get real roots by using Sturm method. Directly search real roots. Much faster, often less accurate
than companion matrix*/
void getRealRootsSturm(const double *const polyy,const int degree,double *const roots,int &numrootss);
/** plane by 3 points routine*/
void plane3points(const double p1[3],const double p2[3],const double p3[3],double w[4],const bool normalize=true);
/** multiply two planes to get a quadratic for from them*/
void planeByplane(const double pl[4],const double p2[4], double Q[3][3],double a[3],double &c);
/** point to plane projection*/
void point2plane(const double p[3],double w[4],double *const dist, double proj[3]);
/** in place inversion of 4x4 matrix*/
void inplace_invert4x4(double M[4][4]);
/** unrolled 4x4 matrix multiply*/
void Matrix4x4MultiplyBy4x4 (const double src1[4][4],const double src2[4][4], double dest[4][4]);
double realCubicSolution(double b, double c, double d);
void quarticEqSolutions(double roots[4], double b, double c, double d, double e, int *num_sol);
/** ray/sphere intersection. ray is o+t*dir */
bool raySphere(const double *const orig,const double *const dir,const double *const sphere_center,const double sphere_radius,double *const t1,double *const t2);
/** ray vs torus intersection. Ray is o+t*dir */
bool rayTorus(int analytical, double invrot[3][3],double torus_center[3],double sphere_center[3],double probe_radius,double major_radius,double radius, int panel,double orig[3],double dir[3],double t[4], int &numInt);
/** get the normal to a sphere*/
void getNormalToSphere(const double *const y,const double *const center,const double radius,double *const normal);
/** project y to sphere and returns the normal*/
void projectToSphere(const double *const y,const double *const center,const double radius,double *const proj,double &dist);

void getMemSpace (double &current_mem_in_MB, double &peak_mem_in_MB);

#if defined(ENABLE_CGAL) && defined(CGAL_LINKED_WITH_TBB)
void initTBB (int num_threads);
#endif


#endif
