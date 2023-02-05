
//---------------------------------------------------------
/*    @file		globals.h
*     @brief	globals.h Includes all the common global parameters
*							    							*/
//---------------------------------------------------------

#ifndef globals_h
#define globals_h

/** @brief enable/disable Microsoft C++ compiler memory leak detection */
//#define DBGMEM_CRT

#ifdef DBGMEM_CRT
	// check vector bounds. To allow that you have to access vector following the
	// here defined access functions
	#define CHECK_BOUNDS
#endif

//if Microsoft C++ is not idenfied memory leak is deactivated
#ifndef _MSC_VER 
	#undef DBGMEM_CRT
#endif

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
	#pragma message(" MESSAGE Microsoft Memory Leaks Detector is Enabled!")
	#include <stdlib.h>
	#include <crtdbg.h>
#endif

#define VERSION "0.7.5"
#define PROGNAME "NanoShaper" 

//////////////////// include section ///////////////
#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include<assert.h>
#include<algorithm>
#include<stdio.h>
#include<cstring>
#include<map>
#include<list>
#include<set>
#include<time.h>
#include<queue>
#include<limits>
#include "ConfigFile.h"
///////////////////////////////////////////////////

using namespace std;

//////////////// Global macros and constants /////////////
#define WARN " <<WARNING>> "
#define ERR " <<ERROR>> "
#define INFO " <<INFO>> "
#define CITE " <<CITATION>> "
#define REMARK " <<REMARK>> "
//////////////////////////////////////////////////////////

// access to EPSMAP,IDEBMAP,STATUSMAP colum-major as in Fortran
// substituted by template access/writing inline functions with optional bounds check

//#define EPSMAP(i,j,k,l,NX,NY,NZ) epsmap[(i)+(NX)*((j)+(NY)*((k)+(NZ)*(l)))]
//#define IDEBMAP(i,j,k,NX,NY) idebmap[(i)+(NX)*((j)+(NY)*(k))]
//#define STATUSMAP(i,j,k,NX,NY) status[(i)+(NX)*((j)+(NY)*(k))]
//#define TEMP_STATUSMAP(i,j,k,NX,NY) tempStatus[(i)+(NX)*((j)+(NY)*(k))]
//#define TEMP_STATUSMAP2(i,j,k,NX,NY) tempStatus2[(i)+(NX)*((j)+(NY)*(k))]

//////////////// Compilation flags ////////////////////////
/** @brief if defined boost threading is enabled */
//#define ENABLE_BOOST_THREADS
//#define ENABLE_BOOST_CHRONO
/** @brief if defined cgal for skin surface is enabled*/
//#define ENABLE_CGAL

//////////////////////////////////////////////////////////////////

#ifdef ENABLE_BOOST_THREADS
	#include <boost/thread/thread.hpp>
	#include <boost/thread/mutex.hpp>
	#include <boost/thread/condition_variable.hpp>
#endif

#ifdef ENABLE_BOOST_CHRONO
	#include <boost/chrono.hpp>
#endif

#ifndef INFINITY
	#define INFINITY 1e100
#endif


///////////////////////////////////////////////////

/** char buffer length */
#define BUFLEN 1024

#define PI 3.14159265
#define HALF_PI 1.57079633
#define QUARTER_PI 0.785398163
#define TWO_PI 6.28318531

////////////////////// arithmetic and R3 vectors tools //////////////
/** absolute value */
#define ABS(x) (x < 0 ? -(x) : (x))

/** cross product */
#define CROSS(dest,v1,v2) \
		dest[0]=(v1[1])*(v2[2])-(v1[2])*(v2[1]); \
		dest[1]=(v1[2])*(v2[0])-(v1[0])*(v2[2]); \
		dest[2]=(v1[0])*(v2[1])-(v1[1])*(v2[0]);

/** dot product */
#define DOT(v1,v2) ((v1[0])*(v2[0])+(v1[1])*(v2[1])+(v1[2])*(v2[2]))

/** normalization routine, t is a temporary. Assumes sqrt from math.h
is available*/
#define NORMALIZE(v,t) \
		t=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); \
		v[0]=v[0]/t; \
		v[1]=v[1]/t; \
		v[2]=v[2]/t;

#define NORMALIZE_PLANE(v,t) \
		t=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); \
		v[0]=v[0]/t; \
		v[1]=v[1]/t; \
		v[2]=v[2]/t; \
		v[3]=v[3]/t;

#define NORMALIZE_S(v,t) \
		t=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); \
		v[0]=v[0]/(t+1e-20); \
		v[1]=v[1]/(t+1e-20); \
		v[2]=v[2]/(t+1e-20);

#define NORMALIZE_S_ASSIGN(w,v,t) \
		t=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); \
		w[0]=v[0]/(t+1e-20); \
		w[1]=v[1]/(t+1e-20); \
		w[2]=v[2]/(t+1e-20);

/** R3 vector copy */
#define ASSIGN(u,v) \
		{ u[0]=v[0]; \
		u[1]=v[1]; \
		u[2]=v[2]; }

/** R4 vector copy */
#define ASSIGN4(u,v) \
		{ u[0]=v[0]; \
		u[1]=v[1]; \
		u[2]=v[2];\
		u[3]=v[3];}

/** invert sign to an R3 vector*/
#define CHANGE_SIGN(u) \
		{	u[0]=-u[0]; \
			u[1]=-u[1]; \
			u[2]=-u[2]; }

/** R3 substraction routine */
#define SUB(dest,v1,v2)\
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2];

/** R3 add routine */
#define ADD(dest,v1,v2)\
          dest[0]=v1[0]+v2[0]; \
          dest[1]=v1[1]+v2[1]; \
          dest[2]=v1[2]+v2[2];

/** R3 mid point routine */
#define MID(dest,v1,v2)\
          dest[0]=(v1[0]+v2[0])*0.5; \
          dest[1]=(v1[1]+v2[1])*0.5; \
          dest[2]=(v1[2]+v2[2])*0.5;

/** R3 mul acc routine */
#define ADD_MUL(dest,v1,v2,a)\
          dest[0]=v1[0]+a*v2[0]; \
          dest[1]=v1[1]+a*v2[1]; \
          dest[2]=v1[2]+a*v2[2];

/** R3 mul acc routine */
#define VEC_MUL(dest,v,a)\
          dest[0]=a*v[0]; \
          dest[1]=a*v[1]; \
          dest[2]=a*v[2];

#define DIST(dist,v1,v2)\
          dist=sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2])); 

#define DIST2(dist,v1,v2)\
          dist=((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2])); 

/** min value */
#define MIN(x, y) ((x)<(y) ? (x):(y))

/** max value */
#define MAX(x, y) ((x)>(y) ? (x):(y))


/*	determinant of matrix
	Computes determinant of matrix m, returning d
 */
#define DETERMINANT_2X2(d,m)						\
{													\
	d = m[0][0] * m[1][1] - m[0][1] * m[1][0];		\
}

/*	determinant of matrix
	Computes determinant of matrix m, returning d
 */

#define DETERMINANT_3X3(d,m)								\
{															\
	d = m[0][0] * (m[1][1]*m[2][2] - m[1][2] * m[2][1]);    \
	d -= m[0][1] * (m[1][0]*m[2][2] - m[1][2] * m[2][0]);   \
	d += m[0][2] * (m[1][0]*m[2][1] - m[1][1] * m[2][0]);   \
}

/* compute adjoint of matrix and scale
 Computes adjoint of matrix m, scales it by s, returning a
 */

#define SCALE_ADJOINT_2X2(a,s,m) \
{                               \
	a[0][0] = (s) * m[1][1];    \
	a[1][0] = - (s) * m[1][0];  \
	a[0][1] = - (s) * m[0][1];  \
	a[1][1] = (s) * m[0][0];    \
}

 
/*	compute adjoint of matrix and scale
	Computes adjoint of matrix m, scales it by s, returning a
 */

#define SCALE_ADJOINT_3X3(a,s,m)	\
{									\
	a[0][0] = (s) * (m[1][1] * m[2][2] - m[1][2] * m[2][1]); \
	a[1][0] = (s) * (m[1][2] * m[2][0] - m[1][0] * m[2][2]); \
	a[2][0] = (s) * (m[1][0] * m[2][1] - m[1][1] * m[2][0]); \
	\
	a[0][1] = (s) * (m[0][2] * m[2][1] - m[0][1] * m[2][2]); \
	a[1][1] = (s) * (m[0][0] * m[2][2] - m[0][2] * m[2][0]); \
	a[2][1] = (s) * (m[0][1] * m[2][0] - m[0][0] * m[2][1]); \
	\
	a[0][2] = (s) * (m[0][1] * m[1][2] - m[0][2] * m[1][1]); \
	a[1][2] = (s) * (m[0][2] * m[1][0] - m[0][0] * m[1][2]); \
	a[2][2] = (s) * (m[0][0] * m[1][1] - m[0][1] * m[1][0]); \
}

/*	inverse of matrix 
	Compute inverse of matrix a, returning determinant m and inverse b
 */

#define INVERT_2X2(b,det,a)			\
{									\
	double tmp;						\
	DETERMINANT_2X2 (det, a);		\
	tmp = 1.0 / (det);              \
	SCALE_ADJOINT_2X2 (b, tmp, a);  \
}

/*	inverse of matrix 
	Compute inverse of matrix a, returning determinant m and inverse b
 */

#define INVERT_3X3(b,det,a)		\
{								\
   double tmp;					\
   DETERMINANT_3X3 (det, a);    \
   tmp = 1.0 / (det);           \
   SCALE_ADJOINT_3X3 (b, tmp, a); \
}

#define PRINT_MAT3(mat)	\
{ \
	cout << endl;					\
	for (int i=0;i<3;i++) \
	{						\
		cout << mat[i][0] << "," << mat[i][1] << "," << mat[i][2];\
		cout << endl;	\
	} \
}

#define PRINT_VEC3(vec)	\
{ \
	cout << endl;					\
	cout << vec[0] << "," << vec[1] << "," << vec[2];\
	cout << endl;	\
	\
}
////////////////////////////////////////////////////////

extern fstream* errorStream;
extern fstream* internals;

class Configuration
{
	public:
	
		double cavVol;
		int numMol;
	
		// grid (DelPhi) params
		double scale;
		double perfill;
		// mol file name
		string molFile;
		// sys name
		string sysName;

		bool multi_diel;
				
		// actions
		bool fillCavities;
		bool buildEpsmaps;
		bool buildStatus;
		bool tri;	
		bool accTri;
		bool smoothing;	
		bool tri2balls;
		bool projBGP;
			
		// save data
		bool saveEpsmaps;
		bool saveIdebmap;
		bool saveBgps;
		bool saveStatusMap;
		bool saveCavities;

		// global parameters
		string operativeMode;		
		int numthd;	
		bool printAvailSurf;	
		int currentSeed;

		// pocket detection
		bool cavAndPockets;
		bool linkPockets;	
		double pocketRadiusBig;		
		double pocketRadiusSmall;		
		double pocketRadiusLink;	
		bool debug;
};

extern Configuration conf;

#endif