
//---------------------------------------------------------
/**    @file		MeshSurface.h
*     @brief	MeshSurface.h is the header for CLASS
*               MeshSurface.cpp								*/
//---------------------------------------------------------

#ifndef MeshSurface_h
#define MeshSurface_h

#include "Surface.h"

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif

#include "./ply/ply.h"

/** tollerance for triangle intersection test and gaussian elimination 1e-5*/
#define EPS 1e-5

// enable or disable culling in intersect_triangle routine 
// #define TEST_CULL


/* #if !defined(OPTIMIZE_GRIDS)
#define GRID_TRIANGLE_MAP_2D(i,j,l,NA,NB) gridTriangleMap2D[ ((j)*(NA) + (j))*MAX_TRIANGLES_2D + (l) ]

#define GRID_TRIANGLE_MAP(i,j,k,l,NX,NY,NZ) gridTriangleMap[ ((k)*((NY)*(NX)) + (j)*(NX) + (i))*MAX_TRIANGLES + (l) ]
#endif */


/** @brief This class represents a converter from an arbitray triangulated mesh surface
to a DelPhi compatible representation. Vertex normals are computed by averaging among
surrounding plane normals. 

@author Sergio Decherchi 
@date 14/10/2011
*/
class MeshSurface: public Surface
{
private:

	/** number of hreads during cells' building, if any. */
	int numThreadDataWrappers = 1;

	unsigned int MAX_TRIANGLES;
	unsigned int AUX_GRID_DIM;
	unsigned int MAX_TRIANGLES_2D;
	unsigned int AUX_GRID_DIM_2D;

	int numVertexes;
	int numTriangles;
	int **faceMatrix;
	VERTEX_TYPE **vertMatrix;
	/* #if !defined(OPTIMIZE_GRIDS)
	// link each auxuliary grid box to the respective triangles list in a 3d grid
	int *gridTriangleMap;
	// link each auxuliary grid box to the respective triangles list in a 2d grid
	int *gridTriangleMap2D;
	#else */
	vector<int> *gridTriangleMap;
	vector<int> *gridTriangleMap2D;
	// #endif
	double scale;
	double side;
	double *x,*y,*z;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	/* #if !defined(OPTIMIZE_GRIDS)
	int ***ind;
	*/
	double **planes;

	/** auxiliary grid sizes*/
	int64_t nx,ny,nz;
	// last row and column sizes of ind_2d matrix for ray casting acceleration
	int64_t last_rows_ind;
	int64_t last_cols_ind;
	/** auxiliary grid sizes for the 2d map*/
	int64_t nx_2d,ny_2d,nz_2d;
	double scale_2d;
	double side_2d;
	double xmin_2d,xmax_2d,ymin_2d,ymax_2d,zmin_2d,zmax_2d;
	/* #if !defined(OPTIMIZE_GRIDS)
	unsigned int **ind_2d;
	#endif */

	/** for each vertex the list of connected triangles */
	vector<int> **vertexTrianglesList;
	/** vertexes normals */
	VERTEX_TYPE **vertNormals;
	/** if enabled vertex normals are computed. This is true for a usual mesh; false
	for MSMS files that provide analytical vertex normals.*/
	//bool computeNormals;
	/** pre process ray casting panel*/
	virtual void preProcessPanel(void);
	/** pre process triangle normals*/
	void preProcessTriangles(void);
	/** build 3d auxiliary grid*/
	bool buildAuxiliaryGrid(void);
	/** intersect a ray into a triangle*/
	int intersect_triangle(double orig[3], double dir[3],double vert0[3], double vert1[3], double vert2[3],double *t, double *u, double *v);
	int intersect_triangle_X(double orig[3], double dir_x, double vert0[3], double vert1[3], double vert2[3], double *t);
	int intersect_triangle_Y(double orig[3], double dir_y, double vert0[3], double vert1[3], double vert2[3], double *t);
	int intersect_triangle_Z(double orig[3], double dir_z, double vert0[3], double vert1[3], double vert2[3], double *t);
	/** project a point in R3 to a triangle. If the projected point is outside
	the triangle the nearest point-segment projection is returned where the segment is 
	the nearest edge of the triangle to the point. The distance and the normal are returned.
	The plane ID (triangle number) and the vertex indexes must be passed.*/
	bool point2triangle(double P[3],double A[3], double B[3], double C[3],double w[4],
						double *proj,double *dist,double *normal,int planeID);
	/** project a point in R3 to a plane*/
	void point2plane(double p[3], double w[4],double *dist, double proj[3]);
	/* flag = inTriangle(P,A,B,C)
	says if the point p is in triangle or not
	flag is +1 if the vertex is in and -1 if it is out
	given the vertices A,B,C*/
	bool inTriangle(double P[3], double A[3], double B[3], double C[3]);
	/** check if a duplicated triangle or vertex is present*/
	bool checkDuplicates(void);
	/** load a mesh in off format*/
	bool loadOFF(char *fileName);
	/** load a mesh in ply format*/
	bool loadPLY(char *fileName);
	
	#if !defined(SINGLE_PASS_RT)
	/** Function used in multi-RT mode only; every patch is pierced by rays from multiple directions (piercing
	is done solely if MINIMIZE_MEMORY is defined in globals.h, otherwise conservative, approximated buffer
	offsets are estimated by approximated patch cells' data. The function is run by multiple threads */
	virtual void getPatchPreIntersectionData (int64_t nxyz[], int panels[2], int thread_id, int potentialIntersections[]);
	#endif
	/** Core patch-based RT function: each patch is pierced through all the needed directions and the pertinent
	buffers are filled. This function is run by multiple threads. */
	virtual void getPatchIntersectionData (int64_t nxyz[], int panels[2], int thread_id, int *netIntersections);

public:
	/** Default constructor*/
	MeshSurface();
	/** set DelPhi environment*/
	MeshSurface(DelPhiShared *ds);
	/** set configuration and DelPhi environment*/
	MeshSurface(ConfigFile *cf, DelPhiShared *ds);

	//////////////////////// INTERFACE MANDATORY METHODS /////////////////////////////////
	virtual int getNumPatches(void);
	virtual bool isPatchBasedRayTracingSupported(void);

	/** Nothing to do in this case*/
	virtual bool build(void)
	{
		return true;
	}

	virtual void postRayCasting(void);
	virtual bool preBoundaryProjection(void)
	{
		// 3d grid is necessary only for boundary grid points projection
		if (projBGP)
			return buildAuxiliaryGrid();
		return false;
	}

	/** Save in off format*/
	virtual bool save(char *fileName);
	/**Load the surface from an OFF/PLY file*/
	virtual bool load(char *fileName);
	/**Load MSMS mesh surface*/
	bool loadMSMS(char *fileName,int numFiles=1);
	/** Print number of vertices and faces*/
	virtual void printSummary(void);
	/** Get a projection of a point on the surface. Return projection and normal*/
	virtual bool getProjection(double p[3],double *proj1,double *proj2,
							   double *proj3,double *normal1,double *normal2,double *normal3);
	/** Get all the intersections of a ray that goes from P1 to P2 over the surface.
	The intersections are returned with increasing distance order.
	the first double in the vector is the t parameter for the intersection of the parametric 
	line and the surface, the double pointer is the normal vector. During ray triangle intersection
	the previously built auxiliary grid is used to speed up computations*/
	virtual void getRayIntersection(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals,int thread_id);
	#if defined(REPORT_FAILED_RAYS)
	virtual void printRayIntersection(double pa[3], double pb[3]);
	#endif
	/** function for the constructor without arguments*/
	virtual void init(void);
	/** functions for the constructor with config file argument*/
	virtual void init(ConfigFile *cf);
	/**function for the denstructor*/
	virtual void clear(void);
	/////////////////////////////////////////////////////////////
	/**for the 3d grid set the max grid size and the maximal number of patches inside a grid cube*/
	void setAuxGrid(unsigned int dim,unsigned int max)
	{
		AUX_GRID_DIM = dim;
		MAX_TRIANGLES = max;
	}

	/**for the 2d grid set the max grid size and the maximal number of patches inside a grid cube.
	The grid cube itself does not exist just a reference, indeed the real quantity is MAX_TRIANGLES_2D
	that is the number of patches along the grid tube*/
	void setAuxGrid2D(unsigned int dim,unsigned int max)
	{
		AUX_GRID_DIM_2D = dim;
		MAX_TRIANGLES_2D = (max*dim);
	}

	virtual ~MeshSurface();
};

// expand it explicitly because Swig is not able to expand it
static class MeshSurfaceRegister{
	static Surface *createSurface(ConfigFile *conf,DelPhiShared *ds)
	{ 
		return new MeshSurface(conf,ds);
	} 
	public: 
		MeshSurfaceRegister() 
		{ 
			surfaceFactory().add("mesh",createSurface); 
		} 
} MeshSurfaceRegisterObject;


//static SurfaceRecorder<MeshSurface> meshRecorder("mesh");

#endif
