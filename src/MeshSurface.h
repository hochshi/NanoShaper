
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
//#define TEST_CULL 

#define GRID_TRIANGLEMAP_2D(i,j,k,NA,NB) gridTriangleMap2D[(k)+(MAX_TRIANGLES_2D-1)*((j)+(NB)*(i))]
#define GRIDTRIANGLEMAP(i,j,k,l,NX,NY,NZ) gridTriangleMap[(l)+(MAX_TRIANGLES-1)*((k)+(NZ)*((j)+(NY)*(i)))]

/** @brief This class represents a converter from an arbitray triangulated mesh surface
to a DelPhi compatible representation. Vertex normals are computed by averaging among
surrounding plane normals. 

@author Sergio Decherchi 
@date 14/10/2011
*/
class MeshSurface: public Surface
{
protected:
	unsigned int MAX_TRIANGLES;
	unsigned int AUX_GRID_DIM;
	unsigned int MAX_TRIANGLES_2D;
	unsigned int AUX_GRID_DIM_2D;
	int numVertexes;
	int numTriangles;
	int** faceMatrix;
	double** vertMatrix;	
	/** link each auxuliary grid box to the respective triangles list in a 3d grid*/
	int* gridTriangleMap;
	/** link each auxuliary grid box to the respective triangles list in a 2d grid*/
	int* gridTriangleMap2D;
	/** auxiliary grid sizes*/
	int nx,ny,nz;
	double scale;
	double side;
	double*x,*y,*z;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	unsigned short*** ind;
	double** planes;

	/** auxiliary grid sizes for the 2d map*/
	int nx_2d,ny_2d,nz_2d;
	double scale_2d;
	double side_2d;
	double xmin_2d,xmax_2d,ymin_2d,ymax_2d,zmin_2d,zmax_2d;
	unsigned int** ind_2d;

	/** for each vertex the list of connected triangles */
	vector<int>** vertexTrianglesList;
	/** vertexes normals */
	double** vertNormals;
	/** if enabled vertex normals are computed. This is true for a usual mesh; false
	for MSMS files that provide analytical vertex normals.*/
	bool computeNormals;
	/** pre process ray casting panel*/
	virtual void preProcessPanel();
	/** pre process triangle normals*/
	void preProcessTriangles();
	/** build 3d auxiliary grid*/
	bool buildAuxiliaryGrid();
	/** intersect a ray into a triangle*/
	int intersect_triangle(double orig[3], double dir[3],double vert0[3], double vert1[3], double vert2[3],double *t, double *u, double *v);
	/** project a point in R3 to a triangle. If the projected point is outside
	the triangle the nearest point-segment projection is returned where the segment is 
	the nearest edge of the triangle to the point. The distance and the normal are returned.
	The plane ID (triangle number) and the vertex indexes must be passed.*/
	bool point2triangle(double P[3],double A[3], double B[3], double C[3],double w[4],
		double* proj,double *dist,double* normal,int planeID,int va,int vb,int vc);
	/** project a point in R3 to a plane*/
	void point2plane(double p[3], double w[4],double* dist, double proj[3]);
	/* flag = inTriangle(P,A,B,C)
	says if the point p is in triangle or not
	flag is +1 if the vertex is in and -1 if it is out
	given the vertices A,B,C*/
	bool inTriangle(double P[3], double A[3], double B[3], double C[3]);
	/** check if a duplicated triangle or vertex is present*/
	bool checkDuplicates();
	/** load a mesh in off format*/
	bool loadOFF(char* fileName);
	/** load a mesh in ply format*/
	bool loadPLY(char* fileName);

public:
	/** Default constructor*/
	MeshSurface();
	/** set DelPhi environment*/
	MeshSurface(DelPhiShared* ds);
	/** set configuration and DelPhi environment*/
	MeshSurface(ConfigFile* cf,DelPhiShared* ds);			

	//////////////////////// INTERFACE MANDATORY METHODS /////////////////////////////////
	/** Nothing to do in this case*/
	virtual bool build()
	{
		return true;
	}

	virtual void postRayCasting()
	{
		// remove 2d grid for ray casting 
		if (gridTriangleMap2D!=NULL)
			deleteVector<int>(gridTriangleMap2D);	
		
		if (ind_2d!=NULL)
			deleteMatrix2D<unsigned int>(last_rows_ind,ind_2d);

		// no more needed
		//gridTriangleMap2D = NULL;
		//ind_2d = NULL;
	}

	virtual bool preBoundaryProjection()
	{
		// 3d grid is necessary only for boundary grid points projection
		if (projBGP)
			return buildAuxiliaryGrid();
		return false;
	}

	/** Save in off format*/
	virtual bool save(char* fileName);
	/**Load the surface from an OFF/PLY file*/
	virtual bool load(char* fileName);
	/**Load MSMS mesh surface*/
	bool loadMSMS(char* fileName,int numFiles=1);
	/** Print number of vertices and faces*/
	virtual void printSummary();		
	/** Get a projection of a point on the surface. Return projection and normal*/
	virtual bool getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3);
	/** Get all the intersections of a ray that goes from P1 to P2 over the surface.
	The interesctions are returned with increasing distance order. 
	the first double in the vector is the t parameter for the intersection of the parametric 
	line and the surface, the double pointer is the normal vector. During ray triangle intersection
	the previously built auxiliary grid is used to speed up computations*/
	virtual void getRayIntersection(double p1[3],double p2[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals);
	/** function for the constructor without arguments*/
	virtual void init();
	/** functions for the constructor with config file argument*/
	virtual void init(ConfigFile* cf);
	/**function for the denstructor*/
	virtual void clear();
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
	static Surface* createSurface(ConfigFile* conf,DelPhiShared* ds) 
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