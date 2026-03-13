//---------------------------------------------------------
/**    @file	SkinSurface.h
*     @brief	SkinSurface.h is the header for CLASS
*               SkinSurface.cpp								*/
//---------------------------------------------------------

#ifndef SkinSurface_h
#define SkinSurface_h

#include "Surface.h"

// #define DEBUG_SKIN

#define DEFAULT_S 0.45 // default shrink factor


/** mixed cell data structure*/
class MixedCell {
public:
	// double quadric[10];
	double *quadric;
	int surface_type;
	vector<double*> mc_points; // these are only pointers, no explicit denstructor

	virtual ~MixedCell()
	{
		if (quadric != NULL)
			deleteVector<double>(quadric);
	}
};


// voronoi facet cell or del edge
class Del1Cell : public MixedCell
{
public:
	vector<double*> planes; // lateral planes

	virtual ~Del1Cell()
	{
		for (unsigned int i=0; i<planes.size(); i++)
			deleteVector<double>(planes[i]);
	}
	// #IMPROVE one could build pointers from here to the planes computed in Del2Cell or viceversa
	double upper[4]; // upper plane
	double lower[4]; // upper plane

	int ids[2]; // the two atoms generating the del edge
};


// reduced voronoi cell
class Del0Cell : public MixedCell
{
public:
	vector<double*> planes; // planes pointer that points to upper/lower of Del1Cell

	int id; // atom index

	virtual ~Del0Cell()
	{
		// do nothing, planes are not managed by this class
	}
};


// del facet cell
class Del2Cell : public MixedCell
{
public:
	double planes[3][4]; // lateral planes
	double upper[4];
	double lower[4];

	int ids[3]; // atom indices

	virtual ~Del2Cell()
	{
		// do nothing, no pointers
	}
};

// reduced thethraedron
class Del3Cell: public MixedCell
{
public:

	// planes of the reduced tethraedron
	double planes[4][4];
	// double planes_or[4][4];
	// double points[4][3];
	double reduced[4][3]; // reduced points in the same order of atom indices
	double vor[3]; // Voronoi center

	int ids[4]; // atom indices

	virtual ~Del3Cell()
	{
		// do nothing, no pointers
	}
};


/*
#if !defined(OPTIMIZE_GRIDS)
#define GRID_MIXED_CELL_MAP_2D(i,j,l,NA,NB) gridMixedCellMap2D[((j)*(NA) + (i))*MAX_MIXED_CELLS_2D + (l)]

#define GRID_MIXED_CELL_MAP(i,j,k,l,NX,NY,NZ) gridMixedCellMap[((k)*((NY)*(NX)) + (j)*(NX) + (i))*MAX_MIXED_CELLS + (l)]
#endif
 */


#ifdef ENABLE_CGAL
//////////////////////// CGAL ///////////////////////////////////////////////////
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
// AV
// #include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// #include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <cassert>
#include <vector>
#include <fstream>
// use to translate in Pov-Ray
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

#ifdef CGAL_USE_BASIC_VIEWER
	#include <CGAL/draw_triangulation_3.h>
#endif
	
////////////////////////////////////////////////////////////////////////////////
#endif

#define DELAUNAY_POINT_CELL		0
#define DELAUNAY_EDGE_CELL		1
#define DELAUNAY_FACET_CELL		2
#define DELAUNAY_TETRA_CELL		3
#define SKIP_CELL				4


/** @brief This class builds and converts to a DelPhi suitable representation the Skin Surface.
All the gathered info is analytically computed both the intersections and the projections. In order
to get an accurate result for the projection routine, as root finding algorithm is used the method of the companion
matrix. The Skin surface was defined in: <i> "H. Edelsbrunner. Deformable smooth surface design. Discrete Comput. Geom., 21:87-115, 1999." </i>

@author Sergio Decherchi 
@date 30/10/2012
*/
class SkinSurface: public Surface
{

#ifdef ENABLE_CGAL
private:

	//AV
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

	typedef CGAL::Triangulation_cell_base_with_info_3< Del3Cell*,K> Cell_with_info;
	typedef CGAL::Regular_triangulation_cell_base_3<K, Cell_with_info> RT_Cell_with_info;

	#if !defined(NO_CGAL_PATCHING)
	typedef CGAL::Triangulation_vertex_base_3<K> Vertex;
	typedef CGAL::Triangulation_data_structure_3<
	CGAL::Regular_triangulation_vertex_base_3<K>,
	RT_Cell_with_info,
	CGAL::Parallel_tag>											tds;
	#else
	// nopatch
	typedef CGAL::Triangulation_vertex_base_with_info_3<int,K> _Vertex_with_info;
	typedef CGAL::Regular_triangulation_vertex_base_3<K,_Vertex_with_info>  Vertex_with_info;

	typedef CGAL::Triangulation_data_structure_3<
	Vertex_with_info,
	RT_Cell_with_info,
	CGAL::Parallel_tag>											tds;
	#endif

	typedef CGAL::Regular_triangulation_3<K, tds>				Rt;
	typedef Rt::Weighted_point									Weighted_point;
	typedef Rt::Vertex_handle									Vertex_handle;

	typedef Rt::Bare_point										Point3;
	typedef Rt::Vertex_iterator									Vertex_iterator;
	typedef Rt::Finite_vertices_iterator						Finite_Vertex_Iterator;
	typedef Rt::Finite_cells_iterator							Finite_Cells_Iterator;
	typedef tds::Cell_circulator								Cell_circulator;
	typedef tds::Cell_handle									Cell_handle;

	typedef Rt::Finite_edges_iterator							Finite_Edges_Iterator;
	typedef Rt::Finite_facets_iterator							Finite_Facets_Iterator;
	typedef Rt::Facet											Facet;

	typedef K::FT												Weight;


	/*
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Regular_triangulation_filtered_traits_3<K>    Traits;
	typedef CGAL::Triangulation_cell_base_with_info_3< Del3Cell*,Traits> Cell_with_info;
	typedef CGAL::Regular_triangulation_cell_base_3<Traits,Cell_with_info> RT_Cell_with_info;
	typedef CGAL::Triangulation_vertex_base_3<Traits> Vertex;
	typedef CGAL::Triangulation_data_structure_3<Vertex,RT_Cell_with_info> tds;
	typedef Traits::RT											Weight;
	typedef Traits::Bare_point									Point3;
	typedef Traits::Weighted_point								Weighted_point;
	typedef CGAL::Regular_triangulation_3<Traits,tds>			Rt;
	typedef Rt::Vertex_iterator									Vertex_iterator;
	typedef Rt::Finite_vertices_iterator						Finite_Vertex_Iiterator;
	typedef Rt::Finite_cells_iterator							Finite_Cells_Iterator;
	typedef Rt::Finite_edges_iterator							Finite_Edges_Iterator;
	typedef Rt::Finite_facets_iterator							Finite_Facets_Iterator;
	typedef Rt::Vertex_handle									Vertex_handle;
	typedef tds::Cell_circulator								Cell_circulator;
	typedef tds::Cell_handle									Cell_handle;
	typedef Rt::Facet											Facet;
	*/


	// used to translate in Pov-Ray
	typedef CGAL::Polyhedron_3<K>								Polyhedron;
	typedef Polyhedron::Halfedge_around_facet_circulator		HF_circulator;

	//a functor computing the plane containing a triangular facet
	struct Plane_from_facet {
		Polyhedron::Plane_3 operator()(Polyhedron::Facet &f) {
		Polyhedron::Halfedge_handle h = f.halfedge();
		return Polyhedron::Plane_3( h->vertex()->point(),
									h->next()->vertex()->point(),
									h->opposite()->vertex()->point());
		}
	};


	#ifdef MULTITHREADED_SKIN_BUILDING
	/** Local per-thread data mostly focused on cell and atom data. */
	struct ThreadDataWrapper
	{
		vector<MixedCell*> mixedComplex;

		#if !defined(NEW_ATOM_PATCHES)
		Del0Cell **atomPatches;
		#endif

		#if !defined(NO_CGAL_PATCHING)
		vector<Weighted_point> l;
		#else
		// nopatch
		vector<std::pair<Weighted_point, int>> l;
		#endif

		int num_cells[4];
	};

	#endif // MULTITHREADED_SKIN_BUILDING

#endif // ENABLE_CGAL


private:

	int numMixedCells, numSkippedCells;
	/** number of threads during cells' building, if any. */
	int numThreadDataWrappers = 1;
	int numThreadDataWrappersPerTask[MAX_TASKS];

	/** number of mixed cells for type.
	type 0 is a Delaunay point + Voronoi cell
	type 1 is a Delaunay edge  + Voronoi facet
	type 2 is a Delaunay facet + Voronoi edge
	type 3 is a Delaunay cell  + Voronoi point*/
	int type[4];
	/** link each auxuliary grid box to the respective mixed cell list */
	/* #if !defined(OPTIMIZE_GRIDS)
	int *gridConnollyCellMap;
	int *gridConnollyCellMap2D;
	#else */
	#if !defined(MULTITHREADED_SKIN_BUILDING)
	vector<int> *gridMixedCellMap;
	vector<int> *gridMixedCellMap2D;
	#else
	vector<pair<int,int>> *gridMixedCellMap;
	vector<pair<int,int>> *gridMixedCellMap2D;
	#endif
	// #endif

	double scale;
	double side;
	double *x,*y,*z;
	double s;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	int ***ind;

	/** auxiliary grid sizes*/
	int64_t nx,ny,nz;
	// last row and column sizes of ind_2d matrix for ray casting acceleration
	int64_t last_rows_ind;
	int64_t last_cols_ind;
	/** auxiliary grid sizes for the 2d map. */
	int64_t nx_2d,ny_2d,nz_2d;
	double scale_2d;
	double side_2d;
	double xmin_2d,xmax_2d,ymin_2d,ymax_2d,zmin_2d,zmax_2d;
	/* #if !defined(OPTIMIZE_GRIDS)
	unsigned int **ind_2d;
	#endif */

	unsigned int AUX_GRID_DIM_SKIN;
	unsigned int MAX_MIXED_CELLS;
	unsigned int AUX_GRID_DIM_SKIN_2D;
	unsigned int MAX_MIXED_CELLS_2D;

	/** for each mixed cell there is a set of planes that define it and a quadric equation. */
	vector<MixedCell*> mixedComplex;
	// MixedCell **mixedComplex;

	/** reduced tet cells that contain usefull reduced points but that does not give a real patch*/
	// vector<Del3Cell*> pendingCells;

	#if defined(MULTITHREADED_SKIN_BUILDING)
	/** local, per-thread data needed in the building stage. */
	// ThreadDataWrapper *thread_data_wrapper;
	ThreadDataWrapper thread_data_wrapper[MAX_TASKS_TIMES_THREADS];
	#endif

	/** compute the skin surface using CGAL regular triangulation and compute all information
	needed by to ray-trace it. */
	bool buildSkin();
	/** map each mixed cell to the auxiliary grid. */
	bool buildAuxiliaryGrid();

	#ifdef ENABLE_CGAL
	/** CGAL build skin */
	bool buildSkinCGAL();

	/** This function builds Delaunay tetrahedra. */
	void BuildDelaunayTetraCells (Rt &rT, vector<MixedCell*> &mixedComplex, int num_cells[]);
	#if !defined(NEW_ATOM_PATCHES)
	/** This function builds Delaunay point cells. */
	void BuildDelaunayPointCells (Rt &rT, vector<MixedCell*> &mixedComplex, Del0Cell **atomPatches, int num_cells[]);
	/** This function builds Delaunay edge cells (prismatic cells w/ parallel polygonal faces on the top/bottom). */
	void BuildDelaunayEdgeCells (Rt &rT, vector<MixedCell*> &mixedComplex, Del0Cell **atomPatches,
								 int num_atoms, int num_cells[], bool *any_error);
	#else
	/** This function builds Delaunay point cells. */
	void BuildDelaunayPointCells (Rt &rT, vector<MixedCell*> &mixedComplex, int atomPatches[], int num_cells[]);
	/** This function builds Delaunay edge cells (prismatic cells w/ parallel polygonal faces on the top/bottom). */
	void BuildDelaunayEdgeCells (Rt &rT, vector<MixedCell*> &mixedComplex, int atomPatches[],
								 int num_atoms, int num_cells[], bool *any_error);
	#endif
	/** This function builds Delaunay facet cells (prismatic cells w/ parallel triangular faces on the top/bottom). */
	void BuildDelaunayFacetCells (Rt &rT, vector<MixedCell*> &mixedComplex,
								  int num_cells[], bool *any_error);

	#if defined(MULTITHREADED_SKIN_BUILDING)
	/** This function builds Delaunay tetrahedra. */
	void BuildDelaunayTetraCells (Rt &rT, ThreadDataWrapper *tdw, int num_cells[]);
	/** This function builds Delaunay point cells. */
	void BuildDelaunayPointCells (Rt &rT, ThreadDataWrapper *tdw, int num_cells[]);
	/** This function builds Delaunay edge cells (prismatic cells w/ parallel polygonal faces on the top/bottom). */
	void BuildDelaunayEdgeCells (Rt &rT, ThreadDataWrapper *tdw, int num_cells[], bool *any_error);
	/** This function builds Delaunay facet cells (prismatic cells w/ parallel triangular faces on the top/bottom). */
	void BuildDelaunayFacetCells (Rt &rT, ThreadDataWrapper *tdw, int num_cells[], bool *any_error);
	void BuildPatches (ThreadDataWrapper *tdw, double minmax[6], bool *any_error);
	#endif // MULTITHREADED_SKIN_BUILDING
	#endif // ENABLE_CGAL

	bool savePovRay;
	bool fastProjection;

public:

	/** Default constructor. */
	SkinSurface();
	/** set DelPhi environment. */
	SkinSurface(DelPhiShared *ds);
	/** set configuration and DelPhi environment. */
	SkinSurface(ConfigFile *cf,DelPhiShared *ds);

	//////////////////////// INTERFACE MANDATORY METHODS /////////////////////////////////
	virtual int getNumPatches(void);
	virtual bool isPatchBasedRayTracingSupported(void);

	/** Compute skin surface. Call it after load*/
	virtual bool build();
	/** Get a projection of a point on the surface. Return projection and normal*/
	virtual bool getProjection(double p[3],double *proj1,double *proj2,
							   double *proj3,double *normal1,double *normal2,double *normal3);
	/** Get all the intersections of a ray that goes from P1 to P2 over the surface.
	The intersections are returned with increasing distance order.
	the first double in the vector is the t parameter for the intersection of the parametric 
	line and the surface, the double pointer is the normal vector. During ray surface intersection
	the previously built auxiliary grid is used to speed up computations*/
	virtual void getRayIntersection(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id);
	#if defined(REPORT_FAILED_RAYS)
	virtual void printRayIntersection(double pa[3], double pb[3]);
	#endif
	virtual void getRayIntersectionX(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id);
	virtual void getRayIntersectionY(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id);
	virtual void getRayIntersectionZ(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id);
	virtual void getRayIntersection(MixedCell *mc, double pa[3], double ray_dir, int id, int *netIntersections, int thread_id);

	/** function for the constructor without arguments. */
	virtual void init();
	/** functions for the constructor with config file argument. */
	virtual void init(ConfigFile *cf);
	/**function for the denstructor. */
	virtual void clear();
	/** pre-process panel to accelerate ray-tracing. */
	virtual void preProcessPanel();
	virtual void postRayCasting();
	virtual bool preBoundaryProjection();

	#if !defined(SINGLE_PASS_RT)
	/** Function used in multi-RT mode only; every patch is pierced by rays from multiple directions (piercing
	is done solely if MINIMIZE_MEMORY is defined in globals.h, otherwise conservative, approximated buffer
	offsets are estimated by approximated patch cells' data. The function is run by multiple threads */
	virtual void getPatchPreIntersectionData (int64_t nxyz[], int panels[2], int thread_id, int potentialIntersections[]);
	#endif
	/** Core patch-based RT function: each patch is pierced through all the needed directions and the pertinent
	buffers are filled. This function is run by multiple threads. */
	virtual void getPatchIntersectionData (int64_t nxyz[], int panels[2], int thread_id, int *netIntersections);
	// Patch normals at intersections can be computed here a posteriori; this can avoid normal data leakage when cleaning.
	// virtual void getPatchNormalsAtIntersections (int64_t nxyz[], int panels[2], int thread_id);
	/** Save it in a simple ASCII format (.skin). */
	virtual bool save(char *fileName);
	/**Load the surface from a file in .skin format. */
	virtual bool load(char *fileName);
	/** Print number of mixed cells and types. */
	virtual void printSummary();

	void setShrinking(double ss) 
	{ 
		double e = 0.05;
		if (ss <= (1.0-e) && ss >= (0.0+e))
		{ 
			s=ss; 
		}
		else
		{
			cout << endl << WARN << "Cannot set " << ss << ". s parameter is in (0+e,1-e], where e is "<< e << ".Setting "<< DEFAULT_S;
			s = DEFAULT_S;
		}
	}
	
	void setFastProjection(bool useFastProjection)
	{
		fastProjection = useFastProjection;
	}

	double getShrinking() 
	{
		return s;
	}

	void setSavePovRay(bool ss) 
	{ 
		savePovRay = ss;
	}
	
	bool getSavePovRay() 
	{
		return savePovRay;
	}

	/**for the 3d grid set the max grid size and the maximal number of patches inside a grid cube*/
	void setAuxGrid(unsigned int dim,unsigned int max)
	{		
		AUX_GRID_DIM_SKIN = dim;
		MAX_MIXED_CELLS = max;
	}

	/**for the 2d grid set the max grid size and the maximal number of patches inside a grid cube.
	The grid cube itself does not exist, it is just a reference; the relevant quantity is MAX_MIXED_CELLS_2D
	that is the number of patches along the grid tube*/
	void setAuxGrid2D(unsigned int dim,unsigned int max)
	{
		AUX_GRID_DIM_SKIN_2D = dim;
		MAX_MIXED_CELLS_2D = (max*dim);
	}
	
	virtual ~SkinSurface();

private:

	#ifdef POINT_METHOD_2
	double halfArea (double x1, double y1, double x2, double y2, double x3, double y3);
	#endif
	/** origin point, direction vector, quadric matrix, intersection parameter values. False
	is returned if no intersection is found. */
	// bool rayQuadricIntersection(double*,double*,double *cache);
	bool rayQuadricIntersection(double*,double*,double*,double *cache);
	bool rayQuadricIntersectionX(double*,double,double*,double*,double *cache);
	bool rayQuadricIntersectionY(double*,double,double*,double*,double *cache);
	bool rayQuadricIntersectionZ(double*,double,double*,double*,double *cache);
	bool rayQuadricIntersectionX(double*,double*,double*);
	bool rayQuadricIntersectionY(double*,double*,double*);
	bool rayQuadricIntersectionZ(double*,double*,double*);
	/** gives true if the point is inside the list of planes*/
	bool isFeasible(MixedCell *mc,double *point);
	/** project a point to a quadric surface defined by Q*/
	void projectToQuadric(double *y,double *Q,int type,double *proj,double *norm,double &dist);
	/** save quadric patch */
	#if !defined(NO_CGAL_PATCHING)
	void saveSkinPatch(ofstream &of,MixedCell *mc,int index,vector<Weighted_point> &la);
	#else
	// nopatch
	void saveSkinPatch(ofstream &of,MixedCell* mc,int index,vector<std::pair<Weighted_point,int>> &la);
	#endif
};


static class SkinSurfaceRegister{
	static Surface *createSurface(ConfigFile *conf,DelPhiShared *ds)
	{ 
		return new SkinSurface(conf,ds);
	} 
	public: 
		SkinSurfaceRegister() 
		{ 
			surfaceFactory().add("skin",createSurface); 
		} 
} SkinSurfaceRegisterObject;

// static SurfaceRecorder<SkinSurface> skinRecorder("skin");

#endif
