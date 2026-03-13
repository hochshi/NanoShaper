
//---------------------------------------------------------
/**    @file		Surface.h
*     @brief	Surface.h is the header for CLASS
*               Surface.cpp								*/
//---------------------------------------------------------

#ifndef Surface_h
#define Surface_h


#include "octree.h"
#include "globals.h"
#include "ConfigFile.h"
#include "SurfaceFactory.h"
#include <stack>


#if defined(USE_VIS_TOOLS)

#define GL_GLEXT_PROTOTYPES

#ifndef NO_OPENGL

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/freeglut.h>
#endif

#ifdef __cplusplus
}
#endif
#endif // NO_OPENGL

#endif // USE_VIS_TOOLS


#ifdef ENABLE_CGAL 
	//////////////////////// CGAL ///////////////////////////////////////////////////
	#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	#include <CGAL/Regular_triangulation_3.h>
	// #include <CGAL/Regular_triangulation_euclidean_traits_3.h>
    // <CGAL/Regular_triangulation_euclidean_traits_3.h
    #include <CGAL/config.h>
	// #include <CGAL/Regular_triangulation_filtered_traits_3.h>
	#include <CGAL/Triangulation_cell_base_with_info_3.h>
	#include <CGAL/Regular_triangulation_cell_base_3.h>
	#include <CGAL/Triangulation_vertex_base_with_info_3.h>
	#include <CGAL/Triangulation_data_structure_3.h>
	////////////////////////////////////////////////////////////////////////////////

	// vor point container
	class VorPoint
	{
	public:
		double vor[3]; // Voronoi center
		bool visited;

		VorPoint()
		{
			visited = false;
		}
	};
#endif


#ifdef COORD_NORM_PACKING
/** grid coordinates + ray dir. + patch coordinates + normal at ray vs patch intersection point */

#if !defined(COMPRESS_INTERSECTION_COORDS)

class coordNormPacket
{
public:

	coordNormPacket(const int x, const int y, const int z,
					VERTEX_TYPE *const v, VERTEX_TYPE *const n,
					const int d)
	{
		if (v != NULL)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}
		else
		{
			vec[0] = 0;
			vec[1] = 0;
			vec[2] = 0;
		}
		if (n != NULL)
		{
			nor[0] = n[0];
			nor[1] = n[1];
			nor[2] = n[2];
		}
		else
		{
			nor[0] = 0;
			nor[1] = 0;
			nor[2] = 0;
		}
		ix = x;
		iy = y;
		iz = z;
		dir = d;
	}
	
	// shallow copy constructor
	coordNormPacket(const coordNormPacket &cnp)
	{
		vec[0] = cnp.vec[0];
		vec[1] = cnp.vec[1];
		vec[2] = cnp.vec[2];
		nor[0] = cnp.nor[0];
		nor[1] = cnp.nor[1];
		nor[2] = cnp.nor[2];
		ix = cnp.ix;
		iy = cnp.iy;
		iz = cnp.iz;
		dir = cnp.dir;
	}
	
	VERTEX_TYPE vec[3], nor[3];
	int ix;
	int iy;
	int iz;
	int dir;
};

#else // COMPRESS_INTERSECTION_COORDS

class compressedCoordNormPacket
{
public:

	compressedCoordNormPacket(const uint64_t x, const uint64_t y, const uint64_t z,
							  VERTEX_TYPE *const v, VERTEX_TYPE *const n, const uint64_t dir)
	{
		vertex_coords = (dir << 60UL) | (z << 40UL) | (y << 20UL) | x;
		if (v != NULL)
		{
			vec[0] = v[0];
			vec[1] = v[1];
			vec[2] = v[2];
		}
		else
		{
			vec[0] = 0;
			vec[1] = 0;
			vec[2] = 0;
		}
		if (n != NULL)
		{
			nor[0] = n[0];
			nor[1] = n[1];
			nor[2] = n[2];
		}
		else
		{
			nor[0] = 0;
			nor[1] = 0;
			nor[2] = 0;
		}
	}

	void getCompressedCoords(int *x, int *y, int *z, int *dir)
	{
		*x  = (int)((this->vertex_coords        ) & ((1UL << 20UL) - 1UL));
		*y  = (int)((this->vertex_coords >> 20UL) & ((1UL << 20UL) - 1UL));
		*z  = (int)((this->vertex_coords >> 40UL) & ((1UL << 20UL) - 1UL));
		*dir = (int)(this->vertex_coords >> 60UL);
	}

	void getCompressedCoords(const uint64_t coords, int *x, int *y, int *z, int *dir)
	{
		*x  = (int)((coords        ) & ((1UL << 20UL) - 1UL));
		*y  = (int)((coords >> 20UL) & ((1UL << 20UL) - 1UL));
		*z  = (int)((coords >> 40UL) & ((1UL << 20UL) - 1UL));
		*dir = (int)(coords >> 60UL);
	}

	VERTEX_TYPE vec[3], nor[3];
	uint64_t vertex_coords;
};
#endif // COMPRESS_INTERSECTION_COORDS

#endif // COORD_NORM_PACKING


/** grid coordinates + ray dir. + patch coordinates or normal at ray vs patch intersection point */
class coordVec
{
public:

	coordVec(const int x, const int y, const int z, VERTEX_TYPE *const v, const int d)
	{
		vec = v;
		ix = x;
		iy = y;
		iz = z;
		dir = d;
	}

	// shallow copy constructor
	coordVec(const coordVec &cv)
	{
		vec = cv.vec;
		ix = cv.ix;
		iy = cv.iy;
		iz = cv.iz;
		dir = cv.dir;
	}

	VERTEX_TYPE *vec;
	int ix;
	int iy;
	int iz;
	int dir;
};


class packet
{
public:

	packet()
	{
		first = nullptr;
		#if !defined(COORD_NORM_PACKING)
		second = nullptr;
		#endif
	}

	packet(const packet &p)
	{
		first = p.first;
		#if !defined(COORD_NORM_PACKING)
		second = p.second;
		#endif
	}

	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *first;
	vector<coordVec> *second;
	#else
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *first;
	#else
	vector<compressedCoordNormPacket> *first;
	#endif
	#endif
};


#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif

#include "tools.h"
#include "DelphiShared.h"

// ids for rays directions
#define X_DIR 0
#define Y_DIR 1
#define Z_DIR 2

// ids for panels
#define PANEL_YZ 0
#define PANEL_XY 1
#define PANEL_XZ 2

#define DEFAULT_VOLUME 11.4
// tolerance on two intersections. If less than EPS_INT then the intersection is the same. def 1e-8
#define EPS_INT 1e-8
// #define MAX_ATOMS_MULTI_GRID 100
// #define DEBUG_SURFACE

#define DELTA 1e-10
#define RAND_DISPLACEMENT 1e-4

// boundary grid point codes
#define INTERNAL_BGP 0
#define EXTERNAL_BGP 1

// molecular surface
#define MOLECULAR_SURFACE 0
// analytical object
#define OBJECT_SURFACE 3
// unknown surface, e.g., mesh
#define GENERIC_SURFACE 2
// 'sum' of two different surface types
#define HYBRID_SURFACE 4

// file format
#define OFF 0
#define OFF_A 1
#define OFF_N 2
#define OFF_N_A 3
#define MSMS_NO_A 4
#define MSMS 5
#define PLY 6

// now using inline templates
// #define GRID_MULTI_MAP(i,j,k,l,NX,NY,NZ) gridMultiMap[( (k)*((NY)*(NX)) + (j)*(NX) + (i) )*MAX_ATOMS_MULTI_GRID + (l)]

#ifdef ENABLE_BOOST_THREADS

#define THREAD_SAFE_SCOPE(SSS) \
	{ \
		boost::mutex::scoped_lock scopedLock(mutex); \
		SSS \
		cout.flush(); \
	}

#else

#define THREAD_SAFE_SCOPE(SSS) \
	{ \
		SSS \
		cout.flush(); \
	}
#endif


#if defined(USE_VIS_TOOLS)
#ifndef NO_OPENGL

#define ROTATE          2
#define ZOOM            3
#define DISTANCE        4
#define QUIT            7

#endif
#endif // USE_VIS_TOOLS





/** @brief Surface class is the general interface that a surface class should have
to be plugged inside DelPhi. Some functions implementations are mandatory
such as load, save, build, etc... Note that the surface is not necessarly a molecular surface.
Build function computes an internal representation of the surface; getSurf translates
that representation in the DelPhi compatible representation.
Note that in order to put a new surface in DelPhi a surface must provide epsmap,
idebmap, computations of the surface area inside a grid cube, identification and projections
of boundary grid points and their surface normalSkinSurfaces; these computations must be done in getSurf while the surface construction
must be performed in build. \n \n

The global variables in the DelPhiShared object that must be filled are: \n
idebmap -> salt (in/out) map \n
epsmap  -> dielectric map \n
ibgp    -> indexes of the boundary grid points \n
scspos  -> coordinates of the projected boundary grid points \n
scsnor  -> coordinates of the outside pointing normal vector over the surface in correspondence of boundary grid points \n
scsarea -> value of the surface area in each cube of the grid where there is a boundary grid point (optional)\n

If a new surface is defined and it only provides the routines of ray intersection and projection
then most information can be still built; in this case one should not overload the getSurf method
that in its base implementation is able to use projections/intersections to build almost all info needed by DelPhi to work; however this mode
of operating has restrictions in fact only two dielectrics can be used (in/out) and the Stern layer (idebmap)
is only in/out; this mode is used in the MeshSurface class where an arbitrary shape is loaded or by the SkinSurface module. Note
that there is no need for the surface (or surface patches) to answer in/out queries; indeed in/out is worked
out by counting the number of times the ray intersects the surface starting from the outside.

To get a full implementation one has to derive the Surface class and re-implement the interface method
load,save,build and getSurf (or ray/projections routine) thus providing all the necessary info to DelPhi solver.


@author Sergio Decherchi
@date 29/06/2013
*/
class Surface
{
#ifdef ENABLE_CGAL 
private:

	//AV
	typedef CGAL::Exact_predicates_inexact_constructions_kernel _K;
 	// Regular T3


	typedef CGAL::Triangulation_cell_base_with_info_3<VorPoint*, _K> _Cell_with_info;
	typedef CGAL::Regular_triangulation_cell_base_3<_K, _Cell_with_info> _RT_Cell_with_info;

	#if !defined(NO_CGAL_PATCHING)
	typedef CGAL::Triangulation_vertex_base_3<_K> _Vertex;

	typedef CGAL::Triangulation_data_structure_3<
	CGAL::Regular_triangulation_vertex_base_3<_K>,
	_RT_Cell_with_info,
	CGAL::Parallel_tag>											_tds;
	#else
	// nopatch
	typedef CGAL::Triangulation_vertex_base_with_info_3<int,_K> _Vertex_with_info;
	typedef CGAL::Regular_triangulation_vertex_base_3<_K, _Vertex_with_info> _RT_Vertex_with_info;

	typedef CGAL::Triangulation_data_structure_3<_RT_Vertex_with_info,
	_RT_Cell_with_info,
	CGAL::Parallel_tag>											  _tds;
	#endif

	typedef CGAL::Regular_triangulation_3<_K, _tds>               _Rt;
	typedef _Rt::Weighted_point                                   _Weighted_point;
	typedef _Rt::Vertex_handle                                    _Vertex_handle;

	typedef _Rt::Bare_point                                       _Point3;
	typedef _Rt::Vertex_iterator                                  _Vertex_iterator;
	typedef _Rt::Finite_vertices_iterator                         _Finite_Vertex_Iterator;
	typedef _Rt::Finite_cells_iterator                            _Finite_Cells_Iterator;
	typedef _tds::Cell_circulator                                 _Cell_circulator;
	typedef _tds::Cell_handle                                     _Cell_handle;

	typedef _Rt::Finite_edges_iterator                            _Finite_Edges_Iterator;
	typedef _Rt::Finite_facets_iterator                           _Finite_Facets_Iterator;
	typedef _Rt::Facet                                            _Facet;

	typedef _K::FT                                                _Weight;

	/*
	typedef CGAL::Triangulation_vertex_base_3<_K> _Vertex;

	typedef CGAL::Triangulation_data_structure_3<
	CGAL::Regular_triangulation_vertex_base_3<_K>,
	_RT_Cell_with_info,
	CGAL::Parallel_tag>											_tds;
	typedef CGAL::Regular_triangulation_3<_K, _tds>				_Rt;
	typedef _Rt::Weighted_point									_Weighted_point;
	typedef _Rt::Vertex_handle									_Vertex_handle;

	typedef _Rt::Bare_point										_Point3;
	typedef _Rt::Vertex_iterator								_Vertex_iterator;
	typedef _Rt::Finite_vertices_iterator						_Finite_Vertex_Iterator;
	typedef _Rt::Finite_cells_iterator							_Finite_Cells_Iterator;
	typedef _tds::Cell_circulator								 _Cell_circulator;
	typedef _tds::Cell_handle									_Cell_handle;

	typedef _Rt::Finite_edges_iterator							_Finite_Edges_Iterator;
	typedef _Rt::Finite_facets_iterator							_Finite_Facets_Iterator;
	typedef _Rt::Facet											_Facet;

	typedef _K::FT												_Weight;
	*/

	/*
	typedef CGAL::Triangulation_cell_base_with_info_3< VorPoint*,_K> _Cell_with_info;
	typedef CGAL::Regular_triangulation_cell_base_3<_K,_Cell_with_info> _RT_Cell_with_info;
	typedef CGAL::Triangulation_vertex_base_3<_K> _Vertex;
	// typedef CGAL::Triangulation_data_structure_3<_Vertex,_RT_Cell_with_info> _tds;

	typedef CGAL::Triangulation_data_structure_3<
            CGAL::Regular_triangulation_vertex_base_3<K>,
            CGAL::Regular_triangulation_cell_base_3<K>,
            CGAL::Parallel_tag>									_tds;


	// typedef CGAL::Triangulation_data_structure_3<
	// CGAL::Regular_triangulation_vertex_base_3<_K>,
	// CGAL::Regular_triangulation_cell_base_3<_K>,
	// CGAL::Parallel_tag>										_tds;
	typedef CGAL::Regular_triangulation_3<_K, _tds>		_Rt;
	typedef _Rt::Bare_point										_Bare_point;
	typedef _Rt::Weighted_point									_Weighted_point;
	typedef _Rt::Vertex_handle									_Vertex_handle;

	typedef _Rt::Bare_point										_Point3;
	typedef _Rt::Vertex_iterator								_Vertex_iterator;
	typedef _Rt::Finite_vertices_iterator					_Finite_Vertex_Iterator;
	typedef _Rt::Finite_cells_iterator						_Finite_Cells_Iterator;
	typedef _tds::Cell_circulator								_Cell_circulator;
	typedef _tds::Cell_handle									_Cell_handle;

	typedef _K::FT													_Weight;
	*/
	/*
	typedef CGAL::Exact_predicates_inexact_constructions_kernel		_K;
	typedef CGAL::Regular_triangulation_filtered_traits_3<_K>		_Traits;
        typedef CGAL::Regular_triangulation_euclidean_traits_3<_K>	_Traits;
	typedef CGAL::Triangulation_cell_base_with_info_3< VorPoint*,_Traits> _Cell_with_info;
	typedef CGAL::Regular_triangulation_cell_base_3<_Traits,_Cell_with_info> _RT_Cell_with_info;
	typedef CGAL::Triangulation_vertex_base_3<_Traits> _Vertex;
	typedef CGAL::Triangulation_data_structure_3<_Vertex,_RT_Cell_with_info> _tds;

	typedef _Traits::RT												_Weight;
	typedef _Traits::Bare_point									_Point3;
	typedef _Traits::Weighted_point								_Weighted_point;

	typedef CGAL::Regular_triangulation_3<_Traits,_tds>	_Rt;
	typedef _Rt::Vertex_iterator									_Vertex_iterator;
	typedef _Rt::Finite_vertices_iterator						_Finite_Vertex_Iterator;
	typedef _Rt::Finite_cells_iterator							_Finite_Cells_Iterator;
	typedef _Rt::Vertex_handle										_Vertex_handle;
	typedef _tds::Cell_circulator									_Cell_circulator;
	typedef _tds::Cell_handle 										_Cell_handle;
	*/

#endif // ENABLE_CGAL


public:
	// current panel under analysis
	int panel;

protected:

	//////////////////////////////////// surf definition variables ///////////////////////
	/** says if the surface represents a molecule, an hybrid system or an object */
	int surfType;
	/** every class at Constructor time must say if it is Ray Casting based or not
	 If it is ray casting based than getSurf will derive the grid based on ray casting
	 if not RC based getSurf will assume that build method has already coloured the grid */
	bool isRCbased;
	/** every surface at constructor time must say if it will provide analytical normals during
	ray tracing or not. */
	bool providesAnalyticalNormals;
	////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////// flags/settings ////////////////////////////////////
	bool projBGP;
	double sternLayer;
	bool accurateTriangulation;
	bool doTriangulation;
	bool fillCavitiesFlag;
	bool computeNormals;
	/** if enabled loaded surface is checked for duplicated vertices */
	bool checkDuplicatedVertices;
	/** if enabled the part of the cavities (or the entire cavities) that are not shaped as a water molecule are
	removed. */
	bool wellShaped;
	/** Probe radius. Currently used in the cavity shape filter */
	double probe_radius;
	bool useLoadBalancing;
	bool vertexAtomsMapFlag;
	bool saveMSMS;
	bool savePLY;
	
	/** Two recent additions to allow choosing patch-based RT algorithm and analythical torus vs ray intersections. */
	bool patchBasedAlgorithm;
	bool analyticalTorusIntersectionAlgorithm;

	/** Collect face intersections or not, even without epsmap */
	bool collectFaceIntersections;

	/** Determine the halo layer thickness in the parallel buildup stage */
	double parallelBuildupHaloThickness;

	/** This allows to force the bulding phase to be serial */
	bool forceSerialBuild;

	/** Flag used for optimizing grids to reduce memory consumption, e.g. with bilevel
	hierarchical grids (instead of full flat uniform grids), or not */
	bool optimizeGrids;

	int maxNumAtoms;
	double domainShrinkage;

	//////////////////////////////////////////////////////////////////////////////////////

	// delphi environment
	DelPhiShared *delphi;

	double **threadPanelVolume;
	
	int *threadFailedRays, *threadTotalRays;
	int panelVolumeFlag[3][2];

	/** how big is the random initial displacement of atoms*/
	double randDisplacement;
	// last nx,ny,nz dimensions seen by Surface class
	int64_t last_nx,last_ny,last_nz;
	#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
	/** 3d matrix of intersections along x rays*/
	Octree<int> *intersectionsMatrixAlongX;
	/** 3d matrix of intersections along y rays*/
	Octree<int> *intersectionsMatrixAlongY;
	/** 3d matrix of intersections along z rays*/
	Octree<int> *intersectionsMatrixAlongZ;
	#if !defined(AVOID_NORMALS_MATRICES)
	/** 3d matrix of normals along x rays*/
	Octree<int> *normalsMatrixAlongX;
	/** 3d matrix of normals along y rays*/
	Octree<int> *normalsMatrixAlongY;
	/** 3d matrix of normals along z rays*/
	Octree<int> *normalsMatrixAlongZ;
	#endif
	#else // OPTIMIZE_INTERSECTIONS_MANAGEMENT
	int **bilevel_intersectionsMatrixAlongX;
	int **bilevel_intersectionsMatrixAlongY;
	int **bilevel_intersectionsMatrixAlongZ;
	#endif // OPTIMIZE_INTERSECTIONS_MANAGEMENT
	

	/** mark wich MC cubes contain triangles to allow fast reject in the second pass. */
	#if !defined(USE_COMPRESSED_GRIDS)
	bool *activeCubes;
	#endif
	// compressed version of activeCubes: in each 32-bit word there are 32 1-bit values
	unsigned int *compressed_activeCubes;
	// compressed buffer useful to skip in triangulationKernel()
	// empty mini grids with 4^3 cubes. This buffer is filled in getVertices().
	unsigned int *compressed_activeMacroCubes;

	#if !defined(USE_COMPRESSED_GRIDS)
	// bool *verticesInsidenessMap;
	bool ***verticesInsidenessMap;
	#endif
	// in each 32-bit word there are 32 1-bit values
	unsigned int *compressed_verticesInsidenessMap;

	double ***scalarField;
	
	// when a scalar field is available scalarField is used instead of verticesInsidenessMap for vertex interpolation
	bool isAvailableScalarField;
	
	////////////////////////////// triangulation data structures ////////////////////////////////
	/** vector of vertex indices for the traingulation obtained by triangulateSurface function. */
	vector<int> triList;
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	/** vector of triangle vertices. */
	vector<VERTEX_TYPE*> vertList;
	/** vector of normal to vertices. */
	vector<VERTEX_TYPE*> normalsList;
	#else
	/** vector of triangle vertices. */
	vector<VERTEX_TYPE> vertList;
	/** vector of normal to vertices. */
	vector<VERTEX_TYPE> normalsList;
	#endif

	/** Buffers employed to access data of vertices and normals with a fairly coherent manner */
	vector<VERTEX_TYPE> verticesBuffers[MAX_TASKS_TIMES_THREADS];
	vector<VERTEX_TYPE> normalsBuffers[MAX_TASKS_TIMES_THREADS];

	/////////////////////////////////////////////////////////////////////////////////////////////

	double delta_accurate_triangulation;
	
	/** type of bgp for each detected bgp **/
	int *bgp_type;
	
	/** grid multi-dielectric map */
	vector<int> *gridMultiMap;
	double gxmin,gymin,gzmin,gside,gscale;
	int64_t ggrid;
	
	/** pointer to a vector which holds the information about the load for each grid slice*/
	int *gridLoad;
	int totalLoad;

	int *vertexAtomsMap;

	/** DelPhi code for inside.*/
	int inside;

	#ifdef ENABLE_BOOST_THREADS
	boost::mutex mutex;
	#endif

	// added keyword to support saving of files 
	string rootFile;

	int numTasks = 1;

	/////////////////////////////////////////////////////////////////////////////////////////////////

	/** 3d cavity detection. */
	void floodFill(int ix,int iy,int iz,int idold,int idnew);
	
	/** Scan line version. */
	void floodFill2(int ix,int iy,int iz,int idold,int idnew);
	void floodFillWithBilevelStatusMap2(int ix,int iy,int iz,int idold,int idnew);
	
	/** Parallel version with scaline */
	void floodFill4(int ix,int iy,int iz,int idold,int idnew);
	
	/** Inner routine for scanline. */
	void floodFill3(	pair<pair<int,int>,int> ind,
						pair<int,int> z_limits, 
						pair<int,int> old_new, 
						pair<queue<pair<pair<int,int>,int>>*, queue<pair<pair<int,int>,int>>*> queues,
						queue<pair<pair<int,int>,int>> *in_queue);
	
	/** This supports inner floodfill inside the 1.4 surface (given status1) to
	understand if two cavities/pockets are status2-linked (check external communication
	between two points by looking at status2 map). Gives true if the cavities/pockets comunicate
	and uses as input two random indices of two cavities/pockets. It stops if a maximal number of
	moves have been reached. A 'move' means moving from one grid point to another one.
	In max moves is returned the number of moves done to get the target. */
	bool innerFloodFill(int *start,int *target,int *status1,int *status2,int &maxMoves,bool debug);
	bool innerFloodFill(int *start,int *target,int **bilevel_status1,int **bilevel_status2,int &maxMoves,bool debug);
	
	/** This gives true if the point is outside vdw surface*/
	bool vdwAccessible(double *p,int &nearest);
	
	/** Ray tracing routine employed to perform partial or full intersections used
	 together with boost threading routines. In order to get a 'robust' ray tracer a
	 particular strategy is adopted during ray-tracing. It can happen that, due to
	 numerical imprecisions of the ray-intersection routine or due to problems of the
	 surface (e.g. degenerate triangles, holes, non manifoldness in general)
	 the number of detected intersections is odd. In this case, the trace of the
	 previous ray is copied to the current one. \n
	 In some cases this strategy can also fill small holes in the surface.
	 
	 SD PB_NEW: added now the package to collect grid and normals intersections. */
	void intersectWithRayBasedAlgorithm(int thread_id,int nb,int start,int end,int iters_block,
										int jump,packet pack, packet gridPack=packet());


	/** This routine is very similar to intersectWithRayBasedAlgorithm(...) but does not compute
	the intersections: it stores vertices' coords and normals using intersection data collected during
	patch-based steps getPatch(Pre)IntersectionData(...), and updates grid data (e.g. the status map),
	if required. */
	void setVerticesAndGridsWithIntersectionData(int thread_id,int nb,int start,int end,int iters_block,
												 int jump,packet pack,packet gridPack=packet());

	/** This function assembles the cross-thread intersection data with the help of octrees
	(if OPTIMIZE_INTERSECTIONS_MANAGEMENT is not defined in globals.h) or bilevel grids. */
	void assembleVerticesList (packet pack, vector<VERTEX_TYPE*> *localVert, vector<VERTEX_TYPE*> *localNormals, int *localIndex);
	void assembleVerticesList (packet pack, int *vertex_index);
	
	/** This function, together with assembleVerticesList(), assembles the cross-thread intersection data */
	void convertLocalGridIndicesToGlobalIndices (packet pack, int indexOffset);

	/** Projector routine, used to perform partial or full intersections with boost
	threading routines. */
	void projector(int start,int end);
	
	/** return +1 if outside and -1 if inside for the vertex indexed by vertInd belonging
	to the grid cube given by i,j,k indexes. */
	char getInsidness(int i, int j, int k, int vertInd);
	
	/** Marching cubes vertex interpolation. */
	void vertexInterp(double isolevel, double *p1, double *p2, double valp1, double valp2, VERTEX_TYPE *p);
	
	/** Build triangles within a Marching Cubes cell. */
	// int getTriangles(double *vertexValues, double isolevel, int **triangles, int ix, int iy, int iz,
	// 				 int NX, int NY, int NZ,int *xedge, int *yedge, int *zedge, int *xedge_down, int *yedge_down);
	int getTriangles(double *vertexValues, double isolevel, int **triangles, int ix, int iy, int iz, int NX, int NY,int NZ);

	/** Return the vertices for a given section on z of the grid. */
	void getVertices(double isolevel,int start_z,int end_z,int jump,vector<coordVec>*,vector<VERTEX_TYPE*>*);

	/** Update in parallel the neighbour vertex lists useful to approximateNormals() */
	void updateVertexTrianglesLists(int **vertexTrianglesList,double **planes,unsigned int vertex_flag[],bool doOnlyList);

	/** Approximate normals based on the surrounding triangle planes normals. */
	void approximateNormals(vector<int> &appNormals,bool doOnlyList);
	
	/** Multi-threaded triangulator. */
	double triangulationKernel(double isolevel,bool revert,int start_z,int end_z,int jump,vector<int> *localTriList,VERTEX_TYPE *localArea);
	
	/** Builds a 3D grid for accelerating nearest atom queries. */
	void buildAtomsMap(void);
	
	/** Deallocate the memory of the 3D nearest atom query */
	void disposeAtomsMap(void);
	
	/** Apply multidielectric correction after grid building. For each internal
	grid point detect the nearest atom and apply its dielectric constant.*/
	void applyMultidielectric(void);
	
	/** swap the state of a point in the epsmap from internal to the nearest atom dielectric*/
	void swap2multi(double gxmin,double gymin,double gzmin,double gside,unsigned int ggrid,vector<int> *gridMultiMap,int i,int j,int k,int l);
	
	/** Build stern layer. */
	void buildSternLayer(void);
	
	/** Clean and alloc intersections. */
	void allocIntersectionsMatrices(int octree_side_size=0);
	void deleteIntersectionsMatrices(void);
	
	#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT) && !defined(AVOID_NORMALS_MATRICES)
	/** Clean and alloc normals. */
	void allocNormalsMatrices(int octree_side_size=0);
	void deleteNormalsMatrices(void);
	#endif
	
	/** Based on provided info, this decides the best format to save the mesh. */
	int deduceFormat(void);

	/** This returns true if the point is completely out: completely
	out means that all its 1-neighbours are out. */
	bool isCompletelyOut(VERTEX_TYPE *pos);

	/** Buffers useful to the patch-based RT. */
	int *num_pixel_intersections;
	pair<VERTEX_TYPE,VERTEX_TYPE*> **pixel_intersections;

	#if !defined(SINGLE_PASS_RT)
	#if defined(MULTITHREADING)
	int *num_patch_intersections, *first_patch_intersection_index;

	pair<VERTEX_TYPE,VERTEX_TYPE*> *temp_intersections_buffer;

	int64_t *intersection_pixel_id;
	#endif
	#else // SINGLE_PASS_RT
	vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> *temp_intersections_buffer;

	vector<int64_t> *intersection_pixel_id;
	#endif // SINGLE_PASS_RT

	#if defined(RAY_VS_CELL_TESTS_CULLING) || defined(RAY_VS_PATCH_TESTS_CULLING)
	int *tagged_pixel_ids;

	bool *screen_buffer;
	#endif

protected:
	
	Surface();
	Surface(ConfigFile *cf);

public:

	// In order to use the set of default surface methods, the following
	// interface methods must be provided.
	// If a fully custom solution is built the mandatory methods can be fake methods
	// and all the work can be made in getSurf();

	//////////////////////// INTERFACE MANDATORY METHODS //////////////////////////////////
	
	virtual int getNumPatches(void) = 0;

	virtual bool isPatchBasedRayTracingSupported(void) = 0;

	/** Build the surface internal representation. Returns if the
	building process was succesful or not. */
	virtual bool build(void) = 0;
	
	/** Save the surface in a specific class dependent format. Return
	if saving was succesful or not. */
	virtual bool save(char *fileName) = 0;
	
	/** Load the surface in a specific class dependent format. Return
	if loading was succesful or not. */
	virtual bool load(char *fileName) = 0;
	
	/** Print a summary of the surface type, status and other stuff. */
	virtual void printSummary(void) = 0;
	
	/** Get a projection of a point on the surface. Return projection and normal. */
	virtual bool getProjection(double p[3],double *proj1,double *proj2,double *proj3,
										double *normal1,double *normal2,double *normal3)=0;
	
	/** Get all the intersections of a ray that goes from P1 to P2 over the surface.
	The intersections are returned with increasing distance order.
	the double in the vector is the t parameter for the intersection of the parametric
	line and the surface, the double pointer is the normal vector. The management of the memory
	of the normal is up to the derived class from surface. */
	virtual void getRayIntersection(double p1[3], double p2[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id) = 0;
	#if defined(REPORT_FAILED_RAYS)
	virtual void printRayIntersection(double p1[3], double p2[3]) = 0;
	#endif

	/** Function for the constructor without arguments. */
	virtual void init(void);

	/** functions for the constructor with config file argument. */
	virtual void init(ConfigFile *cf);

	/** Function for the denstructor. */
	virtual void clear(void);
	/////////////////////////////////////////////////////////////

	//////////////////////// PROVIDED DEFAULT METHODS /////////////////////////////////
	/** Build DelPhi Surface by a specific method. If the defined surface
	 provides at least projection and intersection routines then don't overload
	 this function; it will already provide most of the information needed by DelPhi.
	 The only restrictions are that only two dielectrics are allowed (in/out) and
	 there is no Stern layer. Overload this method if the surface does not provide
	 intersect or project or all information is needed. This mode of
	 operation is used by the MeshSurface class which provides getRayIntersection and
	 getProjection primitives. If requested one can fill cavities by fill flag. In order to
	 use parallel execution ray-tracing is performed using the itersectionFunctor.

	 SD PB_NEW: added the intersectionsInfo input. If different from nullptr than the intersections and normals data is loaded into this vector. */
	virtual bool getSurf(double *surf_vol,bool optimize_grids,bool fill=false,double cav_vol=0,vector<packet> *intersectionsInfo=nullptr);
	
	/** Compute the cavities of the surface by grid flooding. A list of cavities is returned 
	where each list is a list of triplet of indexes. It is up to the caller to free the memory
	of the int *pointers. The input decide wich is the first STATUS code to be checked.*/
	virtual int getCavities(int idStart=STATUS_POINT_TEMPORARY_OUT);
	virtual int getCavitiesWithBilevelStatusMap(int idStart=STATUS_POINT_TEMPORARY_OUT);
	
	/** Fill the previously detected cavities if their volume is bigger that the passed var.
	In the baseline implementation the volume is approximated by the number of grid cubes
	volumes in that cavity. By default all cavities are filled. */
	virtual void fillCavities(double vol=0,bool silent=false);
	
	/** After the cavities are detected and marked this routine filter out the part or the entire
	cavities that are not able to fit the bounding box of a water molecule. That is we filter
	the bad shaped cavities. */
	virtual void filterCavities(bool modStatus=false);
	virtual void filterCavitiesWithBilevelStatusMap(bool modStatus=false);
	
	/** Swap cavity status to STATUS_POINT_TEMPORARY_OUT. */
	virtual void cav2out(void);
	virtual void cav2outWithBilevelStatusMap(void);
	
	/** Triangulate the surface and save it in OFF format. In the Surface class the baseline method
	for triangulation is obtained by employing the marching cube method at each delphi grid cell.
	For each vertex in the grid cube its insideness is computed by voting of the insidness values
	of the incident cubes; thus the scalar field is given by the ensemble of the status map.
	The surface is saved in .off format. */
	virtual double triangulateSurface(bool outputMesh,bool buildAtomsMapHere,double iso=0.0,
									  const char *fileName="triangulatedSurf",bool revert=false);
	
	/** Save mesh in a prescribed format, revert triangles (change plane sign) if requested. */
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	virtual bool savePLYMesh(int format,bool revert,const char *fileName,
							 vector<VERTEX_TYPE*> &vertList,vector<int> &triList,vector<VERTEX_TYPE*> &normalsList);
	virtual bool saveMesh(int format,bool revert,const char *fileName,
						  vector<VERTEX_TYPE*> &vertList,vector<int> &triList,vector<VERTEX_TYPE*> &normalsList);
	virtual bool saveMeshBinary(int format, bool revert, const char *fileName,
								vector<VERTEX_TYPE*> &vertList, vector<int> &triList, vector<VERTEX_TYPE*> &normalsList);
	#else
	virtual bool savePLYMesh(int format,bool revert,const char *fileName,
							 vector<VERTEX_TYPE> &vertList,vector<int> &triList,vector<VERTEX_TYPE> &normalsList);
	virtual bool saveMesh(int format,bool revert,const char *fileName,
						  vector<VERTEX_TYPE> &vertList,vector<int> &triList,vector<VERTEX_TYPE> &normalsList);
	virtual bool saveMeshBinary(int format, bool revert, const char *fileName,
								vector<VERTEX_TYPE> &vertList, vector<int> &triList, vector<VERTEX_TYPE> &normalsList);
	#endif

	void getTempVerticesAndNormals (VERTEX_TYPE tempVertices[],VERTEX_TYPE tempNormals[],int numNeighbours[],int thread_id,int num_threads);
	void getTempVerticesAndNormals (VERTEX_TYPE tempVertices[],VERTEX_TYPE tempNormals[],int numNeighbours[]);

	void getNewVerticesAndNormals (VERTEX_TYPE tempVertices[],VERTEX_TYPE tempNormals[],int numNeighbours[],int thread_id,int num_threads);

	/** Smooth a given mesh and overwrites the given file name. The input/output mesh is in .off format. */
	virtual void smoothSurface(bool outputMesh,bool buildAtomsMapHere,
							   const char *fn="triangulatedSurf",bool revert=false);

	/** This function is called before the ray tracing of the panel. It can be useful if a
	per panel pre-processing step is needed. By default this function does nothing. */
	virtual void preProcessPanel(void)
	{}
	
	/** This function is called after ray casting to perform any post processing, such as memory
	clean-up, after ray casting*/
	virtual void postRayCasting(void)
	{}
	
	/** This function is called before boundary grid projection to perform any pre-processing,
	such as memory setup-up. */
	virtual bool preBoundaryProjection(void)
	{
		return true;
	}

	/** Functions useful for patch-based ray tracing. */
	/** First patch-based RT function; there are buffer allocations and initialisations. */
	void initPatchBasedRayTracing (int64_t nxyz[], int panels[2]);
	#if !defined(SINGLE_PASS_RT)
	/** Function used in multi-RT mode only; every patch is pierced by rays from multiple directions (piercing
	is done solely if MINIMIZE_MEMORY is defined in globals.h, otherwise conservative, approximated buffer
	offsets are estimated by approximated patch cells' data. The function is run by multiple threads */
	virtual void getPatchPreIntersectionData (int64_t nxyz[3], int panels[2], int thread_id, int potentialIntersections[]) {}
	/** Function with other allocations (dependent on the counts obtained in getPatchPreIntersectionData())
	and initialisations; this function, like getPatchPreIntersectionData(), is not used in the
	single-pass version. */
	void allocateIntermediatePatchBuffers (int64_t nxyz[], int panels[2], int potentialIntersections);
	#endif
	/** Core patch-based RT function: each patch is pierced through all the needed directions and the pertinent
	buffers are filled. This function is run by multiple threads. */
	virtual void getPatchIntersectionData (int64_t nxyz[], int panels[2], int thread_id, int *netIntersections) {}
	// Patch normals at intersections can be computed here a posteriori; this can avoid normal data leakage when cleaning.
	// virtual void getPatchNormalsAtIntersections (int64_t nxyz[], int panels[2], int thread_id) {}
	/** This function allocates screen- and ray-centric buffers necessary to setVerticesAndGridsWithIntersectionData(...)
	where vertices, normals and grid data (these, e.g, the status map, many not be required) are stored. */
	void countPixelIntersections (int64_t nxyz[], int panels[2]);
	/** Patch-centric details are converted to ray-centric details via copies realised through the panel+pixel data
	stored in the buffer intersection_pixel_id[]; buffers which are not needed anymore are deallocated. */
	void copyPatchIntersections (int64_t nxyz[], int panels[2]);
	/** Function employed to sort the intersections; this function is run by multiple threads. */
	void reorderPatchIntersections (int pixel_start, int pixel_end);


	/** Translate a triangulation (a set of samples) into a set of balls (ray-casting sampling + power crust). */
	void tri2Balls(void);
	
	/** save current status map in a temporary map*/
	void backupStatus(void);
	
	/** remove temporary status*/
	void removeBackupStatus(void);
	
	/** If between two any cavities there is communication in map st2 and this communication is all internal in map st1 these cavities are merged logically. */
	int linkCavities(int *st1, int *st2);
	int linkCavities(int **st1, int **st2);
	
	/** difference operator: said S1 the first surface (big probe), S2 (small probe) the second surface and S3 the
	 output surface then the difference rule is the following rule: if S1 is not OUT and and S2 is OUT then that's the pocket. */
	bool difference(Surface *);
	bool differenceWithBilevelStatusMap(Surface *);

	/** This is a Connolly regularized difference operator. Difference function is called
	 and then Connolly filter is applied such that the noise can be filtered out. */
	Surface &operator -= (Surface &surf2);

	/** for each cavity/pocket detect the atoms that compose the cavity. */
	void getCavitiesAtoms(void);
	void getCavitiesAtomsWithBilevelStatusMap(void);

	/** Given a triangulation provided by the current surface object, this function returns true
	 if the triangulation points are completely outside or not according to the status map given in surf. */
	void triangulationPointsAreCompletelyOut(Surface *surf,vector<bool> &results);

	/** Save selected points of the triangulation:
	 *  - if filen is provided, write XYZ to file and XYZN to filen (legacy mode);
	 *  - otherwise write a single file with XYZ(+N when available). */
	void saveSelectedPoints(vector<bool> &selection,const char *file,const char *filen = NULL);

	/** Save a subset of the triangulation into triSubset file. The subset is given in a vertex
	wise way by the vector of boolean results. Revert must be true if the mesh comes from the
	triangulation of a pocket or a cavity. The value of the new area is returned. .*/
	double saveTriSubSet(char *triSubset,vector<bool> &results,bool revert);

	/** Access triangulated mesh geometry for external post-processing utilities. */
	int getTriangulatedNumVertices(void) const;
	int getTriangulatedNumTriangles(void) const;
	void getTriangulatedVertex(int idx,double out[3]) const;
	void getTriangulatedTriangle(int triIdx,int &v0,int &v1,int &v2) const;

	/** Nearest-atom query helper and atom-count accessor. */
	bool nearestAtomForPoint(double *pos,int &winner);
	int getNumLoadedAtoms(void) const;
	int getLoadedAtomOutputSerial(int atomIdx) const;

	/** Public wrappers for atoms-map lifecycle used by external orchestration code. */
	void prepareAtomsMap(void);
	void clearAtomsMap(void);

	/** gives MC cube index. If !=-1 then that cube contains triangles. */
	int classifyCube(double *vertexValues,double isolevel);

	///////////////////////////////////////////// getters/setters //////////////////////////////////////////////

	/** Set if surface has to project bgp or not. During cavity detection/surface visualization/triangulation
	bgp projection can be skipped*/
	void setProjBGP(bool flag)
	{
		projBGP = flag;
	}

	bool getProjBGP(void)
	{
		return projBGP;
	}

	double getSternLayer(void)
	{
		return sternLayer;
	}

	void setSternLayer(double l)
	{
		if (l<0)
		{
			cout << endl << WARN << "Cannot set a negative Stern Layer";
			return;
		}
		else
			sternLayer = l;
	}

	void setSaveMSMS(bool m)
	{
		saveMSMS = m;
	}

	bool getSaveMSMS(void)
	{
		return saveMSMS;
	}

	void setSavePLY(bool m)
	{
		savePLY = m;
	}

	bool getSavePLY(void)
	{
		return savePLY;
	}

	void setAccurateTriangulationFlag(bool flag)
	{
		accurateTriangulation = flag;
	}

	bool getAccurateTriangulationFlag(void)
	{
		return accurateTriangulation;
	}

	void setTriangulationFlag(bool flag)
	{
		doTriangulation = flag;
	}

	bool getTriangulationFlag(void)
	{
		return doTriangulation;
	}

	void setVertexAtomsMap(bool f)
	{
		vertexAtomsMapFlag = f;
	}

	bool getVertexAtomsMap(void)
	{
		return vertexAtomsMapFlag;
	}

	void setComputeNormals(bool cn)
	{
		computeNormals = cn;
	}

	bool getComputeNormals(void)
	{
		return computeNormals;
	}

	void setCheckDuplicatedVertices(bool cd)
	{
		checkDuplicatedVertices = cd;
	}

	bool getCheckDuplicatedVertices(void)
	{
		return checkDuplicatedVertices;
	}

	void setKeepWellShapedCavities(bool kwsc)
	{
		wellShaped = kwsc;
	}

	bool getKeepWellShapedCavities(void)
	{
		return wellShaped;
	}

	virtual void setProbeRadius(double probeRadius)
	{
		probe_radius = probeRadius;
	}

	virtual double getProbeRadius(void)
	{ 
		return probe_radius;
	}

	virtual double getRandDisplacement(void)
	{
		return randDisplacement;
	}

	virtual void setRandDisplacement(const double r)
	{
		randDisplacement = r;
	}

	virtual void setLoadBalancing(bool doLoadBalancing)
	{
		useLoadBalancing = doLoadBalancing;
	}

	virtual bool getLoadBalancing(void)
	{
		return useLoadBalancing;
	}

	virtual int getNumTriangles(void)
	{
		return (int)(triList.size() / 3.);
	}

	virtual int getNumVertices(void)
	{
		return (int)vertList.size();
	}

	virtual void setInsideCode(int i)
	{
		inside = i;
	}

	virtual int getInsideCode(void)
	{
		return inside;
	}

	/** Function useful to chose between conventional ray-centric code and newer patch-based version. */
	void setRayTracingAlgorithm (bool is_patch_based)
	{
		patchBasedAlgorithm = is_patch_based;
	}

	void setRayVsTorusIntersectionAlgorithm (bool analytical_torus_intersection)
	{
		analyticalTorusIntersectionAlgorithm = analytical_torus_intersection;
	}

	void setCollectFaceIntersections (bool collect_intersections)
	{
		collectFaceIntersections = collect_intersections;
	}

	void setParallelBuildupHaloThickness (double parallel_buildup_halo_thickness)
	{
		parallelBuildupHaloThickness = parallel_buildup_halo_thickness;
	}

	void setForceSerialBuild (bool serial_build)
	{
		forceSerialBuild = serial_build;
	}

	void setOptimizeGrids (bool optimize_grids)
	{
		optimizeGrids = optimize_grids;
	}

	void setMaxNumAtoms (int max_atoms)
	{
		maxNumAtoms = max_atoms;
	}

	void setDomainShrinkage (double domain_shrinkage)
	{
		domainShrinkage = domain_shrinkage;
	}

	/** Destructor. */
	virtual ~Surface();


	///////////////////////////////////////////// Visualisation functions //////////////////////////////////////////////

	#if defined(USE_VIS_TOOLS)

	void Rotate (double x1, double y1, double z1,
				 double sn1, double cs1, double sn2, double cs2,
				 double *x2, double *y2, double *z2);
	void AntiRotate (double x1, double y1, double z1,
					 double sn1, double cs1, double sn2, double cs2,
					 double *x2, double *y2, double *z2);

	#ifndef NO_OPENGL
	void openWindow (Vis *vis);
	#endif

	void Projection (Vis *vis);

	void Project (double p1[3], double p2[3], Vis *vis);

	#ifndef NO_OPENGL
	void visualizeString (double r, double g, double b, int x, int y, char *string, void *font);

	void visualizeTriangulation (void);

	void initVisualization (int argc, char *argv[], Vis *vis);
	#endif

	#ifndef NO_OPENGL
	void rotateViewpoint (double dx, double dy, Vis *vis);

	void initMenu (Vis *vis);

	#endif // NO_OPENGL

	#endif // USE_VIS_TOOLS

#if defined(USE_VIS_TOOLS)

void Display (void);

#if !defined(NO_OPENGL)
void endVisualization (Vis *vis);

void Reshape (GLsizei w, GLsizei h);

void motionFunction (int x, int y);

void mouseFunction (int button, int state, int x, int y);

void keybordFunction (unsigned char key, int x, int y);

#endif

#endif
};


void processMenuEvents (int option);

void changeView (int option);


#endif
