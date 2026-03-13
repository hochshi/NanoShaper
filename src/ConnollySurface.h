//---------------------------------------------------------
/**    @file	ConnollySurface.h
*     @brief	ConnollySurface.h is the header for CLASS
*               ConnollySurface.cpp								*/
//---------------------------------------------------------

#ifndef ConnollySurface_h
#define ConnollySurface_h

#include "Surface.h"
#include "SurfaceFactory.h"
#include <complex>

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif


// #define DEBUG_CONNOLLY

#define MAX_INCIDENT_PROBES 50

#define DEFAULT_PROBE_RADIUS 1.4 // default probe radius

#if defined(CHECK_BUILDUP_DIFF)
#define EPS_BUILDUP 1E-10
#endif


class ConnollyCell{
public:
	int patch_type;

	#if defined(CHECK_BUILDUP_DIFF)
	// Just a value to be able to order patches and compare multi-thread vs single-thread data
	long int tag;
	#endif

	virtual ~ConnollyCell()
	{}
};


class FacetCell: public ConnollyCell
{
public:
	double planes[4][4]; // trimming tethraedron
	double center[3];
	int plane_indices[3][2]; // every plane depends on a pair of atoms
	int id[3]; // atoms triplet indices
	vector<double*> self_intersection_planes; // additional self intersection planes
	#if defined(OPTIMIZE_CELL_STRUCTURE)
	vector<bool> self_intersection_plane_labels; // labels indicating the signs
	#endif
	FacetCell *mirrorCell;
	bool isSelfIntersecting;

	#if defined(CHECK_BUILDUP_DIFF)
	// tag of facet cell = center[0] + center[1] + center[2]
	#endif

	FacetCell()
	{
		mirrorCell = NULL;
	}

	virtual ~FacetCell()
	{
		for (unsigned int i=0; i<self_intersection_planes.size(); i++)
		{
			#if !defined(OPTIMIZE_CELL_STRUCTURE)
			deleteVector<double>(self_intersection_planes[i]);
			#else
			if (self_intersection_plane_labels[i])
				deleteVector<double>(self_intersection_planes[i]);
			#endif
		}
		self_intersection_planes.clear();
		#if defined(OPTIMIZE_CELL_STRUCTURE)
		self_intersection_plane_labels.clear();
		#endif
	}
};


class EdgeCell: public ConnollyCell
{
public:
	double Rot[3][3]; // rotation matrix of the torus
	double invrot[3][3]; // inverse rotation matrix
	double center[3];
	double clipping_center[3]; // center of the clipping sphere of the torus
	double clipping_radius;
	double major_radius;
	// double u[3], v[3];dir

	#if !defined(OPTIMIZE_CELL_STRUCTURE)
	double cutting_planes[2][4]; // two cutting planes of the torus
	#else
	// pointers to plane data; attention to the sign (see the use)
	double *cutting_planes[2];
	#endif
	// if OPTIMIZE_CELL_STRUCTURE is defined, the vector will be of pointers to existing planes
	// but the sign is changed on demand in the pertinent point of ConnollySurface.cpp
	vector<double*> additional_planes; // additional clipping planes due to tori clipping by singular facets
	vector<bool> flags; // acute / non acute flags for additional planes
	double self_intersection_radius; // radius of the clipping sphere in case of spindle torus
	double rcoi;
	int id[2]; // atom pair indices
	// bool isSelfIntersecting;
	bool acute; // if angle between planes is acute perform "and" side test between planes, otherwise "or" side test is sufficient
	bool isComplex;

	#if defined(CHECK_BUILDUP_DIFF)
	// tag of edge cell = center[0] + center[1] + center[2]
	#endif

	virtual ~EdgeCell()
	{
		#if !defined(OPTIMIZE_CELL_STRUCTURE)
		for (unsigned int i=0; i<additional_planes.size(); i++)
			deleteVector<double>(additional_planes[i]);
		#endif
		additional_planes.clear();

		flags.clear();
	}
};


class PointCell: public ConnollyCell
{
public:
	#if !defined(OPTIMIZE_CELL_STRUCTURE)
	vector<EdgeCell*> neighbours;
	vector<FacetCell*> incidentProbes;
	vector<EdgeCell*> buried_neighbours;
	#else
	vector<double> neighbour_data;
	vector<int> incidentProbes;
	vector<double> buried_neighbour_data;
	#endif

	int id; // atom index

	#if defined(CHECK_BUILDUP_DIFF)
	// tag of point cell = id + RAND_DISPLACEMENT
	#endif

	virtual ~PointCell()
	{
		#if !defined(OPTIMIZE_CELL_STRUCTURE)
		neighbours.clear();

		for (unsigned int i=0; i<buried_neighbours.size(); i++)
			delete buried_neighbours[i];
		buried_neighbours.clear();

		incidentProbes.clear();
		#else
		neighbour_data.clear();
		buried_neighbour_data.clear();
		incidentProbes.clear();
		#endif
	}
};


/*
#if !defined(OPTIMIZE_GRIDS)
// 2d map
#define GRID_CONNOLLY_CELL_MAP_2D(i,j,l,NA,NB) gridConnollyCellMap2D[((j)*NA + (i))*MAX_CONNOLLY_CELLS_2D + (l)]

// 3d map
#define GRID_CONNOLLY_CELL_MAP(i,j,k,l,NX,NY,NZ) gridConnollyCellMap[((k)*((NY)*(NX)) + (j)*(NX) + (i))*MAX_CONNOLLY_CELLS + (l)]

// #define SELF_MAP(i,j,k,l,NX,NY,NZ) gridProbesMap[((k)*((NY)*(NX)) + (j)*(NX) + (i))*MAX_PROBES + (l)]
#endif
*/


#ifdef ENABLE_CGAL
//////////////////////// CGAL ///////////////////////////////////////////////////
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
// #include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/intersections.h>
#include <cassert>
#include <vector>
#include <fstream>
#ifdef CGAL_USE_BASIC_VIEWER
	#include <CGAL/draw_triangulation_3.h>
#endif
////////////////////////////////////////////////////////////////////////////////
#endif

// here singular/regular is in the alpha shape nomenclature

#define POINT_CELL				0
#define REGULAR_EDGE_CELL		1
#define SINGULAR_EDGE_CELL		2
#define REGULAR_FACE_CELL		3
#define SINGULAR_FACE_CELL		4
#define SKIP_CELL				5


/** @brief This class builds and converts to a DelPhi suitable representation the Connolly-Richards Surface.
All the gathered info is analytically computed both the intersections and the projections. First the
alpha shape of the set of atoms is computed, then from that the Connolly surface is computed by build all
the patches.

@author Sergio Decherchi 
@date 15/05/2012
*/
class ConnollySurface: public Surface
{

#ifdef ENABLE_CGAL 
private:

	typedef CGAL::Exact_predicates_inexact_constructions_kernel 		K;
	// typedef CGAL::Regular_triangulation_filtered_traits_3<K>			Traits;
	typedef CGAL::Polyhedron_3<K>										Polyhedron;
	typedef Polyhedron::Halfedge_around_facet_circulator        		HF_circulator;
	typedef CGAL::Vector_3<K>											Vector3;
	typedef CGAL::Plane_3<K>											Plane3;
	typedef CGAL::Ray_3<K>												Ray3;
	typedef CGAL::Line_3<K>												Line3;

	#if !defined(NO_CGAL_PATCHING)
	typedef CGAL::Regular_triangulation_vertex_base_3<K>                    Vertex_with_info;
	#else
    // nopatch
	typedef CGAL::Triangulation_vertex_base_with_info_3<int,K> 				_Vertex_with_info;
	typedef CGAL::Regular_triangulation_vertex_base_3<K,_Vertex_with_info>	Vertex_with_info;
	#endif
	typedef CGAL::Fixed_alpha_shape_vertex_base_3<K,Vertex_with_info>		Vb;

	typedef CGAL::Regular_triangulation_cell_base_3<K> 						RT_Cell_with_info;
	typedef CGAL::Fixed_alpha_shape_cell_base_3<K,RT_Cell_with_info>		Fb;

	typedef CGAL::Triangulation_data_structure_3<Vb,Fb,CGAL::Parallel_tag>	Tds;
	typedef Tds::Cell_circulator											Cell_circulator;

	typedef CGAL::Regular_triangulation_3<K, Tds>							Rt;

	typedef CGAL::Fixed_alpha_shape_3<Rt>								Fixed_alpha_shape_3;
	typedef Fixed_alpha_shape_3::Cell_handle							Alpha_Cell_handle;
	typedef Fixed_alpha_shape_3::Vertex_handle							Alpha_Vertex_handle;
	typedef Fixed_alpha_shape_3::Facet									Alpha_Facet;
	typedef Fixed_alpha_shape_3::Edge									Alpha_Edge;
	typedef Rt::Bare_point												Point3;
	typedef Rt::Weighted_point											Weighted_point;

	typedef Rt::Vertex_iterator											Vertex_iterator;
	typedef Rt::Finite_vertices_iterator                        		Finite_Vertex_Iterator;
	typedef Rt::Finite_cells_iterator									Finite_Cells_Iterator;
	typedef Rt::Finite_edges_iterator									Finite_Edges_Iterator;
	typedef Rt::Finite_facets_iterator									Finite_Facets_Iterator;
	typedef Rt::Vertex_handle											Vertex_handle;
	typedef Rt::Facet													Facet;

	// typedef K::FT													Weight;

	/*
	// OLD
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Regular_triangulation_filtered_traits_3<K>	Traits;
	typedef CGAL::Polyhedron_3<K>				                Polyhedron;
	typedef Polyhedron::Halfedge_around_facet_circulator        HF_circulator;
	typedef CGAL::Vector_3<K>									Vector3;
	typedef CGAL::Plane_3<K>									Plane3;
	typedef CGAL::Ray_3<K>									    Ray3;
	typedef CGAL::Line_3<K>									    Line3;
	typedef CGAL::Triangulation_vertex_base_with_info_3<PointCell*,Traits> Vertex_with_info;
	typedef CGAL::Fixed_alpha_shape_vertex_base_3<Traits,Vertex_with_info>		Vb;
	typedef CGAL::Triangulation_cell_base_with_info_3< FacetCell*[4],Traits> Cell_with_info;
	typedef CGAL::Regular_triangulation_cell_base_3<Traits,Cell_with_info> RT_Cell_with_info;
	typedef CGAL::Fixed_alpha_shape_cell_base_3<Traits,RT_Cell_with_info>         Fb;
	typedef CGAL::Triangulation_data_structure_3<Vb,Fb>			Tds;
	typedef CGAL::Regular_triangulation_3<Traits,Tds>			Rt;
	typedef CGAL::Fixed_alpha_shape_3<Rt>				        Fixed_alpha_shape_3;
	typedef Fixed_alpha_shape_3::Cell_handle					Alpha_Cell_handle;
	typedef Fixed_alpha_shape_3::Vertex_handle					Alpha_Vertex_handle;
	typedef Fixed_alpha_shape_3::Facet							Alpha_Facet;
	typedef Fixed_alpha_shape_3::Edge							Alpha_Edge;
	typedef Traits::Bare_point                                  Point3;
	typedef Traits::Weighted_point                              Weighted_point;

	typedef Rt::Vertex_iterator                                 Vertex_iterator;
	typedef Rt::Finite_vertices_iterator                        Finite_Vertex_Iterator;
	typedef Rt::Finite_cells_iterator							Finite_Cells_Iterator;
	typedef Rt::Finite_edges_iterator							Finite_Edges_Iterator;
	typedef Rt::Finite_facets_iterator							Finite_Facets_Iterator;
	typedef Rt::Vertex_handle                                   Vertex_handle;
	typedef Rt::Facet											Facet;
	*/


	#ifdef MULTITHREADED_SES_BUILDING
	/** Local per-thread data mostly focused on cell and atom data. */
	struct ThreadDataWrapper
	{
		vector<ConnollyCell*> sesComplex;

		#if !defined(NEW_ATOM_PATCHES)
		PointCell **atomPatches;
		#else
		int *atomPatches;
		#endif

		#if !defined(NO_CGAL_PATCHING)
		list<Weighted_point> l;
		#else
		// nopatch
		list<std::pair<Weighted_point, int>> l;
		#endif

		vector<int> exposed;

		int num_cells[5];

		#if defined(OPTIMIZE_BUILDING_MEMORY)
		int my_task_id, my_thread_id;
		#endif
	};


	/** Structure useful to gather data. */
	struct TaggedDataWrapper
	{
		double **minislab_max_y;
		double *slab_max_z;

		int **thread_in_minibin;
		int *task_in_bin;
	};
	#endif

#endif // ENABLE_CGAL


private:

	int numSkippedCells;
	/** number of threads during cells' building, if any. */
	int numThreadDataWrappers = 1;
	int numThreadDataWrappersPerTask[MAX_TASKS];

	/** number of Connolly cells for type */
	int type[5];
	/** link each auxuliary grid box to the respective cell list */
	/* #if !defined(OPTIMIZE_GRIDS)
	int *gridConnollyCellMap;
	int *gridConnollyCellMap2D;
	#else */
	#if !defined(MULTITHREADED_SES_BUILDING)
	vector<int> *gridConnollyCellMap;
	vector<int> *gridConnollyCellMap2D;
	#else
	vector<pair<int,int>> *gridConnollyCellMap;
	vector<pair<int,int>> *gridConnollyCellMap2D;
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
	/** auxiliary grid sizes for the 2d map*/
	int64_t nx_2d,ny_2d,nz_2d;
	double scale_2d;
	double side_2d;
	double xmin_2d,xmax_2d,ymin_2d,ymax_2d,zmin_2d,zmax_2d;
	/* #if !defined(OPTIMIZE_GRIDS)
	unsigned int **ind_2d;
	#endif */

	unsigned int AUX_GRID_DIM_CONNOLLY;
	unsigned int MAX_CONNOLLY_CELLS;
	unsigned int AUX_GRID_DIM_CONNOLLY_2D;
	unsigned int MAX_CONNOLLY_CELLS_2D;

	#if !defined(MULTITHREADED_SES_BUILDING) || defined(CHECK_BUILDUP_DIFF)
	#if !defined(NEW_ATOM_PATCHES)
	/** this vector contains NULL if the atom is not exposed and the point cell if the atom contributes to the surface*/
	PointCell **atomPatches;
	#else
	int *atomPatches;
	#endif

	#if defined(CHECK_BUILDUP_DIFF)
	#if !defined(NEW_ATOM_PATCHES)
	/** this vector contains NULL if the atom is not exposed and the point cell if the atom contributes to the surface*/
	PointCell **st_atomPatches;
	#else
	int *st_atomPatches;
	#endif
	#endif // CHECK_BUILDUP_DIFF
	#endif // !defined(MULTITHREADED_SES_BUILDING) || defined(CHECK_BUILDUP_DIFF)
	
	#if defined(MULTITHREADED_SES_BUILDING)
	/** local, per-thread data needed in the building stage */
	// ThreadDataWrapper *thread_data_wrapper;
	ThreadDataWrapper thread_data_wrapper[MAX_TASKS_TIMES_THREADS];
	#endif

	/** for each cell there is a structure that defines the patch. */
	vector<ConnollyCell*> sesComplex;
	
	/** compute the connolly surface using the CGAL alpha shape module and compute all information
	needed by to ray-trace it*/
	bool buildConnolly(void);
	/** map each cell to the auxiliary grid. */
	bool buildAuxiliaryGrid();

	#ifdef ENABLE_CGAL
	#ifdef MULTITHREADED_SES_BUILDING
	/** Function used to reorder patch cells according to the clipping box radiuses */
	void reorderPatchCells (DelPhiShared *delphi, ThreadDataWrapper *tdw, int cell_start, int cell_end);
	#endif
	/** CGAL build Connolly surface. */
	bool buildConnollyCGAL(void);
	#ifdef MULTITHREADED_SES_BUILDING
	void BuildWeightedPoints (TaggedDataWrapper *tagged_data_wrapper, double slab_pars[], int thread_pars[]);
	#endif
	/** This function builds point cells. */
	#if !defined(NEW_ATOM_PATCHES)
	void BuildPointCells (Fixed_alpha_shape_3 &alpha_shape, vector<ConnollyCell*> &sesComplex, PointCell **atomPatches,
						  vector<int> &exposed, int num_cells[]
						  #if defined(OPTIMIZE_BUILDING_MEMORY)
						  , TaggedDataWrapper *tagged_data_wrapper, int task_id, int thread_id, double grid_pars[]
						  #endif
						  );
	/** This function builds facet, prismatic cells */
	void BuildFacetCells (Fixed_alpha_shape_3 &alpha_shape, vector<ConnollyCell*> &sesComplex, Octree<vector<FacetCell*>> &gridProbesMap,
						  PointCell **atomPatches, int num_cells[], grid_pars);
	/** This function builds edge, prismatic cells */
	void BuildEdgeCells (ofstream &of, Fixed_alpha_shape_3 &alpha_shape, vector<ConnollyCell*> &sesComplex, PointCell **atomPatches,
						 int num_cells[]
						 #if defined(OPTIMIZE_BUILDING_MEMORY)
						 , TaggedDataWrapper *tagged_data_wrapper, int task_id, int thread_id, double grid_pars[]
						 #endif
						 );
	#else // NEW_ATOM_PATCHES
	void BuildPointCells (Fixed_alpha_shape_3 &alpha_shape, vector<ConnollyCell*> &sesComplex, int atomPatches[],
						  vector<int> &exposed, int num_cells[]
						  #if defined(OPTIMIZE_BUILDING_MEMORY)
						  , TaggedDataWrapper *tagged_data_wrapper, int task_id, int thread_id, double grid_pars[]
						  #endif
						  );
	/** This function builds facet, prismatic cells */
	void BuildFacetCells (Fixed_alpha_shape_3 &alpha_shape, vector<ConnollyCell*> &sesComplex, Octree<vector<FacetCell*>> &gridProbesMap, int atomPatches[],
						  int num_cells[], double grid_pars[]);
	/** This function builds edge, prismatic cells */
	void BuildEdgeCells (ofstream &of, Fixed_alpha_shape_3 &alpha_shape, vector<ConnollyCell*> &sesComplex, int atomPatches[],
						 int num_cells[]
						 #if defined(OPTIMIZE_BUILDING_MEMORY)
						 , TaggedDataWrapper *tagged_data_wrapper, int task_id, int thread_id, double grid_pars[]
						 #endif
						 );
	#endif // NEW_ATOM_PATCHES
	/** This function remove self intersections. */
	void RemoveSelfIntersections (vector<ConnollyCell*> &sesComplex, Octree<vector<FacetCell*>> &gridProbesMap, double grid_pars[]);
	#ifdef MULTITHREADED_SES_BUILDING
	void BuildPatches (ThreadDataWrapper *tdw, double grid_pars[]
					   #if defined(OPTIMIZE_BUILDING_MEMORY)
					   , TaggedDataWrapper *tagged_data_wrapper
					   #endif
					   );
	#endif

	bool compPatchTags (ConnollyCell *cc1, ConnollyCell *cc2);

	/** This function checks diverging build-up data obtained by the single-thread exec and by that
	reling on multiple slabs and threads (depending on input setting and eventual
	spatial constraints: if the threads are too many then the number of slabs and
	threads in the build-up is truncated). */
	void checkBuildupDivergences (vector<ConnollyCell*> &st_sesComplex, vector<ConnollyCell*> &sesComplex,
								  int local_divergent_values[16][9], int thread_id, int num_threads);
	#endif

	bool savePovRay;
	// DISMISSED
	// int MAX_PROBES;
	/** self intersections grid perfil. Increase this if you use big probes*/
	double si_perfil;

	ofstream output_file;

public:

	/** Default constructor. */
	ConnollySurface();
	/** set DelPhi environment. */
	ConnollySurface(DelPhiShared *ds);
	/** set configuration and DelPhi environment. */
	ConnollySurface(ConfigFile *cf,DelPhiShared *ds);

	//Fixed_alpha_shape_3 *alpha_shape1;

	//////////////////////// INTERFACE MANDATORY METHODS /////////////////////////////////
	virtual int getNumPatches(void);
	virtual bool isPatchBasedRayTracingSupported(void);

	/** Compute connolly surface. Call it after load*/
	virtual bool build(void);
	/** Get a projection of a point on the surface. Return projection and normal*/
	virtual bool getProjection(double p[3],double *proj1,double *proj2,
							   double *proj3,double *normal1,double *normal2,double *normal3);
	/** Get all the intersections of a ray that goes from P1 to P2 over the surface.
	The intersections are returned with increasing distance order.
	the first double in the vector is the t parameter for the intersection of the parametric 
	line and the surface, the double pointer is the normal vector. During ray surface intersection
	the previously built auxiliary grid is used to speed up computations. */
	virtual void getRayIntersection(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id);
	#if defined(REPORT_FAILED_RAYS)
	virtual void printRayIntersection(double pa[3], double pb[3]);
	#endif
	/** function for the constructor without arguments. */
	virtual void init(void);
	/** functions for the constructor with config file argument. */
	virtual void init(ConfigFile *cf);
	/**function for the denstructor. */
	virtual void clear(void);
	/** pre-process panel to accelerate ray-tracing. */
	virtual void preProcessPanel(void);
	virtual void postRayCasting(void);
	virtual bool preBoundaryProjection(void);
	
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
	/** Save it in a simple ASCII format (.ses)*/
	virtual bool save(char *fileName);
	/**Load the surface from a file in .ses format*/
	virtual bool load(char *fileName);
	/** Print number of cells and types*/
	virtual void printSummary(void);

	void setProbeRadius(double pr) 
	{ 
		double e = 1e-2;
		if (pr > e)
		{ 
			probe_radius=pr; 
		}
		else
		{
			cout << endl << WARN << "Cannot set " << pr << "<=" << e << ". Setting "<< DEFAULT_PROBE_RADIUS;
			probe_radius = DEFAULT_PROBE_RADIUS;
		}
	}
	
	double getProbeRadius(void)
	{
		return probe_radius;
	}

	void setSavePovRay(bool ss) 
	{ 
		savePovRay = ss;
	}
	
	bool getSavePovRay(void)
	{
		return savePovRay;
	}

	/**for the 3d grid set the max grid size and the maximal number of patches inside a grid cube*/
	void setAuxGrid(unsigned int dim,unsigned int max)
	{
		AUX_GRID_DIM_CONNOLLY = dim;
		MAX_CONNOLLY_CELLS = max;
	}

	/**for the 2d grid set the max grid size and the maximal number of patches inside a grid cube.
	The grid cube itself does not exist, it is just a reference; the relevant quantity is MAX_CONNOLLY_CELLS_2D
	that is the number of patches along the grid tube*/
	void setAuxGrid2D(unsigned int dim,unsigned int max)
	{
		AUX_GRID_DIM_CONNOLLY_2D = dim;
		MAX_CONNOLLY_CELLS_2D = (max*dim);
	}

	/*
	void setMaxProbes(int m)
	{
		if (m<=0)
		{
			cout << endl << WARN << "Cannot set max probes <0";
			return;
		}
		MAX_PROBES = m;
	}

	int getMaxProbes(void)
	{
		return MAX_PROBES;
	}
	*/

	void setSIPerfil(double si)
	{
		si_perfil = si;
	}

	double getSIPerfil(void)
	{
		return si_perfil;
	}

	virtual ~ConnollySurface();

private:

	bool rayConnollyCellIntersection(double*,double*,ConnollyCell*,double t[4],int &numInt);
	#if defined(REPORT_FAILED_RAYS)
	bool printRayConnollyCellIntersection(double*,double*,ConnollyCell*,double t[4],int &numInt);
	#endif
	/** This gives true if the point is inside the list of planes*/
	bool isFeasible(ConnollyCell *cc,double *point);
	/** project a point in 3D to a circle in 3D. Input are the point, the radius,center and the
	 p lane where circle belongs and the output is the projection and the distance. Assume that                *
	 the normal to the plane is unitary*/
	void projectToCircle(double *point,double radius,double *center,double *plane,double *proj,double &dist);
	// This intersects a ray with a polyhedron. For now, it handles only the case in which patch_type = REGULAR_EDGE_CELL, but it should be corrected
	bool rayCell(ConnollyCell *cc, double point[3], double dir[3], double t[2]);
	/** project a point to a torus defined by torus_equation*/
	void projectToTorus(double *y,EdgeCell *ec,double *proj,double *norm,double &dist);
	/** given a point on the torus in EdgeCell, it gives the normal to that point without computing the gradient
	explicitly*/
	void getNormalToTorus(double *y,EdgeCell *ec,double *normal);
	/** given a point y compute the normal on that point. This routine does not check that
	the y point really belongs to the surface, this should be assured by the user. If not
	assured the result is meaningless*/
	void getNormal(double *y,ConnollyCell *cc,double *normal);
	/** check the orientation. Assume the planes points toward the visible region of the torus*/
	bool orientation(double *pb_center1,double *pb_center2,double *w1,double *w2);
	/** the aim is to sort probes in clockwise order. As reference the first probe is used*/
	void sortProbes(EdgeCell *ec,FacetCell **fcv,int np,int *sorted);
	void getCoi(double *torus_center,double rcoi,double **sampledPoints,int numPoints,double *u,double *v);
	/** save concave sphere patch in Pov-Ray format*/
	void saveConcaveSpherePatch(ofstream &of,FacetCell *fc,int i);
	void saveSphere(ostream &of,double *center,double radius);
	void saveAtomPatch(ofstream &of,PointCell *pc);
	void saveEdgePatch(ofstream &of,EdgeCell *ec,int size,double bigR,double *u,double *v,double *w,bool isComplex);
};


// expand it explicitly because Swig is not able to expand it
static class ConnollySurfaceRegister{
	static Surface *createSurface(ConfigFile *conf,DelPhiShared *ds)
	{ 
		return new ConnollySurface(conf,ds); 
	} 
	public: 
		ConnollySurfaceRegister() 
		{ 
			surfaceFactory().add("ses",createSurface); 
		} 
} ConnollySurfaceRegisterObject;


// static SurfaceRecorder<ConnollySurface> sesRecorder("ses");

#endif
