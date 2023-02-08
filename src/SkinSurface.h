//---------------------------------------------------------
/**    @file	SkinSurface.h
 *     @brief	SkinSurface.h is the header for CLASS
 *               SkinSurface.cpp
 */
//---------------------------------------------------------

#ifndef SkinSurface_h
#define SkinSurface_h

#include <Surface.h>
#include <spdlog/spdlog.h>

// #define DEBUG_SKIN

#define DEFAULT_S 0.45 // default shrink factor

/** mixed cell data structure*/
class MixedCell {
public:
  double *quadric;
  int surface_type;
  vector<double *> mc_points; // these are only pointer, no explicit denstructor

  virtual ~MixedCell() {
    if (quadric != NULL)
      deleteVector<double>(quadric);
  }
};

// voronoi facet cell or del edge
class Del1Cell : public MixedCell {
public:
  int ids[2];              // the two atoms generating the del edge
  vector<double *> planes; // lateral planes

  virtual ~Del1Cell() {
    for (unsigned int i = 0; i < planes.size(); i++)
      free(planes[i]);
  }

  // #IMPROVE one could build pointers from here to the planes computed in
  // Del2Cell or viceversa
  double upper[4]; // upper plane
  double lower[4]; // upper plane
};

// reduced voronoi cell
class Del0Cell : public MixedCell {
public:
  int id; // atom index
  vector<double *>
      planes; // planes pointer that points to upper/lower of Del1Cell

  virtual ~Del0Cell() {
    // do nothing, planes are not managed by this class
  }
};

// del facet cell
class Del2Cell : public MixedCell {
public:
  int ids[3];          // atom indices
  double planes[3][4]; // lateral planes
  double upper[4];
  double lower[4];

  virtual ~Del2Cell() {
    // do nothing, no pointers
  }
};

// reduced thethraedron
class Del3Cell : public MixedCell {
public:
  int ids[4]; // atom indices
  // planes of the reduced tethraedron
  double planes[4][4];
  double planes_or[4][4];
  double points[4][3];  // reduced points in the same order of atom indices
  double vor[3];        // Voronoi center
  double reduced[4][3]; // reduced points in the same order of atom indices

  virtual ~Del3Cell() {
    // do nothing, no pointers
  }
};

#define GRID_MIXEDCELLMAP_2D(i, j, k, NA, NB)                                  \
  gridMixedCellMap2D[(k) + (MAX_MIXEDCELLS_2D - 1) * ((j) + (NB) * (i))]
// fortran like
// #define GRIDMIXEDCELLMAP(i,j,k,l,NX,NY,NZ)
// gridMixedCellMap[(i)+(NX)*((j)+(NY)*((k)+(NZ)*(l)))]
// c like (NX is ignored)
#define GRIDMIXEDCELLMAP(i, j, k, l, NX, NY, NZ)                               \
  gridMixedCellMap[(l) +                                                       \
                   (MAX_MIXEDCELLS - 1) * ((k) + (NZ) * ((j) + (NY) * (i)))]

#ifdef ENABLE_CGAL
//////////////////////// CGAL
//////////////////////////////////////////////////////
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
// #include <CGAL/Regular_triangulation_euclidean_traits_3.h>
// #include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <cassert>
#include <fstream>
#include <vector>

// use to translate in Pov-Ray
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

////////////////////////////////////////////////////////////////////////////////
#endif

#define DELAUNAY_POINT_CELL 0
#define DELAUNAY_EDGE_CELL 1
#define DELAUNAY_FACET_CELL 2
#define DELAUNAY_TETRA_CELL 3

/** @brief This class builds and converts to a DelPhi suitable representation
the Skin Surface. All the gathered info is analytically computed both the
intersections and the projections. In order to get an accurate result for the
projection routine, as root finding algorithm is used the method of the
companion matrix. The Skin surface was defined in: <i> "H. Edelsbrunner.
Deformable smooth surface design. Discrete Comput. Geom., 21:87-115, 1999." </i>

@author Sergio Decherchi
@date 30/10/2012
*/
class SkinSurface : public Surface {

#ifdef ENABLE_CGAL
private:
  // typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  // typedef CGAL::Regular_triangulation_filtered_traits_3<K>    Traits;
  // typedef CGAL::Triangulation_cell_base_with_info_3< Del3Cell*,Traits>
  // Cell_with_info; typedef
  // CGAL::Regular_triangulation_cell_base_3<Traits,Cell_with_info>
  // RT_Cell_with_info; typedef CGAL::Triangulation_vertex_base_3<Traits>
  // Vertex; typedef
  // CGAL::Triangulation_data_structure_3<Vertex,RT_Cell_with_info> tds; typedef
  // Traits::RT                                          Weight; typedef
  // Traits::Bare_point                                  Point3; typedef
  // Traits::Weighted_point                              Weighted_point; typedef
  // CGAL::Regular_triangulation_3<Traits,tds>           Rt;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::FT Weight;
  typedef K::Point_3 Point;
  typedef K::Weighted_point_3 Weighted_point;
  typedef std::tuple<Weighted_point, int> Indexed_weighted_point;
  typedef CGAL::Regular_triangulation_vertex_base_3<K> Vb0;
  typedef CGAL::Triangulation_vertex_base_with_info_3<int, K, Vb0> Vb;
  typedef CGAL::Regular_triangulation_cell_base_3<K> Cb0;
  typedef CGAL::Triangulation_cell_base_with_info_3<Del3Cell *, Cb0, Cb0> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
  typedef CGAL::Regular_triangulation_3<K, Tds> Rt;
  typedef Rt::Vertex_iterator Vertex_iterator;
  typedef Rt::Finite_vertices_iterator Finite_Vertex_Iterator;
  typedef Rt::Finite_cells_iterator Finite_Cells_Iterator;
  typedef Rt::Finite_edges_iterator Finite_Edges_Iterator;
  typedef Rt::Finite_facets_iterator Finite_Facets_Iterator;
  typedef Rt::Vertex_handle Vertex_handle;
  typedef Tds::Cell_circulator Cell_circulator;
  typedef Tds::Cell_handle Cell_handle;
  typedef Rt::Facet Facet;

  // used to translate in Pov-Ray
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef Polyhedron::Halfedge_around_facet_circulator HF_circulator;

  // a functor computing the plane containing a triangular facet
  struct Plane_from_facet {
    Polyhedron::Plane_3 operator()(Polyhedron::Facet &f) {
      Polyhedron::Halfedge_handle h = f.halfedge();
      return Polyhedron::Plane_3(h->vertex()->point(),
                                 h->next()->vertex()->point(),
                                 h->opposite()->vertex()->point());
    }
  };

#endif

private:
  /** number of mixed cells*/
  int numMixedCells;
  /** number of mixed cells for type.
  type 0 is a Delaunay point + Voronoi cell
  type 1 is a Delaunay edge  + Voronoi facet
  type 2 is a Delaunay facet + Voronoi edge
  type 3 is a Delaunay cell  + Voronoi point*/
  int type[4];
  /** link each auxuliary grid box to the respective mixed cell list */
  int *gridMixedCellMap;
  /** auxiliary grid sizes*/
  int nx, ny, nz;
  double scale;
  double side;
  double *x, *y, *z;
  double s;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  unsigned short ***ind;
  int *gridMixedCellMap2D;

  unsigned int AUX_GRID_DIM_SKIN;
  unsigned int MAX_MIXEDCELLS;
  unsigned int AUX_GRID_DIM_SKIN_2D;
  unsigned int MAX_MIXEDCELLS_2D;

  /** auxiliary grid sizes for the 2d map*/
  int nx_2d, ny_2d, nz_2d;
  double scale_2d;
  double side_2d;
  double xmin_2d, xmax_2d, ymin_2d, ymax_2d, zmin_2d, zmax_2d;
  unsigned int **ind_2d;
  /** reduced tet cells that contain usefull reduced points but that does not
   * give a real patch*/
  vector<Del3Cell *> pendingCells;
  /** for each mixed cell there is a set of planes that define it and a quadric
   * equation*/
  MixedCell **mixedComplex;
  /** compute the skin surface using CGAL regular triangulation and compute all
  information needed by to ray-trace it*/
  bool buildSkin();
  /** map each mixed cell to the auxiliary grid*/
  bool buildAuxiliaryGrid();
#ifdef ENABLE_CGAL
  /** CGAL build skin */
  bool buildSkinCGAL();
#endif

  bool savePovRay;
  bool fastProjection;

public:
  /** Default constructor*/
  SkinSurface();
  /** set DelPhi environment*/
  SkinSurface(DelPhiShared *ds);
  /** set configuration and DelPhi environment*/
  SkinSurface(ConfigFile *cf, DelPhiShared *ds);

  //////////////////////// INTERFACE MANDATORY METHODS
  ////////////////////////////////////
  /** Compute skin surface. Call it after load*/
  virtual bool build();
  /** Save it in a simple ASCII format (.skin)*/
  virtual bool save(char *fileName);
  /**Load the surface from a file in .skin format*/
  virtual bool load(char *fileName);
  /** Print number of mixed cells and types*/
  virtual void printSummary();
  /** Get a projection of a point on the surface. Return projection and normal*/
  virtual bool getProjection(double p[3], double *proj1, double *proj2,
                             double *proj3, double *normal1, double *normal2,
                             double *normal3);
  /** Get all the intersections of a ray that goes from P1 to P2 over the
  surface. The interesctions are returned with increasing distance order. the
  first double in the vector is the t parameter for the intersection of the
  parametric line and the surface, the double pointer is the normal vector.
  During ray surface intersection the previously built auxiliary grid is used to
  speed up computations*/
  virtual void getRayIntersection(double p1[3], double p2[3],
                                  vector<pair<double, double *>> &intersections,
                                  int thdID, bool computeNormals);
  /** function for the constructor without arguments*/
  virtual void init();
  /** functions for the constructor with config file argument*/
  virtual void init(ConfigFile *cf);
  /**function for the denstructor*/
  virtual void clear();
  /** pre-process panel to accelerate ray-tracing*/
  virtual void preProcessPanel();
  virtual void postRayCasting();
  virtual bool preBoundaryProjection();
  /////////////////////////////////////////////////////////////

  void setShrinking(double ss) {
    double e = 0.05;
    if (ss <= (1.0 - e) && ss >= (0.0 + e)) {
      s = ss;
    } else {
      spdlog::warn("Cannot set {}. s parameter is in (0+e,1-e], where e is {}.Setting {}", ss, e, DEFAULT_S);
      s = DEFAULT_S;
    }
  }

  void setFastProjection(bool useFastProjection) {
    fastProjection = useFastProjection;
  }

  double getShrinking() { return s; }

  void setSavePovRay(bool ss) { savePovRay = ss; }

  bool getSavePovRay() { return savePovRay; }

  /**for the 3d grid set the max grid size and the maximal number of patches
   * inside a grid cube*/
  void setAuxGrid(unsigned int dim, unsigned int max) {
    AUX_GRID_DIM_SKIN = dim;
    MAX_MIXEDCELLS = max;
  }

  /**for the 2d grid set the max grid size and the maximal number of patches
  inside a grid cube. The grid cube itself does not exist just a reference,
  indeed the real quantity is MAX_MIXEDCELLS_2D that is the number of patches
  along the grid tube*/
  void setAuxGrid2D(unsigned int dim, unsigned int max) {
    AUX_GRID_DIM_SKIN_2D = dim;
    MAX_MIXEDCELLS_2D = (max * dim);
  }

  virtual ~SkinSurface();

private:
  /** origin point, direction vector, quadric matrix, intersection parameter
  values. False is return if not intersection is found*/
  // bool rayQuadricIntersection(double*,double*,double**,double*,double*,int
  // thdID,double* cache);
  bool rayQuadricIntersection(double *, double *, double *, double *, double *,
                              int thdID, double *cache);
  /** gives true if the point is inside the list of planes*/
  bool isFeasible(MixedCell *mc, double *point);
  /** project a point to a quadric surface defined by Q*/
  void projectToQuadric(double *y, double *Q, int type, double *proj,
                        double *norm, double &dist);
  /** save quadric patch */
  void saveSkinPatch(ofstream &of, MixedCell *mc, int index,
                     vector<Indexed_weighted_point> &la);
};

static class SkinSurfaceRegister {
  static Surface *createSurface(ConfigFile *conf, DelPhiShared *ds) {
    return new SkinSurface(conf, ds);
  }

public:
  SkinSurfaceRegister() { surfaceFactory().add("skin", createSurface); }
} SkinSurfaceRegisterObject;

// static SurfaceRecorder<SkinSurface> skinRecorder("skin");

#endif
