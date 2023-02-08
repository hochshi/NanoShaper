
#include <main_functions.h>
#include <spdlog/spdlog.h>

#ifdef DBGMEM_CRT
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW
#endif

#ifndef DELPHI_BIND
/** @mainpage
        @brief
        NanoShaper is a framework to analyze an arbitrary surface; a particular
   attention is given to molecular surfaces. NanoShaper presents several
   features: <ul> <li> Loading of a triangulated mesh in OFF, PLY or MSMS format
        <li> Analytical build and triangulation of the Skin surface, the Blobby
   surface and the Connolly surface <li> Colour with in/out info of a 3D grid
        <li> Pockets and cavities detection.
        <li> Easy to plug a new surface: it is sufficient to expand Surface
   class, no other modifications are required, just writing the .h and .cpp. The
   surface registration is automatically managed by a template based mechanism.
        <li> It can be compiled as a Standalone, python lib or as a DelPhi
   plug-in.
        </ul>

        History:

        <ul>
        <li> Version 0.3.1:
        <ul>
        <li> Much faster Skin Surface build-up. Acceleration grid 2d also for
   the mesh ray casting.
        </ul>

        <li> Version 0.3.2:
        <ul>
        <li> bug fix on Connolly Surface.
        </ul>

        <li> Version 0.4:
        <ul>
        <li> Using octrees to reduce memory requirements for triangulation and
   using an improved multi-threading strategy. <li> Atom based Multi-dielectric
   is supported (DelPhi). <li> Stern Layer (DelPhi). <li> An efficient,
   ray-casting based algorithm to convert a set of surface samples to a set of
   balls whose envelope approximatess the surface. <li> Atoms in the neighboroud
   of a detected cavity are returned.
        </ul>

        <li> Version 0.5:
        <ul>
        <li> "-=" Filtered Operator to support difference map to detect pockets
   not detectable by cavity detector <li> Operative_Mode keyword. Values =
   normal, difference, membrane <li> Introduced experimental Fast Van Der Waals
   surface. Analytical for FD PDE, not analytical triangulation <li> Introduced
   experimental Coulombic surface <li> From DelPhi atinf and atsurf are
   imported. atsurf is filled  (DelPhi) <li> Fast projection flag to switch
   between fast (approximate) and slow (accurate) SkinSurface bgp projections
        <li> Output for pocket detection in ProShape style added when used with
   DelPhi <li> Improved threads load balancing <li> Every vertex has its own
   reference atoms (OFF+A format defined) <li> Faster flood fill algorithm <li>
   Analytical normals computation or approximation by triangulation <li> OFF+N,
   OFF+N+A, MSMS formats <li> Vertex normals are computed, analytically if
   possible.
        </ul>

        <li> Version 0.6:
        <ul>
        <li> Refactoring using metaprogramming to instantiate surfaces.
                To introduce a new surface, just define that class and register
   it using the registration macro. Now surface constructors access directly the
   configuration file structure. Any new surface thus can define each own
   keywords.
        </ul>

        <li> Version 0.7
        <ul>
        <li> Parallel Marching Cubes: 2 steps MC, cache friendly and trivially
   parallel <li> Dangling vertices now are removed a priori. Dangling vertices
   are those that were sampled, gave an in/out map that is in contrast with
   another ray for numerical reasons. If this happen these vertices, after ray
   tracing, are identified and removed before marching cubes. In this way the
   exact number of MC vertices is known a priori. Another pathological situation
   that may happen on previous versions is that two vertices are false vertices
   on the same edge (let's say an high frequency) now these vertices are
   removed. The previous version of NanoShaper had this added spurious vertices
   in memory; however they are not present in the mesh, so the previous meshes
   are still correct. Moreover now a slightly faster way to color the grid in
   ray tracing is defined. <li> If a normal cannot be computed analytically due
   to any kind of problem, than this is approximated in a subsequent step such
   that any file with normals will have all the normals computed. Most of them
   will be analytical (usually less than 1 over 10000 are not analytical). <li>
   High accuracy timing boost::chrono is optionally available <li> Bug fix in
   blobby surface due to surface factory refactoring <li> Introduced the
   ExampleSurface to let the user better understand how to add a surface <li>
   Using templates instead of macros for memory access/write, check bounds and
   added p=NULL semantic on deletion. from now on vectors of primitive types
   should be allocated, accessed, written, deleted with the given template
   functions <li> Refactoring of main method. <li> Introduced pure virtuals
   init/clear methods to force the Developer to write clean C++ code for the
   Surfaces <li> Optimization on the floodfill. Ping pong to get as much as
   possible cache friendly. <li> Bug fix of ConnollySurface, now Probe_Radius is
   not ignored.
        </ul>

        <li> Version 0.7.4
        <ul>
        <li> Introduced Max_Atoms_Multi_Grid keyword to control atomsmap max
   size. <li> Introduced an algorithm to detected the entrance of a pocket. The
   entrance is saved as a set of points in entranceX.xyz and entranceX.xyzn
   where the first file can be visualized in vmd. <li> If the Connolly surface
   is computed a file name called exposed.xyz saves as Carbons in xyz format
   (VMD) the list of exposed atoms. Additionally a file named exposedIndices.txt
   is the list of 0-based indices of exposed atoms. <li> The area of pockets,
   excluding the mouth is returned. Triangulated files are named
   cav_triX_body.TRI where is triangulation extension and X is the index of the
   pocket <li> Bugfix for MSMS save for blobby

        </ul>

        <li> Version 0.7.5
        <ul>
        <li> Bug fix: Save_Status_map flag was not read. Fixed

        </ul>


        </ul>

        @author Sergio Decherchi
        @date 25-07-2014
        @version 0.7.5
        */
int main(int argc, char *argv[]) {
#ifdef DBGMEM_CRT
  //		 _crtBreakAlloc = 16222;
  _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

  int numargs = argc - 1;
  char confFile[BUFLEN];

  if (!numargs)
    strcpy(confFile, "surfaceConfiguration.prm");
  else
    strcpy(confFile, argv[1]);

  if (numargs > 1)
    cout << endl << WARN << "Ignoring additional non required parameters";

  // check configuration consistency, init error stream, get configuration
  ConfigFileOP cf = load(string(confFile));
  ConfigurationOP conf = parse(cf);

  // operative modes are the NanoShaper lib utilization modes.

  // just build the surface
  if (!conf->operativeMode.compare("normal")) {
    // Set up DelPhi-like environment
    DelPhiShared *dg = new DelPhiShared(conf->scale, conf->perfill,
                                        conf->molFile, conf->buildEpsmaps,
                                        conf->buildStatus, conf->multi_diel);
    // Get surface
    Surface *surf = surfaceFactory().create(cf.get(), dg);
    normalMode(surf, dg, conf);
    cout << endl << INFO_STR << "Cleaning memory...";
    cout.flush();
    delete surf;
    delete dg;
    cout << "ok!";
  }
  // detect pockets
  else if (!conf->operativeMode.compare("pockets")) {
    pocketMode(false, cf, conf);
  } else {
    cout << endl << INFO_STR << "Unknown operative mode";
    cout << endl;
    return -1;
  }

  cite();
  cout << endl << endl;

// Memory leak detection
#ifdef DBGMEM_CRT
  _CrtDumpMemoryLeaks();
#endif

  return 0;
}

#endif

#ifdef DELPHI_BIND

/** This is the DelPhi Fortran-C binding. Returns the number of bgp */
// check visual, borland c++, xlc, watcom, for windows
#if (defined _WIN32) || (defined __WIN32__) || (defined __TOS_WIN__) ||        \
    (defined __WINDOWS__)
extern "C" __declspec(dllexport)
// linux,unix,apple,android
#elif (defined __linux) || (defined __unix) || (defined macintosh) ||          \
    (defined Macintosh) || (defined __APPLE__ && defined __MACH__) ||          \
    (defined __ANDROID__)
extern "C"
#endif

    void epsmakemodule(double xmin, double ymin, double zmin, double xmax,
                       double ymax, double zmax, double c1, double c2,
                       double c3, double rmax, double perf, int *i_epsmap,
                       int igrid, double scale, double *i_scspos,
                       double *i_scsnor, bool *i_idebmap, int *i_ibgp,
                       int inside, int *ibnum, int natom, double *xn1,
                       double *rad, int *d, int maxbgp, double probe,
                       double exrad, char *atinf, int *atsurf) {

  // check configuration consistency, init error stream and get configuration
  ConfigFile *cf = init(string("surfaceConfiguration.prm"));

  // override settings coming from conf with that coming from DelPhi
  if (cf->keyExists("Grid_scale"))
    cf->remove("Grid_scale");
  cf->add<double>("Grid_scale", scale);

  if (cf->keyExists("Grid_perfil"))
    cf->remove("Grid_perfil");
  cf->add<double>("Grid_perfil", perf);

  conf.scale = scale;
  conf.perfill = perf;

  ConfigFile *cf2;
  // get data from custom configuration file
  try {
    cf2 = new ConfigFile("custom.prm");
  } catch (...) {
    cout << endl << ERR << "Cannot read custom.prm";
    cout.flush();
    exit(-1);
  }
  string mode = cf2->read<string>("Surface", "ses");
  delete cf2;

  if (cf->keyExists("Surface"))
    cf->remove("Surface");

  cf->add<string>("Surface", mode);

  if (conf.operativeMode == "normal") {
    cout << endl << INFO_STR << "Binding with DelPhi..";

    // Set up delphi environment
    DelPhiShared *dg = new DelPhiShared();
    dg->DelPhiBinding(xmin, ymin, zmin, xmax, ymax, zmax, c1, c2, c3, rmax,
                      perf, i_epsmap, igrid, scale, i_scspos, i_scsnor,
                      i_idebmap, i_ibgp, maxbgp, conf.buildStatus, atsurf);

    // load atoms info if present
    dg->loadAtoms(natom, xn1, rad, NULL, NULL, atinf);

    // here only means populate epsilon map
    dg->buildEpsmap(true);

    cout << endl << INFO_STR << "DelPhi grid is " << igrid;

    Surface *surf = surfaceFactory().create(cf, dg);
    surf->setProjBGP(true);
    surf->setInsideCode(inside);

    if (exrad > 0) {
      cout << endl << INFO_STR << "Setting Stern Layer -> " << exrad << " Angstrom";
      surf->setSternLayer(exrad);
    }

    // normal protocol
    normalMode(surf, dg);

    // avoid to destroy things that now belong to DelPhi
    dg->finalizeBinding(ibnum);

    cout << endl << INFO_STR << "Cleaning memory...";
    cout.flush();

    delete surf;
    delete dg;

    cout << "ok!";
    cout << endl << endl;
    return;
  }

  // In this mode use the difference on grid operator to infer pockets and
  // cavities (if requested)
  else if (!conf.operativeMode.compare("pockets")) {
    DelPhiShared::parseAtomInfo(atinf, conf.molFile, natom, xn1, rad);

    pocketMode(true, cf);

    cout << endl << INFO_STR << "Brute force exiting to avoid DelPhi solver";
    cout << endl << INFO_STR << "Closing " << PROGNAME << "\n";
    exit(-1);
  } else {
    cout << endl << INFO_STR << "Unknown operative mode: " << conf.operativeMode;
    exit(-1);
  }

  cite();
  cout << endl << endl;
}

#endif
