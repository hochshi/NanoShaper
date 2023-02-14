// Comment to check cpplinter
#include <DelphiShared.h>
#include <Surface.h>
#include <SurfaceFactory.h>
#include <logging.h>
#include <main_functions.h>
#include <memory>

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
using namespace nanoshaper;

int main(int argc, char* argv[]) {
#ifdef DBGMEM_CRT
  //		 _crtBreakAlloc = 16222;
  _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

  int retVal = 0;
  int numargs = argc - 1;
  char confFile[BUFLEN];

  if (!numargs)
    strcpy(confFile, "surfaceConfiguration.prm");
  else
    strcpy(confFile, argv[1]);

  if (numargs > 1)
    logging::log<logging::level::warn>(
        "Ignoring additional non required parameters");

  // check configuration consistency, init error stream, get configuration
  try {
    ConfigFileOP cf = load(string(confFile));
    ConfigurationOP conf = parse(cf);

    // operative modes are the NanoShaper lib utilization modes.

    // just build the surface
    if (!conf->operativeMode.compare("normal")) {
      // Set up DelPhi-like environment
      // DelPhiShared *dg = new DelPhiShared(conf->scale, conf->perfill,
      //                                     conf->molFile, conf->buildEpsmaps,
      //                                     conf->buildStatus,
      //                                     conf->multi_diel);
      DelPhiSharedOP dg = std::make_shared<DelPhiShared>(
          conf->scale, conf->perfill, conf->molFile, conf->buildEpsmaps,
          conf->buildStatus, conf->multi_diel);
      // Get surface
      SurfaceOP surf = SurfaceFactory::getInstance().create(conf, dg);
      if (surf != nullptr) {
        normalMode(surf, dg, conf);
        surf->saveMesh("triangulatedSurf");
      }
      logging::log<logging::level::info>("Cleaning memory...");

      // delete surf;
      // delete dg;
      logging::log<logging::level::info>("ok!");
    }
    // detect pockets
    else if (!conf->operativeMode.compare("pockets")) {
      pocketMode(false, conf, conf);
    } else {
      logging::log<logging::level::info>("Unknown operative mode");
      retVal = -1;
    }
  } catch (const std::exception& e) {
    retVal = -1;
  }

  cite();

// Memory leak detection
#ifdef DBGMEM_CRT
  _CrtDumpMemoryLeaks();
#endif

  return retVal;
}

#endif
