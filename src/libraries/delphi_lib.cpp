// #ifdef DELPHI_BIND
#include <DelphiShared.h>
#include <SurfaceFactory.h>
#include <logging.h>
#include <main_functions.h>
#include <memory>

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
  ConfigFileOP cf = load(string("surfaceConfiguration.prm"));
  ConfigurationOP conf = parse(cf);

  // override settings coming from conf with that coming from DelPhi
  // if (cf->keyExists("Grid_scale"))
  //   cf->remove("Grid_scale");
  // cf->add<double>("Grid_scale", scale);
  // conf->scale = scale;

  // if (cf->keyExists("Grid_perfil"))
  //   cf->remove("Grid_perfil");
  // cf->add<double>("Grid_perfil", perf);
  // conf->perfill = perf;

  conf->scale = scale;
  conf->perfill = perf;

  ConfigFile *cf2;
  // get data from custom configuration file
  try {
    cf2 = new ConfigFile("custom.prm");
  } catch (...) {
    logging::log<logging::level::err>("Cannot read custom.prm");
    exit(-1);
  }
  string mode = cf2->read<string>("Surface", "ses");
  delete cf2;

  if (cf->keyExists("Surface"))
    cf->remove("Surface");

  cf->add<string>("Surface", mode);
  conf->surfName = mode;

  if (conf->operativeMode == "normal") {
    logging::log<logging::level::info>("Binding with DelPhi..");

    // Set up delphi environment
    // DelPhiShared *dg = new DelPhiShared();
    DelPhiSharedOP dg = std::make_shared<DelPhiShared>();
    dg->DelPhiBinding(xmin, ymin, zmin, xmax, ymax, zmax, c1, c2, c3, rmax,
                      perf, i_epsmap, igrid, scale, i_scspos, i_scsnor,
                      i_idebmap, i_ibgp, maxbgp, conf->buildStatus, atsurf);

    // load atoms info if present
    dg->loadAtoms(natom, xn1, rad, NULL, NULL, atinf);

    // here only means populate epsilon map
    dg->buildEpsmap(true);

    logging::log<logging::level::info>("DelPhi grid is {}", igrid);

    SurfaceOP surf = SurfaceFactory::getInstance().create(conf, dg);
    surf->setProjBGP(true);
    surf->setInsideCode(inside);

    if (exrad > 0) {
      logging::log<logging::level::info>("Setting Stern Layer -> {} Angstrom",
                                         exrad);
      surf->setSternLayer(exrad);
    }

    // normal protocol
    normalMode(surf, dg, conf);

    // avoid to destroy things that now belong to DelPhi
    dg->finalizeBinding(ibnum);

    logging::log<logging::level::info>("Cleaning memory...");

    // delete surf;
    // delete dg;

    logging::log<logging::level::info>("ok!");
    return;
  }

  // In this mode use the difference on grid operator to infer pockets and
  // cavities (if requested)
  else if (!conf->operativeMode.compare("pockets")) {
    DelPhiShared::parseAtomInfo(atinf, conf->molFile, natom, xn1, rad);

    pocketMode(true, conf, conf);

    logging::log<logging::level::info>(
        "Brute force exiting to avoid DelPhi solver");
    logging::log<logging::level::info>("Closing {}", PROGNAME);
    exit(-1);
  } else {
    logging::log<logging::level::info>("Unknown operative mode: {}",
                                       conf->operativeMode);
    exit(-1);
  }

  cite();
}

// #endif
