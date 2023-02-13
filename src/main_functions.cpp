#include <ConfigFile.h>
#include <DelphiShared.h>
#include <globals.h>
#include <logging.h>
#include <main_functions.h>
#include <memory>
#include <stdexcept>

// init streams, check configuration file for errors and read variables
ConfigFileOP load(std::string confFile, string delimiter, string comment,
                  string sentry, std::string format) {
  logging::log<logging::level::info>("Starting {} {}", PROGNAME, VERSION);

  ConfigFileOP cf = std::make_shared<ConfigFile>();

  // get data from configuration file
  try {
    cf = std::make_shared<ConfigFile>(confFile.c_str(), delimiter, comment,
                                      sentry, format);
  } catch (...) {
    logging::log<logging::level::err>("Cannot read {}", confFile);
    throw std::invalid_argument("Cannog read config file");
  }

  return cf;
}

ConfigurationOP parse(ConfigFileOP cf) {

  ConfigurationOP conf = std::make_shared<Configuration>();

  conf->buildEpsmaps = cf->read<bool>("Build_epsilon_maps", false);
  conf->saveEpsmaps = cf->read<bool>("Save_eps_maps", false);

  conf->saveStatusMap = cf->read<bool>("Save_Status_map", false);

  conf->buildStatus = cf->read<bool>("Build_status_map", false);
  conf->saveIdebmap = cf->read<bool>("Save_ideb_map", false);
  conf->fillCavities = cf->read<bool>("Cavity_Detection_Filling", false);

  conf->projBGP = cf->read<bool>("Project_boundary_grid_points", false);

  conf->accTri = cf->read<bool>("Accurate_Triangulation", false);
  conf->tri = cf->read<bool>("Triangulation", false);
  conf->operativeMode = cf->read<string>("Operative_Mode", "normal");

  bool dbg = cf->read<bool>("Debug_Internals", false);
  conf->debug = dbg;
  conf->debugStatus = dbg;

  if (cf != NULL) {
    if (!conf->buildEpsmaps && conf->saveEpsmaps) {
      logging::log<logging::level::err>(
          "Asked to save epsmap without builiding it");
      logging::log<logging::level::info>(
          "{} Please set Build_epsilon_maps = true", REMARK);
      throw std::invalid_argument(
          "Asked to save epsmap without building it! "
          "Please set Build_epsilon_maps to true.");
    }

    if (!conf->buildEpsmaps && conf->projBGP) {
      logging::log<logging::level::err>(
          "Cannot project boundary grid points without an epsilon map.");
      logging::log<logging::level::info>(
          "{} Please set Build_epsilon_maps = true", REMARK);
      throw std::invalid_argument(
          "Cannot project boundary grid points without an epsilon map. Please "
          "set Build_epsilon_maps = true");
    }

    if (!conf->accTri && !conf->buildStatus && conf->tri) {
      // status map is needed to deduced in/out vertices
      logging::log<logging::level::err>(
          "If non analytical triangulation is enabled status map is needed.");
      logging::log<logging::level::info>(
          "{} Please set Build_status_map = true", REMARK);
      throw std::invalid_argument(
          "If non analytical triangulation is enabled status map is needed. "
          "Please set Build_status_map = true");
    }

    if (conf->fillCavities && !conf->buildStatus) {
      // status map is needed to search cavities
      logging::log<logging::level::err>(
          "If cavity detection is enabled status map is needed.");
      logging::log<logging::level::info>(
          "{} Please set Build_status_map = true", REMARK);
      throw std::invalid_argument(
          "If cavity detection is enabled status map is needed. Please set "
          "Build_status_map = true");
    }

    if (conf->saveIdebmap && !conf->buildEpsmaps) {
      logging::log<logging::level::err>(
          "Idebmap is computed only if epsilon map is enabled");
      logging::log<logging::level::info>(
          "{} Please set Build_epsilon_maps = true", REMARK);
      throw std::invalid_argument(
          "Idebmap is computed only if epsilon map is enabled. Please set "
          "Build_epsilon_maps = true");
    }

    if (!conf->operativeMode.compare("pockets") && !conf->buildStatus) {
      logging::log<logging::level::warn>(
          "Cannot do pocket detection without status map");
      logging::log<logging::level::info>(
          "{} Please set Build_status_map = true", REMARK);
      throw std::invalid_argument(
          "Cannot do pocket detection without status "
          "map. Please set Build_status_map = true");
    }
  }

  conf->cavVol = cf->read<double>("Conditional_Volume_Filling_Value", 11.4);
  conf->numMol = cf->read<int>("Num_Wat_Pocket", 2);

  // grid (DelPhi) params
  conf->scale = cf->read<double>("Grid_scale", 2.0);
  conf->perfill = cf->read<double>("Grid_perfil", 80.0);
  conf->molFile = cf->read<string>("XYZR_FileName", "temp.xyzr");
  conf->multi_diel = cf->read<bool>("Multi_Dielectric", false);

  // tri
  conf->smoothing = cf->read<bool>("Smooth_Mesh", false);
  conf->tri2balls = cf->read<bool>("Tri2Balls", false);

  // save data
  conf->saveEpsmaps = cf->read<bool>("Save_eps_maps", false);
  conf->saveBgps = cf->read<bool>("Save_bgps", false);
  conf->saveCavities = cf->read<bool>("Save_Cavities", false);

  // globals
  conf->sysName = cf->read<string>("Sys_Name", "mol");
  conf->numthd = cf->read<int>("Number_thread", -1);
  conf->printAvailSurf = cf->read<bool>("Print_Available_Surfaces", false);
  conf->currentSeed = cf->read<int>("Seed", 1);

  // pocket detection
  conf->cavAndPockets = cf->read<bool>("Pockets_And_Cavities", true);
  conf->linkPockets = cf->read<bool>("Link_Pockets", false);
  conf->pocketRadiusBig = cf->read<double>("Pocket_Radius_Big", 3.0);
  conf->pocketRadiusSmall = cf->read<double>("Pocket_Radius_Small", 1.4);
  conf->pocketRadiusLink = cf->read<double>("Pocket_Radius_Link", 1.0);

  // Rest of the values came from the code
  conf->blobby_B = cf->read<double>("Blobbyness", -2.5);

  conf->maxSESDim2D =
      cf->read<unsigned int>("Max_ses_patches_auxiliary_grid_2d_size", 50);
  conf->maxSESPatches2D =
      cf->read<unsigned int>("Max_ses_patches_per_auxiliary_grid_2d_cell", 400);
  conf->maxSESDim =
      cf->read<unsigned int>("Max_ses_patches_auxiliary_grid_size", 100);
  conf->maxSESPatches =
      cf->read<unsigned int>("Max_ses_patches_per_auxiliary_grid_cell", 400);

  conf->mp = cf->read<int>("Max_Probes_Self_Intersections", 100);
  conf->si_perfil =
      cf->read<double>("Self_Intersections_Grid_Coefficient", 1.5);

  conf->radius = cf->read<double>("Example_Surface_Parameter", 1.0);

  conf->sfname = cf->read<string>("Surface_File_Name", "mesh.off");
  conf->maxMeshDim =
      cf->read<unsigned int>("Max_mesh_auxiliary_grid_size", 100);
  conf->maxMeshPatches =
      cf->read<unsigned int>("Max_mesh_patches_per_auxiliary_grid_cell", 250);
  conf->maxMeshDim2D =
      cf->read<unsigned int>("Max_mesh_auxiliary_grid_2d_size", 100);
  conf->maxMeshPatches2D = cf->read<unsigned int>(
      "Max_mesh_patches_per_auxiliary_grid_2d_cell", 250);
  conf->NumMSMSfiles = cf->read<int>("Num_MSMS_files", 1);

  conf->skin_s = cf->read<double>("Skin_Surface_Parameter", 0.45);
  conf->maxSkinDim =
      cf->read<unsigned int>("Max_skin_patches_auxiliary_grid_size", 100);
  conf->maxSkinPatches =
      cf->read<unsigned int>("Max_skin_patches_per_auxiliary_grid_cell", 400);
  conf->maxSkinDim2D =
      cf->read<unsigned int>("Max_skin_patches_auxiliary_grid_2d_size", 50);
  conf->maxSkinPatches2D = cf->read<unsigned int>(
      "Max_skin_patches_per_auxiliary_grid_2d_cell", 400);
  conf->useFastProjection = cf->read<bool>("Skin_Fast_Projection", false);
  conf->savePovRay = cf->read<bool>("Save_PovRay", false);

  conf->checkDuplicatedVertices =
      cf->read<bool>("Check_duplicated_vertices", true);
  conf->wellShaped = cf->read<bool>("Keep_Water_Shaped_Cavities", false);
  conf->probeRadius = cf->read<double>("Probe_Radius", 1.4);
  conf->lb = cf->read<bool>("Load_Balancing", true);
  conf->vaFlag = cf->read<bool>("Vertex_Atom_Info", false);
  conf->computeNormals = cf->read<bool>("Compute_Vertex_Normals", false);
  conf->saveMSMS = cf->read<bool>("Save_Mesh_MSMS_Format", false);
  conf->sternLayer = cf->read<double>("Stern_layer", -1.);
  conf->Max_Atoms_Multi_Grid = cf->read<int>("Max_Atoms_Multi_Grid", 100);
  conf->surfName = cf->read<string>("Surface");

  return conf;
}

void cite() {
  logging::log<logging::level::info>(
      "If you use NanoShaper please cite these works:");
  logging::log<logging::level::info>(
      "\tS. Decherchi, W. Rocchia, \"A general and Robust Ray-Casting-Based "
      "Algorithm for Triangulating Surfaces at the Nanoscale\"; PlosOne");
  logging::log<logging::level::info>(
      "\tlink: "
      "http://www.plosone.org/article/metrics/"
      "info%3Adoi%2F10.1371%2Fjournal.pone.0059744");
}

/** the set of operations in the usual mode of usage. This function is not
responsible for Surface or grid memory. The caller is the responsible.*/
void normalMode(SurfaceOP surf, DelPhiSharedOP dg, ConfigurationOP conf) {
  if (conf->printAvailSurf)
    SurfaceFactory::getInstance().print();

  Timer chrono;
  chrono.start();

  char refName[100];
  strcpy(refName, conf->sysName.c_str());

  // Pre-process surface
  bool outsurf = surf->build();

  if (!outsurf) {
    logging::log<logging::level::err>("Surface construction failed!");
    throw std::runtime_error("Surface construction failed!");
  }

  // Build DelPhi stuff
  surf->getSurf(conf->fillCavities, conf->cavVol);

  double duration = chrono.stop();

  logging::log<logging::level::info>("Surface computation time.. {} [s]",
                                     duration);
  logging::log<logging::level::info>("Estimated volume {} [A^3]",
                                     surf->getVolume());

  if (conf->debugStatus)
    logging::log<logging::level::debug>("volume {}", surf->getVolume());

  if (conf->tri) {
    logging::log<logging::level::info>("Triangulating Surface...");

    Timer chrono;
    chrono.start();

    surf->triangulateSurface();

    if (conf->debugStatus) {
      logging::log<logging::level::debug>("area {}", surf->getArea());
      logging::log<logging::level::debug>("nv {}", surf->getNumVertices());
      logging::log<logging::level::debug>("nt {}", surf->getNumTriangles());
    }

    if (conf->smoothing) {
      logging::log<logging::level::info>("Smoothing surface...");
      surf->smoothSurface();
    }

    double duration = chrono.stop();
    logging::log<logging::level::info>("ok!");

    logging::log<logging::level::info>("Total Triangulation time {} [s]",
                                       duration);
  }

  if (conf->tri2balls) {
    logging::log<logging::level::info>("Converting triangulation to balls...");
    surf->tri2Balls();
    logging::log<logging::level::info>("ok!");
  }

  if (conf->saveEpsmaps) {
    logging::log<logging::level::info>("Saving epsmaps...");
    // Save epsmap
    dg->saveEpsMaps(refName);
    logging::log<logging::level::info>("ok!");
  }

  if (conf->saveBgps) {
    logging::log<logging::level::info>("Saving bgpmap...");
    dg->saveBGP(refName);
    logging::log<logging::level::info>("ok!");
  }

  if (conf->saveStatusMap) {
    logging::log<logging::level::info>("Saving statusmap and cavities...");
    dg->saveStatus(refName);
    logging::log<logging::level::info>("ok!");
  }

  if (conf->saveCavities) {
    logging::log<logging::level::info>("Saving cavities...");
    dg->saveCavities(false);
    logging::log<logging::level::info>("ok!");
  }

  if (conf->saveIdebmap) {
    logging::log<logging::level::info>("Saving idebmap...");
    dg->saveIdebMap(refName);
    logging::log<logging::level::info>("ok!");
  }
}

void pocketMode(bool hasAtomInfo, ConfigurationOP cf, ConfigurationOP conf) {
  double duration = 0;

  bool localEpsMap = false;
  bool localStatusMap = true;
  bool localMulti = false;

  // only debug pockets/cavities
  if (conf->debug) {
    conf->stopDebug();
  }

  Timer chrono;
  chrono.start();

  SurfaceOP surf1, surf2, surf3;
  double areaSurf = 0, volSurf = 0;

  char mol[100] = "temp.txt", refName[100] = "mol";

  /////////////////////////////////////////////////////////////////////////////////////////
  // Set up Surface 1 (fat probe)

  logging::log<logging::level::info>("Step 1 -> fat probe");
  DelPhiSharedOP dg1 = std::make_shared<DelPhiShared>(
      conf->scale, conf->perfill, conf->molFile, localEpsMap, localStatusMap,
      localMulti, hasAtomInfo);
  int natom = dg1->getNumAtoms();
  // cf->remove(string("Surface"));
  // cf->add<string>("Surface", "ses");
  cf->surfName = "ses";
  surf1 = SurfaceFactory::getInstance().create(cf, dg1);
  surf1->setInsideCode(5);
  surf1->setProbeRadius(conf->pocketRadiusBig);
  surf1->setProjBGP(false);
  surf1->setKeepWellShapedCavities(false);

  // Pre-process surface
  bool outsurf = surf1->build();

  if (!outsurf) {
    logging::log<logging::level::err>("Surface 1 construction failed!");
    throw std::runtime_error("Surface 1 construction failed!");
  }

  // fat connolly cancel each cavity (anyway they are smaller than they should
  // be at the end)
  surf1->getSurf(true, INFINITY);

  //////////////////////////////////////////////////////////////////////////////////////////
  // the subsequent surfaces don't need to read atom info

  // Set up Surface 2 (regular probe 1.4)
  // Set up DelPhi grid  2

  logging::log<logging::level::info>("Step 2 -> small probe");

  DelPhiSharedOP dg2 = std::make_shared<DelPhiShared>(
      conf->scale, conf->perfill, conf->molFile, localEpsMap, localStatusMap,
      localMulti, false);
  surf2 = SurfaceFactory::getInstance().create(cf, dg2);
  surf2->setInsideCode(10);
  surf2->setProjBGP(false);
  surf2->setProbeRadius(conf->pocketRadiusSmall);
  surf2->setKeepWellShapedCavities(false);

  /*
  //// to check percolation ///////////////////////
  DelPhiShared* dg2 = new
  DelPhiShared(conf.scale,conf.perfill,conf.mol,localStatusMap,localMulti,false);
  surf2 = SurfaceFactory::getInstance().create(cf,dg2);
  surf2->inside = 10;
  ////////////////////////////////////////////////
  */

  // Pre-process surface
  outsurf = surf2->build();

  if (!outsurf) {
    logging::log<logging::level::err>("Surface 2 construction failed!");
    throw std::runtime_error("Surface 2 construction failed!");
  }

  // if cav and pockets together -> do not perform cavity detection (keep all)
  // if only pockets then perform cavity detection and remove all cavities
  // by default keep both cavities and pockets
  surf2->getSurf(!conf->cavAndPockets, INFINITY);

  // if triangulation is enabled we triangulate the small probe surface, get
  // surface and area
  if (conf->tri) {
    surf2->triangulateSurface();
    surf2->smoothSurface();
    areaSurf = surf2->getArea();
  }

  volSurf = surf2->getVolume();

  // to check percolation ///////////
  /*
        surf1->difference(surf2);
        surf1->triangulateSurface(0.0,"diffmap.off",true);
        #ifdef PYTHON
throw std::exception("");
#else
exit(-1);
#endif
  */
  //////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////
  DelPhiSharedOP dg3;
  if (conf->linkPockets) {
    // Set up Surface 3 (accessibility probe)
    dg3 = std::make_shared<DelPhiShared>(conf->scale, conf->perfill,
                                         conf->molFile, localEpsMap,
                                         localStatusMap, localMulti, false);
    surf3 = SurfaceFactory::getInstance().create(cf, dg3);
    surf2->setProjBGP(false);
    surf2->setProbeRadius(conf->pocketRadiusSmall);
    surf2->setKeepWellShapedCavities(false);
    surf3->setInsideCode(15);

    // Pre-process surface
    outsurf = surf3->build();

    if (!outsurf) {
      logging::log<logging::level::err>("Surface 3 construction failed!");
      throw std::runtime_error("Surface 3 construction failed!");
    }
    // keep original surface
    surf3->getSurf(false);
  }

  logging::log<logging::level::info>("Step 3 -> differential map");

  logging::log<logging::level::info>("Building pockets by difference map...");
  (*surf1) -= (*surf2);

  ////////////////////// recover split cavities links ///////////////
  if (conf->linkPockets) {
    int nr = 1;
    double duration1 = 0;

    Timer chrono1;
    chrono1.start();

    logging::log<logging::level::info>("Step 3.1 -> linking");

    logging::log<logging::level::info>("Linking cavities/pockets...");
    while (nr != 0) {
      // check cavities links, use the link status map as reference map
      // to check accessibility
      nr = surf1->linkCavities(dg2->status, dg3->status);
      logging::log<logging::level::info>("Merged {} cavities", nr);
    }

    duration1 = chrono1.stop();
    logging::log<logging::level::info>("Diff. Step 4 {} [s]", duration1);

    // delete surf3;
    // delete dg3;
  }
  ///////////////////////////////////////////////////////////////////

  logging::log<logging::level::info>("Step 4 -> filtering, envelope building");

  ///////////////////// volume filter ///////////////////////////////
  surf1->fillCavities(11.4 * conf->numMol, false);
  surf1->getCavitiesAtoms();
  logging::log<logging::level::info>("Saving cavities info..");

  if (hasAtomInfo) {
    // save cavities in ProShape format
    dg1->saveCavities2(true, conf->sysName);
  } else {
    // save cavities in NanoShaper format
    bool saveOnlyNonFilled = true;
    dg1->saveCavities(saveOnlyNonFilled);
  }

  ///////////////////// mark support atoms //////////////////////////
  // perform again filtering only to mark the support atoms of each cavity
  surf1->filterCavities(true);

  // write cavities as atoms on files to subsequent triangulation
  int nc = dg1->cavitiesToAtoms(1.4);

  // areas of pockets including the pocket closure
  vector<double> areas;
  // areas of pockets upon removal of the pocket closure surface
  vector<double> areas2;
  vector<double> volumes;
  vector<bool> isPocket;

  if (conf->cavAndPockets) {

    // do cavity detection to recover which is cavity in the reference
    logging::log<logging::level::info>(
        "Recovering cavity/pocket distinction..");
    surf2->getCavities();
    // mark which are pockets and which are cavities
    dg1->markPockets(dg2->status, isPocket);
  }

  logging::log<logging::level::info>(
      "Build the envelope of each cavity/pocket");

  // cf->remove("Surface");
  // cf->add<string>("Surface", "skin");
  cf->surfName = "skin";

  // TODO make it an option
  bool saveEntranceInfo = true;

  for (int i = 0; i < nc; i++) {
    char mol[100], tri_[100];
    snprintf(mol, sizeof(mol), "cav%d.txt", i);
    snprintf(tri_, sizeof(tri_), "cav_tri%d", i);
    DelPhiSharedOP dg_temp =
        std::make_shared<DelPhiShared>(2.0, 50, mol, false, false, false);
    // reset seed
    srand(conf->currentSeed);
    SurfaceOP cs_temp = SurfaceFactory::getInstance().create(cf, dg_temp);
    // relatively big such that non degenerate configurations are avoided
    // due to the high packing of atoms. this has a minimal impact
    // on the estimation of the pocket/surface volume
    cs_temp->setRandDisplacement(0.1);
    cs_temp->setTriangulationFlag(true);
    cs_temp->setCheckDuplicatedVertices(false);
    cs_temp->build();
    cs_temp->getSurf(false);
    volumes.push_back(cs_temp->getVolume());
    if (conf->tri) {
      logging::log<logging::level::info>(
          "Triangulating enveloped cavity/pockets {} ", i);
      cs_temp->triangulateSurface(0.0, tri_);
      areas.push_back(cs_temp->getArea());

      // NanoShaper identifies the set of grid points that are nearest
      // to the given triangulation points and save those points that are
      // completely solvent exposed (i.e. external) in surf2 (the slim surface).
      // Those points represent a point-wise representation of the entrance of
      // the pocket
      if (saveEntranceInfo) {
        if (!isPocket[i]) {
          // logging::log<logging::level::err>( "Cannot find the entrance of a
          // cavity; it must be a pocket.";
          if (!conf->cavAndPockets) {
            logging::log<logging::level::info>(
                "{} You have to enable Cavities and pockets flag to "
                "do a distinction between a pocket and a cavity",
                REMARK);
            throw std::invalid_argument(
                " You have to enalbe Cavities and pockets flag to do a "
                "distinction between a pocket and a cavity.");
          }
        } else {
          // is a pocket, we can identify the entrance
          // true means that it is part of the entrance
          vector<bool> results;
          cs_temp->triangulationPointsAreCompletelyOut(surf2.get(), results);

          // the points that are completely out are those that are represent the
          // entrance together with them we must save the normals

          // surf2 save those normals and and save those vertices selectively of
          // this surface
          char buff[100], buff2[100];
          snprintf(buff, sizeof(buff), "entrance%d.xyz", i);
          snprintf(buff2, sizeof(buff2), "entrance%d.xyzn", i);
          cs_temp->saveSelectedPoints(results, buff, buff2);

          char triSubset[100];
          snprintf(triSubset, sizeof(triSubset), "cav_tri_body%d", i);

          // save triangulation and estimate area removing entrance
          // triangulation points
          vector<bool> results2;
          for (unsigned int i = 0; i < results.size(); i++) {
            results2.push_back(!results[i]);
          }

          bool revert = false;
          double reducedArea =
              cs_temp->saveTriSubSet(triSubset, results2, revert);
          areas2.push_back(reducedArea);
        }
      }
    }
    // delete cs_temp;
    // delete dg_temp;
  }

  duration = chrono.stop();

  if (conf->debug) {
    conf->restartDebug();
  }

  if (hasAtomInfo) {
    char name2[100];
    snprintf(name2, sizeof(name2), "%s.info", conf->sysName.c_str());
    FILE* fp = fopen(name2, "w");
    fprintf(fp, "\n Protein             : ");
    fprintf(fp, "\n Probe radius        : %.2f", conf->pocketRadiusSmall);
    fprintf(fp, "\n Number of atoms : %d", natom);
    fprintf(fp, "\n Total Surface Area  : %.4f", areaSurf);
    fprintf(fp, "\n Total Volume        : %.4f", volSurf);
    fprintf(fp, "\n\n");
    fprintf(fp, "\n Pockets :");
    fprintf(fp, "\n");
    fprintf(fp, "\n Id\t\tN_mth\tSurface\t\t\tVolume\n");

    double sumA = 0, sumV = 0;
    for (int i = 0; i < nc; i++) {
      if (conf->tri) {
        if (conf->cavAndPockets) {
          if (isPocket[i] == true) {
            fprintf(fp, " %d\t\t>0\t\t%.4f\t\t%.4f\n", i + 1, areas[i],
                    volumes[i]);
          } else {
            fprintf(fp, " %d\t\t0\t\t%.4f\t\t%.4f\n", i + 1, areas[i],
                    volumes[i]);
          }
        } else {
          fprintf(fp, " %d\t\t>0\t\t%.4f\t\t%.4f\n", i + 1, areas[i],
                  volumes[i]);
        }

        sumA += areas[i];
        sumV += volumes[i];
      } else {
        if (conf->cavAndPockets) {
          if (isPocket[i] == true) {
            fprintf(fp, "\t%d\t>0\t\t%.4f\t\t%.4f\n", i + 1, 0.0, volumes[i]);
          } else {
            fprintf(fp, "\t%d\t0\t\t%.4f\t\t%.4f\n", i + 1, 0.0, volumes[i]);
          }
        } else {
          fprintf(fp, " %d\t\t>0\t\t%.4f\t\t%.4f\n", i + 1, areas[i],
                  volumes[i]);
        }

        sumV += volumes[i];
      }
    }
    fprintf(fp, "\nTot\t\t%.4f\t\t%.4f", sumA, sumV);
    fclose(fp);
  } else {
    logging::log<logging::level::info>("------------------------------------");
    logging::log<logging::level::info>("     Pocket Detection Summary       ");
    logging::log<logging::level::info>("------------------------------------");
    logging::log<logging::level::info>(
        "Detected a total of {} pockets/cavities having at least the "
        "volume of {} water molecules",
        nc, conf->numMol);
    ;

    for (int i = 0, ii = 0; i < nc; i++) {
      if (conf->cavAndPockets) {
        if (conf->tri) {
          if (isPocket[i]) {
            logging::log<logging::level::info>(
                "Pocket {} vol {} area {} body area {}", i, volumes[i],
                areas[i], areas2[ii]);
            if (conf->debugStatus) {
              logging::log<logging::level::debug>("pocket_vol {}", volumes[i]);
              logging::log<logging::level::debug>(
                  "pocket_area {} pocket_body_area {}", areas[i], areas2[i]);
            }
            ii++;
          } else {
            logging::log<logging::level::info>("Cavity {} vol {} area {}", i,
                                               volumes[i], areas[i]);
            if (conf->debugStatus) {
              logging::log<logging::level::debug>("cav_vol {}", volumes[i]);
              logging::log<logging::level::debug>("cav_area {}", areas[i]);
            }
          }
        } else {
          if (isPocket[i]) {
            logging::log<logging::level::info>("Pocket {} vol {}", i,
                                               volumes[i]);
            if (conf->debugStatus) {
              logging::log<logging::level::debug>("pocket_vol {}", volumes[i]);
            }
          } else {
            logging::log<logging::level::info>("Cavity {} vol {}", i,
                                               volumes[i]);
            if (conf->debugStatus) {
              logging::log<logging::level::debug>("cav_vol {}", volumes[i]);
            }
          }
        }
      } else {
        if (conf->tri) {
          logging::log<logging::level::info>(
              "Pocket {} vol {} area {} body area {}", i, volumes[i], areas[i],
              areas2[ii]);
          if (conf->debugStatus) {
            logging::log<logging::level::debug>("pocket_vol {}", volumes[i]);
            logging::log<logging::level::debug>("pocket_area {}", areas[i]);
          }
          ii++;
        } else {
          logging::log<logging::level::info>("Pocket {} vol {}", i, volumes[i]);
          if (conf->debugStatus) {
            logging::log<logging::level::debug>("pocket_vol {}", volumes[i]);
          }
        }
      }
    }
  }

  logging::log<logging::level::info>("Pocket detection time.. {} [s]",
                                     duration);
  logging::log<logging::level::info>("Cleaning memory...");

  // delete surf1;
  // delete dg1;
  // delete surf2;
  // delete dg2;
  logging::log<logging::level::info>("ok!");
}
