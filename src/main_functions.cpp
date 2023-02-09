#include <DelphiShared.h>
#include <globals.h>
#include <main_functions.h>
#include <ConfigFile.h>
#include <memory>

#ifdef SPDLOG
#include <spdlog/spdlog.h>
#endif

// init streams, check configuration file for errors and read variables
ConfigFileOP load(std::string confFile, string delimiter, string comment, string sentry, std::string format) {
  spdlog::info("Starting {} {}", PROGNAME, VERSION);

  ConfigFileOP cf = std::make_shared<ConfigFile>();

  // get data from configuration file
  try {
    cf = std::make_shared<ConfigFile>(confFile.c_str(), delimiter, comment, sentry, format);
  } catch (...) {
    spdlog::error("Cannot read {}", confFile);
#ifdef PYTHON
    throw std::invalid_argument("");
#else
    exit(-1);
#endif
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
      spdlog::error("Asked to save epsmap without builiding it");
      spdlog::info("{} Please set Build_epsilon_maps = true", REMARK);
#ifdef PYTHON
      throw std::invalid_argument("");
#else
      exit(-1);
#endif
    }

    if (!conf->buildEpsmaps && conf->projBGP) {
      spdlog::error(
          "Cannot project boundary grid points without an epsilon map.");
      spdlog::info("{} Please set Build_epsilon_maps = true", REMARK);
#ifdef PYTHON
      throw std::invalid_argument("");
#else
      exit(-1);
#endif
    }

    if (!conf->accTri && !conf->buildStatus && conf->tri) {
      // status map is needed to deduced in/out vertices
      spdlog::error(
          "If non analytical triangulation is enabled status map is needed.");
      spdlog::info("{} Please set Build_status_map = true", REMARK);
#ifdef PYTHON
      throw std::invalid_argument("");
#else
      exit(-1);
#endif
    }

    if (conf->fillCavities && !conf->buildStatus) {
      // status map is needed to search cavities
      spdlog::error("If cavity detection is enabled status map is needed.");
      spdlog::info("{} Please set Build_status_map = true", REMARK);
#ifdef PYTHON
      throw std::invalid_argument("");
#else
      exit(-1);
#endif
    }

    if (conf->saveIdebmap && !conf->buildEpsmaps) {
      spdlog::error("Idebmap is computed only if epsilon map is enabled");
      spdlog::info("{} Please set Build_epsilon_maps = true", REMARK);
#ifdef PYTHON
      throw std::invalid_argument("");
#else
      exit(-1);
#endif
    }

    if (!conf->operativeMode.compare("pockets") && !conf->buildStatus) {
      spdlog::warn("Cannot do pocket detection without status map");
      spdlog::info("{} Please set Build_status_map = true", REMARK);

#ifdef PYTHON
      throw std::invalid_argument("");
#else
      exit(-1);
#endif
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
  conf->blobby_B = cf->read<double>( "Blobbyness", -2.5 );
  
  conf->maxSESDim2D  = cf->read<unsigned int>("Max_ses_patches_auxiliary_grid_2d_size", 50);
  conf->maxSESPatches2D = cf->read<unsigned int>("Max_ses_patches_per_auxiliary_grid_2d_cell", 400);
  conf->maxSESDim = cf->read<unsigned int>("Max_ses_patches_auxiliary_grid_size", 100);
  conf->maxSESPatches = cf->read<unsigned int>("Max_ses_patches_per_auxiliary_grid_cell", 400);
  
  conf->mp = cf->read<int>("Max_Probes_Self_Intersections", 100);
  conf->si_perfil = cf->read<double>("Self_Intersections_Grid_Coefficient", 1.5);
  
  conf->radius = cf->read<double>( "Example_Surface_Parameter", 1.0);

  conf->sfname = cf->read<string>( "Surface_File_Name", "mesh.off" );
	conf->maxMeshDim = cf->read<unsigned int>( "Max_mesh_auxiliary_grid_size", 100 );
	conf->maxMeshPatches = cf->read<unsigned int>( "Max_mesh_patches_per_auxiliary_grid_cell", 250 );
	conf->maxMeshDim2D = cf->read<unsigned int>( "Max_mesh_auxiliary_grid_2d_size", 100 );
	conf->maxMeshPatches2D = cf->read<unsigned int>( "Max_mesh_patches_per_auxiliary_grid_2d_cell", 250 );
  conf->NumMSMSfiles = cf->read<int>( "Num_MSMS_files", 1 );

  conf->skin_s = cf->read<double>("Skin_Surface_Parameter", 0.45);
  conf->maxSkinDim = cf->read<unsigned int>("Max_skin_patches_auxiliary_grid_size", 100);
  conf->maxSkinPatches = cf->read<unsigned int>("Max_skin_patches_per_auxiliary_grid_cell", 400);
  conf->maxSkinDim2D = cf->read<unsigned int>("Max_skin_patches_auxiliary_grid_2d_size", 50);
  conf->maxSkinPatches2D = cf->read<unsigned int>( "Max_skin_patches_per_auxiliary_grid_2d_cell", 400);
  conf->useFastProjection = cf->read<bool>("Skin_Fast_Projection", false);
  conf->savePovRay = cf->read<bool>("Save_PovRay", false);
  
  conf->checkDuplicatedVertices = cf->read<bool>( "Check_duplicated_vertices", true );
  conf->wellShaped = cf->read<bool>( "Keep_Water_Shaped_Cavities", false );
  conf->probeRadius = cf->read<double>( "Probe_Radius", 1.4 );
  conf->lb = cf->read<bool>( "Load_Balancing",true);
  conf->vaFlag = cf->read<bool>( "Vertex_Atom_Info",false);
  conf->computeNormals = cf->read<bool>( "Compute_Vertex_Normals",false);
  conf->saveMSMS = cf->read<bool>( "Save_Mesh_MSMS_Format",false);
  conf->sternLayer = cf->read<double>( "Stern_layer", -1. );
  conf->Max_Atoms_Multi_Grid = cf->read<int>( "Max_Atoms_Multi_Grid", 100 );
  conf->surfName = cf->read<string>("Surface");

  return conf;
}

void cite() {
  spdlog::info("If you use NanoShaper please cite these works:");
  spdlog::info(
      "\tS. Decherchi, W. Rocchia, \"A general and Robust Ray-Casting-Based "
      "Algorithm for Triangulating Surfaces at the Nanoscale\"; PlosOne");
  spdlog::info("\tlink: "
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
    spdlog::error("Surface construction failed!");
#ifdef PYTHON
    throw std::exception();
#else
    exit(-1);
#endif
  }

  // Build DelPhi stuff
  surf->getSurf(conf->fillCavities, conf->cavVol);

  double duration = chrono.stop();

  spdlog::info("Surface computation time.. {} [s]", duration);
  spdlog::info("Estimated volume {} [A^3]", surf->getVolume());

  if (conf->debugStatus)
#ifdef SPDLOG
    spdlog::debug("volume {}", surf->getVolume());
#endif // SPDLOG

  if (conf->tri) {
    spdlog::info("Triangulating Surface...");

    Timer chrono;
    chrono.start();

    surf->triangulateSurface();

    if (conf->debugStatus) {
#ifdef SPDLOG
      spdlog::debug("area {}", surf->getArea());
      spdlog::debug("nv {}", surf->getNumVertices());
      spdlog::debug("nt {}", surf->getNumTriangles());
#endif // SPDLOG
    }

    if (conf->smoothing) {
      spdlog::info("Smoothing surface...");
      surf->smoothSurface();
    }

    double duration = chrono.stop();
    spdlog::info("ok!");

    spdlog::info("Total Triangulation time {} [s]", duration);
  }

  if (conf->tri2balls) {
    spdlog::info("Converting triangulation to balls...");
    surf->tri2Balls();
    spdlog::info("ok!");
  }

  if (conf->saveEpsmaps) {
    spdlog::info("Saving epsmaps...");
    // Save epsmap
    dg->saveEpsMaps(refName);
    spdlog::info("ok!");
  }

  if (conf->saveBgps) {
    spdlog::info("Saving bgpmap...");
    dg->saveBGP(refName);
    spdlog::info("ok!");
  }

  if (conf->saveStatusMap) {
    spdlog::info("Saving statusmap and cavities...");
    dg->saveStatus(refName);
    spdlog::info("ok!");
  }

  if (conf->saveCavities) {
    spdlog::info("Saving cavities...");
    dg->saveCavities(false);
    spdlog::info("ok!");
  }

  if (conf->saveIdebmap) {
    spdlog::info("Saving idebmap...");
    dg->saveIdebMap(refName);
    spdlog::info("ok!");
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

  spdlog::info("Step 1 -> fat probe");
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
    spdlog::error("Surface 1 construction failed!");
#ifdef PYTHON
    throw std::exception();
#else
    exit(-1);
#endif
  }

  // fat connolly cancel each cavity (anyway they are smaller than they should
  // be at the end)
  surf1->getSurf(true, INFINITY);

  //////////////////////////////////////////////////////////////////////////////////////////
  // the subsequent surfaces don't need to read atom info

  // Set up Surface 2 (regular probe 1.4)
  // Set up DelPhi grid  2

  spdlog::info("Step 2 -> small probe");

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
    spdlog::error("Surface 2 construction failed!");
#ifdef PYTHON
    throw std::exception();
#else
    exit(-1);
#endif
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
      spdlog::error("Surface 3 construction failed!");
#ifdef PYTHON
      throw std::exception();
#else
      exit(-1);
#endif
    }
    // keep original surface
    surf3->getSurf(false);
  }

  spdlog::info("Step 3 -> differential map");

  spdlog::info("Building pockets by difference map...");
  (*surf1) -= (*surf2);

  ////////////////////// recover split cavities links ///////////////
  if (conf->linkPockets) {
    int nr = 1;
    double duration1 = 0;

    Timer chrono1;
    chrono1.start();

    spdlog::info("Step 3.1 -> linking");

    spdlog::info("Linking cavities/pockets...");
    while (nr != 0) {
      // check cavities links, use the link status map as reference map
      // to check accessibility
      nr = surf1->linkCavities(dg2->status, dg3->status);
      spdlog::info("Merged {} cavities", nr);
    }

    duration1 = chrono1.stop();
    spdlog::info("Diff. Step 4 {} [s]", duration1);

    // delete surf3;
    // delete dg3;
  }
  ///////////////////////////////////////////////////////////////////

  spdlog::info("Step 4 -> filtering, envelope building");

  ///////////////////// volume filter ///////////////////////////////
  surf1->fillCavities(11.4 * conf->numMol, false);
  surf1->getCavitiesAtoms();
  spdlog::info("Saving cavities info..");

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
    spdlog::info("Recovering cavity/pocket distinction..");
    surf2->getCavities();
    // mark which are pockets and which are cavities
    dg1->markPockets(dg2->status, isPocket);
  }

  spdlog::info("Build the envelope of each cavity/pocket");

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
      spdlog::info("Triangulating enveloped cavity/pockets {} ", i);
      cs_temp->triangulateSurface(0.0, tri_);
      areas.push_back(cs_temp->getArea());

      // NanoShaper identifies the set of grid points that are nearest
      // to the given triangulation points and save those points that are
      // completely solvent exposed (i.e. external) in surf2 (the slim surface).
      // Those points represent a point-wise representation of the entrance of
      // the pocket
      if (saveEntranceInfo) {
        if (!isPocket[i]) {
          // spdlog::error( "Cannot find the entrance of a cavity; it
          // must be a pocket.";
          if (!conf->cavAndPockets) {
            spdlog::info("{} You have to enable Cavities and pockets flag to "
                         "do a distinction between a pocket and a cavity",
                         REMARK);
#ifdef PYTHON
            throw std::invalid_argument("");
#else
            exit(-1);
#endif
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
    FILE *fp = fopen(name2, "w");
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
    spdlog::info("------------------------------------");
    spdlog::info("     Pocket Detection Summary       ");
    spdlog::info("------------------------------------");
    spdlog::info("Detected a total of {} pockets/cavities having at least the "
                 "volume of {} water molecules",
                 nc, conf->numMol);
    ;

    for (int i = 0, ii = 0; i < nc; i++) {
      if (conf->cavAndPockets) {
        if (conf->tri) {
          if (isPocket[i]) {
            spdlog::info("Pocket {} vol {} area {} body area {}", i, volumes[i],
                         areas[i], areas2[ii]);
            if (conf->debugStatus) {
#ifdef SPDLOG
              spdlog::debug("pocket_vol {}", volumes[i]);
              spdlog::debug("pocket_area {} pocket_body_area {}", areas[i],
                            areas2[i]);
#endif // SPDLOG
            }
            ii++;
          } else {
            spdlog::info("Cavity {} vol {} area {}", i, volumes[i], areas[i]);
            if (conf->debugStatus) {
#ifdef SPDLOG
              spdlog::debug("cav_vol {}", volumes[i]);
              spdlog::debug("cav_area {}", areas[i]);
#endif // SPDLOG
            }
          }
        } else {
          if (isPocket[i]) {
            spdlog::info("Pocket {} vol {}", i, volumes[i]);
            if (conf->debugStatus) {
#ifdef SPDLOG
              spdlog::debug("pocket_vol {}", volumes[i]);
#endif // SPDLOG
            }
          } else {
            spdlog::info("Cavity {} vol {}", i, volumes[i]);
            if (conf->debugStatus) {
#ifdef SPDLOG
              spdlog::debug("cav_vol {}", volumes[i]);
#endif // SPDLOG
            }
          }
        }
      } else {
        if (conf->tri) {
          spdlog::info("Pocket {} vol {} area {} body area {}", i, volumes[i],
                       areas[i], areas2[ii]);
          if (conf->debugStatus) {
#ifdef SPDLOG
            spdlog::debug("pocket_vol {}", volumes[i]);
            spdlog::debug("pocket_area {}", areas[i]);
#endif // SPDLOG
          }
          ii++;
        } else {
          spdlog::info("Pocket {} vol {}", i, volumes[i]);
          if (conf->debugStatus) {
#ifdef SPDLOG
            spdlog::debug("pocket_vol {}", volumes[i]);
#endif // SPDLOG
          }
        }
      }
    }
  }

  spdlog::info("Pocket detection time.. {} [s]", duration);
  spdlog::info("Cleaning memory...");

  // delete surf1;
  // delete dg1;
  // delete surf2;
  // delete dg2;
  spdlog::info("ok!");
}
