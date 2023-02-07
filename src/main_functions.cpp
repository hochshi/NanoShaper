#include "globals.h"
#include <main_functions.h>

// init streams, check configuration file for errors and read variables
ConfigFileOP load(std::string confFile) {
  cout << endl << endl << INFO_STR << "Starting " << PROGNAME << " " << VERSION;

  ConfigFileOP cf = std::make_shared<ConfigFile>();

  // get data from configuration file
  try {
    cf = std::make_shared<ConfigFile>(confFile.c_str());
  } catch (...) {
    cout << endl << ERR << "Cannot read " << confFile;
    exit(-1);
  }

  return cf;
}

ConfigurationOP parse(ConfigFileOP cf) {

  ConfigurationOP conf = std::make_shared<Configuration>();

  conf->errorStream = new fstream("stderror.txt", fstream::out);

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

  if (dbg)
    conf->internals = new fstream("internals.txt", fstream::out);
  else
    conf->internals = NULL;

  if (cf != NULL) {
    if (!conf->buildEpsmaps && conf->saveEpsmaps) {
      cout << endl << ERR << "Asked to save epsmap without builiding it";
      cout << endl << REMARK << "Please set Build_epsilon_maps = true";
      cout << endl;
      exit(-1);
    }

    if (!conf->buildEpsmaps && conf->projBGP) {
      cout << endl
           << ERR
           << "Cannot project boundary grid points without an epsilon map.";
      cout << endl << REMARK << "Please set Build_epsilon_maps = true";
      cout << endl;
      exit(-1);
    }

    if (!conf->accTri && !conf->buildStatus && conf->tri) {
      // status map is needed to deduced in/out vertices
      cout
          << endl
          << ERR
          << "If non analytical triangulation is enabled status map is needed.";
      cout << endl << REMARK << "Please set Build_status_map = true";
      cout << endl;
      exit(-1);
    }

    if (conf->fillCavities && !conf->buildStatus) {
      // status map is needed to search cavities
      cout << endl
           << ERR << "If cavity detection is enabled status map is needed.";
      cout << endl << REMARK << "Please set Build_status_map = true";
      cout << endl;
      exit(-1);
    }

    if (conf->saveIdebmap && !conf->buildEpsmaps) {
      cout << endl
           << ERR << "Idebmap is computed only if epsilon map is enabled";
      cout << endl << REMARK << "Please set Build_epsilon_maps = true";
      cout << endl;
      exit(-1);
    }

    if (!conf->operativeMode.compare("pockets") && !conf->buildStatus) {
      cout << endl << WARN << "Cannot do pocket detection without status map";
      cout << endl << REMARK << "Please set Build_status_map = true";
      cout << endl;
      exit(-1);
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

  return conf;
}

void cite() {
  cout << endl;
  cout << endl << INFO_STR << "If you use NanoShaper please cite these works:";
  cout
      << endl
      << CITE
      << "\tS. Decherchi, W. Rocchia, \"A general and Robust Ray-Casting-Based "
         "Algorithm for Triangulating Surfaces at the Nanoscale\"; PlosOne";
  cout << endl
       << CITE
       << "\tlink: "
          "http://www.plosone.org/article/metrics/"
          "info%3Adoi%2F10.1371%2Fjournal.pone.0059744";
  cout << endl;
}

/** the set of operations in the usual mode of usage. This function is not
responsible for Surface or grid memory. The caller is the responsible.*/
void normalMode(Surface *surf, DelPhiShared *dg, ConfigurationOP conf) {
  if (conf->printAvailSurf)
    surfaceFactory().print();

  Timer chrono;
  chrono.start();

  char refName[100];
  strcpy(refName, conf->sysName.c_str());

  // Pre-process surface
  bool outsurf = surf->build();

  if (!outsurf) {
    cout << endl << ERR << "Surface construction failed!";
    exit(-1);
  }

  // Build DelPhi stuff
  surf->getSurf(conf->fillCavities, conf->cavVol);

  double duration = chrono.stop();

  cout << endl << INFO_STR << "Surface computation time.. " << duration << " [s]";
  cout << endl << INFO_STR << "Estimated volume " << surf->getVolume() << " [A^3]";

  if (conf->debugStatus)
    (*(conf->internals)) << endl << "volume " << surf->getVolume();

  if (conf->tri) {
    cout << endl << INFO_STR << "Triangulating Surface...";

    Timer chrono;
    chrono.start();

    surf->triangulateSurface();

    if (conf->debugStatus) {
      (*(conf->internals)) << endl << "area " << surf->getArea();
      (*(conf->internals)) << endl << "nv " << surf->getNumVertices();
      (*(conf->internals)) << endl << "nt " << surf->getNumTriangles();
    }

    if (conf->smoothing) {
      cout << endl << INFO_STR << "Smoothing surface...";
      surf->smoothSurface();
    }

    double duration = chrono.stop();
    cout << "ok!";

    cout << endl << INFO_STR << "Total Triangulation time " << duration << " [s]";
  }

  if (conf->tri2balls) {
    cout << endl << INFO_STR << "Converting triangulation to balls...";
    surf->tri2Balls();
    cout << "ok!";
  }

  if (conf->saveEpsmaps) {
    cout << endl << INFO_STR << "Saving epsmaps...";
    // Save epsmap
    dg->saveEpsMaps(refName);
    cout << "ok!";
  }

  if (conf->saveBgps) {
    cout << endl << INFO_STR << "Saving bgpmap...";
    dg->saveBGP(refName);
    cout << "ok!";
  }

  if (conf->saveStatusMap) {
    cout << endl << INFO_STR << "Saving statusmap and cavities...";
    dg->saveStatus(refName);
    cout << "ok!";
  }

  if (conf->saveCavities) {
    cout << endl << INFO_STR << "Saving cavities...";
    dg->saveCavities(false);
    cout << "ok!";
  }

  if (conf->saveIdebmap) {
    cout << endl << INFO_STR << "Saving idebmap...";
    dg->saveIdebMap(refName);
    cout << "ok!";
  }
}

void pocketMode(bool hasAtomInfo, ConfigFileOP cf, ConfigurationOP conf) {
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

  Surface *surf1, *surf2, *surf3;
  double areaSurf = 0, volSurf = 0;

  char mol[100] = "temp.txt", refName[100] = "mol";

  /////////////////////////////////////////////////////////////////////////////////////////
  // Set up Surface 1 (fat probe)
  cout << endl;
  cout << endl << INFO_STR << "Step 1 -> fat probe";
  DelPhiShared *dg1 =
      new DelPhiShared(conf->scale, conf->perfill, conf->molFile, localEpsMap,
                       localStatusMap, localMulti, hasAtomInfo);
  int natom = dg1->getNumAtoms();
  cf->remove(string("Surface"));
  cf->add<string>("Surface", "ses");
  surf1 = surfaceFactory().create(cf.get(), dg1);
  surf1->setInsideCode(5);
  surf1->setProbeRadius(conf->pocketRadiusBig);
  surf1->setProjBGP(false);
  surf1->setKeepWellShapedCavities(false);

  // Pre-process surface
  bool outsurf = surf1->build();

  if (!outsurf) {
    cout << endl << ERR << "Surface 1 construction failed!";
    exit(-1);
  }

  // fat connolly cancel each cavity (anyway they are smaller than they should
  // be at the end)
  surf1->getSurf(true, INFINITY);

  //////////////////////////////////////////////////////////////////////////////////////////
  // the subsequent surfaces don't need to read atom info

  // Set up Surface 2 (regular probe 1.4)
  // Set up DelPhi grid  2
  cout << endl;
  cout << endl << INFO_STR << "Step 2 -> small probe";

  DelPhiShared *dg2 =
      new DelPhiShared(conf->scale, conf->perfill, conf->molFile, localEpsMap,
                       localStatusMap, localMulti, false);
  surf2 = surfaceFactory().create(cf.get(), dg2);
  surf2->setInsideCode(10);
  surf2->setProjBGP(false);
  surf2->setProbeRadius(conf->pocketRadiusSmall);
  surf2->setKeepWellShapedCavities(false);

  /*
  //// to check percolation ///////////////////////
  DelPhiShared* dg2 = new
  DelPhiShared(conf.scale,conf.perfill,conf.mol,localStatusMap,localMulti,false);
  surf2 = surfaceFactory().create(cf,dg2);
  surf2->inside = 10;
  ////////////////////////////////////////////////
  */

  // Pre-process surface
  outsurf = surf2->build();

  if (!outsurf) {
    cout << endl << ERR << "Surface 2 construction failed!";
    exit(-1);
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
        exit(-1);
  */
  //////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////
  DelPhiShared *dg3 = NULL;
  if (conf->linkPockets) {
    // Set up Surface 3 (accessibility probe)
    dg3 = new DelPhiShared(conf->scale, conf->perfill, conf->molFile,
                           localEpsMap, localStatusMap, localMulti, false);
    surf3 = surfaceFactory().create(cf.get(), dg3);
    surf2->setProjBGP(false);
    surf2->setProbeRadius(conf->pocketRadiusSmall);
    surf2->setKeepWellShapedCavities(false);
    surf3->setInsideCode(15);

    // Pre-process surface
    outsurf = surf3->build();

    if (!outsurf) {
      cout << endl << ERR << "Surface 3 construction failed!";
      exit(-1);
    }
    // keep original surface
    surf3->getSurf(false);
  }

  cout << endl;
  cout << endl << INFO_STR << "Step 3 -> differential map";

  cout << endl << INFO_STR << "Building pockets by difference map...";
  (*surf1) -= (*surf2);

  ////////////////////// recover split cavities links ///////////////
  if (conf->linkPockets) {
    int nr = 1;
    double duration1 = 0;

    Timer chrono1;
    chrono1.start();

    cout << endl;
    cout << endl << INFO_STR << "Step 3.1 -> linking";

    cout << endl << INFO_STR << "Linking cavities/pockets...";
    cout.flush();
    while (nr != 0) {
      // check cavities links, use the link status map as reference map
      // to check accessibility
      nr = surf1->linkCavities(dg2->status, dg3->status);
      cout << endl << INFO_STR << "Merged " << nr << " cavities";
      cout.flush();
    }

    duration1 = chrono1.stop();
    cout << endl << INFO_STR << "Diff. Step 4 " << duration1 << " [s]";

    delete surf3;
    delete dg3;
  }
  ///////////////////////////////////////////////////////////////////

  cout << endl;
  cout << endl << INFO_STR << "Step 4 -> filtering, envelope building";

  ///////////////////// volume filter ///////////////////////////////
  surf1->fillCavities(11.4 * conf->numMol, false);
  surf1->getCavitiesAtoms();
  cout << endl << INFO_STR << "Saving cavities info..";

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
    cout << endl;
    // do cavity detection to recover which is cavity in the reference
    cout << endl << INFO_STR << "Recovering cavity/pocket distinction..";
    cout.flush();
    surf2->getCavities();
    // mark which are pockets and which are cavities
    dg1->markPockets(dg2->status, isPocket);
  }

  cout << endl << INFO_STR << "Build the envelope of each cavity/pocket";
  cout.flush();

  cf->remove("Surface");
  cf->add<string>("Surface", "skin");

  // TODO make it an option
  bool saveEntranceInfo = true;

  for (int i = 0; i < nc; i++) {
    char mol[100], tri_[100];
    snprintf(mol, sizeof(mol), "cav%d.txt", i);
    snprintf(tri_, sizeof(tri_), "cav_tri%d", i);
    DelPhiShared *dg_temp = new DelPhiShared(2.0, 50, mol, false, false, false);
    // reset seed
    srand(conf->currentSeed);
    Surface *cs_temp = surfaceFactory().create(cf.get(), dg_temp);
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
      cout << endl << INFO_STR << "Triangulating enveloped cavity/pockets " << i;
      cs_temp->triangulateSurface(0.0, tri_);
      areas.push_back(cs_temp->getArea());

      // NanoShaper identifies the set of grid points that are nearest
      // to the given triangulation points and save those points that are
      // completely solvent exposed (i.e. external) in surf2 (the slim surface).
      // Those points represent a point-wise representation of the entrance of
      // the pocket
      if (saveEntranceInfo) {
        if (!isPocket[i]) {
          // cout << endl << ERR << "Cannot find the entrance of a cavity; it
          // must be a pocket.";
          if (!conf->cavAndPockets) {
            cout << endl
                 << REMARK
                 << "You have to enable Cavities and pockets flag to do a "
                    "distinction between a pocket and a cavity";
            exit(-1);
          }
        } else {
          // is a pocket, we can identify the entrance
          // true means that it is part of the entrance
          vector<bool> results;
          cs_temp->triangulationPointsAreCompletelyOut(surf2, results);

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
    delete cs_temp;
    delete dg_temp;
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
    cout << endl << INFO_STR;
    cout << endl << INFO_STR << "------------------------------------";
    cout << endl << INFO_STR << "     Pocket Detection Summary       ";
    cout << endl << INFO_STR << "------------------------------------";
    cout << endl << INFO_STR;
    cout << endl
         << INFO_STR << "Detected a total of " << nc
         << " pockets/cavities having at least the volume of " << conf->numMol
         << " water molecules";

    for (int i = 0, ii = 0; i < nc; i++) {
      if (conf->cavAndPockets) {
        if (conf->tri) {
          if (isPocket[i]) {
            cout << endl
                 << INFO_STR << "Pocket " << i << " vol " << volumes[i] << " area "
                 << areas[i] << " body area " << areas2[ii];
            if (conf->debugStatus) {
              *(conf->internals) << endl << "pocket_vol " << volumes[i];
              *(conf->internals) << endl
                                 << "pocket_area " << areas[i]
                                 << " pocket_body_area " << areas2[i];
            }
            ii++;
          } else {
            cout << endl
                 << INFO_STR << "Cavity " << i << " vol " << volumes[i] << " area "
                 << areas[i];
            if (conf->debugStatus) {
              *(conf->internals) << endl << "cav_vol " << volumes[i];
              *(conf->internals) << endl << "cav_area " << areas[i];
            }
          }
        } else {
          if (isPocket[i]) {
            cout << endl << INFO_STR << "Pocket " << i << " vol " << volumes[i];
            if (conf->debugStatus) {
              *(conf->internals) << endl << "pocket_vol " << volumes[i];
            }
          } else {
            cout << endl << INFO_STR << "Cavity " << i << " vol " << volumes[i];
            if (conf->debugStatus) {
              *(conf->internals) << endl << "cav_vol " << volumes[i];
            }
          }
        }
      } else {
        if (conf->tri) {
          cout << endl
               << INFO_STR << "Pocket " << i << " vol " << volumes[i] << " area "
               << areas[i] << " body area " << areas2[ii];
          if (conf->debugStatus) {
            *(conf->internals) << endl << "pocket_vol " << volumes[i];
            *(conf->internals) << endl << "pocket_area " << areas[i];
          }
          ii++;
        } else {
          cout << endl << INFO_STR << "Pocket " << i << " vol " << volumes[i];
          if (conf->debugStatus) {
            *(conf->internals) << endl << "pocket_vol " << volumes[i];
          }
        }
      }
    }
  }

  cout << endl << INFO_STR << "Pocket detection time.. " << duration << " [s]";
  cout << endl << INFO_STR << "Cleaning memory...";
  cout.flush();

  delete surf1;
  delete dg1;
  delete surf2;
  delete dg2;
  cout << "ok!";
}
