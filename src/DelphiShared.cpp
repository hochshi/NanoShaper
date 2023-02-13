
#include <DelphiShared.h>
#include <logging.h>
#include <tools.h>
#include <stdexcept>

void DelPhiShared::init() {
  epsmap = NULL;
  buildEpsMap = false;
  buildStatus = false;
  x = NULL;
  y = NULL;
  z = NULL;
  idebmap = NULL;
  atoms = NULL;
  scspos = NULL;
  scsnor = NULL;
  scsarea = NULL;
  ibgp = NULL;
  status = NULL;
  cavitiesVec = NULL;
  delphiBinding = false;
  multi_diel = false;
  tempStatus = NULL;
  isAvailableAtomInfo = false;
  atsurf = NULL;
}

DelPhiShared::DelPhiShared() {
  init();
}

DelPhiShared::DelPhiShared(bool map, bool status, bool multi, bool atinfo)
    : buildEpsMap(map),
      buildStatus(status),
      multi_diel(multi),
      isAvailableAtomInfo(atinfo) {}

void DelPhiShared::init(double scale, double perfill, string fn, bool eps_flag,
                        bool stat_flag, bool multi, bool atinfo) {
  buildEpsMap = eps_flag;
  buildStatus = stat_flag;
  multi_diel = multi;
  isAvailableAtomInfo = atinfo;

  init(scale, perfill, fn);
}

void DelPhiShared::init(double scale, double perfill, string fn) {
  logging::log<logging::level::info>("Loading atoms....");

  bool flag = loadAtoms(fn);

  if (!flag) {
    logging::log<logging::level::err>(
        "Missing or corrupt atoms file. Initialization failed");
    throw std::invalid_argument(
        "Missing or corrupt atoms file. Initialization failed");
  }

  flag = buildGrid(scale, perfill);
  if (!flag) {
    logging::log<logging::level::err>("Initialization failed");
    throw std::runtime_error("Initialization failed");
  }
  logging::log<logging::level::info>("Initialization completed");
}

void DelPhiShared::init(double scale, double perfill, const InputData& in) {
  logging::log<logging::level::info>("Loading atoms....");

  int retNum = loadAtoms(in.na, in.x, in.y, in.z, in.r, in.q, in.d, in.ai);

  if (retNum != in.na) {
    logging::log<logging::level::err>(
        "Missing or corrupt atoms input data. Initialization failed");
    throw std::invalid_argument(
        "Missing or corrupt atoms input data. Initialization failed");
  }

  bool flag = buildGrid(scale, perfill);
  if (!flag) {
    logging::log<logging::level::err>("Initialization failed");
    throw std::runtime_error("Initialization failed");
  }
  logging::log<logging::level::info>("Initialization completed");
}

DelPhiShared::DelPhiShared(double scale, double perfill, string fn,
                           bool eps_flag, bool stat_flag, bool multi,
                           bool atinfo) {
  init();
  init(scale, perfill, fn, eps_flag, stat_flag, multi, atinfo);
}

void DelPhiShared::DelPhiBinding(double xmin, double ymin, double zmin,
                                 double xmax, double ymax, double zmax,
                                 double c1, double c2, double c3, double rmax,
                                 double perf, int* local_i_epsmap, int igrid,
                                 double scale, double* local_i_scspos,
                                 double* local_i_scsnor, bool* local_i_idebmap,
                                 int* local_i_ibgp, int maxbgp, bool bstatus,
                                 int* _atsurf) {
  logging::log<logging::level::info>("DelPhi Binding...");
  delphiBinding = true;
  buildStatus = bstatus;
  atsurf = _atsurf;

  if (x != NULL)
    deleteVector<double>(x);
  if (y != NULL)
    deleteVector<double>(y);
  if (z != NULL)
    deleteVector<double>(z);

  if (buildStatus)
    if (status != NULL)
      deleteVector<short>(status);

  this->xmin = xmin;
  this->ymin = ymin;
  this->zmin = zmin;

  this->rmaxdim = rmax;
  this->perfill = perf;
  this->xmax = xmax;
  this->ymax = ymax;
  this->zmax = zmax;

  this->baricenter[0] = c1;
  this->baricenter[1] = c2;
  this->baricenter[2] = c3;

  x = allocateVector<double>(igrid);
  y = allocateVector<double>(igrid);
  z = allocateVector<double>(igrid);

  this->scale = scale;
  side = 1 / scale;
  hside = side / 2;
  A = side * side;
  nx = igrid;
  ny = igrid;
  nz = igrid;

  if (atoms != NULL) {
    for (int i = 0; i < numAtoms; i++)
      delete atoms[i];
    delete[] atoms;
  }

  atoms = NULL;

  // assuming all maps are already allocated before binding
  this->epsmap = local_i_epsmap;
  this->idebmap = local_i_idebmap;
  this->scsnor = local_i_scsnor;
  this->scsarea = NULL;
  this->scspos = local_i_scspos;
  this->ibgp = local_i_ibgp;

  for (int i = 0; i < nx; i++)
    x[i] = xmin + i * side;

  for (int i = 0; i < ny; i++)
    y[i] = ymin + i * side;

  for (int i = 0; i < nz; i++)
    z[i] = zmin + i * side;

  // allocate grid point status memory
  // and initialize all to temporary outside status
  if (buildStatus) {
    /*
    status = allocateMatrix3D<short>(nx,ny,nz);
    for (int i=0;i<nx;i++)
            for (int j=0;j<ny;j++)
                    for (int k=0;k<nz;k++)
                            status[i][j][k]=2;
    */
    int tot = nx * ny * nz;
    status = allocateVector<short>(tot);
    for (int i = 0; i < tot; i++)
      status[i] = STATUS_POINT_TEMPORARY_OUT;
  }
  // set max bgp
  this->maxbgp = maxbgp;

  logging::log<logging::level::info>("ok!");
  logging::log<logging::level::info>(
      "Max number of bgps was set from DelPhi to {}", maxbgp);
}

void DelPhiShared::finalizeBinding(int* ibnum) {
  // update indexes to be conformant to Fortran
  for (int i = 0; i < nbgp; i++) {
    ibgp[3 * i]++;
    ibgp[3 * i + 1]++;
    ibgp[3 * i + 2]++;
  }

  // update indexes to be conformant to Fortran
  if (atsurf != NULL)
    for (int i = 0; i < nbgp; i++)
      atsurf[i]++;

  epsmap = NULL;   // avoid deletion
  scspos = NULL;   // avoid deletion
  scsnor = NULL;   // avoid deletion
  idebmap = NULL;  // avoid deletion
  // scsarea = NULL; //for now not used, so it must be distroied
  ibgp = NULL;  // avoid deletion
  // save ibnum
  (*ibnum) = nbgp;
}

int DelPhiShared::loadAtoms(string fn) {
  // load atoms file
  ifstream fin;
  fin.open(fn.c_str(), ios::in);
  char buffer[BUFLEN];
  double max_rad = 0;

  if (multi_diel)
    logging::log<logging::level::info>("Atom dielectric info is available..");

  if (isAvailableAtomInfo)
    logging::log<logging::level::info>("Atom Info is available..");

  if (fin.fail()) {
    logging::log<logging::level::warn>("Cannot read file {}", fn);
    return 0;
  }

  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<double> r;
  vector<double> q;
  vector<int> d;
  vector<AtomInfo> ai;

  double x_;
  double y_;
  double z_;
  double r_;
  double q_;
  int d_;
  char name[10], resName[10], chain[2];
  int resid;

  while (!fin.eof()) {
    fin.getline(buffer, BUFLEN);
    if (strlen(buffer) > 1) {
      d_ = -1;
      if (!multi_diel && !isAvailableAtomInfo) {
        sscanf(buffer, "%lf %lf %lf %lf", &x_, &y_, &z_, &r_);
      } else if (multi_diel) {
        sscanf(buffer, "%lf %lf %lf %lf %d", &x_, &y_, &z_, &r_, &d_);
      } else if (isAvailableAtomInfo) {
        sscanf(buffer, "%lf %lf %lf %lf %d %s %s %d %s", &x_, &y_, &z_, &r_,
               &d_, name, resName, &resid, chain);
      }

      if (isAvailableAtomInfo || multi_diel) {
        if (-1 == d_) {
          logging::log<logging::level::err>(
              "Cannot get the dielectric value of the atom number {}; please "
              "add it after radius entry to the xyzr input file",
              x.size());
          throw std::runtime_error(
              "Cannot get the dielectric value of an atom number; please "
              "add it after radius entry to the xyzr input file");
        }
        d.push_back(d_);
        ai.emplace_back(string(name), resid, string(resName), string(chain));
      }
      max_rad = MAX(r_, max_rad);
      x.push_back(x_);
      y.push_back(y_);
      z.push_back(z_);
      r.push_back(r_);
    }
  }
  fin.close();

  numAtoms = (int)x.size();
  logging::log<logging::level::info>("Read {} atoms", numAtoms);

  if (numAtoms < 4) {
    logging::log<logging::level::err>(
        "NanoShaper needs at least 4 atoms to work.");
    logging::log<logging::level::info>(
        "{} To emulate 4 atoms you can place dummy atoms with null radius",
        REMARK);
    logging::log<logging::level::info>(
        "{} at the same centers of the real atoms", REMARK);
    throw std::invalid_argument("NanoShaper needs at least 4 atoms to work.");
  }

  if (max_rad == 0) {
    logging::log<logging::level::err>(
        "All null radii? If you are using place holders atoms please "
        "set at least one radius > 0");
    throw std::invalid_argument(
        "All null radii? If you are using place holders atoms please "
        "set at least one radius > 0");
  }

  return loadAtoms(numAtoms, x.data(), y.data(), z.data(), r.data(), q.data(),
                   d.data(), ai.data());
}

int DelPhiShared::loadAtoms(int na, double* pos, double* rad, double* charge,
                            int* dielec, char* atinf) {
  logging::log<logging::level::info>("Read {} atoms", na);

  if (rad == NULL || pos == NULL) {
    logging::log<logging::level::warn>("Atoms info not available!");
    return 0;
  }

  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<double> r;
  vector<double> q;
  vector<int> d;
  vector<AtomInfo> ai;

  double x_;
  double y_;
  double z_;
  double r_;
  double q_;
  int d_;
  string name, resName, chain;
  int resid;
  int lastChar = 0;

  for (int i = 0; i < na; i++) {
    // double qq;
    // int dd;
    if (charge == NULL)
      q_ = 0;
    else
      q_ = q[i];

    if (dielec == NULL)
      d_ = 0;
    else
      d_ = d[i];

    // parse atinf variable and store if available
    if (atinf != NULL) {
      isAvailableAtomInfo = true;
      char buff[100];
      for (int kk = 0; kk < 15; kk++)
        buff[kk] = atinf[lastChar++];
      buff[15] = '\0';
      string str(buff);
      name = str.substr(0, 5);
      // remove white spaces
      name.erase(remove_if(name.begin(), name.end(),
                           static_cast<int (*)(int)>(isspace)),
                 name.end());
      resName = str.substr(6, 3);
      // remove white spaces
      resName.erase(remove_if(resName.begin(), resName.end(),
                              static_cast<int (*)(int)>(isspace)),
                    resName.end());
      resid = atoi(str.substr(11, 4).c_str());
      chain = str.substr(10, 1);
      ai.emplace_back(name, resid, resName, chain);
    }
    x.push_back(pos[i * 3]);
    y.push_back(pos[i * 3 + 1]);
    z.push_back(pos[i * 3 + 2]);
    r.push_back(r[i]);
    q.push_back(q_);
    d.push_back(d_);
  }

  return loadAtoms(numAtoms, x.data(), y.data(), z.data(), r.data(), q.data(),
                   d.data(), ai.data());
}

int DelPhiShared::loadAtoms(int na, const std::vector<double>& x,
                            const std::vector<double>& y,
                            const std::vector<double>& z,
                            const std::vector<double>& r,
                            const std::vector<double>& q,
                            const std::vector<int>& d,
                            const std::vector<AtomInfo>& ai) {
  return loadAtoms(na, x.data(), y.data(), z.data(), r.data(), q.data(),
                   d.data(), ai.data());
}

int DelPhiShared::loadAtoms(int na, const double* x, const double* y,
                            const double* z, const double* r, const double* q,
                            const int* d, const AtomInfo* ai) {

  if (atoms != NULL) {
    for (int i = 0; i < numAtoms; i++)
      delete atoms[i];
    delete[] atoms;
  }

  atoms = new Atom*[na];
  for (int i = 0; i < na; i++) {
    atoms[i] = new Atom(x[i], y[i], z[i], r[i]);
    if (NULL != q) {
      atoms[i]->charge = q[i];
    }
    if (NULL != d) {
      atoms[i]->dielectric = d[i];
    }
    if (NULL != ai) {
      atoms[i]->ai = ai[i];
    }
  }
  numAtoms = na;

  return numAtoms;
}

bool DelPhiShared::buildGrid(double scale, double perfill) {
  if (atoms == NULL) {
    logging::log<logging::level::err>("Cannot build grid with no atoms!");
    return false;
  }

  if (scale < 0) {
    logging::log<logging::level::warn>("Cannot use a <0 scale: setting {}",
                                       DEFAULT_SCALE);
    scale = DEFAULT_SCALE;
  }

  if (perfill < 0 || perfill > 100) {
    logging::log<logging::level::warn>("Perfil is in (0,100). Setting {}",
                                       DEFAULT_PERFIL);
    scale = DEFAULT_PERFIL;
  }

  this->scale = scale;
  this->perfill = perfill;

  double cmin[3] = {INFINITY, INFINITY, INFINITY};
  double cmax[3] = {-INFINITY, -INFINITY, -INFINITY};
  double oldmid[3];

  for (int ia = 0; ia < numAtoms; ia++) {
    cmin[0] = min(cmin[0], atoms[ia]->pos[0] - atoms[ia]->radius);
    cmin[1] = min(cmin[1], atoms[ia]->pos[1] - atoms[ia]->radius);
    cmin[2] = min(cmin[2], atoms[ia]->pos[2] - atoms[ia]->radius);

    cmax[0] = max(cmax[0], atoms[ia]->pos[0] + atoms[ia]->radius);
    cmax[1] = max(cmax[1], atoms[ia]->pos[1] + atoms[ia]->radius);
    cmax[2] = max(cmax[2], atoms[ia]->pos[2] + atoms[ia]->radius);
  }

  oldmid[0] = (cmax[0] + cmin[0]) / 2.;
  oldmid[1] = (cmax[1] + cmin[1]) / 2.;
  oldmid[2] = (cmax[2] + cmin[2]) / 2.;

  logging::log<logging::level::info>("Geometric baricenter ->  {} {} {}",
                                     oldmid[0], oldmid[1], oldmid[2]);

  double v[6];
  v[0] = fabs(cmax[0] - oldmid[0]);
  v[1] = fabs(cmin[0] - oldmid[0]);
  v[2] = fabs(cmax[1] - oldmid[1]);
  v[3] = fabs(cmin[1] - oldmid[1]);
  v[4] = fabs(cmax[2] - oldmid[2]);
  v[5] = fabs(cmin[2] - oldmid[2]);

  rmaxdim = -INFINITY;

  for (int i = 0; i < 6; i++)
    if (rmaxdim < v[i])
      rmaxdim = v[i];

  rmaxdim = 2 * rmaxdim;

  unsigned int igrid = (unsigned int)floor(scale * 100. / perfill * rmaxdim);

  if ((igrid % 2) == 0)
    igrid++;

  logging::log<logging::level::info>("Grid is {}", igrid);

  baricenter[0] = oldmid[0];
  baricenter[1] = oldmid[1];
  baricenter[2] = oldmid[2];

  xmin = oldmid[0] - (igrid - 1) / (2 * scale);
  ymin = oldmid[1] - (igrid - 1) / (2 * scale);
  zmin = oldmid[2] - (igrid - 1) / (2 * scale);

  xmax = oldmid[0] + (igrid - 1) / (2 * scale);
  ymax = oldmid[1] + (igrid - 1) / (2 * scale);
  zmax = oldmid[2] + (igrid - 1) / (2 * scale);

  logging::log<logging::level::info>("MAX {} {} {}", xmax, ymax, zmax);
  logging::log<logging::level::info>("MIN {} {} {}", xmin, ymin, zmin);
  logging::log<logging::level::info>("Perfil {}%", perfill);
  logging::log<logging::level::info>("Rmaxdim {}", rmaxdim);

  logging::log<logging::level::info>("Allocating memory...");

  if (x != NULL)
    deleteVector<double>(x);
  if (y != NULL)
    deleteVector<double>(y);
  if (z != NULL)
    deleteVector<double>(z);

  x = allocateVector<double>(igrid);
  y = allocateVector<double>(igrid);
  z = allocateVector<double>(igrid);

  side = 1 / scale;
  hside = side / 2;
  A = side * side;

  nx = igrid;
  ny = igrid;
  nz = igrid;

  for (int i = 0; i < nx; i++)
    x[i] = xmin + i * side;

  for (int i = 0; i < ny; i++)
    y[i] = ymin + i * side;

  for (int i = 0; i < nz; i++)
    z[i] = zmin + i * side;

  // allocate epsmap memory
  bool flag = true;

  if (buildEpsMap) {
    flag = clearAndAllocEpsMaps();

    if (!flag)
      return false;

    if (idebmap != NULL)
      deleteVector<bool>(idebmap);

    int tot = nx * ny * nz;
    idebmap = allocateVector<bool>(tot);
    for (int i = 0; i < tot; i++)
      idebmap[i] = true;
  }

  // allocate grid point status memory
  // and initialize all to temporary outside status
  if (buildStatus) {
    if (status != NULL)
      deleteVector<short>(status);

    int tot = nx * ny * nz;
    status = allocateVector<short>(tot);
    // int alignementBits = 64;
    // status = allocateAlignedVector<short>(tot,alignementBits);

    if (status == NULL) {
      logging::log<logging::level::err>(
          "Not enough memory to allocate status map");
      throw std::runtime_error("Not enough memory to allocate status map");
    }

    for (int i = 0; i < tot; i++)
      status[i] = STATUS_POINT_TEMPORARY_OUT;
  }

  logging::log<logging::level::info>("ok!");
  return true;
}

bool DelPhiShared::clearAndAllocEpsMaps() {
  if (epsmap != NULL)
    deleteVector<int>(epsmap);

  size_t tot = size_t(nx) * size_t(ny) * size_t(nz) * 3;
  epsmap = allocateVector<int>(tot);

  if (epsmap == NULL) {
    logging::log<logging::level::err>("Error allocating epsmap!");
    return false;
  }

  // init epsmap
  for (size_t i = 0; i < tot; i++)
    epsmap[i] = 0;

  return true;
}

void DelPhiShared::clearEpsMaps() {

  if (epsmap == NULL) {
    logging::log<logging::level::err>("Cannot clear null map!");
    return;
  }

  size_t n = nx * ny * nz * 3;
  for (size_t i = 0; i < n; i++)
    epsmap[i] = 0;
  /*
  for (int i=0;i<nx;i++)
                  for (int j=0;j<ny;j++)
                          for (int k=0;k<nz;k++)
                          {
                                  EPSMAP(i,j,k,0,nx,ny,nz)=0;
                                  EPSMAP(i,j,k,1,nx,ny,nz)=0;
                                  EPSMAP(i,j,k,2,nx,ny,nz)=0;
                          }
  */
  return;
}

void DelPhiShared::saveIdebMap(char* fname) {
  if (idebmap == NULL) {
    logging::log<logging::level::warn>("Cannot save null idebmap!");
    return;
  }

  char f1[BUFLEN];
  // if (debug_stern == -1)
  snprintf(f1, sizeof(f1), "%s.idebmap.txt", fname);
  // else
  //	snprintf(f1, sizeof(f1), "%s.idebmap.stern.txt",fname);
  FILE* fp = fopen(f1, "w");
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++)
        // fprintf(fp,"%d ",(int)IDEBMAP(i,j,k,nx,ny));
        fprintf(fp, "%d ", read3DVector<bool>(idebmap, i, j, k, nx, ny, nz));

      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void DelPhiShared::saveEpsMaps(char* fname) {
  if (epsmap == NULL) {
    logging::log<logging::level::warn>("Cannot save null epsmap!");
    return;
  }

  char f1[BUFLEN], f2[BUFLEN], f3[BUFLEN];
  snprintf(f1, sizeof(f1), "%s.epsmapx.txt", fname);
  snprintf(f2, sizeof(f2), "%s.epsmapy.txt", fname);
  snprintf(f3, sizeof(f3), "%s.epsmapz.txt", fname);

  FILE* fp = fopen(f1, "w");
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++)
        // fprintf(fp,"%d ",EPSMAP(i,j,k,0,nx,ny,nz));
        fprintf(fp, "%d ",
                read4DVector<int>(epsmap, i, j, k, 0, nx, ny, nz, 3));

      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  fp = fopen(f2, "w");
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++)
        // fprintf(fp,"%d ",EPSMAP(i,j,k,1,nx,ny,nz));
        fprintf(fp, "%d ",
                read4DVector<int>(epsmap, i, j, k, 1, nx, ny, nz, 3));

      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  fp = fopen(f3, "w");
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++)
        // fprintf(fp,"%d ",EPSMAP(i,j,k,2,nx,ny,nz));
        fprintf(fp, "%d ",
                read4DVector<int>(epsmap, i, j, k, 2, nx, ny, nz, 3));

      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

DelPhiShared::~DelPhiShared() {
  clear();
}

void DelPhiShared::clear() {
  if (epsmap != NULL)
    deleteVector<int>(epsmap);
  if (ibgp != NULL)
    deleteVector<int>(ibgp);
  if (scspos != NULL)
    deleteVector<double>(scspos);
  if (scsnor != NULL)
    deleteVector<double>(scsnor);
  if (scsarea != NULL)
    deleteVector<double>(scsarea);
  if (x != NULL)
    deleteVector<double>(x);
  if (y != NULL)
    deleteVector<double>(y);
  if (z != NULL)
    deleteVector<double>(z);
  if (status != NULL)
    // deleteMatrix3D<short>(nx,ny,status);
    deleteVector<short>(status);
  // deleteAlignedVector<short>(status);
  if (idebmap != NULL)
    deleteVector<bool>(idebmap);

  if (atoms != NULL) {
    for (int i = 0; i < numAtoms; i++)
      delete atoms[i];
    delete[] atoms;
  }

  // free memory
  if (cavitiesVec != NULL) {
    vector<vector<int*>*>::iterator it;
    for (it = cavitiesVec->begin(); it != cavitiesVec->end(); it++) {
      vector<int*>* inner = (*it);
      vector<int*>::iterator it2;
      for (it2 = inner->begin(); it2 != inner->end(); it2++)
        free((*it2));
      delete inner;
    }
  }
  delete cavitiesVec;
}

void DelPhiShared::saveBGP(char* fname) {
  if (scspos == NULL || scsnor == NULL) {
    logging::log<logging::level::warn>("Cannot save null scspos or scsnor!");
    return;
  }

  char ff[BUFLEN];
  snprintf(ff, sizeof(ff), "%s.projections.txt", fname);
  FILE* fp = fopen(ff, "w");
  FILE* fp2 = fopen("surf.txt", "w");
  fprintf(fp, "%d\n", nbgp);
  for (int iibgp = 0; iibgp < nbgp; iibgp++) {
    fprintf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf", ibgp[iibgp * 3] + 1,
            ibgp[iibgp * 3 + 1] + 1, ibgp[iibgp * 3 + 2] + 1, scspos[iibgp * 3],
            scspos[iibgp * 3 + 1], scspos[iibgp * 3 + 2], scsnor[iibgp * 3],
            scsnor[iibgp * 3 + 1], scsnor[iibgp * 3 + 2]);

    fprintf(fp2, "%lf %lf %lf %lf %lf %lf", scspos[iibgp * 3],
            scspos[iibgp * 3 + 1], scspos[iibgp * 3 + 2], scsnor[iibgp * 3],
            scsnor[iibgp * 3 + 1], scsnor[iibgp * 3 + 2]);

    if (iibgp < nbgp - 1) {
      fprintf(fp, "\n");
      fprintf(fp2, "\n");
    }
  }
  fclose(fp);
  fclose(fp2);
}

void DelPhiShared::saveStatus(char* fname) {
  if (status == NULL) {
    logging::log<logging::level::warn>("Cannot save null status!");
    return;
  }
  char ff[BUFLEN];
  snprintf(ff, sizeof(ff), "%s.status.txt", fname);
  FILE* fp = fopen(ff, "w");

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        // fprintf(fp,"%d ",status[i][j][k]);
        // fprintf(fp,"%d ",STATUSMAP(i,j,k,nx,ny));
        fprintf(fp, "%d ", read3DVector<short>(status, i, j, k, nx, ny, nz));
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  saveCavities(false);
}

void DelPhiShared::initCav2Atoms() {
  int nc = (int)cavitiesVec->size();

  if (cav2atoms.size() != 0)
    clearCav2Atoms();

  cav2atoms.reserve(nc);
  for (int ic = 0; ic < nc; ic++)
    cav2atoms.push_back(new set<int>());
}

void DelPhiShared::clearCav2Atoms() {
  if (cav2atoms.size() != 0) {
    for (unsigned int ic = 0; ic < cav2atoms.size(); ic++)
      delete cav2atoms[ic];
  }
  cav2atoms.clear();
}

void DelPhiShared::markPockets(short* stat, vector<bool>& isPocket) {
  vector<vector<int*>*>::iterator it;
  int num = (int)cavitiesVec->size();

  if (stat == NULL) {
    logging::log<logging::level::warn>(
        "Cannot mark pockets with a reference status map!");
    return;
  }

  if (stat == status) {
    logging::log<logging::level::warn>(
        "In marking pockets the status map passed must be from "
        "another reference object; here the same was passed");
    return;
  }

  short* status = stat;

  if (num == 0)
    return;

  int i = 0, sync_i = 0;

  for (i = 0, sync_i = 0, it = cavitiesVec->begin(); it != cavitiesVec->end();
       it++, i++, sync_i++) {
    if (cavitiesFlag[sync_i] == true) {
      i--;
      continue;
    }

    vector<int*>::iterator it2;
    vector<int*>* vec = (*it);

    bool isPocketFlag = false;

    for (it2 = vec->begin(); it2 != vec->end(); it2++) {
      int* vec = *it2;
      short value =
          read3DVector<short>(status, vec[0], vec[1], vec[2], nx, ny, nz);
      // if at least one was out, then that's a pocket
      // if
      // (STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny)==STATUS_POINT_TEMPORARY_OUT
      // || STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny)==STATUS_POINT_OUT)
      if (value == STATUS_POINT_TEMPORARY_OUT || value == STATUS_POINT_OUT) {
        isPocketFlag = true;
        break;
      }
    }
    isPocket.push_back(isPocketFlag);
  }
  return;
}

void DelPhiShared::parseAtomInfo(char* atinf, string mol, int natom,
                                 double* xn1, double* rad) {
  // get atoms and info from delphi
  FILE* ff = fopen(mol.c_str(), "w");
  int lastChar = 0;

  // parse atinfo into a file
  for (int i = 0; i < natom; i++) {
    char buff[100];
    for (int kk = 0; kk < 15; kk++)
      buff[kk] = atinf[lastChar++];

    buff[15] = '\0';
    string str(buff);
    string name = str.substr(0, 5);
    // avoid disambiguation
    name.erase(
        remove_if(name.begin(), name.end(), static_cast<int (*)(int)>(isspace)),
        name.end());
    string resName = str.substr(6, 3);
    resName.erase(remove_if(resName.begin(), resName.end(),
                            static_cast<int (*)(int)>(isspace)),
                  resName.end());
    fprintf(ff, "%f %f %f %f 0 %s %s %d %s\n", xn1[3 * i], xn1[3 * i + 1],
            xn1[3 * i + 2], rad[i], name.c_str(), resName.c_str(),
            atoi(str.substr(11, 4).c_str()), str.substr(10, 1).c_str());
  }
  fclose(ff);
}

int DelPhiShared::cavitiesToAtoms(double rad) {
  vector<vector<int*>*>::iterator it;
  int num = (int)cavitiesVec->size();

  if (num == 0)
    return 0;

  char buff[100];
  int i = 0, sync_i = 0;
  int nonActive = 0;
  FILE* fp;

  for (i = 0, sync_i = 0, it = cavitiesVec->begin(); it != cavitiesVec->end();
       it++, i++, sync_i++) {
    if (cavitiesFlag[sync_i] == true) {
      nonActive++;
      i--;
      continue;
    }
    snprintf(buff, sizeof(buff), "cav%d.txt", i);
    fp = fopen(buff, "w");
    snprintf(buff, sizeof(buff), "all_cav%d.txt", i);
    FILE* fp2 = fopen(buff, "w");

    vector<int*>::iterator it2;
    vector<int*>* vec = (*it);
    int cavityId = sync_i + STATUS_FIRST_CAV;

    int count = 0;
    double lastx, lasty, lastz;
    for (it2 = vec->begin(); it2 != vec->end(); it2++) {
      int* vec = *it2;
      fprintf(fp2, "%f %f %f %f\n", x[vec[0]], y[vec[1]], z[vec[2]], rad);

      // save only support cavity points
      // if (STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny)==-cavityId)
      if (read3DVector<short>(status, vec[0], vec[1], vec[2], nx, ny, nz) ==
          -cavityId) {
        // printf("\n %d",STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny));
        fprintf(fp, "%f %f %f %f\n", x[vec[0]], y[vec[1]], z[vec[2]], rad);
        lastx = x[vec[0]];
        lasty = y[vec[1]];
        lastz = z[vec[2]];
        count++;
      }
    }

    if (count == 0) {
      logging::log<logging::level::err>("Zero support atoms to save");
      throw std::runtime_error("Zero support atoms to save");
    }

    // clone atoms to let NanoShaper work on this set of points
    if (count < 4) {
      count = 4 - count;
      for (int i = 0; i < count; i++)
        fprintf(fp, "%f %f %f %f\n", lastx, lasty, lastz, rad);
    }
    fclose(fp);
    fclose(fp2);
  }
  fp = fopen("numcav.txt", "w");
  fprintf(fp, "%d", num - nonActive);
  fclose(fp);

  return num - nonActive;
}

void DelPhiShared::saveCavities2(bool onlySaveNonFilled, string sysName) {
  if (isAvailableAtomInfo) {
    FILE* fp;
    char buff[100];
    snprintf(buff, sizeof(buff), "%s.pocket", sysName.c_str());
    fp = fopen(buff, "w");
    int savedIndex = 1;
    for (unsigned int i = 0; i < cav2atoms.size(); i++) {
      set<int>::iterator setIt;
      set<int>* setPt = cav2atoms[i];

      // save only non filled if required, or save all of them
      if ((onlySaveNonFilled && !cavitiesFlag[i]) || (!onlySaveNonFilled)) {
        for (setIt = setPt->begin(); setIt != setPt->end(); setIt++) {
          // printf("\nIndex %d",i);
          fflush(stdout);
          int ii = (*setIt);
          Atom* at = atoms[ii];
          AtomInfo& ai = at->ai;
          if ((ai.getName()).size() == 4)
            fprintf(fp,
                    "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f               "
                    "  %6d POC\n",
                    ii + 1, ai.getName().c_str(), ai.getResName().c_str(),
                    ai.getChain().c_str(), ai.getResNum(), at->pos[0],
                    at->pos[1], at->pos[2], savedIndex);
          else
            fprintf(fp,
                    "ATOM  %5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f              "
                    "   %6d POC\n",
                    ii + 1, ai.getName().c_str(), ai.getResName().c_str(),
                    ai.getChain().c_str(), ai.getResNum(), at->pos[0],
                    at->pos[1], at->pos[2], savedIndex);
        }
        savedIndex++;
      }
    }
    fclose(fp);
  }
}
void DelPhiShared::saveCavities(bool onlySaveNonFilled) {
  if (status == NULL) {
    logging::log<logging::level::warn>(
        "Cannot save cavities with null status!");
    return;
  }

  FILE* fp2 = fopen("cavities.txt", "w");
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        // if (status[i][j][k]>3)
        // if ((STATUSMAP(i,j,k,nx,ny))>STATUS_POINT_OUT)
        if (read3DVector<short>(status, i, j, k, nx, ny, nz) > STATUS_POINT_OUT)
          fprintf(fp2, "\n %f %f %f", x[i], y[j], z[k]);
      }
    }
  }
  fclose(fp2);

  fp2 = fopen("cavitiesSize.txt", "w");
  vector<double>::iterator itCav;
  vector<int*>::iterator itCavVec;
  int i = 0;

  fprintf(fp2, "%f\n", ((float)cavitiesSize.size()));

  // for each cavity save the baricenter and the estimated volume
  for (itCav = cavitiesSize.begin(), i = 0; itCav != cavitiesSize.end();
       itCav++, i++) {
    double barx, bary, barz;
    barx = 0;
    bary = 0;
    barz = 0;
    for (itCavVec = cavitiesVec->at(i)->begin();
         itCavVec != cavitiesVec->at(i)->end(); itCavVec++) {
      barx += x[(*itCavVec)[0]];
      bary += y[(*itCavVec)[1]];
      barz += z[(*itCavVec)[2]];
    }
    double volSize = cavitiesSize.at(i);
    int size = (int)cavitiesVec->at(i)->size();
    barx /= size;
    bary /= size;
    barz /= size;
    int fillFlag;

    if (cavitiesFlag[i])
      fillFlag = 1;
    else
      fillFlag = 0;
    // save baricenter,volume, and filling flag
    fprintf(fp2, "%lf %lf %lf %lf %d\n", barx, bary, barz, volSize, fillFlag);
  }
  fclose(fp2);

  fp2 = fopen("cavAtomsSerials.txt", "w");
  for (unsigned int i = 0; i < cav2atoms.size(); i++) {
    set<int>::iterator setIt;
    set<int>* setPt = cav2atoms[i];

    // save only non filled if required, or save all of them
    if ((onlySaveNonFilled && !cavitiesFlag[i]) || (!onlySaveNonFilled)) {
      for (setIt = setPt->begin(); setIt != setPt->end(); setIt++)
        fprintf(fp2, "%d ", (*setIt) + 1);
      fprintf(fp2, "\n");
    }
  }
  fclose(fp2);

  if (isAvailableAtomInfo) {
    FILE* fp3;
    // logging::log<logging::level::info>("Atoms info is available I am saving
    // also residues");
    fp2 = fopen("residues.txt", "w");
    fp3 = fopen("residues_vmd.txt", "w");

    for (unsigned int i = 0; i < cav2atoms.size(); i++) {
      set<int>::iterator setIt;
      set<int>* setPt = cav2atoms[i];

      // only save the non filled cavities
      if (onlySaveNonFilled && !cavitiesFlag[i]) {
        // detected the involved residues in the cavity
        set<int> involvedResidues;
        for (setIt = setPt->begin(); setIt != setPt->end(); setIt++) {
          int atomIndex = (*setIt);
          involvedResidues.insert(atoms[atomIndex]->ai.getResNum());
        }

        int ll = 0, tot = (int)involvedResidues.size();

        for (setIt = involvedResidues.begin(); setIt != involvedResidues.end();
             setIt++, ll++) {
          int resIndex = (*setIt);
          fprintf(fp2, "%d ", resIndex);
          fprintf(fp3, "(resid %d and not water) ", resIndex);

          if (ll != (tot - 1))
            fprintf(fp3, " or ");
        }
        fprintf(fp2, "\n");
        fprintf(fp3, "\n");
      }
    }
    fclose(fp2);
    fclose(fp3);
  }
}
