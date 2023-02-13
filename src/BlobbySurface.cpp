#include <BlobbySurface.h>
#include <globals.h>
#include <logging.h>

void BlobbySurface::init() {
  B = DEFAULT_BLOBBYNESS;
  cutoff = DEFAULT_CUTOFF;
  surfType = MOLECULAR_SURFACE;
  providesAnalyticalNormals = false;
}

void BlobbySurface::init(ConfigurationOP cf) {
  // double blobby_B = cf->read<double>( "Blobbyness", -2.5 );
  // Set up inside value
  inside = 5;
  setBlobbyness(cf->blobby_B);
  // setBlobbyness(blobby_B);
}
BlobbySurface::BlobbySurface() { init(); }

BlobbySurface::BlobbySurface(DelPhiShared *ds) : MeshSurface(ds) { init(); }

BlobbySurface::BlobbySurface(ConfigurationOP cf, DelPhiShared *ds)
    : MeshSurface(ds) {
  init();
  init(cf);
}

BlobbySurface::~BlobbySurface() { clear(); }

void BlobbySurface::clear() {}

void BlobbySurface::setBlobbyness(double b) {
  if (b >= 0) {
    logging::log<logging::level::warn>(
        "Blobbyness is always a stricly negative real: setting {}",
        DEFAULT_BLOBBYNESS);
    b = DEFAULT_BLOBBYNESS;
  }
  B = b;
}

double BlobbySurface::getBlobbyness() { return B; }

bool BlobbySurface::build() {
  // build blobby in an atom-wise way in order to use a cutoff distance

  // numgrid points that are neighbour of an atom
  int numgrid = (int)(cutoff / delphi->side + 0.5);

  // disable cutoff
  // numgrid = delphi->nx;

  logging::log<logging::level::info>(
      "Using cut-off {} num neighbour grid points {}", cutoff, numgrid);

  if (scalarField != NULL)
    deleteMatrix3D<double>(last_nx, last_ny, scalarField);

  scalarField = allocateMatrix3D<double>(delphi->nx, delphi->ny, delphi->nz);

  if (scalarField == NULL) {
    logging::log<logging::level::err>("Cannot allocate scalar field!");
#ifdef PYTHON
    throw std::exception();
#else
    exit(-1);
#endif
  }

  for (int i = 0; i < delphi->nx; i++)
    for (int j = 0; j < delphi->nx; j++)
      for (int k = 0; k < delphi->nx; k++)
        scalarField[i][j][k] = 0.0;

  Atom **atoms = delphi->atoms;
  printf("\n");
  int na = delphi->numAtoms;
  // fill surface class structures
  for (int i = 0; i < na; i++) {
    printf("\r%sBlobby %.2f%%        ", INFO_STR, ((float)i + 1) / na * 100.0);
    double *pos = atoms[i]->pos;
    double r = atoms[i]->radius;
    double r2 = r * r;
    // get ref grid point
    int ix = (int)rintp((pos[0] - delphi->xmin) / delphi->side);
    int iy = (int)rintp((pos[1] - delphi->ymin) / delphi->side);
    int iz = (int)rintp((pos[2] - delphi->zmin) / delphi->side);

    int start_x = MAX(0, (ix - numgrid));
    int start_y = MAX(0, (iy - numgrid));
    int start_z = MAX(0, (iz - numgrid));

    int end_x = MIN((delphi->nx), (ix + numgrid));
    int end_y = MIN((delphi->ny), (iy + numgrid));
    int end_z = MIN((delphi->nz), (iz + numgrid));

    for (int ii = start_x; ii < end_x; ii++)
      for (int jj = start_y; jj < end_y; jj++)
        for (int kk = start_z; kk < end_z; kk++) {
          double dist2 = 0;
          double p[3];
          p[0] = delphi->x[ii] - delphi->hside;
          p[1] = delphi->y[jj] - delphi->hside;
          p[2] = delphi->z[kk] - delphi->hside;
          DIST2(dist2, p, pos);
          scalarField[ii][jj][kk] += exp(B * (dist2 / r2 - 1));
        }
  }

  bool old = accurateTriangulation;
  bool old_saveMSMS = saveMSMS;
  bool old_vertexAtomsMapFlag = vertexAtomsMapFlag;

  bool old_computeNormals = computeNormals;

  accurateTriangulation = true;
  isAvailableScalarField = true;

  setSaveMSMS(false);
  setVertexAtomsMap(false);
  setComputeNormals(false);

  // allocate empty
  allocIntersectionsMatrices();

  if (computeNormals)
    allocNormalsMatrices();

  triangulateSurface(1.0, "blobby");

  logging::log<logging::level::info>("Smoothing blobby surface...");
  smoothSurface("blobby");

  setSaveMSMS(old_saveMSMS);
  setVertexAtomsMap(old_vertexAtomsMapFlag);
  // setComputeNormals(old_computeNormals);

  accurateTriangulation = old;

  isAvailableScalarField = false;

  logging::log<logging::level::info>("ok!");

  deleteMatrix3D<double>(delphi->nx, delphi->ny, scalarField);

  char buff[BUFLEN];
  // strcpy(buff,"triangulatedSurf_fixed.off");
  strcpy(buff, "blobby.off");

  // the triangulation tools of surface where used. Now they are cleaned

  ////////////// reset surface status ////////////////////////////
  if (triList.size() > 0) {
    vector<int *>::iterator it;
    for (it = triList.begin(); it != triList.end(); it++)
      deleteVector<int>((*it));
  }

  if (vertList.size() > 0) {
    vector<double *>::iterator it;
    for (it = vertList.begin(); it != vertList.end(); it++)
      deleteVector<double>((*it));
  }

  triList.clear();
  vertList.clear();
  totalSurfaceArea = 0;
  totalVolume = 0;
  //////////////////////////////////////////////////////////

  // load its mesh
  load(buff);
  // the rest continue as a mesh surface
  // bool flag = MeshSurface::build();
  MeshSurface::printSummary();
  return true;
}

void BlobbySurface::printSummary() {
  logging::log<logging::level::info>("Blobbyness value {}", getBlobbyness());
  logging::log<logging::level::info>("Cut-off distance {}[A]", cutoff);
}
