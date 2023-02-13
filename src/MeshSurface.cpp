
#include <MeshSurface.h>
#include <globals.h>
#include <logging.h>

void MeshSurface::clear() {
  if (faceMatrix != NULL)
    deleteMatrix2D<int>(numTriangles, faceMatrix);

  if (vertMatrix != NULL)
    deleteMatrix2D<double>(numVertexes, vertMatrix);

  if (vertexTrianglesList != NULL) {
    for (int i = 0; i < numVertexes; i++)
      delete vertexTrianglesList[i];
    delete vertexTrianglesList;
  }

  if (vertNormals != NULL)
    deleteMatrix2D<double>(numVertexes, vertNormals);

  if (gridTriangleMap != NULL)
    free(gridTriangleMap);

  if (gridTriangleMap2D != NULL)
    free(gridTriangleMap2D);

  if (x != NULL)
    deleteVector<double>(x);
  if (y != NULL)
    deleteVector<double>(y);
  if (z != NULL)
    deleteVector<double>(z);
  if (ind != NULL)
    deleteMatrix3D<unsigned short>(nx, ny, ind);
  if (planes != NULL)
    deleteMatrix2D<double>(numTriangles, planes);
}
MeshSurface::~MeshSurface() { clear(); }

void MeshSurface::init() {
  faceMatrix = NULL;
  vertMatrix = NULL;
  vertNormals = NULL;
  gridTriangleMap = NULL;
  gridTriangleMap2D = NULL;
  vertexTrianglesList = NULL;
  planes = NULL;
  providesAnalyticalNormals = false;
  ind = NULL;
  ind_2d = NULL;
  x = NULL;
  y = x;
  z = y;
  MAX_TRIANGLES = 250;
  AUX_GRID_DIM = 100;
  AUX_GRID_DIM_2D = 100;
  MAX_TRIANGLES_2D = (MAX_TRIANGLES * AUX_GRID_DIM_2D);
  surfType = GENERIC_SURFACE;
}
MeshSurface::MeshSurface(DelPhiShared *ds) : Surface() {
  init();
  // set environment
  delphi = ds;
}

void MeshSurface::init(ConfigurationOP cf) {
  string fname =
      cf->sfname; // cf->read<string>( "Surface_File_Name", "mesh.off" );
  unsigned int maxMeshDim = cf->maxMeshDim;
  unsigned int maxMeshPatches = cf->maxMeshPatches;
  unsigned int maxMeshDim2D = cf->maxMeshDim2D;
  unsigned int maxMeshPatches2D = cf->maxMeshPatches2D;
  int numMSMSFiles = cf->NumMSMSfiles;

  setAuxGrid(maxMeshDim, maxMeshPatches);
  setAuxGrid2D(maxMeshDim2D, maxMeshPatches2D);
  // Set up inside value
  inside = 5;

  // get extension
  fname.size();
  string ext("");
  int ind = fname.size();

  while (fname[ind] != '.') {
    ext = fname.substr(ind, fname.size());
    ind--;
  }

  // get the type of surface
  if (!ext.compare("off")) {
    // Load Surface
    bool load_f = load((char *)fname.c_str());
    if (!load_f) {
      logging::log<logging::level::err>("Cannot load {}", fname);
#ifdef PYTHON
      throw std::exception();
#else
      exit(-1);
#endif
    }
  } else if (!ext.compare("vert") || !ext.compare("face")) {
    // Load MSMS surface[s]
    bool load_f = loadMSMS((char *)fname.substr(0, fname.size() - 5).c_str(),
                           numMSMSFiles);
    if (!load_f) {
      logging::log<logging::level::err>("Cannot load {}", fname);
#ifdef PYTHON
      throw std::exception();
#else
      exit(-1);
#endif
    }
    char cc[BUFLEN];
    strcpy(cc, "msms.off");
    save(cc);
  } else {
    logging::log<logging::level::err>(
        "Unknown surface type, please check extension/type");
#ifdef PYTHON
    throw std::exception();
#else
    exit(-1);
#endif
  }
}

MeshSurface::MeshSurface() : Surface() {
  init();
  delphi = NULL;
}

MeshSurface::MeshSurface(ConfigurationOP cf, DelPhiShared *ds) : Surface(cf) {
  init();
  init(cf);
  // set environment
  delphi = ds;
}

void MeshSurface::preProcessPanel() {
  unsigned int igrid = delphi->nx;

  // cannot have an auxiliary grid smaller than that of delphi.
  // all the rest of the code is based on this assumption

  // auxiliary grid is consistent with delphi grid in order to speed-up tracing
  int gridMul = 1;
  while (igrid > AUX_GRID_DIM_2D) {
    gridMul += 2;
    int digrid = delphi->nx;
    while (1) {
      // get nearest odd multiple
      int fixedPoint = (digrid + (gridMul - 1)) / gridMul;
      igrid = fixedPoint * gridMul;
      if (igrid % 2 == 0)
        digrid = igrid + 1;
      else {
        igrid = fixedPoint;
        break;
      }
    }
  }
  // auxiliary scale
  scale_2d = delphi->scale / ((double)gridMul);

  // logging::log<logging::level::info>("Auxiliary grid is " << igrid);

  xmin_2d = delphi->baricenter[0] - (igrid - 1) / (2 * scale_2d);
  ymin_2d = delphi->baricenter[1] - (igrid - 1) / (2 * scale_2d);
  zmin_2d = delphi->baricenter[2] - (igrid - 1) / (2 * scale_2d);

  xmax_2d = delphi->baricenter[0] + (igrid - 1) / (2 * scale_2d);
  ymax_2d = delphi->baricenter[1] + (igrid - 1) / (2 * scale_2d);
  zmax_2d = delphi->baricenter[2] + (igrid - 1) / (2 * scale_2d);

  side_2d = 1 / scale_2d;
  nx_2d = igrid;
  ny_2d = igrid;
  nz_2d = igrid;

  if (gridTriangleMap2D != NULL)
    deleteVector<int>(gridTriangleMap2D);

  if (ind_2d != NULL)
    deleteMatrix2D<unsigned int>(last_rows_ind, ind_2d);

  if (panel == 0) {
    last_rows_ind = ny_2d;
    last_cols_ind = nz_2d;
  } else if (panel == 1) {
    last_rows_ind = nx_2d;
    last_cols_ind = ny_2d;
  } else {
    last_rows_ind = nx_2d;
    last_cols_ind = nz_2d;
  }

  gridTriangleMap2D =
      allocateVector<int>(last_rows_ind * last_cols_ind * MAX_TRIANGLES_2D);

  ind_2d = allocateMatrix2D<unsigned int>(last_rows_ind, last_cols_ind);
  for (int i = 0; i < last_rows_ind; i++)
    for (int j = 0; j < last_cols_ind; j++)
      ind_2d[i][j] = 0;

  // build a bounding box for each triangle and map it to
  // the auxiliary grid
  unsigned int max_t = 0;

  if (vertMatrix == NULL || faceMatrix == NULL) {
    logging::log<logging::level::warn>("Cannot get surface without a mesh!");
    return;
  }

  double *p[3];
  for (int it = 0; it < numTriangles; it++) {
    // triangle points
    p[0] = vertMatrix[faceMatrix[it][0]];
    p[1] = vertMatrix[faceMatrix[it][1]];
    p[2] = vertMatrix[faceMatrix[it][2]];

    // compute bounding box object
    double downx = INFINITY;
    double downy = -INFINITY;
    double downz = INFINITY;

    double upx = -INFINITY;
    double upy = INFINITY;
    double upz = -INFINITY;

    for (unsigned int pind = 0; pind < 3; pind++) {
      double *pp = p[pind];
      downx = MIN(downx, pp[0]);
      downy = MAX(downy, pp[1]);
      downz = MIN(downz, pp[2]);

      upx = MAX(upx, pp[0]);
      upy = MIN(upy, pp[1]);
      upz = MAX(upz, pp[2]);
    }

    // resolve which are the grid cubes that
    // are cut by the bounding box object. These
    // cubes see the mixed cell bounding box
    int ix_start = (int)rintp((downx - xmin_2d) / side_2d);
    int iy_start = (int)rintp((downy - ymin_2d) / side_2d);
    int iz_start = (int)rintp((downz - zmin_2d) / side_2d);

    int ix_end = (int)rintp((upx - xmin_2d) / side_2d);
    int iy_end = (int)rintp((upy - ymin_2d) / side_2d);
    int iz_end = (int)rintp((upz - zmin_2d) / side_2d);

    if (ix_start < 0)
      ix_start = 0;
    if (iy_start < 0)
      iy_start = 0;
    if (iz_start < 0)
      iz_start = 0;

    if (ix_end < 0)
      ix_end = 0;
    if (iy_end < 0)
      iy_end = 0;
    if (iz_end < 0)
      iz_end = 0;

    if (ix_end >= nx_2d)
      ix_end = nx_2d - 1;
    if (iy_end >= ny_2d)
      iy_end = ny_2d - 1;
    if (iz_end >= nz_2d)
      iz_end = nz_2d - 1;

    if (ix_start >= nx_2d)
      ix_start = nx_2d - 1;
    if (iy_start >= ny_2d)
      iy_start = ny_2d - 1;
    if (iz_start >= nz_2d)
      iz_start = nz_2d - 1;

    // map the points into the 2D matrices of each plane

    if (panel == 0) {
      // plane yz
      for (int iz = iz_start; iz <= iz_end; iz++) {
        for (int iy = iy_start; iy >= iy_end; iy--) {
          if (ind_2d[iy][iz] > max_t)
            max_t = ind_2d[iy][iz];

          if (ind_2d[iy][iz] >= MAX_TRIANGLES_2D) {
            logging::log<logging::level::err>(
                "Number of triangles is superior to maximum allowed, please "
                "increase Max_mesh_patches_per_auxiliary_grid_2d_cell");
#ifdef PYTHON
            throw std::exception();
#else
            exit(-1);
#endif
          }
          GRID_TRIANGLEMAP_2D(iy, iz, (ind_2d[iy][iz]), ny_2d, nz_2d) = it;
          ind_2d[iy][iz]++;
        }
      }
    }

    else if (panel == 1) {
      // plane xy
      for (int ix = ix_start; ix <= ix_end; ix++) {
        for (int iy = iy_start; iy >= iy_end; iy--) {
          if (ind_2d[ix][iy] > max_t)
            max_t = ind_2d[ix][iy];

          if (ind_2d[ix][iy] >= MAX_TRIANGLES_2D) {
            logging::log<logging::level::err>(
                "Number of triangles is superior to maximum allowed, please "
                "increase Max_mesh_patches_per_auxiliary_grid_2d_cell");
#ifdef PYTHON
            throw std::exception();
#else
            exit(-1);
#endif
          }
          GRID_TRIANGLEMAP_2D(ix, iy, (ind_2d[ix][iy]), nx_2d, ny_2d) = it;
          ind_2d[ix][iy]++;
        }
      }
    } else {
      // plane xz
      for (int ix = ix_start; ix <= ix_end; ix++) {
        for (int iz = iz_start; iz <= iz_end; iz++) {
          if (ind_2d[ix][iz] > max_t)
            max_t = ind_2d[ix][iz];

          if (ind_2d[ix][iz] >= MAX_TRIANGLES_2D) {
            logging::log<logging::level::err>(
                "Number of triangles is superior to maximum allowed, please "
                "increase Max_mesh_patches_per_auxiliary_grid_2d_cell");
#ifdef PYTHON
            throw std::exception();
#else
            exit(-1);
#endif
          }
          GRID_TRIANGLEMAP_2D(ix, iz, (ind_2d[ix][iz]), nx_2d, nz_2d) = it;
          ind_2d[ix][iz]++;
        }
      }
    }
  }
}

void MeshSurface::preProcessTriangles() {
  if (planes != NULL)
    deleteMatrix2D<double>(numTriangles, planes);

  planes = allocateMatrix2D<double>(numTriangles, 4);

  // compute planes from triangles. Planes are normalized in a further step
  // because if normals estimation is needed it is usually more numerically
  // stable to average the plane normals before normalization

  double *p[3];
  for (int it = 0; it < numTriangles; it++) {
    p[0] = vertMatrix[faceMatrix[it][0]];
    p[1] = vertMatrix[faceMatrix[it][1]];
    p[2] = vertMatrix[faceMatrix[it][2]];

    // compute plane in planes[it]
    plane3points(p[0], p[1], p[2], planes[it]);
  }

  // compute normals for simple triangular meshes. Estimate vertex normal by
  // averaging the normals of connected triangles. Does not compute normals if
  // vertex normals are provided, such as in MSMS files.
  if (computeNormals) {
    if (vertNormals != NULL)
      deleteMatrix2D<double>(numVertexes, vertNormals);

    vertNormals = allocateMatrix2D<double>(numVertexes, 3);

    for (int iv = 0; iv < numVertexes; iv++) {
      // p[0] = vertMatrix[iv];

      double meanNormal[3];

      meanNormal[0] = 0;
      meanNormal[1] = 0;
      meanNormal[2] = 0;
      vector<int>::iterator it;

      for (it = vertexTrianglesList[iv]->begin();
           it != vertexTrianglesList[iv]->end(); it++) {
        int triangleID = *it;
        ADD(meanNormal, meanNormal, planes[triangleID]);
      }

      double vlen = (double)vertexTrianglesList[iv]->size();
      vertNormals[iv][0] = meanNormal[0] / vlen;
      vertNormals[iv][1] = meanNormal[1] / vlen;
      vertNormals[iv][2] = meanNormal[2] / vlen;

      // normalization
      double norm = sqrt(DOT(vertNormals[iv], vertNormals[iv]));
      vertNormals[iv][0] /= norm;
      vertNormals[iv][1] /= norm;
      vertNormals[iv][2] /= norm;
    }
  }

  for (int it = 0; it < numTriangles; it++) {
    p[0] = vertMatrix[faceMatrix[it][0]];
    p[1] = vertMatrix[faceMatrix[it][1]];
    p[2] = vertMatrix[faceMatrix[it][2]];

    // normalization of the plane normal
    double norm = sqrt(DOT(planes[it], planes[it]));
    planes[it][0] /= norm;
    planes[it][1] /= norm;
    planes[it][2] /= norm;
    planes[it][3] = -(p[0][0] * (p[1][1] * p[2][2] - p[2][1] * p[1][2]) +
                      p[1][0] * (p[2][1] * p[0][2] - p[0][1] * p[2][2]) +
                      p[2][0] * (p[0][1] * p[1][2] - p[1][1] * p[0][2]));
    planes[it][3] /= norm;
  }
}
bool MeshSurface::buildAuxiliaryGrid() {
  if (faceMatrix == NULL || vertMatrix == NULL) {
    logging::log<logging::level::warn>(
        "Cannot get surface without a loaded mesh!");
    return false;
  }

  // Perform pre-processing to speed-up intersections and projections.
  // Compute bounding box for each triangle in the mesh, and map each
  // bounding box to the proper grid point.
  // This routine uses an auxiliary grid. Delphi grid is not directly used
  // because it may be too memory consuming. To speed up computations
  // a maximal number of triangles is allowed in each auxiliary grid cell
  // the macro is MAX_TRIANGLES in MeshSurface.h.
  // The macro AUX_GRID_DIM sets the maximally allowed grid size.
  unsigned int igrid = delphi->nx;
  // cannot have an auxiliary grid smaller than that of delphi.
  // all the rest of the code is based on this assumption

  int gridMul = 1;
  while (igrid > AUX_GRID_DIM) {
    gridMul += 2;
    int digrid = delphi->nx;
    while (1) {
      // get nearest odd multiple
      int fixedPoint = (digrid + (gridMul - 1)) / gridMul;
      igrid = fixedPoint * gridMul;
      // printf("\n Tenative Aux grid %d",igrid);
      if (igrid % 2 == 0) {
        digrid = igrid + 1;
        // printf("\n Rejecting grid %d",igrid);
      } else {
        igrid = fixedPoint;
        // printf("\n Accepting grid %d",igrid);
        break;
      }
    }
  }
  // auxiliary scale
  scale = delphi->scale / ((double)gridMul);

  logging::log<logging::level::info>(
      "Auxiliary grid is {} auxiliary scale is {} grid mul {}", igrid, scale,
      gridMul);

  xmin = delphi->baricenter[0] - (igrid - 1) / (2 * scale);
  ymin = delphi->baricenter[1] - (igrid - 1) / (2 * scale);
  zmin = delphi->baricenter[2] - (igrid - 1) / (2 * scale);

  xmax = delphi->baricenter[0] + (igrid - 1) / (2 * scale);
  ymax = delphi->baricenter[1] + (igrid - 1) / (2 * scale);
  zmax = delphi->baricenter[2] + (igrid - 1) / (2 * scale);

  ////////////// Allocate memory for the grid, planes and maps
  /////////////////////////

  if (x != NULL)
    deleteVector<double>(x);
  if (y != NULL)
    deleteVector<double>(y);
  if (z != NULL)
    deleteVector<double>(z);

  x = allocateVector<double>(igrid);
  y = allocateVector<double>(igrid);
  z = allocateVector<double>(igrid);

  if (ind != NULL)
    deleteMatrix3D<unsigned short>(nx, ny, ind);

  side = 1 / scale;
  nx = igrid;
  ny = igrid;
  nz = igrid;

  ind = allocateMatrix3D<unsigned short>(nx, ny, nz);
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        ind[i][j][k] = 0;

  for (int i = 0; i < nx; i++)
    x[i] = xmin + i * side;

  for (int i = 0; i < ny; i++)
    y[i] = ymin + i * side;

  for (int i = 0; i < nz; i++)
    z[i] = zmin + i * side;

  if (gridTriangleMap != NULL)
    free(gridTriangleMap);

  logging::log<logging::level::info>(
      "Allocating {} MB for the auxiliary grid...",
      (nx * ny * nz * MAX_TRIANGLES) * sizeof(int) / 1024.0 / 1024.0);
  gridTriangleMap = (int *)malloc(sizeof(int) * nx * ny * nz * MAX_TRIANGLES);

  //////////////////////////////////////////////////////////////////////////

  double *p[3];
  int max_t = 0;

  // build a bounding box for each triangle and map it to
  // the auxiliary grid
  for (int it = 0; it < numTriangles; it++) {
    p[0] = vertMatrix[faceMatrix[it][0]];
    p[1] = vertMatrix[faceMatrix[it][1]];
    p[2] = vertMatrix[faceMatrix[it][2]];

    // compute bounding box object
    double downx = 1e20;
    double downy = -1e20;
    double downz = 1e20;

    double upx = -1e20;
    double upy = 1e20;
    double upz = -1e20;

    for (unsigned int pind = 0; pind < 3; pind++) {
      double *pp = p[pind];
      downx = MIN(downx, pp[0]);
      downy = MAX(downy, pp[1]);
      downz = MIN(downz, pp[2]);

      upx = MAX(upx, pp[0]);
      upy = MIN(upy, pp[1]);
      upz = MAX(upz, pp[2]);
    }

    // resolve which are the grid cubes that
    // are cut by the bounding object. These
    // cubes see the triangle bounding box
    int ix_start = (int)rintp((downx - xmin) / side);
    int iy_start = (int)rintp((downy - ymin) / side);
    int iz_start = (int)rintp((downz - zmin) / side);

    int ix_end = (int)rintp((upx - xmin) / side);
    int iy_end = (int)rintp((upy - ymin) / side);
    int iz_end = (int)rintp((upz - zmin) / side);

    // printf("\n %f %f %f %f %f %f %d %d %d %d %d
    // %d",downx,downy,downz,upx,upy,upz,ix_start,iy_start,iz_start,ix_end,iy_end,iz_end);
    if (ix_start < 0 || iy_start < 0 || iz_start < 0 || ix_end >= nx ||
        iy_end >= ny || iz_end >= nz) {
      logging::log<logging::level::err>(
          "Triangle bounding box out of grid. Increase perfil");
#ifdef PYTHON
      throw std::exception();
#else
      exit(-1);
#endif
    }

    bool assigned = false;
    for (int iz = iz_start; iz <= iz_end; iz++)
      for (int iy = iy_start; iy >= iy_end; iy--)
        for (int ix = ix_start; ix <= ix_end; ix++) {
          if (ix < 0 || iy < 0 || iz < 0)
            return false;

          if (ind[ix][iy][iz] > max_t)
            max_t = ind[ix][iy][iz];

          if (ind[ix][iy][iz] >= MAX_TRIANGLES) {
            logging::log<logging::level::err>(
                "Number of triangles is superior to maximum allowed, please "
                "increase Max_mesh_patches_per_auxiliary_grid_cell");
#ifdef PYTHON
            throw std::exception();
#else
            exit(-1);
#endif
          }
          // gridTriangleMap[ix][iy][iz][ind[ix][iy][iz]]=it;
          GRIDTRIANGLEMAP(ix, iy, iz, (ind[ix][iy][iz]), nx, ny, nz) = it;
          ind[ix][iy][iz]++;
          assigned = true;
        }

    if (!assigned)
      logging::log<logging::level::err>("Not assigned triangle!");
  }

  logging::log<logging::level::info>("Max triangles per cell -> {}", max_t);

  return true;
}

/** save in OFF format*/
bool MeshSurface::save(char *fileName) {
  ofstream fout;
  fout.open(fileName, ios::out);

  logging::log<logging::level::info>("Writing mesh in OFF file format in {}",
                                     fileName);

  if (fout.fail()) {
    logging::log<logging::level::warn>("Cannot write file {}", fileName);
    return false;
  }

  if (vertMatrix == NULL || faceMatrix == NULL) {
    logging::log<logging::level::warn>("Cannot write null mesh!");
    return false;
  }

  fout << "OFF" << endl;
  fout << "# File created by " << PROGNAME << " version " << VERSION << endl;
  fout << endl;
  fout << numVertexes << " " << numTriangles << " 0" << endl;
  for (int i = 0; i < numVertexes; i++)
    fout << vertMatrix[i][0] << " " << vertMatrix[i][1] << " "
         << vertMatrix[i][2] << endl;
  for (int i = 0; i < numTriangles; i++)
    fout << "3 " << faceMatrix[i][0] << " " << faceMatrix[i][1] << " "
         << faceMatrix[i][2] << endl;

  return true;
}

bool MeshSurface::load(char *fileName) {
  int len = (int)strlen(fileName);
  if (len == 0) {
    logging::log<logging::level::warn>("Cannot load with empty file name!");
    return false;
  }
  string ss(fileName);
  size_t found;
  ss = toLowerCase(ss);
  found = ss.find(".off");

  logging::log<logging::level::info>("Loading mesh in file {}", fileName);

  bool exit = true;

  if (found != string::npos)
    exit = loadOFF(fileName);

  found = ss.find(".ply");
  if (found != string::npos)
    exit = loadPLY(fileName);

  if (!exit)
    return exit;
  else
    preProcessTriangles();

  return true;
}

bool MeshSurface::loadPLY(char *fname) {
  int format = 0, voh, foh, vph, fph;
  int nv, nt, i, j, i1, i2, i3, i4;
  float x, y, z;
  bool triangulate = 0;
  FILE *in;
  char keyword[64], formats[24], version[10];

  if ((in = fopen(fname, "rb")) == NULL) {
    logging::log<logging::level::err>("Can't open input ply file\n");
    return false;
  }

  if (strcmp(readLineFromFile(in), "ply")) {
    logging::log<logging::level::err>("Input doesn't seem a valid ply file.\n");
    return false;
  }
  if (sscanf(readLineFromFile(in), "%7s %24s %10s", keyword, formats, version) <
      3) {
    logging::log<logging::level::err>("Unexpected token or end of file!\n");
    return false;
  }
  if (strcmp(keyword, "format")) {
    logging::log<logging::level::err>("format definition expected!\n");
    return false;
  }

  if (!strcmp(formats, "ascii"))
    format = PLY_FORMAT_ASCII;
  else if (!strcmp(formats, "binary_little_endian"))
    format = PLY_FORMAT_BIN_L;
  else if (!strcmp(formats, "binary_big_endian"))
    format = PLY_FORMAT_BIN_B;
  else {
    logging::log<logging::level::err>("Unrecognized format '%s'\n", formats);
    return false;
  }

  nv = ply_parseElements(in, "vertex");
  vph = ply_getOverhead(in, format, "vertex");
  ply_checkVertexProperties(in);
  voh = ply_getOverhead(in, format, "vertex");
  nt = ply_parseElements(in, "face");
  fph = ply_getOverhead(in, format, "face");
  ply_checkFaceProperties(in);
  foh = ply_getOverhead(in, format, "face");

  if (!sscanf(readLineFromFile(in), "%64s ", keyword)) {
    logging::log<logging::level::err>("Unexpected token or end of file!\n");
    return false;
  }
  while (strcmp(keyword, "end_header"))
    if (!sscanf(readLineFromFile(in), "%64s ", keyword)) {
      logging::log<logging::level::err>("Unexpected token or end of file!\n");
      return false;
    }

  if (vertMatrix != NULL)
    deleteMatrix2D<double>(numVertexes, vertMatrix);

  if (vertexTrianglesList != NULL) {
    for (int i = 0; i < numVertexes; i++)
      delete vertexTrianglesList[i];
    delete vertexTrianglesList;
  }

  numVertexes = nv;
  vertMatrix = allocateMatrix2D<double>(numVertexes, 3);
  // vertexTrianglesList = new vector<int>*[numVertexes];
  vertexTrianglesList = allocateVector<vector<int> *>(numVertexes);

  for (int i = 0; i < numVertexes; i++)
    vertexTrianglesList[i] = new vector<int>();

  for (i = 0; i < nv; i++) {
    ply_readVCoords(in, format, vph, voh, &x, &y, &z);
    vertMatrix[i][0] = x;
    vertMatrix[i][1] = y;
    vertMatrix[i][2] = z;
  }

  if (faceMatrix != NULL)
    deleteMatrix2D<int>(numTriangles, faceMatrix);

  numTriangles = nt;
  faceMatrix = allocateMatrix2D<int>(numTriangles, 3);

  for (i = 0; i < nt; i++) {
    if (ply_readFIndices(in, format, fph, &i4, &i1, &i2, &i3)) {
      if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 3 || i4 >= 4 || i1 > (nv - 1) ||
          i2 > (nv - 1) || i3 > (nv - 1)) {
        logging::log<logging::level::err>(
            "\nloadPLY: Invalid index at face %d!\n", i);
        return false;
      }
      for (j = 3; j <= i4; j++) {
        if (i1 == i2 || i2 == i3 || i3 == i1)
          logging::log<logging::level::warn>(
              "\nloadPLY: Coincident indexes at triangle %d! Skipping.\n", i);
        else {
          faceMatrix[i][0] = i1;
          faceMatrix[i][1] = i2;
          faceMatrix[i][2] = i3;
          vertexTrianglesList[i1]->push_back(i);
          vertexTrianglesList[i2]->push_back(i);
          vertexTrianglesList[i3]->push_back(i);
        }
        i2 = i3;
        if (j < i4) {
          if (!ply_readAnotherFIndex(in, format, &i3)) {
            logging::log<logging::level::err>(
                "\nloadPLY: Couldn't read indexes for face # %d\n", i);
            return false;
          } else
            triangulate = 1;
        } else
          ply_readOverhead(in, format, foh);
      }
    } else
      logging::log<logging::level::err>(
          "\nloadPLY: Couldn't read indexes for face # %d\n", i);
  }

  fclose(in);
  if (checkDuplicatedVertices)
    checkDuplicates();
  printSummary();
  computeNormals = true;
  return true;
}

bool MeshSurface::loadOFF(char *fileName) {
  ifstream fin;
  fin.open(fileName, ios::in);
  int temp;
  char buffer[BUFLEN];

  if (fin.fail()) {
    logging::log<logging::level::warn>("Cannot read file {}", fileName);
    return false;
  }

  int tag = 0;
  int header = 0;

  if (faceMatrix != NULL)
    deleteMatrix2D<int>(numTriangles, faceMatrix);

  if (vertMatrix != NULL)
    deleteMatrix2D<double>(numVertexes, vertMatrix);

  if (vertexTrianglesList != NULL) {
    for (int i = 0; i < numVertexes; i++)
      delete vertexTrianglesList[i];
    deleteVector<vector<int> *>(vertexTrianglesList);
  }

  while (1) {
    fin.getline(buffer, BUFLEN);

    if (fin.eof())
      break;

    // skip empty lines
    if (strlen(buffer) <= 1)
      continue;

    if (tag == 0 && strncmp(buffer, "OFF", 1) == 0) {
      tag = 1;
      continue;
    }

    if (buffer[0] == '#') {
      logging::log<logging::level::info>("Mesh comment line: {}", buffer);
      continue;
    }

    // if OFF is read then the header can be read
    if (header == 0 && tag == 1) {
      sscanf(buffer, "%d %d %d", &numVertexes, &numTriangles, &temp);
      if (numVertexes <= 0 || numTriangles <= 0) {
        logging::log<logging::level::warn>(
            "Number of triangles or vertices <=0 !");
        return false;
      } else {
        faceMatrix = allocateMatrix2D<int>(numTriangles, 3);
        vertMatrix = allocateMatrix2D<double>(numVertexes, 3);
        // vertexTrianglesList = new vector<int>*[numVertexes];
        vertexTrianglesList = allocateVector<vector<int> *>(numVertexes);
        for (int i = 0; i < numVertexes; i++)
          vertexTrianglesList[i] = new vector<int>();
      }
      header = 1;
      break;
    }
  }

  if (header != 1 || tag != 1) {
    logging::log<logging::level::warn>(
        "Cannot read OFF header or number of vertices/triangles, stop reading");
    logging::log<logging::level::warn>("Tag {} Header {}", tag, header);
    return false;
  }
  // read vertices
  for (int i = 0; i < numVertexes; i++) {
    fin.getline(buffer, BUFLEN);
    if (strlen(buffer) <= 1) {
      i--;
      continue;
    }
    float v1, v2, v3;
    sscanf(buffer, "%f %f %f", &v1, &v2, &v3);
    vertMatrix[i][0] = v1;
    vertMatrix[i][1] = v2;
    vertMatrix[i][2] = v3;
  }

  for (int i = 0; i < numTriangles; i++) {
    fin.getline(buffer, BUFLEN);
    if (strlen(buffer) <= 1) {
      i--;
      continue;
    }

    sscanf(buffer, "%d %d %d %d", &temp, &(faceMatrix[i][0]),
           &(faceMatrix[i][1]), &(faceMatrix[i][2]));

    vertexTrianglesList[faceMatrix[i][0]]->push_back(i);
    vertexTrianglesList[faceMatrix[i][1]]->push_back(i);
    vertexTrianglesList[faceMatrix[i][2]]->push_back(i);

    if (temp != 3) {
      logging::log<logging::level::warn>("Non triangular mesh, stop reading");
      return false;
    }
  }

  if (checkDuplicatedVertices)
    checkDuplicates();
  printSummary();
  computeNormals = true;

  return true;
}

bool MeshSurface::loadMSMS(char *fileName, int numFiles) {
  ifstream fin, fin2;

  char buffer[BUFLEN];
  char baseName[BUFLEN];

  char currentFace[BUFLEN];
  char currentVert[BUFLEN];

  strcpy(baseName, fileName);

  int tempNumVertices = 0;
  int tempNumTri = 0;
  numVertexes = 0;
  numTriangles = 0;

  int lastVert = 0;
  int lastTri = 0;

  if (vertMatrix != NULL)
    deleteMatrix2D<double>(numVertexes, vertMatrix);

  if (vertNormals != NULL)
    deleteMatrix2D<double>(numVertexes, vertNormals);

  if (faceMatrix != NULL)
    deleteMatrix2D<int>(numTriangles, faceMatrix);

  // getting the total number of vertices and triangles
  for (int i = 0; i < numFiles; i++) {
    // main component
    if (i == 0) {
      snprintf(currentFace, sizeof(currentFace), "%s.face", baseName);
      snprintf(currentVert, sizeof(currentVert), "%s.vert", baseName);
    }
    // cavities
    else {
      snprintf(currentFace, sizeof(currentFace), "%s_%d.face", baseName, i);
      snprintf(currentVert, sizeof(currentVert), "%s_%d.vert", baseName, i);
    }

    fin.open(currentVert, ios::in);
    fin2.open(currentFace, ios::in);
    logging::log<logging::level::info>(
        "Getting num vertices and triangles in MSMS mesh files {}, {}",
        currentFace, currentVert);
    if (fin.fail() || fin2.fail()) {
      logging::log<logging::level::warn>("One or both MSMS files don't exist");
      return false;
    }

    // comment or empty
    bool comment = true;
    while (comment) {
      fin.getline(buffer, BUFLEN);
      if (buffer[0] != '#' && buffer[0] != '\0')
        comment = false;
    }
    sscanf(buffer, "%d", &tempNumVertices);
    numVertexes += tempNumVertices;
    fin.close();

    // reads num triangles
    comment = true;
    while (comment) {
      fin2.getline(buffer, BUFLEN);
      if (buffer[0] != '#' && buffer[0] != '\0')
        comment = false;
    }
    sscanf(buffer, "%d", &tempNumTri);
    numTriangles += tempNumTri;
    fin2.close();
  }

  logging::log<logging::level::info>("Total number of triangles {}",
                                     numTriangles);
  logging::log<logging::level::info>("Total number of vertices {}",
                                     numVertexes);

  vertMatrix = allocateMatrix2D<double>(numVertexes, 3);
  vertNormals = allocateMatrix2D<double>(numVertexes, 3);
  faceMatrix = allocateMatrix2D<int>(numTriangles, 3);

  // read all
  for (int i = 0; i < numFiles; i++) {
    // main component
    if (i == 0) {
      snprintf(currentFace, sizeof(currentFace), "%s.face", baseName);
      snprintf(currentVert, sizeof(currentVert), "%s.vert", baseName);
    }
    // cavities
    else {
      snprintf(currentFace, sizeof(currentFace), "%s_%d.face", baseName, i);
      snprintf(currentVert, sizeof(currentVert), "%s_%d.vert", baseName, i);
    }

    fin.open(currentVert, ios::in);
    fin2.open(currentFace, ios::in);
    logging::log<logging::level::info>("Loading MSMS mesh files {},{}",
                                       currentFace, currentVert);

    if (fin.fail() || fin2.fail()) {
      logging::log<logging::level::warn>("One or both MSMS files don't exist");
      return false;
    }

    bool comment = true;
    while (comment) {
      fin.getline(buffer, BUFLEN);
      if (buffer[0] != '#' && buffer[0] != '\0')
        comment = false;
    }

    sscanf(buffer, "%d", &tempNumVertices);

    for (int i = lastVert; i < lastVert + tempNumVertices; i++) {
      fin.getline(buffer, BUFLEN);
      // skip empty lines
      if (strlen(buffer) > 1) {
        float a, b, c, d, e, f;
        sscanf(buffer, "%f %f %f %f %f %f", &a, &b, &c, &d, &e, &f);
        vertMatrix[i][0] = a;
        vertMatrix[i][1] = b;
        vertMatrix[i][2] = c;
        vertNormals[i][0] = d;
        vertNormals[i][1] = e;
        vertNormals[i][2] = f;
      } else
        i--;
    }
    fin.close();

    comment = true;
    while (comment) {
      fin2.getline(buffer, BUFLEN);
      if (buffer[0] != '#' && buffer[0] != '\0')
        comment = false;
    }

    sscanf(buffer, "%d", &tempNumTri);

    for (int i = lastTri; i < lastTri + tempNumTri; i++) {
      fin2.getline(buffer, BUFLEN);
      // skip empty lines
      if (strlen(buffer) > 1) {
        sscanf(buffer, "%d %d %d", &(faceMatrix[i][0]), &(faceMatrix[i][1]),
               &(faceMatrix[i][2]));
        faceMatrix[i][0]--;
        faceMatrix[i][1]--;
        faceMatrix[i][2]--;
        faceMatrix[i][0] += lastVert;
        faceMatrix[i][1] += lastVert;
        faceMatrix[i][2] += lastVert;
      } else
        i--;
    }
    fin2.close();

    lastTri += tempNumTri;
    lastVert += tempNumVertices;
  }

  if (checkDuplicatedVertices)
    checkDuplicates();
  printSummary();
  computeNormals = false;

  preProcessTriangles();

  return true;
}

bool MeshSurface::checkDuplicates() {
  if (faceMatrix == NULL || vertMatrix == NULL) {
    logging::log<logging::level::warn>(
        "Cannot check mesh without a loaded mesh!");
    return false;
  }

  bool ret = true;

  // check for duplicated vertices
  set<pair<double, pair<double, double>>> checkV;
  set<pair<double, pair<double, double>>>::iterator stv;

  for (int i = 0; i < numVertexes; i++) {
    pair<double, pair<double, double>> triplet;
    triplet.first = vertMatrix[i][0];
    triplet.second.first = vertMatrix[i][1];
    triplet.second.second = vertMatrix[i][2];
    stv = checkV.find(triplet);
    if (stv == checkV.end()) {
      // ok
    }
    // duplicated vertex
    else {
      logging::log<logging::level::warn>("Duplicated vertex detected! {} {} {}",
                                         triplet.first, triplet.second.first,
                                         triplet.second.second);
      ret = false;
    }

    checkV.insert(triplet);
  }

  // check for duplicated triangles
  set<pair<int, pair<int, int>>> checkT;
  set<pair<int, pair<int, int>>>::iterator st;

  for (int i = 0; i < numTriangles; i++) {
    pair<int, pair<int, int>> triplet;
    triplet.first = faceMatrix[i][0];
    triplet.second.first = faceMatrix[i][1];
    triplet.second.second = faceMatrix[i][2];
    st = checkT.find(triplet);
    if (st == checkT.end()) {
      // ok
    }
    // duplicated triangle
    else {
      logging::log<logging::level::warn>(
          "Duplicated triangle detected! {} {} {}", triplet.first,
          triplet.second.first, triplet.second.second);
      ret = false;
    }

    checkT.insert(triplet);
  }

  return ret;
}

void MeshSurface::printSummary() {
  // logging::log<logging::level::info>("Summary of the triangulated mesh:");
  if (faceMatrix == NULL || vertMatrix == NULL) {
    logging::log<logging::level::warn>("Mesh not loaded!");
  } else {
    logging::log<logging::level::info>("Number of loaded vertices -> {}",
                                       numVertexes);
    logging::log<logging::level::info>("Number of loaded triangles -> {}",
                                       numTriangles);
  }
}

inline void
MeshSurface::getRayIntersection(double pa[3], double pb[3],
                                vector<pair<double, double *>> &intersections,
                                int thdID, bool computeNormals) {
  double dir[3];
  double *p[3];
  double aposx = pa[0], aposy = pa[1], aposz = pa[2];
  dir[0] = pb[0] - aposx;
  dir[1] = pb[1] - aposy;
  dir[2] = pb[2] - aposz;
  double apos[3];
  apos[0] = aposx;
  apos[1] = aposy;
  apos[2] = aposz;

  int ixstart = (int)rintp((pa[0] - delphi->xmin) / delphi->side);
  int iystart = (int)rintp((pa[1] - delphi->ymin) / delphi->side);
  int izstart = (int)rintp((pa[2] - delphi->zmin) / delphi->side);

  int iyaux, izaux, ixaux;

  if (panel == 0) {
    int t1 = (int)rintp((pa[1] - ymin_2d) / side_2d);
    int t2 = (int)rintp((pa[2] - zmin_2d) / side_2d);

    // get auxiliary grid index from delphi grid point (or from starting ray
    // point)
    iyaux = t1;
    izaux = t2;
    int iter;

    int triIndex = ind_2d[iyaux][izaux];
    for (iter = 0; iter < triIndex; iter++) {
      int it = GRID_TRIANGLEMAP_2D(iyaux, izaux, iter, ny_2d, nz_2d);
      p[0] = vertMatrix[faceMatrix[it][0]];
      p[1] = vertMatrix[faceMatrix[it][1]];
      p[2] = vertMatrix[faceMatrix[it][2]];
      double t, u, v;
      int res = 0;
      res = intersect_triangle(apos, dir, p[0], p[1], p[2], &t, &u, &v);

      double proj_point[3];
      if (res) {
        proj_point[0] = apos[0] + dir[0] * t;
        proj_point[1] = apos[1];
        proj_point[2] = apos[2];
        intersections.push_back(pair<double, double *>(t, planes[it]));
      }
    }
  } else if (panel == 1) {
    int t1 = (int)rintp((pa[1] - ymin_2d) / side_2d);
    int t2 = (int)rintp((pa[0] - xmin_2d) / side_2d);

    // get auxiliary grid index
    iyaux = t1;
    ixaux = t2;

    int iter;
    int triIndex = ind_2d[ixaux][iyaux];
    for (iter = 0; iter < triIndex; iter++) {
      int it = GRID_TRIANGLEMAP_2D(ixaux, iyaux, iter, nx_2d, ny_2d);
      p[0] = vertMatrix[faceMatrix[it][0]];
      p[1] = vertMatrix[faceMatrix[it][1]];
      p[2] = vertMatrix[faceMatrix[it][2]];

      double t, u, v;
      int res = intersect_triangle(apos, dir, p[0], p[1], p[2], &t, &u, &v);

      if (res) {
        double proj_point[3];
        proj_point[0] = apos[0];
        proj_point[1] = apos[1];
        proj_point[2] = apos[2] + dir[2] * t;
        intersections.push_back(pair<double, double *>(t, planes[it]));
      }
    }
  } else {

    int t1 = (int)rintp((pa[0] - xmin_2d) / side_2d);
    int t2 = (int)rintp((pa[2] - zmin_2d) / side_2d);

    // get auxiliary grid index
    ixaux = t1;
    izaux = t2;

    int iter;
    int triIndex = ind_2d[ixaux][izaux];
    for (iter = 0; iter < triIndex; iter++) {
      int it = GRID_TRIANGLEMAP_2D(ixaux, izaux, iter, nx_2d, nz_2d);
      p[0] = vertMatrix[faceMatrix[it][0]];
      p[1] = vertMatrix[faceMatrix[it][1]];
      p[2] = vertMatrix[faceMatrix[it][2]];

      double t, u, v;
      int res = intersect_triangle(apos, dir, p[0], p[1], p[2], &t, &u, &v);
      if (res == 1) {
        double proj_point[3];
        proj_point[0] = apos[0];
        proj_point[1] = apos[1] + dir[1] * t;
        proj_point[2] = apos[2];
        intersections.push_back(pair<double, double *>(t, planes[it]));
      }
    }
  }

  if (intersections.size() > 0)
    sort(intersections.begin(), intersections.end(), compKeepIndex);
}

// returns 1 in case of intersection
int MeshSurface::intersect_triangle(double orig[3], double dir[3],
                                    double vert0[3], double vert1[3],
                                    double vert2[3], double *t, double *u,
                                    double *v) {
  double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  double det, inv_det;

  /* find vectors for two edges sharing vert0 */
  SUB(edge1, vert1, vert0);
  SUB(edge2, vert2, vert0);

  /* begin calculating determinant - also used to calculate U parameter */
  CROSS(pvec, dir, edge2);

  /* if determinant is near zero, ray lies in plane of triangle */
  det = DOT(edge1, pvec);

#ifdef TEST_CULL /* define TEST_CULL if culling is desired */
  if (det < EPS)
    return 0;

  /* calculate distance from vert0 to ray origin */
  SUB(tvec, orig, vert0);
  /* calculate U parameter and test bounds */
  *u = DOT(tvec, pvec);
  if (*u < 0.0 || *u > det)
    return 0;

  /* prepare to test V parameter */
  CROSS(qvec, tvec, edge1);

  /* calculate V parameter and test bounds */
  *v = DOT(dir, qvec);
  if (*v < 0.0 || *u + *v > det)
    return 0;

  /* calculate t, scale parameters, ray intersects triangle */
  *t = DOT(edge2, qvec);
  inv_det = 1.0 / det;
  *t *= inv_det;
  *u *= inv_det;
  *v *= inv_det;
#else /* the non-culling branch */
  if (det > -EPS && det < EPS)
    return 0;

  inv_det = 1.0 / det;

  /* calculate distance from vert0 to ray origin */
  SUB(tvec, orig, vert0);
  /* calculate U parameter and test bounds */
  *u = DOT(tvec, pvec) * inv_det;
  if (*u < 0.0 || *u > 1.0)
    return 0;

  /* prepare to test V parameter */
  CROSS(qvec, tvec, edge1);

  /* calculate V parameter and test bounds */
  *v = DOT(dir, qvec) * inv_det;
  if (*v < 0.0 || *u + *v > 1.0)
    return 0;

  /* calculate t, ray intersects triangle */
  *t = DOT(edge2, qvec) * inv_det;
#endif

  return 1;
}

bool MeshSurface::point2triangle(double P[3], double A[3], double B[3],
                                 double C[3], double w[4], double *proj,
                                 double *dist, double *normal, int planeID,
                                 int va, int vb, int vc) {
  bool flag;
  // project to plane
  point2plane(P, w, dist, proj);
  // test triangle
  flag = inTriangle(proj, A, B, C);

  // if not in triangle return the projection
  // on the nearest edge to the query point
  // the normal in this case is estimated as the convex combination
  // of the vertex normals of the vertexes of the edge where the projection
  // lies.
  if (!flag) {
    double t1[3], t2[3], l2, d[3];
    double res1[3], res2[3], res3[3];
    double tt1, tt2, tt3;

    // A,B edge	projection
    SUB(t2, B, A);
    l2 = DOT(t2, t2);
    SUB(t1, proj, A);
    double t = (DOT(t1, t2)) / l2;

    if (t < 0.0) {
      SUB(t1, P, A);
      d[0] = sqrt(DOT(t1, t1)); // Beyond the first end of the segment
      ASSIGN(res1, A);
      // saturation
      t = 0.0;
    } else if (t > 1.0) {
      SUB(t1, P, B);
      d[0] = sqrt(DOT(t1, t1)); // Beyond the second end of the segment
      ASSIGN(res1, B);
      // saturation
      t = 1.0;
    } else {
      ADD_MUL(res1, A, t2, t);
      SUB(t1, res1, P);
      d[0] = sqrt(DOT(t1, t1));
    }
    tt1 = t;

    // A,C edge	projection
    SUB(t2, C, A);
    l2 = DOT(t2, t2);
    SUB(t1, proj, A);
    t = (DOT(t1, t2)) / l2;

    if (t < 0.0) {
      SUB(t1, P, A);
      d[1] = sqrt(DOT(t1, t1)); // Beyond the first end of the segment
      ASSIGN(res2, A);
      t = 0.0;
    } else if (t > 1.0) {
      SUB(t1, P, C);
      d[1] = sqrt(DOT(t1, t1)); // Beyond the second end of the segment
      ASSIGN(res2, C);
      t = 1.0;
    } else {
      ADD_MUL(res2, A, t2, t);
      SUB(t1, res2, P);
      d[1] = sqrt(DOT(t1, t1));
    }
    tt2 = t;

    // B,C edge	projection
    SUB(t2, C, B);
    l2 = DOT(t2, t2);
    SUB(t1, proj, B);
    t = (DOT(t1, t2)) / l2;

    if (t < 0.0) {
      SUB(t1, P, B);
      d[2] = sqrt(DOT(t1, t1)); // Beyond the first end of the segment
      ASSIGN(res3, B);
      t = 0.0;
    } else if (t > 1.0) {
      SUB(t1, P, C);
      d[2] = sqrt(DOT(t1, t1)); // Beyond the second end of the segment
      ASSIGN(res3, C);
      t = 1.0;
    } else {
      ADD_MUL(res3, B, t2, t);
      SUB(t1, res3, P);
      d[2] = sqrt(DOT(t1, t1));
    }

    tt3 = t;

    // get the nearest projection and save the normal
    // as the convex combination of the involved two-vertex normals
    if (d[0] <= d[1] && d[0] <= d[2]) {
      ASSIGN(proj, res1);
      (*dist) = d[0];
      ADD_MUL(normal, A, B, tt1);
      double norm = sqrt(DOT(normal, normal));
      normal[0] /= norm;
      normal[1] /= norm;
      normal[2] /= norm;

    } else if (d[1] <= d[0] && d[1] <= d[2]) {
      ASSIGN(proj, res2);
      (*dist) = d[1];
      ADD_MUL(normal, A, C, tt2);
      double norm = sqrt(DOT(normal, normal));
      normal[0] /= norm;
      normal[1] /= norm;
      normal[2] /= norm;
    } else {
      ASSIGN(proj, res3);
      (*dist) = d[2];
      ADD_MUL(normal, B, C, tt3);
      double norm = sqrt(DOT(normal, normal));
      normal[0] /= norm;
      normal[1] /= norm;
      normal[2] /= norm;
    }
  }
  // assign the normal as the plane normal
  // assumes that the plane normal is normalized
  else
    ASSIGN(normal, planes[planeID]);

  return flag;
}

bool MeshSurface::inTriangle(double P[3], double A[3], double B[3],
                             double C[3]) {
  // Compute vectors
  double v0[3], v1[3], v2[3], dot00, dot01, dot02, dot11, dot12, invDenom, u, v;
  SUB(v0, C, A);
  SUB(v1, B, A);
  SUB(v2, P, A);

  // Compute dot products
  dot00 = DOT(v0, v0);
  dot01 = DOT(v0, v1);
  dot02 = DOT(v0, v2);
  dot11 = DOT(v1, v1);
  dot12 = DOT(v1, v2);

  // Compute barycentric coordinates
  invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  // Check if point is in triangle
  if ((u > 0) && (v > 0) && (u + v < 1))
    return true;
  else
    return false;
}

void MeshSurface::point2plane(double p[3], double w[4], double *dist,
                              double proj[3]) {
  double den = (DOT(w, w));
  double d = sqrt(den);
  double val = DOT(w, p) + w[3];
  proj[0] = p[0] - w[0] * (val / (den));
  proj[1] = p[1] - w[1] * (val / (den));
  proj[2] = p[2] - w[2] * (val / (den));
  (*dist) = fabs(val / d);
}

bool MeshSurface::getProjection(double p[3], double *proj1, double *proj2,
                                double *proj3, double *normal1, double *normal2,
                                double *normal3) {
  // get the triangles that are associated to this grid point
  // by querying the auxiliary grid
  double dist, *pp[3];
  set<int> triangles;
  double hside = delphi->hside;

  double aposx = p[0];
  double aposy = p[1];
  double aposz = p[2];

  // move from delphi grid to auxiliary grid
  int irefx = (int)rintp((aposx - xmin) / side);
  int irefy = (int)rintp((aposy - ymin) / side);
  int irefz = (int)rintp((aposz - zmin) / side);

  int ixaux = irefx;
  int iyaux = irefy;
  int izaux = irefz;

  double epsax = aposx + hside, epsay = aposy + hside, epsaz = aposz + hside;
  double epsaxm = aposx - hside, epsaym = aposy - hside, epsazm = aposz - hside;

  for (int i = 0; i < ind[ixaux][iyaux][izaux]; i++)
    triangles.insert((GRIDTRIANGLEMAP(ixaux, iyaux, izaux, i, nx, ny, nz)));

  // keep the nearest triangle
  set<int>::iterator it;
  double locProj[3], minDist = INFINITY, locNorm[3];
  bool ff = false;

  //	printf("\n\n Dist \n");
  for (it = triangles.begin(); it != triangles.end(); it++) {
    pp[0] = vertMatrix[faceMatrix[*it][0]];
    pp[1] = vertMatrix[faceMatrix[*it][1]];
    pp[2] = vertMatrix[faceMatrix[*it][2]];
    bool flag = point2triangle(p, pp[0], pp[1], pp[2], planes[*it], locProj,
                               &dist, locNorm, *it, faceMatrix[*it][0],
                               faceMatrix[*it][1], faceMatrix[*it][2]);

    //		printf("%lf %lf\n",dist,minDist);

    if (dist < minDist) {
      minDist = dist;
      (*proj1) = locProj[0];
      (*proj2) = locProj[1];
      (*proj3) = locProj[2];

      (*normal1) = locNorm[0];
      (*normal2) = locNorm[1];
      (*normal3) = locNorm[2];

      ff = flag;
    }
  }
  if (triangles.size() == 0) {
    {
#ifdef ENABLE_BOOST_THREADS
      boost::mutex::scoped_lock scopedLock(mutex);
#endif
      logging::log<logging::level::warn>("Approximating bgp with grid point");
    }

    (*proj1) = p[0];
    (*proj2) = p[1];
    (*proj3) = p[2];
    return true;
  }

  // further check that the projection is in the cube
  //	if (!testInCube((*proj1),(*proj2),(*proj3),p[0],p[1],p[2],delphi->side))
  //		logging::log<logging::level::warn>("Out of cube projection in
  // MeshSurface::getProjection!");
  return ff;
}
