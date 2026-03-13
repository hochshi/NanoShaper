
#include "MeshSurface.h"


void MeshSurface::clear()
{
	if (faceMatrix != NULL)
		deleteMatrix2D<int>(numTriangles,3,faceMatrix);

	if (vertMatrix != NULL)
		deleteMatrix2D<VERTEX_TYPE>(numVertexes,3,vertMatrix);

	if (vertexTrianglesList != NULL)
	{
		for (int i=0; i<numVertexes; i++)
			delete vertexTrianglesList[i];
		deleteVector<vector<int>*>(vertexTrianglesList);
	}

	if (vertNormals != NULL)
		deleteMatrix2D<VERTEX_TYPE>(numVertexes,3,vertNormals);

	// remove 2d grid for ray casting
	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridTriangleMap != NULL)
		deleteVector<int>(gridTriangleMap);

	if (ind != NULL)
		deleteMatrix3D<int>(nx,ny,nz,ind);
	
	if (gridTriangleMap2D != NULL)
		deleteVector<int>(gridTriangleMap2D);

	if (ind_2d != NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,last_cols_ind,ind_2d);
	#else */
	if (gridTriangleMap != NULL)
	{
		for (int64_t i=0; i < nx*ny*nz; i++)
		{
			gridTriangleMap[i].clear();
		}
		delete[] gridTriangleMap;
		gridTriangleMap = NULL;
	}
	if (gridTriangleMap2D != NULL)
	{
		for (int64_t i=0; i < last_rows_ind*last_cols_ind; i++)
		{
			gridTriangleMap2D[i].clear();
		}
		delete[] gridTriangleMap2D;
		gridTriangleMap2D = NULL;
	}
	// #endif


	if (x != NULL)
		deleteVector<double>(x);
	if (y != NULL)
		deleteVector<double>(y);
	if (z != NULL)
		deleteVector<double>(z);

	/* #if !defined(OPTIMIZE_GRIDS)
	if (ind != NULL)
		deleteMatrix3D<int>(nx,ny,nz,ind);
	#endif */

	if (planes != NULL)
		deleteMatrix2D<double>(numTriangles,4,planes);


	if (patchBasedAlgorithm && num_pixel_intersections != NULL)
	{
		// remove buffers employed in the patch-based ray tracing

		int64_t NX = delphi->nx;
		int64_t NY = delphi->ny;
		int64_t NZ = delphi->nz;
		int64_t N_MAX = MAX(NX, MAX(NY, NZ));

		int panels[] = {3,3};
		if (!delphi->buildStatus && !delphi->buildEpsMap)
			panels[0] = 0;
		else if (delphi->buildStatus && !delphi->buildEpsMap)
			panels[0] = 1;
		if (!accurateTriangulation || isAvailableScalarField)
			panels[1] = 0;
		int numPanels = panels[0] + panels[1];

		#if !defined(SINGLE_PASS_RT)
		#ifdef MULTITHREADING
		deleteVector<int>(num_patch_intersections);

		deleteVector<int>(first_patch_intersection_index);

		deleteVector<pair<VERTEX_TYPE,VERTEX_TYPE*>>(temp_intersections_buffer);

		deleteVector<int64_t>(intersection_pixel_id);
		#endif // MULTITHREADING
		#else // SINGLE_PASS_RT
		for (int i=0; i < numTriangles*numPanels; i++)
		{
			temp_intersections_buffer[i].clear();
		}
		delete[] temp_intersections_buffer;
		temp_intersections_buffer = NULL;

		for (int i=0; i < numTriangles*numPanels; i++)
		{
			intersection_pixel_id[i].clear();
		}
		delete[] intersection_pixel_id;
		intersection_pixel_id = NULL;
		#endif // SINGLE_PASS_RT

		for (int64_t id=0; id < N_MAX*N_MAX* numPanels; id++)
		{
			if (pixel_intersections[id] != NULL)
				deleteVector<pair<VERTEX_TYPE,VERTEX_TYPE*>>(pixel_intersections[id]);
		}
		deleteVector<pair<VERTEX_TYPE,VERTEX_TYPE*>*>(pixel_intersections);

		deleteVector<int>(num_pixel_intersections);
	}
}


MeshSurface::~MeshSurface()
{
	clear();
}


void MeshSurface::init()
{
	numTriangles = 0;
	faceMatrix = NULL;
	vertMatrix = NULL;
	vertNormals = NULL;
	gridTriangleMap = NULL;
	gridTriangleMap2D = NULL;
	vertexTrianglesList = NULL;
	planes = NULL;
	providesAnalyticalNormals = false;
	/* #if !defined(OPTIMIZE_GRIDS)
	ind = NULL;
	ind_2d = NULL;
	#endif */
	x = NULL;
	y = x;
	z = y;
	MAX_TRIANGLES = 250;
	AUX_GRID_DIM = 100;
	AUX_GRID_DIM_2D = 100;
	MAX_TRIANGLES_2D = (MAX_TRIANGLES*AUX_GRID_DIM_2D);
	surfType = GENERIC_SURFACE;

	if (patchBasedAlgorithm)
	{
		num_pixel_intersections = NULL;
		pixel_intersections = NULL;
	}
}


MeshSurface::MeshSurface(DelPhiShared *ds):Surface()
{
	init();
	// set environment
	delphi = ds;
}


void MeshSurface::init(ConfigFile *cf)
{
	string fname = cf->read<string>( "Surface_File_Name", "mesh.off" );
	unsigned int maxMeshDim = cf->read<unsigned int>( "Max_mesh_auxiliary_grid_size", 100 );
	unsigned int maxMeshPatches = cf->read<unsigned int>( "Max_mesh_patches_per_auxiliary_grid_cell", 250 );
	unsigned int maxMeshDim2D = cf->read<unsigned int>( "Max_mesh_auxiliary_grid_2d_size", 100 );
	unsigned int maxMeshPatches2D = cf->read<unsigned int>( "Max_mesh_patches_per_auxiliary_grid_2d_cell", 250 );
	int numMSMSFiles = cf->read<int>( "Num_MSMS_files", 1 );

	setAuxGrid(maxMeshDim,maxMeshPatches);
	setAuxGrid2D(maxMeshDim2D,maxMeshPatches2D);
	// Set up inside value
	inside = 5;

	// get extension
	string ext("");
	int ind = (int)fname.size();

	while (fname[ind]!='.')
	{
		ext = fname.substr(ind,fname.size());
		ind--;
	}

	// get the type of surface
	if (!ext.compare("off"))
	{
		// Load Surface
		bool load_f = load((char*)fname.c_str());
		if (!load_f)
		{
			cout << endl << ERR << "Cannot load " << fname;
			exit(-1);
		}
	}
	else if (!ext.compare("vert") || !ext.compare("face"))
	{
		// Load MSMS surface[s]
		bool load_f = loadMSMS((char*)fname.substr(0,fname.size()-5).c_str(),numMSMSFiles);
		if (!load_f)
		{
			cout << endl << ERR << "Cannot load " << fname;
			exit(-1);
		}
		char cc[BUFLEN];
		strcpy(cc,"msms.off");
		save(cc);
	}
	else
	{
		cout << endl << ERR << "Unknown surface type, please check extension/type";
		cout << endl;
		exit(-1);
	}
}


MeshSurface::MeshSurface():Surface()
{
	init();
	delphi = NULL;
}


MeshSurface::MeshSurface(ConfigFile *cf,DelPhiShared *ds):Surface(cf)
{
	init();
	init(cf);
	// set environment
	delphi = ds;
}


int MeshSurface::getNumPatches (void)
{
	return numTriangles;
}


bool MeshSurface::isPatchBasedRayTracingSupported (void)
{
	return true;
}


void MeshSurface::postRayCasting()
{
	// remove 2d grid for ray casting
	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridTriangleMap2D != NULL)
		deleteVector<int>(gridTriangleMap2D);

	if (ind_2d != NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,last_cols_ind,ind_2d);
	#else */
	if (gridTriangleMap2D != NULL)
	{
		for (int i=0; i < last_rows_ind*last_cols_ind; i++)
		{
			gridTriangleMap2D[i].clear();
		}
		delete[] gridTriangleMap2D;
		gridTriangleMap2D = NULL;
	}
	// #endif

	if (patchBasedAlgorithm && num_pixel_intersections != NULL)
	{
		// remove buffers employed in the patch-based ray tracing

		int64_t NX = delphi->nx;
		int64_t NY = delphi->ny;
		int64_t NZ = delphi->nz;
		int64_t N_MAX = MAX(NX, MAX(NY, NZ));

		int panels[] = {3,3};
		if (!delphi->buildStatus && !delphi->buildEpsMap)
			panels[0] = 0;
		else if (delphi->buildStatus && !delphi->buildEpsMap)
			panels[0] = 1;
		if (!accurateTriangulation || isAvailableScalarField)
			panels[1] = 0;
		int numPanels = panels[0] + panels[1];

		#if !defined(SINGLE_PASS_RT)
		#ifdef MULTITHREADING
		deleteVector<int>(num_patch_intersections);

		deleteVector<int>(first_patch_intersection_index);

		deleteVector<pair<VERTEX_TYPE,VERTEX_TYPE*>>(temp_intersections_buffer);

		deleteVector<int64_t>(intersection_pixel_id);
		#endif // MULTITHREADING
		#else // SINGLE_PASS_RT
		for (int i=0; i < numTriangles*numPanels; i++)
		{
			temp_intersections_buffer[i].clear();
		}
		delete[] temp_intersections_buffer;
		temp_intersections_buffer = NULL;

		for (int i=0; i < numTriangles*numPanels; i++)
		{
			intersection_pixel_id[i].clear();
		}
		delete[] intersection_pixel_id;
		intersection_pixel_id = NULL;
		#endif // SINGLE_PASS_RT

		for (int id=0; id < N_MAX*N_MAX* numPanels; id++)
		{
			if (pixel_intersections[id] != NULL)
				deleteVector<pair<VERTEX_TYPE,VERTEX_TYPE*>>(pixel_intersections[id]);
		}
		deleteVector<pair <VERTEX_TYPE,VERTEX_TYPE*>*>(pixel_intersections);

		deleteVector<int>(num_pixel_intersections);
	}
}


void MeshSurface::preProcessPanel()
{
	numThreadDataWrappers = 1;

	if (vertMatrix == NULL || faceMatrix == NULL)
	{
		cout << endl << WARN << "Cannot get surface without a mesh";
		return;
	}
	// auxiliary bounding box tree to store the patches (useful to ray tracing acceleration) (?)
	int64_t igrid = delphi->nx;

	//cout << endl << "A->";

	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption

	// auxiliary grid is consistent with delphi grid in order to speed-up tracing
	int gridMul = 1;
	while (igrid > AUX_GRID_DIM_2D)
	{
		gridMul += 2;
		int64_t digrid = delphi->nx;
		while (1)
		{
			// get nearest odd multiple
			int64_t fixedPoint = (digrid+(gridMul-1))/gridMul;
			igrid = fixedPoint*gridMul;
			if (igrid%2 == 0)
				digrid = igrid+1;
			else
			{
				igrid = fixedPoint;
				break;
			}
		}
	}
	// auxiliary scale
	scale_2d = delphi->scale/((double)gridMul);
	side_2d = ((double)gridMul)/delphi->scale;

	//cout << endl << INFO << "Auxiliary grid is " << igrid;

	xmin_2d = delphi->baricenter[0] - (igrid-1)*0.5*side_2d;
	ymin_2d = delphi->baricenter[1] - (igrid-1)*0.5*side_2d;
	zmin_2d = delphi->baricenter[2] - (igrid-1)*0.5*side_2d;

	xmax_2d = delphi->baricenter[0] + (igrid-1)*0.5*side_2d;
	ymax_2d = delphi->baricenter[1] + (igrid-1)*0.5*side_2d;
	zmax_2d = delphi->baricenter[2] + (igrid-1)*0.5*side_2d;

	nx_2d = igrid;
	ny_2d = igrid;
	nz_2d = igrid;

	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridTriangleMap2D != NULL)
		deleteVector<int>(gridTriangleMap2D);

	if (ind_2d != NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,last_cols_ind,ind_2d);
	#else */
	if (gridTriangleMap2D != NULL)
	{
		for (int64_t i=0; i < last_rows_ind*last_cols_ind; i++)
		{
			gridTriangleMap2D[i].clear();
		}
		delete[] gridTriangleMap2D;
		gridTriangleMap2D = NULL;
	}
	// #endif

	if (panel == 0)
	{
		last_rows_ind = nz_2d;
		last_cols_ind = ny_2d;
	}
	else if (panel == 1)
	{
		last_rows_ind = ny_2d;
		last_cols_ind = nx_2d;
	}
	else
	{
		last_rows_ind = nz_2d;
		last_cols_ind = nx_2d;
	}

	/* #if !defined(OPTIMIZE_GRIDS)
	gridTriangleMap2D = allocateVector<int>(last_rows_ind*last_cols_ind*MAX_TRIANGLES_2D);

	ind_2d = allocateMatrix2D<unsigned int>(last_rows_ind,last_cols_ind);
	for (int i=0; i<last_rows_ind; i++)
		for (int j=0; j<last_cols_ind; j++)
			ind_2d[i][j] = 0;

	int max_t = 0;
	#else */
	gridTriangleMap2D = new vector<int> [ last_rows_ind*last_cols_ind ];

	for (int64_t i=0; i<last_rows_ind*last_cols_ind; i++)
		gridTriangleMap2D[i].clear();
	// #endif

	VERTEX_TYPE *p[3];

	// build a bounding box for each triangle and map it to
	// the auxiliary grid
	for (int it=0; it<numTriangles; it++)
	{
		// triangle points
		p[0] = vertMatrix[faceMatrix[it][0]];
		p[1] = vertMatrix[faceMatrix[it][1]];
		p[2] = vertMatrix[faceMatrix[it][2]];

		// compute the bounding box of the object
		double downx = INFINITY;
		double downy = INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = -INFINITY;
		double upz = -INFINITY;

		for (unsigned int pind=0;pind<3;pind++)
		{
			downx = MIN(downx,p[pind][0]);
			downy = MIN(downy,p[pind][1]);
			downz = MIN(downz,p[pind][2]);

			upx = MAX(upx,p[pind][0]);
			upy = MAX(upy,p[pind][1]);
			upz = MAX(upz,p[pind][2]);
		}

		// Determine which are the grid cells that are
		// spanned by the bounding box of the triangle
		int64_t ix_start = (int64_t)rintp(fmax(0., downx-xmin_2d)*scale_2d);
		int64_t iy_start = (int64_t)rintp(fmax(0., downy-ymin_2d)*scale_2d);
		int64_t iz_start = (int64_t)rintp(fmax(0., downz-zmin_2d)*scale_2d);

		int64_t ix_end = (int64_t)rintp(fmax(0., upx-xmin_2d)*scale_2d);
		int64_t iy_end = (int64_t)rintp(fmax(0., upy-ymin_2d)*scale_2d);
		int64_t iz_end = (int64_t)rintp(fmax(0., upz-zmin_2d)*scale_2d);

		if (ix_start >= nx_2d)
			ix_start = nx_2d-1;
		if (iy_start >= ny_2d)
			iy_start = ny_2d-1;
		if (iz_start >= nz_2d)
			iz_start = nz_2d-1;

		if (ix_end >= nx_2d)
			ix_end = nx_2d-1;
		if (iy_end >= ny_2d)
			iy_end = ny_2d-1;
		if (iz_end >= nz_2d)
			iz_end = nz_2d-1;

		// map the points into the 2D matrices of each plane

		if (panel == 0)
		{
			// plane YZ
			for (int64_t iz=iz_start; iz<=iz_end; iz++)
			{
				for (int64_t iy=iy_start; iy<=iy_end; iy++)
				{
					/* #if !defined(OPTIMIZE_GRIDS)
					if (ind_2d[iz][iy] > max_t)
						max_t = ind_2d[iz][iy];

					if (ind_2d[iz][iy] >= MAX_TRIANGLES_2D)
					{
						cout << endl << ERR << "Number of triangles is superior to maximum allowed, please increase Max_mesh_patches_per_auxiliary_grid_2d_cell";
						exit(-1);
					}
					GRID_TRIANGLE_MAP_2D(iy,iz,ind_2d[iz][iy],ny_2d,nz_2d) = it;
					ind_2d[iz][iy]++;
					#else */
					gridTriangleMap2D[ iz*ny_2d + iy ].push_back(it);
					// #endif
				}
			}
		}
		else if (panel == 1)
		{
			// plane XY
			for (int64_t iy=iy_start; iy<=iy_end; iy++)
			{
				for (int64_t ix=ix_start; ix<=ix_end; ix++)
				{
					/* #if !defined(OPTIMIZE_GRIDS)
					if (ind_2d[iy][ix] > max_t)
						max_t = ind_2d[iy][ix];

					if (ind_2d[iy][ix] >= MAX_TRIANGLES_2D)
					{
						cout << endl << ERR << "Number of triangles is superior to maximum allowed, please increase Max_mesh_patches_per_auxiliary_grid_2d_cell";
						exit(-1);
					}
					GRID_TRIANGLE_MAP_2D(ix,iy,ind_2d[iy][ix],nx_2d,ny_2d) = it;
					ind_2d[iy][ix]++;
					#else */
					gridTriangleMap2D[ iy*nx_2d + ix ].push_back(it);
					// #endif
				}
			}
		}
		else
		{
			// plane XZ
			for (int64_t iz=iz_start; iz<=iz_end; iz++)
			{
				for (int64_t ix=ix_start; ix<=ix_end; ix++)
				{
					/* #if !defined(OPTIMIZE_GRIDS)
					if (ind_2d[iz][ix] > max_t)
						max_t = ind_2d[iz][ix];

					if (ind_2d[iz][ix] >= MAX_TRIANGLES_2D)
					{
						cout << endl << ERR << "Number of triangles is superior to maximum allowed, please increase Max_mesh_patches_per_auxiliary_grid_2d_cell";
						exit(-1);
					}
					GRID_TRIANGLE_MAP_2D(ix,iz,ind_2d[iz][ix],nx_2d,nz_2d) = it;
					ind_2d[iz][ix]++;
					#else */
					gridTriangleMap2D[ iz*nx_2d + ix ].push_back(it);
					// #endif
				}
			}
		}
	}
}


void MeshSurface::preProcessTriangles()
{
	#if !defined(FLOAT_VERTICES)
	if (planes != NULL)
		deleteMatrix2D<double>(numTriangles,4,planes);

	planes = allocateMatrix2D<double>(numTriangles,4);

	// compute planes from triangles. Planes are normalized in a further step
	// because if normals estimation is needed it is usually more numerically stable to
	// average the plane normals before normalization

	double *p[3];
	for (int it=0; it<numTriangles; it++)
	{
		p[0] = vertMatrix[faceMatrix[it][0]];
		p[1] = vertMatrix[faceMatrix[it][1]];
		p[2] = vertMatrix[faceMatrix[it][2]];

		// compute plane in planes[it]
		plane3points(p[0],p[1],p[2],planes[it]);
	}
	#else
	if (planes != NULL)
		deleteMatrix2D<double>(numTriangles,4,planes);

	planes = allocateMatrix2D<double>(numTriangles,4);

	// compute planes from triangles. Planes are normalized in a further step
	// because if normals estimation is needed it is usually more numerically stable to
	// average the plane normals before normalization

	double p[3][3];
	for (int it=0; it<numTriangles; it++)
	{
		p[0][0] = vertMatrix[faceMatrix[it][0]][0];
		p[0][1] = vertMatrix[faceMatrix[it][0]][1];
		p[0][2] = vertMatrix[faceMatrix[it][0]][2];

		p[1][0] = vertMatrix[faceMatrix[it][1]][0];
		p[1][1] = vertMatrix[faceMatrix[it][1]][1];
		p[1][2] = vertMatrix[faceMatrix[it][1]][2];

		p[2][0] = vertMatrix[faceMatrix[it][2]][0];
		p[2][1] = vertMatrix[faceMatrix[it][2]][1];
		p[2][2] = vertMatrix[faceMatrix[it][2]][2];

		// compute plane in planes[it]
		plane3points(p[0],p[1],p[2],planes[it]);
	}
	#endif

	// compute normals for simple triangular meshes. Estimate vertex normal by averaging
	// the normals of connected triangles.
	// Does not compute normals if vertex normals are provided, such as in MSMS files.
	if (computeNormals)
	{
		if (vertNormals != NULL)
			deleteMatrix2D<VERTEX_TYPE>(numVertexes,3,vertNormals);

		vertNormals = allocateMatrix2D<VERTEX_TYPE>(numVertexes,3);

		for (int iv=0; iv<numVertexes; iv++)
		{
			VERTEX_TYPE meanNormal[3];

			meanNormal[0] = 0;
			meanNormal[1] = 0;
			meanNormal[2] = 0;
			vector<int>::iterator it;

			for (it = vertexTrianglesList[iv]->begin(); it!=vertexTrianglesList[iv]->end(); it++)
			{
				int triangleID = *it;
				ADD(meanNormal,meanNormal,planes[triangleID]);
			}
			// VERTEX_TYPE vlen = (VERTEX_TYPE)vertexTrianglesList[iv]->size();
			// vertNormals[iv][0] = meanNormal[0]/vlen;
			// vertNormals[iv][1] = meanNormal[1]/vlen;
			// vertNormals[iv][2] = meanNormal[2]/vlen;
			// following normalization does all the needed divisions
			vertNormals[iv][0] = meanNormal[0];
			vertNormals[iv][1] = meanNormal[1];
			vertNormals[iv][2] = meanNormal[2];

			// normalization
			VERTEX_TYPE norm = sqrt(DOT(vertNormals[iv],vertNormals[iv]));
			vertNormals[iv][0] /= norm;
			vertNormals[iv][1] /= norm;
			vertNormals[iv][2] /= norm;
		}
	}

	for (int it=0; it<numTriangles; it++)
	{
		#if !defined(FLOAT_VERTICES)
		p[0] = vertMatrix[faceMatrix[it][0]];
		p[1] = vertMatrix[faceMatrix[it][1]];
		p[2] = vertMatrix[faceMatrix[it][2]];

		double norm = sqrt(DOT(planes[it],planes[it]));
		planes[it][0] /= norm;
		planes[it][1] /= norm;
		planes[it][2] /= norm;
		planes[it][3] = -(p[0][0]*(p[1][1]*p[2][2] - p[2][1]*p[1][2])
						+ p[1][0]*(p[2][1]*p[0][2] - p[0][1]*p[2][2])
						+ p[2][0]*(p[0][1]*p[1][2] - p[1][1]*p[0][2]));
		planes[it][3] /= norm;
		#else
		p[0][0] = vertMatrix[faceMatrix[it][0]][0];
		p[0][1] = vertMatrix[faceMatrix[it][0]][1];
		p[0][2] = vertMatrix[faceMatrix[it][0]][2];

		p[1][0] = vertMatrix[faceMatrix[it][1]][0];
		p[1][1] = vertMatrix[faceMatrix[it][1]][1];
		p[1][2] = vertMatrix[faceMatrix[it][1]][2];

		p[2][0] = vertMatrix[faceMatrix[it][2]][0];
		p[2][1] = vertMatrix[faceMatrix[it][2]][1];
		p[2][2] = vertMatrix[faceMatrix[it][2]][2];

		double norm = sqrt(DOT(planes[it],planes[it]));
		planes[it][0] /= norm;
		planes[it][1] /= norm;
		planes[it][2] /= norm;
		planes[it][3] =  -(p[0][0]*(p[1][1]*p[2][2] - p[2][1]*p[1][2])
						+ p[1][0]*(p[2][1]*p[0][2] - p[0][1]*p[2][2])
						+ p[2][0]*(p[0][1]*p[1][2] - p[1][1]*p[0][2]));
		planes[it][3] /= norm;
		#endif
	}

	double minmax[] = {-100, -62, 30,
		65, 65, 195};

	double mid_x = (minmax[3]+minmax[0])*0.5;
	double mid_y = (minmax[4]+minmax[1])*0.5;
	double mid_z = (minmax[5]+minmax[2])*0.5;

	double dummy_atom[6][3] = {
		{minmax[0],mid_y,mid_z},
		{minmax[3],mid_y,mid_z},
		{mid_x,minmax[1],mid_z},
		{mid_x,minmax[4],mid_z},
		{mid_x,mid_y,minmax[2]},
		{mid_x,mid_y,minmax[5]}
	};

	for (int n; n<6; n++)
		printf (" %f %f %f\n", dummy_atom[n][0], dummy_atom[n][1], dummy_atom[n][2]);
}


bool MeshSurface::buildAuxiliaryGrid()
{
	if (faceMatrix == NULL || vertMatrix == NULL)
	{
		cout << endl << WARN << "Cannot get surface without a loaded mesh!";
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
	int64_t igrid = delphi->nx;
	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption

	int gridMul = 1;
	while (igrid > AUX_GRID_DIM)
	{
		gridMul += 2;
		int64_t digrid = delphi->nx;
		while (1)
		{
			// get nearest odd multiple
			int64_t fixedPoint = (digrid+(gridMul-1))/gridMul;
			igrid = fixedPoint*gridMul;
			
			if (igrid%2 == 0)
			{
				digrid = igrid+1;
			}
			else
			{
				igrid = fixedPoint;
				break;
			}
		}
	}
	// auxiliary scale
	scale = delphi->scale/((double)gridMul);
	side = ((double)gridMul)/delphi->scale;

	cout << endl << INFO << "Auxiliary grid is " << igrid << " auxiliary scale is " << scale << " grid mul " << gridMul;

	xmin = delphi->baricenter[0] - (igrid-1)*0.5*side;
	ymin = delphi->baricenter[1] - (igrid-1)*0.5*side;
	zmin = delphi->baricenter[2] - (igrid-1)*0.5*side;

	xmax = delphi->baricenter[0] + (igrid-1)*0.5*side;
	ymax = delphi->baricenter[1] + (igrid-1)*0.5*side;
	zmax = delphi->baricenter[2] + (igrid-1)*0.5*side;

	////////////// Allocate memory for the grid, planes and maps //////////////////////

	if (x != NULL)
		deleteVector<double>(x);
	if (y != NULL)
		deleteVector<double>(y);
	if (z != NULL)
		deleteVector<double>(z);

	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);

	/* #if !defined(OPTIMIZE_GRIDS)
	if (ind != NULL)
		deleteMatrix3D<int>(nx,ny,nz,ind);
	*/

	nx = igrid;
	ny = igrid;
	nz = igrid;
	
	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridTriangleMap != NULL)
		deleteVector<int>(gridTriangleMap);

	cout << endl << INFO << "Allocating " << (nx*ny*nz*MAX_TRIANGLES)*sizeof(int)/1024.0/1024.0 << " MB" << " for the auxiliary grid...";
	
	gridTriangleMap = allocateVector<int>(nx*ny*nz*MAX_TRIANGLES);

	ind = allocateMatrix3D<int>(nx,ny,nz);
	for (int k=0; k<nz; k++)
		for (int j=0; j<ny; j++)
			for (int i=0; i<nx; i++)
				ind[k][j][i] = 0;
	#else */
	if (gridTriangleMap != NULL)
	{
		for (int64_t i=0; i < nx*ny*nz; i++)
		{
			gridTriangleMap[i].clear();
		}
		delete[] gridTriangleMap;
		gridTriangleMap = NULL;
	}
	gridTriangleMap = new vector<int> [ nx*ny*nz ];

	for (int64_t i=0; i<nx*ny*nz; i++)
		gridTriangleMap[i].clear();
	// #endif

	for (int i=0; i<nx; i++)
		x[i] = xmin + i*side;

	for (int i=0;i<ny;i++)
		y[i] = ymin + i*side;

	for (int i=0;i<nz;i++)
		z[i] = zmin + i*side;

	//////////////////////////////////////////////////////////////////////////

	#if !defined(FLOAT_VERTICES)
	double *p[3];
	#else
	double p[3][3];
	#endif
	int max_t=0;

	// build a bounding box for each triangle and map it to
	// the auxiliary grid
	for (int it=0; it<numTriangles; it++)
	{
		#if !defined(FLOAT_VERTICES)
		p[0] = vertMatrix[faceMatrix[it][0]];
		p[1] = vertMatrix[faceMatrix[it][1]];
		p[2] = vertMatrix[faceMatrix[it][2]];
		#else
		p[0][0] = vertMatrix[faceMatrix[it][0]][0];
		p[0][1] = vertMatrix[faceMatrix[it][0]][1];
		p[0][2] = vertMatrix[faceMatrix[it][0]][2];

		p[1][0] = vertMatrix[faceMatrix[it][1]][0];
		p[1][1] = vertMatrix[faceMatrix[it][1]][1];
		p[1][2] = vertMatrix[faceMatrix[it][1]][2];

		p[2][0] = vertMatrix[faceMatrix[it][2]][0];
		p[2][1] = vertMatrix[faceMatrix[it][2]][1];
		p[2][2] = vertMatrix[faceMatrix[it][2]][2];
		#endif
		// compute the bounding box of the object
		double downx = 1e20;
		double downy = 1e20;
		double downz = 1e20;

		double upx = -1e20;
		double upy = -1e20;
		double upz = -1e20;

		for (unsigned int pind=0; pind<3; pind++)
		{
			downx = MIN(downx,p[pind][0]);
			downy = MIN(downy,p[pind][1]);
			downz = MIN(downz,p[pind][2]);

			upx = MAX(upx,p[pind][0]);
			upy = MAX(upy,p[pind][1]);
			upz = MAX(upz,p[pind][2]);
		}

		// Determine which are the grid cubes that
		// are occupied by the object's bounding box
		int64_t ix_start = (int64_t)rintp((downx-xmin)*scale);
		int64_t iy_start = (int64_t)rintp((downy-ymin)*scale);
		int64_t iz_start = (int64_t)rintp((downz-zmin)*scale);

		int64_t ix_end = (int64_t)rintp((upx-xmin)*scale);
		int64_t iy_end = (int64_t)rintp((upy-ymin)*scale);
		int64_t iz_end = (int64_t)rintp((upz-zmin)*scale);

		if (ix_start<0 || iy_start<0 || iz_start<0 || ix_end>=nx || iy_end>=ny || iz_end>=nz)
		{
			cout << endl << ERR << "Triangle bounding box out of grid. Increase perfil";
			exit(-1);
		}

		for (int64_t iz=iz_start; iz<=iz_end; iz++)
			for (int64_t iy=iy_start; iy<=iy_end; iy++)
				for (int64_t ix=ix_start; ix<=ix_end; ix++)
				{
					/* #if !defined(OPTIMIZE_GRIDS)
					if (ind[iz][iy][ix] > max_t)
						max_t = ind[iz][iy][ix];

					if (ind[iz][iy][ix] >= MAX_TRIANGLES)
					{
						cout << endl << ERR << "Number of triangles is superior to maximum allowed, please increase Max_mesh_patches_per_auxiliary_grid_cell";
						exit(-1);
					}
					// gridTriangleMap[iz][iy][ix][ind[iz][iy][ix]]=it;
					GRID_TRIANGLE_MAP(ix,iy,iz,ind[iz][iy][ix],nx,ny,nz) = it;
					ind[iz][iy][ix]++;
					#else */
					gridTriangleMap[ iz*ny*nx + iy*nx + ix ].push_back(it);
					// #endif
				}
	}
	cout << endl << INFO << "Max triangles per cell -> " << max_t;

	return true;
}


#if !defined(SINGLE_PASS_RT)

#if !defined(MINIMIZE_MEMORY)

void MeshSurface::getPatchPreIntersectionData (int64_t nxyz[3], int panels[2],
											   int thread_id, int potentialIntersections[])
{
	if (vertMatrix == NULL || faceMatrix == NULL)
	{
		return;
	}

	int64_t NX = nxyz[0];
	int64_t NY = nxyz[1];
	int64_t NZ = nxyz[2];
	int64_t N_MAX = MAX(NX, MAX(NY, NZ));

	int numPanels = panels[0] + panels[1];

	double *p[3];


	potentialIntersections[thread_id] = 0;

	// Determine the number of potential intersections per object for allocation purposes
	for (int it = thread_id; it < numTriangles; it += conf.numThreads)
	{
		// triangle points
		p[0] = vertMatrix[faceMatrix[it][0]];
		p[1] = vertMatrix[faceMatrix[it][1]];
		p[2] = vertMatrix[faceMatrix[it][2]];

		// compute the bounding box of the object
		double downx = INFINITY;
		double downy = INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = -INFINITY;
		double upz = -INFINITY;

		for (unsigned int pind=0;pind<3;pind++)
		{
			downx = MIN(downx,p[pind][0]);
			downy = MIN(downy,p[pind][1]);
			downz = MIN(downz,p[pind][2]);

			upx = MAX(upx,p[pind][0]);
			upy = MAX(upy,p[pind][1]);
			upz = MAX(upz,p[pind][2]);
		}

		#ifdef CELL_CULLING
		if (downx > delphi->xmax + EPS_INT)
			continue;
		if (downy > delphi->ymax + EPS_INT)
			continue;
		if (downz > delphi->zmax + EPS_INT)
			continue;
		if (upx < delphi->xmin - (delphi->hside - delta_accurate_triangulation + EPS_INT))
			continue;
		if (upy < delphi->ymin - (delphi->hside - delta_accurate_triangulation + EPS_INT))
			continue;
		if (upz < delphi->zmin - (delphi->hside - delta_accurate_triangulation + EPS_INT))
			continue;
		#endif

		int64_t i_start[3], i_end[3];

		// Determine the object's bounding rectangles
		i_start[0] = (int64_t)rintp(fmax(0., downx-delphi->xmin)*delphi->scale);
		i_start[1] = (int64_t)rintp(fmax(0., downy-delphi->ymin)*delphi->scale);
		i_start[2] = (int64_t)rintp(fmax(0., downz-delphi->zmin)*delphi->scale);

		i_end[0] = (int64_t)rintp(fmax(0., upx-delphi->xmin)*delphi->scale);
		i_end[1] = (int64_t)rintp(fmax(0., upy-delphi->ymin)*delphi->scale);
		i_end[2] = (int64_t)rintp(fmax(0., upz-delphi->zmin)*delphi->scale);

		#if !defined(CELL_CULLING) // if culling was carried out, the following will not be done
		if (i_start[0] >= NX) i_start[0] = NX-1;
		if (i_start[1] >= NY) i_start[1] = NY-1;
		if (i_start[2] >= NZ) i_start[2] = NZ-1;
		#endif
		if (i_end[0] >= NX) i_end[0] = NX-1;
		if (i_end[1] >= NY) i_end[1] = NY-1;
		if (i_end[2] >= NZ) i_end[2] = NZ-1;

		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			if (panels[ int_phase ] == 0) continue;

			double delta = (int_phase == 0) ? 0. : delta_accurate_triangulation - delphi->hside;

			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				int first_dim = 1, last__dim = 2;
				if (panel == 1)
					first_dim = 0, last__dim = 1;
				else if (panel == 2)
					first_dim = 0, last__dim = 2;

				#if defined(MULTITHREADING)
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);
				#endif

				int num_pixels = (1+i_end[last__dim]-i_start[last__dim]) * (1+i_end[first_dim]-i_start[first_dim]);
				#if !defined(MULTITHREADING)
				for (int64_t n=i_start[last__dim]; n<=i_end[last__dim]; n++)
				{
					for (int64_t m=i_start[first_dim]; m<=i_end[first_dim]; m++)
					{
						int64_t id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);
						++num_pixel_intersections[id];
					}
				}
				#else // MULTITHREADING
				num_patch_intersections[ patch_id ] += num_pixels;
				#endif // MULTITHREADING

				potentialIntersections[thread_id] += num_pixels;
			}
		}
	}
}

#else // MINIMIZE_MEMORY

void MeshSurface::getPatchPreIntersectionData (int64_t nxyz[3], int panels[2],
											   int thread_id, int potentialIntersections[])
{
	if (vertMatrix == NULL || faceMatrix == NULL)
	{
		return;
	}

	int64_t NX = nxyz[0];
	int64_t NY = nxyz[1];
	int64_t NZ = nxyz[2];
	int64_t N_MAX = MAX(NX, MAX(NY, NZ));

	int numPanels = panels[0] + panels[1];

	double pa[3];
	#if !defined(FLOAT_VERTICES)
	double *p[3];
	#else
	double p[3][3];
	#endif


	potentialIntersections[thread_id] = 0;

	// Determine the number of potential intersections per object for allocation purposes
	for (int it = thread_id; it < numTriangles; it += conf.numThreads)
	{
		// triangle points
		#if !defined(FLOAT_VERTICES)
		p[0] = vertMatrix[faceMatrix[it][0]];
		p[1] = vertMatrix[faceMatrix[it][1]];
		p[2] = vertMatrix[faceMatrix[it][2]];
		#else
		p[0][0] = vertMatrix[faceMatrix[it][0]][0];
		p[0][1] = vertMatrix[faceMatrix[it][0]][1];
		p[0][2] = vertMatrix[faceMatrix[it][0]][2];

		p[1][0] = vertMatrix[faceMatrix[it][1]][0];
		p[1][1] = vertMatrix[faceMatrix[it][1]][1];
		p[1][2] = vertMatrix[faceMatrix[it][1]][2];

		p[2][0] = vertMatrix[faceMatrix[it][2]][0];
		p[2][1] = vertMatrix[faceMatrix[it][2]][1];
		p[2][2] = vertMatrix[faceMatrix[it][2]][2];
		#endif

		// compute the bounding box of the object
		double downx = INFINITY;
		double downy = INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = -INFINITY;
		double upz = -INFINITY;

		for (unsigned int pind=0;pind<3;pind++)
		{
			downx = MIN(downx,p[pind][0]);
			downy = MIN(downy,p[pind][1]);
			downz = MIN(downz,p[pind][2]);

			upx = MAX(upx,p[pind][0]);
			upy = MAX(upy,p[pind][1]);
			upz = MAX(upz,p[pind][2]);
		}

		#ifdef CELL_CULLING
		if (downx > delphi->xmax + EPS_INT)
			continue;
		if (downy > delphi->ymax + EPS_INT)
			continue;
		if (downz > delphi->zmax + EPS_INT)
			continue;
		if (upx < delphi->xmin - (delphi->hside - delta_accurate_triangulation + EPS_INT))
			continue;
		if (upy < delphi->ymin - (delphi->hside - delta_accurate_triangulation + EPS_INT))
			continue;
		if (upz < delphi->zmin - (delphi->hside - delta_accurate_triangulation + EPS_INT))
			continue;
		#endif

		int64_t i_start[3], i_end[3];

		// Determine the object's bounding rectangles
		i_start[0] = (int64_t)rintp(fmax(0., downx-delphi->xmin)*delphi->scale);
		i_start[1] = (int64_t)rintp(fmax(0., downy-delphi->ymin)*delphi->scale);
		i_start[2] = (int64_t)rintp(fmax(0., downz-delphi->zmin)*delphi->scale);

		i_end[0] = (int64_t)rintp(fmax(0., upx-delphi->xmin)*delphi->scale);
		i_end[1] = (int64_t)rintp(fmax(0., upy-delphi->ymin)*delphi->scale);
		i_end[2] = (int64_t)rintp(fmax(0., upz-delphi->zmin)*delphi->scale);

		#if !defined(CELL_CULLING) // if culling was carried out, the following will not be done
		if (i_start[0]>=NX) i_start[0]=NX-1;
		if (i_start[1]>=NY) i_start[1]=NY-1;
		if (i_start[2]>=NZ) i_start[2]=NZ-1;
		#endif
		if (i_end[0]>=NX) i_end[0]=NX-1;
		if (i_end[1]>=NY) i_end[1]=NY-1;
		if (i_end[2]>=NZ) i_end[2]=NZ-1;

		int64_t pixels[3] = {1+i_end[0]-i_start[0], 1+i_end[1]-i_start[1], 1+i_end[2]-i_start[2]};

		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			if (panels[ int_phase ] == 0) continue;

			double delta = (int_phase == 0) ? 0. : delta_accurate_triangulation - delphi->hside;

			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				double dir[3] = {0., 0., 0.};
				double ray_dir;
				double t,u,v;

				int first_dim, last__dim;

				if (panel == 0)
				{
					first_dim = 1, last__dim = 2;
					pa[0] = delphi->x[0] + delta;
					dir[0] = ray_dir = delphi->x[NX-1] - delphi->x[0];
				}
				else if (panel == 1)
				{
					first_dim = 0, last__dim = 1;
					pa[2] = delphi->z[0] + delta;
					dir[2] = ray_dir = delphi->z[NZ-1] - delphi->z[0];
				}
				else
				{
					first_dim = 0, last__dim = 2;
					pa[1] = delphi->y[0] + delta;
					dir[1] = ray_dir = delphi->y[NY-1] - delphi->y[0];
				}
				#if defined(MULTITHREADING)
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);
				#endif

				for (int64_t rectangle_pixel = 0; rectangle_pixel < pixels[first_dim]*pixels[last__dim]; rectangle_pixel++)
				{
					int64_t n = rectangle_pixel / pixels[first_dim];
					int64_t m = rectangle_pixel - n*pixels[first_dim];

					n += i_start[last__dim];
					m += i_start[first_dim];

					if (panel == 0)
					{
						pa[1] = delphi->y[m] + delta;
						if (pa[1] < downy || pa[1] > upy)
							continue;
						pa[2] = delphi->z[n] + delta;
						if (pa[2] < downz || pa[2] > upz)
							continue;
					}
					else if (panel == 1)
					{
						pa[0] = delphi->x[m] + delta;
						if (pa[0] < downx || pa[0] > upx)
							continue;
						pa[1] = delphi->y[n] + delta;
						if (pa[1] < downy || pa[1] > upy)
							continue;
					}
					else
					{
						pa[0] = delphi->x[m] + delta;
						if (pa[0] < downx || pa[0] > upx)
							continue;
						pa[2] = delphi->z[n] + delta;
						if (pa[2] < downz || pa[2] > upz)
							continue;
					}

					// The ray-triangle intersection could be optimised recalling that many components are nil
					if (intersect_triangle(pa,dir,p[0],p[1],p[2],&t,&u,&v))
					{
						#if !defined(MULTITHREADING)
						int64_t panel_pixel_id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);
						++num_pixel_intersections[ panel_pixel_id ];
						#else
						++num_patch_intersections[ patch_id ];
						#endif

						++potentialIntersections[thread_id];
					}
				}
			}
		}
	}
}

#endif // MINIMIZE_MEMORY

#endif // SINGLE_PASS_RT


void MeshSurface::getPatchIntersectionData (int64_t nxyz[3], int panels[2],
											int thread_id, int *netIntersections)
{
	if (vertMatrix == NULL || faceMatrix == NULL)
	{
		return;
	}

	int64_t NX = nxyz[0];
	int64_t NY = nxyz[1];
	int64_t NZ = nxyz[2];
	int64_t N_MAX = MAX(NX, MAX(NY, NZ));

	int numPanels = panels[0] + panels[1];

	double pa[3];
	#if !defined(FLOAT_VERTICES)
	double *p[3];
	#else
	double p[3][3];
	#endif


	*netIntersections = 0;

	// Perform the per-patch ray casting
	for (int it = thread_id; it < numTriangles; it += conf.numThreads)
	{
		// triangle points
		#if !defined(FLOAT_VERTICES)
		p[0] = vertMatrix[faceMatrix[it][0]];
		p[1] = vertMatrix[faceMatrix[it][1]];
		p[2] = vertMatrix[faceMatrix[it][2]];
		#else
		p[0][0] = vertMatrix[faceMatrix[it][0]][0];
		p[0][1] = vertMatrix[faceMatrix[it][0]][1];
		p[0][2] = vertMatrix[faceMatrix[it][0]][2];

		p[1][0] = vertMatrix[faceMatrix[it][1]][0];
		p[1][1] = vertMatrix[faceMatrix[it][1]][1];
		p[1][2] = vertMatrix[faceMatrix[it][1]][2];

		p[2][0] = vertMatrix[faceMatrix[it][2]][0];
		p[2][1] = vertMatrix[faceMatrix[it][2]][1];
		p[2][2] = vertMatrix[faceMatrix[it][2]][2];
		#endif

		// compute the bounding box of the object
		double downx = INFINITY;
		double downy = INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = -INFINITY;
		double upz = -INFINITY;

		for (unsigned int pind=0;pind<3;pind++)
		{
			downx = MIN(downx,p[pind][0]);
			downy = MIN(downy,p[pind][1]);
			downz = MIN(downz,p[pind][2]);

			upx = MAX(upx,p[pind][0]);
			upy = MAX(upy,p[pind][1]);
			upz = MAX(upz,p[pind][2]);
		}

		#ifdef CELL_CULLING
		if (downx > delphi->xmax + EPS_INT)
			continue;
		if (downy > delphi->ymax + EPS_INT)
			continue;
		if (downz > delphi->zmax + EPS_INT)
			continue;
		if (upx < delphi->xmin - (delphi->hside - delta_accurate_triangulation + EPS_INT))
			continue;
		if (upy < delphi->ymin - (delphi->hside - delta_accurate_triangulation + EPS_INT))
			continue;
		if (upz < delphi->zmin - (delphi->hside - delta_accurate_triangulation + EPS_INT))
			continue;
		#endif

		int64_t i_start[3], i_end[3];

		// Determine the object's bounding rectangles
		i_start[0] = (int64_t)rintp(fmax(0., downx-delphi->xmin)*delphi->scale);
		i_start[1] = (int64_t)rintp(fmax(0., downy-delphi->ymin)*delphi->scale);
		i_start[2] = (int64_t)rintp(fmax(0., downz-delphi->zmin)*delphi->scale);

		i_end[0] = (int64_t)rintp(fmax(0., upx-delphi->xmin)*delphi->scale);
		i_end[1] = (int64_t)rintp(fmax(0., upy-delphi->ymin)*delphi->scale);
		i_end[2] = (int64_t)rintp(fmax(0., upz-delphi->zmin)*delphi->scale);

		#if !defined(CELL_CULLING) // if culling was carried out, the following will not be done
		if (i_start[0] >= NX) i_start[0] = NX-1;
		if (i_start[1] >= NY) i_start[1] = NY-1;
		if (i_start[2] >= NZ) i_start[2] = NZ-1;
		#endif
		if (i_end[0] >= NX) i_end[0] = NX-1;
		if (i_end[1] >= NY) i_end[1] = NY-1;
		if (i_end[2] >= NZ) i_end[2] = NZ-1;

		int64_t pixels[3] = {1+i_end[0]-i_start[0], 1+i_end[1]-i_start[1], 1+i_end[2]-i_start[2]};

		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			if (panels[ int_phase ] == 0) continue;

			double delta = (int_phase == 0) ? 0. : delta_accurate_triangulation - delphi->hside;

			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				double dir[3] = {0., 0., 0.};
				double ray_dir;
				double t,u,v;

				int first_dim, last__dim;

				if (panel == 0)
				{
					first_dim = 1, last__dim = 2;
					pa[0] = delphi->x[0] + delta;
					dir[0] = ray_dir = delphi->x[NX-1] - delphi->x[0];
				}
				else if (panel == 1)
				{
					first_dim = 0, last__dim = 1;
					pa[2] = delphi->z[0] + delta;
					dir[2] = ray_dir = delphi->z[NZ-1] - delphi->z[0];
				}
				else
				{
					first_dim = 0, last__dim = 2;
					pa[1] = delphi->y[0] + delta;
					dir[1] = ray_dir = delphi->y[NY-1] - delphi->y[0];
				}
				#if defined(MULTITHREADING)
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);
				#if !defined(SINGLE_PASS_RT)
				int first_id = first_patch_intersection_index[patch_id];
				#endif
				#endif

				for (int64_t rectangle_pixel = 0; rectangle_pixel < pixels[first_dim]*pixels[last__dim]; rectangle_pixel++)
				{
					int64_t n = rectangle_pixel / pixels[first_dim];
					int64_t m = rectangle_pixel - n*pixels[first_dim];

					n += i_start[last__dim];
					m += i_start[first_dim];

					if (panel == 0)
					{
						pa[1] = delphi->y[m] + delta;
						if (pa[1] < downy || pa[1] > upy)
							continue;
						pa[2] = delphi->z[n] + delta;
						if (pa[2] < downz || pa[2] > upz)
							continue;
					}
					else if (panel == 1)
					{
						pa[0] = delphi->x[m] + delta;
						if (pa[0] < downx || pa[0] > upx)
							continue;
						pa[1] = delphi->y[n] + delta;
						if (pa[1] < downy || pa[1] > upy)
							continue;
					}
					else
					{
						pa[0] = delphi->x[m] + delta;
						if (pa[0] < downx || pa[0] > upx)
							continue;
						pa[2] = delphi->z[n] + delta;
						if (pa[2] < downz || pa[2] > upz)
							continue;
					}

					// The ray-triangle intersection could be optimised recalling that many components are nil
					if (intersect_triangle(pa,dir,p[0],p[1],p[2],&t,&u,&v))
					{
						// the n- and m-dependent integers are rarely computed here since it is very likely
						// that the ray misses the patch
						int64_t panel_pixel_id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);

						#if !defined(SINGLE_PASS_RT)

						#if !defined(MULTITHREADING)
						pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p = &pixel_intersections[panel_pixel_id][ num_pixel_intersections[panel_pixel_id]++ ];
						#else
						pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p = &temp_intersections_buffer[ first_id + num_patch_intersections[patch_id] ];
						intersection_pixel_id[ first_id + num_patch_intersections[patch_id]++ ] = panel_pixel_id;
						#endif

						int_p->first = t;
						int_p->second = nullptr;

						++*netIntersections;

						#else // SINGLE_PASS_RT

						pair<VERTEX_TYPE,VERTEX_TYPE*> new_int;

						new_int.first = t;
						new_int.second = nullptr;

						temp_intersections_buffer[ patch_id ].push_back(new_int);

						intersection_pixel_id[ patch_id ].push_back(panel_pixel_id);

						++*netIntersections;

						#endif // SINGLE_PASS_RT
					}
				}
			}
		}
	}
}


// returns 1 in case of intersection
int MeshSurface::intersect_triangle(double orig[3], double dir[3],
									double vert0[3], double vert1[3], double vert2[3], double *t, double *u, double *v)
{
   	double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
	double det, inv_det;

	/* find vectors for two edges sharing vert0 */
	SUB(edge1, vert1, vert0);
	SUB(edge2, vert2, vert0);

	/* begin calculating determinant - also used to calculate U parameter */
	CROSS(pvec, dir, edge2);

	/* if determinant is near zero, ray lies in plane of triangle */
	det = DOT(edge1, pvec);

	#ifdef TEST_CULL             /* define TEST_CULL if culling is desired */
	if (det < EPS)
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

	/* ray intersects triangle; calculate t, u, v */
	*t = DOT(edge2, qvec) * inv_det;
	#else                    /* the non-culling branch */
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

	/* ray intersects triangle; calculate t, u, v */
	*t = DOT(edge2, qvec) * inv_det;
	#endif

	return 1;
}


// returns 1 in case of intersection, optimised case for ray parallel to axis X
int MeshSurface::intersect_triangle_X(double orig[3], double dir_x,
									  double vert0[3], double vert1[3], double vert2[3], double *t)
{
	double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
	double det, inv_det, u, v;

	/* find vectors for two edges sharing vert0 */
	// SUB(edge1, vert1, vert0);
	edge1[1] = vert1[1] - vert0[1];
	edge1[2] = vert1[2] - vert0[2];
	// SUB(edge2, vert2, vert0);
	edge2[1] = vert2[1] - vert0[1];
	edge2[2] = vert2[2] - vert0[2];

	/* begin calculating determinant - also used to calculate U parameter */
	// CROSS(pvec, dir, edge2);
	pvec[1] = -dir_x * edge2[2];
	pvec[2] =  dir_x * edge2[1];

	/* if determinant is near zero, ray lies in plane of triangle */
	// det = DOT(edge1, pvec);
	det = edge1[1]*pvec[1] + edge1[2]*pvec[2];

	if (det > -EPS && det < EPS)
		return 0;

	inv_det = 1.0 / det;

	/* calculate distance from vert0 to ray origin */
	// SUB(tvec, orig, vert0);
	// tvec[0] = orig[0] - vert0[0];
	tvec[1] = orig[1] - vert0[1];
	tvec[2] = orig[2] - vert0[2];
	/* calculate U parameter and test bounds */
	// u = DOT(tvec, pvec) * inv_det;
	u = (tvec[1]*pvec[1] + tvec[2]*pvec[2]) * inv_det;

	if (u < 0.0 || u > 1.0)
		return 0;

	/* prepare to test V parameter */
	// CROSS(qvec, tvec, edge1);
	qvec[0] = tvec[1]*edge1[2] - tvec[2]*edge1[1];
	// qvec[1] and qvec[2] are eventually calculated later ...

	/* calculate V parameter and test bounds */
	// *v = DOT(dir, qvec) * inv_det;
	v = dir_x * qvec[0] * inv_det;

	if (v < 0.0 || u + v > 1.0)
		return 0;

	// tvec[0], edge1[0] and edge2[0] were not calculated ...
	tvec[0] = orig[0] - vert0[0];
	edge1[0] = vert1[0] - vert0[0];
	edge2[0] = vert2[0] - vert0[0];
	// qvec[1] and qvec[2] were not calculated ...
	qvec[1] = tvec[2]*edge1[0] - tvec[0]*edge1[2];
	qvec[2] = tvec[0]*edge1[1] - tvec[1]*edge1[0];
	/* ray intersects triangle */
	*t = DOT(edge2, qvec) * inv_det;

	return 1;
}


// returns 1 in case of intersection, optimised case for ray parallel to axis Y
int MeshSurface::intersect_triangle_Y(double orig[3], double dir_y,
									  double vert0[3], double vert1[3], double vert2[3], double *t)
{
	double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
	double det, inv_det, u, v;

	/* find vectors for two edges sharing vert0 */
	// SUB(edge1, vert1, vert0);
	edge1[0] = vert1[0] - vert0[0];
	edge1[2] = vert1[2] - vert0[2];
	// SUB(edge2, vert2, vert0);
	edge2[0] = vert2[0] - vert0[0];
	edge2[2] = vert2[2] - vert0[2];

	/* begin calculating determinant - also used to calculate U parameter */
	// CROSS(pvec, dir, edge2);
	pvec[0] =  dir_y * edge2[2];
	pvec[2] = -dir_y * edge2[0];

	/* if determinant is near zero, ray lies in plane of triangle */
	// det = DOT(edge1, pvec);
	det = edge1[0]*pvec[0] + edge1[2]*pvec[2];

	if (det > -EPS && det < EPS)
		return 0;

	inv_det = 1.0 / det;

	/* calculate distance from vert0 to ray origin */
	// SUB(tvec, orig, vert0);
	tvec[0] = orig[0] - vert0[0];
	// tvec[1] = orig[1] - vert0[1];
	tvec[2] = orig[2] - vert0[2];
	/* calculate U parameter and test bounds */
	// u = DOT(tvec, pvec) * inv_det;
	u = (tvec[0]*pvec[0] + tvec[2]*pvec[2]) * inv_det;

	if (u < 0.0 || u > 1.0)
		return 0;

	/* prepare to test V parameter */
	// CROSS(qvec, tvec, edge1);
	qvec[1] = tvec[2]*edge1[0] - tvec[0]*edge1[2];
	// qvec[0] and qvec[2] are eventually calculated later ...

	/* calculate V parameter and test bounds */
	// *v = DOT(dir, qvec) * inv_det;
	v = dir_y * qvec[1] * inv_det;

	if (v < 0.0 || u + v > 1.0)
		return 0;

	// tvec[1], edge1[1] and edge2[1] were not calculated ...
	tvec[1] = orig[1] - vert0[1];
	edge1[1] = vert1[1] - vert0[1];
	edge2[1] = vert2[1] - vert0[1];
	// qvec[0] and qvec[2] were not calculated ...
	qvec[0] = tvec[1]*edge1[2] - tvec[2]*edge1[1];
	qvec[2] = tvec[0]*edge1[1] - tvec[1]*edge1[0];
	/* ray intersects triangle */
	*t = DOT(edge2, qvec) * inv_det;

	return 1;
}


// returns 1 in case of intersection, optimised case for ray parallel to axis Z
int MeshSurface::intersect_triangle_Z(double orig[3], double dir_z,
									  double vert0[3], double vert1[3], double vert2[3], double *t)
{
	double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
	double det, inv_det, u, v;

	/* find vectors for two edges sharing vert0 */
	// SUB(edge1, vert1, vert0);
	edge1[0] = vert1[0] - vert0[0];
	edge1[1] = vert1[1] - vert0[1];
	// SUB(edge2, vert2, vert0);
	edge2[0] = vert2[0] - vert0[0];
	edge2[1] = vert2[1] - vert0[1];

	/* begin calculating determinant - also used to calculate U parameter */
	// CROSS(pvec, dir, edge2);
	pvec[0] = -dir_z * edge2[1];
	pvec[1] =  dir_z * edge2[0];

	/* if determinant is near zero, ray lies in plane of triangle */
	// det = DOT(edge1, pvec);
	det = edge1[0]*pvec[0] + edge1[1]*pvec[1];

	if (det > -EPS && det < EPS)
		return 0;

	inv_det = 1.0 / det;

	/* calculate distance from vert0 to ray origin */
	// SUB(tvec, orig, vert0);
	tvec[0] = orig[0] - vert0[0];
	tvec[1] = orig[1] - vert0[1];
	// tvec[2] = orig[2] - vert0[2];
	/* calculate U parameter and test bounds */
	// u = DOT(tvec, pvec) * inv_det;
	u = (tvec[0]*pvec[0] + tvec[1]*pvec[1]) * inv_det;

	if (u < 0.0 || u > 1.0)
		return 0;

	/* prepare to test V parameter */
	// CROSS(qvec, tvec, edge1);
	qvec[2] = tvec[0]*edge1[1] - tvec[1]*edge1[0];
	// qvec[0] and qvec[1] are eventually calculated later ...

	/* calculate V parameter and test bounds */
	// *v = DOT(dir, qvec) * inv_det;
	v = dir_z * qvec[2] * inv_det;

	if (v < 0.0 || u + v > 1.0)
		return 0;

	// tvec[2], edge1[2] and edge2[2] were not calculated ...
	tvec[2] = orig[2] - vert0[2];
	edge1[2] = vert1[2] - vert0[2];
	edge2[2] = vert2[2] - vert0[2];
	// qvec[0] and qvec[1] were not calculated ...
	qvec[0] = tvec[1]*edge1[2] - tvec[2]*edge1[1];
	qvec[1] = tvec[2]*edge1[0] - tvec[0]*edge1[2];
	/* ray intersects triangle */
	*t = DOT(edge2, qvec) * inv_det;

	return 1;
}


inline void MeshSurface::getRayIntersection(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id)
{
	int first_dim = 1, last__dim = 2;
	int varying_coord = 0;

	if (panel == 1)
	{
		first_dim = 0;
		last__dim = 1;
		varying_coord = 2;
	}
	else if (panel == 2)
	{
		first_dim = 0;
		last__dim = 2;
		varying_coord = 1;
	}

	double min_2d[3];

	min_2d[0] = xmin_2d;
	min_2d[1] = ymin_2d;
	min_2d[2] = zmin_2d;

	int64_t i1 = (int64_t)rintp((pa[first_dim]-min_2d[first_dim])*scale_2d);
	int64_t i2 = (int64_t)rintp((pa[last__dim]-min_2d[last__dim])*scale_2d);

	int64_t n_2d[3];

	n_2d[0] = nx_2d;
	n_2d[1] = ny_2d;
	n_2d[2] = nz_2d;

	int64_t n_2d_first = n_2d[first_dim];
	// int64_t n_2d_last  = n_2d[last__dim];

	/* #if !defined(OPTIMIZE_GRIDS)
	int numTriangles = ind_2d[i2][i1];
	#else */
	int numTriangles = gridTriangleMap2D[ i2*n_2d_first + i1 ].size();
	// #endif

	if (numTriangles == 0) return;

	double dir[3];
	#if !defined(FLOAT_VERTICES)
	double *p[3];
	#else
	double p[3][3];
	#endif

	dir[0] = 0.;
	dir[1] = 0.;
	dir[2] = 0.;
	dir[varying_coord] = pb[varying_coord] - pa[varying_coord];

	for (int iter = 0; iter<numTriangles; iter++)
	{
		/* #if !defined(OPTIMIZE_GRIDS)
		int it = GRID_TRIANGLE_MAP_2D(i1,i2,iter,n_2d_first,n_2d_last);
		#else */
		int it = gridTriangleMap2D[i2*n_2d_first+i1][iter];
		// #endif
		#if !defined(FLOAT_VERTICES)
		p[0] = vertMatrix[faceMatrix[it][0]];
		p[1] = vertMatrix[faceMatrix[it][1]];
		p[2] = vertMatrix[faceMatrix[it][2]];
		#else
		p[0][0] = vertMatrix[faceMatrix[it][0]][0];
		p[0][1] = vertMatrix[faceMatrix[it][0]][1];
		p[0][2] = vertMatrix[faceMatrix[it][0]][2];
		p[1][0] = vertMatrix[faceMatrix[it][1]][0];
		p[1][1] = vertMatrix[faceMatrix[it][1]][1];
		p[1][2] = vertMatrix[faceMatrix[it][1]][2];
		p[2][0] = vertMatrix[faceMatrix[it][2]][0];
		p[2][1] = vertMatrix[faceMatrix[it][2]][1];
		p[2][2] = vertMatrix[faceMatrix[it][2]][2];
		#endif

		double t,u,v;
		// The ray-triangle intersection could be optimised recalling that many components are nil
		if (intersect_triangle(pa,dir,p[0],p[1],p[2],&t,&u,&v))
		{
			#if !defined(FLOAT_VERTICES)
			intersections.push_back(pair<double,double*>(t,planes[it]));
			#else
			float *p = allocateVector<float>(3);
			p[0] = planes[it][0];
			p[1] = planes[it][1];
			p[2] = planes[it][2];
			intersections.push_back(pair<float,float*>(t,p));
			#endif
		}
	}

	if (intersections.size()>0)
		sort(intersections.begin(), intersections.end(), compKeepIndex);
}


#if defined(REPORT_FAILED_RAYS)
inline void MeshSurface::printRayIntersection(double pa[3], double pb[3])
{
	;
}
#endif


bool MeshSurface::point2triangle(double P[3], double A[3], double B[3], double C[3], double w[4], double *proj, double *dist, double *normal, int planeID)
{
	bool flag;
	// project to plane
	point2plane(P,w,dist,proj);
	// test triangle
	flag = inTriangle(proj,A,B,C);

	// if not in triangle return the projection
	// on the nearest edge to the query point
	// the normal in this case is estimated as the convex combination
	// of the vertex normals of the vertexes of the edge where the projection lies.
	if (!flag)
	{
		double t1[3],t2[3],l2,d[3];
		double res1[3],res2[3],res3[3];
		double tt1,tt2,tt3;

		// A,B edge	projection
		SUB(t2,B,A);
		l2 = DOT(t2,t2);
		SUB(t1,proj,A);
		double t = (DOT(t1,t2));

		if (t < 0.0)
		{
			SUB(t1,P,A);
			d[0] = sqrt(DOT(t1,t1));// Beyond the first end of the segment
			ASSIGN(res1,A);
			// saturation
			t = 0.0;
		}
		else if (t > l2)
		{
			SUB(t1,P,B);
			d[0] = sqrt(DOT(t1,t1));// Beyond the second end of the segment
			ASSIGN(res1,B);
			// saturation
			t = 1.0;
		}
		else
		{
			t /= l2;
			ADD_MUL(res1,A,t2,t);
			SUB(t1,res1,P);
			d[0] = sqrt(DOT(t1,t1));
		}
		tt1 = t;

		// A,C edge	projection
		SUB(t2,C,A);
		l2 = DOT(t2,t2);
		SUB(t1,proj,A);
		t = (DOT(t1,t2));

		if (t < 0.0)
		{
			SUB(t1,P,A);
			d[1] = sqrt(DOT(t1,t1));// Beyond the first end of the segment
			ASSIGN(res2,A);
			t = 0.0;
		}
		else if (t > l2)
		{
			SUB(t1,P,C);
			d[1] = sqrt(DOT(t1,t1));// Beyond the second end of the segment
			ASSIGN(res2,C);
			t = 1.0;
		}
		else
		{
			t /= l2;
			ADD_MUL(res2,A,t2,t);
			SUB(t1,res2,P);
			d[1] = sqrt(DOT(t1,t1));
		}
		tt2 = t;

		// B,C edge	projection
		SUB(t2,C,B);
		l2 = DOT(t2,t2);
		SUB(t1,proj,B);
		t = (DOT(t1,t2));

		if (t < 0.0)
		{
			SUB(t1,P,B);
			d[2] = sqrt(DOT(t1,t1));// Beyond the first end of the segment
			ASSIGN(res3,B);
			t = 0.0;
		}
		else if (t > l2)
		{
			SUB(t1,P,C);
			d[2] = sqrt(DOT(t1,t1));// Beyond the second end of the segment
			ASSIGN(res3,C);
			t = 1.0;
		}
		else
		{
			t /= l2;
			ADD_MUL(res3,B,t2,t);
			SUB(t1,res3,P);
			d[2] = sqrt(DOT(t1,t1));
		}
		tt3 = t;

		// get the nearest projection and save the normal
		// as the convex combination of the involved two-vertex normals
		if (d[0] <= d[1] && d[0] <= d[2])
		{
			ASSIGN(proj,res1);
			*dist = d[0];
			ADD_MUL(normal,A,B,tt1);
		}
		else if (d[1] <= d[0] && d[1] <= d[2])
		{
			ASSIGN(proj,res2);
			*dist = d[1];
			ADD_MUL(normal,A,C,tt2);
		}
		else
		{
			ASSIGN(proj,res3);
			*dist = d[2];
			ADD_MUL(normal,B,C,tt3);
		}
		double norm = sqrt(DOT(normal,normal));
		normal[0] /= norm;
		normal[1] /= norm;
		normal[2] /= norm;
	}
	// assign the normal as the plane normal
	// assumes that the plane normal is normalized
	else
		ASSIGN(normal,planes[planeID]);

	return flag;
}


bool MeshSurface::inTriangle(double P[3], double A[3], double B[3], double C[3])
{
	// Compute vectors
	double v0[3],v1[3],v2[3],dot00,dot01,dot02,dot11,dot12,u,v;
	SUB(v0,C,A);
	SUB(v1,B,A);
	SUB(v2,P,A);

	// Compute dot products
	dot00 = DOT(v0,v0);
	dot01 = DOT(v0,v1);
	dot02 = DOT(v0,v2);
	dot11 = DOT(v1,v1);
	dot12 = DOT(v1,v2);

	// Compute barycentric coordinates
	double invDenom = dot00 * dot11 - dot01 * dot01;
	u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	if ((u > 0) && (v > 0) && (u + v < 1))
		return true;
	else
		return false;
}


void MeshSurface::point2plane(double p[3], double w[4], double *dist, double proj[3])
{
	double den = (DOT(w,w));
	double d = sqrt(den);
	double val = DOT(w,p) + w[3];

	proj[0] = p[0] - w[0]*(val/den);
	proj[1] = p[1] - w[1]*(val/den);
	proj[2] = p[2] - w[2]*(val/den);
	*dist = fabs(val/d);
}

bool MeshSurface::getProjection(double p[3], double *proj1, double *proj2,
								double *proj3, double *normal1, double *normal2, double *normal3)
{
	// get the triangles that are associated to this grid point
	// by querying the auxiliary grid
	#if !defined(FLOAT_VERTICES)
	double *pp[3];
	#else
	double pp[3][3];
	#endif
	double dist;
	set<int> triangles;
	// double hside = delphi->hside;

	// move from delphi grid to auxiliary grid
	int64_t irefx = (int64_t)rintp((p[0]-xmin)*scale);
	int64_t irefy = (int64_t)rintp((p[1]-ymin)*scale);
	int64_t irefz = (int64_t)rintp((p[2]-zmin)*scale);

	// double epsax  = p[0]+hside, epsay  = p[1]+hside, epsaz  = p[2]+hside;
	// double epsaxm = p[0]-hside, epsaym = p[1]-hside, epsazm = p[2]-hside;

	/* #if !defined(OPTIMIZE_GRIDS)
	for (int i=0; i<ind[irefz][irefy][irefx]; i++)
		triangles.insert((GRID_TRIANGLE_MAP(irefx,irefy,irefz,i,nx,ny,nz)));
	#else */
	for (int i=0; i<gridTriangleMap[ irefz*ny*nx + irefy*nx + irefx ].size(); i++)
		triangles.insert( gridTriangleMap[ irefz*ny*nx + irefy*nx + irefx ][i] );
	// #endif

	// keep the nearest triangle
	double locProj[3], minDist=INFINITY, locNorm[3];
	bool ff = false;

	// printf("\n\n Dist \n");
	for (set<int>::iterator it = triangles.begin(); it != triangles.end(); it++)
	{
		#if !defined(FLOAT_VERTICES)
		pp[0] = vertMatrix[faceMatrix[*it][0]];
		pp[1] = vertMatrix[faceMatrix[*it][1]];
		pp[2] = vertMatrix[faceMatrix[*it][2]];

		bool flag = point2triangle(p,pp[0],pp[1],pp[2],planes[*it],locProj,&dist,locNorm,*it);
		#else
		pp[0][0] = vertMatrix[faceMatrix[*it][0]][0];
		pp[0][1] = vertMatrix[faceMatrix[*it][0]][1];
		pp[0][2] = vertMatrix[faceMatrix[*it][0]][2];
		pp[1][0] = vertMatrix[faceMatrix[*it][1]][0];
		pp[1][1] = vertMatrix[faceMatrix[*it][1]][1];
		pp[1][2] = vertMatrix[faceMatrix[*it][1]][2];
		pp[2][0] = vertMatrix[faceMatrix[*it][2]][0];
		pp[2][1] = vertMatrix[faceMatrix[*it][2]][1];
		pp[2][2] = vertMatrix[faceMatrix[*it][2]][2];

		double plane[3];
		plane[0] = planes[*it][0];
		plane[1] = planes[*it][1];
		plane[2] = planes[*it][2];

		bool flag = point2triangle(p,pp[0],pp[1],pp[2],plane,locProj,&dist,locNorm,*it);
		#endif

		// printf("%lf %lf\n", dist, minDist);

		if (dist < minDist)
		{
			minDist = dist;
			*proj1 = locProj[0];
			*proj2 = locProj[1];
			*proj3 = locProj[2];

			*normal1 = locNorm[0];
			*normal2 = locNorm[1];
			*normal3 = locNorm[2];

			ff = flag;
		}
	}
	if (triangles.size() == 0)
	{
		{
			#ifdef ENABLE_BOOST_THREADS
			boost::mutex::scoped_lock scopedLock(mutex);
			#endif
			(*errorStream) << endl << WARN << "Approximating bgp with grid point";
		}
		*proj1 = p[0];
		*proj2 = p[1];
		*proj3 = p[2];
		return true;
	}

	// further check that the projection is in the cube
	// if (!testInCube((*proj1),(*proj2),(*proj3),p[0],p[1],p[2],delphi->side))
	//	cout << endl << WARN << "Out of cube projection in MeshSurface::getProjection!";
	return ff;
}


/** save in OFF format*/
bool MeshSurface::save(char *fileName)
{
	ofstream fout;
	fout.open(fileName,ios::out);

	cout << endl << INFO << "Writing mesh in OFF file format in " << fileName << "...";

	if (fout.fail())
	{
		cout << endl << WARN << "Cannot write file " << fileName;
		return false;
	}

	if (vertMatrix == NULL || faceMatrix == NULL)
	{
		cout << endl << WARN << "Cannot write null mesh!";
		return false;
	}

	fout << "OFF" << endl;
	fout << "# File created by "<< PROGNAME << " version " << VERSION << endl ;
	fout << endl;
	fout << numVertexes << " " << numTriangles << " 0" << endl;
	for (int i=0; i<numVertexes; i++)
		fout << vertMatrix[i][0] << " " << vertMatrix[i][1] << " " << vertMatrix[i][2] << endl;
	for (int i=0; i<numTriangles; i++)
		fout << "3 " << faceMatrix[i][0] << " " << faceMatrix[i][1] << " " << faceMatrix[i][2] << endl;

	return true;
}


bool MeshSurface::load(char *fileName)
{
	if ((int)strlen(fileName) == 0)
	{
		cout << WARN << "Cannot load with empty file name!";
		return false;
	}
	string ss(fileName);
	size_t found;
	ss = toLowerCase(ss);
	found = ss.find(".off");

	cout << endl << INFO << "Loading mesh in file " << fileName << "...";

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


bool MeshSurface::loadPLY(char *fname)
{
	int format=0, voh, foh, vph, fph;
	int nv,nt,i,j,i1,i2,i3,i4;
	float x,y,z;
	// bool triangulate = 0;
	FILE *in;
	char keyword[64], formats[24], version[10];

	if ((in = fopen(fname,"rb")) == NULL)
	{
		cout << endl << ERR <<("Can't open input ply file\n");
		return false;
	}

	if (strcmp(readLineFromFile(in),"ply"))
	{
		cout << endl << ERR <<("Input doesn't seem a valid ply file.\n");
		return false;
	}
	if (sscanf(readLineFromFile(in), "%7s %24s %10s", keyword,formats,version) < 3)
	{
		cout << endl << ERR <<("Unexpected token or end of file!\n");
		return false;
	}
	if (strcmp(keyword, "format"))
	{
		cout << endl << ERR <<("format definition expected!\n");
		return false;
	}

	if (!strcmp(formats, "ascii")) format = PLY_FORMAT_ASCII;
	else if (!strcmp(formats, "binary_little_endian")) format = PLY_FORMAT_BIN_L;
	else if (!strcmp(formats, "binary_big_endian")) format = PLY_FORMAT_BIN_B;
	else
	{
		cout << endl << ERR <<("Unrecognized format '%s'\n",formats);
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

	if (!sscanf(readLineFromFile(in), "%64s ",keyword))
	{
		cout << endl << ERR <<("Unexpected token or end of file!\n");
		return false;
	}
	while (strcmp(keyword, "end_header"))
		if (!sscanf(readLineFromFile(in), "%64s ",keyword))
		{
			cout << endl << ERR <<("Unexpected token or end of file!\n");
			return false;
		}

	if (vertMatrix != NULL)
		deleteMatrix2D<VERTEX_TYPE>(numVertexes,3,vertMatrix);

	if (vertexTrianglesList != NULL)
	{
		for (int i=0; i<numVertexes; i++)
			delete vertexTrianglesList[i];
		deleteVector<vector<int>*>(vertexTrianglesList);
	}

	numVertexes = nv;

	vertMatrix = allocateMatrix2D<VERTEX_TYPE>(numVertexes,3);

	vertexTrianglesList = allocateVector<vector<int>*>(numVertexes);

	for (int i=0; i<numVertexes; i++)
		vertexTrianglesList[i] = new vector<int>();

	for (i=0; i<nv; i++)
	{
		ply_readVCoords(in, format, vph, voh, &x, &y, &z);
		vertMatrix[i][0] = x;
		vertMatrix[i][1] = y;
		vertMatrix[i][2] = z;
	}

	if (faceMatrix != NULL)
		deleteMatrix2D<int>(numTriangles,3,faceMatrix);

	numTriangles = nt;
	faceMatrix = allocateMatrix2D<int>(numTriangles,3);

	for (i=0; i<nt; i++)
	{
		if (ply_readFIndices(in, format, fph, &i4, &i1, &i2, &i3))
		{
			if (i1<0 || i2<0 || i3<0 || i4<3 || i4>=4 || i1>(nv-1) || i2>(nv-1) || i3>(nv-1))
			{
				cout << endl << ERR <<("\nloadPLY: Invalid index at face %d!\n",i);
				return false;
			}
			for (j=3; j <= i4; j++)
			{
				if (i1 == i2 || i2 == i3 || i3 == i1) cout << endl << WARN <<("\nloadPLY: Coincident indexes at triangle %d! Skipping.\n",i);
				else
				{
					faceMatrix[i][0] = i1;
					faceMatrix[i][1] = i2;
					faceMatrix[i][2] = i3;
					vertexTrianglesList[i1]->push_back(i);
					vertexTrianglesList[i2]->push_back(i);
					vertexTrianglesList[i3]->push_back(i);
				}
				i2 = i3;
				if (j < i4)
				{
					if (!ply_readAnotherFIndex(in, format, &i3))
					{
						cout << endl << ERR <<("\nloadPLY: Couldn't read indexes for face # %d\n",i);
						return false;
					}
					// else triangulate = 1;
				}
				else ply_readOverhead(in, format, foh);
			}
		}
		else cout << endl << ERR <<("\nloadPLY: Couldn't read indexes for face # %d\n",i);
	}
	fclose(in);

	if (checkDuplicatedVertices)
		checkDuplicates();
	printSummary();
	computeNormals = true;

	return true;
}


bool MeshSurface::loadOFF(char *fileName)
{
	ifstream fin;
	fin.open(fileName,ios::in);
	int temp;
	char buffer[BUFLEN];

	if (fin.fail())
	{
		cout << endl << WARN << "Cannot read file " << fileName;
		return false;
	}

	int tag = 0;
	int header = 0;

	if (faceMatrix != NULL)
		deleteMatrix2D<int>(numTriangles,3,faceMatrix);

	if (vertMatrix != NULL)
		deleteMatrix2D<VERTEX_TYPE>(numVertexes,3,vertMatrix);

	if (vertexTrianglesList != NULL)
	{
		for (int i=0; i<numVertexes; i++)
			delete vertexTrianglesList[i];
		deleteVector<vector<int>*>(vertexTrianglesList);
	}

	while (1)
	{
		fin.getline(buffer,BUFLEN);

		if (fin.eof())
			break;

		//cout << endl << buffer;

		// skip empty lines
		if (strlen(buffer) <= 1)
			continue;

		if (tag == 0 && strncmp(buffer, "OFF", 1) == 0)
		{
			tag = 1;
			continue;
		}

		if (buffer[0] == '#')
		{
			cout << endl << INFO << "Mesh comment line: " << buffer;
			continue;
		}

		// if OFF is read then the header can be read
		if (header == 0 && tag == 1)
		{
			sscanf(buffer,"%d %d %d",&numVertexes,&numTriangles,&temp);
			if (numVertexes <=0 || numTriangles <= 0)
			{
				cout << endl << WARN << "Number of triangles or vertices <=0 !";
				return false;
			}
			else
			{
				faceMatrix = allocateMatrix2D<int>(numTriangles,3);

				vertMatrix = allocateMatrix2D<VERTEX_TYPE>(numVertexes,3);

				vertexTrianglesList = allocateVector<vector<int>*>(numVertexes);
				for (int i=0; i<numVertexes; i++)
					vertexTrianglesList[i] = new vector<int>();
			}
			header = 1;
			break;
		}
	}

	if (header !=1 || tag !=1)
	{
		cout << endl << WARN << "Cannot read OFF header or number of vertices/triangles, stop reading";
		cout << endl << WARN << "Tag " << tag << " Header " << header;
		return false;
	}
	// read vertices
	for (int i=0; i<numVertexes; i++)
	{
		fin.getline(buffer,BUFLEN);
		if (strlen(buffer) <= 1)
		{
			i--;
			continue;
		}
		float v1,v2,v3;
		sscanf(buffer, "%f %f %f", &v1,&v2,&v3);
		vertMatrix[i][0] = v1;
		vertMatrix[i][1] = v2;
		vertMatrix[i][2] = v3;
	}

	for (int i=0; i<numTriangles; i++)
	{
		fin.getline(buffer,BUFLEN);
		if (strlen(buffer) <= 1)
		{
			i--;
			continue;
		}

		sscanf(buffer, "%d %d %d %d", &temp,&faceMatrix[i][0],&faceMatrix[i][1],&faceMatrix[i][2]);

		vertexTrianglesList[faceMatrix[i][0]]->push_back(i);
		vertexTrianglesList[faceMatrix[i][1]]->push_back(i);
		vertexTrianglesList[faceMatrix[i][2]]->push_back(i);

		if (temp != 3)
		{
			cout << endl << WARN << "Non triangular mesh, stop reading";
			return false;
		}
	}

	if (checkDuplicatedVertices)
		checkDuplicates();
	printSummary();
	computeNormals = true;

	return true;
}


bool MeshSurface::loadMSMS(char *fileName,int numFiles)
{
	ifstream fin,fin2;

	char buffer[BUFLEN];
	char baseName[BUFLEN];

	char currentFace[2*BUFLEN];
	char currentVert[2*BUFLEN];

	strcpy(baseName,fileName);

	int tempNumVertices = 0;
	int tempNumTri = 0;
	numVertexes = 0;
	numTriangles = 0;

	int lastVert = 0;
	int lastTri = 0;

	#if !defined(FLOAT_VERTICES)
	if (vertMatrix != NULL)
		deleteMatrix2D<double>(numVertexes,3,vertMatrix);

	if (vertNormals!=NULL)
		deleteMatrix2D<double>(numVertexes,3,vertNormals);
	#else
	if (vertMatrix != NULL)
		deleteMatrix2D<float>(numVertexes,3,vertMatrix);

	if (vertNormals!=NULL)
		deleteMatrix2D<float>(numVertexes,3,vertNormals);
	#endif

	if (faceMatrix!=NULL)
		deleteMatrix2D<int>(numTriangles,3,faceMatrix);

	// getting the total number of vertices and triangles
	for (int i=0; i<numFiles; i++)
	{
		// main component
		if (i == 0)
		{
			sprintf(currentFace, "%s.face", baseName);
			sprintf(currentVert, "%s.vert", baseName);
		}
		// cavities
		else
		{
			sprintf(currentFace, "%s_%d.face", baseName,i);
			sprintf(currentVert, "%s_%d.vert", baseName,i);
		}

		fin.open(currentVert,ios::in);
		fin2.open(currentFace,ios::in);
		cout << endl << INFO << "Getting num vertices and triangles in MSMS mesh files " << currentFace << "," << currentVert << "...";

		if (fin.fail() || fin2.fail())
		{
			cout << WARN << "One or both MSMS files don't exist";
			return false;
		}

		// comment or empty
		bool comment = true;
		while (comment)
		{
			fin.getline(buffer,BUFLEN);
			if (buffer[0] != '#' && buffer[0] != '\0')
				comment = false;
		}
		sscanf(buffer, "%d", &tempNumVertices);
		numVertexes += tempNumVertices;
		fin.close();

		// reads num triangles
		comment = true;
		while (comment)
		{
			fin2.getline(buffer,BUFLEN);
			if (buffer[0] != '#' && buffer[0] != '\0')
				comment = false;
		}
		sscanf(buffer, "%d", &tempNumTri);
		numTriangles += tempNumTri;
		fin2.close();
	}

	cout << endl << INFO << "Total number of triangles " << numTriangles;
	cout << endl << INFO << "Total number of vertices " << numVertexes;

	vertMatrix = allocateMatrix2D<VERTEX_TYPE>(numVertexes,3);
	vertNormals = allocateMatrix2D<VERTEX_TYPE>(numVertexes,3);
	faceMatrix = allocateMatrix2D<int>(numTriangles,3);

	// read all
	for (int i=0; i<numFiles; i++)
	{
		// main component
		if (i == 0)
		{
			sprintf(currentFace, "%s.face", baseName);
			sprintf(currentVert, "%s.vert", baseName);
		}
		// cavities
		else
		{
			sprintf(currentFace, "%s_%d.face", baseName,i);
			sprintf(currentVert, "%s_%d.vert", baseName,i);
		}

		fin.open(currentVert,ios::in);
		fin2.open(currentFace,ios::in);
		cout << endl << INFO << "Loading MSMS mesh files " << currentFace << "," << currentVert << "...";

		if (fin.fail() || fin2.fail())
		{
			cout << WARN << "One or both MSMS files don't exist";
			return false;
		}

		bool comment = true;
		while (comment)
		{
			fin.getline(buffer,BUFLEN);
			if (buffer[0] != '#' && buffer[0] != '\0')
				comment = false;
		}

		sscanf(buffer, "%d", &tempNumVertices);

		for (int i=lastVert; i < lastVert + tempNumVertices; i++)
		{
			fin.getline(buffer,BUFLEN);
			// skip empty lines
			if (strlen(buffer) > 1)
			{
				float a,b,c,d,e,f;
				sscanf(buffer, "%f %f %f %f %f %f", &a,&b,&c,&d,&e,&f);
				vertMatrix[i][0] = a;
				vertMatrix[i][1] = b;
				vertMatrix[i][2] = c;
				vertNormals[i][0] = d;
				vertNormals[i][1] = e;
				vertNormals[i][2] = f;
			}
			else
				i--;
		}
		fin.close();

		comment = true;
		while (comment)
		{
			fin2.getline(buffer,BUFLEN);
			if (buffer[0]!='#' && buffer[0]!='\0')
				comment = false;
		}
		sscanf(buffer,"%d",&tempNumTri);

		for (int i=lastTri; i < lastTri + tempNumTri; i++)
		{
			fin2.getline(buffer,BUFLEN);
			// skip empty lines
			if (strlen(buffer) > 1)
			{
				sscanf(buffer, "%d %d %d",
					   &(faceMatrix[i][0]),&(faceMatrix[i][1]),&(faceMatrix[i][2]));
				faceMatrix[i][0]--;
				faceMatrix[i][1]--;
				faceMatrix[i][2]--;
				faceMatrix[i][0] += lastVert;
				faceMatrix[i][1] += lastVert;
				faceMatrix[i][2] += lastVert;
			}
			else
				i--;
		}
		fin2.close();

		lastTri += tempNumTri;
		lastVert += tempNumVertices;
	}
	if (checkDuplicatedVertices)
		checkDuplicates();
	printSummary();

	// MDM: bug fixed, computeNormals restored to the original value
	bool tempNormals = computeNormals;
	computeNormals = false;

	preProcessTriangles();

	computeNormals = tempNormals;

	return true;
}


bool MeshSurface::checkDuplicates()
{
	if (faceMatrix == NULL || vertMatrix == NULL)
	{
		cout << endl << WARN << "Cannot check mesh without a loaded mesh!";
		return false;
	}

	bool ret = true;

	// check for duplicated vertices
	set<pair<VERTEX_TYPE,pair<VERTEX_TYPE,VERTEX_TYPE>>> checkV;
	set<pair<VERTEX_TYPE,pair<VERTEX_TYPE,VERTEX_TYPE>>>::iterator stv;

	for (int i=0; i<numVertexes; i++)
	{
		pair<VERTEX_TYPE,pair<VERTEX_TYPE,VERTEX_TYPE>> triplet;

		triplet.first = vertMatrix[i][0];
		triplet.second.first = vertMatrix[i][1];
		triplet.second.second = vertMatrix[i][2];
		stv = checkV.find(triplet);

		if (stv == checkV.end())
		{
			// ok
		}
		// duplicated vertex
		else
		{
			cout << endl << WARN << "Duplicated vertex detected! " << triplet.first << " "
			<< triplet.second.first << " " << triplet.second.second;
			ret = false;
		}

		checkV.insert(triplet);
	}

	// check for duplicated triangles
	set<pair<int,pair<int,int>>> checkT;
	set<pair<int,pair<int,int>>>::iterator st;

	for (int i=0; i<numTriangles; i++)
	{
		pair<int,pair<int,int>> triplet;
		triplet.first = faceMatrix[i][0];
		triplet.second.first = faceMatrix[i][1];
		triplet.second.second = faceMatrix[i][2];
		st = checkT.find(triplet);

		if (st == checkT.end())
		{
			// ok
		}
		// duplicated triangle
		else
		{
			cout << endl << WARN << "Duplicated triangle detected! " << triplet.first << " "
			<< triplet.second.first << " " << triplet.second.second;
			ret = false;
		}

		checkT.insert(triplet);
	}

	return ret;
}


void MeshSurface::printSummary()
{
	//cout << endl << INFO << "Summary of the triangulated mesh:";
	if (faceMatrix == NULL || vertMatrix == NULL)
	{
		cout << endl << WARN << "Mesh not loaded!";
	}
	else
	{
		cout << endl << INFO << "Number of loaded vertices -> " << numVertexes;
		cout << endl << INFO << "Number of loaded triangles -> " << numTriangles;
	}
}
