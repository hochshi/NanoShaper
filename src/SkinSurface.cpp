#include "SkinSurface.h"


void SkinSurface::clear()
{
	// remove 2d grid for ray casting
	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridMixedCellMap != NULL)
		deleteVector<int>(gridMixedCellMap);
	
	if (ind != NULL)
		deleteMatrix3D<int>(nx,ny,nz,ind);
	
	if (gridMixedCellMap2D != NULL)
		deleteVector<int>(gridMixedCellMap2D);

	if (ind_2d != NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,last_cols_ind,ind_2d);
	#else */
	if (gridMixedCellMap != NULL)
	{
		for (int64_t i=0; i < nx*ny*nz; i++)
		{
			gridMixedCellMap[i].clear();
		}
		delete[] gridMixedCellMap;
		gridMixedCellMap = NULL;
	}
	if (gridMixedCellMap2D != NULL)
	{
		for (int64_t i=0; i < last_rows_ind*last_cols_ind; i++)
		{
			gridMixedCellMap2D[i].clear();
		}
		delete[] gridMixedCellMap2D;
		gridMixedCellMap2D = NULL;
	}
	// #endif

	if (x != NULL)
		deleteVector<double>(x);
	if (y != NULL)
		deleteVector<double>(y);
	if (z != NULL)
		deleteVector<double>(z);


	#if defined(MULTITHREADED_SKIN_BUILDING)

	// if (thread_data_wrapper != NULL)
	{
		for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
		{
			ThreadDataWrapper *tdw = &thread_data_wrapper[thread_id];

			for (int i=0; i<tdw->mixedComplex.size(); i++)
			{
				if (tdw->mixedComplex[i] == NULL)
					continue;

				if (tdw->mixedComplex[i]->quadric != NULL)
					deleteVector<double>(tdw->mixedComplex[i]->quadric);
			}
			for (int i=0; i<tdw->mixedComplex.size(); i++)
			{
				delete tdw->mixedComplex[i];
				tdw->mixedComplex[i] = NULL;
			}
			tdw->mixedComplex.clear();
		}
		// delete thread_data_wrapper;
		// thread_data_wrapper = NULL;
	}

	#endif

	for (int i=0; i<numMixedCells; i++)
	{
		mixedComplex[i]->mc_points.clear();

		delete mixedComplex[i];
		mixedComplex[i] = NULL;
	}

	mixedComplex.clear();

	if (patchBasedAlgorithm && num_pixel_intersections != NULL)
	{
		// remove bufferse employed in the patch-based ray tracer

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
		for (int i=0; i < numMixedCells*numPanels; i++)
		{
			temp_intersections_buffer[i].clear();
		}
		delete[] temp_intersections_buffer;
		temp_intersections_buffer = NULL;

		for (int i=0; i < numMixedCells*numPanels; i++)
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


SkinSurface::~SkinSurface()
{
	clear();
}


void SkinSurface::init()
{
	gridMixedCellMap = NULL;
	gridMixedCellMap2D = NULL;

	#ifdef MULTITHREADED_SKIN_BUILDING
	// thread_data_wrapper = NULL;
	#endif
	// mixedComplex = nullptr;
	mixedComplex.clear();

	numMixedCells = 0;
	/* #if !defined(OPTIMIZE_GRIDS)
	ind = NULL;
	ind_2d = NULL;
	#endif */
	x = NULL;
	y = x;
	z = y;
	s = 0.5;
	savePovRay = false;
	// aggressive settings can be used in a 64 bit machine with 6 GB of memory
	AUX_GRID_DIM_SKIN = 100;
	MAX_MIXED_CELLS = 400;
	AUX_GRID_DIM_SKIN_2D = 50; // aggressive setting is 150
	MAX_MIXED_CELLS_2D = (400*AUX_GRID_DIM_SKIN_2D); // aggressive setting is (200*AUX_GRID_DIM_SKIN_2D)
	surfType = MOLECULAR_SURFACE;
	fastProjection = false;
	providesAnalyticalNormals = true;
	
	if (patchBasedAlgorithm)
	{
		num_pixel_intersections = NULL;
		pixel_intersections = NULL;
	}
}


void SkinSurface::init(ConfigFile *cf)
{
	double skin_s = cf->read<double>( "Skin_Surface_Parameter", 0.45 );
	unsigned int maxSkinDim = cf->read<unsigned int>( "Max_skin_patches_auxiliary_grid_size", 100 );
	unsigned int maxSkinPatches = cf->read<unsigned int>( "Max_skin_patches_per_auxiliary_grid_cell", 400 );
	unsigned int maxSkinDim2D = cf->read<unsigned int>( "Max_skin_patches_auxiliary_grid_2d_size", 50 );
	unsigned int maxSkinPatches2D = cf->read<unsigned int>( "Max_skin_patches_per_auxiliary_grid_2d_cell", 400 );
	bool useFastProjection = cf->read<bool>("Skin_Fast_Projection",false);
	bool savePovRay = cf->read<bool>("Save_PovRay",false);

	setAuxGrid(maxSkinDim,maxSkinPatches);
	setAuxGrid2D(maxSkinDim2D,maxSkinPatches2D);
	setFastProjection(useFastProjection);
	inside = 5;
	setShrinking(skin_s);
	setSavePovRay(savePovRay);
}


SkinSurface::SkinSurface(DelPhiShared *ds):Surface()
{
	init();
	// set environment
	delphi = ds;
}


SkinSurface::SkinSurface():Surface()
{
	delphi = NULL;
	init();
}


SkinSurface::SkinSurface(ConfigFile *cf,DelPhiShared *ds):Surface(cf)
{
	init();
	init(cf);
	// set environment
	delphi = ds;
}


void SkinSurface::postRayCasting()
{
	// remove 2d grid for ray casting
	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridMixedCellMap2D != NULL)
		deleteVector<int>(gridMixedCellMap2D);

	if (ind_2d != NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,last_cols_ind,ind_2d);
	#else */
	if (gridMixedCellMap2D != NULL)
	{
		for (int64_t i=0; i < last_rows_ind*last_cols_ind; i++)
		{
			gridMixedCellMap2D[i].clear();
		}
		delete[] gridMixedCellMap2D;
		gridMixedCellMap2D = NULL;
	}
	// #endif
	

	if (num_pixel_intersections != NULL)
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
		for (int i=0; i < numMixedCells*numPanels; i++)
		{
			temp_intersections_buffer[i].clear();
		}
		delete[] temp_intersections_buffer;
		temp_intersections_buffer = NULL;

		for (int i=0; i < numMixedCells*numPanels; i++)
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


bool SkinSurface::preBoundaryProjection()
{
	// 3d grid is necessary only for boundary grid points projection
	if (projBGP)
		return buildAuxiliaryGrid();
	return false;
}


int SkinSurface::getNumPatches (void)
{
	return numMixedCells;
}


bool SkinSurface::isPatchBasedRayTracingSupported (void)
{
	return true;
}


bool SkinSurface::build()
{
	bool f = false;
	#ifdef ENABLE_CGAL
		f = buildSkinCGAL();
	#else
		cout << endl << ERR << "Skin surface cannot be used, please install CGAL and rebuild\n";
		exit(-1);
	#endif
	if (!f)
	{
		cout << endl << ERR << "Error during skin build-up";
	}
	return f;
}


#ifdef ENABLE_CGAL

void SkinSurface::BuildDelaunayTetraCells (Rt &rT, vector<MixedCell*> &mixedComplex, int num_cells[])
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// Compute Voronoi points and reduced Delaunay solids (Skin) ///////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	#if defined(QUADRIC_CULLING)
	double checkWeights[4];
	#endif

	double dummy_atom_weight = -1.0;

	int currentCell = 0;

	for (Finite_Cells_Iterator fcit = rT.finite_cells_begin();
		 fcit != rT.finite_cells_end(); fcit++)
	{
		const Point3 &p = rT.geom_traits().construct_weighted_circumcenter_3_object()(fcit->vertex(0)->point(),fcit->vertex(1)->point(),fcit->vertex(2)->point(),fcit->vertex(3)->point());

		Del3Cell *&mc = (Del3Cell*&)fcit->info();

		mc = new Del3Cell();
		mc->surface_type = DELAUNAY_TETRA_CELL;
		mc->vor[0] = p.x();
		mc->vor[1] = p.y();
		mc->vor[2] = p.z();

		mixedComplex.push_back(mc);

		double *vor = mc->vor;

		// minkowsky sum of tethraedron vertices and voronoi point. Record the atom from which
		// it was computed
		// double points[4][3];
		double a, b, c;
		double radius2 = INFINITY;

		for (int i=0; i<4; i++)
		{
			Weighted_point &wp = fcit->vertex(i)->point();
			#if defined(QUADRIC_CULLING)
			checkWeights[i] = wp.weight();
			#endif

			// nopatch
			mc->ids[i] = fcit->vertex(i)->info();

			mc->reduced[i][0] = wp.x();
			mc->reduced[i][1] = wp.y();
			mc->reduced[i][2] = wp.z();

			// points[i][0] = wp.x();
			// points[i][1] = wp.y();
			// points[i][2] = wp.z();

			a = vor[0] - mc->reduced[i][0];
			mc->reduced[i][0] += s*a;

			b = vor[1] - mc->reduced[i][1];
			mc->reduced[i][1] += s*b;

			c = vor[2] - mc->reduced[i][2];
			mc->reduced[i][2] += s*c;

			if (i == 0)
			{
				double norm2 = a*a+b*b+c*c;
				radius2 = -(wp.weight()-norm2);
			}
		}

		// avoid storing non sense patch, but keep its voronoi point
		if (radius2 < 0)
		{
			mc->surface_type = SKIP_CELL;
			continue;
		}
		#if defined(QUADRIC_CULLING)
		// avoid storing nonsense patches, but it keeps the reduced points
		if (checkWeights[0] == dummy_atom_weight ||
			checkWeights[1] == dummy_atom_weight ||
			checkWeights[2] == dummy_atom_weight ||
			checkWeights[3] == dummy_atom_weight)
		{
			mc->surface_type = SKIP_CELL;
			continue;
		}
		#endif
		////////////////////////// store equation and simplex ////////////////////////////////////////////

		plane3points(mc->reduced[2],mc->reduced[1],mc->reduced[0],mc->planes[0],false);
		plane3points(mc->reduced[3],mc->reduced[2],mc->reduced[0],mc->planes[1],false);
		plane3points(mc->reduced[0],mc->reduced[1],mc->reduced[3],mc->planes[3],false);
		plane3points(mc->reduced[3],mc->reduced[1],mc->reduced[2],mc->planes[2],false);

		// plane3points(points[2],points[1],points[0],mc->planes_or[0],false);
		// plane3points(points[3],points[2],points[0],mc->planes_or[1],false);
		// plane3points(points[0],points[1],points[3],mc->planes_or[3],false);
		// plane3points(points[3],points[1],points[2],mc->planes_or[2],false);

		mc->mc_points.push_back(mc->reduced[0]);
		mc->mc_points.push_back(mc->reduced[1]);
		mc->mc_points.push_back(mc->reduced[2]);
		mc->mc_points.push_back(mc->reduced[3]);

		// NOTE:
		// these planes are already correctly oriented because the choosen
		// order of points in producing the planes leverages the combinatoric properties CGAL internal data structures.
		// In a few words, using points on a tethraedron in that order one is assured to get correctly oriented planes
		// In a non CGAL implementation these planes must be oriented

		double *QQ = allocateVector<double>(10);
		mc->quadric = QQ;

		QQ[0] = 1;
		QQ[1] = 1;
		QQ[2] = 1;
		QQ[3] = -radius2*(1-s)+vor[0]*vor[0]+vor[1]*vor[1]+vor[2]*vor[2];
		QQ[4] = 0;
		QQ[7] = 0;
		QQ[9] = -vor[0];
		QQ[5] = 0;
		QQ[8] = -vor[1];
		QQ[6] = -vor[2];

		////////////////////////// end store equation and simplex ////////////////////////////////////////////

		// mixedComplex.push_back(mc);
		++currentCell;
	}
	num_cells[DELAUNAY_TETRA_CELL] = currentCell;
}


void SkinSurface::BuildDelaunayPointCells (Rt &rT, vector<MixedCell*> &mixedComplex,
										   #if !defined(NEW_ATOM_PATCHES)
										   Del0Cell **atomPatches,
										   #else
										   int atomPatches[],
										   #endif
										   int num_cells[])
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// Compute Voronoi (VDW) or reduced Voronoi solids (Skin) //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	vector<Cell_handle> cells;
	cells.reserve(1000);

	double dummy_atom_weight = -1.0;

	int cellCount = 0;

	for (Finite_Vertex_Iterator fvit = rT.finite_vertices_begin(); fvit != rT.finite_vertices_end(); fvit++)
	{
		const Weight &w0 = fvit->point().weight();

		// skip virtual atoms
		if (w0 == dummy_atom_weight)
			continue;

		if (rT.is_infinite(fvit->cell()))
			continue;

		// current atom id around which we are moving
		// nopatch
		const int &currentId = fvit->info();

		cells.clear();
		rT.incident_cells(fvit,std::back_inserter(cells));

		bool infiniteCell = false;

		// start check if all thetraedra are feasible
		for (std::vector<Cell_handle>::iterator it = cells.begin(); it != cells.end(); it++)
		{
			if (rT.is_infinite(*it))
			{
				infiniteCell = true;
				break;
			}
		}
		if (infiniteCell)
			continue;

		// this cell is finite, it is acceptable
		Del0Cell *mc = new Del0Cell();
		mc->id = currentId;
		mc->surface_type = DELAUNAY_POINT_CELL;

		#if !defined(NEW_ATOM_PATCHES)
		atomPatches[currentId] = mc;
		#else
		// Retrieve index of element within mixedComplex[*] for later purposes
		// (ixedComplex.size() without -1, since mixedComplex.push_back(mc) appears later)
		atomPatches[currentId] = mixedComplex.size();
		#endif

		double *QQ = allocateVector<double>(10);
		mc->quadric = QQ;

		QQ[0] = 1;
		QQ[1] = 1;
		QQ[2] = 1;
		QQ[3] = -s*(fvit->point().weight())+fvit->point().x()*fvit->point().x()+fvit->point().y()*fvit->point().y()+fvit->point().z()*fvit->point().z();
		QQ[4] = 0;
		QQ[7] = 0;
		QQ[9] = -fvit->point().x();
		QQ[5] = 0;
		QQ[8] = -fvit->point().y();
		QQ[6] = -fvit->point().z();

		////////////////////////// end store equation and simplex ////////////////////////////////////////////
		mixedComplex.push_back(mc);
		++cellCount;
	}
	num_cells[DELAUNAY_POINT_CELL] = cellCount;
}


void SkinSurface::BuildDelaunayEdgeCells (Rt &rT, vector<MixedCell*> &mixedComplex,
										  #if !defined(NEW_ATOM_PATCHES)
										  Del0Cell **atomPatches,
										  #else
										  int atomPatches[],
										  #endif
										  int num_atoms, int num_cells[], bool *any_error)
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// Compute Delaunay Edges + Voronoi facets /////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double *vorFace[3];

	double RM_[4][4];
	double T_[4][4];
	double TEMP_[4][4];
	// double C_[4][4];
	double Cinv_[4][4];
	double M_[4][4];
	double Q_[4][4];
	double sign;

	// Skin Surface needs the remaining solids
	int num_d1_v2_patches = 0;

	*any_error = true;

	double dummy_atom_weight = -1.0;

	// Delaunay Edge vs Voronoi facets

	vector<double*> upperPoints, lowerPoints;

	upperPoints.reserve(50);
	lowerPoints.reserve(50);

	for (Finite_Edges_Iterator fei = rT.finite_edges_begin(); fei != rT.finite_edges_end(); fei++)
	{
		// up and down edge vertices
		const Vertex_handle &v1 = fei->first->vertex(fei->second);
		const Vertex_handle &v2 = fei->first->vertex(fei->third);

		// connecting one or two non atoms
		if (v1->point().weight() == dummy_atom_weight ||
			v2->point().weight() == dummy_atom_weight)
			continue;

		Cell_circulator cit = rT.incident_cells(*fei);
		int num = (int)circulator_size(cit);

		#if !defined(NO_CGAL_PATCHING)
		const int refAtom1 = v1->point().index();
		const int refAtom2 = v2->point().index();
		#else
        // nopatch
		const int refAtom1 = v1->info();
		const int refAtom2 = v2->info();
		#endif

		double v1p[3], v2p[3];

		v2p[0] = v2->point().x();
		v2p[1] = v2->point().y();
		v2p[2] = v2->point().z();

		v1p[0] = v1->point().x();
		v1p[1] = v1->point().y();
		v1p[2] = v1->point().z();

		upperPoints.clear();
		lowerPoints.clear();

		Del1Cell *mc = new Del1Cell();
		mc->ids[0] = refAtom1;
		mc->ids[1] = refAtom2;
		mc->surface_type = DELAUNAY_EDGE_CELL;

		// collect all the points
		for (int i=0; i<num; i++)
		{
			const Cell_handle &ch = cit;
			const Del3Cell *ci = ch->info();

			double *p1 = NULL;
			double *p2 = NULL;

			for (int j=0; j<4; j++)
			{
				if (p1 == NULL && ci->ids[j] == refAtom1) {
					p1 = (double *)ci->reduced[j];
				}
				if (p2 == NULL && ci->ids[j] == refAtom2) {
					p2 = (double *)ci->reduced[j];
				}
			}

			if (p1 == NULL || p2 == NULL)
			{
				for (int j=0; j<4; j++)
					cout << endl << ERR << ci->ids[j];

				cout << endl << ERR << "Index not found";
				if (p1 == NULL) cout << endl << ERR << refAtom1;
				if (p2 == NULL) cout << endl << ERR << refAtom2;
				return;
			}
			// store the pointer to the reduced point just identified
			mc->mc_points.push_back(p1);
			mc->mc_points.push_back(p2);
			upperPoints.push_back(p1);
			lowerPoints.push_back(p2);

			// printf("\n vor %x", ci->vor);
			// store the vornoi points to infer voronoi original plane
			// this is later used to compute the focus of the cell
			if (i < 3) vorFace[i] = (double *)ci->vor;

			++cit;
		}
		// we have the three points to compute the upper and lower facets;
		// recall that the first point is up, the second down, the third up and so on
		plane3points(upperPoints[0],upperPoints[1],upperPoints[2],mc->upper);
		plane3points(lowerPoints[0],lowerPoints[1],lowerPoints[2],mc->lower);

		// now the corresponding reduced voronoi cells must link to these planes
		#if !defined(NEW_ATOM_PATCHES)
		if (refAtom1 < num_atoms)
			atomPatches[refAtom1]->planes.push_back(mc->upper);
		if (refAtom2 < num_atoms)
			atomPatches[refAtom2]->planes.push_back(mc->lower);
		#else
		if (refAtom1 < num_atoms)
		{
			Del0Cell *dc = (Del0Cell *)mixedComplex[ atomPatches[refAtom1] ];
			dc->planes.push_back(mc->upper);
		}
		if (refAtom2 < num_atoms)
		{
			Del0Cell *dc = (Del0Cell *)mixedComplex[ atomPatches[refAtom2] ];
			dc->planes.push_back(mc->lower);
		}
		#endif

		// plane line intersection point to get focus
		double mat[3][3];
		double invmat[3][3];

		mat[0][0] = v2p[0] - v1p[0];
		mat[1][0] = v2p[1] - v1p[1];
		mat[2][0] = v2p[2] - v1p[2];

		mat[0][1] =	vorFace[1][0] - vorFace[0][0];
		mat[1][1] = vorFace[1][1] - vorFace[0][1];
		mat[2][1] = vorFace[1][2] - vorFace[0][2];

		mat[0][2] =	vorFace[2][0] - vorFace[0][0];
		mat[1][2] = vorFace[2][1] - vorFace[0][1];
		mat[2][2] = vorFace[2][2] - vorFace[0][2];

		double det;

		INVERT_3X3(invmat,det,mat)

		double v[3];

		v[0] = v2p[0] - vorFace[0][0];
		v[1] = v2p[1] - vorFace[0][1];
		v[2] = v2p[2] - vorFace[0][2];

		double res = 0;
		for (int j=0; j<3; j++)
			res += invmat[0][j] * v[j];

		double focus[3];

		focus[0] = v2p[0] + res*(v1p[0]-v2p[0]);
		focus[1] = v2p[1] + res*(v1p[1]-v2p[1]);
		focus[2] = v2p[2] + res*(v1p[2]-v2p[2]);

		// generate a reference inner point
		double inner[3];

		inner[0] = (upperPoints[0][0]+upperPoints[1][0]+upperPoints[2][0]+lowerPoints[0][0]+lowerPoints[1][0]+lowerPoints[2][0])/6.;
		inner[1] = (upperPoints[0][1]+upperPoints[1][1]+upperPoints[2][1]+lowerPoints[0][1]+lowerPoints[1][1]+lowerPoints[2][1])/6.;
		inner[2] = (upperPoints[0][2]+upperPoints[1][2]+upperPoints[2][2]+lowerPoints[0][2]+lowerPoints[1][2]+lowerPoints[2][2])/6.;

		// decide orientation based on the position of the focus of the computed planes
		sign = DOT(mc->upper,inner) + mc->upper[3];
		if (sign>0) ASSIGN4(mc->upper, -mc->upper)

		sign = DOT(mc->lower,inner) + mc->lower[3];
		if (sign>0) ASSIGN4(mc->lower, -mc->lower)

		// build lateral planes
		for (int i=0; i<num-1; i++)
		{
			double *plane = allocateVector<double>(4);
			plane3points(upperPoints[i],lowerPoints[i],upperPoints[i+1],plane);
			// orient
			sign = DOT(plane,inner) + plane[3];
			if (sign>0) ASSIGN4(plane, -plane)
			mc->planes.push_back(plane);
		}

		// add the last one
		double *plane = allocateVector<double>(4);
		plane3points(upperPoints[num-1],lowerPoints[num-1],upperPoints[0],plane);
		// orient

		sign = DOT(plane,inner) + plane[3];
		if (sign>0) ASSIGN4(plane, -plane)

		mc->planes.push_back(plane);

		// compute the reduced cell radius
		double a = v1p[0] - focus[0];
		double b = v1p[1] - focus[1];
		double c = v1p[2] - focus[2];
		double radius2 = v1->point().weight() - (a*a+b*b+c*c);

		double o2[3], o1[3], o3[3], ttt;

		// vector parallel to voronoi face
		SUB(o2,vorFace[0],focus)
		NORMALIZE(o2,ttt)

		// vector parallel to delaunay edge
		SUB(o1,v2p,v1p)
		NORMALIZE(o1,ttt)

		CROSS(o3,o1,o2)
		NORMALIZE(o3,ttt)

		// o1,o2,o3 is an orthonormal set

		for (int ii=0;ii<4;ii++)
			for (int jj=0;jj<4;jj++)
			{
				// C_[ii][jj]		=	0;
				Cinv_[ii][jj]	=	0;
				M_[ii][jj]		=	0;
				Q_[ii][jj]		=	0;
				TEMP_[ii][jj]	=	0;
				RM_[ii][jj]		=	0;
				T_[ii][jj]		=	0;
			}
		RM_[0][0] = o1[0];
		RM_[1][0] = o1[1];
		RM_[2][0] = o1[2];

		RM_[0][1] = o3[0];
		RM_[1][1] = o3[1];
		RM_[2][1] = o3[2];

		RM_[0][2] = o2[0];
		RM_[1][2] = o2[1];
		RM_[2][2] = o2[2];

		RM_[3][3] = 1.;

		T_[0][0] = 1.;
		T_[1][1] = 1.;
		T_[2][2] = 1.;
		T_[0][3] = focus[0];
		T_[1][3] = focus[1];
		T_[2][3] = focus[2];
		T_[3][3] = 1.;

		Matrix4x4MultiplyBy4x4(T_,RM_,Cinv_);

		//for (int ii=0;ii<4;ii++)
		//	for (int jj=0;jj<4;jj++)
		//		C_[ii][jj]=Cinv_[ii][jj];

		inplace_invert4x4(Cinv_);

		M_[0][0] = -1./(1.-s);
		M_[1][1] = 1./s;
		M_[2][2] = 1./s;
		M_[3][3] = -radius2;

		Matrix4x4MultiplyBy4x4(M_,Cinv_,TEMP_);

		for (int ii = 0; ii<4; ii++)
			for (int jj = ii+1; jj<4; jj++)
			{
				double temp = Cinv_[ii][jj];
				Cinv_[ii][jj] = Cinv_[jj][ii];
				Cinv_[jj][ii] = temp;
			}
		Matrix4x4MultiplyBy4x4(Cinv_,TEMP_,Q_);

		////////////////////////// store equation and simplex ////////////////////////////////////////////

		double *QQ = allocateVector<double>(10);
		mc->quadric = QQ;

		// Symmetric matrix coding
		//	Q0	Q4	Q7	Q9
		//		Q1	Q5	Q8
		//			Q2	Q6
		//				Q3

		QQ[0] = Q_[0][0];
		QQ[1] = Q_[1][1];
		QQ[2] = Q_[2][2];
		QQ[3] = Q_[3][3];
		QQ[4] = Q_[0][1];
		QQ[7] = Q_[0][2];
		QQ[9] = Q_[0][3];
		QQ[5] = Q_[1][2];
		QQ[8] = Q_[1][3];
		QQ[6] = Q_[2][3];

		////////////////////////// end store equation and simplex ////////////////////////////////////////////

		mixedComplex.push_back(mc);
		++num_d1_v2_patches;
	}
	num_cells[DELAUNAY_EDGE_CELL] = num_d1_v2_patches;

	*any_error = false;
}


void SkinSurface::BuildDelaunayFacetCells (Rt &rT, vector<MixedCell*> &mixedComplex, int num_cells[], bool *any_error)
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////// Compute Delaunay Facets + Voronoi Edges /////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double RM_[4][4];
	double T_[4][4];
	double TEMP_[4][4];
	// double C_[4][4];
	double Cinv_[4][4];
	double M_[4][4];
	double Q_[4][4];
	// plane orientation test
	double sign;

	#if defined(QUADRIC_CULLING)
	double checkWeights[4];
	#endif

	double dummy_atom_weight = -1.0;

	int num_d2_v1_patches = 0;

	*any_error = true;

	for (Finite_Facets_Iterator ffi = rT.finite_facets_begin(); ffi != rT.finite_facets_end(); ffi++)
	{
		const Cell_handle &ch = ffi->first;
		// get the mirror facet/cell go on if infinite
		const Facet &mirrorFace = rT.mirror_facet(*ffi);
		const Cell_handle &ch2 = mirrorFace.first;

		if (rT.is_infinite(ch) || rT.is_infinite(ch2))
			continue;

		Del2Cell *mc = new Del2Cell();
		mc->surface_type = DELAUNAY_FACET_CELL;

		mixedComplex.push_back(mc);

		Del3Cell *ci = ch->info();
		Del3Cell *ci2 = ch2->info();

		double *upperPoints[3], *lowerPoints[3];

		upperPoints[0] = NULL;
		upperPoints[1] = NULL;
		upperPoints[2] = NULL;

		lowerPoints[0] = NULL;
		lowerPoints[1] = NULL;
		lowerPoints[2] = NULL;

		double delFace[3][3];
		double weight;

		// it is the local index of the tethraedron vertex
		int f_index1 = ffi->second;
		// int f_index2 = mirrorFace.second;

		int ind = 0;

		// recover the three reduced common vertices
		for (int i=0; i<4; i++)
		{
			// this vertex is opposite and cannot be shared, so skip
			if (i == f_index1)
				continue;

			#if defined(QUADRIC_CULLING)
			checkWeights[ind] = ch->vertex(i)->point().weight();
			#endif

			// nopatch
			const int index1 = ch->vertex(i)->info();

			mc->ids[ind] = index1;
			delFace[ind][0] = ch->vertex(i)->point().x();
			delFace[ind][1] = ch->vertex(i)->point().y();
			delFace[ind][2] = ch->vertex(i)->point().z();

			// store the weight of the del point that will be used to compute the focus
			if (ind == 0)
				weight = ch->vertex(i)->point().weight();

			double *p1 = NULL;
			double *p2 = NULL;

			for (int j=0; j<4; j++)
			{
				// search the corresponding reduced point
				if (p1 == NULL && ci->ids[j] == index1) {
					p1 = ci->reduced[j];
				}
				if (p2 == NULL && ci2->ids[j] == index1) {
					p2 = ci2->reduced[j];
				}
			}
			// weird..
			if (p1 == NULL || p2 == NULL)
			{
				return;
			}
			// store the pointer to the reduced point just identified
			mc->mc_points.push_back(p1);
			mc->mc_points.push_back(p2);

			upperPoints[ind] = p1;
			lowerPoints[ind] = p2;

			++ind;
		}
		#ifdef QUADRIC_CULLING
		// avoid storing nonsense patches, but it keeps the reduced points;
		if (checkWeights[0] == dummy_atom_weight ||
			checkWeights[1] == dummy_atom_weight ||
			checkWeights[2] == dummy_atom_weight)
		{
			// delete mc;
			// mc = NULL;
			mc->surface_type = SKIP_CELL;
			continue;
		}
		#endif

		double *vor1 = ci->vor;
		double *vor2 = ci2->vor;

		// plane line intersection point to get focus
		double mat[3][3];
		double invmat[3][3];

		mat[0][0] = vor2[0] - vor1[0];
		mat[1][0] = vor2[1] - vor1[1];
		mat[2][0] = vor2[2] - vor1[2];

		mat[0][1] =	delFace[1][0] - delFace[0][0];
		mat[1][1] = delFace[1][1] - delFace[0][1];
		mat[2][1] = delFace[1][2] - delFace[0][2];

		mat[0][2] =	delFace[2][0] - delFace[0][0];
		mat[1][2] = delFace[2][1] - delFace[0][1];
		mat[2][2] = delFace[2][2] - delFace[0][2];

		double det;
		INVERT_3X3(invmat,det,mat)

		double v[3];

		v[0] = vor2[0] - delFace[0][0];
		v[1] = vor2[1] - delFace[0][1];
		v[2] = vor2[2] - delFace[0][2];

		double res = 0;
		for (int j=0; j<3; j++)
			res += invmat[0][j]*v[j];

		double focus[3];

		focus[0] = vor2[0] + res*(vor1[0]-vor2[0]);
		focus[1] = vor2[1] + res*(vor1[1]-vor2[1]);
		focus[2] = vor2[2] + res*(vor1[2]-vor2[2]);

		// compute the reduced cell radius
		double a = delFace[0][0]-focus[0];
		double b = delFace[0][1]-focus[1];
		double c = delFace[0][2]-focus[2];
		double radius2 = weight - (a*a+b*b+c*c);

		// generate a reference inner point
		double inner[3];

		inner[0] = (upperPoints[0][0]+upperPoints[1][0]+upperPoints[2][0]+lowerPoints[0][0]+lowerPoints[1][0]+lowerPoints[2][0])/6.;
		inner[1] = (upperPoints[0][1]+upperPoints[1][1]+upperPoints[2][1]+lowerPoints[0][1]+lowerPoints[1][1]+lowerPoints[2][1])/6.;
		inner[2] = (upperPoints[0][2]+upperPoints[1][2]+upperPoints[2][2]+lowerPoints[0][2]+lowerPoints[1][2]+lowerPoints[2][2])/6.;

		// compute the del facets planes and orient
		plane3points(mc->mc_points[0],mc->mc_points[2],mc->mc_points[4],mc->upper);
		// orient
		sign = DOT(mc->upper,inner) + mc->upper[3];
		if (sign>0) ASSIGN4(mc->upper, -mc->upper)

		plane3points(mc->mc_points[1],mc->mc_points[3],mc->mc_points[5],mc->lower);
		// orient
		sign = DOT(mc->lower,inner) + mc->lower[3];
		if (sign>0) ASSIGN4(mc->lower, -mc->lower)

		// compute the lateral planes and orient
		double *plane;
		plane = mc->planes[0];
		plane3points(mc->mc_points[0],mc->mc_points[1],mc->mc_points[2],plane);
		// orient
		sign = DOT(plane,inner) + plane[3];
		if (sign>0) ASSIGN4(plane, -plane)

		plane = mc->planes[1];
		plane3points(mc->mc_points[2],mc->mc_points[3],mc->mc_points[4],plane);
		// orient
		sign = DOT(plane,inner) + plane[3];
		if (sign>0) ASSIGN4(plane, -plane)

		plane = mc->planes[2];
		plane3points(mc->mc_points[4],mc->mc_points[5],mc->mc_points[0],plane);
		// orient
		sign = DOT(plane,inner) + plane[3];
		if (sign>0) ASSIGN4(plane, -plane)

		double o1[3], o2[3], o3[3], ttt;

		// vector parallel to voronoi edge
		SUB(o1,vor1,vor2)
		NORMALIZE(o1,ttt)

		// vector parallel to delaunay face
		SUB(o2,delFace[0],delFace[1])
		NORMALIZE(o2,ttt)

		// cross product
		CROSS(o3,o1,o2)
		NORMALIZE(o3,ttt)

		// o1,o2,o3 is an orthonormal set
		// compute equation

		for (int ii=0;ii<4;ii++)
			for (int jj=0;jj<4;jj++)
			{
				// C_[ii][jj]		=	0;
				Cinv_[ii][jj]	=	0;
				M_[ii][jj]		=	0;
				Q_[ii][jj]		=	0;
				TEMP_[ii][jj]	=	0;
				RM_[ii][jj]		=	0;
				T_[ii][jj]		=	0;
			}
		RM_[0][0] = o2[0];
		RM_[1][0] = o2[1];
		RM_[2][0] = o2[2];

		RM_[0][1] = o3[0];
		RM_[1][1] = o3[1];
		RM_[2][1] = o3[2];

		RM_[0][2] = o1[0];
		RM_[1][2] = o1[1];
		RM_[2][2] = o1[2];

		RM_[3][3] = 1.;

		T_[0][0] = 1.;
		T_[1][1] = 1.;
		T_[2][2] = 1.;
		T_[0][3] = focus[0];
		T_[1][3] = focus[1];
		T_[2][3] = focus[2];
		T_[3][3] = 1.;

		Matrix4x4MultiplyBy4x4(T_,RM_,Cinv_);

		//for (int ii=0;ii<4;ii++)
		//	for (int jj=0;jj<4;jj++)
		//		C_[ii][jj] = Cinv_[ii][jj];

		inplace_invert4x4(Cinv_);

		M_[0][0] = -1./(1.-s);
		M_[1][1] = -1./(1.-s);
		M_[2][2] = 1./s;
		M_[3][3] = -radius2;

		Matrix4x4MultiplyBy4x4(M_,Cinv_,TEMP_);

		for (int ii = 0; ii<4; ii++)
			for (int jj = ii+1; jj<4; jj++)
			{
				double temp = Cinv_[ii][jj];
				Cinv_[ii][jj] = Cinv_[jj][ii];
				Cinv_[jj][ii] = temp;
			}
		Matrix4x4MultiplyBy4x4(Cinv_,TEMP_,Q_);

		////////////////////////// store equation and simplex ////////////////////////////////////////////

		double *QQ = allocateVector<double>(10);
		mc->quadric = QQ;

		// Symmetric matrix coding
		//	Q0	Q4	Q7	Q9
		//		Q1	Q5	Q8
		//			Q2	Q6
		//				Q3

		QQ[0] = Q_[0][0];
		QQ[1] = Q_[1][1];
		QQ[2] = Q_[2][2];
		QQ[3] = Q_[3][3];
		QQ[4] = Q_[0][1];
		QQ[7] = Q_[0][2];
		QQ[9] = Q_[0][3];
		QQ[5] = Q_[1][2];
		QQ[8] = Q_[1][3];
		QQ[6] = Q_[2][3];

		////////////////////////// end store equation and simplex ////////////////////////////////////////////

		// ok this mc really exists and it is acceptable
		// mixedComplex.push_back(mc);
		++num_d2_v1_patches;
	}
	num_cells[DELAUNAY_FACET_CELL] = num_d2_v1_patches;

	*any_error = false;
}


#ifdef MULTITHREADED_SKIN_BUILDING

void SkinSurface::BuildPatches (ThreadDataWrapper *tdw, double minmax[6], bool *any_error)
{
	#ifdef CGAL_LINKED_WITH_TBB
	Rt rT;

	// Construct the locking data-structure, using the bounding-box of the points
	Rt::Lock_data_structure locking_ds(CGAL::Bbox_3(minmax[0], minmax[1], minmax[2], minmax[3], minmax[4], minmax[5]), 50);
	rT.set_lock_data_structure(&locking_ds);
	rT.insert (tdw->l.begin(), tdw->l.end());
	#else
	rT.insert (tdw->l.begin(), tdw->l.end());
	#endif
	
	#ifdef CGAL_USE_BASIC_VIEWER
	// CGAL::draw(rT);
	#endif

	assert( rT.is_valid() );
	assert( rT.dimension() == 3 );


	#if defined(NEW_ATOM_PATCHES)
	int *atomPatches = allocateVector<int>(delphi->atoms.size());
	#endif


	#if !defined(NEW_ATOM_PATCHES)
	tdw->atomPatches = allocateVector<Del0Cell*>(delphi->atoms.size());
	#endif

	BuildDelaunayTetraCells (rT, tdw->mixedComplex, tdw->num_cells);

	bool error;

	#if !defined(NEW_ATOM_PATCHES)
	BuildDelaunayPointCells (rT, tdw->mixedComplex, tdw->atomPatches, tdw->num_cells);

	BuildDelaunayEdgeCells (rT, tdw->mixedComplex, tdw->atomPatches, delphi->atoms.size(), tdw->num_cells, &error);
	#else
	BuildDelaunayPointCells (rT, tdw->mixedComplex, atomPatches, tdw->num_cells);

	BuildDelaunayEdgeCells (rT, tdw->mixedComplex, atomPatches, delphi->atoms.size(), tdw->num_cells, &error);
	#endif
	*any_error |= error;

	BuildDelaunayFacetCells (rT, tdw->mixedComplex, tdw->num_cells, &error);
	*any_error |= error;


	#if !defined(NEW_ATOM_PATCHES)
	deleteVector<Del0Cell*>(tdw->atomPatches);
	#else
	deleteVector<int>(atomPatches);
	#endif
}

#endif


bool SkinSurface::buildSkinCGAL()
{
	auto chrono_start = chrono::high_resolution_clock::now();

	#if !defined(NO_CGAL_PATCHING)
	vector<Weighted_point> l;
	#else
	// nopatch
	vector<std::pair<Weighted_point, int>> l;
	#endif

	double x, y, z, r;
	double max_x=-1e6,min_x=1e6,max_y=-1e6,min_y=1e6,max_z=-1e6,min_z=1e6;
	
	cout << endl << INFO << "Building skin surface..";


	if (numMixedCells != 0)
	{
		for (int i=0; i < numMixedCells; i++)
			delete mixedComplex[i];

		// deleteVector<MixedCell*>(mixedComplex);
		mixedComplex.clear();
	}

	int dummy_atoms = 6;

	l.reserve(delphi->atoms.size() + dummy_atoms);

	// for (int i=0; i<delphi->numAtoms; i++)
	for (int i=0; i<delphi->atoms.size(); i++)
	{
		{
			#if defined(ENABLE_BOOST_THREADS) && defined(MULTITHREADED_POCKET_LOOP)
			if (conf.parallelPocketLoop)
				boost::mutex::scoped_lock scopedLock(mutex);
			#endif
			delphi->atoms[i].pos[0] += randDisplacement*(randnum()-0.5);
			delphi->atoms[i].pos[1] += randDisplacement*(randnum()-0.5);
			delphi->atoms[i].pos[2] += randDisplacement*(randnum()-0.5);
		}
		x = delphi->atoms[i].pos[0];
		y = delphi->atoms[i].pos[1];
		z = delphi->atoms[i].pos[2];
		r = delphi->atoms[i].radius;

		#if !defined(NO_CGAL_PATCHING)
		l.push_back(Weighted_point(Point3(x,y,z), (r*r)/s, i));
		#else
		// nopatch
		l.push_back(std::make_pair(Weighted_point(Point3(x,y,z), (r*r)/s), i));
		#endif

		max_x = max(max_x,x+r);
		max_y = max(max_y,y+r);
		max_z = max(max_z,z+r);
		
		min_x = min(min_x,x-r);
		min_y = min(min_y,y-r);
		min_z = min(min_z,z-r);
	}
	double mid_x = (max_x+min_x)*0.5;
	double mid_y = (max_y+min_y)*0.5;
	double mid_z = (max_z+min_z)*0.5;

	min_x -= fabs(mid_x-min_x)*2.;
	min_y -= fabs(mid_y-min_y)*2.;
	min_z -= fabs(mid_z-min_z)*2.;

	max_x += fabs(mid_x-max_x)*2.;
	max_y += fabs(mid_y-max_y)*2.;
	max_z += fabs(mid_z-max_z)*2.;

	double dummy_atom[6][3] = {
		{min_x,mid_y,mid_z},
		{max_x,mid_y,mid_z},
		{mid_x,min_y,mid_z},
		{mid_x,max_y,mid_z},
		{mid_x,mid_y,min_z},
		{mid_x,mid_y,max_z}
	};

	double dummy = -1.0;

	// add bounding box using nil/negative weights to let easy detection of these virtual atoms
	#if !defined(NO_CGAL_PATCHING)
	for (int i=0; i < dummy_atoms; i++)
		l.push_back(Weighted_point(Point3(dummy_atom[i][0],dummy_atom[i][1],dummy_atom[i][2]),dummy,delphi->atoms.size()+i));
	#else
	for (int i=0; i < dummy_atoms; i++)
		l.push_back(std::make_pair(Weighted_point(Point3(dummy_atom[i][0],dummy_atom[i][1],dummy_atom[i][2]),dummy),delphi->atoms.size()+i));
	#endif

	cout << endl << INFO << "Regular triangulation and patch computations ...";
	cout.flush();


	int num_cells[] = {0, 0, 0, 0};


	#if !defined(MULTITHREADED_SKIN_BUILDING)
	numThreadDataWrappers = 1;
	#else
	numThreadDataWrappers = conf.numThreads;
	#endif


	#ifdef MULTITHREADED_SKIN_BUILDING

	// thread_data_wrapper = new ThreadDataWrapper[ numTasks * numThreadDataWrappers ]();
	
	// allocation of private, per-thread data
	for (int t_id=0; t_id<numTasks; t_id++)
	{
		for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
		{
			ThreadDataWrapper *tdw = &thread_data_wrapper[t_id*numThreadDataWrappers+thd_id];

			tdw->mixedComplex.reserve( MAX(10, delphi->atoms.size()*20 / (double)(numTasks*numThreadDataWrappers) );

			tdw->l.reserve(delphi->atoms.size() / (2.*numTasks*numThreadDataWrappers));

			for (int i=0; i<4; i++)
				tdw->num_cells[i] = 0;
		}
	}

	// here there is a code that subdivides the domain along the Z-axis in N=numThreads=numThreadDataWrappers
	// slabs whose thickness varies according to the population of atoms;
	// a large number of bins/slabs (numThreads * bin_factor) stores the number of atoms
	// in each bin
	double slab_max_z[numThreadDataWrappers];
	
	// Thread-specific, Z-slab dependent assignments with halo layers of 9 Angstrom
	double halo_layer_size = 9;

	int bin_factor = 1000;
	int num_bins = numThreadDataWrappers * bin_factor;
	int num_bin_atoms[num_bins];
	int tagged_thread_in_bin[num_bins];

	memset(num_bin_atoms, 0, num_bins * sizeof(int));

	double bin_slab_size = (max_z - min_z) / num_bins;
	double inv_bin_slab_size = num_bins / (max_z - min_z);

	for (int i=0; i<delphi->atoms.size(); i++)
	{
		int bin_id = (int)((delphi->atoms[i].pos[2] - min_z) * inv_bin_slab_size);
		++num_bin_atoms[bin_id];
	}
	int target_num_atoms_per_thread = (int)((double)delphi->atoms.size() / (double)numThreadDataWrappers);
	int current_atoms_sum = 0;
	int previous_atoms_sum = 0;
	int thread_id = 0;

	for (int k=0; k < num_bins; k++)
	{
		current_atoms_sum += num_bin_atoms[k];

		double upper_slab_edge = min_z + (k+1) * bin_slab_size;
		// be sure that the current slab's thickness is at least halo_layer_size

		tagged_thread_in_bin[k] = thread_id;

		if (thread_id > 0)
		{
			if (upper_slab_edge - slab_max_z[ thread_id-1 ] < halo_layer_size)
				continue;
		}
		if (current_atoms_sum - previous_atoms_sum >= target_num_atoms_per_thread)
		{
			previous_atoms_sum = current_atoms_sum;
			slab_max_z[thread_id] = upper_slab_edge;
			++thread_id;
		}
		if (thread_id == numThreadDataWrappers - 1)
		{
			slab_max_z[ thread_id ] = max_z;

			for (++k; k < num_bins; k++)
				tagged_thread_in_bin[k] = thread_id;
			break;
		}
	}
	int max_numThreadDataWrappers = numThreadDataWrappers;
	int last_k;

	for (int k=0; k < numThreadDataWrappers-1; k++)
	{
		last_k = k;

		if (slab_max_z[ k+1 ] - slab_max_z[ k ] < halo_layer_size)
		{
			// the slab is enlarged so that its border will be the grid border
			slab_max_z[ k ] = min_z;

			if (max_numThreadDataWrappers > k+1)
			{
				max_numThreadDataWrappers = k+1;
				break;
			}
		}
	}
	for (int k = last_k + 1; k < numThreadDataWrappers; k++)
		slab_max_z[ k ] = max_z;

	if (max_numThreadDataWrappers < numThreadDataWrappers)
	{
		for (int k=0; k < num_bins; k++)
		{
			if (tagged_thread_in_bin[ k ] > max_numThreadDataWrappers - 1)
				tagged_thread_in_bin[ k ] = max_numThreadDataWrappers - 1;
		}
		numThreadDataWrappers = max_numThreadDataWrappers;
	}
	if (numThreadDataWrappers > 1)
	{
		if (slab_max_z[ numThreadDataWrappers-1 ] - slab_max_z[ numThreadDataWrappers-2 ] < halo_layer_size)
		{
			cout << endl << ERR << "Too small subdomain for multithreaded cells' building!";
			return false;
		}
	}

	/*
	printf (" min_z, max_z = %f %f\n", min_z, max_z);
	for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
		printf ("thread ID, slab_max_z = %i %.14e\n", thread_id, slab_max_z[thread_id]);
	*/

	// code to check the atoms' population across bins
	int slab_atoms[numThreadDataWrappers];

	for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
		slab_atoms[thread_id] = 0;

	for (int i=0; i<delphi->atoms.size(); i++)
	{
		int bin_id = (int)((delphi->atoms[i].pos[2] - min_z) * inv_bin_slab_size);
		++slab_atoms[ tagged_thread_in_bin[bin_id] ];
	}

	for (int i=0; i<delphi->atoms.size(); i++)
	{
		x = delphi->atoms[i].pos[0];
		y = delphi->atoms[i].pos[1];
		z = delphi->atoms[i].pos[2];
		r = delphi->atoms[i].radius;

		int bin_id = (int)((z - min_z) * inv_bin_slab_size);
		int thread_slab_id = tagged_thread_in_bin[ bin_id ];
		int involved_thread[3];
		int num_involved_threads = 1;

		involved_thread[0] = thread_slab_id;

		// If the reference thread is not the last one, check if the above
		// thread is also involved because of the halo space
		if (thread_slab_id < numThreadDataWrappers - 1)
		{
			// check if the atom is in the halo space of the above thread
			if (z >= slab_max_z[ thread_slab_id ] - halo_layer_size)
				involved_thread[ num_involved_threads++ ] = thread_slab_id + 1;
		}
		// If the reference thread is not the first one, check if the below thread is also involved because of the halo space
		if (thread_slab_id > 0)
		{
			// check if the atom is in the halo space of the below thread
			if (z <= slab_max_z[ thread_slab_id-1 ] + halo_layer_size)
				involved_thread[ num_involved_threads++ ] = thread_slab_id - 1;
		}
		// Atom data is loaded by 1-3 threads
		#if !defined(NO_CGAL_PATCHING)
		for (int m = 0; m < num_involved_threads; m++)
			thread_data_wrapper[ involved_thread[m] ].l.push_back( l[i] );
		#else
		for (int m = 0; m < num_involved_threads; m++)
			thread_data_wrapper[ involved_thread[m] ].l.push_back( std::make_pair(Weighted_point(Point3(x,y,z), (r*r)/s), i) );
		#endif
	}
	for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
	{
		for (int i=delphi->atoms.size(); i<delphi->atoms.size() + dummy_atoms; i++)
		{
			#if !defined(NO_CGAL_PATCHING)
			thread_data_wrapper[ thread_id ].l.push_back( l[i] );
			#else
			thread_data_wrapper[ thread_id ].l.push_back( std::make_pair(Weighted_point(Point3(x,y,z), (r*r)/s), i) );
			#endif
		}
	}


	bool *any_error = new bool [numThreadDataWrappers];

	for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
		any_error[thread_id] = 0;
	
	double minmax[] = {min_x, min_y, min_z, max_x, max_y, max_z};
	
	#ifdef ENABLE_BOOST_THREADS
	boost::thread_group thdGroup;
	#endif
	for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
	{
		ThreadDataWrapper *tdw = &thread_data_wrapper[thread_id];

		#if defined(ENABLE_BOOST_THREADS)
		thdGroup.create_thread(boost::bind(&SkinSurface::BuildPatches, this, tdw, minmax, &any_error[thread_id]));
		#else
		BuildPatches (tdw, minmax, &any_error[thread_id]);
		#endif
	}
	#ifdef ENABLE_BOOST_THREAS
	thdGroup.join_all();
	#endif
	

	for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
	{
		thread_data_wrapper[thread_id].l.clear();
	}

	bool error = 0;
	
	for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
		error |= any_error[thread_id];
	
	delete any_error;


	if (error)
	{
		for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
		{
			ThreadDataWrapper *tdw = &thread_data_wrapper[thread_id];

			for (unsigned int i=0; i<tdw->mixedComplex.size(); i++)
				delete tdw->mixedComplex[i];
			
			tdw->mixedComplex.clear();
		}
		return false;
	}


	// assign pointers, consolidate mixed complex
	mixedComplex.reserve(delphi->atoms.size()*20);
	
	for (int thread_id = 0; thread_id < numThreadDataWrappers; thread_id++)
	{
		ThreadDataWrapper *tdw = &thread_data_wrapper[thread_id];

		for (int m = 0; m < tdw->mixedComplex.size(); m++)
		{
			MixedCell *mc = tdw->mixedComplex[m];

			double downz = INFINITY;

			if (mc->surface_type == DELAUNAY_POINT_CELL)
			{
				int index = ((Del0Cell*)mc)->id;
				
				// avoid bounding box virtual atom
				if (index >= delphi->atoms.size())
					continue;

				double radius = delphi->atoms[index].radius;
				double sphere_center_z = delphi->atoms[index].pos[2];
				
				downz = sphere_center_z - radius;
			}
			else
			{
				for (unsigned int i=0; i<mc->mc_points.size(); i++)
					downz = MIN(downz, mc->mc_points[i][2]);
			}
			// if the bottom border of the patch resides within the slab of another thread the cell is discarded ...
			double slab_border;

			// lower slab border
			if (thread_id == 0)
				slab_border = min_z;
			else
				slab_border = slab_max_z[ thread_id-1 ];
			// ... in particular,
			// check if the lower border of the patch is inside another slab/subdomain;
			// if so, the patch will be processed by the thread responsible for the neighbour slab
			if (downz < slab_border)
			{
				--tdw->num_cells[ mc->surface_type ];
				mc->surface_type == SKIP_CELL;
				continue;
			}
			// If the lower border of the patch is above the upper border of the current slab
			// then the patch will be processed by the thread responsible for the upper slab
			if (downz > slab_max_z[ thread_id ])
			{
				--tdw->num_cells[ mc->surface_type ];
				mc->surface_type = SKIP_CELL;
				continue;
			}
			mixedComplex.push_back(mc);
		}
	}


	for (int thread_id=0; thread_id<numThreadDataWrappers; thread_id++)
	{
		for (int i=0; i<4; i++)
			num_cells[i] += thread_data_wrapper[thread_id].num_cells[i];
	}

	#else // MULTITHREADED_SKIN_BUILDING

	Rt rT;


	#ifdef CGAL_LINKED_WITH_TBB
	// Construct the locking data-structure, using the bounding-box of the points
	Rt::Lock_data_structure locking_ds(CGAL::Bbox_3(min_x, min_y, min_z, max_x, max_y, max_z), 50);
	rT.set_lock_data_structure(&locking_ds);
	rT.insert (l.begin(), l.end());
	#else
	rT.insert (l.begin(), l.end());
	#endif

	#ifdef CGAL_USE_BASIC_VIEWER
	CGAL::draw(rT);
	#endif

	assert( rT.is_valid() );
	assert( rT.dimension() == 3 );

	mixedComplex.reserve(delphi->atoms.size()*20);


	#if !defined(NEW_ATOM_PATCHES)
	// Del0Cell **atomPatches = allocateVector<Del0Cell*>(delphi->numAtoms);
    Del0Cell **atomPatches = allocateVector<Del0Cell*>(delphi->atoms.size());
	#else
	int *atomPatches = allocateVector<int>(delphi->atoms.size());
	#endif


	bool any_error = 0;
	bool error;


	BuildDelaunayTetraCells (rT, mixedComplex, num_cells);

	BuildDelaunayPointCells (rT, mixedComplex, atomPatches, num_cells);

	BuildDelaunayEdgeCells (rT, mixedComplex, atomPatches, delphi->atoms.size(), num_cells, &error);
	any_error |= error;

	BuildDelaunayFacetCells (rT, mixedComplex, num_cells, &error);
	any_error |= error;

	#if defined(NEW_ATOM_PATCHES)
	deleteVector<int>(atomPatches);
	#endif

	numMixedCells = (int)mixedComplex.size();

	#endif // MULTITHREADED_SKIN_BUILDING


	type[DELAUNAY_TETRA_CELL] = num_cells[DELAUNAY_TETRA_CELL];
	type[DELAUNAY_POINT_CELL] = num_cells[DELAUNAY_POINT_CELL];
	type[DELAUNAY_EDGE_CELL]  = num_cells[DELAUNAY_EDGE_CELL];
	type[DELAUNAY_FACET_CELL] = num_cells[DELAUNAY_FACET_CELL];


	auto chrono_end = chrono::high_resolution_clock::now();
	
	chrono::duration<double> build_time = chrono_end - chrono_start;
	cout << endl << INFO << "Regular triangulation and patches' build-up computing time.. ";
	printf ("%.4e [s]", build_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	if (!conf.parallelPocketLoop)
	{
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after build-up is " << current_mem_in_MB << " MB";
	}
	#endif

	// Be sure of be out of pocket mode, otherwise p
	if (savePovRay)
	{
		ofstream of;
		cout << endl << INFO << "Saving surface in Pov-Ray in skin.pov...";
		cout.flush();
		of.open("skin.pov");
		of << "#include \"shapes.inc\" ";
		of << "\n#include \"colors.inc\" ";
		of << "\nglobal_settings {max_trace_level 3}";
		of << "\nbackground { color Black }";
		of << "\ncamera {";
		of << "\n\tlocation  <" << mid_x << ","  << max_y <<"," << mid_z <<">";
		of << "\nlook_at  <" << mid_x << "," << mid_y << "," << mid_z << "> translate-<" << mid_x << "," << mid_y << "," << mid_z << ">" << " rotate<0,0,clock> "<< "translate<" << mid_x << "," << mid_y << "," << mid_z << ">}";
		of << "\nlight_source {<" << mid_x << ","  << max_y <<"," << mid_z <<">" << " color White "<< "translate-<" << mid_x << "," << mid_y << "," << mid_z << ">" << " rotate<0,0,clock> "<< "translate<" << mid_x << "," << mid_y << "," << mid_z << ">}";

		#if defined(MULTITHREADED_POCKET_LOOP)
		if (!conf.parallelPocketLoop)
		{
			cout << endl << WARN << "Cannot write skin patches in pocket mode with pocket-level multithreading with thread-safe randnum() without mutexes";
		#endif
			for (int k=0; k < numMixedCells; k++)
			{
				MixedCell *mc = mixedComplex[k];

				if (mc->surface_type == SKIP_CELL)
					continue;

				saveSkinPatch(of,mc,k,l);
			}
		#if defined(MULTITHREADED_POCKET_LOOP)
		}
		#endif

		cout << "ok!";
		of.close();
	}

	l.clear();

	printSummary();

	return true;
}

#endif // ENABLE_CGAL


#ifdef ENABLE_CGAL

#if !defined(NO_CGAL_PATCHING)
void SkinSurface::saveSkinPatch(ofstream &of, MixedCell *mc,int index,vector<Weighted_point> &la)
#else
// nopatch
void SkinSurface::saveSkinPatch(ofstream &of, MixedCell *mc,int index,vector<std::pair<Weighted_point, int>> &la)
#endif
{
	char clipplanes[1000],buff2[1000];

	if (mc->surface_type == DELAUNAY_EDGE_CELL)
	{	
		Del1Cell *ec = (Del1Cell*)mc;

		of << "\n// ---------------------------------------------//";
		of << "\n// ------------- Delunay edge cell -------------//";
		of << "\n// ---------------------------------------------//";
		of << "\n// atoms " << ec->ids[0] << "," << ec->ids[1];

		sprintf(clipplanes,"Clipping_Mesh%d", index);
		of << "\n#declare " << clipplanes << "=\n mesh{";

		// this is slow but allows Pov-Ray to work well
		Polyhedron poly;
		// little noise to avoid problems and to let povray work
		vector<Point3> local;

		for (unsigned int i=0; i<ec->mc_points.size(); i++)
		{
			double *p = ec->mc_points[i];
			local.push_back(Point3(p[0]+1e-3*(randnum()-0.5),p[1]+1e-3*(randnum()-0.5),p[2]+1e-3*(randnum()-0.5)));
		}

		// compute convex hull to recover mesh triangles
		CGAL::convex_hull_3(local.begin(), local.end(), poly);	

		// save the mesh 
		for (Polyhedron::Facet_iterator fit = poly.facets_begin(); fit != poly.facets_end(); fit++)
		{	
			of << "\n triangle{\n";
			int g = 0;
			for (HF_circulator h = fit->facet_begin();g<3;h++)
			{
				of << "<" << h->vertex()->point().x() << "," << h->vertex()->point().y() << "," << h->vertex()->point().z() << ">";
				if (g < 2) of << ",";
				g++;
			}
			of << "}";
		}

		// generate a reference inner point
		double inner[3],t1[3],t2[3];
		
		t1[0] = (ec->mc_points[0][0]+ec->mc_points[2][0]+ec->mc_points[4][0])/3.;
		t1[1] = (ec->mc_points[0][1]+ec->mc_points[2][1]+ec->mc_points[4][1])/3.;
		t1[2] = (ec->mc_points[0][2]+ec->mc_points[2][2]+ec->mc_points[4][2])/3.;

		t2[0] = (ec->mc_points[1][0]+ec->mc_points[3][0]+ec->mc_points[5][0])/3.;
		t2[1] = (ec->mc_points[1][1]+ec->mc_points[3][1]+ec->mc_points[5][1])/3.;
		t2[2] = (ec->mc_points[1][2]+ec->mc_points[3][2]+ec->mc_points[5][2])/3.;

		inner[0] = (t1[0]+t2[0])*0.5;
		inner[1] = (t1[1]+t2[1])*0.5;
		inner[2] = (t1[2]+t2[2])*0.5;

		of  << "\n inside_vector<" << inner[0] << "," << inner[1] << "," << inner[2] << ">";
		of << "}";

			
		// this is formally correct but Pov-Ray seems not to be able to manage this code for 
		// more than very small molecules
		/*
		// add clipping planes
		of << "\n\n #declare " << clipplanes << " =" ;
		sprintf(buff3,"\n\n // clipping planes \n\n intersection { \n ");
		of << buff3;
		
		for (unsigned int l=0; l<ec->planes.size(); l++)
		{
			sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",ec->planes[l][0],ec->planes[l][1],ec->planes[l][2],-(ec->planes[l][3])/(sqrt((DOT(ec->planes[l],ec->planes[l])))));
			of << buff3;
		}

		sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",ec->lower[0],ec->lower[1],ec->lower[2],-(ec->lower[3])/(sqrt((DOT(ec->lower,ec->lower)))));
		of << buff3;

		sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",ec->upper[0],ec->upper[1],ec->upper[2],-(ec->upper[3])/(sqrt((DOT(ec->upper,ec->upper)))));
		of << buff3;
		of << "\n}";
		*/
		
		double A,B,C,D,E,F,G,H,I,J;
		A = mc->quadric[0];
		B = mc->quadric[1];
		C = mc->quadric[2];
		J = mc->quadric[3];
		D = 2*mc->quadric[4];
		E = 2*mc->quadric[7];
		G = 2*mc->quadric[9];
		F = 2*mc->quadric[5];
		H = 2*mc->quadric[8];
		I = 2*mc->quadric[6];

		sprintf(buff2, "\n\n quadric { \n <%f,%f,%f>,<%f,%f,%f>,<%f,%f,%f>,%f", A,B,C,D,E,F,G,H,I,J);
		of << buff2;
		of << "\n pigment{color Yellow}";
		of << "\n bounded_by { " << clipplanes << " } \n clipped_by{ bounded_by}}"; 
		of << "\n";
	}
	else if (mc->surface_type == DELAUNAY_POINT_CELL)
	{
		Del0Cell *pc = (Del0Cell*)mc;

		of << "\n// ---------------------------------------------//";
		of << "\n// ------------- Delunay point cell ------------//";
		of << "\n// ---------------------------------------------//";
		of << "\n// atom " << pc->id;

		#if !defined(NO_CGAL_PATCHING)
		sprintf(buff2, "\n\n sphere { \n <%f,%f,%f>,%f", la[pc->id].x(),la[pc->id].y(),la[pc->id].z(),sqrt(la[pc->id].weight()*s));
		#else
		sprintf(buff2, "\n\n sphere { \n <%f,%f,%f>,%f", la[pc->id].first.x(),la[pc->id].first.y(),la[pc->id].first.z(),sqrt(la[pc->id].first.weight()*s));
		#endif

		of << buff2;
		of << "\n pigment{color Green}}";
	}
	else if (mc->surface_type == DELAUNAY_FACET_CELL)
	{
		Del2Cell *fc = (Del2Cell*)mc;

		of << "\n// ----------------------------------------------//";
		of << "\n// ------------- Delunay facet cell -------------//";
		of << "\n// ----------------------------------------------//";
		of << "\n// atoms " << fc->ids[0] << "," << fc->ids[1] << "," << fc->ids[2];

		sprintf(clipplanes, "Clipping_Mesh%d", index);
		of << "\n#declare " << clipplanes << "=\n mesh{";

		// this is slow but allows Pov-Ray to work well
		Polyhedron poly;
		// little noise to avoid problems
		vector<Point3> local;

		for (unsigned int i=0; i<fc->mc_points.size(); i++)
		{
			double *p = fc->mc_points[i];
			local.push_back(Point3(p[0]+1e-3*(randnum()-0.5),p[1]+1e-3*(randnum()-0.5),p[2]+1e-3*(randnum()-0.5)));
		}

		// compute convex hull to recover mesh triangles
		CGAL::convex_hull_3(local.begin(), local.end(), poly);

		// save the mesh 
		for (Polyhedron::Facet_iterator fit = poly.facets_begin(); fit != poly.facets_end(); fit++)
		{
			of << "\n triangle{\n";
			int g = 0;
			for (HF_circulator h = fit->facet_begin(); g<3; h++)
			{
				of << "<" << h->vertex()->point().x() << "," << h->vertex()->point().y() << "," << h->vertex()->point().z() << ">";
				if (g<2) of << ",";
				g++;
			}
			of << "}";
		}
		// generate a reference inner point
		double inner[3], t1[3], t2[3];
		
		t1[0] = (fc->mc_points[0][0]+fc->mc_points[2][0]+fc->mc_points[4][0])/3.;
		t1[1] = (fc->mc_points[0][1]+fc->mc_points[2][1]+fc->mc_points[4][1])/3.;
		t1[2] = (fc->mc_points[0][2]+fc->mc_points[2][2]+fc->mc_points[4][2])/3.;

		t2[0] = (fc->mc_points[1][0]+fc->mc_points[3][0]+fc->mc_points[5][0])/3.;
		t2[1] = (fc->mc_points[1][1]+fc->mc_points[3][1]+fc->mc_points[5][1])/3.;
		t2[2] = (fc->mc_points[1][2]+fc->mc_points[3][2]+fc->mc_points[5][2])/3.;

		inner[0] = (t1[0]+t2[0])*0.5;
		inner[1] = (t1[1]+t2[1])*0.5;
		inner[2] = (t1[2]+t2[2])*0.5;

		of  << "\n inside_vector<" << inner[0] << "," << inner[1] << "," << inner[2] << ">";
		of << "}";

		// this is formally correct but Pov-Ray seems not to be able to manage this code for 
		// more than very small molecules
		/*
		// add clipping planes
		of << "\n\n #declare " << clipplanes << " =" ;
		sprintf(buff3,"\n\n // clipping planes \n\n intersection { \n ");
		of << buff3;
		
		for (int l=0;l<3;l++)
		{
			sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",fc->planes[l][0],fc->planes[l][1],fc->planes[l][2],-(fc->planes[l][3])/(sqrt((DOT(fc->planes[l],fc->planes[l])))));
			of << buff3;
		}
		
		sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",fc->lower[0],fc->lower[1],fc->lower[2],-(fc->lower[3])/(sqrt((DOT(fc->lower,fc->lower)))));
		of << buff3;

		sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",fc->upper[0],fc->upper[1],fc->upper[2],-(fc->upper[3])/(sqrt((DOT(fc->upper,fc->upper)))));
		of << buff3;
		of << "\n}";
		*/

		double A,B,C,D,E,F,G,H,I,J;
		A = mc->quadric[0];
		B = mc->quadric[1];
		C = mc->quadric[2];
		J = mc->quadric[3];
		D = 2*mc->quadric[4];
		E = 2*mc->quadric[7];
		G = 2*mc->quadric[9];
		F = 2*mc->quadric[5];
		H = 2*mc->quadric[8];
		I = 2*mc->quadric[6];

		sprintf(buff2,"\n\n quadric { \n <%f,%f,%f>,<%f,%f,%f>,<%f,%f,%f>,%f",A,B,C,D,E,F,G,H,I,J);
		of << buff2;
		of << "\n pigment{color Magenta}";
		of << "\n bounded_by { " << clipplanes << " } \n clipped_by{ bounded_by}}"; 
		of << "\n";

	}
	else if (mc->surface_type == DELAUNAY_TETRA_CELL)
	{
		Del3Cell *tc = (Del3Cell*)mc;

		of << "\n// ----------------------------------------------------//";
		of << "\n// ------------- Delunay tethraedron cell -------------//";
		of << "\n// ----------------------------------------------------//";
		of << "\n// atoms " << tc->ids[0] << "," << tc->ids[1] << "," << tc->ids[2] << "," << tc->ids[3];
		
		sprintf(clipplanes,"Clipping_Mesh%d",index);
		of << "\n#declare " << clipplanes << "=\n mesh{";

		// this is slow but allows Pov-Ray to work well
		Polyhedron poly;
		// little noise to avoid problems
		vector<Point3> local;
		double pp[3];
		for (unsigned int i=0; i<4; i++)
		{
			pp[0] = tc->reduced[i][0];
			pp[1] = tc->reduced[i][1];
			pp[2] = tc->reduced[i][2];
			local.push_back(Point3(pp[0]+1e-3*(randnum()-0.5),pp[1]+1e-3*(randnum()-0.5),pp[2]+1e-3*(randnum()-0.5)));
		}

		// compute convex hull to recover mesh triangles
		CGAL::convex_hull_3(local.begin(), local.end(), poly);

		// save the mesh 
		for (Polyhedron::Facet_iterator fit = poly.facets_begin(); fit != poly.facets_end(); fit++)
		{
			of << "\n triangle{\n";
			int g = 0;
			for (HF_circulator h = fit->facet_begin();g<3;h++)
			{
				of << "<" << h->vertex()->point().x() << "," << h->vertex()->point().y() << "," << h->vertex()->point().z() << ">";
				if (g<2) of << ",";
				g++;
			}
			of << "}";
		}

		of  << "\n inside_vector<" << tc->vor[0] << "," << tc->vor[1] << "," << tc->vor[2] << ">";
		of << "}";

		// this is formally correct but Pov-Ray seems not to be able to manage this code for 
		// more than very small molecules
		/*
		sprintf(clipplanes,"Clipping_Planes%d",index);
		of << "\n\n #declare " << clipplanes << " =" ;
		of << "\n intersection {";
		for (int l=0;l<4;l++)
		{
			sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",tc->planes[l][0],tc->planes[l][1],tc->planes[l][2],-(tc->planes[l][3])/(sqrt((DOT(tc->planes[l],tc->planes[l])))));
			of << buff3;
		}
		of << "\n}";
		*/
		double radius2,v[3],norm2,w,radius;
		#if !defined(NO_CGAL_PATCHING)
		v[0] = la[tc->ids[0]].x() - tc->vor[0];
		v[1] = la[tc->ids[0]].y() - tc->vor[1];
		v[2] = la[tc->ids[0]].z() - tc->vor[2];
		#else
		// nopatch
		v[0] = la[tc->ids[0]].first.x() - tc->vor[0];
		v[1] = la[tc->ids[0]].first.y() - tc->vor[1];
		v[2] = la[tc->ids[0]].first.z() - tc->vor[2];
		#endif

		norm2 = DOT(v,v);
		#if !defined(NO_CGAL_PATCHING)
		w = la[tc->ids[0]].weight();
		#else
		w = la[tc->ids[0]].first.weight();
		#endif

		radius2 = -(w-norm2);
		radius = sqrt(radius2*(1-s));
		sprintf(buff2,"\n\n sphere { \n <%f,%f,%f>,%f",tc->vor[0],tc->vor[1],tc->vor[2],radius);
		of << buff2;
		of << "\n pigment{color Red}";
		of << "\n bounded_by { " << clipplanes << " } \n clipped_by{bounded_by}}";
	}
}

#endif // ENABLE_CGAL


void SkinSurface::preProcessPanel()
{
	if (numMixedCells == 0)
	{
		cout << endl << WARN << "Cannot get surface without a computed mixed complex!";
		return;
	}

	int64_t igrid = delphi->nx;
	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption

	// auxiliary grid is consistent with delphi grid in order to speed-up tracing
	int64_t gridMul = 1;
	while (igrid > AUX_GRID_DIM_SKIN_2D)
	{
		gridMul += 2;
		int64_t digrid = delphi->nx;
		while (1)
		{
			// get nearest odd multiple
			int64_t fixedPoint = (digrid+(gridMul-1))/gridMul;
			igrid = fixedPoint*gridMul;
			if (igrid%2 == 0)
				digrid=igrid+1;
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
	
	xmin_2d = delphi->baricenter[0]-(igrid-1)*0.5*side_2d;
	ymin_2d = delphi->baricenter[1]-(igrid-1)*0.5*side_2d;
	zmin_2d = delphi->baricenter[2]-(igrid-1)*0.5*side_2d;

	xmax_2d = delphi->baricenter[0]+(igrid-1)*0.5*side_2d;
	ymax_2d = delphi->baricenter[1]+(igrid-1)*0.5*side_2d;
	zmax_2d = delphi->baricenter[2]+(igrid-1)*0.5*side_2d;
	
	nx_2d = igrid;
	ny_2d = igrid;
	nz_2d = igrid;

	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridMixedCellMap2D != NULL)
		deleteVector<int>(gridMixedCellMap2D);

	if (ind_2d != NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,last_cols_ind,ind_2d);
	#else */
	if (gridMixedCellMap2D != NULL)
	{
		for (int64_t i=0; i < last_rows_ind*last_cols_ind; i++)
		{
			gridMixedCellMap2D[i].clear();
		}
		delete[] gridMixedCellMap2D;
		gridMixedCellMap2D = NULL;
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
	gridMixedCellMap2D = allocateVector<int>(last_rows_ind*last_cols_ind*MAX_MIXED_CELLS_2D);

	// look up, it has been deallocated
	ind_2d = allocateMatrix2D<unsigned int>(last_rows_ind,last_cols_ind);
	for (int i=0;i<last_rows_ind;i++)
		for (int j=0;j<last_cols_ind;j++)
			ind_2d[i][j] = 0;

	int max_t = 0;
	#else */
	#if !defined(MULTITHREADED_SKIN_BUILDING)
	gridMixedCellMap2D = new vector<int> [ last_rows_ind*last_cols_ind ];
	#else
	gridMixedCellMap2D = new vector<pair<int,int>> [ last_rows_ind*last_cols_ind ];
	#endif
	// #endif

	int max_num_cells_per_pixel = 0;

	// build a bounding box for cell and map it to the auxiliary grid
	#if !defined(MULTITHREADED_SKIN_BUILDING)
	for (unsigned int it=0; it<getNumPatches(); it++)
	{
		{
			MixedCell *mc = mixedComplex[it];

			if (mc->surface_type == SKIP_CELL)
				continue;
	#else
	for (int thd_id=0; thd_id<numThreadDataWrappers; thd_id++)
	{
		for (unsigned int it=0; it < thread_data_wrapper[thd_id].mixedComplex.size(); it++)
		{
			MixedCell *mc = thread_data_wrapper[thd_id].mixedComplex[it];

			if (mc->surface_type == SKIP_CELL)
				continue;
	#endif
			/*
			if (ll.size() == 0)
			{
				cout << endl << ERR << "Empty mixed cell!";
				exit(-1);
			}
			*/

			// compute the bounding box of the object
			double downx = INFINITY;
			double downy = INFINITY;
			double downz = INFINITY;

			double upx = -INFINITY;
			double upy = -INFINITY;
			double upz = -INFINITY;

			// use the cube of the atom as bounding box object for voronoi cells
			if (mc->surface_type == DELAUNAY_POINT_CELL)
			{
				int index = ((Del0Cell*)mc)->id;

				// avoid mapping a non atom (bounding box virtual atom)
				if (index >= delphi->atoms.size())
				{
					continue;
				}
				// double radius = delphi->atoms[index]->radius;
				// double *sphere_center = delphi->atoms[index]->pos;
				double radius = delphi->atoms[index].radius;
				double *sphere_center = delphi->atoms[index].pos;

				downx = sphere_center[0]-radius;
				downy = sphere_center[1]-radius;
				downz = sphere_center[2]-radius;

				upx = sphere_center[0]+radius;
				upy = sphere_center[1]+radius;
				upz = sphere_center[2]+radius;
			}
			else
			{
				// mixed cell points
				vector<double*> &ll = mc->mc_points;

				for (unsigned int i=0; i<ll.size(); i++)
				{
					downx = MIN(downx,ll[i][0]);
					downy = MIN(downy,ll[i][1]);
					downz = MIN(downz,ll[i][2]);

					upx = MAX(upx,ll[i][0]);
					upy = MAX(upy,ll[i][1]);
					upz = MAX(upz,ll[i][2]);
				}
			}
			// Determine which are the grid cells that are
			// spanned by the bounding box of the object
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

						if (ind_2d[iz][iy] >= MAX_MIXED_CELLS_2D)
						{
							cout << endl << ERR << "Number of mixed cells is superior to maximum allowed, please increase Max_skin_patches_per_auxiliary_grid_2d_cell";
							exit(-1);
						}
						GRID_MIXED_CELL_MAP_2D(iy,iz,ind_2d[iz][iy],ny_2d,nz_2d) = it;
						ind_2d[iz][iy]++;
						#else */
						max_num_cells_per_pixel = max(max_num_cells_per_pixel, (int)gridMixedCellMap2D[ iz*ny_2d + iy ].size());
						#if !defined(MULTITHREADED_SKIN_BUILDING)
						gridMixedCellMap2D[ iz*ny_2d + iy ].push_back(it);
						#else
						gridMixedCellMap2D[ iz*ny_2d + iy ].push_back( pair<int,int>(thd_id, it) );
						#endif
						// #endif
					}
				}
			}
			else if (panel == 1)
			{
				// plane XY
				for (int64_t iy=iy_start; iy <= iy_end; iy++)
				{
					for (int64_t ix=ix_start; ix <= ix_end; ix++)
					{
						/* #if !defined(OPTIMIZE_GRIDS)
						if (ind_2d[iy][ix] > max_t)
							max_t = ind_2d[iy][ix];

						if (ind_2d[iy][ix] >= MAX_MIXED_CELLS_2D)
						{
							cout << endl << ERR << "Number of mixed cells is superior to maximum allowed, please increase Max_skin_patches_per_auxiliary_grid_2d_cell";
							exit(-1);
						}
						GRID_MIXED_CELL_MAP_2D(ix,iy,ind_2d[iy][ix],nx_2d,ny_2d) = it;
						ind_2d[iy][ix]++;
						#else */
						max_num_cells_per_pixel = max(max_num_cells_per_pixel, (int)gridMixedCellMap2D[ iy*nx_2d + ix ].size());
						#if !defined(MULTITHREADED_SKIN_BUILDING)
						gridMixedCellMap2D[ iy*nx_2d + ix ].push_back(it);
						#else
						gridMixedCellMap2D[ iy*nx_2d + ix ].push_back( pair<int,int>(thd_id, it) );
						#endif
						// #endif
					}
				}
			}
			else
			{
				// plane XZ
				for (int64_t iz=iz_start; iz <= iz_end; iz++)
				{
					for (int64_t ix=ix_start; ix <= ix_end; ix++)
					{
						/* #if !defined(OPTIMIZE_GRIDS)
						if (ind_2d[iz][ix] > max_t)
							max_t = ind_2d[iz][ix];

						if (ind_2d[iz][ix] >= MAX_MIXED_CELLS_2D)
						{
							cout << endl << ERR << "Number of mixed cells is superior to maximum allowed, please increase Max_skin_patches_per_auxiliary_grid_2d_cell";
							exit(-1);
						}
						GRID_MIXED_CELL_MAP_2D(ix,iz,ind_2d[iz][ix],nx_2d,nz_2d) = it;
						ind_2d[iz][ix]++;
						#else */
						max_num_cells_per_pixel = max(max_num_cells_per_pixel, (int)gridMixedCellMap2D[ iz*nx_2d + ix ].size());
						#if !defined(MULTITHREADED_SKIN_BUILDING)
						gridMixedCellMap2D[ iz*nx_2d + ix ].push_back(it);
						#else
						gridMixedCellMap2D[ iz*nx_2d + ix ].push_back( pair<int,int>(thd_id, it) );
						#endif
						// #endif
					}
				}
			}
		}
	}
	cout << endl << INFO << "Max number of Skin cells per panel pixel: " << max_num_cells_per_pixel;
}


bool SkinSurface::buildAuxiliaryGrid()
{
	if (getNumPatches() == 0)
	{
		cout << endl << WARN << "Cannot get surface without a computed mixed complex!";
		return false;
	}

	// Perform pre-processing to speed-up intersections and projections.
	// Compute bounding box for each mixed cell in the surface, and map each
	// bounding box to the proper grid point.
	// This routine uses an auxiliary grid. Delphi grid is not directly used
	// because it may be too memory consuming. To speed up computations
	// a maximal number of mixed cells is allowed in each auxiliary grid cell
	// the macro is MAX_MIXED_CELLS in SkinSurface.h.
	// The macro AUX_GRID_DIM_SKIN sets the maximally allowed grid size.
	int64_t igrid = delphi->nx;
	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption

	/*
	double layer = 0.0;

	if (igrid > AUX_GRID_DIM_SKIN)
		igrid = AUX_GRID_DIM_SKIN;

	scale = ((double)(igrid-1))*delphi->perfill/(100.*(delphi->rmaxdim+layer));*/

	// auxiliary grid is consistent with delphi grid in order to speed-up tracing
	int64_t gridMul = 1;
	while (igrid > AUX_GRID_DIM_SKIN)
	{
		gridMul += 2;
		int64_t digrid = delphi->nx;
		while (1)
		{
			// get nearest odd multiple
			int64_t fixedPoint = (digrid+(gridMul-1))/gridMul;
			igrid = fixedPoint*gridMul;
			if (igrid%2 == 0)
				digrid=igrid+1;
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

	cout << endl << INFO << "Auxiliary grid is " << igrid;
	
	xmin = delphi->baricenter[0]-(igrid-1)*0.5*side;
	ymin = delphi->baricenter[1]-(igrid-1)*0.5*side;
	zmin = delphi->baricenter[2]-(igrid-1)*0.5*side;

	xmax = delphi->baricenter[0]+(igrid-1)*0.5*side;
	ymax = delphi->baricenter[1]+(igrid-1)*0.5*side;
	zmax = delphi->baricenter[2]+(igrid-1)*0.5*side;

	////////////// Allocate memory for the grid, planes and maps //////////////////////

	if (x != NULL)
		deleteVector<double>(x);
	if (y != NULL)
		deleteVector<double>(y);
	if (z != NULL)
		deleteVector<double>(z);

	if (ind != NULL)
		deleteMatrix3D<int>(nx,ny,nz,ind);

	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);

	nx = igrid;
	ny = igrid;
	nz = igrid;

	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridMixedCellMap != NULL)
		deleteVector<int>(gridMixedCellMap);

	cout << endl << INFO << "Allocating " << (nx*ny*nz*MAX_MIXED_CELLS)*sizeof(int)/1024.0/1024.0 << " MB" << " for the auxiliary grid...";
	
	gridMixedCellMap = allocateVector<int>(nx*ny*nz*MAX_MIXED_CELLS);

	ind = allocateMatrix3D<int>(nx,ny,nz);
	for (int k=0; k<nz; k++)
		for (int j=0; j<ny; j++)
			for (int i=0; i<nx; i++)
				ind[k][j][i] = 0;
	#else */
	if (gridMixedCellMap != NULL)
	{
		for (int64_t i=0; i < nx*ny*nz; i++)
		{
			gridMixedCellMap[i].clear();
		}
		delete[] gridMixedCellMap;
		gridMixedCellMap = NULL;
	}
	#if !defined(MULTITHREADED_SKIN_BUILDING)
	gridMixedCellMap = new vector<int> [ nx*ny*nz ];
	#else
	gridMixedCellMap = new vector<pair<int,int>> [ nx*ny*nz ];
	#endif
	// #endif
	
	for (int i=0;i<nx;i++)
		x[i] = xmin + i*side;
	
	for (int i=0;i<ny;i++)
		y[i] = ymin + i*side;

	for (int i=0;i<nz;i++)
		z[i] = zmin + i*side;
	
	cout << "ok!";
	//////////////////////////////////////////////////////////////////////////

	// build a bounding box for each mixed cell and map it to
	// the auxiliary grid
	int max_t = 0;

	cout << endl << INFO << "Mapping auxiliary grid...";
	
	#if !defined(MULTITHREADED_SKIN_BUILDING)
	for (int it=0; it<numMixedCells; it++)
	{
		{
			MixedCell *mc = mixedComplex[it];

			if (mc->surface_type == SKIP_CELL)
				continue;
	#else
	for (int thd_id=0; thd_id<numThreadDataWrappers; thd_id++)
	{
		for (unsigned int it=0; it < thread_data_wrapper[thd_id].mixedComplex.size(); it++)
		{
			MixedCell *mc = thread_data_wrapper[thd_id].mixedComplex[it];

			if (mc->surface_type == SKIP_CELL)
				continue;
	#endif
			/*
			if (mc->mc_points.size() == 0)
			{
				cout << endl << ERR << "Empty mixed cell!";
				exit(-1);
			}
			*/

			// compute the bounding box of the object
			double downx = INFINITY;
			double downy = INFINITY;
			double downz = INFINITY;

			double upx = -INFINITY;
			double upy = -INFINITY;
			double upz = -INFINITY;

			// use the cube of the atom as bounding box object for voronoi cells
			if (mc->surface_type == DELAUNAY_POINT_CELL)
			{
				int index = ((Del0Cell*)mc)->id;

				// avoid mapping a non atom (bounding box virtual atom)
				// if (index >= delphi->numAtoms)
				if (index >= delphi->atoms.size())
					continue;

				// double radius = delphi->atoms[index]->radius;
				// double *sphere_center = delphi->atoms[index]->pos;
				double radius = delphi->atoms[index].radius;
				double *sphere_center = delphi->atoms[index].pos;

				downx = sphere_center[0]-radius;
				downy = sphere_center[1]-radius;
				downz = sphere_center[2]-radius;

				upx = sphere_center[0]+radius;
				upy = sphere_center[1]+radius;
				upz = sphere_center[2]+radius;
			}
			else
			{
				// mixed cell points
				vector<double*> &ll = mc->mc_points;

				for (unsigned int i=0; i<ll.size(); i++)
				{
					downx = MIN(downx,ll[i][0]);
					downy = MIN(downy,ll[i][1]);
					downz = MIN(downz,ll[i][2]);

					upx = MAX(upx,ll[i][0]);
					upy = MAX(upy,ll[i][1]);
					upz = MAX(upz,ll[i][2]);
				}
			}
			// resolve which are the grid cubes that
			// are cut by the bounding box object. These
			// cubes see the mixed cell bounding box
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

			for (int64_t iz=iz_start; iz<=iz_end; iz++)
			{
				for (int64_t iy=iy_start; iy<=iy_end; iy++)
				{
					for (int64_t ix=ix_start; ix<=ix_end; ix++)
					{
						/* #if !defined(OPTIMIZE_GRIDS)
						if (ind[iz][iy][ix] > max_t)
							max_t = ind[iz][iy][ix];

						if (ind[iz][iy][ix] >= MAX_MIXED_CELLS)
						{
							cout << endl << ERR << "Number of mixed cells is superior to maximum allowed, please increase Max_skin_patches_per_auxiliary_grid_cell";
							exit(-1);
						}
						GRID_MIXED_CELL_MAP(ix,iy,iz,ind[iz][iy][ix],nx,ny,nz) = it;
						ind[iz][iy][ix]++;
						#else */
						#if !defined(MULTITHREADED_SKIN_BUILDING)
						gridMixedCellMap[ iz*ny*nx + iy*nx + ix ].push_back(it);
						#else
						gridMixedCellMap[ iz*ny*nx + iy*nx + ix ].push_back( pair<int,int>(thd_id, it) );
						#endif
						// #endif
					}
				}
			}
		}
	}

	cout << "ok!";
	cout << endl << INFO << "Max mixed cells per auxiliary cell -> " << max_t;

	return true;
}


#ifdef POINT_METHOD_2
inline double SkinSurface::halfArea (double x1, double y1, double x2, double y2, double x3, double y3)
{
   return fabs(x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2));
}
#endif


#if !defined(SINGLE_PASS_RT)

#if !defined(MINIMIZE_MEMORY)

void SkinSurface::getPatchPreIntersectionData (int64_t nxyz[3], int panels[2],
											   int thread_id, int potentialIntersections[])
{
	if (getNumPatches() == 0)
	{
		return;
	}

	int64_t NX = nxyz[0];
	int64_t NY = nxyz[1];
	int64_t NZ = nxyz[2];
	int64_t N_MAX = MAX(NX, MAX(NY, NZ));

	int numPanels = panels[0] + panels[1];


	// assumption: max number of intersections per patch = 2;
	int max_intersections_per_ray = 2;

	potentialIntersections[thread_id] = 0;

	// Determine the number of potential intersections per object for allocation purposes
	for (int it = thread_id; it < numMixedCells; it += conf.numThreads)
	{
		MixedCell *mc = mixedComplex[it];

		if (mc->surface_type == SKIP_CELL)
			continue;

		#ifdef RAY_VS_CELL_TESTS_CULLING
		vector<double*> &ll = mc->mc_points;

		int projected_cell_vertices_coords[3][ll.size()];
		#endif

		// compute the bounding box of the object
		double downx = INFINITY;
		double downy = INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = -INFINITY;
		double upz = -INFINITY;

		double *sphere_center;
		double radius;

		if (mc->surface_type == DELAUNAY_POINT_CELL)
		{
			int index = ((Del0Cell*)mc)->id;

			// avoid bounding box virtual atom
			if (index >= delphi->atoms.size())
			{
				continue;
			}
			sphere_center = delphi->atoms[index].pos;
			radius = delphi->atoms[index].radius;

			downx = sphere_center[0]-radius;
			downy = sphere_center[1]-radius;
			downz = sphere_center[2]-radius;

			upx = sphere_center[0]+radius;
			upy = sphere_center[1]+radius;
			upz = sphere_center[2]+radius;
		}
		else
		{
			// mixed cell points
			vector<double*> &ll = mc->mc_points;

			for (unsigned int i=0; i<ll.size(); i++)
			{
				downx = MIN(downx,ll[i][0]);
				downy = MIN(downy,ll[i][1]);
				downz = MIN(downz,ll[i][2]);

				upx = MAX(upx,ll[i][0]);
				upy = MAX(upy,ll[i][1]);
				upz = MAX(upz,ll[i][2]);
			}
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

			#ifdef RAY_VS_CELL_TESTS_CULLING
			if (mc->surface_type == DELAUNAY_FACET_CELL ||
				mc->surface_type == DELAUNAY_EDGE_CELL)
			{
				for (unsigned int i=0; i<ll.size(); i++)
				{
					// storing the floating point projections of the vertices of the current cell
					projected_cell_vertices_coords[0][i] = (int)rintp((ll[i][0] - (delta + delphi->xmin))*delphi->scale);
					projected_cell_vertices_coords[1][i] = (int)rintp((ll[i][1] - (delta + delphi->ymin))*delphi->scale);
					projected_cell_vertices_coords[2][i] = (int)rintp((ll[i][2] - (delta + delphi->zmin))*delphi->scale);
				}
			}
			#endif
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

				#if !defined(RAY_VS_CELL_TESTS_CULLING)
				#if !defined(EQ_CULLING)
				int64_t num_pixels = pixels[last__dim] * pixels[first_dim];
				#if !defined(MULTITHREADING)
				for (int64_t n=i_start[last__dim]; n<=i_end[last__dim]; n++)
				{
					for (int64_t m=i_start[first_dim]; m<=i_end[first_dim]; m++)
					{
						int64_t id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);
						num_pixel_intersections[id] += max_intersections_per_ray;
					}
				}
				#else // MULTITHREADING
				num_patch_intersections[ patch_id ] += max_intersections_per_ray * num_pixels;
				#endif // MULTITHREADING

				potentialIntersections[thread_id] += max_intersections_per_ray * num_pixels;
				#else // EQ_CULLING

				if (mc->surface_type != DELAUNAY_POINT_CELL)
				{
					int64_t num_pixels = pixels[last__dim] * pixels[first_dim];
					#if !defined(MULTITHREADING)
					for (int64_t n=i_start[last__dim]; n<=i_end[last__dim]; n++)
					{
						for (int64_t m=i_start[first_dim]; m<=i_end[first_dim]; m++)
						{
							int64_t id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);
							num_pixel_intersections[id] += max_intersections_per_ray;
						}
					}
					#else // MULTITHREADING
					num_patch_intersections[ patch_id ] += max_intersections_per_ray * num_pixels;
					#endif // MULTITHREADING

					potentialIntersections[thread_id] += max_intersections_per_ray * num_pixels;
				}
				else
				{
					for (int64_t rectangle_pixel = 0; rectangle_pixel < pixels[first_dim]*pixels[last__dim]; rectangle_pixel++)
					{
						int64_t n = rectangle_pixel / pixels[first_dim];
						int64_t m = rectangle_pixel - n*pixels[first_dim];

						n += i_start[last__dim];
						m += i_start[first_dim];

						double orig[2];

						if (panel == 0) {
							orig[0] = delphi->y[m] + delta;
							orig[1] = delphi->z[n] + delta;
						}
						else if (panel == 1) {
							orig[0] = delphi->x[m] + delta;
							orig[1] = delphi->y[n] + delta;
						}
						else {
							orig[0] = delphi->x[m] + delta;
							orig[1] = delphi->z[n] + delta;
						}

						double dx = sphere_center[first_dim] - orig[0];
						double dy = sphere_center[last__dim] - orig[1];

						if (dx*dx + dy*dy > radius*radius)
							continue;

						#if !defined(MULTITHREADING)
						int64_t id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);
						num_pixel_intersections[id] += max_intersections_per_ray;
						#else
						num_patch_intersections[ patch_id ] += max_intersections_per_ray;
						#endif

						potentialIntersections[thread_id] += max_intersections_per_ray;
					}
				}
				#endif // EQ_CULLING
				#else // RAY_VS_CELL_TESTS_CULLING
				int num_tagged_pixels = 0;

				int v_start[2], v_end[2];
				int proj_v_pos[3];

				if (mc->surface_type == DELAUNAY_POINT_CELL)
				{
					#if !defined(EQ_CULLING)
					int64_t num_pixels = pixels[last__dim] * pixels[first_dim];
					#if !defined(MULTITHREADING)
					for (int64_t n=i_start[last__dim]; n<=i_end[last__dim]; n++)
					{
						for (int64_t m=i_start[first_dim]; m<=i_end[first_dim]; m++)
						{
							int64_t id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);
							num_pixel_intersections[id] += max_intersections_per_ray;
						}
					}
					#else // MULTITHREADING
					num_patch_intersections[ patch_id ] += max_intersections_per_ray * num_pixels;
					#endif // MULTITHREADING

					potentialIntersections[thread_id] += max_intersections_per_ray * num_pixels;
					#else // EQ_CULLING

					for (int64_t rectangle_pixel = 0; rectangle_pixel < pixels[first_dim]*pixels[last__dim]; rectangle_pixel++)
					{
						int64_t n = rectangle_pixel / pixels[first_dim];
						int64_t m = rectangle_pixel - n*pixels[first_dim];

						n += i_start[last__dim];
						m += i_start[first_dim];

						double orig[2];

						if (panel == 0) {
							orig[0] = delphi->y[m] + delta;
							orig[1] = delphi->z[n] + delta;
						}
						else if (panel == 1) {
							orig[0] = delphi->x[m] + delta;
							orig[1] = delphi->y[n] + delta;
						}
						else {
							orig[0] = delphi->x[m] + delta;
							orig[1] = delphi->z[n] + delta;
						}

						double dx = sphere_center[first_dim] - orig[0];
						double dy = sphere_center[last__dim] - orig[1];

						if (dx*dx + dy*dy > radius*radius)
							continue;

						#if !defined(MULTITHREADING)
						int64_t id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);
						num_pixel_intersections[id] += max_intersections_per_ray;
						#else
						num_patch_intersections[ patch_id ] += max_intersections_per_ray;
						#endif

						potentialIntersections[thread_id] += max_intersections_per_ray;
					}
					#endif // EQ_CULLING
				}
				else if (mc->surface_type == DELAUNAY_FACET_CELL || mc->surface_type == DELAUNAY_EDGE_CELL)
				{
					// face_group_code = 1 for lateral rectangles 1 and 2 for top and bottom polygons (and each face is subdivided into triangles)
					for (int face_group_code = 1; face_group_code <= 2; face_group_code++)
					{
						int loop_pars[] = {face_group_code, face_group_code<<1, face_group_code, face_group_code};

						// for each k-th cell's face
						for (int k = 0; k < loop_pars[0]; k++)
						{
							// decomposition of the cell's polygon faces into triangles
							for (int j = k; j < (int)ll.size() - loop_pars[1]; j += loop_pars[2])
							{
								int first = (face_group_code == 1) ? j : k;
								proj_v_pos[0] = projected_cell_vertices_coords[first_dim][first];
								proj_v_pos[1] = projected_cell_vertices_coords[first_dim][j+loop_pars[3]];
								proj_v_pos[2] = projected_cell_vertices_coords[first_dim][j+loop_pars[3]*2];

								// eventual culling of bounding rectangle of triangle vertices' projections (first axis)
								if ((v_end[0] = MAX(proj_v_pos[0], MAX(proj_v_pos[1], proj_v_pos[2]))) < 0)
									continue;
								if ((v_start[0] = MIN(proj_v_pos[0], MIN(proj_v_pos[1], proj_v_pos[2]))) >= nxyz[first_dim])
									continue;

								proj_v_pos[0] = projected_cell_vertices_coords[last__dim][first];
								proj_v_pos[1] = projected_cell_vertices_coords[last__dim][j+loop_pars[3]];
								proj_v_pos[2] = projected_cell_vertices_coords[last__dim][j+loop_pars[3]*2];

								// eventual culling of bounding rectangle of triangle vertices' projections (second axis)
								if ((v_end[1] = MAX(proj_v_pos[0], MAX(proj_v_pos[1], proj_v_pos[2]))) < 0)
									continue;
								if ((v_start[1] = MIN(proj_v_pos[0], MIN(proj_v_pos[1], proj_v_pos[2]))) >= nxyz[last__dim])
									continue;

								// eventual clipping of vertices of projected triangle
								if (v_start[0] < 0               ) v_start[0] = 0;
								if (v_end[0]   >= nxyz[first_dim]) v_end[0]   = nxyz[first_dim] - 1;
								if (v_start[1] < 0               ) v_start[1] = 0;
								if (v_end[1]   >= nxyz[last__dim]) v_end[1]   = nxyz[last__dim] - 1;

								double v_pos[2][3] = {
									{ll[first][first_dim], ll[j+loop_pars[3]][first_dim], ll[j+loop_pars[3]*2][first_dim]},
									{ll[first][last__dim], ll[j+loop_pars[3]][last__dim], ll[j+loop_pars[3]*2][last__dim]}
								};

								double A = halfArea (v_pos[0][0], v_pos[1][0], v_pos[0][1], v_pos[1][1], v_pos[0][2], v_pos[1][2]);

								for (int n=v_start[1]; n<=v_end[1]; n++)
								{
									double y = delphi->z[n] + delta;
									if (panel == 1) y = delphi->y[n] + delta;

									for (int m=v_start[0]; m<=v_end[0]; m++)
									{
										int64_t pixel_id = n*N_MAX + m;

										if (screen_buffer[ thread_id*N_MAX*N_MAX + pixel_id ])
											continue;

										double x = delphi->x[m] + delta;
										if (panel == 0) x = delphi->y[m] + delta;

										#if !defined(POINT_METHOD_2)
										double d = (v_pos[1][1]-v_pos[1][2]) * (v_pos[0][0]-v_pos[0][2]) + (v_pos[0][2]-v_pos[0][1]) * (v_pos[1][0]-v_pos[1][2]);
										double a = (v_pos[1][1]-v_pos[1][2]) * (          x-v_pos[0][2]) + (v_pos[0][2]-v_pos[0][1]) * (          y-v_pos[1][2]);

										// negative barycentric coordinate? if so, (x,y) is not inside the 2D triangle. The test is conservative, also for single prec.
										if (a < -0.000001*d) continue;

										double b = (v_pos[1][2]-v_pos[1][0]) * (          x-v_pos[0][2]) + (v_pos[0][0]-v_pos[0][2]) * (          y-v_pos[1][2]);
										// negative barycentric coordinate(s)? if so, (x,y) is not inside the 2D triangle. The test is conservative, also for single prec.
										if (b < -0.000001*d || a + b > 1.000001*d) continue;
										#else // POINT_METHOD_2

										double A1 = halfArea (          x,           y, v_pos[0][1], v_pos[1][1], v_pos[0][2], v_pos[1][2]);
										double A2 = halfArea (v_pos[0][0], v_pos[1][0],           x,           y, v_pos[0][2], v_pos[1][2]);
										double A3 = halfArea (v_pos[0][0], v_pos[1][0], v_pos[0][1], v_pos[1][1],           x,           y);

										if ((A * 0.999999 > A1+A2+A3) || (A * 1.000001 < A1+A2+A3)) continue;
										#endif // POINT_METHOD_2

										// element of the boolean mask is marked
										screen_buffer[ thread_id*N_MAX*N_MAX + pixel_id ] = 1;

										#if !defined(MULTITHREADING)
										int64_t id = (pixel_id)*numPanels + (panels[0]*int_phase + panel);
										num_pixel_intersections[id] += max_intersections_per_ray;
										#else
										num_patch_intersections[ patch_id ] += max_intersections_per_ray;
										#endif

										potentialIntersections[thread_id] += max_intersections_per_ray;

										// storage of the labelled elements of the boolean mask
										tagged_pixel_ids[ thread_id*N_MAX*N_MAX + num_tagged_pixels++ ] = pixel_id;
									}
								}
							}
						}
					}
				}
				else
				{
					char tetra_vtx_ids[4][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};

					for (int k=0; k<4; k++)
					{
						proj_v_pos[0] = projected_cell_vertices_coords[first_dim][ tetra_vtx_ids[k][0] ];
						proj_v_pos[1] = projected_cell_vertices_coords[first_dim][ tetra_vtx_ids[k][1] ];
						proj_v_pos[2] = projected_cell_vertices_coords[first_dim][ tetra_vtx_ids[k][2] ];

						// eventual culling of bounding rectangle of triangle vertices' projections (first axis)
						if ((v_end[0] = MAX(proj_v_pos[0], MAX(proj_v_pos[1], proj_v_pos[2]))) < 0)
							continue;
						if ((v_start[0] = MIN(proj_v_pos[0], MIN(proj_v_pos[1], proj_v_pos[2]))) >= nxyz[first_dim])
							continue;

						proj_v_pos[0] = projected_cell_vertices_coords[last__dim][ tetra_vtx_ids[k][0] ];
						proj_v_pos[1] = projected_cell_vertices_coords[last__dim][ tetra_vtx_ids[k][1] ];
						proj_v_pos[2] = projected_cell_vertices_coords[last__dim][ tetra_vtx_ids[k][2] ];

						// eventual culling of bounding rectangle of triangle vertices' projections (second axis)
						if ((v_end[1] = MAX(proj_v_pos[0], MAX(proj_v_pos[1], proj_v_pos[2]))) < 0)
							continue;
						if ((v_start[1] = MIN(proj_v_pos[0], MIN(proj_v_pos[1], proj_v_pos[2]))) >= nxyz[last__dim])
							continue;

						// eventual clipping of vertices of projected triangle
						if (v_start[0] < 0               ) v_start[0] = 0;
						if (v_end[0]   >= nxyz[first_dim]) v_end[0]   = nxyz[first_dim] - 1;
						if (v_start[1] < 0               ) v_start[1] = 0;
						if (v_end[1]   >= nxyz[last__dim]) v_end[1]   = nxyz[last__dim] - 1;

						double v_pos[2][3] = {
							{ll[ tetra_vtx_ids[k][0] ][first_dim], ll[ tetra_vtx_ids[k][1] ][first_dim], ll[ tetra_vtx_ids[k][2] ][first_dim]},
							{ll[ tetra_vtx_ids[k][0] ][last__dim], ll[ tetra_vtx_ids[k][1] ][last__dim], ll[ tetra_vtx_ids[k][2] ][last__dim]}
						};

						double A = halfArea (v_pos[0][0], v_pos[1][0], v_pos[0][1], v_pos[1][1], v_pos[0][2], v_pos[1][2]);

						for (int n=v_start[1]; n<=v_end[1]; n++)
						{
							double y = delphi->z[n] + delta;
							if (panel == 1) y = delphi->y[n] + delta;

							for (int m=v_start[0]; m<=v_end[0]; m++)
							{
								int64_t pixel_id = n*N_MAX + m;

								if (screen_buffer[ thread_id*N_MAX*N_MAX + pixel_id ])
									continue;

								double x = delphi->x[m] + delta;
								if (panel == 0) x = delphi->y[m] + delta;

								#if !defined(POINT_METHOD_2)
								double d = (v_pos[1][1]-v_pos[1][2]) * (v_pos[0][0]-v_pos[0][2]) + (v_pos[0][2]-v_pos[0][1]) * (v_pos[1][0]-v_pos[1][2]);
								double a = (v_pos[1][1]-v_pos[1][2]) * (          x-v_pos[0][2]) + (v_pos[0][2]-v_pos[0][1]) * (          y-v_pos[1][2]);

								// negative barycentric coordinate? if so, (x,y) is not inside the 2D triangle. The test is conservative, also for single prec.
								if (a < -0.000001*d) continue;

								double b = (v_pos[1][2]-v_pos[1][0]) * (          x-v_pos[0][2]) + (v_pos[0][0]-v_pos[0][2]) * (          y-v_pos[1][2]);
								// negative barycentric coordinate(s)? if so, (x,y) is not inside the 2D triangle. The test is conservative, also for single prec.
								if (b < -0.000001*d || a + b > 1.000001*d) continue;
								#else // POINT_METHOD_2

								double A1 = halfArea (          x,           y, v_pos[0][1], v_pos[1][1], v_pos[0][2], v_pos[1][2]);
								double A2 = halfArea (v_pos[0][0], v_pos[1][0],           x,           y, v_pos[0][2], v_pos[1][2]);
								double A3 = halfArea (v_pos[0][0], v_pos[1][0], v_pos[0][1], v_pos[1][1],           x,           y);

								if ((A * 0.999999 > A1+A2+A3) || (A * 1.000001 < A1+A2+A3)) continue;
								#endif // POINT_METHOD_2

								// element of the boolean mask is marked
								screen_buffer[ thread_id*N_MAX*N_MAX + pixel_id ] = 1;

								#if !defined(MULTITHREADING)
								int64_t id = (pixel_id)*numPanels + (panels[0]*int_phase + panel);
								num_pixel_intersections[id] += max_intersections_per_ray;
								#else
								num_patch_intersections[ patch_id ] += max_intersections_per_ray;
								#endif

								potentialIntersections[thread_id] += max_intersections_per_ray;

								// storage of the labelled elements of the boolean mask
								tagged_pixel_ids[ thread_id*N_MAX*N_MAX + num_tagged_pixels++ ] = pixel_id;
							}
						}
					}
				}
				for (int k = 0; k < num_tagged_pixels; k++)
				{
					// resetting the elements of the labelling buffer
					screen_buffer[ thread_id*N_MAX*N_MAX + tagged_pixel_ids[k] ] = 0;
				}
				#endif // RAY_VS_CELL_TESTS_CULLING
			}
		}
	}
}

#else // MINIMIZE_MEMORY

void SkinSurface::getPatchPreIntersectionData (int64_t nxyz[3], int panels[2],
											   int thread_id, int potentialIntersections[])
{
	if (getNumPatches() == 0)
	{
		return;
	}

	int64_t NX = nxyz[0];
	int64_t NY = nxyz[1];
	int64_t NZ = nxyz[2];
	int64_t N_MAX = MAX(NX, MAX(NY, NZ));

	int numPanels = panels[0] + panels[1];

	double pa[3];


	potentialIntersections[thread_id] = 0;

	// Perform the per-patch ray casting
	for (int64_t it = thread_id; it < numMixedCells; it += conf.numThreads)
	{
		MixedCell *mc = mixedComplex[it];

		if (mc->surface_type == SKIP_CELL)
			continue;

		double *Q = mc->quadric;

		// compute the bounding box of the object
		double downx = INFINITY;
		double downy = INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = -INFINITY;
		double upz = -INFINITY;

		double *sphere_center;
		double radius;

		if (mc->surface_type == DELAUNAY_POINT_CELL)
		{
			int index = ((Del0Cell*)mc)->id;

			// avoid mapping a non atom (bounding box virtual atom)
			if (index >= delphi->atoms.size())
			{
				continue;
			}
			sphere_center = delphi->atoms[index].pos;
			radius = delphi->atoms[index].radius;

			downx = sphere_center[0]-radius;
			downy = sphere_center[1]-radius;
			downz = sphere_center[2]-radius;

			upx = sphere_center[0]+radius;
			upy = sphere_center[1]+radius;
			upz = sphere_center[2]+radius;
		}
		else
		{
			// mixed cell points
			vector<double*> &ll = mc->mc_points;

			for (unsigned int i=0; i<ll.size(); i++)
			{
				downx = MIN(downx,ll[i][0]);
				downy = MIN(downy,ll[i][1]);
				downz = MIN(downz,ll[i][2]);

				upx = MAX(upx,ll[i][0]);
				upy = MAX(upy,ll[i][1]);
				upz = MAX(upz,ll[i][2]);
			}
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
				int first_dim, last__dim;
				int varying_coord;

				if (panel == 0)
				{
					first_dim = 1, last__dim = 2;
					pa[0] = delphi->xmin + delta;
					varying_coord = 0;
				}
				else if (panel == 1)
				{
					first_dim = 0, last__dim = 1;
					pa[2] = delphi->zmin + delta;
					varying_coord = 2;
				}
				else
				{
					first_dim = 0, last__dim = 2;
					pa[1] = delphi->ymin + delta;
					varying_coord = 1;
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

					// double cache[9] = {o3+o3,o3*o3,o2+o2,2*o2*o3,o2*o2,o1+o1,2*o1*o3,2*o1*o2,o1*o1};
					double t[2];

					bool ff;
					if (panel == 0) ff = rayQuadricIntersectionX(Q, pa, t);
					if (panel == 1) ff = rayQuadricIntersectionZ(Q, pa, t);
					if (panel == 2) ff = rayQuadricIntersectionY(Q, pa, t);
					if (!ff) continue;

					for (int i=0; i<2; i++)
					{
						double intPoint[3] = {pa[0], pa[1], pa[2]};
						intPoint[varying_coord] += t[i];

						if (!isFeasible(mc,intPoint)) continue;

						#if !defined(MULTITHREADING)
						int64_t panel_pixel_id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);
						++num_pixel_intersections[panel_pixel_id];
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


void SkinSurface::getPatchIntersectionData (int64_t nxyz[3], int panels[2],
											int thread_id, int *netIntersections)
{
	if (getNumPatches() == 0)
	{
		return;
	}

	int64_t NX = nxyz[0];
	int64_t NY = nxyz[1];
	int64_t NZ = nxyz[2];
	int64_t N_MAX = MAX(NX, MAX(NY, NZ));

	int numPanels = panels[0] + panels[1];

	double pa[3];


	*netIntersections = 0;

	// Perform the per-patch ray casting
	for (int it = thread_id; it < numMixedCells; it += conf.numThreads)
	{
		MixedCell *mc = mixedComplex[it];

		if (mc->surface_type == SKIP_CELL)
			continue;

		double *Q = mc->quadric;

		#ifdef RAY_VS_CELL_TESTS_CULLING
		vector<double*> &ll = mc->mc_points;

		int projected_cell_vertices_coords[3][ll.size()];

		int thd_id = thread_id;
		#endif

		double double_normal_sign = (mc->surface_type == DELAUNAY_TETRA_CELL) ? -2.0 : 2.0;

		// compute the bounding box of the object
		double downx = INFINITY;
		double downy = INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = -INFINITY;
		double upz = -INFINITY;

		double *sphere_center;
		double radius;

		if (mc->surface_type == DELAUNAY_POINT_CELL)
		{
			int index = ((Del0Cell*)mc)->id;

			// avoid mapping a non atom (bounding box virtual atom)
			if (index >= delphi->atoms.size())
			{
				continue;
			}
			sphere_center = delphi->atoms[index].pos;
			radius = delphi->atoms[index].radius;

			downx = sphere_center[0]-radius;
			downy = sphere_center[1]-radius;
			downz = sphere_center[2]-radius;

			upx = sphere_center[0]+radius;
			upy = sphere_center[1]+radius;
			upz = sphere_center[2]+radius;
		}
		else
		{
			vector<double*> &ll = mc->mc_points;

			for (unsigned int i=0; i<ll.size(); i++)
			{
				downx = MIN(downx,ll[i][0]);
				downy = MIN(downy,ll[i][1]);
				downz = MIN(downz,ll[i][2]);

				upx = MAX(upx,ll[i][0]);
				upy = MAX(upy,ll[i][1]);
				upz = MAX(upz,ll[i][2]);
			}
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

		i_start[0] = (int64_t)rintp(fmax(0., (downx-delphi->xmin))*delphi->scale);
		i_start[1] = (int64_t)rintp(fmax(0., (downy-delphi->ymin))*delphi->scale);
		i_start[2] = (int64_t)rintp(fmax(0., (downz-delphi->zmin))*delphi->scale);

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

			#ifdef RAY_VS_PATCH_TESTS_CULLING
			if (mc->surface_type == DELAUNAY_FACET_CELL ||
				mc->surface_type == DELAUNAY_EDGE_CELL  ||
				mc->surface_type == DELAUNAY_TETRA_CELL)
			{
				for (unsigned int i=0; i<ll.size(); i++)
				{
					// storing the floating point projections of the vertices of the current cell
					projected_cell_vertices_coords[0][i] = (int)rintp((ll[i][0] - (delta + delphi->xmin))*delphi->scale);
					projected_cell_vertices_coords[1][i] = (int)rintp((ll[i][1] - (delta + delphi->ymin))*delphi->scale);
					projected_cell_vertices_coords[2][i] = (int)rintp((ll[i][2] - (delta + delphi->zmin))*delphi->scale);
				}
			}
			#endif
			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				double ray_dir;

				int first_dim, last__dim;
				int varying_coord;

				if (panel == 0)
				{
					first_dim = 1, last__dim = 2;
					pa[0] = delphi->xmin + delta;
					ray_dir = delphi->xmax - delphi->xmin;
					varying_coord = 0;
				}
				else if (panel == 1)
				{
					first_dim = 0, last__dim = 1;
					pa[2] = delphi->zmin + delta;
					ray_dir = delphi->zmax - delphi->zmin;
					varying_coord = 2;
				}
				else
				{
					first_dim = 0, last__dim = 2;
					pa[1] = delphi->ymin + delta;
					ray_dir = delphi->ymax - delphi->ymin;
					varying_coord = 1;
				}
				#if defined(MULTITHREADING)
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);
				#if !defined(SINGLE_PASS_RT)
				int first_id = first_patch_intersection_index[patch_id];
				#endif
				#endif

				#if !defined(RAY_VS_PATCH_TESTS_CULLING)
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
					#ifdef EQ_CULLING
					if (mc->surface_type == DELAUNAY_POINT_CELL)
					{
						double dx = sphere_center[first_dim] - pa[first_dim];
						double dy = sphere_center[last__dim] - pa[last__dim];

						if (dx*dx + dy*dy > radius*radius)
							continue;
					}
					#endif

					// double cache[9] = {o3+o3,o3*o3,o2+o2,2*o2*o3,o2*o2,o1+o1,2*o1*o3,2*o1*o2,o1*o1};
					double t[2];

					bool ff;
					if (panel == 0) ff = rayQuadricIntersectionX(Q, pa, t);
					if (panel == 1) ff = rayQuadricIntersectionZ(Q, pa, t);
					if (panel == 2) ff = rayQuadricIntersectionY(Q, pa, t);
					if (!ff) continue;

					for (int i=0; i<2; i++)
					{
						double intPoint[3] = {pa[0], pa[1], pa[2]};
						intPoint[varying_coord] += t[i];

						if (!isFeasible(mc,intPoint)) continue;

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

						int_p->first = t[i] / ray_dir;

						VERTEX_TYPE *normal = nullptr;
						if (computeNormals)
						{
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);

							normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

							normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
							normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
							normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
						}
						int_p->second = normal;

						++*netIntersections;

						#else // SINGLE_PASS_RT

						VERTEX_TYPE *normal = nullptr;
						if (computeNormals)
						{
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);

							normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

							normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
							normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
							normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
						}
						temp_intersections_buffer[ patch_id ].push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i]/ray_dir, normal));

						intersection_pixel_id[ patch_id ].push_back(panel_pixel_id);

						++*netIntersections;

						#endif // SINGLE_PASS_RT
					}
				}
				#else // RAY_VS_PATCH_TESTS_CULLING
				int num_tagged_pixels = 0;

				int v_start[2], v_end[2];

				int proj_v_pos[3];

				if (mc->surface_type == DELAUNAY_POINT_CELL)
				{
					for (int64_t rectangle_pixel = 0; rectangle_pixel < pixels[first_dim]*pixels[last__dim]; rectangle_pixel++)
					{
						int64_t n = rectangle_pixel / pixels[first_dim];
						int64_t m = rectangle_pixel - n*pixels[first_dim];

						n += i_start[last__dim];
						m += i_start[first_dim];

						if (panel == 0) {
							pa[1] = delphi->y[m] + delta;
							pa[2] = delphi->z[n] + delta;
						}
						else if (panel == 1) {
							pa[0] = delphi->x[m] + delta;
							pa[1] = delphi->y[n] + delta;
						}
						else {
							pa[0] = delphi->x[m] + delta;
							pa[2] = delphi->z[n] + delta;
						}
						#ifdef EQ_CULLING
						double dx = sphere_center[first_dim] - pa[first_dim];
						double dy = sphere_center[last__dim] - pa[last__dim];

						if (dx*dx + dy*dy > radius*radius)
							continue;
						#endif

						// double cache[9] = {o3+o3,o3*o3,o2+o2,2*o2*o3,o2*o2,o1+o1,2*o1*o3,2*o1*o2,o1*o1};
						double t[2]

						bool ff;
						if (panel == 0) ff = rayQuadricIntersectionX(Q, pa, t);
						if (panel == 1) ff = rayQuadricIntersectionZ(Q, pa, t);
						if (panel == 2) ff = rayQuadricIntersectionY(Q, pa, t);
						if (!ff) continue;

						for (int i=0; i<2; i++)
						{
							double intPoint[3] = {o1, o2, o3};
							intPoint[varying_coord] += t[i];

							if (!isFeasible(mc,intPoint)) continue;

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

							int_p->first = t[i] / ray_dir;

							VERTEX_TYPE *normal = nullptr;
							if (computeNormals)
							{
								normalsBuffers[thread_id].push_back(0.);
								normalsBuffers[thread_id].push_back(0.);
								normalsBuffers[thread_id].push_back(0.);

								normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

								normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
								normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
								normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
							}
							int_p->second = normal;

							++*netIntersections;

							#else // SINGLE_PASS_RT

							VERTEX_TYPE *normal = nullptr;
							if (computeNormals)
							{
								normalsBuffers[thread_id].push_back(0.);
								normalsBuffers[thread_id].push_back(0.);
								normalsBuffers[thread_id].push_back(0.);

								normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

								normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
								normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
								normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
							}
							temp_intersections_buffer[ patch_id ].push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i]/ray_dir, normal));

							intersection_pixel_id[ patch_id ].push_back(panel_pixel_id);

							++*netIntersections;

							#endif // SINGLE_PASS_RT
						}
					}
				}
				else if (mc->surface_type == DELAUNAY_FACET_CELL || mc->surface_type == DELAUNAY_EDGE_CELL)
				{
					// face_group_code = 1 for lateral rectangles 1 and 2 for top and bottom polygons (and each face is subdivided into triangles)
					for (int face_group_code = 1; face_group_code <= 2; face_group_code++)
					{
						int loop_pars[] = {face_group_code, face_group_code<<1, face_group_code, face_group_code};

						// for each k-th cell's face
						for (int k = 0; k < loop_pars[0]; k++)
						{
							// decomposition of the cell's polygon faces into triangles
							for (int j = k; j < (int)ll.size() - loop_pars[1]; j += loop_pars[2])
							{
								int first = (face_group_code == 1) ? j : k;
								proj_v_pos[0] = projected_cell_vertices_coords[first_dim][first];
								proj_v_pos[1] = projected_cell_vertices_coords[first_dim][j+loop_pars[3]];
								proj_v_pos[2] = projected_cell_vertices_coords[first_dim][j+loop_pars[3]*2];

								// eventual culling of bounding rectangle of triangle vertices' projections (first axis)
								if ((v_end[0] = MAX(proj_v_pos[0], MAX(proj_v_pos[1], proj_v_pos[2]))) < 0)
									continue;
								if ((v_start[0] = MIN(proj_v_pos[0], MIN(proj_v_pos[1], proj_v_pos[2]))) >= nxyz[first_dim])
									continue;

								proj_v_pos[0] = projected_cell_vertices_coords[last__dim][first];
								proj_v_pos[1] = projected_cell_vertices_coords[last__dim][j+loop_pars[3]];
								proj_v_pos[2] = projected_cell_vertices_coords[last__dim][j+loop_pars[3]*2];

								// eventual culling of bounding rectangle of triangle vertices' projections (second axis)
								if ((v_end[1] = MAX(proj_v_pos[0], MAX(proj_v_pos[1], proj_v_pos[2]))) < 0)
									continue;
								if ((v_start[1] = MIN(proj_v_pos[0], MIN(proj_v_pos[1], proj_v_pos[2]))) >= nxyz[last__dim])
									continue;

								// eventual clipping of vertices of projected triangle
								if (v_start[0] < 0               ) v_start[0] = 0;
								if (v_end[0]   >= nxyz[first_dim]) v_end[0]   = nxyz[first_dim] - 1;
								if (v_start[1] < 0               ) v_start[1] = 0;
								if (v_end[1]   >= nxyz[last__dim]) v_end[1]   = nxyz[last__dim] - 1;

								double v_pos[2][3] = {
									{ll[first][first_dim], ll[j+loop_pars[3]][first_dim], ll[j+loop_pars[3]*2][first_dim]},
									{ll[first][last__dim], ll[j+loop_pars[3]][last__dim], ll[j+loop_pars[3]*2][last__dim]}
								};

								double A = halfArea (v_pos[0][0], v_pos[1][0], v_pos[0][1], v_pos[1][1], v_pos[0][2], v_pos[1][2]);

								for (int n=v_start[1]; n<=v_end[1]; n++)
								{
									double y;

									if (panel == 0)
										y = pa[2] = delphi->z[n] + delta;
									else if (panel == 1)
										y = pa[1] = delphi->y[n] + delta;
									else
										y = pa[2] = delphi->z[n] + delta;

									for (int m=v_start[0]; m<=v_end[0]; m++)
									{
										int64_t pixel_id = n*N_MAX + m;

										if (screen_buffer[ thd_id*N_MAX*N_MAX + pixel_id ])
											continue;

										double x;

										if (panel == 0)
											x = pa[1] = delphi->y[m] + delta;
										else if (panel == 1)
											x = pa[0] = delphi->x[m] + delta;
										else
											x = pa[0] = delphi->x[m] + delta;

										#if !defined(POINT_METHOD_2)
										double d = (v_pos[1][1]-v_pos[1][2]) * (v_pos[0][0]-v_pos[0][2]) + (v_pos[0][2]-v_pos[0][1]) * (v_pos[1][0]-v_pos[1][2]);
										double a = (v_pos[1][1]-v_pos[1][2]) * (          x-v_pos[0][2]) + (v_pos[0][2]-v_pos[0][1]) * (          y-v_pos[1][2]);

										// negative barycentric coordinate? if so, (x,y) is not inside the 2D triangle. The test is conservative, also for single prec.
										if (a < -0.000001*d) continue;

										double b = (v_pos[1][2]-v_pos[1][0]) * (          x-v_pos[0][2]) + (v_pos[0][0]-v_pos[0][2]) * (          y-v_pos[1][2]);
										// negative barycentric coordinate(s)? if so, (x,y) is not inside the 2D triangle. The test is conservative, also for single prec.
										if (b < -0.000001*d || a + b > 1.000001*d) continue;
										#else // POINT_METHOD_2

										double A1 = halfArea (          x,           y, v_pos[0][1], v_pos[1][1], v_pos[0][2], v_pos[1][2]);
										double A2 = halfArea (v_pos[0][0], v_pos[1][0],           x,           y, v_pos[0][2], v_pos[1][2]);
										double A3 = halfArea (v_pos[0][0], v_pos[1][0], v_pos[0][1], v_pos[1][1],           x,           y);

										if ((A * 0.999999 > A1+A2+A3) || (A * 1.000001 < A1+A2+A3)) continue;
										#endif // POINT_METHOD_2

										// element of the boolean mask is marked
										screen_buffer[ thd_id*N_MAX*N_MAX + pixel_id ] = 1;

										// storage of the labelled elements of the boolean mask
										tagged_pixel_ids[ thd_id*N_MAX*N_MAX + num_tagged_pixels++ ] = pixel_id;

										int old_numIntersections = *netIntersections;

										// double cache[12] = {o3+o3,o3*o3,o2+o2,2*o2*o3,o2*o2,o1+o1,2*o1*o3,2*o1*o2,o1*o1};
										double t[2];

										bool ff;
										if (panel == 0) ff = rayQuadricIntersectionX(Q, pa, t);
										if (panel == 1) ff = rayQuadricIntersectionZ(Q, pa, t);
										if (panel == 2) ff = rayQuadricIntersectionY(Q, pa, t);
										if (!ff) continue;

										for (int i=0; i<2; i++)
										{
											double intPoint[3] = {o1, o2, o3};
											intPoint[varying_coord] += t[i];

											if (!isFeasible(mc,intPoint)) continue;

											// the n- and m-dependent integers are rarely computed here since it is very likely
											// that the ray misses the patch
											int64_t panel_pixel_id = (pixel_id)*numPanels + (panels[0]*int_phase + panel);

											#if !defined(SINGLE_PASS_RT)

											#if !defined(MULTITHREADING)
											pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p = &pixel_intersections[panel_pixel_id][ num_pixel_intersections[panel_pixel_id]++ ];
											#else
											pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p = &temp_intersections_buffer[ first_id + num_patch_intersections[patch_id] ];
											intersection_pixel_id[ first_id + num_patch_intersections[patch_id]++ ] = panel_pixel_id;
											#endif

											int_p->first = t[i] / ray_dir;

											VERTEX_TYPE *normal = nullptr;
											if (computeNormals)
											{
												normalsBuffers[thread_id].push_back(0.);
												normalsBuffers[thread_id].push_back(0.);
												normalsBuffers[thread_id].push_back(0.);

												normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

												normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
												normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
												normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
											}
											int_p->second = normal;

											++*netIntersections;

											#else // SINGLE_PASS_RT

											VERTEX_TYPE *normal = nullptr;
											if (computeNormals)
											{
												normalsBuffers[thread_id].push_back(0.);
												normalsBuffers[thread_id].push_back(0.);
												normalsBuffers[thread_id].push_back(0.);

												normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

												normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
												normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
												normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
											}
											temp_intersections_buffer[ patch_id ].push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i]/ray_dir, normal));

											intersection_pixel_id[ patch_id ].push_back(panel_pixel_id);

											++*netIntersections;

											#endif // SINGLE_PASS_RT
										}
									}
								}
							}
						}
					}
				}
				else
				{
					char tetra_vtx_ids[4][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};

					for (int k=0; k<4; k++)
					{
						proj_v_pos[0] = projected_cell_vertices_coords[first_dim][ tetra_vtx_ids[k][0] ];
						proj_v_pos[1] = projected_cell_vertices_coords[first_dim][ tetra_vtx_ids[k][1] ];
						proj_v_pos[2] = projected_cell_vertices_coords[first_dim][ tetra_vtx_ids[k][2] ];

						// eventual culling of bounding rectangle of triangle vertices' projections (first axis)
						if ((v_end[0] = MAX(proj_v_pos[0], MAX(proj_v_pos[1], proj_v_pos[2]))) < 0)
							continue;
						if ((v_start[0] = MIN(proj_v_pos[0], MIN(proj_v_pos[1], proj_v_pos[2]))) >= nxyz[first_dim])
							continue;

						proj_v_pos[0] = projected_cell_vertices_coords[last__dim][ tetra_vtx_ids[k][0] ];
						proj_v_pos[1] = projected_cell_vertices_coords[last__dim][ tetra_vtx_ids[k][1] ];
						proj_v_pos[2] = projected_cell_vertices_coords[last__dim][ tetra_vtx_ids[k][2] ];

						// eventual culling of bounding rectangle of triangle vertices' projections (second axis)
						if ((v_end[1] = MAX(proj_v_pos[0], MAX(proj_v_pos[1], proj_v_pos[2]))) < 0)
							continue;
						if ((v_start[1] = MIN(proj_v_pos[0], MIN(proj_v_pos[1], proj_v_pos[2]))) >= nxyz[last__dim])
							continue;

						// eventual clipping of vertices of projected triangle
						if (v_start[0] < 0               ) v_start[0] = 0;
						if (v_end[0]   >= nxyz[first_dim]) v_end[0]   = nxyz[first_dim] - 1;
						if (v_start[1] < 0               ) v_start[1] = 0;
						if (v_end[1]   >= nxyz[last__dim]) v_end[1]   = nxyz[last__dim] - 1;

						double v_pos[2][3] = {
							{ll[ tetra_vtx_ids[k][0] ][first_dim], ll[ tetra_vtx_ids[k][1] ][first_dim], ll[ tetra_vtx_ids[k][2] ][first_dim]},
							{ll[ tetra_vtx_ids[k][0] ][last__dim], ll[ tetra_vtx_ids[k][1] ][last__dim], ll[ tetra_vtx_ids[k][2] ][last__dim]}
						};

						double A = halfArea (v_pos[0][0], v_pos[1][0], v_pos[0][1], v_pos[1][1], v_pos[0][2], v_pos[1][2]);

						for (int n=v_start[1]; n<=v_end[1]; n++)
						{
							double y;

							if (panel == 0)
								y = pa[2] = delphi->z[n] + delta;
							else if (panel == 1)
								y = pa[1] = delphi->y[n] + delta;
							else
								y = pa[2] = delphi->z[n] + delta;

							for (int m=v_start[0]; m<=v_end[0]; m++)
							{
								int64_t pixel_id = n*N_MAX + m;

								if (screen_buffer[ thd_id*N_MAX*N_MAX + pixel_id ])
									continue;

								double x;

								if (panel == 0)
									x = pa[1] = delphi->y[m] + delta;
								else if (panel == 1)
									x = pa[0] = delphi->x[m] + delta;
								else
									x = pa[0] = delphi->x[m] + delta;

								#if !defined(POINT_METHOD_2)
								double d = (v_pos[1][1]-v_pos[1][2]) * (v_pos[0][0]-v_pos[0][2]) + (v_pos[0][2]-v_pos[0][1]) * (v_pos[1][0]-v_pos[1][2]);
								double a = (v_pos[1][1]-v_pos[1][2]) * (          x-v_pos[0][2]) + (v_pos[0][2]-v_pos[0][1]) * (          y-v_pos[1][2]);

								// negative barycentric coordinate? if so, (x,y) is not inside the 2D triangle. The test is conservative, also for single prec.
								if (a < -0.000001*d) continue;

								double b = (v_pos[1][2]-v_pos[1][0]) * (          x-v_pos[0][2]) + (v_pos[0][0]-v_pos[0][2]) * (          y-v_pos[1][2]);
								// negative barycentric coordinate(s)? if so, (x,y) is not inside the 2D triangle. The test is conservative, also for single prec.
								if (b < -0.000001*d || a + b > 1.000001*d) continue;
								#else // POINT_METHOD_2

								double A1 = halfArea (          x,           y, v_pos[0][1], v_pos[1][1], v_pos[0][2], v_pos[1][2]);
								double A2 = halfArea (v_pos[0][0], v_pos[1][0],           x,           y, v_pos[0][2], v_pos[1][2]);
								double A3 = halfArea (v_pos[0][0], v_pos[1][0], v_pos[0][1], v_pos[1][1],           x,           y);

								if ((A * 0.999999 > A1+A2+A3) || (A * 1.000001 < A1+A2+A3)) continue;
								#endif // POINT_METHOD_2

								// element of the boolean mask is marked
								screen_buffer[ thd_id*N_MAX*N_MAX + pixel_id ] = 1;

								// storage of the labelled elements of the boolean mask
								tagged_pixel_ids[ thd_id*N_MAX*N_MAX + num_tagged_pixels++ ] = pixel_id;

								// double cache[12] = {o3+o3,o3*o3,o2+o2,2*o2*o3,o2*o2,o1+o1,2*o1*o3,2*o1*o2,o1*o1};
								double t[2];

								bool ff;
								if (panel == 0) ff = rayQuadricIntersectionX(Q, pa, t);
								if (panel == 1) ff = rayQuadricIntersectionZ(Q, pa, t);
								if (panel == 2) ff = rayQuadricIntersectionY(Q, pa, t);
								if (!ff) continue;

								for (int i=0; i<2; i++)
								{
									double intPoint[3] = {o1, o2, o3};
									intPoint[varying_coord] += t[i];

									if (!isFeasible(mc,intPoint)) continue;

									// the n- and m-dependent integers are rarely computed here since it is very likely
									// that the ray misses the patch
									int64_t panel_pixel_id = (pixel_id)*numPanels + (panels[0]*int_phase + panel);

									#if !defined(SINGLE_PASS_RT)

									#if !defined(MULTITHREADING)
									pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p = &pixel_intersections[panel_pixel_id][ num_pixel_intersections[panel_pixel_id]++ ];
									#else
									pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p = &temp_intersections_buffer[ first_id + num_patch_intersections[patch_id] ];
									intersection_pixel_id[ first_id + num_patch_intersections[patch_id]++ ] = panel_pixel_id;
									#endif

									int_p->first = t[i] / ray_dir;

									VERTEX_TYPE *normal = nullptr;
									if (computeNormals)
									{
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);

										normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

										normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
										normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
										normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
									}
									int_p->second = normal;

									++*netIntersections;

									#else // SINGLE_PASS_RT

									VERTEX_TYPE *normal = nullptr;
									if (computeNormals)
									{
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);

										normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

										normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
										normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
										normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
									}
									temp_intersections_buffer[ patch_id ].push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i]/ray_dir, normal));

									intersection_pixel_id[ patch_id ].push_back(panel_pixel_id);

									++*netIntersections;

									#endif // SINGLE_PASS_RT
								}
							}
						}
					}
				}
				for (int k = 0; k < num_tagged_pixels; k++)
				{
					// resetting the elements of the labelling buffer
					screen_buffer[ thd_id*N_MAX*N_MAX + tagged_pixel_ids[k] ] = 0;
				}
				#endif // RAY_VS_PATCH_TESTS_CULLING
			}
		}
	}
}


/*
#if !defined(SINGLE_PASS_RT)

// This can be used if the normals' buffer common to the threads is filled after calculating and storing the intersections
void SkinSurface::getPatchNormalsAtIntersections (int64_t nxyz[3], int panels[2], int thread_id)
{
	if (getNumPatches() == 0 || !computeNormals)
	{
		return;
	}

	int NX = nxyz[0];
	int NY = nxyz[1];
	int NZ = nxyz[2];
	int N_MAX = MAX(NX, MAX(NY, NZ));

	int numPanels = panels[0] + panels[1];


	for (int it = thread_id; it < numMixedCells; it += conf.numThreads)
	{
		MixedCell *mc = mixedComplex[it];

		if (mc->surface_type == SKIP_CELL)
			continue;

		double *Q = mc->quadric;

		double double_normal_sign = (mc->surface_type == DELAUNAY_TETRA_CELL) ? -2.0 : 2.0;

		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			double delta = (int_phase == 0) ? 0. : delta_accurate_triangulation - delphi->hside;

			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);
				int first_id = first_patch_intersection_index[patch_id];

				for (int i=0; i<num_patch_intersections[patch_id]; i++)
				{
					int64_t panel_pixel_id = intersection_pixel_id[ first_id + i ];

					pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p = &temp_intersections_buffer[ first_id + i ];

					// reconstruction of intersection point intPoint[] so as to calculate and store the normal
					double intPoint[3];
					double ray_dir;

					int64_t pixel_id = (panel_pixel_id - (panels[0]*int_phase + panel)) / numPanels;
					int n = pixel_id / N_MAX;
					int m = pixel_id - n*N_MAX;

					if (panel == 0)
					{
						intPoint[0] = delphi->x[0] + delta;
						intPoint[1] = delphi->y[m] + delta;
						intPoint[2] = delphi->z[n] + delta;
						ray_dir = delphi->x[NX-1] - delphi->x[0];
						intPoint[0] += int_p->first * ray_dir;
					}
					else if (panel == 1)
					{
						intPoint[0] = delphi->x[m] + delta;
						intPoint[1] = delphi->y[n] + delta;
						intPoint[2] = delphi->z[0] + delta;
						ray_dir = delphi->z[NZ-1] - delphi->z[0];
						intPoint[2] += int_p->first * ray_dir;
					}
					else
					{
						intPoint[0] = delphi->x[m] + delta;
						intPoint[1] = delphi->y[0] + delta;
						intPoint[2] = delphi->z[n] + delta;
						ray_dir = delphi->y[NY-1] - delphi->y[0];
						intPoint[1] += int_p->first * ray_dir;
					}
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);

					VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

					normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
					normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
					normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
					int_p->second = normal;
				}
			}
		}
	}
}

#else // SINGLE_PASS_RT

// This can be used if the normals' buffer common to the threads is filled after calculating and storing the intersections
void SkinSurface::getPatchNormalsAtIntersections (int64_t nxyz[3], int panels[2], int thread_id)
{
	if (getNumPatches() == 0 || !computeNormals)
	{
		return;
	}

	int NX = nxyz[0];
	int NY = nxyz[1];
	int NZ = nxyz[2];
	int N_MAX = MAX(NX, MAX(NY, NZ));

	int numPanels = panels[0] + panels[1];


	for (int it = thread_id; it < numMixedCells; it += conf.numThreads)
	{
		MixedCell *mc = mixedComplex[it];

		if (mc->surface_type == SKIP_CELL)
			continue;

		double *Q = mc->quadric;

		double double_normal_sign = (mc->surface_type == DELAUNAY_TETRA_CELL) ? -2.0 : 2.0;

		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			double delta = (int_phase == 0) ? 0. : delta_accurate_triangulation - delphi->hside;

			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);

				for (int i=0; i<intersection_pixel_id[patch_id].size(); i++)
				{
					int64_t panel_pixel_id = intersection_pixel_id[ patch_id ][i];

					pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p = &temp_intersections_buffer[patch_id][i];

					// reconstruction of intersection point intPoint[] so as to calculate and store the normal
					double intPoint[3];
					double ray_dir;

					int64_t pixel_id = (panel_pixel_id - (panels[0]*int_phase + panel)) / numPanels;
					int n = pixel_id / N_MAX;
					int m = pixel_id - n*N_MAX;

					if (panel == 0)
					{
						intPoint[0] = delphi->x[0] + delta;
						intPoint[1] = delphi->y[m] + delta;
						intPoint[2] = delphi->z[n] + delta;
						ray_dir = delphi->x[NX-1] - delphi->x[0];
						intPoint[0] += int_p->first * ray_dir;
					}
					else if (panel == 1)
					{
						intPoint[0] = delphi->x[m] + delta;
						intPoint[1] = delphi->y[n] + delta;
						intPoint[2] = delphi->z[0] + delta;
						ray_dir = delphi->z[NZ-1] - delphi->z[0];
						intPoint[2] += int_p->first * ray_dir;
					}
					else
					{
						intPoint[0] = delphi->x[m] + delta;
						intPoint[1] = delphi->y[0] + delta;
						intPoint[2] = delphi->z[n] + delta;
						ray_dir = delphi->y[NY-1] - delphi->y[0];
						intPoint[1] += int_p->first * ray_dir;
					}
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);

					VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

					normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
					normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
					normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
					int_p->second = normal;
				}
			}
		}
	}
}

#endif // SINGLE_PASS_RT
*/


bool SkinSurface::isFeasible(MixedCell *mc, double point[3])
{
	if (mc->surface_type == DELAUNAY_POINT_CELL)
	{
		Del0Cell *mc0 = (Del0Cell*)mc;

		for (int it=0; it < (int)mc0->planes.size(); it++)
		{
			double *plane = mc0->planes[it];
			double val = DOT(plane,point) + plane[3];
			if (val < 0)
				return false;
		}
		return true;
	}
	else if (mc->surface_type == DELAUNAY_EDGE_CELL)
	{
		Del1Cell *mc1 = (Del1Cell*)mc;

		for (int it = 0; it < (int)mc1->planes.size(); it++)
		{
			double *plane = mc1->planes[it];
			double val = DOT(plane,point) + plane[3];
			if (val > 0)
				return false;
		}
		double val = DOT(mc1->upper,point) + mc1->upper[3];
		if (val > 0)
			return false;
		
		val = DOT(mc1->lower,point) + mc1->lower[3];
		if (val > 0)
			return false;
		
		return true;
	}
	else if (mc->surface_type == DELAUNAY_FACET_CELL)
	{
		Del2Cell *mc2 = (Del2Cell*)mc;
		for (int it = 0; it<3; it++)
		{
			double *plane = mc2->planes[it];
			double val = DOT(plane,point) + plane[3];
			if (val > 0)
				return false;
		}
		double val = DOT(mc2->upper,point) + mc2->upper[3];
		if (val > 0)
			return false;
		
		val = DOT(mc2->lower,point) + mc2->lower[3];
		if (val > 0)
			return false;
		
		return true;
	}
	else if (mc->surface_type == DELAUNAY_TETRA_CELL)
	{
		Del3Cell *mc3 = (Del3Cell*)mc;
		for (int it = 0; it<4; it++)
		{
			double *plane = mc3->planes[it];
			double val = DOT(plane,point) + plane[3];
			if (val > 0)
				return false;
		}
		return true;
	}
	else
	{
		cout << endl << ERR << "Unrecognized patch type!";
		return false;
	}
}


/*
inline bool SkinSurface::rayQuadricIntersection(double *Q, double *i1, double *i2, double *cache)
{
	double a=0,b=0,c=0;
	
	//for (int i=0;i<4;i++)
	//	for (int j=0;j<4;j++)
	//		a += dir[i]*Q[i][j]*dir[j];
	
	//inline using symmetry and caching
	//double q44=Q[3][3],q33=Q[2][2],q22=Q[1][1],q11=Q[0][0],q34=Q[2][3],q24=Q[1][3],q23=Q[1][2],q14=Q[0][3],q13=Q[0][2],q12=Q[0][1];
	double q44=Q[3],q33=Q[2],q22=Q[1],q11=Q[0],q34=Q[6],q24=Q[8],q23=Q[5],q14=Q[9],q13=Q[7],q12=Q[4];
 
	a = cache[0]*q44+cache[1]*q34+cache[2]*q33+cache[3]*q24+cache[4]*q23+cache[5]*q22+cache[6]*q14+cache[7]*q13+cache[8]*q12+cache[9]*q11;
 
	// bug fix (extremely rare)
	if (a == 0)
		return false;
 
	//for (int i=0;i<4;i++)
	//	for (int j=0;j<4;j++)
	//		b += orig[i]*Q[i][j]*dir[j];
 
	b = cache[20]*q44+cache[21]*q34+cache[22]*q33+cache[23]*q24+cache[24]*q23+cache[25]*q22+cache[26]*q14+cache[27]*q13+cache[28]*q12+cache[29]*q11;
 
	//for (int i=0;i<4;i++)
	//	for (int j=0;j<4;j++)
	//		c += orig[i]*Q[i][j]*orig[j];
 
	c = cache[10]*q44+cache[11]*q34+cache[12]*q33+cache[13]*q24+cache[14]*q23+cache[15]*q22+cache[16]*q14+cache[17]*q13+cache[18]*q12+cache[19]*q11;
 
	double delta = b*b-a*c;
 
	if (delta<0)
		return false;
	else
	{
		delta = sqrt(delta);
		(*i1) = (-b-delta)/a;
		(*i2) = (-b+delta)/a;
	}
	return true;
}
*/


inline bool SkinSurface::rayQuadricIntersection(double *Q, double *i1, double *i2, double *cache)
{
	double a, b, c;
	
	//a = 0;
	// Alert: orig and dir are currently 3-element arrays ...
	//for (int i=0;i<4;i++)
	//	for (int j=0;j<4;j++)
	//		a += dir[i]*Q[i][j]*dir[j];
	
	//inline using symmetry and caching
	//double q44=Q[3][3],q33=Q[2][2],q22=Q[1][1],q11=Q[0][0],q34=Q[2][3],q24=Q[1][3],q23=Q[1][2],q14=Q[0][3],q13=Q[0][2],q12=Q[0][1];
	double q44=Q[3],q33=Q[2],q22=Q[1],q11=Q[0],q34=Q[6],q24=Q[8],q23=Q[5],q14=Q[9],q13=Q[7],q12=Q[4];
	
	a = cache[0]*q33+cache[1]*q23+cache[2]*q22+cache[3]*q13+cache[4]*q12+cache[5]*q11;
	
	// bug fix (extremely rare)
	if (a == 0)
		return false;
	
	// Alert: orig and dir are currently 3-element arrays ...
	//for (int i=0;i<4;i++)
	//	for (int j=0;j<4;j++)
	//		b += orig[i]*Q[i][j]*dir[j];
	
	b = cache[15]*q34+cache[16]*q33+cache[17]*q24+cache[18]*q23+cache[19]*q22+cache[20]*q14+cache[21]*q13+cache[22]*q12+cache[23]*q11;
	
	//for (int i=0;i<4;i++)
	//	for (int j=0;j<4;j++)
	//		c += orig[i]*Q[i][j]*orig[j];
	
	c = q44+cache[6]*q34+cache[7]*q33+cache[8]*q24+cache[9]*q23+cache[10]*q22+cache[11]*q14+cache[12]*q13+cache[13]*q12+cache[14]*q11;
	
	double delta = b*b-a*c;
	
	if (delta<0)
		return false;
	
	delta = sqrt(delta);
	(*i1) = (-b-delta)/a;
	(*i2) = (-b+delta)/a;
	
	return true;
}


inline bool SkinSurface::rayQuadricIntersectionX(double *Q, double dir_x, double *i1, double *i2, double *cache)
{
	// optimised X-specific algorithm
	double a, b, c;

	double q44=Q[3],q33=Q[2],q22=Q[1],q11=Q[0],q34=Q[6],q24=Q[8],q23=Q[5],q14=Q[9],q13=Q[7],q12=Q[4];

	a = dir_x*dir_x*q11;

	// bug fix (extremely rare)
	if (a == 0)
		return false;

	b = dir_x*q14+cache[9]*q13+cache[10]*q12+cache[11]*q11;

	c = q44+cache[0]*q34+cache[1]*q33+cache[2]*q24+cache[3]*q23+cache[4]*q22+cache[5]*q14+cache[6]*q13+cache[7]*q12+cache[8]*q11;

	double delta = b*b-a*c;

	if (delta<0)
		return false;

	delta = sqrt(delta);
	(*i1) = (-b-delta)/a;
	(*i2) = (-b+delta)/a;

	return true;
}


inline bool SkinSurface::rayQuadricIntersectionY(double *Q, double dir_y, double *i1, double *i2, double *cache)
{
	// optimised Y-specific algorithm
	double a, b, c;

	double q44=Q[3],q33=Q[2],q22=Q[1],q11=Q[0],q34=Q[6],q24=Q[8],q23=Q[5],q14=Q[9],q13=Q[7],q12=Q[4];

	a = dir_y*dir_y*q22;

	// bug fix (extremely rare)
	if (a == 0)
		return false;

	b = dir_y*q24+cache[9]*q23+cache[10]*q22+cache[11]*q12;

	c = q44+cache[0]*q34+cache[1]*q33+cache[2]*q24+cache[3]*q23+cache[4]*q22+cache[5]*q14+cache[6]*q13+cache[7]*q12+cache[8]*q11;

	double delta = b*b-a*c;

	if (delta<0)
		return false;

	delta = sqrt(delta);
	(*i1) = (-b-delta)/a;
	(*i2) = (-b+delta)/a;

	return true;
}


inline bool SkinSurface::rayQuadricIntersectionZ(double *Q, double dir_z, double *i1, double *i2, double *cache)
{
	// optimised Z-specific algorithm
	double a, b, c;

	double q44=Q[3],q33=Q[2],q22=Q[1],q11=Q[0],q34=Q[6],q24=Q[8],q23=Q[5],q14=Q[9],q13=Q[7],q12=Q[4];

	a = dir_z*dir_z*q33;

	// bug fix (extremely rare)
	if (a == 0)
		return false;

	b = dir_z*q34+cache[9]*q33+cache[10]*q23+cache[11]*q13;

	c = q44+cache[0]*q34+cache[1]*q33+cache[2]*q24+cache[3]*q23+cache[4]*q22+cache[5]*q14+cache[6]*q13+cache[7]*q12+cache[8]*q11;

	double delta = b*b-a*c;

	if (delta<0)
		return false;

	delta = sqrt(delta);
	(*i1) = (-b-delta)/a;
	(*i2) = (-b+delta)/a;

	return true;
}


inline bool SkinSurface::rayQuadricIntersectionX(double *Q, double o[3], double t[2])
{
	// optimised X-specific algorithm with dir_x = 1
	double q44=Q[3],q33=Q[2],q22=Q[1],q11=Q[0],q34=Q[6],q24=Q[8],q23=Q[5],q14=Q[9],q13=Q[7],q12=Q[4];

	// bug fix (extremely rare)
	if (q11 == 0)
		return false;

	double b = q14+o[2]*q13+o[1]*q12+o[0]*q11;

	double c = q44+(o[2]+o[2])*q34+(o[2]*o[2])*q33+(o[1]+o[1])*q24+(2*o[1]*o[2])*q23+(o[1]*o[1])*q22+(o[0]+o[0])*q14+(2*o[0]*o[2])*q13+(2*o[0]*o[1])*q12+(o[0]*o[0])*q11;

	double delta = b*b-q11*c;

	if (delta<0)
		return false;

	delta = sqrt(delta);
	t[0] = (-b-delta)/q11;
	t[1] = (-b+delta)/q11;

	return true;
}


inline bool SkinSurface::rayQuadricIntersectionY(double *Q, double o[3], double t[2])
{
	// optimised Y-specific algorithm with dir_y = 1
	double q44=Q[3],q33=Q[2],q22=Q[1],q11=Q[0],q34=Q[6],q24=Q[8],q23=Q[5],q14=Q[9],q13=Q[7],q12=Q[4];

	// bug fix (extremely rare)
	if (q22 == 0)
		return false;

	double b = q24+o[2]*q23+o[1]*q22+o[0]*q12;

	double c = q44+(o[2]+o[2])*q34+(o[2]*o[2])*q33+(o[1]+o[1])*q24+(2*o[1]*o[2])*q23+(o[1]*o[1])*q22+(o[0]+o[0])*q14+(2*o[0]*o[2])*q13+(2*o[0]*o[1])*q12+(o[0]*o[0])*q11;

	double delta = b*b-q22*c;

	if (delta<0)
		return false;

	delta = sqrt(delta);
	t[0] = (-b-delta)/q22;
	t[1] = (-b+delta)/q22;

	return true;
}


inline bool SkinSurface::rayQuadricIntersectionZ(double *Q, double o[3], double t[2])
{
	// optimised Z-specific algorithm with dir_z = 1
	double q44=Q[3],q33=Q[2],q22=Q[1],q11=Q[0],q34=Q[6],q24=Q[8],q23=Q[5],q14=Q[9],q13=Q[7],q12=Q[4];

	// bug fix (extremely rare)
	if (q33 == 0)
		return false;

	double b = q34+o[2]*q33+o[1]*q23+o[0]*q13;

	double c = q44+(o[2]+o[2])*q34+(o[2]*o[2])*q33+(o[1]+o[1])*q24+(2*o[1]*o[2])*q23+(o[1]*o[1])*q22+(o[0]+o[0])*q14+(2*o[0]*o[2])*q13+(2*o[0]*o[1])*q12+(o[0]*o[0])*q11;

	double delta = b*b-q33*c;

	if (delta<0)
		return false;

	delta = sqrt(delta);
	t[0] = (-b-delta)/q33;
	t[1] = (-b+delta)/q33;

	return true;
}


void SkinSurface::getRayIntersection(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id)
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
	int numQuadrics = ind_2d[i2][i1];
	#else */
	int numQuadrics = gridMixedCellMap2D[ i2*n_2d_first + i1 ].size();
	// #endif

	if (numQuadrics == 0) return;

	// this cache exploits the particular structure of skin ray-tracing in order to improve speed.
	// all intersections are with quadrics so the ray invariant part is precomputed once and for all.
	double cache[24];

	double dir[3];

	dir[0] = 0.;
	dir[1] = 0.;
	dir[2] = 0.;
	dir[varying_coord] = pb[varying_coord] - pa[varying_coord];

	// direction and origin are constant
	// double pat[4],dirt[4];
	// ASSIGN(pat,pa);
	// ASSIGN(dirt,dir);

	// load cache with constant-ray info
	// this is the quadric-ray intersection that is constant
	// given the current ray
	// double d1 = dirt[0],d2 = dirt[1],d3 = dirt[2],d4 = dirt[3];
	double d1 = dir[0], d2 = dir[1], d3 = dir[2]; // d4 is dirt[3] = 0 ...
	double o1 = pa[0], o2 = pa[1], o3 = pa[2];

	// cache[0] = 0.; // cache[0] = d4*d4;
	// cache[1] = 0.; // cache[1] = 2*d3*d4;
	// cache[2] = d3*d3;
	// cache[3] = 0.; // cache[3] = 2*d2*d4;
	// cache[4] = 2*d2*d3;
	// cache[5] = d2*d2;
	// cache[6] = 0.; // cache[6] = 2*d1*d4;
	// cache[7] = 2*d1*d3;
	// cache[8] = 2*d1*d2;
	// cache[9] = d1*d1;
	// cache[10] = 1.; // cache[10] = o4*o4;
	// cache[11] = 2*o3; // cache[11] = 2*o3*o4;
	// cache[12] = o3*o3;
	// cache[13] = 2*o2; // cache[13] = 2*o2*o4;
	// cache[14] = 2*o2*o3;
	// cache[15] = o2*o2;
	// cache[16] = 2*o1; // cache[16] = 2*o1*o4;
	// cache[17] = 2*o1*o3;
	// cache[18] = 2*o1*o2;
	// cache[19] = o1*o1;
	// cache[20] = 0.; // cache[20] = d4*o4;
	// cache[21] = d3; // cache[21] = d3*o4+d4*o3;
	// cache[22] = d3*o3;
	// cache[23] = d2; // cache[23] = d2*o4+d4*o2;
	// cache[24] = d2*o3+d3*o2;
	// cache[25] = d2*o2;
	// cache[26] = d1; // cache[26] = d1*o4+d4*o1;
	// cache[27] = d1*o3+d3*o1;
	// cache[28] = d1*o2+d2*o1;
	// cache[29] = d1*o1;

	// More compact cache than in the previous version
	cache[0] = d3*d3;
	cache[1] = 2*d2*d3;
	cache[2] = d2*d2;
	cache[3] = 2*d1*d3;
	cache[4] = 2*d1*d2;
	cache[5] = d1*d1;
	cache[6] = 2*o3;
	cache[7] = o3*o3;
	cache[8] = 2*o2;
	cache[9] = 2*o2*o3;
	cache[10] = o2*o2;
	cache[11] = 2*o1;
	cache[12] = 2*o1*o3;
	cache[13] = 2*o1*o2;
	cache[14] = o1*o1;
	cache[15] = d3;
	cache[16] = d3*o3;
	cache[17] = d2;
	cache[18] = d2*o3+d3*o2;
	cache[19] = d2*o2;
	cache[20] = d1;
	cache[21] = d1*o3+d3*o1;
	cache[22] = d1*o2+d2*o1;
	cache[23] = d1*o1;

	for (int iter = 0; iter<numQuadrics; iter++)
	{
		/* #if !defined(OPTIMIZE_GRIDS)
		int it = GRID_MIXED_CELL_MAP_2D(i1,i2,iter,n_2d_first,n_2d_last);
		#else */
		#if !defined(MULTITHREADED_SKIN_BUILDING)
		int it = gridMixedCellMap2D[i2*n_2d_first+i1][iter];

		MixedCell *mc = mixedComplex[it];
		#else
		pair<int,int> it = gridMixedCellMap2D[i2*n_2d_first+i1][iter];

		MixedCell *mc = thread_data_wrapper[ it.first ].mixedComplex[ it.second ];
		#endif
		// #endif

		// // ray-object culling for an object's bounding rectangle which is not intersected by the ray
		// // compute the bounding box of the object
		// int panel_dims[3][2] = {{1,2}, {0,1}, {0,2}};

		// double down1 = INFINITY;
		// double down2 = INFINITY;

		// double up1 = -INFINITY;
		// double up2 = -INFINITY;

		// int m = panel_dims[panel][0];
		// int n = panel_dims[panel][1];

		// // use the cube of the atom as bounding box object for voronoi cells
		// if (mc->surface_type == DELAUNAY_POINT_CELL)
		// {
		// 	int index = ((Del0Cell*)mc)->id;

		// 	// avoid mapping a non atom (bounding box virtual atom)
		// 	if (index >= delphi->numAtoms)
		// 	{
		// 		continue;
		// 	}
		// 	double radius = delphi->atoms[index]->radius;
		// 	double *sphere_center = delphi->atoms[index]->pos;
		//
		// 	down1 = sphere_center[m]-radius;
		// 	down2 = sphere_center[n]-radius;

		// 	up1 = sphere_center[m]+radius;
		// 	up2 = sphere_center[n]+radius;
		// }
		// else
		// {
		// 	// mixed cell points
		// 	vector<double*> &ll = mc->mc_points;

		// 	for (unsigned int i=0; i<ll.size(); i++)
		// 	{
		// 		down1 = MIN(down1,ll[i][m]);
		// 		down2 = MIN(down2,ll[i][n]);
		//
		// 		up1 = MAX(up1,ll[i][m]);
		// 		up2 = MAX(up2,ll[i][n]);
		// 	}
		// }
		// if (pa[m] < down1 || pa[m] > up1) continue;
		// if (pa[n] < down2 || pa[n] > up2) continue;
		// // end of ray-object culling algorithm
		// // The above cached buffers can eventually be calculated once here
		// // (at the eventual first ray-rectangle intersection), and not above

		double intPoint[3];
		double t[2];

		double *Q = mc->quadric;

		// ray quadric intersection
		bool ff = rayQuadricIntersection(Q,&t[0],&t[1],cache);

		// no intersection
		if (!ff)
			continue;

		double double_normal_sign = (mc->surface_type == DELAUNAY_TETRA_CELL) ? -2. : 2.;

		for (int i=0; i<2; i++)
		{
			intPoint[0] = pa[0];
			intPoint[1] = pa[1];
			intPoint[2] = pa[2];
			intPoint[varying_coord] += t[i] * dir[varying_coord];

			// feasibility test: check if the intersection is inside the mixed complex cell
			if (!isFeasible(mc,intPoint))
				continue;

			// compute normal
			if (computeNormals)
			{
				normalsBuffers[thread_id].push_back(0.);
				normalsBuffers[thread_id].push_back(0.);
				normalsBuffers[thread_id].push_back(0.);

				VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

				normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
				normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
				normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
				intersections.push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i], normal));
			}
			else
			{
				intersections.push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i], nullptr));
			}
		}
	}

	if (intersections.size()>0)
		sort(intersections.begin(), intersections.end(), compKeepIndex);
}


#if defined(REPORT_FAILED_RAYS)
inline void SkinSurface::printRayIntersection(double pa[3], double pb[3])
{
	;
}
#endif


inline void SkinSurface::getRayIntersectionX(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id)
{
	int64_t i1 = (int64_t)rintp((pa[1]-ymin_2d)*scale_2d);
	int64_t i2 = (int64_t)rintp((pa[2]-zmin_2d)*scale_2d);

	/* #if !defined(OPTIMIZE_GRIDS)
	int numQuadrics = ind_2d[i2][i1];
	#else */
	int numQuadrics = gridMixedCellMap2D[ i2*ny_2d + i1 ].size();
	// #endif
	
	if (numQuadrics == 0) return;

	// this cache exploits the particular structure of skin ray-tracing in order to improve speed.
	// all intersections are with quadrics so the ray invariant part is precomputed once and for all.
	double cache[12];
	
	double dir_x = pb[0] - pa[0];
	double o1 = pa[0],o2 = pa[1],o3 = pa[2];

	cache[0] = o3+o3;
	cache[1] = o3*o3;
	cache[2] = o2+o2;
	cache[3] = 2*o2*o3;
	cache[4] = o2*o2;
	cache[5] = o1+o1;
	cache[6] = 2*o1*o3;
	cache[7] = 2*o1*o2;
	cache[8] = o1*o1;
	cache[9] = dir_x*o3;
	cache[10] = dir_x*o2;
	cache[11] = dir_x*o1;

	for (int iter=0; iter<numQuadrics; iter++)
	{
		/* #if !defined(OPTIMIZE_GRIDS)
		int it = GRID_MIXED_CELL_MAP_2D(i1,i2,iter,ny_2d,nz_2d);
		#else */
		#if !defined(MULTITHREADED_SKIN_BUILDING)
		int it = gridMixedCellMap2D[i2*ny_2d+i1][iter];

		MixedCell *mc = mixedComplex[it];
		#else
		pair<int,int> it = gridMixedCellMap2D[i2*ny_2d+i1][iter];

		MixedCell *mc = thread_data_wrapper[ it.first ].mixedComplex[ it.second ];
		#endif
		// #endif
		
		double intPoint[3];
		double t[2];
		
		double *Q = mc->quadric;

		// ray quadric intersection
		bool ff = rayQuadricIntersectionX(Q, dir_x, &t[0], &t[1], cache);
		
		// no intersection
		if (!ff)
			continue;
		
		double double_normal_sign = (mc->surface_type == DELAUNAY_TETRA_CELL) ? -2. : 2.;
		
		for (int i=0; i<2; i++)
		{
			intPoint[0] = pa[0] + t[i] * dir_x;
			intPoint[1] = pa[1];
			intPoint[2] = pa[2];
			
			// feasibility test: check if the intersection is inside the mixed complex cell
			if (!isFeasible(mc,intPoint))
				continue;

			// compute normal
			if (computeNormals)
			{
				normalsBuffers[thread_id].push_back(0.);
				normalsBuffers[thread_id].push_back(0.);
				normalsBuffers[thread_id].push_back(0.);

				VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

				normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
				normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
				normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
				intersections.push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i], normal));
			}
			else
			{
				intersections.push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i], nullptr));
			}
		}
	}
	
	if (intersections.size()>0)
		sort(intersections.begin(), intersections.end(), compKeepIndex);
}


inline void SkinSurface::getRayIntersectionY(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id)
{
	int64_t i1 = (int64_t)rintp((pa[0]-xmin_2d)*scale_2d);
	int64_t i2 = (int64_t)rintp((pa[2]-zmin_2d)*scale_2d);

	/* #if !defined(OPTIMIZE_GRIDS)
	int numQuadrics = ind_2d[i2][i1];
	#else */
	int numQuadrics = gridMixedCellMap2D[ i2*nx_2d + i1 ].size();
	// #endif

	if (numQuadrics == 0) return;

	// this cache exploits the particular structure of skin ray-tracing in order to improve speed.
	// all intersections are with quadrics so the ray invariant part is precomputed once and for all.
	double cache[12];

	double dir_y = pb[1] - pa[1];
	double o1 = pa[0],o2 = pa[1],o3 = pa[2];

	// Minimal cache
	cache[0] = o3+o3;
	cache[1] = o3*o3;
	cache[2] = o2+o2;
	cache[3] = 2*o2*o3;
	cache[4] = o2*o2;
	cache[5] = o1+o1;
	cache[6] = 2*o1*o3;
	cache[7] = 2*o1*o2;
	cache[8] = o1*o1;
	cache[9] = dir_y*o3;
	cache[10] = dir_y*o2;
	cache[11] = dir_y*o1;

	for (int iter=0; iter<numQuadrics; iter++)
	{
		/* #if !defined(OPTIMIZE_GRIDS)
		int it = GRID_MIXED_CELL_MAP_2D(i1,i2,iter,nx_2d,nz_2d);
		#else */
		#if !defined(MULTITHREADED_SKIN_BUILDING)
		int it = gridMixedCellMap2D[i2*nx_2d+i1][iter];

		MixedCell *mc = mixedComplex[it];
		#else
		pair<int,int> it = gridMixedCellMap2D[i2*nx_2d+i1][iter];

		MixedCell *mc = thread_data_wrapper[ it.first ].mixedComplex[ it.second ];
		#endif
		// #endif

		double intPoint[3];
		double t[2];

		double *Q = mc->quadric;

		// ray quadric intersection
		bool ff = rayQuadricIntersectionY(Q, dir_y, &t[0], &t[1], cache);

		// no intersection
		if (!ff)
			continue;

		double double_normal_sign = (mc->surface_type == DELAUNAY_TETRA_CELL) ? -2. : 2.;

		for (int i=0; i<2; i++)
		{
			intPoint[0] = pa[0];
			intPoint[1] = pa[1] + t[i] * dir_y;
			intPoint[2] = pa[2];

			// feasibility test: check if the intersection is inside the mixed complex cell
			if (!isFeasible(mc,intPoint))
				continue;

			// compute normal
			if (computeNormals)
			{
				normalsBuffers[thread_id].push_back(0.);
				normalsBuffers[thread_id].push_back(0.);
				normalsBuffers[thread_id].push_back(0.);

				VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

				normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
				normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
				normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
				intersections.push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i], normal));
			}
			else
			{
				intersections.push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i], nullptr));
			}
		}
	}

	if (intersections.size()>0)
		sort(intersections.begin(), intersections.end(), compKeepIndex);
}


inline void SkinSurface::getRayIntersectionZ(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id)
{
	int64_t i1 = (int64_t)rintp((pa[0]-xmin_2d)*scale_2d);
	int64_t i2 = (int64_t)rintp((pa[1]-ymin_2d)*scale_2d);

	/* #if !defined(OPTIMIZE_GRIDS)
	int numQuadrics = ind_2d[i2][i1];
	#else */
	int numQuadrics = gridMixedCellMap2D[i2*nx_2d+i1].size();
	// #endif

	if (numQuadrics == 0) return;

	// this cache exploits the particular structure of skin ray-tracing in order to improve speed.
	// all intersections are with quadrics so the ray invariant part is precomputed once and for all.
	double cache[12];

	double dir_z = pb[2] - pa[2];
	double o1 = pa[0],o2 = pa[1],o3 = pa[2];

	// More compact cache than in the previous version
	cache[0] = o3+o3;
	cache[1] = o3*o3;
	cache[2] = o2+o2;
	cache[3] = 2*o2*o3;
	cache[4] = o2*o2;
	cache[5] = o1+o1;
	cache[6] = 2*o1*o3;
	cache[7] = 2*o1*o2;
	cache[8] = o1*o1;
	cache[9] = dir_z*o3;
	cache[10] = dir_z*o2;
	cache[11] = dir_z*o1;

	for (int iter=0; iter<numQuadrics; iter++)
	{
		/* #if !defined(OPTIMIZE_GRIDS)
		int it = GRID_MIXED_CELL_MAP_2D(i1,i2,iter,nx_2d,ny_2d);
		#else */
		#if !defined(MULTITHREADED_SKIN_BUILDING)
		int it = gridMixedCellMap2D[i2*nx_2d+i1][iter];

		MixedCell *mc = mixedComplex[it];
		#else
		pair<int,int> it = gridMixedCellMap2D[i2*nx_2d+i1][iter];

		MixedCell *mc = thread_data_wrapper[ it.first ].mixedComplex[ it.second ];
		#endif
		// #endif

		double intPoint[3];
		double t[2];

		double *Q = mc->quadric;

		// ray quadric intersection
		bool ff = rayQuadricIntersectionZ(Q, dir_z, &t[0], &t[1], cache);

		// no intersection
		if (!ff)
			continue;

		double double_normal_sign = (mc->surface_type == DELAUNAY_TETRA_CELL) ? -2. : 2.;

		for (int i=0; i<2; i++)
		{
			intPoint[0] = pa[0];
			intPoint[1] = pa[1];
			intPoint[2] = pa[2] + t[i] * dir_z;

			// feasibility test: check if the intersection is inside the mixed complex cell
			if (!isFeasible(mc,intPoint))
				continue;

			// compute normal
			if (computeNormals)
			{
				normalsBuffers[thread_id].push_back(0.);
				normalsBuffers[thread_id].push_back(0.);
				normalsBuffers[thread_id].push_back(0.);

				VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

				normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
				normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
				normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
				intersections.push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i], normal));
			}
			else
			{
				intersections.push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i], nullptr));
			}
		}
	}

	if (intersections.size()>0)
		sort(intersections.begin(), intersections.end(), compKeepIndex);
}


inline void SkinSurface::getRayIntersection(MixedCell *mc, double pa[3], double ray_dir, int id, int *netIntersections, int thread_id)
{
	double t[2];
	double o1 = pa[0], o2 = pa[1], o3 = pa[2];

	int varying_coord = 0;
	if (panel == 1)
		varying_coord = 2;
	else if (panel == 2)
		varying_coord = 1;

	double cache[12] = {o3+o3,o3*o3,o2+o2,2*o2*o3,o2*o2,o1+o1,2*o1*o3,2*o1*o2,o1*o1,ray_dir*o3,ray_dir*o2,ray_dir*o1};
	double *Q = mc->quadric;

	bool ff;
	if (panel == 0) ff = rayQuadricIntersectionX(Q, ray_dir, &t[0],&t[1], cache);
	if (panel == 1) ff = rayQuadricIntersectionZ(Q, ray_dir, &t[0],&t[1], cache);
	if (panel == 2) ff = rayQuadricIntersectionY(Q, ray_dir, &t[0],&t[1], cache);
	if (!ff) return;

	double double_normal_sign = (mc->surface_type != DELAUNAY_TETRA_CELL) ? -2.0 : 2.0;

	for (int i=0; i<2; i++)
	{
		double intPoint[3] = {o1, o2, o3};
		intPoint[varying_coord] += t[i] * ray_dir;

		if (!isFeasible(mc,intPoint)) continue;

		pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p = &pixel_intersections[id][ num_pixel_intersections[id]++ ];

		VERTEX_TYPE *normal = nullptr;
		if (computeNormals)
		{
			normalsBuffers[thread_id].push_back(0.);
			normalsBuffers[thread_id].push_back(0.);
			normalsBuffers[thread_id].push_back(0.);

			normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

			normal[0] = double_normal_sign*(Q[0]*intPoint[0]+Q[4]*intPoint[1]+Q[7]*intPoint[2]+Q[9]);
			normal[1] = double_normal_sign*(Q[4]*intPoint[0]+Q[1]*intPoint[1]+Q[5]*intPoint[2]+Q[8]);
			normal[2] = double_normal_sign*(Q[7]*intPoint[0]+Q[5]*intPoint[1]+Q[2]*intPoint[2]+Q[6]);
		}
		int_p->first = t[i];
		int_p->second = normal;
		++*netIntersections;
	}
}


bool SkinSurface::getProjection(double p[3], double *proj1, double *proj2, double *proj3,
								double *normal1, double *normal2, double *normal3)
{
	// get the mixed cells that are associated to this grid point
	// by querying the auxiliary grid
	double dist;
	#if !defined(MULTITHREADED_SKIN_BUILDING)
	set<int>cells;
	#else
	set<pair<int,int>>cells;
	#endif

	// move from delphi grid to auxiliary grid
	int irefx = (int)rintp((p[0]-xmin)*scale);
	int irefy = (int)rintp((p[1]-ymin)*scale);
	int irefz = (int)rintp((p[2]-zmin)*scale);

	/* #if !defined(OPTIMIZE_GRIDS)
	for (int i=0; i<ind[irefz][irefy][irefx]; i++)
		cells.insert((GRID_MIXED_CELL_MAP(irefx,irefy,irefz,i,nx,ny,nz)));
	#else */
	for (int i=0; i<gridMixedCellMap[ irefz*ny*nx + irefy*nx + irefx ].size(); i++)
	{
		cells.insert( gridMixedCellMap[ irefz*ny*nx + irefy*nx + irefx ][i] );
	}
	// #endif


	// keep the nearest patch
	double locProj[3], locNorm[3];
	double minDist = INFINITY;

	dist      =  0;
	locProj[0] = 0;
	locProj[1] = 0;
	locProj[2] = 0;

	locNorm[0] = 0;
	locNorm[1] = 0;
	locNorm[2] = 0;

	#if !defined(MULTITHREADED_SKIN_BUILDING)
	for (set<int>::iterator it = cells.begin(); it != cells.end(); it++)
	#else
	for (set<pair<int,int>>::iterator it = cells.begin(); it != cells.end(); it++)
	#endif
	{
		#if !defined(MULTITHREADED_SKIN_BUILDING)
		MixedCell *mc = mixedComplex[*it];
		#else
		MixedCell *mc = thread_data_wrapper[ it->first ].mixedComplex[ it->second ];
		#endif

		if (mc->surface_type == SKIP_CELL)
			continue;

		double *Q = mc->quadric;

		projectToQuadric(p,Q,mc->surface_type,locProj,locNorm,dist);
		
		if (dist == INFINITY)
			continue;

		bool flag = isFeasible(mc,locProj);

		if (!flag)
			continue;

		if (dist<minDist)
		{
			minDist = dist;
			(*proj1) = locProj[0];
			(*proj2) = locProj[1];
			(*proj3) = locProj[2];

			(*normal1) = locNorm[0];
			(*normal2) = locNorm[1];
			(*normal3) = locNorm[2];
		}
	}
	if (cells.size() == 0 || minDist == INFINITY)
	{
		{  
			#ifdef ENABLE_BOOST_THREADS
			boost::mutex::scoped_lock scopedLock(mutex);
			#endif
			(*errorStream) << endl << WARN << "Approximating bgp with grid point";
		}
		(*proj1) = p[0];
		(*proj2) = p[1];
		(*proj3) = p[2];
		return true;
	}
	// further check that the projection is in the cube	
	//if (!testInCube((*proj1),(*proj2),(*proj3),p[0],p[1],p[2],delphi->side))
	//	cout << endl << WARN << "Out of cube projection in SkinSurface::getProjection!";
	
	return true;
}


void SkinSurface::projectToQuadric(double *y,double *Q,int type,double *proj,double *normal,double &dist)
{
	// sphere case performs fast projection
	if (type == 0 || type == 3)
	{
		double cc[3];
		cc[0] = -Q[9]/Q[0];
		cc[1] = -Q[8]/Q[0];
		cc[2] = -Q[6]/Q[0];
		double k = Q[3]/Q[0];
		double r2 = -k+cc[0]*cc[0]+cc[1]*cc[1]+cc[2]*cc[2];

		// immaginary sphere, skip
		if (r2 < 0)
		{
			dist = INFINITY;
			return;
		}
		double r = sqrt(r2);
		
		// existing sphere is projected
		double pp[3];
		double ttt;

		SUB(pp,y,cc);
		NORMALIZE_S(pp,ttt);
		ADD_MUL(proj,cc,pp,r);

		// get the normal depending on the insidness value
		double insideness = 0;
		double yi,yj;

		for (int i=0; i<4; i++)
		{
			if (i == 3)
				yi = 1;
			else
				yi = y[i];

			for (int j=0; j<4; j++)
			{
				if (j == 3)
					yj = 1;
				else
					yj = y[j];

				double qij = 0;
				
				if (i==0 && j==0)
					qij = Q[0];
				else if (i==1 && j==1)
					qij = Q[1];
				else if (i==2 && j==2)
					qij = Q[2];
				else if (i==3 && j==3)
					qij = Q[3];
				if ((i==1 && j==0) || (i==0 && j==1))
					qij = Q[4];
				else if ((i==2 && j==0) || (i==0 && j==2))
					qij = Q[7];
				else if ((i==3 && j==0) || (i==0 && j==3))
					qij = Q[9];
				else if ((i==1 && j==2) || (i==2 && j==1))
					qij = Q[5];
				else if ((i==2 && j==3) || (i==3 && j==2))
					qij = Q[6];
				else if ((i==1 && j==3) || (i==3 && j==1))
					qij = Q[8];
				
				insideness += yi*yj*qij;
			}
		}
     
		if (insideness > 0)
			insideness = 1;
		else
			insideness = -1;
	     
		if (type == 3)
			insideness = -insideness;
			
		DIST(dist,proj,y);
		SUB(normal,y,proj);

		if (insideness < 0)
		{
			normal[0] = -normal[0];
			normal[1] = -normal[1];
			normal[2] = -normal[2];
		}
		double t;
		NORMALIZE(normal,t);
		/*
		printf("\n r2 = %f\n",r2);
		printf("\ny %f %f %f\n",y[0],y[1],y[2]);
		printf("\nc %f %f %f r2 %f\n",cc[0],cc[1],cc[2],r2);
		printf("\nproj %f %f %f \n",proj[0],proj[1],proj[2]);
		printf("\ndist %f r %f\n",dist,r);
	    getchar();
		*/
		return;
	}

	// generic quadric projection

	double b[3];
	double poly[7];

	/*
	 * double c = Q[3][3];
	b[0] = 2*Q[3][0];
	b[1] = 2*Q[3][1];
	b[2] = 2*Q[3][2];

	// compute sixth degree polynomial
	double a11 = Q[0][0];
	double a12 = Q[0][1];
	double a13 = Q[0][2];
	double a22 = Q[1][1];
	double a33 = Q[2][2];
	double a23 = Q[1][2];
	*/

	double c = Q[3];
	b[0] = 2*Q[9];
	b[1] = 2*Q[8];
	b[2] = 2*Q[6];

	// compute sixth degree polynomial
	double a11 = Q[0];
	double a12 = Q[4];
	double a13 = Q[7];
	double a22 = Q[1];
	double a33 = Q[2];
	double a23 = Q[5];

	double a112 = a11*a11;
	double a222 = a22*a22;
	double a332 = a33*a33;

	double a232 = a23*a23;
	double a132 = a13*a13;
	double a122 = a12*a12;

	double b02 = b[0]*b[0];
	double b12 = b[1]*b[1];
	double b22 = b[2]*b[2];

	double y02 = y[0]*y[0];
	double y12 = y[1]*y[1];
	double y22 = y[2]*y[2];

	double a234 = a232*a232;
	double a134 = a132*a132;
	double a123 = a12*a12*a12;
	double a124 = a123*a12;

	double a233 = a232*a23;
	double a133 = a132*a13;

	// t^6 coefficient
	poly[6] = 16* (a11*a22*a33-a122*a33-a11*a232+2*a12*a13*a23-a132*a22)*
					(4*a11*a22*a33*c-4*a122*a33*c-4*a11*a232*c+8*a12*a13*a23*c
					-4*a132*a22*c-b02*a22*a33+2*b[0]*
					b[1]*a12*a33-b12*a11*a33+b02*a232-2*b[0]*b[1]*a13*a23-
					2*b[0]*b[2]*a12*a23+2*b[1]*b[2]*a11*a23+2*b[0]*b[2]*a13*a22
					-b22*a11*a22+b12*
					a132-2*b[1]*b[2]*a12*a13+b22*a122);   

	// t^5 coefficient
	poly[5] = 16*(a22*a33+a11*a33-a232+a11*a22-a132-a122)*
					(4*a11*a22*a33*c-4*a122*a33*c-4*a11*a232*c+8*a12*a13*a23*c
					-4*a132*a22*c-b02*a22*a33+2*b[0]*b[1]*a12*a33-
					b12*a11*a33+b02*a232-2*b[0]*b[1]*a13*a23-
					2*b[0]*b[2]*a12*a23+2*b[1]*b[2]*a11*a23+2*b[0]*b[2]*a13*a22
					-b22*a11*a22+b12*a132-2*b[1]*
					b[2]*a12*a13+b22*a122);

	// t^4 coefficient
	poly[4] = +(((16*a222+64*a11*a22-32*a122+16*a112)*a332+
				((-32*a22-64*a11)*a232+64*a12*a13*a23+64*a11*a222+
				(-64*a132-64*a122+64*a112)*a22-32*a11*a132-64*a11*a122)*a33+16*a234
				+(-64*a11*a22+32*a132+32*a122-32*a112)*a232+
				(64*a12*a13*a22+64*a11*a12*a13)*a23+(16*a112-32*a132)*a222
				+(-64*a11*a132-32*a11*a122)*a22+16*a134+32*a122*a132+16*a124)*c+((16*y02*a11+16*b[0]*y[0])*a222+
				(-16*y02*a122+(-32*y[0]*y[1]*a11-16*b[0]*y[1]-16*y[0]*b[1])*
				a12+16*y12*a112-4*b12-16*b02)*a22+32*y[0]*y[1]*a123+(-16*y12*a11+16*b[1]*y[1]+16*b[0]*y[0])*a122+
				((-16*b[0]*y[1]-16*y[0]*b[1])*a11+24*b[0]*b[1])*a12+16*b[1]*y[1]*a112+(-16*b12-4*b02)*a11)*a332+(
				((-32*y02*a11-32*b[0]*y[0])*a22+16*y02*a122+(32*y[0]*y[1]*a11+16*b[0]*y[1]+
				16*y[0]*b[1])*a12-16*y12*a112+4*b12+16*b02)*a232+(
				((32*y02*a12+32*y[0]*y[1]*a11+16*b[0]*y[1]+16*y[0]*b[1])*a13+
				(32*y[0]*y[2]*a11+16*b[0]*y[2]+16*y[0]*b[2])*a12-32*y[1]*y[2]*a112+8*b[1]*b[2])*a22+
				(-96*y[0]*y[1]*a122+(32*y12*a11-32*b[1]*y[1]-32*b[0]*y[0])*a12+(16*b[0]*y[1]+
				16*y[0]*b[1])*a11-24*b[0]*b[1])*a13-32*y[0]*y[2]*a123+
				(32*y[1]*y[2]*a11-16*b[1]*y[2]-16*y[1]*b[2])*a122+((16*b[0]*y[2]+16*y[0]*b[2])*a11-
				24*b[0]*b[2])*a12+(-16*b[1]*y[2]-16*y[1]*b[2])*a112+32*b[1]*b[2]*a11)*a23+
				(-16*y02*a132+(-32*y[0]*y[2]*a11-16*b[0]*y[2]-16*y[0]*b[2])*a13+16*y22*a112
				-4*b22-16*b02)*a222+((32*y[0]*y[1]*a12-32*y12*a11)*a132+
				(32*y[0]*y[2]*a122+(32*y[1]*y[2]*a11+16*b[1]*y[2]+16*y[1]*b[2])*a12+32*b[0]*b[2])*a13-
				32*y22*a11*a122+32*b[0]*b[1]*a12+(-16*b22-16*b12-16*b02)*a11)*a22+
				(16*y12*a122+(16*b[0]*y[1]+16*y[0]*b[1])*a12-32*b[1]*y[1]*a11+16*b12+4*b02)*a132+
				(-32*y[1]*y[2]*a123+(-16*b[0]*y[2]-16*y[0]*b[2])*a122+((16*b[1]*y[2]+
				16*y[1]*b[2])*a11-24*b[1]*b[2])*a12+8*b[0]*b[2]*a11)*a13+16*y22*a124+8*b22*a122+32*b[0]*b[1]*a11* 
				a12+(-4*b22-16*b12)*a112)*a33+(16*y02*a11+16*b[0]*y[0])*a234+
				((-32*y02*a12-32*y[0]*y[1]*a11-16*b[0]*y[1]-16*y[0]*b[1])*a13+(-32*y[0]*y[2]*a11-
				16*b[0]*y[2]-16*y[0]*b[2])*a12+32*y[1]*y[2]*a112-8*b[1]*b[2])*a233+(
				(16*y02*a132+(32*y[0]*y[2]*a11+16*b[0]*y[2]+16*y[0]*b[2])*a13-16*y22*a112+4*b22+
				16*b02)*a22+(64*y[0]*y[1]*a12+16*y12*a11+16*b[1]*y[1]+16*b[0]*y[0])*a132+
				(64*y[0]*y[2]*a122+(-96*y[1]*y[2]*a11+16*b[1]*y[2]+16*y[1]*b[2])*a12+(-16*b[0]*y[2]-
				16*y[0]*b[2])*a11-8*b[0]*b[2])*a13+(16*y22*a11+16*b[2]*y[2]+16*b[0]*y[0])*a122+
				((-16*b[0]*y[1]-16*y[0]*b[1])*a11-8*b[0]*b[1])*a12+(16*b[2]*y[2]+16*b[1]*y[1])*a112+
				8*b02*a11)*a232+((-32*y[0]*y[1]*a133+
				(-96*y[0]*y[2]*a12+32*y[1]*y[2]*a11-16*b[1]*y[2]-16*y[1]*b[2])*a132+(
				(32*y22*a11-32*b[2]*y[2]-32*b[0]*y[0])*a12+(16*b[0]*y[1]+16*y[0]*b[1])*a11-24*b[0]*b[1])*a13+
				((16*b[0]*y[2]+16*y[0]*b[2])*a11-24*b[0]*b[2])*a12+(-16*b[1]*y[2]-16*y[1]*b[2])*a112+32*b[1]*b[2]*a11)*a22
				+(-32*y12*a12-16*b[0]*y[1]-16*y[0]*b[1])*a133+
				(64*y[1]*y[2]*a122+(16*b[0]*y[2]+16*y[0]*b[2])*a12+(16*b[1]*y[2]+16*y[1]*b[2])*a11-8*b[1]*b[2])*a132+
				(-32*y22*a123+(16*b[0]*y[1]+16*y[0]*b[1])*a122+((-32*b[2]*y[2]-32*b[1]*y[1])*a11+8*b22+8*b12+
				8*b02)*a12-24*b[0]*b[1]*a11)*a13+(-16*b[0]*y[2]-16*y[0]*b[2])*a123+
				((16*b[1]*y[2]+16*y[1]*b[2])*a11-8*b[1]*b[2])*a122-24*b[0]*b[2]*a11*a12+24*b[1]*b[2]*a112)*a23+
				(32*y[0]*y[2]*a133+(-16*y22*a11+16*b[2]*y[2]+16*b[0]*y[0])*a132+((-16*b[0]*y[2]-
				16*y[0]*b[2])*a11+24*b[0]*b[2])*a13+16*b[2]*y[2]*a112+(-16*b22-4*b02)*a11)*a222+(16*
				y12*a134-32*y[1]*y[2]*a12*a133+(16*y22*a122+(-16*b[0]*y[1]-16*y[0]*b[1])*a12+8*b12)*a132+
				((16*b[0]*y[2]+16*y[0]*b[2])*a122+((16*b[1]*y[2]+16*y[1]*b[2])*a11-24*b[1]*b[2])*a12+32*b[0]*b[2]*a11)*a13
				+(-32*b[2]*y[2]*a11+16*b22+4*b02)*a122+8*b[0]*b[1]*a11*a12+
				(-16*b22-4*b12)*a112)*a22+16*b[1]*y[1]*a134+((-16*b[1]*y[2]-16*y[1]*b[2])*a12-8*b[0]*b[2])*a133
				+((16*b[2]*y[2]+16*b[1]*y[1])*a122-8*b[0]*b[1]*a12+(4*b22+16*b12)*a11)*
				a132+((-16*b[1]*y[2]-16*y[1]*b[2])*a123-8*b[0]*b[2]*a122-24*b[1]*b[2]*a11*a12)*a13+16*b[2]*y[2]*a124
				-8*b[0]*b[1]*a123+(16*b22+4*b12)*a11*a122);

	// t^3 coefficient
	poly[3] = +(((16*a22+16*a11)*
				a332+(-16*a232+16*a222+64*a11*a22-16*a132-32*a122+16*a112)*a33
				+(-16*a22-32*a11)*a232+32*a12*a13*a23+16*a11*a222+(-32*a132-16*a122+16*a112)*a22-16*a11*a132-16*a11*
				a122)*c+(((16*y12+16*y02)*a11+16*b[0]*y[0])*a22+(-16*y12-16*y02)*a122+(-16*b[0]*y[1]
				-16*y[0]*b[1])*a12+16*b[1]*y[1]*a11-4*b12-4*b02)*a332+(
				((-16*y12-16*y02)*a11-16*b[0]*y[0])*a232+(-32*y[1]*y[2]*a11*a22+
				((32*y12+32*y02)*a12+16*b[0]*y[1]+16*y[0]*b[1])*a13+32*y[1]*y[2]*a122+(16*b[0]*y[2]+16*y[0]*b[2])*a12
				+(-16*b[1]*y[2]-16*y[1]*b[2])*a11+8*b[1]*b[2])*a23+((16*y22+16*y02)*a11+16*b[0]*y[0])*a222+
				((-16*y12-16*y02)*a132+(-32*y[0]*y[2]*a11-16*b[0]*y[2]-16*y[0]*b[2])*a13+
				(-16*y22-16*y02)*a122+(-32*y[0]*y[1]*a11-16*b[0]*y[1]-16*y[0]*b[1])*a12+(16*y22+16*y12)*a112
				-4*b22-4*b12-16*b02)*a22-16*b[1]*y[1]*a132+
				(32*y[0]*y[2]*a122+(16*b[1]*y[2]+16*y[1]*b[2])*a12+8*b[0]*b[2])*a13+32*y[0]*y[1]*a123+
				((-16*y22-16*y12)*a11+16*b[1]*y[1]+16*b[0]*y[0])*a122+
				((-16*b[0]*y[1]-16*y[0]*b[1])*a11+24*b[0]*b[1])*a12+16*b[1]*y[1]*a112+(-4*b22-
				16*b12-4*b02)*a11)*a33+32*y[1]*y[2]*a11*a233+
				(((-16*y22-16*y02)*a11-16*b[0]*y[0])*a22+(32*y[0]*y[2]*a11-64*y[1]*y[2]*a12)*a13+32*y[0]*y[1]*a11*a12
				+(-16*y22-16*y12)*a112+(16*b[2]*y[2]+16*b[1]*y[1])*a11+8*b02)* 
				a232+((32*y[1]*y[2]*a132+((32*y22+32*y02)*a12+16*b[0]*y[1]+16*y[0]*b[1])*a13+(16*b[0]*y[2]
				+16*y[0]*b[2])*a12+(-16*b[1]*y[2]-16*y[1]*b[2])*a11+8*b[1]*b[2])*a22-64*y[0]* 
				y[2]*a12*a132+(-64*y[0]*y[1]*a122+((32*y22+32*y12)*a11-32*b[2]*y[2]-32*b[1]*y[1]
				-32*b[0]*y[0])*a12+(16*b[0]*y[1]+16*y[0]*b[1])*a11-16*b[0]*b[1])*a13+
				((16*b[0]*y[2]+16*y[0]*b[2])*a11-16*b[0]*b[2])*a12+(-16*b[1]*y[2]-16*y[1]*b[2])*a112+24*b[1]*b[2]*a11)*a23+
				((-16*y22-16*y02)*a132+(-16*b[0]*y[2]-16*y[0]*b[2])*a13+16*b[2]*y[2]*a11-4*b22-4*b02)*a222+(32*y[0]*y[2]*a133+
				(32*y[0]*y[1]*a12+(-16*y22-16*y12)*a11+16*b[2]*y[2]+16*b[0]*y[0])*a132+((16*b[1]*y[2]+16*y[1]*b[2])*a12
				+(-16*b[0]*y[2]-16*y[0]*b[2])*a11+24*b[0]*b[2])*a13-16*b[2]*y[2]*
				a122+8*b[0]*b[1]*a12+16*b[2]*y[2]*a112+(-16*b22-4*b12-4*b02)*a11)*a22+(8*b12-16*b[1]*y[1]*a11)*a132+
				(((16*b[1]*y[2]+16*y[1]*b[2])*a11-16*b[1]*b[2])*a12+8*b[0]*b[2]*a11)*a13+(8*b22-16*b[2]*y[2]*a11)*a122
				+8*b[0]*b[1]*a11*a12+(-4*b22-4*b12)*a112);

	// t^2 coefficient
	poly[2] =	+((4*a332+(16*a22+16*a11)*a33-8*a232+4*a222+16*a11*a22-8*a132-8*a122+4*a112)*c
				+(4*y12*a22+8*y[0]*y[1]*a12+4*y02*a11+4*b[1]*y[1]+4*b[0]*y[0])*a332+(-4*y12*a232+
				(-8*y[1]*y[2]*a22-8*y[0]*y[1]*a13-8*y[0]*y[2]*a12-4*b[1]*y[2]-4*y[1]*b[2])*a23+4*y22*a222+
				((16*y22+16*y12+16*y02)*a11+16*b[0]*y[0])*a22-4*y02*a132+
				(-8*y[1]*y[2]*a12-8*y[0]*y[2]*a11-4*b[0]*y[2]-4*y[0]*b[2])*a13+(-8*y22-16*y12-16*y02)*a122
				+(-16*b[0]*y[1]-16*y[0]*b[1])*a12+4*y22*a112+16*b[1]*y[1]*a11-b22-4*b12 
				-4*b02)*a33+8*y[1]*y[2]*a233+(-4*y22*a22+8*y[0]*y[2]*a13+8*y[0]*y[1]*a12+(-16*y22-16*y12-
				8*y02)*a11+4*b[2]*y[2]+4*b[1]*y[1]-8*b[0]*y[0])*a232+(
				(-8*y[0]*y[1]*a13-8*y[0]*y[2]*a12-4*b[1]*y[2]-4*y[1]*b[2])*a22+8*y[1]*y[2]*a132+
				((24*y22+24*y12+24*y02)*a12-8*y[0]*y[1]*a11+12*b[0]*y[1]+12*y[0]*b[1])*a13+8*y[1]*
				y[2]*a122+(-8*y[0]*y[2]*a11+12*b[0]*y[2]+12*y[0]*b[2])*a12+8*y[1]*y[2]*a112
				+(-16*b[1]*y[2]-16*y[1]*b[2])*a11+6*b[1]*b[2])*a23+
				(8*y[0]*y[2]*a13+4*y02*a11+4*b[2]*y[2]+4*b[0]*y[0])*a222+((-16*y22-8*y12-16*y02)*a132
				+(-8*y[1]*y[2]*a12-16*b[0]*y[2]-16*y[0]*b[2])*a13-4*y02*a122+
				(-8*y[0]*y[1]*a11-4*b[0]*y[1]-4*y[0]*b[1])*a12+4*y12*a112+16*b[2]*y[2]*a11-4*b22-b12-4*b02)*a22+8*y[0]*y[2]*a133+
				(8*y[0]*y[1]*a12-4*y22*a11+4*b[2]*y[2]-8*b[1]*y[1]+4*b[0]*y[0])*a132+
				(8*y[0]*y[2]*a122+(-8*y[1]*y[2]*a11+12*b[1]*y[2]+12*y[1]*b[2])*a12+(-4*b[0]*y[2]
				-4*y[0]*b[2])*a11+6*b[0]*b[2])*a13+8*y[0]*y[1]*a123+
				(-4*y12*a11-8*b[2]*y[2]+4*b[1]*y[1]+4*b[0]*y[0])*a122+((-4*b[0]*y[1]-4*y[0]*b[1])*a11+6*b[0]*b[1])*a12+
				(4*b[2]*y[2]+4*b[1]*y[1])*a112+(-4*b22-4*b12-b02)*a11);

	// t coefficient
	poly[1] = ((4*a33+4*a22+4*a11)*c+((4*y22+4*y12)*a22+8*y[0]*y[1]*a12+(4*y22+4*y02)*a11+
				4*b[1]*y[1]+4*b[0]*y[0])*a33+(-4*y22-4*y12)*a232+
				(-8*y[0]*y[1]*a13-8*y[0]*y[2]*a12+8*y[1]*y[2]*a11-4*b[1]*y[2]-4*y[1]*b[2])*a23
				+(8*y[0]*y[2]*a13+(4*y12+4*y02)*a11+4*b[2]*y[2]+4*b[0]*y[0])*a22+(-4*y22-4*y02)*a132+
				(-8*y[1]*y[2]*a12-4*b[0]*y[2]-4*y[0]*b[2])*a13+(-4*y12-4*y02)*a122+(-4*b[0]*y[1]
				-4*y[0]*b[1])*a12+(4*b[2]*y[2]+4*b[1]*y[1])*a11-b22-b12-b02);

	// t^0 coefficient
	poly[0] = +c+y22*a33+2*y[1]*y[2]*a23+y12*a22+2*y[0]*y[2]*a13+2*y[0]*y[1]*a12+y02*a11+b[2]*y[2]+b[1]*y[1]+b[0]*y[0];

	double roots[6];
	int num = 6;
	
	if (fastProjection)
		getRealRootsSturm(poly,6,roots,num);
	else
		getRealRootsCompanion(poly,6,roots,num);

	// skip non properly projected patch (typical case of distant patch of no interest)
	if (num == 0)
	{
		// printf("\n %lf %lf %lf %lf %lf %lf %lf",poly[6],poly[5],poly[4],poly[3],poly[2],poly[1],poly[0]);
		dist = INFINITY;
		// getchar();
		return;
	}

	double projt[6][3];
	int opt = 7;
	double minDist = INFINITY;

	// get solution for each root and choose the nearest point
	for (int i=0; i<num; i++)
	{
		double jj = y[0]-roots[i]*b[0];
		double kk = y[1]-roots[i]*b[1];
		double ll = y[2]-roots[i]*b[2];

		/*
		double aa = Q[0][0]*(2*roots[i])+1;
		double ee = Q[1][1]*(2*roots[i])+1;
		double ii = Q[2][2]*(2*roots[i])+1;
		
		double dd = Q[0][1]*(2*roots[i]);
		double gg = Q[0][2]*(2*roots[i]);
		double ff = Q[1][2]*(2*roots[i]);
		*/

		// solve a 3x3 linear system by Cramer rule
		// perform all inline and exploit symmetry		
		double aa = Q[0]*(2*roots[i])+1;
		double ee = Q[1]*(2*roots[i])+1;
		double ii = Q[2]*(2*roots[i])+1;
		
		double dd = Q[4]*(2*roots[i]);
		double gg = Q[7]*(2*roots[i]);
		double ff = Q[5]*(2*roots[i]);

		double ff2 = ff*ff;
		double fg = ff*gg;

		double den = aa*ee*ii+2*dd*fg-gg*gg*ee-dd*dd*ii-aa*ff2;
		double num1 = jj*ee*ii+dd*ff*ll+kk*fg-gg*ee*ll-dd*kk*ii-jj*ff2;
		double num2 = aa*kk*ii+jj*fg+gg*dd*ll-gg*gg*kk-jj*dd*ii-aa*ff*ll;
		double num3 = aa*ee*ll+dd*kk*gg+jj*dd*ff-jj*ee*gg-dd*dd*ll-aa*kk*ff;
		
		if (den != 0.)
		{
			projt[i][0] = num1/den;
			projt[i][1] = num2/den;
			projt[i][2] = num3/den;
		}
		else
		{
			projt[i][0] = num1/(den+1e-20);
			projt[i][1] = num2/(den+1e-20);
			projt[i][2] = num3/(den+1e-20);
		}
		// printf("\n %f %f %f\n",projt[i][0],projt[i][1],projt[i][2]);
		double a = (projt[i][0]-y[0]);
		double b = (projt[i][1]-y[1]);
		double c = (projt[i][2]-y[2]);
		double distl = a*a+b*b+c*c;

		if (distl < minDist)
		{
			minDist = distl;
			opt = i;
		}
	}
	proj[0] = projt[opt][0];
	proj[1] = projt[opt][1];
	proj[2] = projt[opt][2];
	dist = sqrt(minDist);
	// printf("\nopt %f %f %f\n",proj[0],proj[1],proj[2]);

	// get the normal depending on the insidness value
	double insideness = 0;
	double yi,yj;

	for (int i=0; i<4; i++)
	{
		if (i == 3)
			yi = 1;
		else
			yi = y[i];

		for (int j=0; j<4; j++)
		{
			if (j == 3)
				yj = 1;
			else
				yj = y[j];

			double qij = 0;
			
			if (i==0 && j==0)
				qij = Q[0];
			else if (i==1 && j==1)
				qij = Q[1];
			else if (i==2 && j==2)
				qij = Q[2];
			else if (i==3 && j==3)
				qij = Q[3];
			if ((i==1 && j==0) || (i==0 && j==1))
				qij = Q[4];
			else if ((i==2 && j==0) || (i==0 && j==2))
				qij = Q[7];
			else if ((i==3 && j==0) || (i==0 && j==3))
				qij = Q[9];
			else if ((i==1 && j==2) || (i==2 && j==1))
				qij = Q[5];
			else if ((i==2 && j==3) || (i==3 && j==2))
				qij = Q[6];
			else if ((i==1 && j==3) || (i==3 && j==1))
				qij = Q[8];
			
			//insideness += yi*yj*Q[i][j];
			insideness += yi*yj*qij;
		}
	}
	// check the correctness of the projection
	double test = 0;
	for (int i=0; i<4; i++)
	{
		if (i==3)
			yi = 1;
		else
			yi = proj[i];

		for (int j=0; j<4; j++)
		{
			if (j==3)
				yj = 1;
			else
				yj = proj[j];

			double qij = 0;
			
			if (i==0 && j==0)
				qij = Q[0];
			else if (i==1 && j==1)
				qij = Q[1];
			else if (i==2 && j==2)
				qij = Q[2];
			else if (i==3 && j==3)
				qij = Q[3];
			if ((i==1 && j==0) || (i==0 && j==1))
				qij = Q[4];
			else if ((i==2 && j==0) || (i==0 && j==2))
				qij = Q[7];
			else if ((i==3 && j==0) || (i==0 && j==3))
				qij = Q[9];
			else if ((i==1 && j==2) || (i==2 && j==1))
				qij = Q[5];
			else if ((i==2 && j==3) || (i==3 && j==2))
				qij = Q[6];
			else if ((i==1 && j==3) || (i==3 && j==1))
				qij = Q[8];

			//test += yi*yj*Q[i][j];
			test += yi*yj*qij;
		}
	}	
	if (fabs(test) > 1e-3)
	{
		//Projected point does not belong to quadric
		//Unreliable projection detected is discarded (probably a distant patch)
		dist = INFINITY;
		return;
	}

	// concave sphere
	if (type == 3)
		insideness = -insideness;

	SUB(normal,y,proj);

	if (insideness < 0)
	{
		normal[0] = -normal[0];
		normal[1] = -normal[1];
		normal[2] = -normal[2];
	}
	double t;
	NORMALIZE(normal,t);
}


/** save in .skin format*/
bool SkinSurface::save(char *fileName)
{
	ofstream fout;
	fout.open(fileName,ios::out);
	
	cout << endl << INFO << "Writing skin in .skin file format in " << fileName << "...";
	
	if (fout.fail())
	{
		cout << endl << WARN << "Cannot write file " << fileName;
		return false;
	}
	
	if (getNumPatches() == 0)
	{
		cout << endl << WARN << "Cannot write null mesh!";
		return false;
	}
	
	// save TODO
	return true;
}


bool SkinSurface::load(char *fileName)
{
	int len = (int)strlen(fileName);
	if (len == 0)
	{
		cout << WARN << "Cannot load with empty file name!";
		return false;
	}
	
	//cout << endl << INFO << "Loading Skin Surface in file " << fileName << "...";
	
	// TODO LOAD
	return true;
}


void SkinSurface::printSummary()
{
	cout << endl << INFO << "Shrinking value " << getShrinking();

	if (numMixedCells == 0)
	{
		cout << endl << WARN << "Skin surface not loaded!";
	}
	else
	{
		numSkippedCells = 0;

		for (int it = 0; it < numMixedCells; it++)
		{
			if (mixedComplex[it]->surface_type == SKIP_CELL)
				++numSkippedCells;
		}
		cout << endl << INFO << "Number of mixed cells -> " << numMixedCells - numSkippedCells;
		cout << endl << INFO << "Number of del_point/vor_cell -> " << type[DELAUNAY_POINT_CELL];
		cout << endl << INFO << "Number of del_edge/vor_facet -> " << type[DELAUNAY_EDGE_CELL];
		cout << endl << INFO << "Number of del_facet/vor_edge -> " << type[DELAUNAY_FACET_CELL];
		cout << endl << INFO << "Number of del_cell/vor_point -> " << type[DELAUNAY_TETRA_CELL];
		/*
		if (internals != NULL)
		{
			*internals << endl << "mixedcells " << numMixedCells - numSkippedCells;
			*internals << endl << "del_point " << type[0];
			*internals << endl << "del_edge " << type[1];
			*internals << endl << "del_facet " << type[2];
			*internals << endl << "del_cell " << type[3];
		}
		*/
	}
}
