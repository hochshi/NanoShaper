
#include "ConnollySurface.h"


void ConnollySurface::clear()
{
	// remove 2d grid for ray casting
	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridConnollyCellMap != NULL)
		#if !defined(MULTITHREADED_SES_BUILDING)
		deleteVector<int>(gridConnollyCellMap);
		#else
		deleteVector<pair<int,int>>(gridConnollyCellMap);
		#endif

	if (ind != NULL)
		deleteMatrix3D<int>(nx,ny,nz,ind);
	
	if (gridConnollyCellMap2D != NULL)
		deleteVector<int>(gridConnollyCellMap2D);

	if (ind_2d != NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,last_cols_ind,ind_2d);
	#else */
	if (gridConnollyCellMap != NULL)
	{
		for (int64_t i=0; i < nx*ny*nz; i++)
		{
			gridConnollyCellMap[i].clear();
		}
		delete[] gridConnollyCellMap;
		gridConnollyCellMap = NULL;
	}
	if (gridConnollyCellMap2D != NULL)
	{
		for (int64_t i=0; i < last_rows_ind*last_cols_ind; i++)
		{
			gridConnollyCellMap2D[i].clear();
		}
		delete[] gridConnollyCellMap2D;
		gridConnollyCellMap2D = NULL;
	}
	// #endif

	if (x != NULL)
		deleteVector<double>(x);
	if (y != NULL)
		deleteVector<double>(y);
	if (z != NULL)
		deleteVector<double>(z);


	#if !defined(MULTITHREADED_SES_BUILDING)

	#if !defined(NEW_ATOM_PATCHES)
	if (atomPatches != NULL)
	{
		deleteVector<PointCell*>(atomPatches);
	}
	#else
	deleteVector<int>(atomPatches);
	#endif

	#else // MULTITHREADED_SES_BUILDING

	// if (thread_data_wrapper != NULL)
	{
		for (int task_id=0; task_id<numTasks; task_id++)
		{
			for (int thread_id=0; thread_id<numThreadDataWrappersPerTask[task_id]; thread_id++)
			{
				ThreadDataWrapper *tdw = &thread_data_wrapper[task_id*numThreadDataWrappers+thread_id];

				#if !defined(NEW_ATOM_PATCHES)
				deleteVector<PointCell*>(tdw->atomPatches);
				#else
				deleteVector<int>(tdw->atomPatches);
				#endif

				for (int i=0; i<tdw->sesComplex.size(); i++)
				{
					delete tdw->sesComplex[i];
					tdw->sesComplex[i] = NULL;
				}
				tdw->sesComplex.clear();
			}
		}
		// delete thread_data_wrapper;
		// thread_data_wrapper = NULL;
	}

	#endif // MULTITHREADED_SES_BUILDING

	#if defined(CHECK_BUILDUP_DIFF)
	#if !defined(NEW_ATOM_PATCHES)
	if (st_atomPatches != NULL)
	{
		deleteVector<PointCell*>(st_atomPatches);
	}
	#else
	deleteVector<int>(st_atomPatches);
	#endif
	#endif // CHECK_BUILDUP_DIFF

	sesComplex.clear();

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
		for (int i=0; i < sesComplex.size()*numPanels; i++)
		{
			temp_intersections_buffer[i].clear();
		}
		delete[] temp_intersections_buffer;
		temp_intersections_buffer = NULL;

		for (int i=0; i < sesComplex.size()*numPanels; i++)
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


ConnollySurface::~ConnollySurface()
{
	clear();
}


void ConnollySurface::init()
{
	gridConnollyCellMap = NULL;
	gridConnollyCellMap2D = NULL;
	#if !defined(MULTITHREADED_SES_BUILDING)
	atomPatches = NULL;
	#else
	// thread_data_wrapper = NULL;
	#endif
	sesComplex.clear();
	x = NULL;
	y = x;
	z = y;
	savePovRay = false;
	// default grid settings
	AUX_GRID_DIM_CONNOLLY = 100;
	AUX_GRID_DIM_CONNOLLY_2D = 50;
	// max connolly cells for each acc grid cube
	MAX_CONNOLLY_CELLS = 400;
	MAX_CONNOLLY_CELLS_2D = (400*AUX_GRID_DIM_CONNOLLY_2D);
	surfType = MOLECULAR_SURFACE;
	providesAnalyticalNormals = true;
	// MAX_PROBES = 100;
	si_perfil = 1.5;
	
	if (patchBasedAlgorithm)
	{
		num_pixel_intersections = NULL;
		pixel_intersections = NULL;
	}
}


void ConnollySurface::init(ConfigFile *cf)
{
	unsigned int maxSESDim2D = cf->read<unsigned int>( "Max_ses_patches_auxiliary_grid_2d_size", 100 );
	unsigned int maxSESPatches2D = cf->read<unsigned int>( "Max_ses_patches_per_auxiliary_grid_2d_cell", 400 );
	unsigned int maxSESDim = cf->read<unsigned int>( "Max_ses_patches_auxiliary_grid_size", 50 );
	unsigned int maxSESPatches = cf->read<unsigned int>( "Max_ses_patches_per_auxiliary_grid_cell", 400 );
	bool savePovRay = cf->read<bool>( "Save_PovRay", false );
	// int mp = cf->read<int>( "Max_Probes_Self_Intersections",200);
	double si_perfil = cf->read<double>( "Self_Intersections_Grid_Coefficient", 1.5);

	setAuxGrid(maxSESDim,maxSESPatches);
	setAuxGrid2D(maxSESDim2D,maxSESPatches2D);
	inside = 5;
	// setMaxProbes(mp);
	setSIPerfil(si_perfil);
	setSavePovRay(savePovRay);
}


ConnollySurface::ConnollySurface():Surface()
{
	init();
}


ConnollySurface::ConnollySurface(DelPhiShared *ds):Surface()
{
	init();
	// set environment
	delphi = ds;
}


ConnollySurface::ConnollySurface(ConfigFile *cf, DelPhiShared *ds):Surface(cf)
{
	init();
	// set environment
	delphi = ds;
	init(cf);
}


int ConnollySurface::getNumPatches (void)
{
	return sesComplex.size();
}


bool ConnollySurface::isPatchBasedRayTracingSupported (void)
{
	return true;
}


bool ConnollySurface::build()
{
	bool flag;
	#ifdef ENABLE_CGAL 
		flag = buildConnollyCGAL();
	#else
		cout << endl << WARN << "Connolly surface requires CGAL lib";
		return false;
	#endif

    if (!flag)
	{
		cout << endl << ERR << "Connolly surface construction failed";
		return flag;
	}
	return flag;
}


void ConnollySurface::postRayCasting()
{
	// remove 2d grid for ray casting
	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridConnollyCellMap2D != NULL)
		#if !defined(MULTITHREADED_SES_BUILDING)
		deleteVector<int>(gridConnollyCellMap2D);
		#else
		deleteVector<pair<int,int>>(gridConnollyCellMap2D);
		#endif

	if (ind_2d != NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,last_cols_ind,ind_2d);
	#else */
	if (gridConnollyCellMap != NULL)
	{
		for (int64_t i=0; i < nx*ny*nz; i++)
		{
			gridConnollyCellMap[i].clear();
		}
		delete[] gridConnollyCellMap;
		gridConnollyCellMap = NULL;
	}
	if (gridConnollyCellMap2D != NULL)
	{
		for (int64_t i=0; i < last_rows_ind*last_cols_ind; i++)
		{
			gridConnollyCellMap2D[i].clear();
		}
		delete[] gridConnollyCellMap2D;
		gridConnollyCellMap2D = NULL;
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
		for (int i=0; i < sesComplex.size()*numPanels; i++)
		{
			temp_intersections_buffer[i].clear();
		}
		delete[] temp_intersections_buffer;
		temp_intersections_buffer = NULL;

		for (int i=0; i < sesComplex.size()*numPanels; i++)
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


bool ConnollySurface::preBoundaryProjection()
{
	bool flag;
	// useful to project bgp
	if (projBGP)
	{
		flag = buildAuxiliaryGrid();
		if (!flag)
		{
			cout << endl << ERR << "Cannot build 3D auxiliary grid";
			return flag; 
		}
	}
	return true;
}


#ifdef MULTITHREADED_SES_BUILDING

void ConnollySurface::BuildWeightedPoints (TaggedDataWrapper *tdw, double slab_pars[], int thread_pars[])
{
	double x, y, z, r2;

	double halo_layer_size = slab_pars[0];
	double inv_bin_slab_size = slab_pars[1];
	double inv_minibin_slab_size = slab_pars[2];

	for (int i = thread_pars[0]; i < thread_pars[1]; i += thread_pars[2])
	{
		x = delphi->atoms[i].pos[0];
		y = delphi->atoms[i].pos[1];
		z = delphi->atoms[i].pos[2];
		r2 = delphi->atoms[i].radius;
		r2 += probe_radius;
		r2 *= r2;

		int bin_id = (int)((z - (delphi->zmin - delphi->hside)) * inv_bin_slab_size);
		int task_slab_id = tdw->task_in_bin[ bin_id ];

		int minibin_id = (int)((y - (delphi->ymin - delphi->hside)) * inv_minibin_slab_size);

		int involved_task[3];
		int num_involved_tasks = 1;

		involved_task[ 0 ] = task_slab_id;

		// If the reference task is not the last one, check if the
		// above task is also involved because of the halo space
		if (task_slab_id < numTasks - 1)
		{
			// check if the atom is in the halo space of the above task
			if (z >= tdw->slab_max_z[ task_slab_id ] - halo_layer_size)
				involved_task[ num_involved_tasks++ ] = task_slab_id + 1;
		}
		// If the reference task is not the first one, check if the
		// below task is also involved because of the halo space
		if (task_slab_id > 0)
		{
			// check if the atom is in the halo space of the below task
			if (z <= tdw->slab_max_z[ task_slab_id-1 ] + halo_layer_size)
				involved_task[ num_involved_tasks++ ] = task_slab_id - 1;
		}

		for (int k=0; k < num_involved_tasks; k++)
		{
			task_slab_id = involved_task[k];

			int involved_thread[3];
			int num_involved_threads = 1;

			int minislab_id = tdw->thread_in_minibin[ task_slab_id ][ minibin_id ];

			involved_thread[ 0 ] = minislab_id;

			// If the reference thread is not the last one, check if the
			// right thread is also involved because of the halo space
			if (minislab_id < numThreadDataWrappersPerTask[task_slab_id] - 1)
			{
				// check if the atom is in the halo space of the right thread
				if (y >= tdw->minislab_max_y[ task_slab_id ][ minislab_id ] - halo_layer_size)
					involved_thread[ num_involved_threads++ ] = minislab_id + 1;
			}
			// If the reference thread is not the first one, check if the
			// left task is also involved because of the halo space
			if (minislab_id > 0)
			{
				// check if the atom is in the halo space of the below thread
				if (y <= tdw->minislab_max_y[ task_slab_id ][ minislab_id-1 ] + halo_layer_size)
					involved_thread[ num_involved_threads++ ] = minislab_id - 1;
			}
			// For each task, atom data is loaded by 1-3 threads
			for (int j = 0; j < num_involved_threads; j++)
			{
				#ifdef ENABLE_BOOST_THREADS
				// boost::mutex::scoped_lock scopedLock(mutex);
				#endif
				#if !defined(NO_CGAL_PATCHING)
				thread_data_wrapper[ task_slab_id*numThreadDataWrappers + involved_thread[j] ].l.push_front(Weighted_point(Point3(x,y,z), r2, i));
				#else
				// nopatch
				thread_data_wrapper[ task_slab_id*numThreadDataWrappers + involved_thread[j] ].l.push_front(std::make_pair(Weighted_point(Point3(x,y,z), r2), i));
				#endif
			}
		}
	}
}

#endif


void ConnollySurface::BuildPointCells (Fixed_alpha_shape_3 &alpha_shape, vector<ConnollyCell*> &sesComplex,
									   #if !defined(NEW_ATOM_PATCHES)
									   PointCell **atomPatches,
									   #else
									   int atomPatches[],
									   #endif
									   vector<int> &exposed, int num_cells[]
									   #if defined(OPTIMIZE_BUILDING_MEMORY)
									   , TaggedDataWrapper *tagged_data_wrapper,
									   int task_id, int thread_id, double grid_pars[]
									   #endif
									   )
{
	for (Finite_Vertex_Iterator fvit = alpha_shape.finite_vertices_begin(); alpha_shape.finite_vertices_end() != fvit; fvit++)
	{
		int ctype = alpha_shape.classify(fvit);

		// vertex classification is equal to face classification
		// if you go here, the atom is exposed and this a point cell
		if (ctype == Fixed_alpha_shape_3::REGULAR || ctype == Fixed_alpha_shape_3::SINGULAR)
		{
			// retrieve atom index
			#if !defined(NO_CGAL_PATCHING)
			int atom_id = fvit->point().index();
			#else
			// nopatch
			int atom_id = fvit->info();
			#endif

			#if defined(OPTIMIZE_BUILDING_MEMORY)
			// For some reasons, the execution is incorrect if this code is used

			double downz = delphi->atoms[ atom_id ].pos[2] - delphi->atoms[ atom_id ].radius;
			double slab_border;

			if (task_id == 0)
				/// slab_border = gzmax + 0.5*gside;
				slab_border = delphi->zmax + delphi->hside;
			else
				slab_border = tagged_data_wrapper->slab_max_z[ task_id-1 ];

			if (downz < slab_border)
				continue;

			if (downz > tagged_data_wrapper->slab_max_z[ task_id ])
				continue;

			double downy = delphi->atoms[ atom_id ].pos[1] - delphi->atoms[ atom_id ].radius;
			double minislab_border;

			if (thread_id == 0)
				minislab_border = delphi->ymax + delphi->hside;
			else
				minislab_border = tagged_data_wrapper->minislab_max_y[ task_id ][ thread_id-1 ];

			if (downy < minislab_border)
				continue;

			if (downy > tagged_data_wrapper->minislab_max_y[ task_id ][ thread_id ])
				continue;
			#endif // OPTIMIZE_BUILDING_MEMORY

			PointCell *pc = new PointCell();
			sesComplex.push_back(pc);
			pc->id = atom_id;
			pc->patch_type = POINT_CELL;
			#if !defined(NEW_ATOM_PATCHES)
			atomPatches[pc->id] = pc;
			#else
			// Retrieve index of element within sesComplex[*] for later purposes
			atomPatches[pc->id] = sesComplex.size() - 1;
			#endif
			++num_cells[POINT_CELL];

			#if defined(CHECK_BUILDUP_DIFF)
			pc->tag = pc->id;
			#endif

			exposed.push_back(pc->id);

			// get all buried neighbours
			vector<Alpha_Edge> ll;
			alpha_shape.incident_edges(fvit,std::back_inserter(ll));

			for (unsigned int i=0; i<ll.size(); i++)
			{
				int ctype = alpha_shape.classify(ll[i]);

				// connection to buried atoms. Get the hidden clipping sphere
				if (ctype == Fixed_alpha_shape_3::INTERIOR)
				{
					const Vertex_handle &v1 = ll[i].first->vertex(ll[i].second);
					const Vertex_handle &v2 = ll[i].first->vertex(ll[i].third);
					double diff[3];

					diff[0] = v2->point().x() - v1->point().x();
					diff[1] = v2->point().y() - v1->point().y();
					diff[2] = v2->point().z() - v1->point().z();

					double dij2 = DOT(diff,diff);
					// double dij = sqrt(dij2);
					// ratio = (...) / (2.0*dij)
					double ratio = (v1->point().weight() + dij2 - v2->point().weight()) * 0.5;

					#if !defined(OPTIMIZE_CELL_STRUCTURE)
					// tt2 = .... - ratio*ratio
					double tt2 = v1->point().weight() - ratio*ratio / dij2;
					double bigR = sqrt(tt2);

					// nopatch
					int ind1 = v1->info();
					int ind2 = v2->info();

					double r1 = delphi->atoms[ind1].radius;
					// double r2 = delphi->atoms[ind2].radius;

					EdgeCell *ec = new EdgeCell();
					pc->buried_neighbours.push_back(ec);

					ec->id[0] = ind1;
					ec->id[1] = ind2;

					ec->major_radius = bigR;

					double rj2[3], cof;
					// cof = ratio / dij
					cof = ratio / dij2;
					VEC_MUL(rj2,diff,cof)

					double torus_center[3];
					torus_center[0] = rj2[0] + v1->point().x();
					torus_center[1] = rj2[1] + v1->point().y();
					torus_center[2] = rj2[2] + v1->point().z();

					double torus_center_v1_dist2 = DOT(rj2,rj2);

					// 2) get clipping sphere
					// 2.1) get an arbitrary probe position
					// get by Pitagora theorem the distance from the torus center to the probe
					double torus_center_probe_dist2 = v1->point().weight() - torus_center_v1_dist2;
					double rcoi = sqrt(torus_center_probe_dist2);
					ec->rcoi = rcoi;
					// get an arbitrary probe position
					// 2.2) get the radius and center of the visibility sphere
					double w[3],ty;
					NORMALIZE_S_ASSIGN(w,diff,ty)

					// double d = -DOT(torus_center,w);

					double va[3], u[3], ttt;
					va[0] = 0;
					va[1] = w[2];
					va[2] = -w[1];

					NORMALIZE_S_ASSIGN(u,va,ttt)

					// get a random probe position
					double RandomProbe[3];
					ADD_MUL(RandomProbe,torus_center,u,rcoi)

					double temp1[3];
					temp1[0] = RandomProbe[0] - v1->point().x();
					temp1[1] = RandomProbe[1] - v1->point().y();
					temp1[2] = RandomProbe[2] - v1->point().z();

					double pai = sqrt(DOT(temp1,temp1));

					double temp2[3];
					temp2[0] = RandomProbe[0] - v2->point().x();
					temp2[1] = RandomProbe[1] - v2->point().y();
					temp2[2] = RandomProbe[2] - v2->point().z();

					double paj = sqrt(DOT(temp2,temp2));

					double xx[3];
					xx[0] = temp1[0]*(r1/pai);
					xx[1] = temp1[1]*(r1/pai);
					xx[2] = temp1[2]*(r1/pai);

					double aiaj[3];
					aiaj[0] = v2->point().x() - v1->point().x();
					aiaj[1] = v2->point().y() - v1->point().y();
					aiaj[2] = v2->point().z() - v1->point().z();

					double ctemp[3], coeff;
					coeff = pai/(pai+paj);
					VEC_MUL(ctemp,aiaj,coeff)

					double clipping_center[3];
					clipping_center[0] = ctemp[0] + v1->point().x();
					clipping_center[1] = ctemp[1] + v1->point().y();
					clipping_center[2] = ctemp[2] + v1->point().z();

					double ttemp[3];
					SUB(ttemp,xx,ctemp)
					double r_clipping = sqrt(DOT(ttemp,ttemp));

					ASSIGN(ec->clipping_center,clipping_center)
					ec->clipping_radius = r_clipping;

					#else // OPTIMIZE_CELL_STRUCTURE

					// tt2 = .... - ratio*ratio
					// double tt2 = v1->point().weight() - ratio*ratio / dij2;
					// double bigR = sqrt(tt2);

					// nopatch
					int ind1 = v1->info();
					// int ind2 = v2->info();

					double r1 = delphi->atoms[ind1].radius;
					// double r2 = delphi->atoms[ind2].radius;

					// ec->id[0] = ind1;
					// ec->id[1] = ind2;

					// ec->major_radius = bigR;

					double rj2[3], cof;
					// cof = ratio / dij
					cof = ratio / dij2;
					VEC_MUL(rj2,diff,cof)

					double torus_center[3];
					torus_center[0] = rj2[0] + v1->point().x();
					torus_center[1] = rj2[1] + v1->point().y();
					torus_center[2] = rj2[2] + v1->point().z();

					double torus_center_v1_dist2 = DOT(rj2,rj2);

					// 2) get clipping sphere
					// 2.1) get an arbitrary probe position
					// get by Pitagora theorem the distance from the torus center to the probe
					double torus_center_probe_dist2 = v1->point().weight() - torus_center_v1_dist2;
					double rcoi = sqrt(torus_center_probe_dist2);
					// ec->rcoi = rcoi;
					// get an arbitrary probe position
					// 2.2) get the radius and center of the visibility sphere
					double w[3],ty;
					NORMALIZE_S_ASSIGN(w,diff,ty)

					// double d = -DOT(torus_center,w);

					double va[3], u[3], ttt;
					va[0] = 0;
					va[1] = w[2];
					va[2] = -w[1];

					NORMALIZE_S_ASSIGN(u,va,ttt)

					// get a random probe position
					double RandomProbe[3];
					ADD_MUL(RandomProbe,torus_center,u,rcoi)

					double temp1[3];
					temp1[0] = RandomProbe[0] - v1->point().x();
					temp1[1] = RandomProbe[1] - v1->point().y();
					temp1[2] = RandomProbe[2] - v1->point().z();

					double pai = sqrt(DOT(temp1,temp1));

					double temp2[3];
					temp2[0] = RandomProbe[0] - v2->point().x();
					temp2[1] = RandomProbe[1] - v2->point().y();
					temp2[2] = RandomProbe[2] - v2->point().z();

					double paj = sqrt(DOT(temp2,temp2));

					double xx[3];
					xx[0] = temp1[0]*(r1/pai);
					xx[1] = temp1[1]*(r1/pai);
					xx[2] = temp1[2]*(r1/pai);

					double aiaj[3];
					aiaj[0] = v2->point().x() - v1->point().x();
					aiaj[1] = v2->point().y() - v1->point().y();
					aiaj[2] = v2->point().z() - v1->point().z();

					double ctemp[3], coeff;
					coeff = pai/(pai+paj);
					VEC_MUL(ctemp,aiaj,coeff)

					double clipping_center[3];
					clipping_center[0] = ctemp[0] + v1->point().x();
					clipping_center[1] = ctemp[1] + v1->point().y();
					clipping_center[2] = ctemp[2] + v1->point().z();

					double ttemp[3];
					SUB(ttemp,xx,ctemp)
					double r2_clipping = DOT(ttemp,ttemp);

					pc->buried_neighbour_data.push_back(clipping_center[0]);
					pc->buried_neighbour_data.push_back(clipping_center[1]);
					pc->buried_neighbour_data.push_back(clipping_center[2]);
					pc->buried_neighbour_data.push_back(r2_clipping);
					#endif // OPTIMIZE_CELL_STRUCTURE
				}
			}

			#if !defined(OPTIMIZE_CELL_STRUCTURE)
			// cout << endl << pc->buried_neighbours.size();
			#else
			// cout << endl << ((int)(pc->buried_neighbour_data.size()) >> 2);
			#endif
		}
		// else
		//	cout << endl << "buried";
	}
}


void ConnollySurface::BuildFacetCells (Fixed_alpha_shape_3 &alpha_shape, vector<ConnollyCell*> &sesComplex,
									   Octree<vector<FacetCell*>> &gridProbesMap,
									   #if !defined(NEW_ATOM_PATCHES)
									   PointCell **atomPatches,
									   #else
									   int atomPatches[],
									   #endif
									   int num_cells[], double grid_pars[])
{
	set<FacetCell*> checkList;

	double gxmin = grid_pars[0];
	double gymin = grid_pars[1];
	double gzmin = grid_pars[2];
	double gscale = grid_pars[6];
	int ggrid = (int)(grid_pars[7] + 1.E-7);

	for (Finite_Facets_Iterator ffit = alpha_shape.finite_facets_begin(); ffit != alpha_shape.finite_facets_end(); ffit++)
	{
		int ctype = alpha_shape.classify(*ffit);

		if (ctype == Fixed_alpha_shape_3::REGULAR || ctype == Fixed_alpha_shape_3::SINGULAR)
		{
			bool isRegular = true;
			// probe touches both sides
			if (ctype == Fixed_alpha_shape_3::SINGULAR)
			{
				isRegular = false;
			}

			double P1[3], P2[3], P3[3], P12[3], P13[3], P23[3], n[3];

			// int index0;
			int index1, index2, index3;
			// double w0;
			double w1, w2, w3;

			// get the mirror facet such that the outwards normal is correct
			// depending on the classification of the cell
			// which gives the same facet viewed from the other adjacent cells
			if (alpha_shape.classify(ffit->first) == Fixed_alpha_shape_3::INTERIOR)
			{
				const Alpha_Facet &facetp = alpha_shape.mirror_facet(*ffit);
				// opposite point
				// const Weighted_point &p0 = facetp.first->vertex((facetp.second)&3)->point();
				// get points on facet
				const Weighted_point &p1 = facetp.first->vertex((facetp.second+1)&3)->point();
				const Weighted_point &p2 = facetp.first->vertex((facetp.second+2)&3)->point();
				const Weighted_point &p3 = facetp.first->vertex((facetp.second+3)&3)->point();
				// get the proper normal
				// const Vector3 nn = (facetp.second % 2 == 1) ? CGAL::normal(p1, p2, p3) : CGAL::normal(p1, p3, p2);
				const Vector3 nn = (facetp.second % 2 == 1) ? CGAL::normal((Point3)p1, (Point3)p2, (Point3)p3) : CGAL::normal((Point3)p1, (Point3)p3, (Point3)p2);

				n[0] = nn.x();
				n[1] = nn.y();
				n[2] = nn.z();

				P1[0] = p1.x();
				P1[1] = p1.y();
				P1[2] = p1.z();

				P2[0] = p2.x();
				P2[1] = p2.y();
				P2[2] = p2.z();

				P3[0] = p3.x();
				P3[1] = p3.y();
				P3[2] = p3.z();

				// nopatch
				// index0 = facetp.first->vertex((facetp.second+0)&3)->info();
				index1 = facetp.first->vertex((facetp.second+1)&3)->info();
				index2 = facetp.first->vertex((facetp.second+2)&3)->info();
				index3 = facetp.first->vertex((facetp.second+3)&3)->info();

				// w0 = p0.weight();
				w1 = p1.weight();
				w2 = p2.weight();
				w3 = p3.weight();
			}
			else
			{
				const Alpha_Facet *facetp = &(*ffit);
				// opposite point
				// const Weighted_point &p0 = facetp->first->vertex((facetp->second)&3)->point();
				// get points on facet
				const Weighted_point &p1 = facetp->first->vertex((facetp->second+1)&3)->point();
				const Weighted_point &p2 = facetp->first->vertex((facetp->second+2)&3)->point();
				const Weighted_point &p3 = facetp->first->vertex((facetp->second+3)&3)->point();
				// get the proper normal
				// const Vector3 nn = (facetp->second % 2 == 1) ? CGAL::normal(p1, p2, p3) : CGAL::normal(p1, p3, p2);
				const Vector3 nn = (facetp->second % 2 == 1) ? CGAL::normal((Point3)p1, (Point3)p2, (Point3)p3) : CGAL::normal((Point3)p1, (Point3)p3, (Point3)p2);

				n[0] = nn.x();
				n[1] = nn.y();
				n[2] = nn.z();

				P1[0] = p1.x();
				P1[1] = p1.y();
				P1[2] = p1.z();

				P2[0] = p2.x();
				P2[1] = p2.y();
				P2[2] = p2.z();

				P3[0] = p3.x();
				P3[1] = p3.y();
				P3[2] = p3.z();

				// nopatch
				// index0 = facetp->first->vertex((facetp->second+0)&3)->info();
				index1 = facetp->first->vertex((facetp->second+1)&3)->info();
				index2 = facetp->first->vertex((facetp->second+2)&3)->info();
				index3 = facetp->first->vertex((facetp->second+3)&3)->info();

				// w0 = p0.weight();
				w1 = p1.weight();
				w2 = p2.weight();
				w3 = p3.weight();
			}
			SUB(P12,P2,P1)
			SUB(P13,P3,P1)
			SUB(P23,P3,P2)

			// get probe center
			// const double A = P12.squared_length();
			// const double B = P23.squared_length();
			// const double C = P13.squared_length();
			const double A = DOT(P12,P12);
			const double B = DOT(P23,P23);
			const double C = DOT(P13,P13);
			const double D = w1;
			const double E = w2;
			const double F = w3;
			const double X = D+E-A;
			const double Y = E+F-B;
			const double Z = D+F-C;
			const double Varg = 4.*D*E*F-F*X*X-D*Y*Y-E*Z*Z+X*Y*Z;

			if (Varg < 0.)
			{
				cout << endl << WARN << "Negative sqrt argument of V=1/12*sqrt(...) " << Varg;
				exit(-1);
			}
			const double V = 1./12.*sqrt(Varg);

			// const double G = P12*P13;
			const double G = DOT(P12,P13);
			// const double H = n*n;
			const double H = DOT(n,n);
			// const double I = P2*P2-P1*P1;
			const double I = (DOT(P2,P2))-(DOT(P1,P1));
			// const double J = P3*P3-P1*P1;
			const double J = (DOT(P3,P3))-(DOT(P1,P1));
			// const double K = n*P1;
			const double K = DOT(n,P1);
			const double L = I-E+D;
			const double M = J-F+D;
			const double S = C*L-G*M;
			const double T = A*M-G*L;

			// Vector3 Q = (K*n+0.5*(S*P12+T*P13))/H;
			double Q[3];

			// Q[0] = (K/H)*n[0]+(S/(2*H))*P12[0]+(T/(2*H))*P13[0];
			// Q[1] = (K/H)*n[1]+(S/(2*H))*P12[1]+(T/(2*H))*P13[1];
			// Q[2] = (K/H)*n[2]+(S/(2*H))*P12[2]+(T/(2*H))*P13[2];
			Q[0] = (K*n[0]+0.5*(S*P12[0]+T*P13[0]))/H;
			Q[1] = (K*n[1]+0.5*(S*P12[1]+T*P13[1]))/H;
			Q[2] = (K*n[2]+0.5*(S*P12[2]+T*P13[2]))/H;
			// Vector3 Probe = Q + ((6*V)/(H))*n;
			double probe[3], coeff;

			coeff = (6.0*V) / H;

			probe[0] = Q[0] + coeff*n[0];
			probe[1] = Q[1] + coeff*n[1];
			probe[2] = Q[2] + coeff*n[2];

			FacetCell *fc1 = new FacetCell();
			if (isRegular)
				fc1->patch_type = REGULAR_FACE_CELL;
			else
				fc1->patch_type = SINGULAR_FACE_CELL;

			sesComplex.push_back(fc1);
			++num_cells[ fc1->patch_type ];

			// fc1->center[0] = Probe.x();
			// fc1->center[1] = Probe.y();
			// fc1->center[2] = Probe.z();

			fc1->center[0] = probe[0];
			fc1->center[1] = probe[1];
			fc1->center[2] = probe[2];

			#if defined(CHECK_BUILDUP_DIFF)
			fc1->tag = fc1->center[0] + fc1->center[1] + fc1->center[2];
			#endif

			fc1->id[0] = index1;
			fc1->id[1] = index2;
			fc1->id[2] = index3;

			#if !defined(OPTIMIZE_CELL_STRUCTURE)
			#if !defined(NEW_ATOM_PATCHES)
			atomPatches[index1]->incidentProbes.push_back(fc1);
			atomPatches[index2]->incidentProbes.push_back(fc1);
			atomPatches[index3]->incidentProbes.push_back(fc1);
			#else
			PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[index1] ];
			PointCell *pc2 = (PointCell *)sesComplex[ atomPatches[index2] ];
			PointCell *pc3 = (PointCell *)sesComplex[ atomPatches[index3] ];

			pc1->incidentProbes.push_back(fc1);
			pc2->incidentProbes.push_back(fc1);
			pc3->incidentProbes.push_back(fc1);
			#endif
			#else // OPTIMIZE_CELL_STRUCTURE
			PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[index1] ];
			PointCell *pc2 = (PointCell *)sesComplex[ atomPatches[index2] ];
			PointCell *pc3 = (PointCell *)sesComplex[ atomPatches[index3] ];

			pc1->incidentProbes.push_back( sesComplex.size()-1 );
			pc2->incidentProbes.push_back( sesComplex.size()-1 );
			pc3->incidentProbes.push_back( sesComplex.size()-1 );
			#endif // OPTIMIZE_CELL_STRUCTURE

			// Vector3 PP1, PP2, PP3;
			double PP1[3], PP2[3], PP3[3], YY[3];

			// if it is singular this automatically clips the sphere without adding planes
			if (!isRegular)
			{
				// tetrahedron points
				// PP1 = P1;
				ASSIGN(PP1,P1);
				// PP2 = P2;
				ASSIGN(PP2,P2);
				// PP3 = P3;
				ASSIGN(PP3,P3);
				fc1->isSelfIntersecting = true;
			}
			// avoid unecessary self intersection clipping
			else
			{
				// tetrahedron points
				// PP1 = Probe + 50000.*(P1-Probe);
				SUB(YY,P1,probe)
				ADD_MUL(PP1,probe,YY,50000.)
				// PP2 = Probe + 50000.*(P2-Probe);
				SUB(YY,P2,probe)
				ADD_MUL(PP2,probe,YY,50000.)
				// PP3 = Probe + 50000.*(P3-Probe);
				SUB(YY,P3,probe)
				ADD_MUL(PP3,probe,YY,50000.)
				fc1->isSelfIntersecting = false;
			}
			// 0 meaning that we have no additional planes.
			// If numSelfIntersections == 0 and isSelfIntersecting == 0 then this means that
			// the only self intersections is due to mirror facet

			// double t1[3], t2[3], t3[3], probe[3];
			double *t1, *t2, *t3;
			t1 = PP1;
			t2 = PP2;
			t3 = PP3;

			// t1[0] = PP1.x();
			// t1[1] = PP1.y();
			// t1[2] = PP1.z();

			// t2[0] = PP2.x();
			// t2[1] = PP2.y();
			// t2[2] = PP2.z();

			// t3[0] = PP3.x();
			// t3[1] = PP3.y();
			// t3[2] = PP3.z();

			// probe[0] = Probe.x();
			// probe[1] = Probe.y();
			// probe[2] = Probe.z();

			plane3points(t1,t2,t3,fc1->planes[0]);
			plane3points(t1,t2,probe,fc1->planes[1]);
			plane3points(t2,t3,probe,fc1->planes[2]);
			plane3points(t1,t3,probe,fc1->planes[3]);

			fc1->plane_indices[0][0] = fc1->id[0];
			fc1->plane_indices[0][1] = fc1->id[1];

			fc1->plane_indices[1][0] = fc1->id[1];
			fc1->plane_indices[1][1] = fc1->id[2];

			fc1->plane_indices[2][0] = fc1->id[0];
			fc1->plane_indices[2][1] = fc1->id[2];

			// fix orientation

			// interior reference point
			// Vector3 inner = ((PP1+PP2+PP3)/3.0+Probe)*0.5;
			double inner2[3];

			// inner2[0] = inner.x();
			// inner2[1] = inner.y();
			// inner2[2] = inner.z();

			inner2[0] = ((1./3.)*(PP1[0]+PP2[0]+PP3[0])+probe[0])*0.5;
			inner2[1] = ((1./3.)*(PP1[1]+PP2[1]+PP3[1])+probe[1])*0.5;
			inner2[2] = ((1./3.)*(PP1[2]+PP2[2]+PP3[2])+probe[2])*0.5;

			for (int ff=0; ff<4; ff++)
			{
				double orient = DOT(inner2,fc1->planes[ff]);
				// ok
				if (orient<-(fc1->planes[ff][3]))
				{}
				// swap
				else
				{
					fc1->planes[ff][0] = -fc1->planes[ff][0];
					fc1->planes[ff][1] = -fc1->planes[ff][1];
					fc1->planes[ff][2] = -fc1->planes[ff][2];
					fc1->planes[ff][3] = -fc1->planes[ff][3];
				}
			}

			if (ctype == Fixed_alpha_shape_3::SINGULAR)
			{
				// Vector3 Probe = Q - ((6*V)/(H))*n;
				double probe[3];

				probe[0] = Q[0] - coeff*n[0];
				probe[1] = Q[1] - coeff*n[1];
				probe[2] = Q[2] - coeff*n[2];

				FacetCell *fc2 = new FacetCell();
				fc2->mirrorCell = fc1;
				fc1->mirrorCell = fc2;
				fc2->patch_type = SINGULAR_FACE_CELL;
				fc2->isSelfIntersecting = true;

				sesComplex.push_back(fc2);
				++num_cells[SINGULAR_FACE_CELL];

				// fc2->center[0] = Probe.x();
				// fc2->center[1] = Probe.y();
				// fc2->center[2] = Probe.z();

				fc2->center[0] = probe[0];
				fc2->center[1] = probe[1];
				fc2->center[2] = probe[2];

				#if defined(CHECK_BUILDUP_DIFF)
				fc2->tag = fc2->center[0] + fc2->center[1] + fc2->center[2];
				#endif

				fc2->id[0] = index1;
				fc2->id[1] = index2;
				fc2->id[2] = index3;

				#if !defined(OPTIMIZE_CELL_STRUCTURE)
				#if !defined(NEW_ATOM_PATCHES)
				atomPatches[index1]->incidentProbes.push_back(fc2);
				atomPatches[index2]->incidentProbes.push_back(fc2);
				atomPatches[index3]->incidentProbes.push_back(fc2);
				#else
				PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[index1] ];
				PointCell *pc2 = (PointCell *)sesComplex[ atomPatches[index2] ];
				PointCell *pc3 = (PointCell *)sesComplex[ atomPatches[index3] ];

				pc1->incidentProbes.push_back(fc2);
				pc2->incidentProbes.push_back(fc2);
				pc3->incidentProbes.push_back(fc2);
				#endif
				#else // OPTIMIZE_CELL_STRUCTURE
				PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[index1] ];
				PointCell *pc2 = (PointCell *)sesComplex[ atomPatches[index2] ];
				PointCell *pc3 = (PointCell *)sesComplex[ atomPatches[index3] ];

				pc1->incidentProbes.push_back( sesComplex.size()-1 );
				pc2->incidentProbes.push_back( sesComplex.size()-1 );
				pc3->incidentProbes.push_back( sesComplex.size()-1 );
				#endif // OPTIMIZE_CELL_STRUCTURE

				// probe[0] = Probe.x();
				// probe[1] = Probe.y();
				// probe[2] = Probe.z();

				plane3points(t1,t2,t3,fc2->planes[0]);
				plane3points(t1,t2,probe,fc2->planes[1]);
				plane3points(t2,t3,probe,fc2->planes[2]);
				plane3points(t1,t3,probe,fc2->planes[3]);

				fc2->plane_indices[0][0] = fc2->id[0];
				fc2->plane_indices[0][1] = fc2->id[1];

				fc2->plane_indices[1][0] = fc2->id[1];
				fc2->plane_indices[1][1] = fc2->id[2];

				fc2->plane_indices[2][0] = fc2->id[0];
				fc2->plane_indices[2][1] = fc2->id[2];

				// interior reference point

				// Vector3 inner = ((PP1+PP2+PP3)/3.0+Probe)*0.5;
				// inner2[0] = inner.x();
				// inner2[1] = inner.y();
				// inner2[2] = inner.z();

				inner2[0] = ((1./3.)*(PP1[0]+PP2[0]+PP3[0])+probe[0])*0.5;
				inner2[1] = ((1./3.)*(PP1[1]+PP2[1]+PP3[1])+probe[1])*0.5;
				inner2[2] = ((1./3.)*(PP1[2]+PP2[2]+PP3[2])+probe[2])*0.5;

				for (int ff=0; ff<4; ff++)
				{
					double orient = DOT(inner2,fc2->planes[ff]);
					// ok
					if (orient >= -fc2->planes[ff][3])
					{
						// swap
						fc2->planes[ff][0] = -fc2->planes[ff][0];
						fc2->planes[ff][1] = -fc2->planes[ff][1];
						fc2->planes[ff][2] = -fc2->planes[ff][2];
						fc2->planes[ff][3] = -fc2->planes[ff][3];
					}
				}
			}
			// no mirror patch
			else
				fc1->mirrorCell = NULL;

			// check if self-intersection must be tested or not
			double plane_dist, proj[3], plane[4];

			plane3points(P1,P2,P3,plane);
			point2plane(fc1->center,plane,&plane_dist,proj);

			if (plane_dist <= probe_radius && ctype != Fixed_alpha_shape_3::SINGULAR)
			{
				set<FacetCell*>::iterator it = checkList.find(fc1);

				if (it != checkList.end())
					continue;

				int ix = (int)rintp((fc1->center[0]-gxmin)*gscale);
				int iy = (int)rintp((fc1->center[1]-gymin)*gscale);
				int iz = (int)rintp((fc1->center[2]-gzmin)*gscale);

				if (ix >= ggrid || iy >= ggrid || iz >= ggrid || ix < 0 || iy < 0 || iz < 0)
				{
					cout << endl << ERR << "Too big probe, to support this probe increase Self_Intersections_Grid_Coefficient\n";
					exit(-1);
				}
				checkList.insert(fc1);

				// int64_t max_ind = (int64_t)SELF_MAP(ix,iy,iz,0,ggrid,ggrid,ggrid);

				// max_ind++;
				// if ((int)max_ind >= MAX_PROBES)
				// {
				// 	cout << endl << ERR << "Increase MAX_PROBES";
				// 	exit(-1);
				// }
				// SELF_MAP(ix,iy,iz,0,ggrid,ggrid,ggrid) = (FacetCell*)max_ind;
				// SELF_MAP(ix,iy,iz,max_ind,ggrid,ggrid,ggrid) = fc1;

				vector<FacetCell*> &vec = gridProbesMap(ix,iy,iz);
				if (vec.size() == 0)
					vec.reserve(10);
				vec.push_back(fc1);
			}
		}
	}
	// cout << endl << INFO << "Checking " << checkList.size() << " probes for self intersections...";
	// cout.flush();
}


void ConnollySurface::BuildEdgeCells (ofstream &of, Fixed_alpha_shape_3 &alpha_shape,
									  vector<ConnollyCell*> &sesComplex,
									  #if !defined(NEW_ATOM_PATCHES)
									  PointCell **atomPatches,
									  #else
									  int atomPatches[],
									  #endif
									  int num_cells[]
									  #if defined(OPTIMIZE_BUILDING_MEMORY)
									  , TaggedDataWrapper *tagged_data_wrapper,
									  int task_id, int thread_id, double grid_pars[]
									  #endif
									  )
{
	#if defined(CHECK_MAX_RADII)
	static double max_major_radius = 0.0;
	static double max_clipping_radius = 0.0;
	#endif

	for (Finite_Edges_Iterator feit = alpha_shape.finite_edges_begin(); feit != alpha_shape.finite_edges_end(); feit++)
	{
		int ctype = alpha_shape.classify(*feit);

		// exposed edge
		if (ctype == Fixed_alpha_shape_3::REGULAR || ctype == Fixed_alpha_shape_3::SINGULAR)
		{
			// up and down edge vertices
			const Vertex_handle &v1 = feit->first->vertex(feit->second);
			const Vertex_handle &v2 = feit->first->vertex(feit->third);

			// Vector3 diff = v2->point()-v1->point();
			double diff[3];
			diff[0] = v2->point().x() - v1->point().x();
			diff[1] = v2->point().y() - v1->point().y();
			diff[2] = v2->point().z() - v1->point().z();

			// double dij2 = diff.x()*diff.x()+diff.y()*diff.y()+diff.z()*diff.z();
			double dij2 = DOT(diff,diff);
			// double dij = sqrt(dij2);
			// ratio = (...) / (2.0*dij)
			double ratio = (v1->point().weight() + dij2 - v2->point().weight()) * 0.5;
			// tt2 = .... - ratio*ratio
			double tt2 = v1->point().weight() - ratio*ratio / dij2;
			// compute major torus radius
			double bigR = sqrt(tt2);

            // nopatch
			int ind1 = v1->info();
			int ind2 = v2->info();

			double r1 = delphi->atoms[ind1].radius;
			// double r2 = delphi->atoms[ind2].radius;

			// ratio2 = (...) / (2*dij)
			double ratio2 = (v2->point().weight() + dij2 - v1->point().weight()) * 0.5;
			// bigR2 = sqrt(... - ratio2*ratio2)
			double bigR2 = sqrt(v2->point().weight() - ratio2*ratio2 / dij2);

			if (fabs(bigR-bigR2) > 1e-6)
				cout << endl << WARN << "Non precise torus inner radius " << bigR << " "<< bigR2;

			// compute torus clipping sphere:
			// 1) get torus center
			double rj2[3], cof;
			// ratio / dij;
			cof = ratio / dij2;
			VEC_MUL(rj2,diff,cof)

			double torus_center[3];
			torus_center[0] = rj2[0] + v1->point().x();
			torus_center[1] = rj2[1] + v1->point().y();
			torus_center[2] = rj2[2] + v1->point().z();

			double torus_center_v1_dist2 = DOT(rj2,rj2);

			// 2) get clipping sphere
			// 2.1) get an arbitrary probe position
			// get by Pitagora theorem the distance from the torus center to the probe
			double torus_center_probe_dist2 = v1->point().weight() - torus_center_v1_dist2;
			double rcoi = sqrt(torus_center_probe_dist2);

			// get an arbitrary probe position
			// 2.2) get the radius and center of the visibility sphere
			double w[3],ty;
			NORMALIZE_S_ASSIGN(w,diff,ty)

			// double d = -DOT(torus_center,w);

			// plane of the COI is w*x+d=0. Torus center lies on this plane
			// get a random point on the COI (Circle Of Intersection). The COI can be analytically computed.
			double cross_temp[3], va[3], vb[3], u[3], v[3], ttt;
			va[0] = 0;
			va[1] = w[2];
			va[2] = -w[1];
			vb[0] = w[2];
			vb[1] = 0;
			vb[2] = -w[0];

			NORMALIZE_S_ASSIGN(u,va,ttt)
			CROSS(cross_temp,va,vb)
			CROSS(v,cross_temp,va)
			NORMALIZE_S(v,ttt)

			// sample coi for debug
			// double tt = 0;
			// for (double tcoi=0.0; tcoi < TWO_PI; tcoi += 0.1)
			// {
			// 	// sample the equation of the coi
			// 	Vector3 P = torus_center+rcoi*(cos(tcoi)*U+sin(tcoi)*V);
			// 	char buff2[BUFLEN];
			// 	sprintf(buff2,"\n\n sphere { \n <%lf,%lf,%lf>,%lf",(P.x()),(P.y()),(P.z()),0.1);
			// 	of << buff2;
			// 	of << "\n pigment{ color Yellow}}";
			// }

			// get a random probe position
			double RandomProbe[3];
			ADD_MUL(RandomProbe,torus_center,u,rcoi)

			double temp1[3];
			temp1[0] = RandomProbe[0] - v1->point().x();
			temp1[1] = RandomProbe[1] - v1->point().y();
			temp1[2] = RandomProbe[2] - v1->point().z();

			double pai = sqrt(DOT(temp1,temp1));

			double temp2[3];
			temp2[0] = RandomProbe[0] - v2->point().x();
			temp2[1] = RandomProbe[1] - v2->point().y();
			temp2[2] = RandomProbe[2] - v2->point().z();

			double paj = sqrt(DOT(temp2,temp2));

			double aiaj[3];
			aiaj[0] = v2->point().x() - v1->point().x();
			aiaj[1] = v2->point().y() - v1->point().y();
			aiaj[2] = v2->point().z() - v1->point().z();

			double ctemp[3], coeff;
			coeff = pai / (pai+paj);
			VEC_MUL(ctemp,aiaj,coeff)

			double clipping_center[3];
			clipping_center[0] = ctemp[0] + v1->point().x();
			clipping_center[1] = ctemp[1] + v1->point().y();
			clipping_center[2] = ctemp[2] + v1->point().z();

			double xx[3];
			xx[0] = temp1[0]*(r1/pai);
			xx[1] = temp1[1]*(r1/pai);
			xx[2] = temp1[2]*(r1/pai);

			double ttemp[3];
			SUB(ttemp,xx,ctemp)
			double r2_clipping = DOT(ttemp,ttemp);
			double clipping_radius = sqrt(r2_clipping);

			#if defined(CHECK_MAX_RADII)
			max_clipping_radius = fmax(max_clipping_radius, clipping_radius);
			#endif

			#if defined(OPTIMIZE_BUILDING_MEMORY)
			// For some reasons, the execution is incorrect if this code is used

			double downz = clipping_center[2] - clipping_radius;
			double slab_border;

			if (task_id == 0)
				slab_border = delphi->zmax + delphi->hside;
			else
				slab_border = tagged_data_wrapper->slab_max_z[ task_id-1 ];

			if (downz < slab_border)
				continue;

			if (downz > tagged_data_wrapper->slab_max_z[ task_id ])
				continue;

			double downy = clipping_center[1] - clipping_radius;
			double minislab_border;

			if (thread_id == 0)
				minislab_border = delphi->ymax + delphi->hside;
			else
				minislab_border = tagged_data_wrapper->minislab_max_y[ task_id ][ thread_id-1 ];

			if (downy < minislab_border)
				continue;

			if (downy > tagged_data_wrapper->minislab_max_y[ task_id ][ thread_id ])
				continue;
			#endif // OPTIMIZE_BUILDING_MEMORY

			EdgeCell *ec = new EdgeCell();
			sesComplex.push_back(ec);

			ec->id[0] = ind1;
			ec->id[1] = ind2;

			ec->major_radius = bigR;

			#if defined(CHECK_MAX_RADII)
			max_major_radius = fmax(max_major_radius, bigR);
			#endif

			// classify the type of the edge
			if ((ctype == Fixed_alpha_shape_3::REGULAR))
			{
				ec->patch_type = REGULAR_EDGE_CELL;
			}
			else
			{
				ec->patch_type = SINGULAR_EDGE_CELL;
			}
			++num_cells[ ec->patch_type ];

			ec->rcoi = rcoi;

			// ec->u[0] = u[0];
			// ec->u[1] = u[1];
			// ec->u[2] = u[2];

			// ec->v[0] = v[0];
			// ec->v[1] = v[1];
			// ec->v[2] = v[2];

			ASSIGN(ec->center,torus_center)

			#if defined(CHECK_BUILDUP_DIFF)
			ec->tag = ec->center[0] + ec->center[1] + ec->center[2];
			#endif

			//  set of incident probes to a given edge
			FacetCell *fcv[MAX_INCIDENT_PROBES];

			// if it is singular no clipping plane is needed
			ec->isComplex = false;
			int index = 0;

			// if it is regular needs clipping
			// now get the two triangle patches if the edge is regular
			if (ec->patch_type == REGULAR_EDGE_CELL)
			{
				#if !defined(OPTIMIZE_CELL_STRUCTURE)
				#if !defined(NEW_ATOM_PATCHES)
				vector<FacetCell*> &f1 = atomPatches[ind1]->incidentProbes;
				#else
				PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[ind1] ];
				vector<FacetCell*> &f1 = pc1->incidentProbes;
				#endif

				unsigned int np = (unsigned int)f1.size();

				// get all the incident probes. Most of the cases will be 2.
				for (unsigned int indf=0; indf<np; indf++)
				{
					int i1 = (f1)[indf]->id[0];
					int i2 = (f1)[indf]->id[1];
					int i3 = (f1)[indf]->id[2];

					if ((i1 == ind1 && i2 == ind2) ||
						(i2 == ind1 && i1 == ind2) ||
						(i1 == ind1 && i3 == ind2) ||
						(i3 == ind1 && i1 == ind2) ||
						(i2 == ind1 && i3 == ind2) ||
						(i3 == ind1 && i2 == ind2))
					{
						fcv[index++] = (f1)[indf];
					}
				}
				#else // OPTIMIZE_CELL_STRUCTURE
				PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[ind1] ];

				int np = pc1->incidentProbes.size();

				// get all the incident probes. Most of the cases will be 2.
				for (int indf=0; indf<np; indf++)
				{
					FacetCell *f1 = (FacetCell *)sesComplex[ pc1->incidentProbes[indf] ];

					int i1 = f1->id[0];
					int i2 = f1->id[1];
					int i3 = f1->id[2];

					if ((i1 == ind1 && i2 == ind2) ||
						(i2 == ind1 && i1 == ind2) ||
						(i1 == ind1 && i3 == ind2) ||
						(i3 == ind1 && i1 == ind2) ||
						(i2 == ind1 && i3 == ind2) ||
						(i3 == ind1 && i2 == ind2))
					{
						fcv[index++] = f1;
					}
				}
				#endif // OPTIMIZE_CELL_STRUCTURE

				if (fcv[0] == NULL || fcv[1] == NULL)
				{
					cout << endl << ERR << "Error at atoms " << ind1 << "," << ind2;
					cout << endl << ERR << "Regular edge with no probe stations?";
					exit(-1);
				}

				// now get the two planes that clip the torus and decide
				// if the clipping is in "AND" or "OR" mode. That is the angle is respectively  <pi or >pi between planes

				// two border probes: this is the most common situation
				// recover the two reference planes and deduce if the "AND" or "OR" configuaration is needed
				if (index == 2)
				{
					bool found = false;
					// get first plane
					for (int ii=0; ii<3; ii++)
					{
						int i1 = fcv[0]->plane_indices[ii][0];
						int i2 = fcv[0]->plane_indices[ii][1];

						// printf("\n(1) %d,%d (%d,%d)",i1,i2,ind1,ind2);

						if (((i1 == ind1) && (i2 == ind2)) || ((i1 == ind2) && (i2 == ind1)))
						{
							found = true;
							#if !defined(OPTIMIZE_CELL_STRUCTURE)
							ec->cutting_planes[0][0] = -fcv[0]->planes[ii+1][0];
							ec->cutting_planes[0][1] = -fcv[0]->planes[ii+1][1];
							ec->cutting_planes[0][2] = -fcv[0]->planes[ii+1][2];
							ec->cutting_planes[0][3] = -fcv[0]->planes[ii+1][3];
							#else // OPTIMIZE_CELL_STRUCTURE
							// the sign is the opposite; the later use takes care of this
							ec->cutting_planes[0] = fcv[0]->planes[ii+1];
							#endif // OPTIMIZE_CELL_STRUCTURE
							// cout << endl << -fcv[0]->planes[ii+1][0] << "," << -fcv[0]->planes[ii+1][1] << "," <<-fcv[0]->planes[ii+1][2] << "," << -fcv[0]->planes[ii+1][3];
							break;
						}
					}
					if (!found)
					{
						cout << endl << ERR << "Cannot detect correct plane!";
						exit(-1);
					}
					found = false;

					// get second plane
					for (int ii=0; ii<3; ii++)
					{
						int i1 = fcv[1]->plane_indices[ii][0];
						int i2 = fcv[1]->plane_indices[ii][1];

						// printf("\n(2) %d,%d (%d,%d)",i1,i2,ind1,ind2);

						if (((i1 == ind1) && (i2 == ind2)) || ((i1 == ind2) && (i2 == ind1)))
						{
							found = true;
							#if !defined(OPTIMIZE_CELL_STRUCTURE)
							ec->cutting_planes[1][0] = -fcv[1]->planes[ii+1][0];
							ec->cutting_planes[1][1] = -fcv[1]->planes[ii+1][1];
							ec->cutting_planes[1][2] = -fcv[1]->planes[ii+1][2];
							ec->cutting_planes[1][3] = -fcv[1]->planes[ii+1][3];
							#else // OPTIMIZE_CELL_STRUCTURE
							// the sign is the opposite; the later use takes care of this
							ec->cutting_planes[1] = fcv[1]->planes[ii+1];
							#endif // OPTIMIZE_CELL_STRUCTURE
							// cout << endl << -fcv[1]->planes[ii+1][0] << "," << -fcv[1]->planes[ii+1][1] << "," <<-fcv[1]->planes[ii+1][2] << "," << -fcv[1]->planes[ii+1][3];
							break;
						}
					}

					if (!found)
					{
						cout << endl << ERR << "Cannot detect correct plane!";
						exit(-1);
					}

					// get the two reference vectors of the planes.
					// double vref1[3],vref2[3];

					// SUB(vref1,fcv[0]->center,ec->center)
					// SUB(vref2,fcv[1]->center,ec->center)

					// decide AND or TEST based on orientation

					ec->acute = orientation(fcv[0]->center, fcv[1]->center, ec->cutting_planes[0], ec->cutting_planes[1]);
					#if defined(OPTIMIZE_CELL_STRUCTURE)
					// the sign has to be changed because if OPTIMIZE_CELL_STRUCTURE is defined
					// existing planes are used but the sign of the elements is wrong
					ec->acute = !ec->acute;
					#endif
				}
				// more than 2 incident probes: complicate situation ...
				else if (index > 2)
				{
					// indices of the sorted set of probes
					int sorted[MAX_INCIDENT_PROBES];

					if ((index%2) != 0)
					{
						cout << endl << ERR << "Torus circulator gives an odd number of probes!";
						exit(-1);
					}

					// manage this situation by explictly trasversing the set of stations of the probes along the current edge
					// sort the probe positions along the COI in order to get a pair of planes every time.
					sortProbes(ec,fcv,index,sorted);

					// as first index find a singular face. Choose the direction by
					// picking a non singular face
					int start_index = 0;
					int direction = 0;

					// very probably at least one facet is singular: do fast detection of the first probe
					for (int u=0; u<index; u++)
					{
						if (fcv[sorted[u]]->patch_type == SINGULAR_FACE_CELL)
						{
							start_index = u;
							int next = (u+1)%index;
							int prev = (index+u-1)%index;

							if (fcv[sorted[u]]->mirrorCell == fcv[sorted[next]])
								direction = -2;
							else if (fcv[sorted[u]]->mirrorCell == fcv[sorted[prev]])
								direction = +2;
							else
							{
								cout << endl << ERR << "Inconsistency in torus circulator";
								exit(-1);
							}
							break;
						}
					}

					// no singular facets, use other, slower criterion
					if (direction == 0)
					{
						ec->isComplex = true;
						double *c1 = fcv[sorted[0]]->center;
						double *c2 = fcv[sorted[1]]->center;
						double mid2[3],plane[3];

						MID(mid2,c1,c2)
						SUB(mid2,mid2,ec->center)

						bool found = false;
						// get refence plane for the test
						for (int ii=0; ii<3; ii++)
						{
							int i1 = fcv[sorted[0]]->plane_indices[ii][0];
							int i2 = fcv[sorted[0]]->plane_indices[ii][1];

							if (((i1 == ind1) && (i2 == ind2)) || ((i1 == ind2) && (i2 == ind1)))
							{
								found = true;

								plane[0] = -fcv[sorted[0]]->planes[ii+1][0];
								plane[1] = -fcv[sorted[0]]->planes[ii+1][1];
								plane[2] = -fcv[sorted[0]]->planes[ii+1][2];
								#if !defined(OPTIMIZE_CELL_STRUCTURE)
								// cout << endl << ec->cutting_planes[0][0] << "," << ec->cutting_planes[0][1] << "," << ec->cutting_planes[0][2] << "," << ec->cutting_planes[0][3];
								#else
								// cout << endl << -ec->cutting_planes[0][0] << "," << -ec->cutting_planes[0][1] << "," << -ec->cutting_planes[0][2] << "," << -ec->cutting_planes[0][3];
								#endif
								break;
							}
						}
						if (!found)
						{
							cout << endl <<ERR<< "Cannot identify plane!";
							exit(-1);
						}

						double test = -DOT(mid2,plane);
						start_index = 0;
						if (test > 0)
							direction = +2;
						else
							direction = -2;
					}

					int currentProbe1 = start_index;
					int currentProbe2 = (index + start_index + (direction>>1))%index;

					// cout << endl << "Choosen direction "<<direction;

					// get all needed planes
					while (1)
					{
						FacetCell *fc1 = fcv[sorted[currentProbe1]];
						FacetCell *fc2 = fcv[sorted[currentProbe2]];

						bool found = false;
						double *plane1, *plane2;

						// get first plane
						for (int ii=0; ii<3; ii++)
						{
							int i1 = fc1->plane_indices[ii][0];
							int i2 = fc1->plane_indices[ii][1];

							// printf("\n(1) %d,%d (%d,%d)",i1,i2,ind1,ind2);

							if (((i1 == ind1) && (i2 == ind2)) || ((i1 == ind2) && (i2 == ind1)))
							{
								found = true;
								#if !defined(OPTIMIZE_CELL_STRUCTURE)
								double *plane = allocateVector<double>(4);
								plane1 = plane;
								plane[0] = -fc1->planes[ii+1][0];
								plane[1] = -fc1->planes[ii+1][1];
								plane[2] = -fc1->planes[ii+1][2];
								plane[3] = -fc1->planes[ii+1][3];
								ec->additional_planes.push_back(plane);
								// cout << endl << ec->cutting_planes[0][0] << "," << ec->cutting_planes[0][1] << "," << ec->cutting_planes[0][2] << "," << ec->cutting_planes[0][3];
								#else
								plane1 = fc1->planes[ii+1];
								ec->additional_planes.push_back(fc1->planes[ii+1]);
								// cout << endl << -ec->cutting_planes[0][0] << "," << -ec->cutting_planes[0][1] << "," << -ec->cutting_planes[0][2] << "," << -ec->cutting_planes[0][3];
								#endif
								break;
							}
						}

						if (!found)
						{
							cout << endl << ERR << "Cannot detect correct plane!";
							exit(-1);
						}

						found = false;

						// get second plane
						for (int ii=0; ii<3; ii++)
						{
							int i1 = fc2->plane_indices[ii][0];
							int i2 = fc2->plane_indices[ii][1];

							// printf("\n(2) %d,%d (%d,%d)", i1,i2,ind1,ind2);

							if ((i1 == ind1 && i2 == ind2) || (i1 == ind2 && i2 == ind1))
							{
								found = true;
								#if !defined(OPTIMIZE_CELL_STRUCTURE)
								double *plane = allocateVector<double>(4);
								plane2 = plane;
								plane[0] = -fc2->planes[ii+1][0];
								plane[1] = -fc2->planes[ii+1][1];
								plane[2] = -fc2->planes[ii+1][2];
								plane[3] = -fc2->planes[ii+1][3];
								ec->additional_planes.push_back(plane);
								// cout << endl << ec->cutting_planes[0][0] << "," << ec->cutting_planes[0][1] << "," << ec->cutting_planes[0][2] << "," << ec->cutting_planes[0][3];
								#else
								plane2 = fc2->planes[ii+1];
								ec->additional_planes.push_back(fc2->planes[ii+1]);
								// cout << endl << -ec->cutting_planes[0][0] << "," << -ec->cutting_planes[0][1] << "," << -ec->cutting_planes[0][2] << "," << -ec->cutting_planes[0][3];
								#endif
								break;
							}
						}
						if (!found)
						{
							cout << endl << ERR << "Cannot detect correct plane!";
							exit(-1);
						}
						// get the two reference vectors of the planes.
						// double vref1[3],vref2[3];

						// SUB(vref1,fc1->center,ec->center)
						// SUB(vref2,fc2->center,ec->center)

						// decide AND or TEST based on orientation
						// bool flag = orientation(vref1,vref2,ec,fc1->center,fc2->center,plane1,plane2);
						bool flag = orientation(fc1->center,fc2->center,plane1,plane2);
						#if defined(OPTIMIZE_CELL_STRUCTURE)
						flag = !flag;
						#endif
						// cout << endl << "Direction " << direction;
						// cout << endl << "Orientation " << flag;
						ec->flags.push_back(flag);

						// go the next pair using the correct direction
						currentProbe1 += direction;
						currentProbe1 = (index + currentProbe1)%index;
						currentProbe2 = (index + currentProbe1 + (direction>>1))%index;
						// cycle is restarting, stop
						if (currentProbe1 == start_index)
							break;
					}
				}
			}
			// formal "OR". Not executed in practice
			else
				ec->acute = false;

			if (bigR >= probe_radius)
			{
				// ec->isSelfIntersecting = false;
				ec->self_intersection_radius = -1.0; // the sign of the radius incubates the boolean isSelfIntersecting
			}
			else
			{
				// ec->isSelfIntersecting = true;
				ec->self_intersection_radius = sqrt(probe_radius*probe_radius - bigR*bigR);

				// if the edge is singular and self intersecting now probes pair must be collected
			}

			ec->clipping_center[0] = clipping_center[0];
			ec->clipping_center[1] = clipping_center[1];
			ec->clipping_center[2] = clipping_center[2];

			ec->clipping_radius = clipping_radius;

			#if defined(CHECK_MAX_RADII)
			max_clipping_radius = fmax(max_clipping_radius, clipping_radius);
			#endif

			ec->Rot[0][0] = v[0];
			ec->Rot[0][1] = u[0];
			ec->Rot[0][2] = w[0];
			ec->Rot[1][0] = v[1];
			ec->Rot[1][1] = u[1];
			ec->Rot[1][2] = w[1];
			ec->Rot[2][0] = v[2];
			ec->Rot[2][1] = u[2];
			ec->Rot[2][2] = w[2];

			double det2;
			INVERT_3X3(ec->invrot,det2,ec->Rot)

			#if !defined(OPTIMIZE_CELL_STRUCTURE)
			#if !defined(NEW_ATOM_PATCHES)
			atomPatches[ind1]->neighbours.push_back(ec);
			atomPatches[ind2]->neighbours.push_back(ec);
			#else
			PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[ind1] ];
			PointCell *pc2 = (PointCell *)sesComplex[ atomPatches[ind2] ];

			pc1->neighbours.push_back(ec);
			pc2->neighbours.push_back(ec);
			#endif
			#else // OPTIMIZE_CELL_STRUCTURE
			PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[ind1] ];
			PointCell *pc2 = (PointCell *)sesComplex[ atomPatches[ind2] ];

			pc1->neighbour_data.push_back(clipping_center[0]);
			pc1->neighbour_data.push_back(clipping_center[1]);
			pc1->neighbour_data.push_back(clipping_center[2]);
			pc1->neighbour_data.push_back(r2_clipping);
			pc2->neighbour_data.push_back(clipping_center[0]);
			pc2->neighbour_data.push_back(clipping_center[1]);
			pc2->neighbour_data.push_back(clipping_center[2]);
			pc2->neighbour_data.push_back(r2_clipping);
			#endif // OPTIMIZE_CELL_STRUCTURE
		}
	}
	#if defined(CHECK_MAX_RADII)
	printf (" \n Max major radius = %.3f\n", max_major_radius);
	printf (" Max clipping radius = %.3f\n", max_clipping_radius);
	#endif
}


void ConnollySurface::RemoveSelfIntersections (vector<ConnollyCell*> &sesComplex, Octree<vector<FacetCell*>> &gridProbesMap, double grid_pars[])
{
	double gxmin = grid_pars[0];
	double gymin = grid_pars[1];
	double gzmin = grid_pars[2];
	double gscale = grid_pars[6];
	int ggrid = (int)(grid_pars[7] + 1.E-7);

	double probe_radius2 = probe_radius*probe_radius;

	// remove self intersections
	for (unsigned int i=0; i<sesComplex.size(); i++)
	{
		ConnollyCell *cc = sesComplex[i];

		if (cc->patch_type != REGULAR_FACE_CELL && cc->patch_type != SINGULAR_FACE_CELL)
			continue;

		FacetCell *fc1 = (FacetCell*)cc;

		int ix = (int)rintp((fc1->center[0]-gxmin)*gscale);
		int iy = (int)rintp((fc1->center[1]-gymin)*gscale);
		int iz = (int)rintp((fc1->center[2]-gzmin)*gscale);

		// cout << endl << ix << " " << iy << " " << iz;

		for (int k=0; k<27; k++)
		{
			int cx = ix+shift_map[k][0];
			int cy = iy+shift_map[k][1];
			int cz = iz+shift_map[k][2];

			if (cx >= ggrid || cy >= ggrid || cz >= ggrid || cx < 0 || cy < 0 || cz < 0)
				continue;

			vector<FacetCell*> &vec = gridProbesMap(cx,cy,cz);

			for (int64_t j=0; j<vec.size(); j++)
			{
				// FacetCell *fc2 = SELF_MAP(cx,cy,cz,j+1,ggrid,ggrid,ggrid);
				FacetCell *fc2 = vec[j];

				if (fc1 == fc2)
					continue;

				double *c1 = fc1->center;
				double *c2 = fc2->center;
				double dist2;

				DIST2(dist2,c1,c2)
				// cout << endl << dist2 << " " << 4*probe_radius2;
				// there is a self intersection, must clip
				if (dist2 < 4*probe_radius2)
				{
					// a reference point is the mid point
					double ref_point[3],dir[3],temp;
					ADD(ref_point,c2,c1)
					VEC_MUL(ref_point,ref_point,0.5)

					fc1->isSelfIntersecting = true;
					fc2->isSelfIntersecting = true;

					// get plane direction
					SUB(dir,c1,c2)
					NORMALIZE_S(dir,temp)

					double bias = -DOT(dir,ref_point);

					double *plane = allocateVector<double>(4);
					plane[0] = dir[0];
					plane[1] = dir[1];
					plane[2] = dir[2];
					plane[3] = bias;

					double orient = DOT(fc1->center,plane);
					// ok
					if (orient < -plane[3])
					{}
					// swap
					else
					{
						plane[0] = -dir[0];
						plane[1] = -dir[1];
						plane[2] = -dir[2];
						plane[3] = -bias;
					}

					fc1->self_intersection_planes.push_back(plane);
					#if !defined(OPTIMIZE_CELL_STRUCTURE)
					double *plane2 = allocateVector<double>(4);
					plane2[0] = -plane[0];
					plane2[1] = -plane[1];
					plane2[2] = -plane[2];
					plane2[3] = -plane[3];
					fc2->self_intersection_planes.push_back(plane2);
					#else
					fc1->self_intersection_plane_labels.push_back(true);
					// the following shares the same above data but the sign is recorded via the following flag
					fc2->self_intersection_planes.push_back(plane);
					fc2->self_intersection_plane_labels.push_back(false);
					#endif
				}
			}
		}
	}
}


#ifdef MULTITHREADED_SES_BUILDING

void ConnollySurface::BuildPatches (ThreadDataWrapper *tdw, double grid_pars[]
									#if defined(OPTIMIZE_BUILDING_MEMORY)
									, TaggedDataWrapper *tagged_data_wrapper
									#endif
									)
{
	double gxmin = grid_pars[0];
	double gymin = grid_pars[1];
	double gzmin = grid_pars[2];
	double gscale = grid_pars[6];
	int ggrid = (int)(grid_pars[7] + 1.E-7);

	Octree<vector<FacetCell*>> gridProbesMap(ggrid,vector<FacetCell*>());

	Fixed_alpha_shape_3 alpha_shape(tdw->l.begin(), tdw->l.end(), 0.0);

	////////////////////////////////// atom cells ///////////////////////////////////////////////
	#if !defined(OPTIMIZE_BUILDING_MEMORY)
	BuildPointCells (alpha_shape, tdw->sesComplex, tdw->atomPatches, tdw->exposed, tdw->num_cells);
	#else
	BuildPointCells (alpha_shape, tdw->sesComplex, tdw->atomPatches, tdw->exposed, tdw->num_cells,
					 tagged_data_wrapper, tdw->my_task_id, tdw->my_thread_id, grid_pars);
	#endif

	////////////////////////////////// facet cells ///////////////////////////////////////////////
	BuildFacetCells (alpha_shape, tdw->sesComplex, gridProbesMap, tdw->atomPatches,
					 tdw->num_cells, grid_pars);

	////////////////////////////////// edge cells ///////////////////////////////////////////////
	#if !defined(OPTIMIZE_BUILDING_MEMORY)
	BuildEdgeCells (output_file, alpha_shape, tdw->sesComplex, tdw->atomPatches, tdw->num_cells);
	#else
	BuildEdgeCells (output_file, alpha_shape, tdw->sesComplex, tdw->atomPatches, tdw->num_cells,
					tagged_data_wrapper, tdw->my_task_id, tdw->my_thread_id, grid_pars);
	#endif

	////////////////////////////////// remove self intersections //////////////////////////////////////////////
	RemoveSelfIntersections (tdw->sesComplex, gridProbesMap, grid_pars);

	tdw->l.clear();
}

#endif


#ifdef MULTITHREADED_SES_BUILDING

void ConnollySurface::reorderPatchCells (DelPhiShared *delphi, ThreadDataWrapper *tdw, int cell_start, int cell_end)
{
	auto compCellSizes = [delphi] (ConnollyCell *cc1, ConnollyCell *cc2)
	{
		if (cc1->patch_type == REGULAR_FACE_CELL || cc1->patch_type == SINGULAR_FACE_CELL)
		{
			// In this case the radiuses are all the same and equal to the probe radius
			return true;
		}
		double radius1, radius2;

		if (cc1->patch_type == SINGULAR_EDGE_CELL || cc1->patch_type == REGULAR_EDGE_CELL)
		{
			radius1 = ((EdgeCell*)&cc1)->clipping_radius;
			radius2 = ((EdgeCell*)&cc2)->clipping_radius;
		}
		else if (cc1->patch_type == POINT_CELL)
		{
			radius1 = delphi->atoms[ ((PointCell*)&cc1)->id ].radius;
			radius2 = delphi->atoms[ ((PointCell*)&cc2)->id ].radius;
		}
		return radius1 < radius2;
	};

	sort (&tdw->sesComplex[cell_start], &tdw->sesComplex[cell_end], compCellSizes);
}

#endif // MULTITHREADED_SES_BUILDING


#ifdef CHECK_BUILDUP_DIFF

bool ConnollySurface::compPatchTags (ConnollyCell *cc1, ConnollyCell *cc2)
{
	return cc1->tag < cc2->tag;
}



void ConnollySurface::checkBuildupDivergences (vector<ConnollyCell*> &st_sesComplex, vector<ConnollyCell*> &sesComplex,
											   int local_divergent_values[16][9], int thread_id, int num_threads)
{
	for (int n = thread_id; n < sesComplex.size(); n += num_threads)
	{
		ConnollyCell *cc = sesComplex[n];

		FacetCell *fc, *st_fc;
		EdgeCell *ec, *st_ec;
		PointCell *pc, *st_pc;

		double *sphere_center, *st_sphere_center;
		double radius, st_radius;

		int found = 0;

		if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
		{
			found = -1;

			fc = (FacetCell*)cc;
			sphere_center = fc->center;
			radius = probe_radius;
		}
		else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
		{
			found = -2;

			ec = (EdgeCell*)cc;
			sphere_center = ec->clipping_center;
			radius = ec->clipping_radius;
		}
		else if (cc->patch_type == POINT_CELL)
		{
			found = -3;

			pc = (PointCell*)cc;
			radius = delphi->atoms[pc->id].radius;
			sphere_center = delphi->atoms[pc->id].pos;
		}

		for (int m = 0; m < st_sesComplex.size(); m++)
		{
			ConnollyCell *st_cc = st_sesComplex[m];

			if (st_cc->patch_type == REGULAR_FACE_CELL || st_cc->patch_type == SINGULAR_FACE_CELL)
			{
				st_fc = (FacetCell*)st_cc;

				if (st_cc->patch_type == cc->patch_type && (
					(st_fc->id[0] == fc->id[0] && st_fc->id[1] == fc->id[1] && st_fc->id[2] == fc->id[2]) ||
					(st_fc->id[0] == fc->id[0] && st_fc->id[1] == fc->id[2] && st_fc->id[2] == fc->id[1]) ||
					(st_fc->id[0] == fc->id[1] && st_fc->id[1] == fc->id[2] && st_fc->id[2] == fc->id[0]) ||
					(st_fc->id[0] == fc->id[1] && st_fc->id[1] == fc->id[0] && st_fc->id[2] == fc->id[2]) ||
					(st_fc->id[0] == fc->id[2] && st_fc->id[1] == fc->id[0] && st_fc->id[2] == fc->id[1]) ||
					(st_fc->id[0] == fc->id[2] && st_fc->id[1] == fc->id[1] && st_fc->id[2] == fc->id[0])))
				{
					found = 1;

					st_sphere_center = st_fc->center;
					st_radius = probe_radius;
				}
			}
			else if (st_cc->patch_type == SINGULAR_EDGE_CELL || st_cc->patch_type == REGULAR_EDGE_CELL)
			{
				st_ec = (EdgeCell*)st_cc;

				if (st_cc->patch_type == cc->patch_type && (
					(st_ec->id[0] == ec->id[0] && st_ec->id[1] == ec->id[1]) ||
					(st_ec->id[0] == ec->id[1] && st_ec->id[1] == ec->id[0])))
				{
					found = 2;

					st_sphere_center = st_ec->clipping_center;
					st_radius = st_ec->clipping_radius;
				}
			}
			else if (st_cc->patch_type == POINT_CELL)
			{
				st_pc = (PointCell*)st_cc;

				if (st_cc->patch_type == cc->patch_type &&
					st_pc->id == pc->id)
				{
					found = 3;

					st_radius = delphi->atoms[st_pc->id].radius;
					st_sphere_center = delphi->atoms[st_pc->id].pos;
				}
			}
		}
		if (found <= 0)
		{
			if (found == -1)
				++local_divergent_values[thread_id][0];
			else if (found == -3)
				++local_divergent_values[thread_id][1];
			else
				++local_divergent_values[thread_id][2];
			continue;
		}
		if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
		{
			if (fabs((st_fc->center[0]-fc->center[0])/st_fc->center[0]) > EPS_BUILDUP ||
				fabs((st_fc->center[1]-fc->center[1])/st_fc->center[1]) > EPS_BUILDUP ||
				fabs((st_fc->center[2]-fc->center[2])/st_fc->center[2]) > EPS_BUILDUP)
			{
				// This relies on the assumption that patches matches if and only if the centers are also identical
				found = -1;
				continue;
			}

			for (int j=0; j<4; j++)
			{
				for (int i=0; i<4; i++)
				{
					double rel_diff = fabs((st_fc->planes[j][i]-fc->planes[j][i])/st_fc->planes[j][i]);

					if (rel_diff > EPS_BUILDUP)
					{
						printf ("\n plane[j][i] rel. diff. = %.10e\n", rel_diff);
						printf (" values: %.10e, %.10e\n", st_fc->planes[j][i], fc->planes[j][i]);
						++local_divergent_values[thread_id][3];
					}
				}
			}
		}
		else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
		{
			if (fabs((st_ec->center[0]-ec->center[0])/st_ec->center[0]) > EPS_BUILDUP ||
				fabs((st_ec->center[1]-ec->center[1])/st_ec->center[1]) > EPS_BUILDUP ||
				fabs((st_ec->center[2]-ec->center[2])/st_ec->center[2]) > EPS_BUILDUP)
			{
				// This relies on the assumption that patches matches if and only if the centers are also identical
				found = -2;
				continue;
			}

			for (int j=0; j<3; j++)
			{
				for (int i=0; i<3; i++)
				{
					double rel_diff = fabs((fabs(st_ec->invrot[j][i])-fabs(ec->invrot[j][i]))/st_ec->invrot[j][i]);

					if (rel_diff > EPS_BUILDUP)
					{
						printf ("\n invrot[j][i] rel. diff. = %.10e\n", rel_diff);
						printf (" abs values: %.10e, %.10e\n", fabs(st_ec->invrot[j][i]), fabs(ec->invrot[j][i]));
						++local_divergent_values[thread_id][4];
					}
				}
			}
			double rel_diff = fabs((st_ec->clipping_center[0]-ec->clipping_center[0])/st_ec->clipping_center[0]);

			if (rel_diff > EPS_BUILDUP)
			{
				printf ("\n clipping center[0] rel. diff. = %.10e\n", rel_diff);
				++local_divergent_values[thread_id][5];
			}
			rel_diff = fabs((st_ec->clipping_center[1]-ec->clipping_center[1])/st_ec->clipping_center[1]);

			if (rel_diff > EPS_BUILDUP)
			{
				printf ("\n clipping center[1] rel. diff. = %.10e\n", rel_diff);
				++local_divergent_values[thread_id][5];
			}
			rel_diff = fabs((st_ec->clipping_center[2]-ec->clipping_center[2])/st_ec->clipping_center[2]);

			if (rel_diff > EPS_BUILDUP)
			{
				printf ("\n clipping center[2] rel. diff. = %.10e\n", rel_diff);
				++local_divergent_values[thread_id][5];
			}
			rel_diff = fabs((st_ec->clipping_radius-ec->clipping_radius)/st_ec->clipping_radius);

			if (rel_diff > EPS_BUILDUP)
			{
				printf (" clipping radius rel. diff. = %.10e\n", rel_diff);
				++local_divergent_values[thread_id][6];
			}
			rel_diff = fabs((st_ec->major_radius-ec->major_radius)/st_ec->major_radius);

			if (rel_diff > EPS_BUILDUP)
			{
				printf ("\n major radius rel. diff. = %.10e\n", rel_diff);
				++local_divergent_values[thread_id][7];
			}
			rel_diff = fabs((st_ec->rcoi-ec->rcoi)/st_ec->rcoi);

			if (rel_diff > EPS_BUILDUP)
			{
				printf ("\n rcoi rel. diff. = %.10e\n", rel_diff);
				++local_divergent_values[thread_id][8];
			}
		}
	}
}


#endif // CHECK_BUILDUP_DIFF


bool ConnollySurface::buildConnollyCGAL()
{
	auto chrono_start = chrono::high_resolution_clock::now();

	char outputFile[BUFLEN];
	strcpy(outputFile, "connolly.pov");

	double x, y, z, r;
	double max_x=-1e6, min_x=1e6, max_y=-1e6, min_y=1e6, max_z=-1e6, min_z=1e6;

	// setup self intersections grid of the probe
	double gside;
	double gscale = 0.5/probe_radius;
	
	int ggrid = 1;
	
	cout << endl << INFO << "Adjusting self intersection grid ";

	while (1)
	{
		ggrid = (int)floor(gscale*si_perfil*delphi->rmaxdim);
		// cout << ggrid << ",";
		
		if (ggrid <= 100)
		{
			if (ggrid <= 0)
				gscale *= 2.0;
			else
				break;
		}
		else
		{
			gscale *= 0.5;
		}
	}
	
	cout << endl << INFO << "Self intersection grid is (before) " << ggrid;

	int p2 = 1;
	// get nearest power of two
	while (1)
	{
		p2 <<= 1;
		if (p2 > ggrid)
		{
			// the only safe thing in powers of two
			// not so efficient but it works
			ggrid = p2 >> 1;
			break;
		}
	}
	gscale = (double)ggrid/(si_perfil*delphi->rmaxdim);
	gside = (si_perfil*delphi->rmaxdim)/(double)ggrid;
	
	cout << endl << INFO << "Self intersection grid is " << ggrid;
	double gxmin = delphi->baricenter[0] - (ggrid-1)*0.5*gside;
	double gymin = delphi->baricenter[1] - (ggrid-1)*0.5*gside;
	double gzmin = delphi->baricenter[2] - (ggrid-1)*0.5*gside;
	
	double gxmax = delphi->baricenter[0] + (ggrid-1)*0.5*gside;
	double gymax = delphi->baricenter[1] + (ggrid-1)*0.5*gside;
	double gzmax = delphi->baricenter[2] + (ggrid-1)*0.5*gside;

	#if !defined(MULTITHREADED_SES_BUILDING)
	numThreadDataWrappers = 1;
	#else
	numThreadDataWrappers = min(conf.numThreads, (int)delphi->atoms.size());

	double halo_layer_size;

	if (!forceSerialBuild)
	{
		if (!conf.parallelBuildupHaloThicknessInput)
		{
			// Task-specific, Z-slab dependent assignments with halo layers of 9-18 Angstrom

			if (probe_radius <= 2.)
				halo_layer_size = 12.0;
			else if (probe_radius > 2 && probe_radius <= 3)
				halo_layer_size = 11;
			else if (probe_radius > 3 && probe_radius <= 6)
				halo_layer_size = 18;
			else
			{
				// force 1 thread for debug
				numThreadDataWrappers = 1;

				cout << endl << WARN << "The number of threads for building the surface has been reduced to 1 for probes > 6 Angstroms";
			}
		}
		else
		{
			halo_layer_size = conf.parallelBuildupHaloThickness;
		}
	}
	else
	{
		numThreadDataWrappers = 1;
	}
	#endif


	for (int task_id = 0; task_id < numTasks; task_id++)
	{
		// This is the intial setting; it might be changed later if there are too few atoms ...
		numThreadDataWrappersPerTask[task_id] = numThreadDataWrappers;
	}

	#ifdef MULTITHREADED_SES_BUILDING

	// thread_data_wrapper = new ThreadDataWrapper[ numTasks * numThreadDataWrappers ]();

	// here there is a code that subdivides the domain along the Z-axis in N=numTasks slabs
	// and each of these slabs in M=numThreads=numThreadDataWrappers mini-slabs along the Y-axis;
	// the thickness of the slabs and mini-slabs varies according to the population of atoms;
	// a large number of bins/slabs (numTasks/numThreads * bin_factor) stores the number of
	// atoms in each bin
	double slab_max_z[ numTasks ];
	double minislab_max_y[ numTasks ][ numThreadDataWrappers ];


	int bin_factor = 1000;
	int num_bins = numTasks * bin_factor;
	int num_minibins = conf.numThreads * bin_factor;
	
	int tagged_task_in_bin[ num_bins ];
	int tagged_thread_in_minibin[ numTasks ][ num_minibins ];

	int slab_atoms[numTasks];

	double bin_slab_size     = (delphi->zmax - delphi->zmin + delphi->side) / num_bins;
	double inv_bin_slab_size = num_bins / (delphi->zmax - delphi->zmin + delphi->side);

	double minibin_slab_size     = (delphi->ymax - delphi->ymin + delphi->side) / num_minibins;
	double inv_minibin_slab_size = num_minibins / (delphi->ymax - delphi->ymin + delphi->side);


	for (int task_id = 0; task_id < numTasks; task_id++)
	{
		slab_max_z[ task_id ] = delphi->zmax + delphi->hside;

		for (int thread_id = 0; thread_id < numThreadDataWrappers; thread_id++)
			minislab_max_y[ task_id ][ thread_id ] = delphi->ymax + delphi->hside;
	}

	#endif // MULTITHREADED_SES_BUILDING


	#if !defined(MULTITHREADED_SES_BUILDING) || defined(CHECK_BUILDUP_DIFF)
	#if !defined(NO_CGAL_PATCHING)
	list<Weighted_point> l;
	#else
	// nopatch
	list<std::pair<Weighted_point, int>> l;
	#endif
	#endif

	// for (int i=0; i<delphi->numAtoms; i++)
	for (int i=0; i<delphi->atoms.size(); i++)
	{
		delphi->atoms[i].pos[0] += randDisplacement*(randnum()-0.5);
		delphi->atoms[i].pos[1] += randDisplacement*(randnum()-0.5);
		delphi->atoms[i].pos[2] += randDisplacement*(randnum()-0.5);

		x = delphi->atoms[i].pos[0];
		y = delphi->atoms[i].pos[1];
		z = delphi->atoms[i].pos[2];
		r = delphi->atoms[i].radius;

		max_x = max(max_x,x+r);
		max_y = max(max_y,y+r);
		max_z = max(max_z,z+r);

		min_x = std::min(min_x,x-r);
		min_y = std::min(min_y,y-r);
		min_z = std::min(min_z,z-r);

		#if !defined(MULTITHREADED_SES_BUILDING) || defined(CHECK_BUILDUP_DIFF)
		r += probe_radius;
		#if !defined(NO_CGAL_PATCHING)
		l.push_front(Weighted_point(Point3(x,y,z), (r*r), i));
		#else
		// nopatch
		l.push_front(std::make_pair(Weighted_point(Point3(x,y,z), (r*r)), i));
		#endif
		#endif
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



	#ifdef MULTITHREADED_SES_BUILDING

	if (numTasks > 1)
	{
		int num_bin_atoms[num_bins];
		
		
		memset(num_bin_atoms, 0, num_bins * sizeof(int));
		
		// Count the number of atoms in each bin
		for (int i=0; i<delphi->atoms.size(); i++)
		{
			// count no. of atoms in each task bin
			int bin_id = (int)((z - (delphi->zmin - delphi->hside)) * inv_bin_slab_size);
			++num_bin_atoms[bin_id];
		}

		int target_num_atoms_per_task = (int)((double)delphi->atoms.size() / (double)numTasks);
		int current_atoms_sum = 0;
		int previous_atoms_sum = 0;
		int task_id = 0;
		
		// Determine the edge of each task slab
		for (int k=0; k < num_bins; k++)
		{
			current_atoms_sum += num_bin_atoms[k];

			double upper_slab_edge = (delphi->zmin - delphi->hside) + (k+1) * bin_slab_size;

			tagged_task_in_bin[k] = task_id;

			if (task_id > 0)
			{
				// be sure that the current slab's thickness is at least halo_layer_size
				if (upper_slab_edge - slab_max_z[ task_id-1 ] < halo_layer_size)
					continue;
			}
			if (current_atoms_sum - previous_atoms_sum >= target_num_atoms_per_task)
			{
				previous_atoms_sum = current_atoms_sum;
				slab_max_z[task_id] = upper_slab_edge;
				++task_id;
			}
			if (task_id == numTasks - 1)
			{
				slab_max_z[ task_id ] = delphi->zmax + delphi->hside;

				for (++k; k < num_bins; k++)
					tagged_task_in_bin[k] = task_id;
				break;
			}
		}
		if (slab_max_z[ numTasks-1 ] - slab_max_z[ numTasks-2 ] < halo_layer_size)
		{
			cout << endl << ERR << "Too small subdomain for multi-task cells' building!";
			return false;
		}

		// code to determine the atoms' population across slabs
		for (int task_id=0; task_id<numTasks; task_id++)
			slab_atoms[task_id] = 0;

		for (int i=0; i<delphi->atoms.size(); i++)
		{
			int bin_id = (int)((delphi->atoms[i].pos[2] - (delphi->zmin - delphi->hside)) * inv_bin_slab_size);
			++slab_atoms[ tagged_task_in_bin[bin_id] ];
		}
	}
	else
	{
		for (int k=0; k < num_bins; k++)
			tagged_task_in_bin[k] = 0;

		slab_atoms[0] = delphi->atoms.size();
	}

	if (numThreadDataWrappers > 1)
	{
		int num_minibin_atoms[ numTasks*num_minibins ];

		memset(num_minibin_atoms, 0, numTasks*num_minibins * sizeof(int));

		// Count the number of atoms in each minibin
		for (int i=0; i<delphi->atoms.size(); i++)
		{
			// First, determine the task bin and slab
			int bin_id = (int)((delphi->atoms[i].pos[2] - (delphi->zmin - delphi->hside)) * inv_bin_slab_size);
			int task_id = tagged_task_in_bin[bin_id];

			// Second, determine the thread minibin and count no. of atoms in each minibin
			int minibin_id = (int)((delphi->atoms[i].pos[1] - (delphi->ymin - delphi->hside)) * inv_minibin_slab_size);
			++num_minibin_atoms[ task_id*num_minibins + minibin_id ];
		}

		for (int task_id=0; task_id < numTasks; task_id++)
		{
			int target_num_atoms_per_thread = (int)((double)slab_atoms[task_id] / (double)numThreadDataWrappers);

			int current_atoms_sum = 0;
			int previous_atoms_sum = 0;

			int thread_id = 0;
			
			// Determine the edge of each thread minislab
			for (int j=0; j < num_minibins; j++)
			{
				current_atoms_sum += num_minibin_atoms[ task_id*num_minibins + j ];

				double right_minislab_edge = (delphi->ymin - delphi->hside) + (j+1) * minibin_slab_size;

				tagged_thread_in_minibin[ task_id ][ j ] = thread_id;

				if (thread_id > 0)
				{
					// be sure that the current minislab's thickness is at least halo_layer_size
					if (right_minislab_edge - minislab_max_y[ task_id ][ thread_id-1 ] < halo_layer_size)
						continue;
				}
				if (current_atoms_sum - previous_atoms_sum >= target_num_atoms_per_thread)
				{
					previous_atoms_sum = current_atoms_sum;
					minislab_max_y[ task_id ][ thread_id ] = right_minislab_edge;
					++thread_id;
				}
				if (thread_id == numThreadDataWrappers - 1)
				{
					minislab_max_y[ task_id ][ thread_id ] = delphi->ymax + delphi->hside;

					for (++j; j < num_minibins; j++)
						tagged_thread_in_minibin[ task_id ][j] = thread_id;
					break;
				}
			}
			int max_numThreadDataWrappers = numThreadDataWrappers;
			int last_j;

			for (int j=0; j < numThreadDataWrappers-1; j++)
			{
				last_j = j;

				if (minislab_max_y[task_id][ j+1 ] - minislab_max_y[task_id][ j ] < halo_layer_size)
				{
					// the slab is enlarged so that its border will be the grid border
					minislab_max_y[task_id][ j ] = delphi->ymax + delphi->hside;

					if (max_numThreadDataWrappers > j+1)
					{
						max_numThreadDataWrappers = j+1;
						break;
					}
				}
			}
			for (int j = last_j + 1; j < numThreadDataWrappers; j++)
				minislab_max_y[task_id][ j ] = delphi->ymax + delphi->hside;

			if (max_numThreadDataWrappers < numThreadDataWrappers)
			{
				for (int j=0; j < num_minibins; j++)
				{
					if (tagged_thread_in_minibin[ task_id ][ j ] > max_numThreadDataWrappers - 1)
						tagged_thread_in_minibin[ task_id ][ j ] = max_numThreadDataWrappers - 1;
				}
				numThreadDataWrappersPerTask[task_id] = max_numThreadDataWrappers;

				cout << endl << WARN << "Imposed limitation: number of threads in parallel build-up stage: " << max_numThreadDataWrappers;
			}
		}
	}
	else
	{
		for (int task_id=0; task_id < numTasks; task_id++)
		{
			for (int k=0; k < num_minibins; k++)
				tagged_thread_in_minibin[ task_id ][ k ] = 0;
		}
	}

	double slab_pars[] = {halo_layer_size, inv_bin_slab_size, inv_minibin_slab_size};

	TaggedDataWrapper tagged_data_wrapper;

	tagged_data_wrapper.slab_max_z = slab_max_z;
	tagged_data_wrapper.task_in_bin = tagged_task_in_bin;

	tagged_data_wrapper.minislab_max_y = allocateVector<double*>(numTasks);
	tagged_data_wrapper.thread_in_minibin = allocateVector<int*>(numTasks);

	for (int task_id=0; task_id<numTasks; task_id++)
	{
		tagged_data_wrapper.minislab_max_y[ task_id ] = minislab_max_y[ task_id ];
		tagged_data_wrapper.thread_in_minibin[ task_id ] = tagged_thread_in_minibin[ task_id ];
	}

	int thread_pars[] = {0, (int)delphi->atoms.size(), 1};
	BuildWeightedPoints (&tagged_data_wrapper, slab_pars, thread_pars);

	/*
	// the bellow parallel version does not perfectly reproduce the serial version when BOOST is used, for some reasons
	// (recall that now in BuildWeightedPoints (...) the line including the mutex is commented because the above call
	//  does not use BOOST)

	#ifdef ENABLE_BOOST_THREADS
	boost::thread_group thdGroup;
	#endif

	for (int thd_id=0; thd_id<numThreadDataWrappers; thd_id++)
	{
		int thread_pars[] = {thd_id, (int)delphi->atoms.size(), numThreadDataWrappers};

		#if defined(ENABLE_BOOST_THREADS)
		thdGroup.create_thread(boost::bind(&ConnollySurface::BuildWeightedPoints, this, &tagged_data_wrapper, slab_pars, thread_pars));
		#else
		BuildWeightedPoints (&tagged_data_wrapper, slab_pars, thread_pars);
		#endif
	}
	#ifdef ENABLE_BOOST_THREADS
	thdGroup.join_all();
	#endif
	*/


	cout << endl << INFO << "Allocating self intersection grid and computing alpha shape complex ...";
	cout.flush();

	vector<int> exposed;
	exposed.reserve(10);


	// Allocate buffers for each task and thread
	for (int task_id=0; task_id<numTasks; task_id++)
	{
		for (int thd_id=0; thd_id<numThreadDataWrappersPerTask[task_id]; thd_id++)
		{
			ThreadDataWrapper *tdw = &thread_data_wrapper[task_id*numThreadDataWrappers+thd_id];

			tdw->sesComplex.reserve( MAX(10, delphi->atoms.size()*4/(2.*numTasks*numThreadDataWrappers)) );

			#if !defined(NEW_ATOM_PATCHES)
			tdw->atomPatches = allocateVector<PointCell*>(delphi->atoms.size());

			for (int i=0; i<delphi->atoms.size(); i++)
				tdw->atomPatches[i] = NULL;
			#else
			tdw->atomPatches = allocateVector<int>(delphi->atoms.size());
			#endif
			
			for (int i=0; i<5; i++)
				tdw->num_cells[i] = 0;

			#if defined(OPTIMIZE_BUILDING_MEMORY)
			tdw->my_task_id   = task_id;
			tdw->my_thread_id = thd_id;
			#endif
		}
	}

	double grid_pars[] = {gxmin, gymin, gzmin, gxmax, gymax, gzmax, gscale, (double)ggrid};

	for (int task_id=0; task_id<numTasks; task_id++)
	{
		#ifdef ENABLE_BOOST_THREADS
		boost::thread_group thdGroup;
		#endif
		for (int thd_id=0; thd_id<numThreadDataWrappersPerTask[task_id]; thd_id++)
		{
			ThreadDataWrapper *tdw = &thread_data_wrapper[task_id*numThreadDataWrappers+thd_id];

			#if !defined(OPTIMIZE_BUILDING_MEMORY)
			#if defined(ENABLE_BOOST_THREADS)
			thdGroup.create_thread(boost::bind(&ConnollySurface::BuildPatches, this,
											   tdw, grid_pars));
			#else
			BuildPatches (tdw, grid_pars);
			#endif

			#else // OPTIMIZE_BUILDING_MEMORY

			#if defined(ENABLE_BOOST_THREADS)
			thdGroup.create_thread(boost::bind(&ConnollySurface::BuildPatches, this,
											   tdw, grid_pars, &tagged_data_wrapper));
			#else
			BuildPatches (tdw, grid_pars, &tagged_data_wrapper);
			#endif
			#endif // OPTIMIZE_BUILDING_MEMORY
		}
		#ifdef ENABLE_BOOST_THREADS
		thdGroup.join_all();
		#endif
	}

	/*
	for (int task_id=0; task_id<numTasks; task_id++)
	{
		for (int thd_id=0; thd_id<numThreadDataWrappersPerTask[task_id]; thd_id++)
		{
			ThreadDataWrapper *tdw = &thread_data_wrapper[task_id*numThreadDataWrappers+thd_id];

			reorderPatchCells (delphi, tdw, 0, tdw->sesComplex.size());
		}
	}
	*/

	deleteVector<int*>(tagged_data_wrapper.thread_in_minibin);
	deleteVector<double*>(tagged_data_wrapper.minislab_max_y);


	sesComplex.reserve(delphi->atoms.size()*4);

	// The following removes dubplicated patch cells

	#if !defined(OPTIMIZE_BUILDING_MEMORY)
	for (int task_id=0; task_id<numTasks; task_id++)
	{
		for (int thread_id=0; thread_id<numThreadDataWrappersPerTask[task_id]; thread_id++)
		{
			ThreadDataWrapper *tdw = &thread_data_wrapper[task_id*numThreadDataWrappers+thread_id];

			for (int m = 0; m < tdw->sesComplex.size(); m++)
			{
				ConnollyCell *cc = tdw->sesComplex[m];
				
				double *sphere_center;
				double radius;
				
				if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
				{
					FacetCell *fc = (FacetCell*)cc;
					sphere_center = fc->center;
					radius = probe_radius;
				}
				else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
				{
					EdgeCell *ec = (EdgeCell*)cc;
					sphere_center = ec->clipping_center;
					radius = ec->clipping_radius;
				}
				else if (cc->patch_type == POINT_CELL)
				{
					PointCell *pc = (PointCell*)cc;
					radius = delphi->atoms[pc->id].radius;
					sphere_center = delphi->atoms[pc->id].pos;
				}
				// if the bottom border of the patch resides within the slab of another task the cell is discarded ...
				double downz = sphere_center[2] - radius;
				double slab_border;
				
				// lower slab border
				if (task_id == 0)
					slab_border = delphi->zmin - delphi->hside;
				else
					slab_border = slab_max_z[ task_id-1 ];
				// ... in particular,
				// check if the lower border of the patch is inside another slab;
				// if so, the patch will be processed by the task responsible for the neighbour slab
				if (downz < slab_border)
				{
					--tdw->num_cells[ cc->patch_type ];
					cc->patch_type = SKIP_CELL;
					continue;
				}
				// If the lower border of the patch is above the upper border of the current slab
				// then the patch will be processed by the task responsible for the upper slab
				if (downz > slab_max_z[ task_id ])
				{
					--tdw->num_cells[ cc->patch_type ];
					cc->patch_type = SKIP_CELL;
					continue;
				}
				
				// if the left border of the patch resides within the minislab of another thread the cell is discarded ...
				double downy = sphere_center[1] - radius;
				double minislab_border;
				
				// left slab border
				if (thread_id == 0)
					minislab_border = delphi->ymin - delphi->hside;
				else
					minislab_border = minislab_max_y[ task_id ][ thread_id-1 ];
				// ... in particular,
				// check if the left border of the patch is inside another minislab;
				// if so, the patch will be processed by the thread responsible for the neighbour slab
				if (downy < minislab_border)
				{
					--tdw->num_cells[ cc->patch_type ];
					cc->patch_type = SKIP_CELL;
					continue;
				}
				// If the lower border of the patch is above the upper border of the current slab
				// then the patch will be processed by the task responsible for the upper slab
				if (downy > minislab_max_y[ task_id ][ thread_id ])
				{
					--tdw->num_cells[ cc->patch_type ];
					cc->patch_type = SKIP_CELL;
					continue;
				}
				sesComplex.push_back(cc);
			}
		}
	}
	#else // OPTIMIZE_BUILDING_MEMORY
	for (int task_id=0; task_id<numTasks; task_id++)
	{
		for (int thread_id=0; thread_id<numThreadDataWrappersPerTask[task_id]; thread_id++)
		{
			ThreadDataWrapper *tdw = &thread_data_wrapper[task_id*numThreadDataWrappers+thread_id];

			for (int m = 0; m < tdw->sesComplex.size(); m++)
			{
				ConnollyCell *cc = tdw->sesComplex[m];

				double *sphere_center;
				double radius;

				if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
				{
					FacetCell *fc = (FacetCell*)cc;
					sphere_center = fc->center;
					radius = probe_radius;
				}
				else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
				{
					EdgeCell *ec = (EdgeCell*)cc;
					sphere_center = ec->clipping_center;
					radius = ec->clipping_radius;
				}
				else if (cc->patch_type == POINT_CELL)
				{
					PointCell *pc = (PointCell*)cc;
					radius = delphi->atoms[pc->id].radius;
					sphere_center = delphi->atoms[pc->id].pos;
				}
				// if the bottom border of the patch resides within the slab of another task the cell is discarded ...
				double downz = sphere_center[2] - radius;
				double slab_border;

				// lower slab border
				if (task_id == 0)
					slab_border = delphi->zmin - delphi->hside;
				else
					slab_border = slab_max_z[ task_id-1 ];
				// ... in particular,
				// check if the lower border of the patch is inside another slab;
				// if so, the patch will be processed by the task responsible for the neighbour slab
				if (downz < slab_border)
				{
					--tdw->num_cells[ cc->patch_type ];
					cc->patch_type = SKIP_CELL;
					continue;
				}
				// If the lower border of the patch is above the upper border of the current slab
				// then the patch will be processed by the task responsible for the upper slab
				if (downz > slab_max_z[ task_id ])
				{
					--tdw->num_cells[ cc->patch_type ];
					cc->patch_type = SKIP_CELL;
					continue;
				}

				// if the left border of the patch resides within the minislab of another thread the cell is discarded ...
				double downy = sphere_center[1] - radius;
				double minislab_border;

				// left slab border
				if (thread_id == 0)
					minislab_border = gymin;
				else
					minislab_border = minislab_max_y[ task_id ][ thread_id-1 ];
				// ... in particular,
				// check if the left border of the patch is inside another minislab;
				// if so, the patch will be processed by the thread responsible for the neighbour slab
				if (downy < minislab_border)
				{
					--tdw->num_cells[ cc->patch_type ];
					cc->patch_type = SKIP_CELL;
					continue;
				}
				// If the lower border of the patch is above the upper border of the current slab
				// then the patch will be processed by the task responsible for the upper slab
				if (downy > minislab_max_y[ task_id ][ thread_id ])
				{
					--tdw->num_cells[ cc->patch_type ];
					cc->patch_type = SKIP_CELL;
					continue;
				}
				sesComplex.push_back(cc);
			}
		}
	}
	#endif // OPTIMIZE_BUILDING_MEMORY

	int num_cells[5] = {0, 0, 0, 0, 0};

	for (int task_id=0; task_id<numTasks; task_id++)
	{
		for (int thd_id=0; thd_id<numThreadDataWrappersPerTask[task_id]; thd_id++)
		{
			for (int i=0; i<5; i++)
				num_cells[i] += thread_data_wrapper[task_id*numThreadDataWrappers+thd_id].num_cells[i];
			
			for (int i=0; i<thread_data_wrapper[task_id*numThreadDataWrappers+thd_id].exposed.size(); i++)
				exposed.push_back( thread_data_wrapper[task_id*numThreadDataWrappers+thd_id].exposed[i] );
			
			thread_data_wrapper[task_id*numThreadDataWrappers+thd_id].exposed.clear();
		}
	}

	#endif // MULTITHREADED_SES_BUILDING


	#if !defined(MULTITHREADED_SES_BUILDING)
	cout << endl << INFO << "Allocating self intersection grid and computing alpha shape complex....";
	cout.flush();

	vector<int> exposed;
	exposed.reserve(10);

	sesComplex.reserve(delphi->atoms.size()*4);

	#if !defined(NEW_ATOM_PATCHES)
	atomPatches = allocateVector<PointCell*>(delphi->atoms.size());

	for (int i=0; i<delphi->atoms.size(); i++)
		atomPatches[i] = NULL;
	#else
	atomPatches = allocateVector<int>(delphi->atoms.size());
	#endif
	#endif // MULTITHREADED_SES_BUILDING


	#if defined(CHECK_BUILDUP_DIFF)
	vector<int> st_exposed;
	st_exposed.reserve(10);

	vector<ConnollyCell*> st_sesComplex;
	#endif


	#if !defined(MULTITHREADED_SES_BUILDING) || defined(CHECK_BUILDUP_DIFF)
	Octree<vector<FacetCell*>> gridProbesMap(ggrid,vector<FacetCell*>());

	Fixed_alpha_shape_3 alpha_shape(l.begin(), l.end(), 0.0);
	#endif


	#if !defined(MULTITHREADED_SES_BUILDING)

	int num_cells[] = {0, 0, 0, 0, 0};

	////////////////////////////////// atom cells ///////////////////////////////////////////////
	BuildPointCells (alpha_shape, sesComplex, atomPatches, exposed, num_cells);

	////////////////////////////////// facet cells ///////////////////////////////////////////////
	BuildFacetCells (alpha_shape, sesComplex, gridProbesMap, atomPatches,
					 num_cells, gxmin, gymin, gzmin, gscale, ggrid);

	////////////////////////////////// edge cells ///////////////////////////////////////////////
	BuildEdgeCells (output_file, alpha_shape, sesComplex, atomPatches, num_cells);

	////////////////////////////////// remove self intersections //////////////////////////////////////////////
	RemoveSelfIntersections (sesComplex, gridProbesMap, grid_pars);

	#endif // MULTITHREADED_SES_BUILDING


	#if defined(CHECK_BUILDUP_DIFF)

	st_sesComplex.reserve(delphi->atoms.size()*4);

	#if !defined(NEW_ATOM_PATCHES)
	st_atomPatches = allocateVector<PointCell*>(delphi->atoms.size());

	for (int i=0; i<delphi->atoms.size(); i++)
		st_atomPatches[i] = NULL;
	#else
	st_atomPatches = allocateVector<int>(delphi->atoms.size());
	#endif

	int st_num_cells[] = {0, 0, 0, 0, 0};

	////////////////////////////////// atom cells ///////////////////////////////////////////////
	BuildPointCells (alpha_shape, st_sesComplex, st_atomPatches, st_exposed, st_num_cells);

	////////////////////////////////// facet cells ///////////////////////////////////////////////
	BuildFacetCells (alpha_shape, st_sesComplex, gridProbesMap, st_atomPatches,
					 st_num_cells, grid_pars);

	////////////////////////////////// edge cells ///////////////////////////////////////////////
	BuildEdgeCells (output_file, alpha_shape, st_sesComplex, st_atomPatches, st_num_cells);

	////////////////////////////////// remove self intersections //////////////////////////////////////////////
	RemoveSelfIntersections (st_sesComplex, gridProbesMap, grid_pars);


	int divergent_face_patches     = 0;
	int divergent_edge_patches     = 0;
	int divergent_point_patches    = 0;
	int divergent_planes           = 0;
	int divergent_invrot_matrices  = 0;
	int digergent_clipping_centers = 0;
	int digergent_clipping_radii   = 0;
	int divergent_major_radii      = 0;
	int divergent_rcoi             = 0;

	int local_divergent_values[16][9];

	for (int thd_id=0; thd_id<16; thd_id++)
	{
		for (int m=0; m<9; m++)
		{
			local_divergent_values[thd_id][m] = 0;
		}
	}

	#ifdef ENABLE_BOOST_THREADS
	boost::thread_group thdGroup;
	#endif
	for (int thd_id=0; thd_id<16; thd_id++)
	{
		#if defined(ENABLE_BOOST_THREADS)
		thdGroup.create_thread(boost::bind(&ConnollySurface::checkBuildupDivergences, this,
										   st_sesComplex, sesComplex, local_divergent_values, thd_id, 16));
		#else
		checkBuildupDivergences (st_sesComplex, sesComplex, local_divergent_values, thd_id, 16);
		#endif
	}
	#ifdef ENABLE_BOOST_THREADS
	thdGroup.join_all();
	#endif

	for (int thd_id=0; thd_id<16; thd_id++)
	{
		divergent_face_patches 		+= local_divergent_values[thd_id][0];
		divergent_edge_patches 		+= local_divergent_values[thd_id][1];
		divergent_point_patches 	+= local_divergent_values[thd_id][2];
		divergent_planes 			+= local_divergent_values[thd_id][3];
		divergent_invrot_matrices 	+= local_divergent_values[thd_id][4];
		digergent_clipping_centers 	+= local_divergent_values[thd_id][5];
		digergent_clipping_radii 	+= local_divergent_values[thd_id][6];
		divergent_major_radii 		+= local_divergent_values[thd_id][7];
		divergent_rcoi 				+= local_divergent_values[thd_id][8];
	}
	int global_counts[] = {
		st_num_cells[REGULAR_FACE_CELL] + st_num_cells[SINGULAR_FACE_CELL],
		st_num_cells[REGULAR_EDGE_CELL] + st_num_cells[SINGULAR_EDGE_CELL],
		st_num_cells[POINT_CELL]};

	printf ("\n");
	printf ("Divergent face patches = %i out of %i\n"				, divergent_face_patches, global_counts[0]);
	printf ("Divergent edge patches = %i out of %i\n"				, divergent_edge_patches, global_counts[1]);
	printf ("Divergent point patches = %i out of %i\n"				, divergent_point_patches, global_counts[2]);
	printf ("Divergent planes = %i out of %i\n"						, divergent_planes, global_counts[0] * 16);
	printf ("Divergent invrot matrices' elements = %i out of %i\n"	, divergent_invrot_matrices, global_counts[1] * 9);
	printf ("Divergent clipping centers = %i out of %i\n"			, digergent_clipping_centers, global_counts[1]);
	printf ("Divergent clipping radii = %i out of %i\n"				, digergent_clipping_radii, global_counts[1]);
	printf ("Divergent major radii = %i out of %i\n"				, divergent_major_radii, global_counts[1]);
	printf ("Divergent rcoi = %i out of %i\n"						, divergent_rcoi, global_counts[1]);


	st_sesComplex.clear();

	#endif // CHECK_BUILDUP_DIFF


	#if !defined(MULTITHREADED_SES_BUILDING) || defined(CHECK_BUILDUP_DIFF)
	l.clear();
	#endif


	type[POINT_CELL]         = num_cells[POINT_CELL];
	type[REGULAR_EDGE_CELL]  = num_cells[REGULAR_EDGE_CELL];
	type[SINGULAR_EDGE_CELL] = num_cells[SINGULAR_EDGE_CELL];
	type[REGULAR_FACE_CELL]  = num_cells[REGULAR_FACE_CELL];
	type[SINGULAR_FACE_CELL] = num_cells[SINGULAR_FACE_CELL];


	auto chrono_end = chrono::high_resolution_clock::now();

	chrono::duration<double> build_time = chrono_end - chrono_start;
	cout << endl << INFO << "Surface build-up time.. ";
	printf ("%.4e [s]", build_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	if (!conf.parallelPocketLoop)
	{
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after build-up is " << current_mem_in_MB << " MB";
	}
	#endif


	if (savePovRay)
	{
		output_file.open(outputFile);

		output_file << "#include \"shapes.inc\" ";
		output_file << "\n#include \"colors.inc\" ";
		output_file << "\nglobal_settings {max_trace_level 3}";
		output_file << "\nbackground { color Black }";
		output_file << "\ncamera {";
		output_file << "\n\tlocation  <" << mid_x << ","  << max_y <<"," << mid_z <<">";
		output_file << "\fnlook_at  <" << mid_x << "," << mid_y << "," << mid_z << "> translate-<" << mid_x << "," << mid_y << "," << mid_z << ">" << " rotate<0,0,clock> "<< "translate<" << mid_x << "," << mid_y << "," << mid_z << ">}";
		output_file << "\nlight_source {<" << mid_x << ","  << max_y <<"," << mid_z <<">" << " color White "<< "translate-<" << mid_x << "," << mid_y << "," << mid_z << ">" << " rotate<0,0,clock> "<< "translate<" << mid_x << "," << mid_y << "," << mid_z << ">}";

		for (unsigned int it=0; it<sesComplex.size(); it++)
		{
			ConnollyCell *cc = sesComplex[it];

			if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
			{
				EdgeCell *ec = (EdgeCell*)cc;

				double u[3], v[3], w[3];

				v[0] = ec->Rot[0][0];
				u[0] = ec->Rot[0][1];
				w[0] = ec->Rot[0][2];
				v[1] = ec->Rot[1][0];
				u[1] = ec->Rot[1][1];
				w[1] = ec->Rot[1][2];
				v[2] = ec->Rot[2][0];
				u[2] = ec->Rot[2][1];
				w[2] = ec->Rot[2][2];

				if (ec->additional_planes.size() == 0)
				{
					// 2 is the number of cutting planes
					saveEdgePatch(output_file, ec, 2,
								  ec->major_radius, u, v, w, ec->isComplex);
				}
				else
				{
					saveEdgePatch(output_file, ec,
								  ec->additional_planes.size(), ec->major_radius, u, v, w, ec->isComplex);
				}
			}
		}

		for (unsigned int i=0; i<sesComplex.size(); i++)
		{
			ConnollyCell *cp = sesComplex[i];

			if (cp->patch_type == REGULAR_FACE_CELL || cp->patch_type == SINGULAR_FACE_CELL)
			{
				FacetCell *fc = (FacetCell*)cp;

				saveConcaveSpherePatch(output_file,fc, i);
			}
		}


		for (unsigned int i=0;i<sesComplex.size();i++)
		{
			ConnollyCell *cp = sesComplex[i];

			if (cp->patch_type == POINT_CELL)
			{
				PointCell *pc = (PointCell*)cp;

				saveAtomPatch(output_file, pc);
			}
		}

		output_file.close();
	}

	char fileName[1024];
	sprintf(fileName,"%sexposed.xyz",rootFile.c_str());
	
	FILE *fpExp = fopen(fileName, "w");

	fprintf(fpExp,"%zu\n", exposed.size());
	fprintf(fpExp,"Exposed_atoms\n");
	for (unsigned int i=0; i<exposed.size(); i++)
		fprintf(fpExp, "C %f %f %f\n",
				delphi->atoms[exposed[i]].pos[0], delphi->atoms[exposed[i]].pos[1], delphi->atoms[exposed[i]].pos[2]);
	fclose(fpExp);

	sprintf(fileName, "%sexposedIndices.txt", rootFile.c_str());
	fpExp = fopen(fileName, "w");
	for (unsigned int i=0; i<exposed.size(); i++)
		fprintf(fpExp,"%d\n", exposed[i]);
	fclose(fpExp);


	printSummary();

	return true;
}


void ConnollySurface::preProcessPanel()
{
	if (sesComplex.size() == 0)
	{
		cout << endl << WARN << "Cannot get surface with an empty SES complex";
		return;
	}

	int64_t igrid = delphi->nx;
	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption

	// auxiliary grid is consistent with delphi grid in order to speed-up tracing
	int64_t gridMul = 1;
	while (igrid > AUX_GRID_DIM_CONNOLLY_2D)
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

	// cout << endl << INFO << "Auxiliary grid is " << igrid;
	
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
	if (gridConnollyCellMap2D != NULL)
		#if !defined(MULTITHREADED_SES_BUILDING)
		deleteVector<int>(gridConnollyCellMap2D);
		#else
		deleteVector<pair<int,int>>(gridConnollyCellMap2D);
		#endif
	
	if (ind_2d != NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,last_cols_ind,ind_2d);
	#else */
	if (gridConnollyCellMap2D != NULL)
	{
		for (int64_t i=0; i < last_rows_ind*last_cols_ind; i++)
		{
			gridConnollyCellMap2D[i].clear();
		}
		delete[] gridConnollyCellMap2D;
		gridConnollyCellMap2D = NULL;
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
	#if !defined(MULTITHREADED_SES_BUILDING)
	gridConnollyCellMap2D = allocateVector<int>(last_rows_ind*last_cols_ind*MAX_CONNOLLY_CELLS_2D);
	#else
	gridConnollyCellMap2D = allocateVector<pair<int,int>>(last_rows_ind*last_cols_ind*MAX_CONNOLLY_CELLS_2D);
	#endif

	// look up, it has been deallocated
	ind_2d = allocateMatrix2D<unsigned int>(last_rows_ind,last_cols_ind);
	for (int i=0;i<last_rows_ind;i++)
		for (int j=0;j<last_cols_ind;j++)
			ind_2d[i][j] = 0;

	int max_t = 0;
	#else */
	#if !defined(MULTITHREADED_SES_BUILDING)
	gridConnollyCellMap2D = new vector<int> [ last_rows_ind*last_cols_ind ];
	#else
	gridConnollyCellMap2D = new vector<pair<int,int>> [ last_rows_ind*last_cols_ind ];
	#endif
	// #endif

	int max_num_cells_per_pixel = 0;

	// build a bounding box for cell and map it to the auxiliary grid
	#if !defined(MULTITHREADED_SES_BUILDING)
	for (unsigned int it=0; it<sesComplex.size(); it++)
	{
		{
			{
				// map connolly cell
				ConnollyCell *cc = sesComplex[it];

				if (cc->patch_type == SKIP_CELL)
					continue;
	#else
	for (int task_id=0; task_id<numTasks; task_id++)
	{
		for (int thd_id=0; thd_id<numThreadDataWrappersPerTask[task_id]; thd_id++)
		{
			int tdw_id = task_id * numThreadDataWrappers + thd_id;

			for (unsigned int it=0; it < thread_data_wrapper[ tdw_id ].sesComplex.size(); it++)
			{
				ConnollyCell *cc = thread_data_wrapper[ tdw_id ].sesComplex[it];

				if (cc->patch_type == SKIP_CELL)
					continue;
	#endif
				double *sphere_center;
				double radius;

				if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
				{
					FacetCell *fc = (FacetCell*)cc;
					sphere_center = fc->center;
					radius = probe_radius;
				}
				else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
				{
					EdgeCell *ec = (EdgeCell*)cc;
					sphere_center = ec->clipping_center;
					radius = ec->clipping_radius;
				}
				else if (cc->patch_type == POINT_CELL)
				{
					PointCell *pc = (PointCell*)cc;
					radius = delphi->atoms[pc->id].radius;
					sphere_center = delphi->atoms[pc->id].pos;
				}

				// compute the bounding box of the object
				double downx = sphere_center[0]-radius;
				double downy = sphere_center[1]-radius;
				double downz = sphere_center[2]-radius;

				double upx = sphere_center[0]+radius;
				double upy = sphere_center[1]+radius;
				double upz = sphere_center[2]+radius;
	          
				// Determine which are the grid cells that are
				// occupied by the bounding box of the object
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

							if (ind_2d[iz][iy] >= MAX_CONNOLLY_CELLS_2D)
							{
								cout << endl << ERR << "Number of connolly cells is superior to maximum allowed, please increase Max_ses_patches_per_auxiliary_grid_2d_cell";
								exit(-1);
							}
							GRID_CONNOLLY_CELL_MAP_2D(iy,iz,ind_2d[iz][iy],ny_2d,nz_2d) = it;
							ind_2d[iz][iy]++;
							#else */
							max_num_cells_per_pixel = max(max_num_cells_per_pixel, (int)gridConnollyCellMap2D[ iz*ny_2d + iy ].size());
							#if !defined(MULTITHREADED_SES_BUILDING)
							gridConnollyCellMap2D[ iz*ny_2d + iy ].push_back(it);
							#else
							gridConnollyCellMap2D[ iz*ny_2d + iy ].push_back( pair<int,int>(tdw_id, it) );
							#endif
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

							if (ind_2d[iy][ix] >= MAX_CONNOLLY_CELLS_2D)
							{
								cout << endl << ERR << "Number of connolly cells is superior to maximum allowed, please increase Max_ses_patches_per_auxiliary_grid_2d_cell";
								exit(-1);
							}
							GRID_CONNOLLY_CELL_MAP_2D(ix,iy,ind_2d[iy][ix],nx_2d,ny_2d) = it;
							ind_2d[iy][ix]++;
							#else */
							max_num_cells_per_pixel = max(max_num_cells_per_pixel, (int)gridConnollyCellMap2D[ iy*nx_2d + ix ].size());
							#if !defined(MULTITHREADED_SES_BUILDING)
							gridConnollyCellMap2D[ iy*nx_2d + ix ].push_back(it);
							#else
							gridConnollyCellMap2D[ iy*nx_2d + ix ].push_back( pair<int,int>(tdw_id, it) );
							#endif
							// #endif
						}
					}
				}
				else
				{
					// plane XZ
					for (int64_t iz=iz_start;iz <= iz_end;iz++)
					{
						for (int64_t ix=ix_start;ix <= ix_end;ix++)
						{
							/* #if !defined(OPTIMIZE_GRIDS)
							if (ind_2d[iz][ix] > max_t)
								max_t = ind_2d[iz][ix];
				
							if (ind_2d[iz][ix] >= MAX_CONNOLLY_CELLS_2D)
							{
								cout << endl << ERR << "Number of connolly cells is superior to maximum allowed, please increase Max_ses_patches_per_auxiliary_grid_cell";
								exit(-1);
							}
							GRID_CONNOLLY_CELL_MAP_2D(ix,iz,ind_2d[iz][ix],nx_2d,nz_2d) = it;
							ind_2d[iz][ix]++;
							#else */
							max_num_cells_per_pixel = max(max_num_cells_per_pixel, (int)gridConnollyCellMap2D[ iz*nx_2d + ix ].size());
							#if !defined(MULTITHREADED_SES_BUILDING)
							gridConnollyCellMap2D[ iz*nx_2d + ix ].push_back(it);
							#else
							gridConnollyCellMap2D[ iz*nx_2d + ix ].push_back( pair<int,int>(tdw_id, it) );
							#endif
							// #endif
						}
					}
				}
			}
		}
	}
	cout << endl << INFO << "Max number of SES cells per panel pixel: " << max_num_cells_per_pixel;

	/*
	if (gridLoad != NULL)
		deleteVector<int>(gridLoad);

	// printf("\n Load balancer...");
	// the load is splitted along the dimension defined by surface class
	// query the auxiliary on how many patches are encountered by that ray
	if (panel == 0)
	{
		totalLoad = 0;
		gridLoad = allocateVector<int>(delphi->nz);
		for (int k=0; k<delphi->nz; k++)
		{
			gridLoad[k] = 0;
			int t1 = (int)rintp((delphi->z[k]-zmin_2d)*scale_2d);
			// along Y
			for (int j=0;j<delphi->ny;j++)
			{
				int t2 = (int)rintp((delphi->y[j]-ymin_2d)*scale_2d);
				gridLoad[k] += ind_2d[t1][t2];
			}
			totalLoad += gridLoad[k];
		}
	}
	else if (panel == 1)
	{
		totalLoad = 0;
		gridLoad = allocateVector<int>(delphi->ny);
		for (int j=0; j<delphi->ny; j++)
		{
			gridLoad[j] = 0;
			int t1 = (int)rintp((delphi->y[j]-ymin_2d)*scale_2d);
			// along X
			for (int i=0;i<delphi->nx;i++)
			{
				int t2 = (int)rintp((delphi->x[i]-xmin_2d)*scale_2d);
				gridLoad[j] += ind_2d[t1][t2];
			}
			totalLoad += gridLoad[j];
		}
	}
	else 
	{
		totalLoad = 0;
		gridLoad = allocateVector<int>(delphi->nz);
		for (int k=0;k<delphi->nz;k++)
		{
			gridLoad[k] = 0;
			int t1 = (int)rintp((delphi->z[k]-zmin_2d)*scale_2d);
			// along X
			for (int i=0; i<delphi->nx; i++)
			{
				int t2 = (int)rintp((delphi->x[i]-xmin_2d)*scale_2d);
				gridLoad[k] += ind_2d[t1][t2];
			}
			totalLoad += gridLoad[k];
		}
	}
	*/
}


bool ConnollySurface::buildAuxiliaryGrid()
{
	if (sesComplex.size() == 0)
	{
		cout << endl << WARN << "Cannot get surface with an empty SES complex";
		return false;
	}
	
	// Perform pre-processing to speed-up intersections and projections.
	// Compute bounding box for each connoly cell in the surface, and map each
	// bounding box to the proper grid point.
	// This routine uses an auxiliary grid. Delphi grid is not directly used
	// because it may be too memory consuming. To speed up computations
	// a maximal number of connolly cells is allowed in each auxiliary grid cell
	// the macro is MAX_CONNOLLY_CELLS in ConnollySurface.h.
	// The macro AUX_GRID_DIM_CONNOLLY sets the maximally allowed grid size.
	int64_t igrid = delphi->nx;
	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption
	
	// auxiliary grid is consistent with delphi grid in order to speed-up tracing
	int64_t gridMul = 1;
	while (igrid > AUX_GRID_DIM_CONNOLLY)
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
	
	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);
	
	/* #if !defined(OPTIMIZE_GRIDS)
	// delete before the new dimensions are computed
	if (ind != NULL)
		deleteMatrix3D<int>(nx,ny,nz,ind);
	#endif */
	
	nx = igrid;
	ny = igrid;
	nz = igrid;
	
	/* #if !defined(OPTIMIZE_GRIDS)
	if (gridConnollyCellMap != NULL)
		deleteVector<int>(gridConnollyCellMap);
	
	cout << endl << INFO << "Allocating " << (nx*ny*nz*MAX_CONNOLLY_CELLS)*sizeof(int)/1024.0/1024.0 << " MB" << " for the auxiliary grid...";
	
	gridConnollyCellMap = allocateVector<int>(nx*ny*nz*MAX_CONNOLLY_CELLS);
	
	// look up, it has been deallocated!
	ind = allocateMatrix3D<int>(nx,ny,nz);
	
	for (int k=0; k<nz; k++)
		for (int j=0; j<ny; j++)
			for (int i=0; i<nx; i++)
				ind[k][j][i] = 0;
	#else */
	if (gridConnollyCellMap != NULL)
	{
		for (int64_t i=0; i < nx*ny*nz; i++)
		{
			gridConnollyCellMap[i].clear();
		}
		delete[] gridConnollyCellMap;
		gridConnollyCellMap = NULL;
	}
	#if !defined(MULTITHREADED_SES_BUILDING)
	gridConnollyCellMap = new vector<int> [ nx*ny*nz ];
	#else
	gridConnollyCellMap = new vector<pair<int,int>> [ nx*ny*nz ];
	#endif
	// #endif
	
	for (int i=0; i<nx; i++)
		x[i] = xmin + i*side;
	
	for (int i=0; i<ny; i++)
		y[i] = ymin + i*side;
	
	for (int i=0; i<nz; i++)
		z[i] = zmin + i*side;
	
	cout << "ok!";
	
	// build a bounding box for cell and map it to
	// the auxiliary grid
	int max_t = 0;
	
	cout << endl << INFO << "Mapping auxiliary grid...";
	
	#if !defined(MULTITHREADED_SES_BUILDING)
	for (unsigned int it=0; it<sesComplex.size(); it++)
	{
		{
			{
				// map connolly cell
				ConnollyCell *cc = sesComplex[it];

				if (cc->patch_type == SKIP_CELL)
					continue;
	#else
	for (int task_id=0; task_id<numTasks; task_id++)
	{
		for (int thd_id=0; thd_id<numThreadDataWrappersPerTask[task_id]; thd_id++)
		{
			int tdw_id = task_id * numThreadDataWrappers + thd_id;

			for (unsigned int it=0; it < thread_data_wrapper[ tdw_id ].sesComplex.size(); it++)
			{
				ConnollyCell *cc = thread_data_wrapper[ tdw_id ].sesComplex[it];

				if (cc->patch_type == SKIP_CELL)
					continue;
	#endif
				double *sphere_center;
				double radius;

				if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
				{
					FacetCell *fc = (FacetCell*)cc;
					sphere_center = fc->center;
					radius = probe_radius;
				}
				else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
				{
					EdgeCell *ec = (EdgeCell*)cc;
					sphere_center = ec->clipping_center;
					radius = ec->clipping_radius;
				}
				else if (cc->patch_type == POINT_CELL)
				{
					PointCell *pc = (PointCell*)cc;
					radius = delphi->atoms[pc->id].radius;
					sphere_center = delphi->atoms[pc->id].pos;
				}

				double downx = sphere_center[0]-radius;
				double downy = sphere_center[1]-radius;
				double downz = sphere_center[2]-radius;

				double upx = sphere_center[0]+radius;
				double upy = sphere_center[1]+radius;
				double upz = sphere_center[2]+radius;

				// Determine which are the grid cubes that
				// are occupied by the object's bounding box
				int64_t ix_start = (int64_t)rintp(fmax(0., downx-xmin)*scale);
				int64_t iy_start = (int64_t)rintp(fmax(0., downy-ymin)*scale);
				int64_t iz_start = (int64_t)rintp(fmax(0., downz-zmin)*scale);

				int64_t ix_end = (int64_t)rintp(fmax(0., upx-xmin)*scale);
				int64_t iy_end = (int64_t)rintp(fmax(0., upy-ymin)*scale);
				int64_t iz_end = (int64_t)rintp(fmax(0., upz-zmin)*scale);
		
				if (ix_start >= nx)
					ix_start = nx-1;
				if (iy_start >= ny)
					iy_start = ny-1;
				if (iz_start >= nz)
					iz_start = nz-1;

				if (ix_end >= nx)
					ix_end = nx-1;
				if (iy_end >= ny)
					iy_end = ny-1;
				if (iz_end >= nz)
					iz_end = nz-1;

				for (int64_t iz=iz_start; iz<=iz_end; iz++)
				{
					for (int64_t iy=iy_start; iy<=iy_end; iy++)
					{
						for (int64_t ix=ix_start; ix<=ix_end; ix++)
						{
							/* #if !defined(OPTIMIZE_GRIDS)
							if (ind[iz][iy][ix] > max_t)
								max_t = ind[iz][iy][ix];

							if (ind[iz][iy][ix] >= MAX_CONNOLLY_CELLS)
							{
								cout << endl << ERR << "Number of connolly cells is superior to maximum allowed, please increase Max_ses_patches_per_auxiliary_grid_2d_cell";
								exit(-1);
							}
							GRID_CONNOLLY_CELL_MAP(ix,iy,iz,ind[iz][iy][ix],nx,ny,nz) = it;
							ind[iz][iy][ix]++;
							#else */
							#if !defined(MULTITHREADED_SES_BUILDING)
							gridConnollyCellMap[ iz*ny*nx + iy*nx + ix ].push_back(it);
							#else
							gridConnollyCellMap[ iz*ny*nx + iy*nx + ix ].push_back( pair<int,int>(tdw_id, it) );
							#endif
							// #endif
						}
					}
				}
			}
		}
	}
	
	cout << "ok!";
	cout << endl << INFO << "Max Connolly cells per auxiliary cell -> " << max_t;
	
	return true;
}


#if !defined(SINGLE_PASS_RT)

#if !defined(MINIMIZE_MEMORY)

void ConnollySurface::getPatchPreIntersectionData (int64_t nxyz[3], int panels[2],
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
	for (int it = thread_id; it < sesComplex.size(); it += conf.numThreads)
	{
		// map connolly cell
		ConnollyCell *cc = sesComplex[it];

		if (cc->patch_type == SKIP_CELL)
			continue;

		double *sphere_center;
		double radius;

		if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
		{
			FacetCell *fc = (FacetCell*)cc;
			sphere_center = fc->center;
			radius = probe_radius;
		}
		else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
		{
			EdgeCell *ec = (EdgeCell*)cc;
			sphere_center = ec->clipping_center;
			radius = ec->clipping_radius;
		}
		else if (cc->patch_type == POINT_CELL)
		{
			PointCell *pc = (PointCell*)cc;
			sphere_center = delphi->atoms[pc->id].pos;
			radius = delphi->atoms[pc->id].radius;
		}
		double downx = sphere_center[0]-radius;
		double downy = sphere_center[1]-radius;
		double downz = sphere_center[2]-radius;

		double upx = sphere_center[0]+radius;
		double upy = sphere_center[1]+radius;
		double upz = sphere_center[2]+radius;

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
				int first_dim = 1, last__dim = 2;
				if (panel == 1)
					first_dim = 0, last__dim = 1;
				else if (panel == 2)
					first_dim = 0, last__dim = 2;

				#if defined(MULTITHREADING)
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);
				#endif

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
		}
	}
}

#else // MINIMIZE_MEMORY

void ConnollySurface::getPatchPreIntersectionData (int64_t nxyz[3], int panels[2],
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
	for (int it = thread_id; it < sesComplex.size(); it += conf.numThreads)
	{
		// map connolly cell
		ConnollyCell *cc = sesComplex[it];

		if (cc->patch_type == SKIP_CELL)
			continue;

		double *sphere_center;
		double *torus_center;
		EdgeCell *ec;
		double radius;

		bool doTorus = false;

		if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
		{
			FacetCell *fc = (FacetCell*)cc;
			sphere_center = fc->center;
			radius = probe_radius;
		}
		else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
		{
			ec = (EdgeCell*)cc;
			sphere_center = ec->clipping_center;
			radius = ec->clipping_radius;
			torus_center = ec->center;
			doTorus = true;
		}
		else if (cc->patch_type == POINT_CELL)
		{
			PointCell *pc = (PointCell*)cc;
			sphere_center = delphi->atoms[pc->id].pos;
			radius = delphi->atoms[pc->id].radius;
		}
		double downx = sphere_center[0]-radius;
		double downy = sphere_center[1]-radius;
		double downz = sphere_center[2]-radius;

		double upx = sphere_center[0]+radius;
		double upy = sphere_center[1]+radius;
		double upz = sphere_center[2]+radius;

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
				double t[4];
				double dir[3] = {0., 0., 0.};
				double ray_dir;

				int first_dim, last__dim;
				int varying_coord;

				if (panel == 0)
				{
					first_dim = 1, last__dim = 2;
					pa[0] = delphi->x[0] + delta;
					dir[0] = ray_dir = delphi->x[NX-1] - delphi->x[0];
					varying_coord = 0;
				}
				else if (panel == 1)
				{
					first_dim = 0, last__dim = 1;
					pa[2] = delphi->z[0] + delta;
					dir[2] = ray_dir = delphi->z[NZ-1] - delphi->z[0];
					varying_coord = 2;
				}
				else
				{
					first_dim = 0, last__dim = 2;
					pa[1] = delphi->y[0] + delta;
					dir[1] = ray_dir = delphi->y[NY-1] - delphi->y[0];
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
					#ifdef EQ_CULLING
					double dx = sphere_center[first_dim] - pa[first_dim];
					double dy = sphere_center[last__dim] - pa[last__dim];

					if (dx*dx + dy*dy > radius*radius)
						continue;
					#endif

					bool det = raySphere(pa,dir,sphere_center,radius,&t[0],&t[1]);

					if (!det)
						continue;

					int numInt = 0;
					// if it is a sphere we have finished
					if (!doTorus)
					{
						// there are 2 ray-sphere intersections
						numInt = 2;
					}
					else
					{
						rayTorus (analyticalTorusIntersectionAlgorithm, ec->invrot, torus_center, sphere_center,
								  probe_radius, ec->major_radius, radius, panel, pa, dir, t, numInt);
					}

					for (int i=0; i<numInt; i++)
					{
						double intPoint[3] = {pa[0], pa[1], pa[2]};
						intPoint[varying_coord] += t[i] * ray_dir;

						if (!isFeasible(cc,intPoint))
							continue;

						#if !defined(MULTITHREADING)
						// the n- and m-dependent integers are rarely computed here since it is very likely
						// that the ray misses the patch
						int64_t panel_pixel_id = (n*N_MAX+m)*numPanels + (panels[0]*int_phase + panel);
						++num_pixel_intersections[panel_pixel_id];
						#else
						++num_patch_intersections[patch_id];
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


void ConnollySurface::getPatchIntersectionData (int64_t nxyz[3], int panels[2],
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
	for (int it = thread_id; it < sesComplex.size(); it += conf.numThreads)
	{
		// map connolly cell
		ConnollyCell *cc = sesComplex[it];

		if (cc->patch_type == SKIP_CELL)
			continue;

		double *sphere_center;
		double *torus_center;
		EdgeCell *ec;
		double radius;

		bool doTorus = false;

		if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
		{
			FacetCell *fc = (FacetCell*)cc;
			sphere_center = fc->center;
			radius = probe_radius;
		}
		else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
		{
			ec = (EdgeCell*)cc;
			sphere_center = ec->clipping_center;
			radius = ec->clipping_radius;
			torus_center = ec->center;
			doTorus = true;
		}
		else if (cc->patch_type == POINT_CELL)
		{
			PointCell *pc = (PointCell*)cc;
			sphere_center = delphi->atoms[pc->id].pos;
			radius = delphi->atoms[pc->id].radius;
		}
		double downx = sphere_center[0]-radius;
		double downy = sphere_center[1]-radius;
		double downz = sphere_center[2]-radius;

		double upx = sphere_center[0]+radius;
		double upy = sphere_center[1]+radius;
		double upz = sphere_center[2]+radius;

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
				double t[4];
				double dir[3] = {0., 0., 0.};
				double ray_dir;

				int first_dim, last__dim;
				int varying_coord;

				if (panel == 0)
				{
					first_dim = 1, last__dim = 2;
					pa[0] = delphi->x[0] + delta;
					dir[0] = ray_dir = delphi->x[NX-1] - delphi->x[0];
					varying_coord = 0;
				}
				else if (panel == 1)
				{
					first_dim = 0, last__dim = 1;
					pa[2] = delphi->z[0] + delta;
					dir[2] = ray_dir = delphi->z[NZ-1] - delphi->z[0];
					varying_coord = 2;
				}
				else
				{
					first_dim = 0, last__dim = 2;
					pa[1] = delphi->y[0] + delta;
					dir[1] = ray_dir = delphi->y[NY-1] - delphi->y[0];
					varying_coord = 1;
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
					#ifdef EQ_CULLING
					double dx = sphere_center[first_dim] - pa[first_dim];
					double dy = sphere_center[last__dim] - pa[last__dim];

					if (dx*dx + dy*dy > radius*radius)
						continue;
					#endif

					bool det = raySphere(pa,dir,sphere_center,radius,&t[0],&t[1]);

					if (!det)
						continue;

					int numInt = 0;
					// if it is a sphere we have finished
					if (!doTorus)
					{
						// there are 2 ray-sphere intersections
						numInt = 2;
					}
					else
					{
						// if (!rayCell(cc, pa, dir, t))
						// 	continue;

						rayTorus (analyticalTorusIntersectionAlgorithm, ec->invrot, torus_center, sphere_center,
								  probe_radius, ec->major_radius, radius, panel, pa, dir, t, numInt);
					}

					for (int i=0; i<numInt; i++)
					{
						double intPoint[3] = {pa[0], pa[1], pa[2]};
						intPoint[varying_coord] += t[i] * ray_dir;

						if (!isFeasible(cc, intPoint))
							continue;

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

						int_p->first = t[i];

						VERTEX_TYPE *normal = nullptr;
						if (computeNormals)
						{
							double n[3];
							getNormal(intPoint,cc,n);

							normalsBuffers[thread_id].push_back(n[0]);
							normalsBuffers[thread_id].push_back(n[1]);
							normalsBuffers[thread_id].push_back(n[2]);

							normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];
						}
						int_p->second = normal;

						++*netIntersections;

						#else // SINGLE_PASS_RT

						VERTEX_TYPE *normal = nullptr;
						if (computeNormals)
						{
							double n[3];
							getNormal(intPoint,cc,n);

							normalsBuffers[thread_id].push_back(n[0]);
							normalsBuffers[thread_id].push_back(n[1]);
							normalsBuffers[thread_id].push_back(n[2]);

							normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];
						}
						temp_intersections_buffer[ patch_id ].push_back(pair<VERTEX_TYPE,VERTEX_TYPE*>(t[i], normal));

						intersection_pixel_id[ patch_id ].push_back(panel_pixel_id);

						++*netIntersections;

						#endif // SINGLE_PASS_RT
					}
				}
			}
		}
	}
}


/*
#if !defined(SINGLE_PASS_RT)

// This can be used if the normals' buffer common to the threads is filled after calculating and storing the intersections
void ConnollySurface::getPatchNormalsAtIntersections (int64_t nxyz[3], int panels[2], int thread_id)
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


	for (int it = thread_id; it < sesComplex.size(); it += conf.numThreads)
	{
		ConnollyCell *cc = sesComplex[it];

		if (cc->patch_type == SKIP_CELL)
			continue;

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
					double n[3];
					getNormal(intPoint,cc,n);

					normalsBuffers[thread_id].push_back(n[0]);
					normalsBuffers[thread_id].push_back(n[1]);
					normalsBuffers[thread_id].push_back(n[2]);

					int_p->second = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];
				}
			}
		}
	}
}

#else // SINGLE_PASS_RT

// This can be used if the normals' buffer common to the threads is filled after calculating and storing the intersections
void ConnollySurface::getPatchNormalsAtIntersections (int64_t nxyz[3], int panels[2], int thread_id)
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


	for (int it = thread_id; it < sesComplex.size(); it += conf.numThreads)
	{
		ConnollyCell *cc = sesComplex[it];

		if (cc->patch_type == SKIP_CELL)
			continue;

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
					int64_t n = pixel_id / N_MAX;
					int64_t m = pixel_id - n*N_MAX;

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
					double n[3];
					getNormal(intPoint,cc,n);

					normalsBuffers[thread_id].push_back(n[0]);
					normalsBuffers[thread_id].push_back(n[1]);
					normalsBuffers[thread_id].push_back(n[2]);

					int_p->second = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];
				}
			}
		}
	}
}

#endif
*/


bool ConnollySurface::isFeasible(ConnollyCell *cc, double point[3])
{
	if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
	{
		FacetCell *fc = (FacetCell*)cc;
		int start;
		
		if (cc->patch_type == SINGULAR_FACE_CELL)
			start = 0;
		else
			start = 1;
		
		for (int i=start; i<4; i++)
		{
			double tst = DOT(fc->planes[i],point) + fc->planes[i][3];
			if (tst > 0)
				return false;
		}
		
		for (unsigned int i=0; i<fc->self_intersection_planes.size(); i++)
		{
			double tst = DOT(fc->self_intersection_planes[i],point) + fc->self_intersection_planes[i][3];

			#if defined(OPTIMIZE_CELL_STRUCTURE)
			if (!fc->self_intersection_plane_labels[i])
				tst = -tst;
			#endif

			if (tst > 0)
				return false;
		}
		
		return true;
	}
	else if (cc->patch_type == SINGULAR_EDGE_CELL)
	{
		EdgeCell *ec = (EdgeCell*)cc;
		
		double r = ec->self_intersection_radius;

		// if (ec->isSelfIntersecting)
		if (r >= 0.)
		{
			// check for self intersection
			double dd;

			DIST2(dd,ec->center,point)

			if (dd < r*r)
				return false;
		}
		return true;
	}
	else if (cc->patch_type == REGULAR_EDGE_CELL)
	{
		EdgeCell *ec = (EdgeCell*)cc;
		
		double r = ec->self_intersection_radius;

		// if (ec->isSelfIntersecting)
		if (r >= 0.)
		{
			// check for self intersection
			double dd;

			DIST2(dd,ec->center,point)

			if (dd < r*r)
				return false;
		}
		
		// it is inside clipping sphere, check the clipping planes
		if (ec->additional_planes.size() != 0)
		{
			int plane_index = 0;
			for (unsigned int i = 0; i<ec->flags.size(); i++, plane_index += 2)
			{
				bool acute = ec->flags[i];
				
				#if !defined(OPTIMIZE_CELL_STRUCTURE)
				if (acute)
				{
					double tst = DOT(ec->additional_planes[plane_index],point) + ec->additional_planes[plane_index][3];
					if (tst > 0)
						continue;
					
					tst = DOT(ec->additional_planes[plane_index+1],point) + ec->additional_planes[plane_index+1][3];
					if (tst > 0)
						continue;
					
					return true;
				}
				else
				{
					double tst = DOT(ec->additional_planes[plane_index],point) + ec->additional_planes[plane_index][3];
					if (tst <= 0)
						return true;
					
					tst = DOT(ec->additional_planes[plane_index+1],point) + ec->additional_planes[plane_index+1][3];
					if (tst <= 0)
						return true;
				}
				#else // OPTIMIZE_CELL_STRUCTURE
				if (acute)
				{
					double tst = -(DOT(ec->additional_planes[plane_index],point) + ec->additional_planes[plane_index][3]);
					if (tst > 0)
						continue;

					tst = -(DOT(ec->additional_planes[plane_index+1],point) + ec->additional_planes[plane_index+1][3]);
					if (tst > 0)
						continue;

					return true;
				}
				else
				{
					double tst = -(DOT(ec->additional_planes[plane_index],point) + ec->additional_planes[plane_index][3]);
					if (tst <= 0)
						return true;

					tst = -(DOT(ec->additional_planes[plane_index+1],point) + ec->additional_planes[plane_index+1][3]);
					if (tst <= 0)
						return true;
				}
				#endif // OPTIMIZE_CELL_STRUCTURE
			}
			// no planes pair gives a positive result
			return false;
		}
		// fast path
		else
		{
			#if !defined(OPTIMIZE_CELL_STRUCTURE)
			if (ec->acute)
			{
				double tst = DOT(ec->cutting_planes[0],point) + ec->cutting_planes[0][3];
				if (tst > 0)
					return false;
				
				tst = DOT(ec->cutting_planes[1],point) + ec->cutting_planes[1][3];
				if (tst > 0)
					return false;
				
				return true;
			}
			else
			{
				double tst = DOT(ec->cutting_planes[0],point) + ec->cutting_planes[0][3];
				if (tst <= 0)
					return true;
				
				tst = DOT(ec->cutting_planes[1],point) + ec->cutting_planes[1][3];
				if (tst <= 0)
					return true;
				
				return false;
			}
			#else // OPTIMIZE_CELL_STRUCTURE
			if (ec->acute)
			{
				double tst = -(DOT(ec->cutting_planes[0],point) + ec->cutting_planes[0][3]);
				if (tst > 0)
					return false;

				tst = -(DOT(ec->cutting_planes[1],point) + ec->cutting_planes[1][3]);
				if (tst > 0)
					return false;

				return true;
			}
			else
			{
				double tst = -(DOT(ec->cutting_planes[0],point) + ec->cutting_planes[0][3]);
				if (tst <= 0)
					return true;

				tst = -(DOT(ec->cutting_planes[1],point) + ec->cutting_planes[1][3]);
				if (tst <= 0)
					return true;

				return false;
			}
			#endif // OPTIMIZE_CELL_STRUCTURE
		}
		
		return false;
	}
	else if (cc->patch_type == POINT_CELL)
	{
		PointCell *pc = (PointCell*)cc;
		
		#if !defined(OPTIMIZE_CELL_STRUCTURE)
		// filter on tori clipping spheres (shifted voronoi planes of exposed atoms)
		for (unsigned int i=0; i<pc->neighbours.size(); i++)
		{
			double *center = pc->neighbours[i]->clipping_center;
			double radius = pc->neighbours[i]->clipping_radius;
			double dd;

			DIST2(dd,center,point)
			dd -= radius*radius;

			if (dd < 0)
				return false;
		}
		// filter on buried clipping tori (shifted voronoi planes of buried atoms)
		for (unsigned int i=0; i<pc->buried_neighbours.size(); i++)
		{
			double *center = pc->buried_neighbours[i]->clipping_center;
			double radius = pc->buried_neighbours[i]->clipping_radius;
			double dd;

			DIST2(dd,center,point)
			dd -= radius*radius;

			if (dd < 0)
				return false;
		}
		#else // OPTIMIZE_CELL_STRUCTURE
		// filter on tori clipping spheres (shifted voronoi planes of exposed atoms)
		for (unsigned int i=0; i<pc->neighbour_data.size(); i += 4)
		{
			double *center = &pc->neighbour_data[i];
			double radius2 = pc->neighbour_data[i+3];
			double dd;

			DIST2(dd,center,point)
			dd -= radius2;

			if (dd < 0)
				return false;
		}
		for (unsigned int i=0; i<pc->buried_neighbour_data.size(); i += 4)
		{
			double *center = &pc->buried_neighbour_data[i];
			double radius2 = pc->buried_neighbour_data[i+3];
			double dd;

			DIST2(dd,center,point)
			dd -= radius2;

			if (dd < 0)
				return false;
		}
		#endif
		return true;
	}
	cout << endl << ERR << "Cannot get an answer in feasibility test";
	return false;
}


/*
// Incomplete and incorrect experimental function
bool ConnollySurface::rayCell(ConnollyCell *cc, double orig[3], double dir[3], double t[2])
{
	if (cc->patch_type != REGULAR_EDGE_CELL)
	{
		return true;
	}

	EdgeCell *ec = (EdgeCell*)cc;

	// t[*] are the two segment parameters given by the intersection between the ray point+t*dir and the sphere of the cell domain

	double entry_point[3], exit_point[3];
	double entry_value, exit_value;
	double p_new[3], t_new;

	int num_intersections = 0;

	entry_point[0] = orig[0];
	entry_point[1] = orig[1];
	entry_point[2] = orig[2];

	exit_point[0] = orig[0] + dir[0];
	exit_point[1] = orig[1] + dir[1];
	exit_point[2] = orig[2] + dir[2];

	// it is inside clipping sphere, check the clipping planes
	if (ec->additional_planes.size() != 0)
	{
		for (unsigned int i = 0; i<ec->flags.size(); i++)
		{
			int plane_index = i << 1;

			bool acute = ec->flags[i];

			#if !defined(OPTIMIZE_CELL_STRUCTURE)
			if (acute)
			{
				entry_value = DOT(ec->additional_planes[plane_index],entry_point) + ec->additional_planes[plane_index][3];

				exit_value = DOT(ec->additional_planes[plane_index],exit_point) + ec->additional_planes[plane_index][3];

				if (entry_value * exit_value <= 0.)
				{
					// determine the intersection parameter t given by the plane vs segmnt
					t_new = -entry_value / (DOT(ec->additional_planes[plane_index], dir));

					p_new[0] = entry_point[0] + t_new*dir[0];
					p_new[1] = entry_point[1] + t_new*dir[1];
					p_new[2] = entry_point[2] + t_new*dir[2];

					if (isFeasible(cc, p_new))
					{
						++num_intersections;
					}
				}
				entry_value = DOT(ec->additional_planes[plane_index+1],entry_point) + ec->additional_planes[plane_index+1][3];

				exit_value = DOT(ec->additional_planes[plane_index+1],exit_point) + ec->additional_planes[plane_index+1][3];

				if (entry_value * exit_value <= 0.)
				{
					// determine the intersection parameter t given by the plane vs segmnt
					t_new = -entry_value / (DOT(ec->additional_planes[plane_index+1], dir));

					p_new[0] = entry_point[0] + t_new*dir[0];
					p_new[1] = entry_point[1] + t_new*dir[1];
					p_new[2] = entry_point[2] + t_new*dir[2];

					if (isFeasible(cc, p_new))
					{
						++num_intersections;
					}
				}
			}
			else
			{
				// to be completed
				return true;
			}
			#else // OPTIMIZE_CELL_STRUCTURE
			if (acute)
			{
				entry_value = -(DOT(ec->additional_planes[plane_index],entry_point) + ec->additional_planes[plane_index][3]);

				exit_value = -(DOT(ec->additional_planes[plane_index],exit_point) + ec->additional_planes[plane_index][3]);

				if (entry_value * exit_value <= 0.)
				{
					// determine the intersection parameter t given by the plane vs segmnt
					t_new = entry_value / (DOT(ec->additional_planes[plane_index], dir));

					p_new[0] = entry_point[0] + t_new*dir[0];
					p_new[1] = entry_point[1] + t_new*dir[1];
					p_new[2] = entry_point[2] + t_new*dir[2];

					if (fabs(DOT(ec->additional_planes[plane_index],p_new) + ec->additional_planes[plane_index][3]) > 1.e-6) printf ("?");

					if (isFeasible(cc, p_new))
					{
						++num_intersections;
					}
				}
				entry_value = -(DOT(ec->additional_planes[plane_index+1],entry_point) + ec->additional_planes[plane_index+1][3]);

				exit_value = -(DOT(ec->additional_planes[plane_index+1],exit_point) + ec->additional_planes[plane_index+1][3]);

				if (entry_value * exit_value <= 0.)
				{
					// determine the intersection parameter t given by the plane vs segmnt
					t_new = entry_value / (DOT(ec->additional_planes[plane_index+1], dir));

					p_new[0] = entry_point[0] + t_new*dir[0];
					p_new[1] = entry_point[1] + t_new*dir[1];
					p_new[2] = entry_point[2] + t_new*dir[2];

					if (isFeasible(cc, p_new))
					{
						++num_intersections;
					}
				}
			}
			else
			{
				// to be completed
				return true;
			}
			#endif // OPTIMIZE_CELL_STRUCTURE
		}
	}
	// fast path
	else
	{
		#if !defined(OPTIMIZE_CELL_STRUCTURE)
		if (ec->acute)
		{
			entry_value = DOT(ec->cutting_planes[0],entry_point) + ec->cutting_planes[0][3];

			exit_value = DOT(ec->cutting_planes[0],exit_point) + ec->cutting_planes[0][3];

			if (entry_value * exit_value <= 0.)
			{
				// determine the intersection parameter t given by the plane vs segmnt
				t_new = -entry_value / DOT(ec->cutting_planes[0], dir);

				p_new[0] = entry_point[0] + t_new*dir[0];
				p_new[1] = entry_point[1] + t_new*dir[1];
				p_new[2] = entry_point[2] + t_new*dir[2];

				if (isFeasible(cc, p_new))
				{
					++num_intersections;
				}
			}
			entry_value = DOT(ec->cutting_planes[1],entry_point) + ec->cutting_planes[1][3];

			exit_value = DOT(ec->cutting_planes[1],exit_point) + ec->cutting_planes[1][3];

			if (entry_value * exit_value <= 0.)
			{
				// determine the intersection parameter t given by the plane vs segmnt
				t_new = -entry_value / DOT(ec->cutting_planes[1], dir);

				p_new[0] = entry_point[0] + t_new*dir[0];
				p_new[1] = entry_point[1] + t_new*dir[1];
				p_new[2] = entry_point[2] + t_new*dir[2];

				if (isFeasible(cc, p_new))
				{
					++num_intersections;
				}
			}
		}
		else
		{
			// to be completed
			return true;
		}
		#else // OPTIMIZE_CELL_STRUCTURE
		if (ec->acute)
		{
			entry_value = DOT(ec->cutting_planes[0],entry_point) + ec->cutting_planes[0][3];

			exit_value = DOT(ec->cutting_planes[0],exit_point) + ec->cutting_planes[0][3];

			if (entry_value * exit_value <= 0.)
			{
				// determine the intersection parameter t given by the plane vs segmnt
				t_new = -entry_value / DOT(ec->cutting_planes[0], dir);

				p_new[0] = entry_point[0] + t_new*dir[0];
				p_new[1] = entry_point[1] + t_new*dir[1];
				p_new[2] = entry_point[2] + t_new*dir[2];

				if (isFeasible(cc, p_new))
				{
					++num_intersections;
				}
			}
			entry_value = DOT(ec->cutting_planes[1],entry_point) + ec->cutting_planes[1][3];

			exit_value = DOT(ec->cutting_planes[1],exit_point) + ec->cutting_planes[1][3];

			if (entry_value * exit_value <= 0.)
			{
				// determine the intersection parameter t given by the plane vs segmnt
				t_new = -entry_value / DOT(ec->cutting_planes[1], dir);

				p_new[0] = entry_point[0] + t_new*dir[0];
				p_new[1] = entry_point[1] + t_new*dir[1];
				p_new[2] = entry_point[2] + t_new*dir[2];

				if (isFeasible(cc, p_new))
				{
					++num_intersections;
				}
			}
		}
		else
		{
			// to be completed
			return true;
		}
		#endif // OPTIMIZE_CELL_STRUCTURE
	}

	if (num_intersections > 0)
	{
		return true;
	}
	return false;
}
*/


bool ConnollySurface::rayConnollyCellIntersection(double *orig, double *dir, ConnollyCell *cc, double t[4], int &numInt)
{
	double *sphere_center;
	double *torus_center;
	EdgeCell *ec;
	double radius;

	bool doTorus = false;

	if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
	{
		FacetCell *fc = (FacetCell*)cc;
		sphere_center = fc->center;
		radius        = probe_radius;
	}
	else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
	{
		ec = (EdgeCell*)cc;
		sphere_center = ec->clipping_center;
		radius        = ec->clipping_radius;
		torus_center  = ec->center;
		doTorus = true;
	}
	else if (cc->patch_type == POINT_CELL)
	{
		PointCell *pc = (PointCell*)cc;
		sphere_center = delphi->atoms[pc->id].pos;
		radius        = delphi->atoms[pc->id].radius;
	}
	else
		cout << endl << ERR << "Cannot get patch type in intersection test";

	#ifdef EQ_CULLING
	int panel_dims[3][2] = {{1,2}, {0,1}, {0,2}};

	int n = panel_dims[panel][0];
	int m = panel_dims[panel][1];

	// culling if the ray does not pierce the bounding box of the patch
	if (orig[m] < sphere_center[m] - radius) return false;
	if (orig[m] > sphere_center[m] + radius) return false;
	if (orig[n] < sphere_center[n] - radius) return false;
	if (orig[n] > sphere_center[n] + radius) return false;
	#endif

	bool det = raySphere(orig,dir,sphere_center,radius,&t[0],&t[1]);

	if (!det)
		return false;

	numInt = 2;

	// if it is a sphere we have finished
	if (!doTorus)
	{
		return true;
	}

	return rayTorus (analyticalTorusIntersectionAlgorithm, ec->invrot, torus_center, sphere_center,
					 probe_radius, ec->major_radius, radius, panel, orig, dir, t, numInt);
}


#if defined(REPORT_FAILED_RAYS)
bool ConnollySurface::printRayConnollyCellIntersection(double *orig, double *dir, ConnollyCell *cc, double t[4], int &numInt)
{
	double *sphere_center;
	double *torus_center;
	EdgeCell *ec;
	double radius;

	bool doTorus = false;

	if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
	{
		FacetCell *fc = (FacetCell*)cc;
		sphere_center = fc->center;
		radius        = probe_radius;
	}
	else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
	{
		ec = (EdgeCell*)cc;
		sphere_center = ec->clipping_center;
		radius        = ec->clipping_radius;
		torus_center  = ec->center;
		doTorus = true;
	}
	else if (cc->patch_type == POINT_CELL)
	{
		PointCell *pc = (PointCell*)cc;
		sphere_center = delphi->atoms[pc->id].pos;
		radius        = delphi->atoms[pc->id].radius;
	}
	else
		cout << endl << ERR << "Cannot get patch type in intersection test";

	#ifdef EQ_CULLING
	int panel_dims[3][2] = {{1,2}, {0,1}, {0,2}};

	int n = panel_dims[panel][0];
	int m = panel_dims[panel][1];

	// culling if the ray does not pierce the bounding box of the patch
	if (orig[m] < sphere_center[m] - radius) return false;
	if (orig[m] > sphere_center[m] + radius) return false;
	if (orig[n] < sphere_center[n] - radius) return false;
	if (orig[n] > sphere_center[n] + radius) return false;
	#endif

	bool det = raySphere(orig,dir,sphere_center,radius,&t[0],&t[1]);

	if (!det)
		return false;

	numInt = 2;

	// if it is a sphere we have finished
	if (!doTorus)
	{
		return true;
	}

	double roots[4],pp[3],temp[3];
	double orig1[3];
	double ddir[3];
	double dir_norm;

	// start the ray from the clipping sphere intersection
	// ADD_MUL(orig1,orig,dir,t[0])
	// ADD_MUL(orig2,orig,dir,t[1])
	// SUB(tdir,orig2,orig1)

	// some optimised code ...
	int varying_coord = (panel == 0) ? 0 : ((panel == 1) ? 2 : 1);
	double ray_dir = dir[varying_coord];

	ASSIGN(orig1,orig);
	orig1[varying_coord] += t[0] * ray_dir;
	dir_norm = (t[1] - t[0]) * ray_dir;

	// build torus roots equation
	SUB(temp,orig1,torus_center)

	for (int i=0;i<3;i++)
	{
		// roto-translate origin point
		pp[i] = ec->invrot[i][0]*temp[0] + ec->invrot[i][1]*temp[1] + ec->invrot[i][2]*temp[2];
		// rotate ray direction
		ddir[i] = ec->invrot[i][varying_coord]*dir_norm;
	}
	double r2,R2;
	r2 = probe_radius*probe_radius;
	R2 = ec->major_radius*ec->major_radius;

	double pd,p2,coeff;

	// pd = DOT(pp,ddir);
	pd = pp[0]*ec->invrot[0][varying_coord] + pp[1]*ec->invrot[1][varying_coord] + pp[2]*ec->invrot[2][varying_coord];
	p2 = DOT(pp,pp);
	coeff = p2-r2-R2;

	double B = 4*pd;
	double C = 2*coeff + 4*pd*pd + 4*R2*ec->invrot[2][varying_coord]*ec->invrot[2][varying_coord];
	double D = 4*pd*coeff + 8*R2*pp[2]*ec->invrot[2][varying_coord];
	double E = coeff*coeff - 4*R2*(r2-pp[2]*pp[2]);

	int numroots;

	if (analyticalTorusIntersectionAlgorithm)
	{
		quarticEqSolutions(roots, B, C, D, E, &numroots);

		printf ("%.1f %.3f %.3f ", probe_radius, ec->major_radius, radius);
		printf ("%.14e %.14e %.14e ", pp[0], pp[1], pp[2]);
		printf ("%.14e %.14e %.14e ", ddir[0], ddir[1], ddir[2]);
		printf ("%.14e %.14e %.14e %.14e ", B, C, D, E);
		printf ("%.14e %.14e %.14e ", t[0], t[1], ray_dir);
		printf ("%i ", numroots);
		for (int i=0; i<numroots; i++)
			printf ("%.14e ", roots[i]);

		printf ("\n");
	}
	else
	{
		double polyy[5];

		numroots = 4;

		polyy[0] = E;
		polyy[1] = D;
		polyy[2] = C;
		polyy[3] = B;
		polyy[4] = 1.0;

		int order = 4;

		getRealRootsSturm(polyy,order,roots,numroots);
	}

	if (numroots == 0)
		return false;

	double t0 = t[0];
	numInt = 0;

	for (int i=0; i<numroots; i++)
	{
		if (roots[i] < 0)
			continue;

		double point[3], dd;

		// roots[i] /= dir_norm;
		// ADD_MUL(point,orig1,tdir,roots[i])
		ASSIGN(point,orig1)
		point[varying_coord] += roots[i]; // roots[i]*dir_norm
		DIST2(dd,sphere_center,point)
		// acceptable
		if (dd < radius*radius)
		{
			// t[ numInt++ ] = (point[v_coord]-orig[v_coord]) / ray_dir;
			//               = (orig1[v_coord] (=pa[v_coord]+t[0]*ray_dir) + roots[i] - pa[v_coord]) / ray_dir =
			//               = t0 + roots[i]/ray_dir
			t[ numInt++ ] = t0 + roots[i]/ray_dir;
		}
	}

	if (numInt == 0)
		return false;

	return true;
}
#endif


void ConnollySurface::getRayIntersection(double pa[3], double pb[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id)
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
	int numCells = ind_2d[i2][i1];
	#else */
	int numCells = gridConnollyCellMap2D[ i2*n_2d_first + i1 ].size();
	// #endif
	
	if (numCells == 0) return;
	
	double dir[3];
	
	dir[0] = 0.;
	dir[1] = 0.;
	dir[2] = 0.;
	dir[varying_coord] = pb[varying_coord] - pa[varying_coord];
	
	for (unsigned int iter = 0; iter<numCells; iter++)
	{
		double t[4];
		double intPoint[3];
		int ni;

		/* #if !defined(OPTIMIZE_GRIDS)
		int it = GRID_CONNOLLY_CELL_MAP_2D(i1,i2,iter,n_2d_first,n_2d_last);
		#else */
		#if !defined(MULTITHREADED_SES_BUILDING)
		int it = gridConnollyCellMap2D[i2*n_2d_first+i1][iter];

		ConnollyCell *cc = sesComplex[it];
		#else
		pair<int,int> it = gridConnollyCellMap2D[i2*n_2d_first+i1][iter];

		ThreadDataWrapper *tdw = &thread_data_wrapper[ it.first ];

		ConnollyCell *cc = tdw->sesComplex[ it.second ];
		#endif
		// #endif

		bool ff = rayConnollyCellIntersection(pa,dir,cc,t,ni);

		// no intersection
		if (!ff)
			continue;

		for (int i=0; i<ni; i++)
		{
			intPoint[0] = pa[0];
			intPoint[1] = pa[1];
			intPoint[2] = pa[2];
			intPoint[varying_coord] += t[i] * dir[varying_coord];
			
			// feasibility test: check if the intersection is inside the cell
			if (!isFeasible(cc, intPoint))
				continue;

			// compute normal
			if (computeNormals)
			{
				double n[3];
				getNormal(intPoint,cc,n);

				normalsBuffers[thread_id].push_back(n[0]);
				normalsBuffers[thread_id].push_back(n[1]);
				normalsBuffers[thread_id].push_back(n[2]);

				VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

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
void ConnollySurface::printRayIntersection(double pa[3], double pb[3])
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
	int numCells = ind_2d[i2][i1];
	#else */
	int numCells = gridConnollyCellMap2D[ i2*n_2d_first + i1 ].size();
	// #endif

	if (numCells == 0) return;

	double dir[3];

	dir[0] = 0.;
	dir[1] = 0.;
	dir[2] = 0.;
	dir[varying_coord] = pb[varying_coord] - pa[varying_coord];

	printf ("Ray start\n");

	for (unsigned int iter = 0; iter<numCells; iter++)
	{
		double t[4];
		double intPoint[3];
		int ni;

		/* #if !defined(OPTIMIZE_GRIDS)
		int it = GRID_CONNOLLY_CELL_MAP_2D(i1,i2,iter,n_2d_first,n_2d_last);
		#else */
		#if !defined(MULTITHREADED_SES_BUILDING)
		int it = gridConnollyCellMap2D[i2*n_2d_first+i1][iter];

		bool ff = printRayConnollyCellIntersection(pa,dir,sesComplex[it],t,ni);
		#else
		pair<int,int> it = gridConnollyCellMap2D[i2*n_2d_first+i1][iter];

		ThreadDataWrapper *tdw = &thread_data_wrapper[ it.first ];

		ConnollyCell *cc = tdw->sesComplex[ it.second ];

		bool ff = printRayConnollyCellIntersection(pa,dir,cc,t,ni);
		#endif
		// #endif
	}
	printf ("Ray end\n");
}
#endif // REPORT_FAILED_RAYS


void ConnollySurface::projectToCircle(double *point, double radius, double *center, double *plane, double *proj, double &dist)
{
	double pdist;

	point2plane(point, plane, &pdist, proj);
	SUB(proj,proj,center)
	radius /= sqrt((DOT(proj,proj)));
	ADD_MUL(proj,center,proj,radius)
	DIST(dist,point,proj)
}


void ConnollySurface::projectToTorus(double *y, EdgeCell *ec, double *proj, double *norm, double &dist)
{
	double temp[3], pp[3], torus_plane[4] = {0,0,+1,0}, proj_major[3], minor_plane[4], torus_center[3]={0,0,0};

	// roto-translate the point into the torus reference
	SUB(temp,y,ec->center)
	for (int i=0; i<3; i++)
	{
		pp[i] = ec->invrot[i][0]*temp[0] + ec->invrot[i][1]*temp[1] + ec->invrot[i][2]*temp[2];
	}
	// first project to the big circle
	projectToCircle(pp,ec->major_radius,torus_center,torus_plane,proj_major,dist);

	// need 2 more points to get the plane of the previously identified circle. 
	// Two good points are along the axis together with center of the circle.
	double p2[3]={0,0,+1}, p3[3]={0,0,-1};
	plane3points(proj_major,p2,p3,minor_plane);
	
	// second project to the small circle
	projectToCircle(pp,probe_radius,proj_major,minor_plane,proj,dist);

	// deroto-translate
	for (int i=0; i<3; i++)
	{
		pp[i] = ec->Rot[i][0]*proj[0] + ec->Rot[i][1]*proj[1] + ec->Rot[i][2]*proj[2];
	}
	ADD(proj,pp,ec->center)
	DIST(dist,proj,y)
}


void ConnollySurface::getNormal(double *y, ConnollyCell *cc, double *normal)
{
	double *sphere_center,sphere_radius;

	if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
	{
		FacetCell *fc = (FacetCell*)cc;
		sphere_center = fc->center;
		sphere_radius = probe_radius;
		SUB(normal,sphere_center,y)
		normal[0] /= sphere_radius;
		normal[1] /= sphere_radius;
		normal[2] /= sphere_radius;
	}
	else if (cc->patch_type == SINGULAR_EDGE_CELL || cc->patch_type == REGULAR_EDGE_CELL)
	{
		EdgeCell *ec = (EdgeCell*)cc;
		getNormalToTorus(y,ec,normal);
	}
	else if (cc->patch_type == POINT_CELL)
	{
		PointCell *pc = (PointCell*)cc;
		sphere_center = delphi->atoms[pc->id].pos;
		sphere_radius = delphi->atoms[pc->id].radius;
		SUB(normal,y,sphere_center)
		normal[0] /= sphere_radius;
		normal[1] /= sphere_radius;
		normal[2] /= sphere_radius;
	}
}


/*
void ConnollySurface::getNormalToTorus(double *y, EdgeCell *ec, double *normal)
{
	double temp[3], pp[3], torus_plane[4]={0,0,+1,0}, torus_center[3]={0,0,0}, dist, C[3], aa;
			
	// roto-translate the point into the torus reference 	
	SUB(temp,y,ec->center)
	for (int i=0; i<3; i++)
	{
		pp[i] = ec->invrot[i][0]*temp[0] + ec->invrot[i][1]*temp[1] + ec->invrot[i][2]*temp[2];
	}

	projectToCircle(pp,ec->major_radius,torus_center,torus_plane,C,dist);
	SUB(C,C,pp)
	NORMALIZE(C,aa)

	// de-rotate the normal vector
	for (int i=0; i<3; i++)
	{
		normal[i] = ec->Rot[i][0]*C[0] + ec->Rot[i][1]*C[1] + ec->Rot[i][2]*C[2];
	}
}
*/


// Optimised version of the above
void ConnollySurface::getNormalToTorus(double *y, EdgeCell *ec, double *normal)
{
	double temp[3], C[3];

	// roto-translate the point into the torus reference
	SUB(temp,y,ec->center)

	for (int i=0; i<3; i++)
	{
		C[i] = ec->invrot[i][0]*temp[0] + ec->invrot[i][1]*temp[1] + ec->invrot[i][2]*temp[2];
	}
	double d = C[0]*C[0] + C[1]*C[1];
	double radius_minus_1 = ec->major_radius/sqrt(d) - 1.;

	d = sqrt(d*radius_minus_1*radius_minus_1 + C[2]*C[2]);

	C[0] *= radius_minus_1 / d;
	C[1] *= radius_minus_1 / d;
	C[2] /= -d;

	// de-rotate the normal vector
	for (int i=0; i<3; i++)
	{
		normal[i] = ec->Rot[i][0]*C[0] + ec->Rot[i][1]*C[1] + ec->Rot[i][2]*C[2];
	}
}


bool ConnollySurface::getProjection(double p[3], double *proj1, double *proj2,
									double *proj3, double *normal1, double *normal2, double *normal3)
{
	// get the cells that are associated to this grid point
	// by querying the auxiliary grid
	double dist;
	#if !defined(MULTITHREADED_SES_BUILDING)
	set<int>cells;
	#else
	set<pair<int,int>>cells;
	#endif
	
	// move from delphi grid to auxiliary grid
	int64_t irefx = (int64_t)rintp((p[0]-xmin)*scale);
	int64_t irefy = (int64_t)rintp((p[1]-ymin)*scale);
	int64_t irefz = (int64_t)rintp((p[2]-zmin)*scale);
	
	/* #if !defined(OPTIMIZE_GRIDS)
	for (int i=0; i < ind[irefz][irefy][irefx]; i++)
		cells.insert((GRID_CONNOLLY_CELL_MAP(irefx,irefy,irefz,i,nx,ny,nz)));
	#else */
	for (int i=0; i < gridConnollyCellMap[ irefz*ny*nx + irefy*nx + irefx ].size(); i++)
	{
		cells.insert( gridConnollyCellMap[ irefz*ny*nx + irefy*nx + irefx ][i] );
	}
	// #endif
	
	// keep the nearest patch
	double locProj[3], locNorm[3];
	double minDist = INFINITY;
	
	dist = 0;
	locProj[0] = 0;
	locProj[1] = 0;
	locProj[2] = 0;
	
	locNorm[0] = 0;
	locNorm[1] = 0;
	locNorm[2] = 0;
	
	#if !defined(MULTITHREADED_SES_BUILDING)
	for (set<int>::iterator it = cells.begin(); it != cells.end(); it++)
	#else
	for (set<pair<int,int>>::iterator it = cells.begin(); it != cells.end(); it++)
	#endif
	{
		#if !defined(MULTITHREADED_SES_BUILDING)
		ConnollyCell *cc = sesComplex[*it];
		#else
		ConnollyCell *cc = thread_data_wrapper[ it->first ].sesComplex[ it->second ];

		if (cc->patch_type == SKIP_CELL)
			continue;
		#endif
		
		if (cc->patch_type == POINT_CELL)
		{
			PointCell *pc = (PointCell*)cc;
			projectToSphere(p,delphi->atoms[pc->id].pos,delphi->atoms[pc->id].radius,locProj,dist);
		}
		else if (cc->patch_type == REGULAR_FACE_CELL || cc->patch_type == SINGULAR_FACE_CELL)
		{
			FacetCell *fc = (FacetCell*)cc;
			projectToSphere(p,fc->center,probe_radius,locProj,dist);
		}
		else
		{
			EdgeCell *ec = (EdgeCell*)cc;
			projectToTorus(p,ec,locProj,locNorm,dist);
		}
		
		if (!isFeasible(cc,locProj))
			continue;
		
		if (dist < minDist)
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
	
	if (cells.size() == 0)
	{
		cout << endl << WARN << "Empty cell in getProjection!";
		return false;
	}
	
	// Connolly surface is not that nice. Projections can't be analytical everywhere.
	// We are trying to project a point into a singularity? Old DelPhi style.
	if (minDist == INFINITY)
	{
		bool fixed = false;
		
		double **sampledPoints;
		int coiNum = 500;
		sampledPoints = allocateMatrix2D<double>(coiNum,3);
		
		#if !defined(MULTITHREADED_SES_BUILDING)
		for (set<int>::iterator it = cells.begin(); it != cells.end(); it++)
		#else
		for (set<pair<int,int>>::iterator it = cells.begin(); it != cells.end(); it++)
		#endif
		{
			#if !defined(MULTITHREADED_SES_BUILDING)
			ConnollyCell *cc = sesComplex[*it];
			#else
			ConnollyCell *cc = thread_data_wrapper[ it->first ].sesComplex[ it->second ];
			#endif

			if (cc->patch_type == SINGULAR_EDGE_CELL ||  cc->patch_type == REGULAR_EDGE_CELL)
			{
				EdgeCell *ec = (EdgeCell*)cc;
				
				// try to manage singularity. Explicitly project to probes and keep the nearest one
				// if (ec->isSelfIntersecting)
				if (ec->self_intersection_radius >= 0.)
				{
					// identify the nearest the common probes and check the feasibility of the projection
					// for at least one of them.
					
					//  set of incident probes to a given edge
					FacetCell *fcv[MAX_INCIDENT_PROBES];
					int index = 0;
					
					// manage SINGULAR EDGE?
					if (ec->patch_type == SINGULAR_EDGE_CELL)
					{
						{
							#ifdef ENABLE_BOOST_THREADS
							boost::mutex::scoped_lock scopedLock(mutex);
							#endif
							(*errorStream) << endl << WARN << "Singular edge projection";
						}
					}
					
					// get all the sourrounding probe stations
					if (ec->patch_type == REGULAR_EDGE_CELL)
					{
						#if !defined(OPTIMIZE_CELL_STRUCTURE)

						#if !defined(MULTITHREADED_SES_BUILDING)
						#if !defined(NEW_ATOM_PATCHES)
						vector<FacetCell*> &f1 = atomPatches[ec->id[0]]->incidentProbes;
						#else
						PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[ec->id[0]] ];
						vector<FacetCell*> &f1 = pc1->incidentProbes;
						#endif
						#else // MULTITHREADED_SES_BUILDING
						ThreadDataWrapper *tdw = &thread_data_wrapper[ it->first ];
						#if !defined(NEW_ATOM_PATCHES)
						vector<FacetCell*> &f1 = tdw->atomPatches[ec->id[0]]->incidentProbes;
						#else
						int cell_id = tdw->atomPatches[ec->id[0]];
						PointCell *pc1 = (PointCell *)tdw->sesComplex[ cell_id ];
						vector<FacetCell*> &f1 = pc1->incidentProbes;
						#endif
						#endif // MULTITHREADED_SES_BUILDING
						
						int ref1 = ec->id[0];
						int ref2 = ec->id[1];
						int np = (int)f1.size();
						
						// get all the incident probes. Most of the cases will be 2.
						for (unsigned int indf=0; indf<np; indf++)
						{
							int i1 = (f1)[indf]->id[0];
							int i2 = (f1)[indf]->id[1];
							int i3 = (f1)[indf]->id[2];
							
							if ((i1 == ref1 && i2 == ref2) ||
								(i2 == ref1 && i1 == ref2) ||
								(i1 == ref1 && i3 == ref2) ||
								(i3 == ref1 && i1 == ref2) ||
								(i2 == ref1 && i3 == ref2) ||
								(i3 == ref1 && i2 == ref2))
							{
								fcv[index++] = (f1)[indf];
							}
						}
						#else // OPTIMIZE_CELL_STRUCTURE

						#if !defined(MULTITHREADED_SES_BUILDING)
						PointCell *pc1 = (PointCell *)sesComplex[ atomPatches[ec->id[0]] ];
						#else // MULTITHREADED_SES_BUILDING
						ThreadDataWrapper *tdw = &thread_data_wrapper[ it->first ];

						int cell_id = tdw->atomPatches[ec->id[0]];
						PointCell *pc1 = (PointCell *)tdw->sesComplex[ cell_id ];
						#endif // MULTITHREADED_SES_BUILDING

						int ref1 = ec->id[0];
						int ref2 = ec->id[1];
						int np = pc1->incidentProbes.size();

						// get all the incident probes. Most of the cases will be 2.
						for (int indf=0; indf<np; indf++)
						{
							#if !defined(MULTITHREADED_SES_BUILDING)
							FacetCell *f1 = (FacetCell *)sesComplex[ pc1->incidentProbes[indf] ];
							#else
							FacetCell *f1 = (FacetCell *)tdw->sesComplex[ pc1->incidentProbes[indf] ];
							#endif

							int i1 = f1->id[0];
							int i2 = f1->id[1];
							int i3 = f1->id[2];

							if ((i1 == ref1 && i2 == ref2) ||
								(i2 == ref1 && i1 == ref2) ||
								(i1 == ref1 && i3 == ref2) ||
								(i3 == ref1 && i1 == ref2) ||
								(i2 == ref1 && i3 == ref2) ||
								(i3 == ref1 && i2 == ref2))
							{
								fcv[index++] = f1;
							}
						}
						#endif // OPTIMIZE_CELL_STRUCTURE
					}
					double u[3], v[3];

					v[0] = ec->Rot[0][0];
					u[0] = ec->Rot[0][1];
					v[1] = ec->Rot[1][0];
					u[1] = ec->Rot[1][1];
					v[2] = ec->Rot[2][0];
					u[2] = ec->Rot[2][1];
					getCoi(ec->center,ec->rcoi,sampledPoints,coiNum,u,v);
					
					// for all the points get projection. Keep the nearest feasible
					for (int kk=0; kk<coiNum; kk++)
					{
						projectToSphere(p,sampledPoints[kk],probe_radius,locProj,dist);

						for (int ll=0; ll<index; ll++)
						{
							if (!isFeasible(fcv[ll],locProj))
								continue;
							
							if (dist < minDist)
							{
								minDist = dist;
								(*proj1) = locProj[0];
								(*proj2) = locProj[1];
								(*proj3) = locProj[2];
								(*normal1) = 0;
								(*normal2) = 0;
								(*normal3) = 0;
								fixed = true;
							}
						}
					}
				}
			}
		} // end singular projection
		deleteMatrix2D<double>(coiNum,3,sampledPoints);
		
		if (!fixed)
		{
			(*proj1) = p[0];
			(*proj2) = p[1];
			(*proj3) = p[2];
			{
				#ifdef ENABLE_BOOST_THREADS
				boost::mutex::scoped_lock scopedLock(mutex);
				#endif
				(*errorStream) << endl << WARN << "Approximating bgp with grid point ";
			}
		}
		else
		{
			// non tangent continuity has been properly managed.
			// the obtained projection is the best thing allowed by Connolly definition.
		}
	}
	return true;
}


/** check the orientation. Assume the planes points toward the visible region of the torus*/
bool ConnollySurface::orientation(double *pb_center1, double *pb_center2, double *w1, double *w2)
{
	// intersect the lines starting at the probes position and following the normal orientations
	// double a = DOT(w1,w1);
	double b = DOT(w1,w2);
	double c = DOT(w2,w2);
	double w0[3];

	SUB(w0,pb_center1,pb_center2)

	double d = DOT(w1,w0);
	double e = DOT(w2,w0);

	// double sc = (b*e-c*d)/(a*c-b*b);
	// only sign
	// double sc = (b*e-c*d);

	// cout << endl << "sc param " << sc;
	if (b*e < c*d)
		return true;
	else 
		return false;
}


void ConnollySurface::sortProbes(EdgeCell *ec, FacetCell **fcv, int np, int *sorted)
{
	vector<indexed_double> angles;
	double ref[3];
	double tt;
	double ws[3],temp[3],vs[3];

	SUB(ref,fcv[0]->center,ec->center)
	NORMALIZE_S(ref,tt);
	SUB(temp,fcv[1]->center,ec->center)
	NORMALIZE_S(temp,tt);
	CROSS(vs,ref,temp)
	CROSS(ws,vs,ref)

	for (int i=1; i<np; i++)
	{
		double temp[3];
		SUB(temp,fcv[i]->center,ec->center)
		NORMALIZE_S(temp,tt)
		double a1 = DOT(ws,temp);
		double a2 = DOT(ref,temp);
		double angle = atan2(a1,a2);
		// [0,2pi]
		if (angle < 0)
			angle = TWO_PI+angle;
		// cout << endl << "Angle " << angle;
		angles.push_back(indexed_double(angle,i));
	}
	// remember indexing and sort
	sort(angles.begin(),angles.end(),index_double_comparator);
	// save
	sorted[0] = 0;
	// cout << endl << "Ordered";
	for (int i=1; i<np; i++)
	{
		sorted[i] = angles[i-1].second;
		// cout << endl << "\t" << sorted[i];
	}
	return;
}


void ConnollySurface::getCoi(double *torus_center,double rcoi,double **sampledPoints,int numPoints,double *u,double *v)
{
	// FILE *fp;
	// fp = fopen("tempp.txt","a");

	double step = TWO_PI/numPoints;
	double tcoi = 0.;

	for (int i=0; i<numPoints; i++)
	{
		sampledPoints[i][0] = torus_center[0] + rcoi*(cos(tcoi)*u[0] + sin(tcoi)*v[0]);
		sampledPoints[i][1] = torus_center[1] + rcoi*(cos(tcoi)*u[1] + sin(tcoi)*v[1]);
		sampledPoints[i][2] = torus_center[2] + rcoi*(cos(tcoi)*u[2] + sin(tcoi)*v[2]);
		tcoi += step;

		// fprintf(fp, "\nsphere {");
		// fprintf(fp, "\n<%f,%f,%f>,0.1\n", sampledPoints[i][0],sampledPoints[i][1],sampledPoints[i][2]);
		// fprintf(fp, "\npigment{ color White}}");
	}
	// fclose(fp);
}


void ConnollySurface::saveConcaveSpherePatch(ofstream &of, FacetCell *fc, int i)
{
	char buff2[BUFLEN],buff[BUFLEN],temp[BUFLEN];
	sprintf(buff,"Mesh_CS%d",i);
	of << "\n\n #declare " << buff << "=" ;

	sprintf(buff2,"\n\n intersection { \n ");
	of << endl << buff2;

	for (int l=0;l<4;l++)
	{
		sprintf(temp, "plane{<%f,%f,%f>,%f}", fc->planes[l][0],fc->planes[l][1],fc->planes[l][2],-(fc->planes[l][3])/(sqrt((DOT(fc->planes[l],fc->planes[l])))));
		of << endl << temp;
	}

	// cout << endl << "Num " << fc->numSelfIntersections;
	if (fc->isSelfIntersecting)
	{
		of << endl << "// self intersection planes" << endl;
		for (unsigned int l=0; l<fc->self_intersection_planes.size(); l++)
		{
			#if !defined(OPTIMIZE_CELL_STRUCTURE)
			// cout << endl << fc->self_intersection_planes[l][0] << " " << fc->self_intersection_planes[l][1] <<" "<< fc->self_intersection_planes[l][2]
			sprintf(temp, "plane{<%f,%f,%f>,%f}", fc->self_intersection_planes[l][0],fc->self_intersection_planes[l][1],fc->self_intersection_planes[l][2],-(fc->self_intersection_planes[l][3])/(sqrt((DOT(fc->self_intersection_planes[l],fc->self_intersection_planes[l])))));
			#else
			if (fc->self_intersection_plane_labels[l])
			{
				// cout << endl << fc->self_intersection_planes[l][0] << " " << fc->self_intersection_planes[l][1] <<" "<< fc->self_intersection_planes[l][2]
				sprintf(temp, "plane{<%f,%f,%f>,%f}", fc->self_intersection_planes[l][0],fc->self_intersection_planes[l][1],fc->self_intersection_planes[l][2],-(fc->self_intersection_planes[l][3])/(sqrt((DOT(fc->self_intersection_planes[l],fc->self_intersection_planes[l])))));
			}
			else
			{
				// cout << endl << -fc->self_intersection_planes[l][0] << " " << -fc->self_intersection_planes[l][1] <<" "<< -fc->self_intersection_planes[l][2]
				sprintf(temp, "plane{<%f,%f,%f>,%f}", -fc->self_intersection_planes[l][0],-fc->self_intersection_planes[l][1],-fc->self_intersection_planes[l][2],+(fc->self_intersection_planes[l][3])/(sqrt((DOT(fc->self_intersection_planes[l],fc->self_intersection_planes[l])))));
			}
			#endif
			of << endl << temp;
		}
	}
	// cout << endl << "-------";

	of << "\n}";
	// save the sphere with its radius
	sprintf(buff2, "\n\n sphere { \n <%lf,%lf,%lf>,%lf", fc->center[0],fc->center[1],fc->center[2],probe_radius);
	of << buff2;
	// if (fc->patch_type == SINGULAR_FACE_CELL)
	// 	of << "\n pigment{ color Orange}";
	// else
	of << "\n pigment{ color Red}";
	of << "\n bounded_by { " << buff << " } \n clipped_by{ bounded_by}}";
	}


	void ConnollySurface::saveSphere(ostream &of, double *center, double radius)
	{
		char buff2[BUFLEN];
		sprintf(buff2, "\n\n sphere { \n <%lf,%lf,%lf>,%lf  pigment{color Green}}", center[0],center[1],center[2],radius);
		of << buff2;
		}


		void ConnollySurface::saveAtomPatch(ofstream &of, PointCell *pc)
		{
			char buff2[BUFLEN],buff[BUFLEN];

			of << "\n// -------------- ATOM " << pc->id << "-------------------";
			sprintf(buff, "AtomClip%d", pc->id);
			// add
			of << "\n\n #declare " << buff << " =" ;
			of << "\n union {";
			double *sphere;
			double radius;

			#if !defined(OPTIMIZE_CELL_STRUCTURE)
			for (unsigned int k=0; k<pc->neighbours.size(); k++)
			{
				sphere = pc->neighbours[k]->clipping_center;
				radius = pc->neighbours[k]->clipping_radius;
				sprintf(buff2, "\n\n sphere { \n <%lf,%lf,%lf>,%lf}", sphere[0],sphere[1],sphere[2],radius);
				of << buff2;
			}
			for (unsigned int k=0; k<pc->buried_neighbours.size(); k++)
			{
				sphere = pc->buried_neighbours[k]->clipping_center;
				radius = pc->buried_neighbours[k]->clipping_radius;
				sprintf(buff2, "\n\n sphere { \n <%lf,%lf,%lf>,%lf}", sphere[0],sphere[1],sphere[2],radius);
				of << buff2;
			}
			#else
			for (unsigned int k=0; k<pc->neighbour_data.size(); k += 4)
			{
				sphere = &pc->neighbour_data[k];
				radius = sqrt(pc->neighbour_data[k+3]);
				sprintf(buff2, "\n\n sphere { \n <%lf,%lf,%lf>,%lf}", sphere[0],sphere[1],sphere[2],radius);
				of << buff2;
			}
			for (unsigned int k=0; k<pc->buried_neighbour_data.size(); k += 4)
			{
				sphere = &pc->buried_neighbour_data[k];
				radius = sqrt(pc->buried_neighbour_data[k+3]);
				sprintf(buff2, "\n\n sphere { \n <%lf,%lf,%lf>,%lf}", sphere[0],sphere[1],sphere[2],radius);
				of << buff2;
			}
			#endif

			of << "\n}";

			of << "\ndifference{";
			// sprintf(buff2, "\n\n sphere { \n <%lf,%lf,%lf>,%lf  pigment{color Green}}", delphi->atoms[pc->id]->pos[0],delphi->atoms[pc->id]->pos[1],delphi->atoms[pc->id]->pos[2],delphi->atoms[pc->id]->radius);
			sprintf(buff2, "\n\n sphere { \n <%lf,%lf,%lf>,%lf  pigment{color Green}}", delphi->atoms[pc->id].pos[0],delphi->atoms[pc->id].pos[1],delphi->atoms[pc->id].pos[2],delphi->atoms[pc->id].radius);
			of << buff2;
			of << "\nobject{" << buff << "}}";
			}


			void ConnollySurface::saveEdgePatch(ofstream &of, EdgeCell *ec, int size, double bigR, double *u, double *v, double *w, bool isComplex)
			{
				char buff2[3*BUFLEN],buff[3*BUFLEN],buff3[3*BUFLEN],clipplanes[BUFLEN],clip[BUFLEN];

				of << "\n// ------------- torus -------------//";
				of << "\n// pair " << ec->id[0] << "," << ec->id[1];
				of << "\n\n";
				sprintf(clipplanes,"Clipping_Planes%d", size);

				// singular edge cell has not clipping planes
				if (ec->patch_type != SINGULAR_EDGE_CELL)
				{
					// add clipping planes
					of << "\n\n #declare " << clipplanes << " =" ;

					if (ec->additional_planes.size() == 0)
					{
						if (ec->acute)
							sprintf(buff3, "\n\n // clipping planes \n\n intersection { \n ");
						else
							sprintf(buff3, "\n\n // clipping planes \n\n union { \n ");

						of << buff3;
						for (int l=0; l<2; l++)
						{
							#if !defined(OPTIMIZE_CELL_STRUCTURE)
							sprintf(buff3, "\nplane{<%f,%f,%f>,%f}", ec->cutting_planes[l][0],ec->cutting_planes[l][1],ec->cutting_planes[l][2],-(ec->cutting_planes[l][3])/(sqrt((DOT(ec->cutting_planes[l],ec->cutting_planes[l])))));
							of << buff3;
							#else
							sprintf(buff3, "\nplane{<%f,%f,%f>,%f}", -ec->cutting_planes[l][0],-ec->cutting_planes[l][1],-ec->cutting_planes[l][2],ec->cutting_planes[l][3]/(sqrt((DOT(ec->cutting_planes[l],ec->cutting_planes[l])))));
							of << buff3;
							#endif
						}
						of << "\n}";
					}
					else
					{
						of << "\nunion{\n";
						int jj=0;
						for (unsigned int i=0; i<ec->additional_planes.size(); i+=2)
						{
							if (ec->flags[jj])
								sprintf(buff3,"\n\n // clipping planes \n\n intersection { \n ");
							else
								sprintf(buff3,"\n\n // clipping planes \n\n union { \n ");

							jj++;

							of << buff3;
							#if !defined(OPTIMIZE_CELL_STRUCTURE)
							sprintf(buff3, "\nplane{<%f,%f,%f>,%f}", ec->additional_planes[i][0],ec->additional_planes[i][1],ec->additional_planes[i][2],-(ec->additional_planes[i][3])/(sqrt((DOT(ec->additional_planes[i],ec->additional_planes[i])))));
							of << buff3;
							sprintf(buff3, "\nplane{<%f,%f,%f>,%f}", ec->additional_planes[i+1][0],ec->additional_planes[i+1][1],ec->additional_planes[i+1][2],-(ec->additional_planes[i+1][3])/(sqrt((DOT(ec->additional_planes[i+1],ec->additional_planes[i+1])))));
							#else
							sprintf(buff3, "\nplane{<%f,%f,%f>,%f}", -ec->additional_planes[i][0],-ec->additional_planes[i][1],-ec->additional_planes[i][2],ec->additional_planes[i][3]/(sqrt((DOT(ec->additional_planes[i],-ec->additional_planes[i])))));
							of << buff3;
							sprintf(buff3, "\nplane{<%f,%f,%f>,%f}", -ec->additional_planes[i+1][0],-ec->additional_planes[i+1][1],-ec->additional_planes[i+1][2],ec->additional_planes[i+1][3]/(sqrt((DOT(ec->additional_planes[i+1],-ec->additional_planes[i+1])))));
							#endif
							of << buff3;

							of << "\n}";
						}
						of << "\n}";
					}
				}

				sprintf(clip,"Clip%d",size);
				of << "\n\n #declare " << clip << " =" ;

				// manage self intersection if necessary
				if (bigR < probe_radius)
				{
					sprintf(buff2, "\n\n difference { \n ""sphere {  <%lf,%lf,%lf>,%lf }  ""\n sphere {  <%lf,%lf,%lf>,%lf } }",
							// (clipping_center.x()),(clipping_center.y()),(clipping_center.z()),r_clipping,
							//torus_center.x(),torus_center.y(),torus_center.z(),sqrt(-bigR*bigR+probe_radius*probe_radius));
							(ec->clipping_center[0]),(ec->clipping_center[1]),(ec->clipping_center[2]),ec->clipping_radius,
							ec->center[0],ec->center[1],ec->center[2],sqrt(-bigR*bigR+probe_radius*probe_radius));
					of << buff2;
				}
				else
				{
					sprintf(buff2, "\nsphere {  <%lf,%lf,%lf>,%lf }  ",
							// (clipping_center.x()),(clipping_center.y()),(clipping_center.z()),r_clipping);
							(ec->clipping_center[0]),(ec->clipping_center[1]),(ec->clipping_center[2]),ec->clipping_radius);
					of << buff2;
				}

				// global clipping object
				sprintf(buff, "Clipping_Object%d", getNumPatches());

				// add all
				of << "\n\n #declare " << buff << " =" ;
				// torus,clipping spheres, and also clipping planes
				if (ec->patch_type != SINGULAR_EDGE_CELL)
					sprintf(buff2, "\n\n intersection { \n"
					"object{%s}\nobject{%s}}", clipplanes,clip);
				// only torus with clipping spheres
				else
					sprintf(buff2, "\n\nobject{%s}", clip);

				of << buff2;


				sprintf(buff2, "\n\n torus { \n %lf,%lf", bigR,probe_radius);
				of << buff2;
				// if (!isComplex)
				// 	of << "\n pigment{ color White}";
				// else
				// 	of << "\n pigment{ color Yellow}";

				// sprintf(buff2,"\n matrix <%lf,%lf,%lf,\n%lf,%lf,%lf,\n%lf,%lf,%lf,\n0,0,0>", U.x(),U.y(),U.z(), w.x(),w.y(),w.z(), V.x(),V.y(),V.z());

				sprintf(buff2, "\n matrix <%lf,%lf,%lf,\n%lf,%lf,%lf,\n%lf,%lf,%lf,\n0,0,0>",
						u[0],u[1],u[2], w[0],w[1],w[2], v[0],v[1],v[2]);
				of << buff2;
				// of << "\n translate"<< "<" << torus_center.x() << "," << torus_center.y() << "," << torus_center.z() << ">";
				of << "\n translate"<< "<" << ec->center[0] << "," << ec->center[1] << "," << ec->center[2] << ">";
				of << "\nsturm";
				of << "\n bounded_by { " << buff << " } \n clipped_by{ bounded_by}}";

				// to draw the clipping sphere
				// sprintf(buff2, "\nsphere {  <%lf,%lf,%lf>,%lf   \n pigment{ color Yellow filter 0.7}}",
				// (ec->clipping_center[0]),(ec->clipping_center[1]),(ec->clipping_center[2]),ec->clipping_radius);
				// of << buff2;
				}


bool ConnollySurface::save(char *fileName)
{
	return false;
}


bool ConnollySurface::load(char *fileName)
{
	return false;
}


void ConnollySurface::printSummary()
{
	cout << endl << INFO << "Probe Radius value " << getProbeRadius();

	if (getNumPatches() == 0)
	{
		cout << endl << WARN << "Connolly surface not loaded!";
	}
	else
	{
		cout << endl << INFO << "Number of ses cells -> " << getNumPatches();
		cout << endl << INFO << "Number of del_point cells -> " << type[POINT_CELL];
		cout << endl << INFO << "Number of regular del_edge cells -> " << type[REGULAR_EDGE_CELL];
		cout << endl << INFO << "Number of singular del_edge cells -> " << type[SINGULAR_EDGE_CELL];
		cout << endl << INFO << "Number of regular del_facet cells -> " << type[REGULAR_FACE_CELL];
		cout << endl << INFO << "Number of singular del_facet cells -> " << type[SINGULAR_FACE_CELL];
		/*
		if (internals != NULL)
		{
			*internals << endl << "cells " << getNumPatches();
			*internals << endl << "del_point " << type[POINT_CELL];
			*internals << endl << "rdel_edge " << type[REGULAR_EDGE_CELL];
			*internals << endl << "sdel_edge " << type[SINGULAR_EDGE_CELL];
			*internals << endl << "rdel_facet " << type[REGULAR_FACE_CELL];
			*internals << endl << "sdel_facet " << type[SINGULAR_FACE_CELL];
		}
		*/
	}
}
