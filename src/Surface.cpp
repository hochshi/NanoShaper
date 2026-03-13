
//---------------------------------------------------------
/**    @file		Surface.cpp
*     @brief	Surface.cpp is the Surface CLASS
*															*/
//---------------------------------------------------------

#include "Surface.h"
#include <fstream>

#ifdef PLY_ENABLED
#include "tinyply/tinyply.h"
#endif

static inline VERTEX_TYPE *sanitizeNormalPtrForThread(
	vector<VERTEX_TYPE> (&normalsBuffers)[MAX_TASKS_TIMES_THREADS],
	int thread_id,
	bool computeNormals,
	bool providesAnalyticalNormals,
	VERTEX_TYPE *normal)
{
	if (normal == NULL || !(computeNormals && providesAnalyticalNormals))
		return NULL;
	if (normalsBuffers[thread_id].empty())
		return NULL;
	VERTEX_TYPE *begin = &normalsBuffers[thread_id][0];
	VERTEX_TYPE *end = begin + normalsBuffers[thread_id].size();
	if (normal < begin || (normal + 2) >= end)
		return NULL;
	return normal;
}

void Surface::init()
{
	delphi = NULL;
	panel = 0;
	inside = 1;
	accurateTriangulation = true;
	fillCavitiesFlag = false;
	isAvailableScalarField = false;
	projBGP = false;
	scalarField = NULL;
	delta_accurate_triangulation = DELTA;
	checkDuplicatedVertices = false;
	wellShaped = false;
	probe_radius = 1.4;

	#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
	intersectionsMatrixAlongX = NULL;
	intersectionsMatrixAlongY = NULL;
	intersectionsMatrixAlongZ = NULL;
	#if !defined(AVOID_NORMALS_MATRICES)
	normalsMatrixAlongX = NULL;
	normalsMatrixAlongY = NULL;
	normalsMatrixAlongZ = NULL;
	#endif
	#else
	bilevel_intersectionsMatrixAlongX = NULL;
	bilevel_intersectionsMatrixAlongY = NULL;
	bilevel_intersectionsMatrixAlongZ = NULL;
	#endif

	bgp_type = NULL;
	gridMultiMap = NULL;
	#if !defined(USE_COMPRESSED_GRIDS)
	verticesInsidenessMap = NULL;
	#endif
	compressed_verticesInsidenessMap = NULL;
	sternLayer = -1;
	isRCbased = true;
	randDisplacement = RAND_DISPLACEMENT;
	gridLoad = NULL;
	useLoadBalancing = true;
	totalLoad = 0;
	vertexAtomsMap = NULL;
	vertexAtomsMapFlag = false;
	computeNormals = false;
	saveMSMS = false;
	savePLY = false;
	providesAnalyticalNormals = false;
	#if !defined(USE_COMPRESSED_GRIDS)
	activeCubes = NULL;
	#endif
	compressed_activeCubes = NULL;
	patchBasedAlgorithm = true;
	analyticalTorusIntersectionAlgorithm = true;
	
	if (patchBasedAlgorithm)
	{
		num_pixel_intersections = NULL;
		pixel_intersections = NULL;
	}
}


void Surface::init(ConfigFile *cf)
{
	bool projBGP = cf->read<bool>("Project_boundary_grid_points", false);
	bool accTri = cf->read<bool>("Accurate_Triangulation", false);
	bool doTri = cf->read<bool>("Triangulation", false);
	bool checkDuplicatedVertices = cf->read<bool>( "Check_duplicated_vertices", true);
	bool wellShaped = cf->read<bool>("Keep_Water_Shaped_Cavities", false);
	double probeRadius = cf->read<double>("Probe_Radius", 1.4);
	bool lb = cf->read<bool>("Load_Balancing", true);
	bool vaFlag = cf->read<bool>( "Vertex_Atom_Info", false);
	bool computeNormals = cf->read<bool>("Compute_Vertex_Normals", false);
	bool saveMSMS = cf->read<bool>("Save_Mesh_MSMS_Format", false);
	bool savePLY = cf->read<bool>("Save_Mesh_PLY_Format", false);
	double sternLayer = cf->read<double>("Stern_layer", -1.);
	rootFile = cf->read<string>("Root_FileName","");
	bool patch_based = cf->read<bool>("Patch_Based_Algorithm", true);
	bool analytical_torus_intersection = cf->read<bool>("Analytical_Ray_Vs_Torus_Intersection", true);
	bool collect_face_intersections = cf->read<bool>("Collect_Face_Intersections", false);
	double parallel_buildup_halo_thickness = cf->read<double>("Parallel_Buildup_Halo_Thickness", 0.0);
	bool serial_build = cf->read<bool>("Force_Serial_Build", false);
	int max_atoms = cf->read<int>("Max_Num_Atoms", -1);
	double domain_shrinkage = cf->read<double>("Domain_Shrinkage", 0.);

	setProjBGP(projBGP);
	setAccurateTriangulationFlag(accTri);
	setTriangulationFlag(doTri);
	setCheckDuplicatedVertices(checkDuplicatedVertices);
	setKeepWellShapedCavities(wellShaped);
	setProbeRadius(probeRadius);
	setLoadBalancing(lb);
	setVertexAtomsMap(vaFlag);
	setComputeNormals(computeNormals);
	setSaveMSMS(saveMSMS);
	setSavePLY(savePLY);
	setRayTracingAlgorithm(patch_based);
	setRayVsTorusIntersectionAlgorithm(analytical_torus_intersection);
	setCollectFaceIntersections(collect_face_intersections);
	setParallelBuildupHaloThickness(parallel_buildup_halo_thickness);
	setForceSerialBuild (serial_build);
	setMaxNumAtoms (max_atoms);
	setDomainShrinkage (domain_shrinkage);
	
	// if >0 enable stern layer, else disabled by default
	if (sternLayer > 0)
		setSternLayer(sternLayer);

	if (parallel_buildup_halo_thickness != 0.0)
		conf.parallelBuildupHaloThicknessInput = true;
	else
		conf.parallelBuildupHaloThicknessInput = false;
}


Surface::Surface()
{
	init();
}


Surface::Surface(ConfigFile *cf)
{
	init();
	init(cf);
}


void Surface::clear()
{
	if (triList.size() > 0)
	{
		triList.clear();
	}
	if (vertList.size() > 0)
	{
		vertList.clear();
	}
	if (normalsList.size() > 0)
	{
		normalsList.clear();
	}
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	for (int i=0; i < conf.numThreads; i++)
	{
		verticesBuffers[i].clear();
		normalsBuffers[i].clear();
	}
	#endif

	deleteIntersectionsMatrices();
	
	#if !defined(USE_COMPRESSED_GRIDS)
	if (!optimizeGrids)
	{
		if (verticesInsidenessMap != NULL)
			deleteMatrix3D<bool>(delphi->nx, delphi->ny, delphi->nz, verticesInsidenessMap);
	}
	else
	#endif
	{
		if (compressed_verticesInsidenessMap != NULL)
			deleteVector<unsigned int>(compressed_verticesInsidenessMap);
	}

	if (scalarField != NULL)
		deleteMatrix3D<double>(delphi->nx,delphi->ny,delphi->nz,scalarField);

	if (gridLoad != NULL)
		deleteVector<int>(gridLoad);

	if (vertexAtomsMap != NULL)
		deleteVector<int>(vertexAtomsMap);

	init();
}


Surface::~Surface()
{
	clear();
}


int Surface::getNumPatches (void)
{
	cout << endl << WARN << "Number of patches not available!";
	return 0;
}


bool Surface::isPatchBasedRayTracingSupported (void)
{
	cout << endl << ERR << "This surface type does not support patch based ray-tracing";
	exit(-1);
}


bool Surface::build()
{
	cout << endl << WARN << "Build surface not supported!";
	return false;
}


bool Surface::save(char *fileName)
{
	cout << endl << WARN << "Save surface not supported!";
	return false;
}


bool Surface::load(char *fileName)
{
	cout << endl << WARN << "Load surface not supported";
	return false;
}


void Surface::printSummary()
{
	cout << endl << WARN << "Print summary not supported!";
}


void Surface::allocIntersectionsMatrices(int octree_side_size)
{
	deleteIntersectionsMatrices();

	#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)

	intersectionsMatrixAlongX = new Octree<int>(octree_side_size,-1);

	intersectionsMatrixAlongY = new Octree<int>(octree_side_size,-1);

	intersectionsMatrixAlongZ = new Octree<int>(octree_side_size,-1);

	#else // OPTIMIZE_INTERSECTIONS_MANAGEMENT

	bilevel_intersectionsMatrixAlongX = allocateBilevelGridCells<int>(delphi->nx, delphi->ny, delphi->nz);

	bilevel_intersectionsMatrixAlongY = allocateBilevelGridCells<int>(delphi->nx, delphi->ny, delphi->nz);

	bilevel_intersectionsMatrixAlongZ = allocateBilevelGridCells<int>(delphi->nx, delphi->ny, delphi->nz);

	#endif // OPTIMIZE_INTERSECTIONS_MANAGEMENT
}


void Surface::deleteIntersectionsMatrices(void)
{
	#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)

	if (intersectionsMatrixAlongX != NULL)
	{
		for (int k=0; k<delphi->nz; k++)
			for (int j=0; j<delphi->ny; j++)
				for (int i=0; i<delphi->nx; i++)
					if (intersectionsMatrixAlongX->at(i,j,k) != -1)
						intersectionsMatrixAlongX->erase(i,j,k);
		delete intersectionsMatrixAlongX;
		intersectionsMatrixAlongX = NULL;
	}
	if (intersectionsMatrixAlongY != NULL)
	{
		for (int k=0; k<delphi->nz; k++)
			for (int j=0; j<delphi->ny; j++)
				for (int i=0; i<delphi->nx; i++)
					if (intersectionsMatrixAlongY->at(i,j,k) != -1)
						intersectionsMatrixAlongY->erase(i,j,k);
		delete intersectionsMatrixAlongY;
		intersectionsMatrixAlongY = NULL;
	}
	if (intersectionsMatrixAlongZ != NULL)
	{
		for (int k=0; k<delphi->nz; k++)
			for (int j=0; j<delphi->ny; j++)
				for (int i=0; i<delphi->nx; i++)
					if (intersectionsMatrixAlongZ->at(i,j,k) != -1)
						intersectionsMatrixAlongZ->erase(i,j,k);
		delete intersectionsMatrixAlongZ;
		intersectionsMatrixAlongZ = NULL;
	}

	#else // OPTIMIZE_INTERSECTIONS_MANAGEMENT

	if (bilevel_intersectionsMatrixAlongX != NULL)
	{
		deleteBilevelGridCells<int>(bilevel_intersectionsMatrixAlongX, delphi->nx, delphi->ny, delphi->nz);
		deleteVector<int *>(bilevel_intersectionsMatrixAlongX);
	}
	if (bilevel_intersectionsMatrixAlongY != NULL)
	{
		deleteBilevelGridCells<int>(bilevel_intersectionsMatrixAlongY, delphi->nx, delphi->ny, delphi->nz);
		deleteVector<int *>(bilevel_intersectionsMatrixAlongY);
	}
	if (bilevel_intersectionsMatrixAlongZ != NULL)
	{
		deleteBilevelGridCells<int>(bilevel_intersectionsMatrixAlongZ, delphi->nx, delphi->ny, delphi->nz);
		deleteVector<int *>(bilevel_intersectionsMatrixAlongZ);
	}
	#endif // OPTIMIZE_INTERSECTIONS_MANAGEMENT
}


#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT) && !defined(AVOID_NORMALS_MATRICES)

void Surface::allocNormalsMatrices(int octree_side_size)
{
	deleteNormalsMatrices();

	normalsMatrixAlongX = new Octree<int>(octree_side_size,-1);

	normalsMatrixAlongY = new Octree<int>(octree_side_size,-1);

	normalsMatrixAlongZ = new Octree<int>(octree_side_size,-1);
}


void Surface::deleteNormalsMatrices(void)
{
	if (normalsMatrixAlongX != NULL)
	{
		for (int k=0; k<delphi->nz; k++)
			for (int j=0; j<delphi->ny; j++)
				for (int i=0; i<delphi->nx; i++)
					if (normalsMatrixAlongX->at(i,j,k) != -1)
						normalsMatrixAlongX->erase(i,j,k);
		delete normalsMatrixAlongX;
		normalsMatrixAlongX = NULL;
	}
	if (normalsMatrixAlongY != NULL)
	{
		for (int k=0; k<delphi->nz; k++)
			for (int j=0; j<delphi->ny; j++)
				for (int i=0; i<delphi->nx; i++)
					if (normalsMatrixAlongY->at(i,j,k) != -1)
						normalsMatrixAlongY->erase(i,j,k);
		delete normalsMatrixAlongY;
		normalsMatrixAlongY = NULL;
	}
	if (normalsMatrixAlongZ != NULL)
	{
		for (int k=0; k<delphi->nz; k++)
			for (int j=0; j<delphi->ny; j++)
				for (int i=0; i<delphi->nx; i++)
					if (normalsMatrixAlongZ->at(i,j,k) != -1)
						normalsMatrixAlongZ->erase(i,j,k);
		delete normalsMatrixAlongZ;
		normalsMatrixAlongZ = NULL;
	}
}

#endif // OPTIMIZE_INTERSECTIONS_MANAGEMENT && AVOID_NORMALS_MATRICES


/** build epsmap and compute boundary grid points by intersection
 and projection routines. The minimal number of projections is performed.
 The ray tracing part is parallelized with boost threading if enabled.
 Ray tracing is implemented in two manners: the conventional, ray-centric
 way and the newer patch-based one, in which each patch is processed by
 one thread which pierces it with the rays spanned by its bounding box.
 
 // SD: PB_NEW: modified to allow returning intersections. */
bool Surface::getSurf(double *surf_volume, bool optimize_grids, bool fillCav, double cav_vol, vector<packet> *intersectionsInfo)
{
	int num_threads = conf.numThreads;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	double volPanel[3] = {0,0,0};

	// These will label the taken conditionals during RT and will be useful for volume calculations
	panelVolumeFlag[0][0] = 0;
	panelVolumeFlag[0][1] = 0;
	panelVolumeFlag[1][0] = 0;
	panelVolumeFlag[1][1] = 0;
	panelVolumeFlag[2][0] = 0;
	panelVolumeFlag[2][1] = 0;

	setOptimizeGrids (optimize_grids);

	fillCavitiesFlag = fillCav;

	if (delphi == NULL)
	{
		cout << endl << WARN << "Cannot get surface without DelPhi environment!";
		return false;
	}

	// per thread independent data structures. This avoids resorting to
	// a mutex on the octree; thus, it is faster and simpler

	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *buffersIntersections;
	vector<coordVec> *buffersNormals;
	#else
	// intersection and normal details are packed in one structure
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *buffersIntersections;
	#else
	vector<compressedCoordNormPacket> *buffersIntersections;
	#endif
	#endif


	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *buffersIntersections_grid;
	vector<coordVec> *buffersNormals_grid;
	#else
	// intersection and normal details are packed in one structure
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *buffersIntersections_grid;
	#else
	vector<compressedCoordNormPacket> *buffersIntersections_grid;
	#endif
	#endif

	int *netInts_per_thd;


	// if Ray casting based perform parallel ray casting
	// if does not assume that grids colouring took place into the build phase
	if (isRCbased)
	{
		// clean if required
		triList.clear();
		vertList.clear();
		normalsList.clear();

		#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
		int octree_side_size = 2;
		bool found_label = 0;

		while (!found_label)
		{
			if (octree_side_size < NX ||
				octree_side_size < NY ||
				octree_side_size < NZ)
			{
				octree_side_size *= 2;
			}
			else
			{
				found_label = 1;
			}
		}
		// this is always needed. Vertices are stored here in any case
		allocIntersectionsMatrices(octree_side_size);

		#if !defined(AVOID_NORMALS_MATRICES)
		// this is used only if requested and are analytically computed
		// if not analytically computed, the octree will not be used
		if (computeNormals && providesAnalyticalNormals)
			allocNormalsMatrices(octree_side_size);
		#endif

		#else // OPTIMIZE_INTERSECTIONS_MANAGEMENT

		// Vertices' indexes are stored with a bilevel grid
		allocIntersectionsMatrices();

		#endif // OPTIMIZE_INTERSECTIONS_MANAGEMENT

		
		#if !defined(USE_COMPRESSED_GRIDS)
		if (!optimizeGrids)
		{
			if (accurateTriangulation)
			{
				if (verticesInsidenessMap != NULL)
					deleteMatrix3D<bool>(last_nx, last_ny, last_nz, verticesInsidenessMap);
					// deleteVector<bool>(verticesInsidenessMap);

				// verticesInsidenessMap = allocateVector<bool>(NX*NY*NZ);
				verticesInsidenessMap = allocateMatrix3D<bool>(NX, NY, NZ);
			
				/*
				int64_t tot = NX*NY*NZ;
				for (int64_t i=0;i<tot;i++)
					verticesInsidenessMap[i] = true;
				*/

				for (int k=0; k<NZ; k++)
					for (int j=0; j<NY; j++)
						for (int i=0; i<NX; i++)
							verticesInsidenessMap[k][j][i] = true;
			}
			else
				verticesInsidenessMap = NULL;
		}
		else
		#endif
		{
			// the vertices' insideness map has a compressed version: 32 one-bit values are store in one 32-bit word
			if (accurateTriangulation)
			{
				// compressed version of the vertices insideness grid: 32 consecutive values are stored in one element
				if (compressed_verticesInsidenessMap != NULL)
					deleteVector<unsigned int>(compressed_verticesInsidenessMap);

				compressed_verticesInsidenessMap = allocate32xCompressedGrid(true, NX, NY, NZ);
			}
			else
				compressed_verticesInsidenessMap = NULL;
		}

		last_nx = NX;
		last_ny = NY;
		last_nz = NZ;
		
		int na, nb;

		if (!delphi->getMultiDiel())


		threadPanelVolume = allocateMatrix2D<double>(3, num_threads);
		threadFailedRays = allocateVector<int>(num_threads);
		threadTotalRays = allocateVector<int>(num_threads);

		// per thread independent data structures. This avoids resorting to
		// a mutex on the hierarchical data structure; thus, it is faster and simpler

		#if !defined(COORD_NORM_PACKING)
		buffersIntersections = new vector<coordVec> [num_threads];
		buffersNormals = new vector<coordVec> [num_threads];
		#else
		// intersection and normal details are packed in one structure
		#if !defined(COMPRESS_INTERSECTION_COORDS)
		buffersIntersections = new vector<coordNormPacket> [num_threads];
		#else
		buffersIntersections = new vector<compressedCoordNormPacket> [num_threads];
		#endif
		#endif

		// SD PB_NEW
		if (intersectionsInfo != nullptr)
		{
			#if !defined(COORD_NORM_PACKING)
			buffersIntersections_grid = new vector<coordVec> [num_threads];
			buffersNormals_grid = new vector<coordVec> [num_threads];
			#else
			#if !defined(COMPRESS_INTERSECTION_COORDS)
			buffersIntersections_grid = new vector<coordNormPacket> [num_threads];
			#else
			buffersIntersections_grid = new vector<compressedCoordNormPacket> [num_threads];
			#endif
			#endif
		}

		netInts_per_thd = allocateVector<int>(num_threads);

		for (int j=0; j < num_threads; j++)
			netInts_per_thd[j] = 0;

		for (int j=0; j < num_threads; j++)
		{
			verticesBuffers[j].reserve( max(10, (int)(max(30., 15.*conf.scale) * 3.*max(1000,getNumPatches()) / (double)num_threads)) );

			if (computeNormals && providesAnalyticalNormals)
			{
				normalsBuffers[j].reserve( max(10, (int)(max(30., 15.*conf.scale) * 3.*max(1000,getNumPatches()) / (double)num_threads)) );
			}
		}
		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		if (computeNormals && !providesAnalyticalNormals)
		{
			// Only the first thread manages the normals in approximateNormals()
			normalsBuffers[0].reserve( max(10, (int)(max(30., 15.*conf.scale) * 3.*getNumPatches())) );
		}
		#endif

		auto chrono_start = chrono::high_resolution_clock::now();

		int numFails = 0;
		int numTotalRays = 0;

		if (!patchBasedAlgorithm)
		{
			cout.flush();
			cout << endl << INFO << "Conventional ray-based ray-tracing ";

			for (int i=0; i<num_threads; i++)
			{
				buffersIntersections[i].reserve(1000);
				#if !defined(COORD_NORM_PACKING)
				buffersNormals[i].reserve(1000);
				#endif
			}

			if (intersectionsInfo != nullptr)
			{
				for (int i=0; i<num_threads; i++)
				{
					buffersIntersections_grid[i].reserve(1000);
					#if !defined(COORD_NORM_PACKING)
					buffersNormals_grid[i].reserve(1000);
					#endif
				}
			}
		}
		else
		{
			cout.flush();
			cout << endl << INFO << "Patch-based ray-tracing ";

			if (isPatchBasedRayTracingSupported() == false)
			{
				cout << endl << ERR << "This surface type does not support patch based ray-tracing";
				cout << endl << REMARK << "Please set Patch_Based_Algorithm = false";
				cout << endl;
				exit(-1);
			}

			// panels[0] is the number of ray tracing steps with rays centred in one way
			// panels[1] is the number of ray tracing steps with rays centred in another one
			int panels[] = {3,3};
			if (!collectFaceIntersections && !delphi->buildStatus && !delphi->buildEpsMap)
				panels[0] = 0;
			else if (!collectFaceIntersections && delphi->buildStatus && !delphi->buildEpsMap)
				panels[0] = 1;
			if (!accurateTriangulation || isAvailableScalarField)
				panels[1] = 0;

			int64_t nxyz[3] = {NX, NY, NZ};


			#if !defined(SINGLE_PASS_RT)
			int *potInts_per_thd;

			int potentialIntersections = 0;
			#endif
			int netIntersections = 0;

			#if !defined(SINGLE_PASS_RT)
			potInts_per_thd = allocateVector<int>(num_threads);
			#endif

			// initPatchBasedRayTracing() is the first of multiple stages devoted to build the intersections'
			// and normals' buffers needed for the RT routine setVerticesAndGridsWithIntersectionData();
			// see Surface.h for info
			initPatchBasedRayTracing (nxyz, panels);

			int old_num_threads = num_threads;
			num_threads = MIN(getNumPatches(), num_threads);

			#ifdef ENABLE_BOOST_THREADS
			boost::thread_group thdGroup;
			#endif

			#if !defined(SINGLE_PASS_RT)

			#if defined(MULTITHREADING)
			for (int i=0; i<num_threads; i++)
			{
				#ifdef ENABLE_BOOST_THREADS
				thdGroup.create_thread(boost::bind(&Surface::getPatchPreIntersectionData, this, nxyz, panels, i, potInts_per_thd));
				#else
				getPatchPreIntersectionData (nxyz, panels, i, potInts_per_thd);
				#endif
			}
			#ifdef ENABLE_BOOST_THREADS
			thdGroup.join_all();
			#endif
			#else
			for (int i=0; i<num_threads; i++)
			{
				getPatchPreIntersectionData (nxyz, panels, i, potInts_per_thd);
			}
			#endif
			for (int i=0; i<num_threads; i++)
				potentialIntersections += potInts_per_thd[i];

			allocateIntermediatePatchBuffers (nxyz, panels, potentialIntersections);

			#endif // SINGLE_PASS_RT

			#if defined(MULTITHREADING)
			for (int i=0; i<num_threads; i++)
			{
				#ifdef ENABLE_BOOST_THREADS
				thdGroup.create_thread(boost::bind(&Surface::getPatchIntersectionData, this, nxyz, panels, i, &netInts_per_thd[i]));
				#else
				getPatchIntersectionData (nxyz, panels, i, &netInts_per_thd[i]);
				#endif
			}
			#ifdef ENABLE_BOOST_THREADS
			thdGroup.join_all();
			#endif
			#else
			for (int i=0; i<num_threads; i++)
			{
				getPatchIntersectionData (nxyz, panels, i, &netInts_per_thd[i]);
			}
			#endif
			for (int i=0; i<num_threads; i++) {
				netIntersections += netInts_per_thd[i];
			}

			countPixelIntersections (nxyz, panels);

			/*
			if (computeNormals && providesAnalyticalNormals)
			{
				for (int i=0; i<num_threads; i++)
				{
					#ifdef ENABLE_BOOST_THREADS
					thdGroup.create_thread(boost::bind(&Surface::getPatchNormalsAtIntersections, this, nxyz, panels, i));
					#else
					getPatchNormalsAtIntersections (nxyz, panels, i);
					#endif
				}
				#ifdef ENABLE_BOOST_THREADS
				thdGroup.join_all();
				#endif
				#else
				for (int i=0; i<num_threads; i++)
				{
					getPatchNormalsAtIntersections (nxyz, panels, i);
				}
			}
			*/

			copyPatchIntersections (nxyz, panels);


			auto reordering_chrono_start = chrono::high_resolution_clock::now();

			for (int i=0; i<num_threads; i++)
			{
				int64_t N_MAX = MAX(NX, MAX(NY, NZ));

				#ifdef ENABLE_BOOST_THREADS
				thdGroup.create_thread(boost::bind(&Surface::reorderPatchIntersections, this, i, N_MAX*N_MAX*(panels[0]+panels[1])));
				#else
				reorderPatchIntersections (i, N_MAX*N_MAX*(panels[0]+panels[1]);
				#endif
			}
			#ifdef ENABLE_BOOST_THREADS
			thdGroup.join_all();
			#endif

			auto reordering_chrono_end = chrono::high_resolution_clock::now();

			chrono::duration<double> reordering_time = reordering_chrono_end - reordering_chrono_start;
			cout << endl << INFO << "Intersections' reordering phase time... ";
			printf ("%.4e [s]", reordering_time.count());


			// #if !defined(SINGLE_PASS_RT)
			// cout << "number of potential intersections: " << potentialIntersections << endl;
			// #endif
			// cout << "net number of intersections: " << netIntersections << endl;

			for (int i=0; i<num_threads; i++)
			{
				buffersIntersections[i].reserve(netInts_per_thd[i]);
				#if !defined(COORD_NORM_PACKING)
				buffersNormals[i].reserve(netInts_per_thd[i]);
				#endif

				#if defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				// intersections store pointers into verticesBuffers/normalsBuffers.
				// At high resolution, reallocation here would invalidate those pointers.
				int nints = netInts_per_thd[i];
				#if !defined(SINGLE_PASS_RT)
				// potential intersections are an upper-bound from the pre-pass;
				// use it as safer reserve to avoid pointer invalidation.
				if (potInts_per_thd != NULL && potInts_per_thd[i] > nints)
					nints = potInts_per_thd[i];
				#endif
				if (nints < 0)
					nints = 0;
				size_t requiredCoords = 3ULL * (size_t)nints;
				if (verticesBuffers[i].capacity() < requiredCoords)
					verticesBuffers[i].reserve(requiredCoords);
				if (computeNormals && providesAnalyticalNormals && normalsBuffers[i].capacity() < requiredCoords)
					normalsBuffers[i].reserve(requiredCoords);
				#endif
			}

			if (intersectionsInfo != nullptr)
			{
				for (int i=0; i<num_threads; i++)
				{
					buffersIntersections_grid[i].reserve(netInts_per_thd[i]);
					#if !defined(COORD_NORM_PACKING)
					buffersNormals_grid[i].reserve(netInts_per_thd[i]);
					#endif
				}
			}

			#if !defined(SINGLE_PASS_RT)
			deleteVector<int>(potInts_per_thd);
			#endif

			num_threads = old_num_threads;
		}

		// numPanels is the number of ray tracing steps with rays centred in one way or the other one
		int numPanels = 3;
		if (!collectFaceIntersections && !delphi->buildStatus && !delphi->buildEpsMap)
			numPanels = 0;
		else if (!collectFaceIntersections && delphi->buildStatus && !delphi->buildEpsMap)
			numPanels = 1;
		if (accurateTriangulation && !isAvailableScalarField)
			numPanels = 3;

		// Phase 1, ray trace from each coordinate plane if necessary
		for (panel=0; panel < numPanels; panel++)
		{
			for (int l=0; l<num_threads; l++)
			{
				threadFailedRays[l] = 0;
				threadTotalRays[l] = 0;
				threadPanelVolume[panel][l] = 0;
			}

			if (!patchBasedAlgorithm)
			{
				preProcessPanel();
			}

			cout.flush();
			cout << endl << INFO << "Ray-tracing panel " << panel << "...";

			// left YZ panel
			if (panel == 0)
			{
				na = NZ;
				nb = NY;
			}
			// front XY panel
			else if (panel == 1)
			{
				na = NY;
				nb = NX;
			}
			// bottom XZ panel
			else
			{
				na = NZ;
				nb = NX;
			}

			// trivial split
			int chunk = (int)((double)na / (double)num_threads);
			int rem = na%num_threads;
			#if defined(ENABLE_BOOST_THREADS)
			/*
			if (optimizeGrids && delphi->buildStatus)
			{
				// status map is a coarse grid of 4^3 voxels, so that conflict-free cross-thread writes
				// are possible if each thread proceeds launching ray slabs which are multiple of 4*nb rays
				chunk = (int)((double)na / (double)(4*num_threads));
				if (chunk == 0)
					chunk = 1;
				chunk *= 4;
				rem = na - chunk * num_threads;
				if (rem > 0)
					chunk += 4;
			}
			*/
			#endif

			#ifdef ENABLE_BOOST_THREADS
			boost::thread_group thdGroup;
			#endif
			
			if (!optimizeGrids || !delphi->buildStatus)
			{
				// setup split
				int start, stop;

				if (!useLoadBalancing)
				{
					for (int j=0; j<num_threads; j++)
					{
						if (j == 0)
						{
							start = 0;
							stop = chunk;
						}
						else
						{
							start = stop;
							stop = start+chunk;
						}
						if (j < rem)
							stop++;

						packet pack;
						pack.first = &buffersIntersections[j];
						#if !defined(COORD_NORM_PACKING)
						pack.second = &buffersNormals[j];
						#endif

						// SD PB_NEW
						packet pack_grid;
						if (intersectionsInfo != nullptr)
						{
							pack_grid.first = &buffersIntersections_grid[j];
							#if !defined(COORD_NORM_PACKING)
							pack_grid.second = &buffersNormals_grid[j];
							#endif
						}

						if (patchBasedAlgorithm)
						{
							#if defined(ENABLE_BOOST_THREADS)
							thdGroup.create_thread(boost::bind(&Surface::setVerticesAndGridsWithIntersectionData, this,
															   j, nb, start, stop, 1, 1, pack, pack_grid));
							#else
							setVerticesAndGridsWithIntersectionData (j, nb, start, stop, 1, 1, pack, pack_grid);
							#endif
						}
						else
						{
							#if defined(ENABLE_BOOST_THREADS)
							thdGroup.create_thread(boost::bind(&Surface::intersectWithRayBasedAlgorithm, this,
															   j, nb, start, stop, 1, 1, pack, pack_grid));
							#else
							intersectWithRayBasedAlgorithm (j, nb, start, stop, 1, 1, pack, pack_grid);
							#endif
						}

						// SD PB_NEW loading the intersections info
						if (intersectionsInfo != nullptr)
						{
							if (panel == 0)
								intersectionsInfo->push_back(pack_grid);
						}
					}
				}
				else
				{
					int jump = num_threads;
					for (int i=0; i<num_threads; i++)
					{
						packet pack;
						pack.first = &buffersIntersections[i];
						#if !defined(COORD_NORM_PACKING)
						pack.second = &buffersNormals[i];
						#endif

						// SD PB_NEW
						packet pack_grid;
						if (intersectionsInfo != nullptr)
						{
							pack_grid.first = &buffersIntersections_grid[i];
							#if !defined(COORD_NORM_PACKING)
							pack_grid.second = &buffersNormals_grid[i];
							#endif
						}

						if (patchBasedAlgorithm)
						{
							#if defined(ENABLE_BOOST_THREADS)
							thdGroup.create_thread(boost::bind(&Surface::setVerticesAndGridsWithIntersectionData, this,
															   i, nb, i, na, 1, jump, pack, pack_grid));
							#else
							setVerticesAndGridsWithIntersectionData (i, nb, i, na, 1, jump, pack, pack_grid);
							#endif
						}
						else
						{
							#if defined(ENABLE_BOOST_THREADS)
							thdGroup.create_thread(boost::bind(&Surface::intersectWithRayBasedAlgorithm, this,
															   i, nb, i, na, 1, jump, pack, pack_grid));
							#else
							intersectWithRayBasedAlgorithm (i, nb, i, na, 1, jump, pack, pack_grid);
							#endif
						}

						// SD PB_NEW loading the intersections info
						if (intersectionsInfo != nullptr)
						{
							if (panel == 0)
								intersectionsInfo->push_back(pack_grid);
						}
					}
				}
			}
			else
			{
				// each thread deals with a XY ray slab which has a multiple of 4*nb rays
				// status map is a coarse grid of 4^3 voxels, so that conflict-free cross-thread writes
				// are possible if each thread proceeds launching ray slabs which are multiple of 4*nb rays
				int fine_grid_size = 4;
				int jump = num_threads * fine_grid_size;

				for (int j=0; j<num_threads; j++)
				{
					int start = j*fine_grid_size;
					int stop = na;

					packet pack;
					pack.first = &buffersIntersections[j];
					#if !defined(COORD_NORM_PACKING)
					pack.second = &buffersNormals[j];
					#endif

					// SD PB_NEW
					packet pack_grid;
					if (intersectionsInfo != nullptr)
					{
						pack_grid.first = &buffersIntersections_grid[j];
						#if !defined(COORD_NORM_PACKING)
						pack_grid.second = &buffersNormals_grid[j];
						#endif
					}

					if (patchBasedAlgorithm)
					{
						#if defined(ENABLE_BOOST_THREADS)
						thdGroup.create_thread(boost::bind(&Surface::setVerticesAndGridsWithIntersectionData, this,
														   j, nb, start, stop, fine_grid_size, jump, pack, pack_grid));
						#else
						setVerticesAndGridsWithIntersectionData (j, nb, start, stop, fine_grid_size, jump, pack, pack_grid);
						#endif
					}
					else
					{
						#if defined(ENABLE_BOOST_THREADS)
						thdGroup.create_thread(boost::bind(&Surface::intersectWithRayBasedAlgorithm, this,
														   j, nb, start, stop, fine_grid_size, jump, pack, pack_grid));
						#else
						intersectWithRayBasedAlgorithm (j, nb, start, stop, fine_grid_size, jump, pack, pack_grid);
						#endif
					}

					// SD PB_NEW loading the intersections info
					if (intersectionsInfo != nullptr)
					{
						if (panel == 0)
							intersectionsInfo->push_back(pack_grid);
					}
				}
			}

			// end setupnumThreads
			#if defined(ENABLE_BOOST_THREADS)
			// join; final part of the computation of the volume within the surface
			thdGroup.join_all();
			#endif

			// reduce
			for (int j=0; j<num_threads; j++)
			{
				volPanel[panel] += threadPanelVolume[panel][j];
				numFails += threadFailedRays[j];
				numTotalRays += threadTotalRays[j];
			}
			cout << "ok!";
		}

		auto chrono_end = chrono::high_resolution_clock::now();

		chrono::duration<double> ray_tracing_time = chrono_end - chrono_start;
		cout << endl << INFO << "Ray-tracing computation time... ";
		printf ("%.4e [s]", ray_tracing_time.count());

		#if !defined(AVOID_MEM_CHECKS)
		if (!conf.parallelPocketLoop)
		{
			double current_mem_in_MB, peak_mem_in_MB;
			getMemSpace (current_mem_in_MB, peak_mem_in_MB);
			cout << endl << INFO << "Memory required after ray tracing is " << current_mem_in_MB << " MB";
		}
		#endif

		printf("\n%sApproximated %d rays (%.5f %%)", INFO, numFails, 100. * numFails / numTotalRays);
	}
	
	// before bgp identification and vertices storage, cavities are filled if requested
	if (fillCav)
	{
		auto chrono_cav_start = chrono::high_resolution_clock::now();

		cout << endl << INFO << "Performing cavity detection and conditional filling...";
		cout.flush();

		int cav;
		if (!optimizeGrids)
			cav = getCavities();
		else
			cav = getCavitiesWithBilevelStatusMap();

		cout << "ok!";
		cout << endl << INFO << "Detected " << cav << " cavitiy[ies]";
		cout.flush();

		fillCavities(cav_vol);

		if (wellShaped)
		{
			cout << endl << INFO << "Performing cavities shape filtering...";
			cout.flush();

			if (!optimizeGrids)
				filterCavities();
			else
				filterCavitiesWithBilevelStatusMap();
		}
		cout << endl << INFO << "Recovering cavities atoms...";
		cout.flush();

		if (!optimizeGrids)
			getCavitiesAtoms();
		else
			getCavitiesAtomsWithBilevelStatusMap();

		cout << "ok!";

		auto chrono_cav_end = chrono::high_resolution_clock::now();

		chrono::duration<double> chrono_cav_time = chrono_cav_end - chrono_cav_start;
		cout << endl << INFO << "Cavity detection time is ";
		printf ("%.4e [s]", chrono_cav_time.count());

		#if !defined(AVOID_MEM_CHECKS)
		if (!conf.parallelPocketLoop)
		{
			double current_mem_in_MB, peak_mem_in_MB;
			getMemSpace (current_mem_in_MB, peak_mem_in_MB);
			cout << endl << INFO << "Memory required after cavity detection is " << current_mem_in_MB << " MB";
		}
		#endif
	}

	// The following code removes the switched off vertices and assembles the others thanks to assembleVerticesList()
	if (isRCbased)
	{
		auto chrono_assembling_start = chrono::high_resolution_clock::now();

		#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
		cout << endl << INFO << "Assembling octrees...";
		#else
		cout << endl << INFO << "Assembling vertex and normal data...";
		#endif


		int approx_num_vertices = 0;

		for (int j=0; j<num_threads; j++)
		{
			approx_num_vertices += netInts_per_thd[j];
		}
		vertList.reserve(3 * approx_num_vertices);

		if (computeNormals && providesAnalyticalNormals)
		{
			normalsList.reserve(3 * approx_num_vertices);
		}
		deleteVector<int>(netInts_per_thd);

		int vertex_index = 0;

		for (int j=0; j<num_threads; j++)
		{
			packet pack;
			pack.first = &buffersIntersections[j];
			#if !defined(COORD_NORM_PACKING)
			pack.second = &buffersNormals[j];
			#endif

			assembleVerticesList (pack, &vertex_index);

			#if defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			verticesBuffers[j].clear();
			#endif

			if (computeNormals && providesAnalyticalNormals)
			{
				#if defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				normalsBuffers[j].clear();
				#endif
			}
		}

		/*
		// Parallel assembling version; thread safety is not guaranteed in some cases ...

		vector<vector<VERTEX_TYPE*>*> localVert;
		vector<vector<VERTEX_TYPE*>*> localNormals;

		localVert.reserve(num_threads);
		localNormals.reserve(num_threads);

		for (int i=0; i<num_threads; i++)
		{
			localVert.push_back(new vector<VERTEX_TYPE*>());
			localNormals.push_back(new vector<VERTEX_TYPE*>());

			localVert[i]->reserve(max(1000, netInts_per_thd[i]));
			localNormals[i]->reserve(max(1000, netInts_per_thd[i]));
		}
		deleteVector<int>(netInts_per_thd);


		#ifdef ENABLE_BOOST_THREADS
		boost::thread_group thdGroup;
		#endif

		// 2-step stage devoted tossembling final octrees/bilevel grids and vertices/normals matrix
		// It is multithreaded and the second step converts local indices to global ones thanks to
		// index offsets calculated between the two steps
		int *localIndices = allocateVector<int>(num_threads);
		int *indexOffsets = allocateVector<int>(num_threads);

		for (int j=0; j<num_threads; j++)
			localIndices[j] = 0;

		for (int j=0; j<num_threads; j++)
		{
			packet pack;
			pack.first = &buffersIntersections[j];
			#if !defined(COORD_NORM_PACKING)
			pack.second = &buffersNormals[j];
			#endif

			#ifdef ENABLE_BOOST_THREADS
			thdGroup.create_thread(boost::bind(&Surface::assembleVerticesList, this, pack, localVert[j], localNormals[j], &localIndices[j]));
			#else
			assembleVerticesList (pack, localVert[j], localNormals[j], &localIndices[j]);
			#endif
		}
		#if defined(ENABLE_BOOST_THREADS)
		thdGroup.join_all();
		#endif

		// The local indices store in the octrees/bilevel grids are converted to global indices with the following procedure
		indexOffsets[0] = 0;
		for (int j=1; j<num_threads; j++)
			indexOffsets[j] = indexOffsets[j-1] + localIndices[j-1];

		for (int j=0; j<num_threads; j++)
		{
			packet pack;
			pack.first = &buffersIntersections[j];
			#if !defined(COORD_NORM_PACKING)
			pack.second = &buffersNormals[j];
			#endif

			#ifdef ENABLE_BOOST_THREADS
			thdGroup.create_thread(boost::bind(&Surface::convertLocalGridIndicesToGlobalIndices, this, pack, indexOffsets[j]));
			#else
			convertLocalGridIndicesToGlobalIndices (pack, indexOffsets[j]);
			#endif
		}
		#if defined(ENABLE_BOOST_THREADS)
		thdGroup.join_all();
		#endif


		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		int approx_num_vertices = 0;

		for (int j=0; j<num_threads; j++)
		{
			approx_num_vertices += localVert[j]->size();
		}
		vertList.reserve(approx_num_vertices);

		if (computeNormals && providesAnalyticalNormals)
		{
			normalsList.reserve(approx_num_vertices);
		}
		#else
		// Do not allocate a long buffer before using and deleting local buffers (see below)
		vertList.reserve(3 * localVert[0]->size());

		if (computeNormals && providesAnalyticalNormals)
		{
			normalsList.reserve(3 * localNormals[0]->size());
		}
		#endif

		for (int j=0; j<num_threads; j++)
		{
			vector<VERTEX_TYPE*> *lv = localVert[j];
			vector<VERTEX_TYPE*>::iterator it2;

			vector<VERTEX_TYPE*> *ln;
			vector<VERTEX_TYPE*>::iterator it3;

			if (computeNormals && providesAnalyticalNormals)
			{
				ln = localNormals[j];
				it3 = ln->begin();
			}

			for (it2 = lv->begin(); it2 != lv->end(); it2++)
			{
				VERTEX_TYPE *v = *it2;

				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertList.push_back(v);
				#else
				vertList.push_back(v[0]);
				vertList.push_back(v[1]);
				vertList.push_back(v[2]);
				#endif

				if (computeNormals && providesAnalyticalNormals)
				{
					VERTEX_TYPE *vn = *it3;

					#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
					normalsList.push_back(vn);
					#else
					normalsList.push_back(vn[0]);
					normalsList.push_back(vn[1]);
					normalsList.push_back(vn[2]);
					#endif

					it3++;
				}
			}
			#if defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			verticesBuffers[j].clear();
			#endif
			localVert[j]->clear();
			delete localVert[j];

			if (computeNormals && providesAnalyticalNormals)
			{
				#if defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				normalsBuffers[j].clear();
				#endif
				localNormals[j]->clear();
				delete localNormals[j];
			}
		}
		localVert.clear();
		localNormals.clear();

		deleteVector<int>(localIndices);
		deleteVector<int>(indexOffsets);
		*/

		for (int i=0; i<num_threads; i++)
		{
			buffersIntersections[i].clear();
			#if !defined(COORD_NORM_PACKING)
			buffersNormals[i].clear();
			#endif
		}
		delete[] buffersIntersections;
		#if !defined(COORD_NORM_PACKING)
		delete[] buffersNormals;
		#endif

		cout << "ok!";

		auto chrono_assembling_end = chrono::high_resolution_clock::now();

		chrono::duration<double> chrono_assembling_time = chrono_assembling_end - chrono_assembling_start;
		#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
		cout << endl << INFO << "Assembling vertex/normal data time with octrees is ";
		#else
		cout << endl << INFO << "Assembling vertex/normal data time with bilevel grids is ";
		#endif
		printf ("%.4e [s]", chrono_assembling_time.count());

		deleteMatrix2D<double>(3,num_threads,threadPanelVolume);
		deleteVector<int>(threadFailedRays);
		deleteVector<int>(threadTotalRays);

		#if !defined(AVOID_MEM_CHECKS)
		if (!conf.parallelPocketLoop)
		{
			double current_mem_in_MB, peak_mem_in_MB;
			getMemSpace (current_mem_in_MB, peak_mem_in_MB);
			cout << endl << INFO << "Memory required after assembling is " << current_mem_in_MB << " MB";
		}
		#endif
	}

	// bugfix:idebmap has been already fixed by fillcavities, this step would restore
	// the original STATUS_MAP that does not store the final in/out information

	// Stern Layer
	// build it if enabled and if a epsilon map flag was enabled
	// TODO if (sternLayer > 0 && (delphi->buildEpsMap or getDelphiBinding()))
	if (sternLayer > 0 && delphi->buildEpsMap)
	{
		cout << endl << INFO << "Computing Stern Layer...";
		buildSternLayer();
		cout << "ok!";
	}

	// rest bgp count
	delphi->nbgp = 0;


	auto chrono_post_ray_casting_start = chrono::high_resolution_clock::now();

	// custom clean up
	postRayCasting();

	auto chrono_post_ray_casting_end = chrono::high_resolution_clock::now();

	chrono::duration<double> chrono_post_ray_casting_time = chrono_post_ray_casting_end - chrono_post_ray_casting_start;
	cout << endl << INFO << "Post ray tracing time is ";
	printf ("%.4e [s]", chrono_post_ray_casting_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	if (!conf.parallelPocketLoop)
	{
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after post ray tracing is " << current_mem_in_MB << " MB";
	}
	#endif


	// Phase 1.5. If requested correct epsilon map to take into account multi-dielectric
	if (delphi->getMultiDiel())
	{
		if (surfType != MOLECULAR_SURFACE)
		{
			cout << endl << WARN << "Cannot use multi-dielectric in a non molecular surface";
		}
		else if (!delphi->buildEpsMap)
		{
			cout << endl << WARN << "Cannot apply multi-dielectric correction without epsilon map";
		}
		else
		{
			cout << endl << INFO << "Applying multiple dielectric correction...";
			buildAtomsMap();
			applyMultidielectric();
			cout << "ok!";
		}
	}

	int internal_bgps = 0;
	int external_bgps = 0;

	// Phase 2. boundary grid points are located and projected if requested
	if (projBGP)
	{
		vector<int*> bgp;
		vector<int> bgp_type_temp;

		preBoundaryProjection();

		cout << endl << INFO << "Detecting boundary grid points...";

		for (int64_t iz=1; iz<NZ; iz++)
		{
			for (int64_t iy=1; iy<NY; iy++)
			{
				for (int64_t ix=1; ix<NX; ix++)
				{
					// int kost = delphi->EPSMAP(ix,iy,iz,0,NX,NY,NZ);
					const int konst = read4DVector<int>(delphi->epsmap,ix,iy,iz,0,NX,NY,NZ,3);
					/*
					if (konst != delphi->EPSMAP(ix,iy,iz,1,NX,NY,NZ) ||
						konst != delphi->EPSMAP(ix,iy,iz,2,NX,NY,NZ) ||
						konst != delphi->EPSMAP((ix-1),iy,iz,0,NX,NY,NZ) ||
						konst != delphi->EPSMAP(ix,(iy-1),iz,1,NX,NY,NZ) ||
						konst != delphi->EPSMAP(ix,iy,(iz-1),2,NX,NY,NZ) )
					{*/

					if (konst != read4DVector<int>(delphi->epsmap,ix,iy,iz,1,NX,NY,NZ,3) ||
						konst != read4DVector<int>(delphi->epsmap,ix,iy,iz,2,NX,NY,NZ,3) ||
						konst != read4DVector<int>(delphi->epsmap,ix-1,iy,iz,0,NX,NY,NZ,3) ||
						konst != read4DVector<int>(delphi->epsmap,ix,iy-1,iz,1,NX,NY,NZ,3) ||
						konst != read4DVector<int>(delphi->epsmap,ix,iy,iz-1,2,NX,NY,NZ,3) )
					{
						int *temp = allocateVector<int>(3);
						temp[0] = ix;
						temp[1] = iy;
						temp[2] = iz;
						bgp.push_back(temp);

						// if only one point is external than the bgp is external (surface bgp)
						// else it is internal bgp (different dielectrics)
						/*
						if (delphi->EPSMAP(ix,iy,iz,0,NX,NY,NZ) == 0 ||
							delphi->EPSMAP(ix,iy,iz,1,NX,NY,NZ) == 0 ||
							delphi->EPSMAP(ix,iy,iz,2,NX,NY,NZ) == 0 ||
							delphi->EPSMAP((ix-1),iy,iz,0,NX,NY,NZ) == 0 ||
							delphi->EPSMAP(ix,(iy-1),iz,1,NX,NY,NZ) == 0 ||
							delphi->EPSMAP(ix,iy,(iz-1),2,NX,NY,NZ) == 0)
						{*/
						if (
						0 == read4DVector<int>(delphi->epsmap,ix,iy,iz,0,NX,NY,NZ,3) ||
						0 == read4DVector<int>(delphi->epsmap,ix,iy,iz,1,NX,NY,NZ,3) ||
						0 == read4DVector<int>(delphi->epsmap,ix,iy,iz,2,NX,NY,NZ,3) ||
						0 == read4DVector<int>(delphi->epsmap,ix-1,iy,iz,0,NX,NY,NZ,3) ||
						0 == read4DVector<int>(delphi->epsmap,ix,iy-1,iz,1,NX,NY,NZ,3) ||
						0 == read4DVector<int>(delphi->epsmap,ix,iy,iz-1,2,NX,NY,NZ,3) )
						{
							bgp_type_temp.push_back(EXTERNAL_BGP);
							external_bgps++;
						}
						else
						{
							bgp_type_temp.push_back(INTERNAL_BGP);
							internal_bgps++;
						}
					}
				}
			}
		}
		delphi->nbgp = (int)bgp.size();

		if (!delphi->getDelphiBinding())
		{
			delphi->scspos	= allocateVector<double>(3*delphi->nbgp);
			delphi->scsnor	= allocateVector<double>(3*delphi->nbgp);
			delphi->scsarea	= allocateVector<double>(delphi->nbgp);
			delphi->ibgp	= allocateVector<int>(3*delphi->nbgp);
		}

		if (delphi->getDelphiBinding() && bgp.size() >= (unsigned)delphi->maxbgp)
		{
			cout << endl << ERR << "Number of bgp is " << bgp.size() << " and the maximum allowed is " << delphi->maxbgp;
			cout << endl << ERR << "Please increase ibmx in DelPhi and recompile";
			exit(-1);
		}

		if (bgp_type != NULL)
			deleteVector<int>(bgp_type);

		bgp_type = allocateVector<int>(delphi->nbgp);

		// freeze vector
		int i=0;

		for (vector<int*>::iterator it=bgp.begin(); it != bgp.end(); it++)
		{
			//delphi->ibgp[i] = *it;
			int *vv = *it;
			delphi->ibgp[3*i  ] = vv[0];
			delphi->ibgp[3*i+1] = vv[1];
			delphi->ibgp[3*i+2] = vv[2];
			deleteVector<int>(vv);
			i++;
		}

		for (unsigned int l=0; l < bgp_type_temp.size(); l++)
			bgp_type[l] = bgp_type_temp[l];

		cout << "ok!";
		cout << endl << INFO << "Detected " << delphi->nbgp << " boundary grid points (bgp)";
		cout << endl << INFO << "Detected " << external_bgps << " external bgps, and " << internal_bgps << " internal bgps";
		cout << endl << INFO << "Scaling bgps...";


		#if !defined(MULTITHREADED_POCKET_LOOP)

		#ifdef ENABLE_BOOST_THREADS
		// the overhead of thread creation is quite low, so it is
		// convenient to dispatch as many threads as possible so as
		// to saturate all cores and to minimize their unbalancing
		num_threads = MIN(64, MAX(1, delphi->nbgp));
		// cout << endl << INFO << "Automatically dispatching " << num_threads << " threads";
		#endif

		#ifdef ENABLE_BOOST_THREADS
		boost::thread_group thdGroup;
		#endif

		int chunk = delphi->nbgp/num_threads;
		int rem = delphi->nbgp%num_threads;

		// setup split
		int start, stop;

		for (int j=0; j<num_threads; j++)
		{
			if (j == 0)
			{
				start = 0;
				stop = chunk;
			}
			else
			{
				start = stop;
				stop = start+chunk;
			}

			if (j < rem)
				stop++;

			#ifdef ENABLE_BOOST_THREADS
			thdGroup.create_thread(boost::bind(&Surface::projector, this, start, stop));
			#else
			projector(start, stop);
			#endif
		}
		// end setup
		
		#ifdef ENABLE_BOOST_THREADS
		thdGroup.join_all();
		#endif

		#else // MULTITHREADED_POCKET_LOOP

		projector(0, delphi->nbgp);

		#endif // MULTITHREADED_POCKET_LOOP

		cout << "ok!";

		deleteVector<int>(bgp_type);
		bgp_type = NULL;

		int *atsurf = delphi->atsurf;

		// setting atsurf = NULL it is a way to avoid this computation
		if (delphi->getDelphiBinding() && atsurf != NULL)
		{
			cout << endl << INFO << "Linking boundary grid points to nearest atom";
			// if multi-diel is enabled the atoms map is already available
			// if not we have to build it
			if (!delphi->getMultiDiel())
				buildAtomsMap();

			double *l_scspos = delphi->scspos;

			for (int i=0; i<delphi->nbgp; i++)
			{
				int nearestAtom = -1;
				double *v = &l_scspos[3*i];
				vdwAccessible(v, nearestAtom);

				if (nearestAtom == -1)
					cout << endl << WARN << "Cannot detect nearest atom for bgp index " << i;

				atsurf[i] = nearestAtom;
			}
		}
	}
	
	// if multi-dielectric or stern layer was computed then dispose atoms map
	// in any case if delphi is bound and bgps have been computed we have to drop it
	if (delphi->getMultiDiel() || (delphi->getDelphiBinding() && delphi->atsurf != NULL && projBGP))
		disposeAtomsMap();

	int flag_sum;
	flag_sum  = panelVolumeFlag[0][0] + panelVolumeFlag[1][0] + panelVolumeFlag[2][0];
	flag_sum += panelVolumeFlag[0][1] + panelVolumeFlag[1][1] + panelVolumeFlag[2][1];

	*surf_volume = delphi->A * (volPanel[0]+volPanel[1]+volPanel[2]) / flag_sum;

	// if (delphi->buildStatus)
	// 	*surf_volume = delphi->A * volPanel[0];
	// else
	// 	*surf_volume = delphi->A * (volPanel[0]+volPanel[1]+volPanel[2]) / 3.;

	return true;
}


void Surface::buildAtomsMap()
{
	// Build a 3D accelaration grid for the atoms
	double rmax = 0;
	// get the biggest atom 
	// for (int i=0; i<delphi->numAtoms; i++)
	for (int i=0; i<delphi->atoms.size(); i++)
		rmax = MAX((delphi->atoms[i].radius),rmax);
	
	gscale = 0.5/rmax;
	ggrid = 1;

	while (1)
	{
		ggrid = (unsigned int)floor(gscale*1.3*delphi->rmaxdim);

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
	gside = 1./gscale;
	gxmin = delphi->baricenter[0] - (ggrid-1)*0.5*gside;
	gymin = delphi->baricenter[1] - (ggrid-1)*0.5*gside;
	gzmin = delphi->baricenter[2] - (ggrid-1)*0.5*gside;

	if (gridMultiMap != NULL)
	{
		for (int64_t i=0; i < ggrid*ggrid*ggrid; i++)
		{
			gridMultiMap[i].clear();
		}
		delete[] gridMultiMap;
		gridMultiMap = NULL;
	}
	gridMultiMap = new vector<int> [ ggrid*ggrid*ggrid ];

	// fill 3D acceleration grid with atoms
	// for (int i=0; i<delphi->numAtoms; i++)
	for (int i=0; i<delphi->atoms.size(); i++)
	{
		int64_t ix = (int64_t)rintp((delphi->atoms[i].pos[0] - gxmin)*gscale);
		int64_t iy = (int64_t)rintp((delphi->atoms[i].pos[1] - gymin)*gscale);
		int64_t iz = (int64_t)rintp((delphi->atoms[i].pos[2] - gzmin)*gscale);

		gridMultiMap[ iz*(ggrid*ggrid) + iy*ggrid + ix ].push_back(i);
	}
}


void Surface::disposeAtomsMap()
{
	// remove acceleration grid
	if (gridMultiMap != NULL)
	{
		for (int64_t i=0; i < ggrid*ggrid*ggrid; i++)
		{
			gridMultiMap[i].clear();
		}
		delete[] gridMultiMap;
		gridMultiMap = NULL;
	}
}


void Surface::applyMultidielectric()
{
	// for each internal point in the eps map gets the nearest atom
	// and fix the epsilon map accordingly.

	if (gridMultiMap == NULL)
	{
		cout << endl << WARN << "Cannot apply multi-dielectric correction without atoms map";
		return;
	}

	// for each internal point get nearest atom and change eps value accordingly
	for (int64_t k=0; k<delphi->nz; k++)
	{
		for (int64_t j=0; j<delphi->ny; j++)
		{
			for (int64_t i=0; i<delphi->nx; i++)
			{
				// for each cube face fix the epsmap value
				//int value = delphi->EPSMAP(i,j,k,0,delphi->nx,delphi->ny,delphi->nz);
				int value = read4DVector<int>(delphi->epsmap,i,j,k,0,delphi->nx,delphi->ny,delphi->nz,3);

				if (value == inside)
					swap2multi(gxmin,gymin,gzmin,gside,ggrid,gridMultiMap,i,j,k,0);

				//value = delphi->EPSMAPint *(i,j,k,1,(delphi->nx),(delphi->ny),(delphi->nz));
				value = read4DVector<int>(delphi->epsmap,i,j,k,1,delphi->nx,delphi->ny,delphi->nz,3);

				if (value == inside)
					swap2multi(gxmin,gymin,gzmin,gside,ggrid,gridMultiMap,i,j,k,1);
				
				//value = delphi->EPSMAP(i,j,k,2,(delphi->nx),(delphi->ny),(delphi->nz));
				value = read4DVector<int>(delphi->epsmap,i,j,k,2,delphi->nx,delphi->ny,delphi->nz,3);
				
				if (value == inside)
					swap2multi(gxmin,gymin,gzmin,gside,ggrid,gridMultiMap,i,j,k,2);
			}
		}
	}
}


void Surface::swap2multi(double gxmin, double gymin, double gzmin, double gside, unsigned int ggrid,
						 vector<int> *gridMultiMap, int i, int j, int k, int l)
{
	// generate point position
	double pos[3];

	pos[0] = delphi->x[i];
	pos[1] = delphi->y[j];
	pos[2] = delphi->z[k];

	int ori=i, orj=j, ork=k, orl=l;
	// position the point on the corrseponding cube side. Pos now is the position of the midpoint
	pos[l] += delphi->hside;

	// move from midpoint of delphi grid to auxiliary atoms grid
	int ix = (int)rintp((pos[0] - gxmin)*gscale);
	int iy = (int)rintp((pos[1] - gymin)*gscale);
	int iz = (int)rintp((pos[2] - gzmin)*gscale);

	double minDist = INFINITY;
	int winner = -1;
	
	// get the nearest atom to set the dielectric constant
	// the dielectric constant is mapped according to the additively weighted voronoi diagram
	// that is the signed distance from the point p is ||p-c||^2-r^2 where c is the center
	// of the atom and r is the radius. The minimum signed distance wins.
	for (k=0; k<SHIFT_MAP; k++)
	{
		int64_t cx = ix+shift_map[k][0];
		int64_t cy = iy+shift_map[k][1];
		int64_t cz = iz+shift_map[k][2];

		if (cx>=ggrid || cy>=ggrid || cz>=ggrid || cx<0 || cy<0 || cz<0)
			continue;

		int num_atoms = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ].size();

		for (int j=0; j<num_atoms; j++)
		{
			int atom_index = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ][j];

			double signed_dist = 0;
			DIST2(signed_dist,delphi->atoms[atom_index].pos,pos)

			double rad = delphi->atoms[atom_index].radius;
			signed_dist -= rad*rad;

			if (signed_dist < minDist)
			{
				minDist = signed_dist;
				winner = atom_index;
			}
		}
	}
	
	if (winner == -1)
	{
		cout << endl << WARN << "No winner atom found!";
		return;
	}

	// store winner
	// delphi->EPSMAP(ori,orj,ork,orl,delphi->nx,delphi->ny,delphi->nz) = delphi->atoms[winner]->dielectric;
	write4DVector<int>(delphi->epsmap,delphi->atoms[winner].dielectric,ori,orj,ork,orl,delphi->nx,delphi->ny,delphi->nz,3);
	return;
}


int Surface::getCavities(int idStart)
{
	// search cavities in status vector
	int id = idStart;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int *status = delphi->status;


	if (status == NULL)
	{
		cout << endl << ERR << "Cannot do cavity detection without a status map";
		cout << endl << REMARK << "Please set Build_status_map = true";
		cout << endl;
		exit(-1);
	}

	// free memory 
	if (delphi->cavitiesVec != NULL)
	{
		for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin();
			  it != delphi->cavitiesVec->end(); it++)
		{
			vector<int*> *inner = (*it);
			for (vector<int*>::iterator it2 = inner->begin(); it2 != inner->end(); it2++)
				deleteVector<int>(*it2);
			delete inner;
		}
		delete delphi->cavitiesVec;
	}
	
	// clean memory
	delphi->cavitiesSize.clear();
	delphi->cavitiesFlag.clear();

	// allocate empty vector
	delphi->cavitiesVec = new vector<vector<int*>*>();
	delphi->cavitiesVec->reserve(50);
	
	int times = 0;
	int startZ = 0;
	int endZ = NZ-1;
	int stepZ = 1;
	int sig = 1;
	
	while (1)
	{
		int64_t i,j,k;
		bool cavities = false;
		
		if (times%2 == 0)
		{
			startZ = 0;
			endZ = NZ-1;
			stepZ = 1;
			sig = 1;
		}
		else
		{
			startZ = NZ-1;
			endZ = 0;
			stepZ = -1;
			sig = -1;
		}
		
		times++;
		
		// get starting point
		for (k = startZ; sig*k <= endZ; k += stepZ)
		{
			for (j=0; j<NY; j++)
			{
				for (i=0; i<NX; i++)
				{
					// if (STATUSMAP(i,j,k,NX,NY) == STATUS_POINT_TEMPORARY_OUT)
					if (read3DVector<int>(status,i,j,k,NX,NY,NZ) == STATUS_POINT_TEMPORARY_OUT)
					{
						cavities = true;
						break;
					}
				}
				if (cavities)
					break;
			}
			if (cavities)
				break;
		}

		if (!cavities)
			break;

		// new status
		id++;

		// mark from temporary outside to cavity index.
		// if idStart = STATUS_POINT_TEMPORARY_OUT in the first pass STATUS_POINT_TEMPORARY_OUT -> STATUS_POINT_OUT
		// from that moment on, if there are still cavities they are all marked with STATUS_POINT_TEMPORARY_OUT
		// and they are detected one by one and associated to a cavity index
		// floodFill(i,j,k,STATUS_POINT_TEMPORARY_OUT,id);
		// floodFill4(i,j,k,STATUS_POINT_TEMPORARY_OUT,id);
		floodFill2(i,j,k,STATUS_POINT_TEMPORARY_OUT,id);
	}

	int numCavities = id - STATUS_POINT_OUT;
	delphi->cavitiesSize.resize(numCavities);
	delphi->cavitiesFlag.resize(numCavities);

	for (int i=0; i<numCavities; i++)
		delphi->cavitiesFlag[i] = false;

	return numCavities;
}


int Surface::getCavitiesWithBilevelStatusMap(int idStart)
{
	// search cavities in status vector
	int id = idStart;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;
	int64_t coarse_nx = getCoarseN(NX);
	int64_t coarse_ny = getCoarseN(NY);

	int **bilevel_status = delphi->bilevel_status;


	if (bilevel_status == NULL)
	{
		cout << endl << ERR << "Cannot do cavity detection without a bilevel status map";
		cout << endl << REMARK << "Please set Build_status_map = true";
		cout << endl;
		exit(-1);
	}

	// free memory
	if (delphi->cavitiesVec != NULL)
	{
		for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin();
			  it != delphi->cavitiesVec->end(); it++)
		{
			vector<int*> *inner = (*it);
			for (vector<int*>::iterator it2 = inner->begin(); it2 != inner->end(); it2++)
				deleteVector<int>(*it2);
			delete inner;
		}
		delete delphi->cavitiesVec;
	}

	// clean memory
	delphi->cavitiesSize.clear();
	delphi->cavitiesFlag.clear();

	// allocate empty vector
	delphi->cavitiesVec = new vector<vector<int*>*>();
	delphi->cavitiesVec->reserve(50);

	int times = 0;
	int startZ = 0;
	int endZ = NZ-1;
	int stepZ = 1;
	int sig = 1;

	while (1)
	{
		int64_t i,j,k;
		bool cavities = false;

		if (times%2 == 0)
		{
			startZ = 0;
			endZ = NZ-1;
			stepZ = 1;
			sig = 1;
		}
		else
		{
			startZ = NZ-1;
			endZ = 0;
			stepZ = -1;
			sig = -1;
		}

		times++;

		// get starting point
		for (k = startZ; sig*k <= endZ; k += stepZ)
		{
			int64_t coarse_k = getCoarseID(NZ, k);
			int64_t fine_k = getFineID(coarse_k, k);

			for (j=0; j<NY; j++)
			{
				int64_t coarse_j = getCoarseID(NY, j);
				int64_t fine_j = getFineID(coarse_j, j);

				for (i=0; i<NX; i++)
				{
					int64_t coarse_i = getCoarseID(NX, i);
					int64_t coarse_index = coarse_k*(coarse_ny*coarse_nx) + coarse_j*coarse_nx + coarse_i;

					if (bilevel_status[coarse_index] == NULL)
					{
						cavities = true;
						break;
					}
					int64_t fine_i = getFineID(coarse_i, i);
					int64_t fine_index = getUnrolledFineID(fine_i, fine_j, fine_k);

					if (bilevel_status[coarse_index][fine_index] == STATUS_POINT_TEMPORARY_OUT)
					{
						cavities = true;
						break;
					}
				}
				if (cavities)
					break;
			}
			if (cavities)
				break;
		}

		if (!cavities)
			break;

		// new status
		id++;

		// mark from temporary outside to cavity index.
		// if idStart = STATUS_POINT_TEMPORARY_OUT in the first pass STATUS_POINT_TEMPORARY_OUT -> STATUS_POINT_OUT
		// from that moment on, if there are still cavities they are all marked with STATUS_POINT_TEMPORARY_OUT
		// and they are detected one by one and associated to a cavity index
		// floodFillWithBilevelStatusMap(i,j,k,STATUS_POINT_TEMPORARY_OUT,id);
		// floodFillWithBilevelStatusMap4(i,j,k,STATUS_POINT_TEMPORARY_OUT,id);
		floodFillWithBilevelStatusMap2(i,j,k,STATUS_POINT_TEMPORARY_OUT,id);
	}

	int numCavities = id - STATUS_POINT_OUT;
	delphi->cavitiesSize.resize(numCavities);
	delphi->cavitiesFlag.resize(numCavities);

	for (int i=0; i<numCavities; i++)
		delphi->cavitiesFlag[i] = false;

	return numCavities;
}


void Surface::getCavitiesAtoms()
{
	delphi->initCav2Atoms();
	vector<set<int>*> &cav2atoms = delphi->cav2atoms;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int *status = delphi->status;


	if (delphi->buildStatus == false)
	{
		cout << endl << WARN << "Cannot get cavity atoms without a status map";
		return;
	}

	if (delphi->cavitiesVec->size() == 0)
		return;

	buildAtomsMap();
	int i=0;

	// cout << endl << cav2atoms.size() << " " << delphi->cavitiesVec->size();

	for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin();
		 it != delphi->cavitiesVec->end(); it++)
	{
		vector<int*> *vec = (*it);

		// Get bgps for each cavity. And obtain the nearest atom
		// in order to avoid using epsmap a slightly changed notion of bgp is employed
		for (vector<int*>::iterator it2 = vec->begin(); it2 != vec->end(); it2++)
		{
			int *v = (*it2);
			int64_t ix = v[0], iy = v[1], iz = v[2];

			// cout << endl << ix << " " << iy << " " << iz;

			// int kost = STATUSMAP(ix,iy,iz,NX,NY);
			const int kost = read3DVector<int>(status,ix,iy,iz,NX,NY,NZ);

			// is it pseudo-bgp? if yes continue else skip.
			/*
			if (kost != STATUSMAP((ix+1),iy,iz,NX,NY) ||
				kost != STATUSMAP((ix-1),iy,iz,NX,NY) ||
				kost != STATUSMAP(ix,(iy+1),iz,NX,NY) ||
				kost != STATUSMAP(ix,(iy-1),iz,NX,NY) ||
				kost != STATUSMAP(ix,iy,(iz+1),NX,NY) ||
				kost != STATUSMAP(ix,iy,(iz-1),NX,NY))
			*/
			if (kost != read3DVector<int>(status,ix+1,iy,iz,NX,NY,NZ) ||
				kost != read3DVector<int>(status,ix-1,iy,iz,NX,NY,NZ) ||
				kost != read3DVector<int>(status,ix,iy+1,iz,NX,NY,NZ) ||
				kost != read3DVector<int>(status,ix,iy-1,iz,NX,NY,NZ) ||
				kost != read3DVector<int>(status,ix,iy,iz+1,NX,NY,NZ) ||
				kost != read3DVector<int>(status,ix,iy,iz-1,NX,NY,NZ))
			{
				;
			}
			else
			{
				continue;
			}

			double gridPoint[3];

			gridPoint[0] = delphi->x[ix];
			gridPoint[1] = delphi->y[iy];
			gridPoint[2] = delphi->z[iz];

			// move from grid point of delphi grid to auxiliary atoms grid
			ix = (int64_t)rintp((gridPoint[0]-gxmin)*gscale);
			iy = (int64_t)rintp((gridPoint[1]-gymin)*gscale);
			iz = (int64_t)rintp((gridPoint[2]-gzmin)*gscale);

			double minDist = INFINITY;
			int first = -1;

			// get the nearest atom
			for (int k=0; k<SHIFT_MAP; k++)
			{
				int64_t cx = ix+shift_map[k][0];
				int64_t cy = iy+shift_map[k][1];
				int64_t cz = iz+shift_map[k][2];

				if (cx >= ggrid || cy >= ggrid || cz >= ggrid || cx < 0 || cy < 0 || cz < 0)
					continue;

				int num_atoms = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ].size();

				for (int j=0; j<num_atoms; j++)
				{
					int atom_index = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ][j];

					double signed_dist = 0;

					DIST2(signed_dist,delphi->atoms[atom_index].pos,gridPoint)
					double radius2 = delphi->atoms[atom_index].radius2;
					signed_dist -= radius2;

					if (signed_dist < minDist)
					{
						minDist = signed_dist;
						first = atom_index;
					}
				}
			}
			if (first == -1)
				cout << endl << WARN << "No nearest atom in cavity/pocket!";

			cav2atoms[i]->insert(first);
		}
		i++;
	}
	disposeAtomsMap();
}


void Surface::getCavitiesAtomsWithBilevelStatusMap()
{
	delphi->initCav2Atoms();
	vector<set<int>*> &cav2atoms = delphi->cav2atoms;
	
	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int **bilevel_status = delphi->bilevel_status;


	if (delphi->buildStatus == false)
	{
		cout << endl << WARN << "Cannot get cavity atoms without a bilevel status map";
		return;
	}

	if (delphi->cavitiesVec->size() == 0)
		return;
	
	buildAtomsMap();
	int i=0;

	// cout << endl << cav2atoms.size() << " " << delphi->cavitiesVec->size();

	for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin();
		  it != delphi->cavitiesVec->end(); it++)
	{	
		vector<int*> *vec = (*it);
		
		// Get bgps for each cavity. And obtain the nearest atom
		// in order to avoid using epsmap a slightly changed notion of bgp is employed
		for (vector<int*>::iterator it2 = vec->begin(); it2 != vec->end(); it2++)
		{
			int *v = (*it2);
			int64_t ix = v[0], iy = v[1], iz = v[2];

			// cout << endl << ix << " " << iy << " " << iz;
			
			const int kost = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz,NX,NY,NZ);
			
			if (kost == readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix+1,iy,iz,NX,NY,NZ) &&
				kost == readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix-1,iy,iz,NX,NY,NZ) &&
				kost == readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy+1,iz,NX,NY,NZ) &&
				kost == readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy-1,iz,NX,NY,NZ) &&
				kost == readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz+1,NX,NY,NZ) &&
				kost == readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz-1,NX,NY,NZ))
			{
				continue;
			}

			double gridPoint[3];

			gridPoint[0] = delphi->x[ix];
			gridPoint[1] = delphi->y[iy];
			gridPoint[2] = delphi->z[iz];

			// move from grid point of delphi grid to auxiliary atoms grid
			ix = (int64_t)rintp((gridPoint[0]-gxmin)*gscale);
			iy = (int64_t)rintp((gridPoint[1]-gymin)*gscale);
			iz = (int64_t)rintp((gridPoint[2]-gzmin)*gscale);
		
			double minDist = INFINITY;
			int first = -1;
			
			// get the nearest atom
			for (int k=0; k<SHIFT_MAP; k++)
			{
				int64_t cx = ix+shift_map[k][0];
				int64_t cy = iy+shift_map[k][1];
				int64_t cz = iz+shift_map[k][2];

				if (cx >= ggrid || cy >= ggrid || cz >= ggrid || cx < 0 || cy < 0 || cz < 0)
					continue;

				int num_atoms = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ].size();

				for (int j=0; j<num_atoms; j++)
				{
					int atom_index = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ][j];

					double signed_dist = 0;

					DIST2(signed_dist,delphi->atoms[atom_index].pos,gridPoint)
					double radius2 = delphi->atoms[atom_index].radius2;
					signed_dist -= radius2;
					
					if (signed_dist < minDist)
					{
						minDist = signed_dist;
						first = atom_index;
					}
				}
			}
			if (first == -1)
				cout << endl << WARN << "No nearest atom in cavity/pocket!";

			cav2atoms[i]->insert(first);
		}
		i++;
	}
	disposeAtomsMap();
}


void Surface::fillCavities(double vol, bool silent)
{
	if (vol < 0)
	{
		cout << endl << WARN << "Cannot fill with a negative volume. Setting " << DEFAULT_VOLUME;
		vol = DEFAULT_VOLUME;
	}
	// analyze each cavity and fill if required
	int i = 0;
	double cubeVol = delphi->side*delphi->side*delphi->side;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	if (!silent)
	{
		cout << endl << INFO << "Threshold volume is " << vol;
		cout << endl << INFO << "Tot num cavities is " << delphi->cavitiesVec->size();
	}
	for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++)
	{
		if (!silent)
			cout << endl << INFO << "Cavity " << i;
		
		double cavVol = (*it)->size()*cubeVol;
		if (!silent)
			printf("\tvol is %.4lf [A^3] \t", cavVol);

		if (internals != NULL)
		{
			*internals << endl << "cav " << cavVol;
		}

		delphi->cavitiesSize[i] = cavVol;
		delphi->cavitiesFlag[i] = false;

		// we don't fill cavities on STATUS_MAP because we want
		// to retain the full memory of cavities on this structure

		if (cavVol <= vol)
		{
			delphi->cavitiesFlag[i] = true;
			vector<int*> *vec = (*it);

			// filling eps map
			for (vector<int*>::iterator it2 = vec->begin(); it2 != vec->end(); it2++)
			{
				int *v = (*it2);
				// that grid point is filled

				// apply correction only if epsmap is used
				if (delphi->buildEpsMap)
				{
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2],0,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2],1,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2],2,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0]-1,v[1],v[2],0,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1]-1,v[2],1,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2]-1,2,NX,NY,NZ,3);

					//delphi->IDEBMAP((v[0]),(v[1]),(v[2]),NX,NY) = false;
					write3DVector<bool>(delphi->idebmap,false,v[0],v[1],v[2],NX,NY,NZ);
				}

				if (accurateTriangulation && !isAvailableScalarField)
				{
					#if !defined(USE_COMPRESSED_GRIDS)
					if (!optimizeGrids)
					{
						verticesInsidenessMap[v[2]][v[1]][v[0]] = false;
						verticesInsidenessMap[v[2]][v[1]][v[0]+1] = false;
						verticesInsidenessMap[v[2]][v[1]+1][v[0]] = false;
						verticesInsidenessMap[v[2]][v[1]+1][v[0]+1] = false;

						verticesInsidenessMap[v[2]+1][v[1]][v[0]] = false;
						verticesInsidenessMap[v[2]+1][v[1]][v[0]+1] = false;
						verticesInsidenessMap[v[2]+1][v[1]+1][v[0]] = false;
						verticesInsidenessMap[v[2]+1][v[1]+1][v[0]+1] = false;
						/*
						write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2],NX,NY,NZ);
						write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2],NX,NY,NZ);
						write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2],NX,NY,NZ);
						write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2],NX,NY,NZ);

						write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2]+1,NX,NY,NZ);
						write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2]+1,NX,NY,NZ);
						write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2]+1,NX,NY,NZ);
						write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2]+1,NX,NY,NZ);
						*/
					}
					else
					#endif
					{
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1],v[2],NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1],v[2],NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1]+1,v[2],NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2],NX,NY,NZ);

						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1],v[2]+1,NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1],v[2]+1,NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1]+1,v[2]+1,NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2]+1,NX,NY,NZ);
					}
				}
				// update the scalar filed consistently with cavity detection
				if (accurateTriangulation && isAvailableScalarField)
				{
					scalarField[v[2]][v[1]][v[0]] = +INFINITY;
					scalarField[v[2]][v[1]][v[0]+1] = +INFINITY;
					scalarField[v[2]][v[1]+1][v[0]] = +INFINITY;
					scalarField[v[2]][v[1]+1][v[0]+1] = +INFINITY;

					scalarField[v[2]+1][v[1]][v[0]] = +INFINITY;
					scalarField[v[2]+1][v[1]][v[0]+1] = +INFINITY;
					scalarField[v[2]+1][v[1]+1][v[0]] = +INFINITY;
					scalarField[v[2]+1][v[1]+1][v[0]+1] = +INFINITY;
				}
			}
			if (!silent)
				cout << "filled";
		}
		else
		{
			if (!silent)
				cout << "non filled";
		}
		i++;
	}
}


/*
// conservative filter

void Surface::filterCavities()
{	
	int ngrid = (int)rintp(2.*probe_radius*delphi->scale);
	if ((ngrid%2) == 0)
		ngrid++;

	int k = (ngrid-1)/2;
	// cout << " Probe radius is " << probe_radius << " [A], grid points " << k;
	int cavityId;

	double dist2ref = probe_radius*probe_radius;
	int i = 0;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	// analyze each cavity
	for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin();
		  it != delphi->cavitiesVec->end(); it++, i++)
	{
		cout << endl << INFO << "Cavity " << i ;
		if (delphi->cavitiesFlag[i])
		{
			cout << " already filled";
			continue;
		}
		else
			cout << " Under analysis...";

		vector<int*> *currentCavity = (*it);

		// mark each cavity point as temporary STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK code
		for (vector<int*>::iterator it2 = currentCavity->begin();
			  it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;
			cavityId = delphi->status[v[0]][v[1]][v[2]];
			delphi->status[v[0]][v[1]][v[2]] = STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK;
		}
		// for each cavity point 'march' around it in order to understand if the 
		// water molecule can fit in the cavity
		
		// a cavity point is filled if a probe sphere cannot stay into the cavity.

		bool keepCavity = false;

		for (vector<int*>::iterator it2 = currentCavity->begin();
			  it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;
			int start_ix = v[0]-k;
			int start_iy = v[1]-k;
			int start_iz = v[2]-k;

			int end_ix = v[0]+k;
			int end_iy = v[1]+k;
			int end_iz = v[2]+k;
		
			bool oneProbeFits = true;

			for (int iz = start_iz; iz<end_iz; iz++)
				for (int iy = start_iy; iy<end_iy; iy++)
					for (int ix = start_ix; ix<end_ix; ix++)
					{
						double a = (delphi->x[ix]-delphi->x[v[0]])*(delphi->x[ix]-delphi->x[v[0]]);
						double b = (delphi->y[iy]-delphi->y[v[1]])*(delphi->y[iy]-delphi->y[v[1]]);
						double c = (delphi->z[iz]-delphi->z[v[2]])*(delphi->z[iz]-delphi->z[v[2]]);

						double dist2 = a+b+c;
						if (dist2 > dist2ref)
							continue;

						// there is no space we went outside the cavity, go on to the next grid point
						if (delphi->status[iz][iy][ix] != STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK)
						{
							oneProbeFits = false;
							break;
						}
					}
					if (!oneProbeFits)
						break;
				}
				if (!oneProbeFits)
					break;
			}

			// water molecule fits, stop, keep cavity
			if (oneProbeFits)
			{
				keepCavity = true;
				break;
			}
		}

		for (vector<int*>::iterator it2 = currentCavity->begin();
			  it2 != currentCavity->end(); it2++)
		{	
			int *v = *it2;
			// clean cavity
			if (!keepCavity)
			{
				if (delphi->buildEpsMap)
				{
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),0,NX,NY,NZ) = inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),1,NX,NY,NZ) = inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),2,NX,NY,NZ) = inside;
					delphi->EPSMAP((v[0]-1),(v[1]),(v[2]),0,NX,NY,NZ) = inside;
					delphi->EPSMAP((v[0]),(v[1]-1),(v[2]),1,NX,NY,NZ) = inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]-1),2,NX,NY,NZ) = inside;
					
					delphi->IDEBMAP((v[0]),(v[1]),(v[2]),NX,NY) = false;
				}

				if (accurateTriangulation && !isAvailableScalarField)
				{
					verticesInsidenessMap[v[2]][v[1]][v[0]] = false;
					verticesInsidenessMap[v[2]][v[1]][v[0]+1] = false;
					verticesInsidenessMap[v[2]][v[1]+1][v[0]] = false;
					verticesInsidenessMap[v[2]][v[1]+1][v[0]+1] = false;

					verticesInsidenessMap[v[2]+1][v[1]][v[0]] = false;
					verticesInsidenessMap[v[2]+1][v[1]][v[0]+1] = false;
					verticesInsidenessMap[v[2]+1][v[1]+1][v[0]] = false;
					verticesInsidenessMap[v[2]+1][v[1]+1][v[0]+1] = false;
				}
				
				// update the scalar filed consistently with cavity shape detection
				if (accurateTriangulation && isAvailableScalarField)
				{
					scalarField[v[2]][v[1]][v[0]] = +INFINITY;
					scalarField[v[2]][v[1]][v[0]+1] = +INFINITY;
					scalarField[v[2]][v[1]+1][v[0]] = +INFINITY;
					scalarField[v[2]][v[1]+1][v[0]+1] = +INFINITY;

					scalarField[v[2]+1][v[1]][v[0]] = +INFINITY;
					scalarField[v[2]+1][v[1]][v[0]+1] = +INFINITY;
					scalarField[v[2]+1][v[1]+1][v[0]] = +INFINITY;
					scalarField[v[2]+1][v[1]+1][v[0]+1] = +INFINITY;
				}
			}
			
			// restore orginal value in any case
			delphi->status[v[0]][v[1]][v[2]] = cavityId;
		} // end for each cavity point

		if (!keepCavity)
		{
			delphi->cavitiesFlag[i] = true;
			cout << "Filled";
		}
		else
			cout << "Non filled";

	} // end for each cavity
}
*/


void Surface::cav2out()
{
	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int *status = delphi->status;

	int i = 0;

	// for each cavity if it is not filled revert to OUT
	// if it is filled revert to IN

	if (delphi->cavitiesVec == NULL)
	{
		cout << endl << "WARN" << "Cannot convert cavities to out flags if cavities are absent!";
		return;
	}

	for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++, i++)
	{
		vector<int*> *currentCavity = (*it);

		for (vector<int*>::iterator it2 = currentCavity->begin(); it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;

			// if a confirmed cavity (it can be STATUS_POINT_INISDE if filtered) then switch to out
			// if it is a point that is a support of the cavity, switch it too
			/*
			if (STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) >= STATUS_FIRST_CAV ||
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) <= STATUS_FIRST_SUPPORT_CAV)
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = STATUS_POINT_TEMPORARY_OUT;
			*/
			int val = read3DVector<int>(status,v[0],v[1],v[2],NX,NY,NZ);
			if (val >= STATUS_FIRST_CAV || val <= STATUS_FIRST_SUPPORT_CAV)
				write3DVector<int>(status,STATUS_POINT_TEMPORARY_OUT,v[0],v[1],v[2],NX,NY,NZ);
		}
	}
}


void Surface::cav2outWithBilevelStatusMap()
{
	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;
	int64_t coarse_nx = getCoarseN(NX);
	int64_t coarse_ny = getCoarseN(NY);

	int **bilevel_status = delphi->bilevel_status;

	int i = 0;

	// for each cavity if it is not filled revert to OUT
	// if it is filled revert to IN

	if (delphi->cavitiesVec == NULL)
	{
		cout << endl << "WARN" << "Cannot convert cavities to out flags if cavities are absent!";
		return;
	}

	for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++, i++)
	{
		vector<int*> *currentCavity = (*it);

		for (vector<int*>::iterator it2 = currentCavity->begin(); it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;

			int64_t coarse_i = getCoarseID(NX, v[0]);
			int64_t coarse_j = getCoarseID(NY, v[1]);
			int64_t coarse_k = getCoarseID(NZ, v[2]);
			int64_t fine_i = getFineID(coarse_i, v[0]);
			int64_t fine_j = getFineID(coarse_j, v[1]);
			int64_t fine_k = getFineID(coarse_k, v[2]);
			int64_t coarse_index = coarse_k*(coarse_ny*coarse_nx) + coarse_j*coarse_nx + coarse_i;
			int64_t fine_index = getUnrolledFineID(fine_i, fine_j, fine_k);

			// if a confirmed cavity (it can be STATUS_POINT_INISDE if filtered) then switch to out
			// if it is a point that is a support of the cavity, switch it too
			int val;

			if (bilevel_status[coarse_index] == NULL)
				val = STATUS_POINT_TEMPORARY_OUT;
			else
				val = bilevel_status[coarse_index][fine_index];

			if (val >= STATUS_FIRST_CAV || val <= STATUS_FIRST_SUPPORT_CAV)
			{
				if (bilevel_status[coarse_index] == NULL)
				{
					bilevel_status[coarse_index] = allocateBilevelMinigridCells<int>(STATUS_POINT_TEMPORARY_OUT);
				}
				bilevel_status[coarse_index][fine_index] = STATUS_POINT_TEMPORARY_OUT;
			}
		}
	}
}


void Surface::filterCavities(bool modStatus)
{
	// cout << " Probe radius to filter cavities is " << probe_radius << " [A]";
	int cavityId;
	double dist2ref = probe_radius*probe_radius;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int *status = delphi->status;

	int i = 0;


	// Connolly like filter
	for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++, i++)
	{
		if (delphi->cavitiesFlag[i])
			continue;

		vector<int*> *currentCavity = (*it);
		cavityId = i + STATUS_FIRST_CAV;

		// mark each cavity point as temporary STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK code
		for (vector<int*>::iterator it2 = currentCavity->begin(); it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;

			write3DVector<int>(status,STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK,v[0],v[1],v[2],NX,NY,NZ);
		}

		// for each cavity point 'march' around it in order to understand if the
		// water molecule can fit in around the current cavity point
		for (vector<int*>::iterator it2 = currentCavity->begin(); it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;

			double downx = delphi->x[v[0]] - probe_radius;
			double downy = delphi->y[v[1]] - probe_radius;
			double downz = delphi->z[v[2]] - probe_radius;

			double upx = delphi->x[v[0]] + probe_radius;
			double upy = delphi->y[v[1]] + probe_radius;
			double upz = delphi->z[v[2]] + probe_radius;

			// Determine which are the grid cells that
			// are occupied by the bounding box of the object
			int64_t ix_start = (int64_t)rintp((downx-delphi->xmin)*delphi->scale);
			int64_t iy_start = (int64_t)rintp((downy-delphi->ymin)*delphi->scale);
			int64_t iz_start = (int64_t)rintp((downz-delphi->zmin)*delphi->scale);

			int64_t ix_end = (int64_t)rintp((upx-delphi->xmin)*delphi->scale);
			int64_t iy_end = (int64_t)rintp((upy-delphi->ymin)*delphi->scale);
			int64_t iz_end = (int64_t)rintp((upz-delphi->zmin)*delphi->scale);

			if (ix_start < 0) ix_start = 0;
			if (iy_start < 0) iy_start = 0;
			if (iz_start < 0) iz_start = 0;

			if (ix_end < 0) ix_end = 0;
			if (iy_end < 0) iy_end = 0;
			if (iz_end < 0) iz_end = 0;

			if (ix_end >= NX) ix_end = NX-1;
			if (iy_end >= NY) iy_end = NY-1;
			if (iz_end >= NZ) iz_end = NZ-1;

			if (ix_start >= NX) ix_start = NX-1;
			if (iy_start >= NY) iy_start = NY-1;
			if (iz_start >= NZ) iz_start = NZ-1;

			bool go_out = false;

			for (int64_t iz = iz_start; iz<=iz_end; iz++)
			{
				for (int64_t iy = iy_start; iy<=iy_end; iy++)
				{
					for (int64_t ix = ix_start; ix<=ix_end; ix++)
					{
						double a = (delphi->x[ix] - delphi->x[v[0]]);
						double b = (delphi->y[iy] - delphi->y[v[1]]);
						double c = (delphi->z[iz] - delphi->z[v[2]]);

						double dist2 = a*a + b*b + c*c;
						if (dist2 > dist2ref)
							continue;

						int val = read3DVector<int>(status,ix,iy,iz,NX,NY,NZ);

						// this is inside the probe and it is an out point, stop immediately, go to next cavity point
						/*
						 *					if (STATUSMAP(ix,iy,iz,NX,NY) != STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK &&
						 *						STATUSMAP(ix,iy,iz,NX,NY) != STATUS_CAVITY_POINT_SHAPE_OK &&
						 *						STATUSMAP(ix,iy,iz,NX,NY) != -cavityId)
						 */

						if (val != STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK &&
							val != STATUS_CAVITY_POINT_SHAPE_OK &&
							val != -cavityId)
						{
							go_out = true;
							break;
						}
					}
					if (go_out)
						break;
				}
				if (go_out)
					break;
			}

			// Record that a probe can stay completely inside
			// assign all points inside accordingly
			if (!go_out)
			{
				for (int64_t iz = iz_start; iz <= iz_end; iz++)
				{
					for (int64_t iy = iy_start; iy <= iy_end; iy++)
					{
						for (int64_t ix = ix_start; ix <= ix_end; ix++)
						{
							double a = (delphi->x[ix]-delphi->x[v[0]]);
							double b = (delphi->y[iy]-delphi->y[v[1]]);
							double c = (delphi->z[iz]-delphi->z[v[2]]);

							double dist2 = a*a + b*b + c*c;
							if (dist2 > dist2ref)
								continue;

							// don't overwrite support cavity points
							if (read3DVector<int>(status,ix,iy,iz,NX,NY,NZ) != -cavityId)
								write3DVector<int>(status,STATUS_CAVITY_POINT_SHAPE_OK,ix,iy,iz,NX,NY,NZ);
						}
					}
				}
				// the point is marked as a support of the cavity
				write3DVector<int>(status,-cavityId,v[0],v[1],v[2],NX,NY,NZ);
			}
		}

		// number of grid points removed
		int numRemoved = 0;

		// the parts of the cavity that remained unmarked are due to the
		// non possibility to fit on it the probe
		// thus these points are remarked as inside.
		for (vector<int*>::iterator it2 = currentCavity->begin(); it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;

			// here the probe does not fit thus put inside. this point is no more in cavity
			if (read3DVector<int>(status,v[0],v[1],v[2],NX,NY,NZ) == STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK)
			{
				numRemoved++;

				if (delphi->buildEpsMap)
				{
					/*
					 *				delphi->EPSMAP((v[0]),(v[1]),(v[2]),0,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]),(v[1]),(v[2]),1,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]),(v[1]),(v[2]),2,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]-1),(v[1]),(v[2]),0,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]),(v[1]-1),(v[2]),1,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]),(v[1]),(v[2]-1),2,NX,NY,NZ) = inside;
					 *
					 *				delphi->IDEBMAP((v[0]),(v[1]),(v[2]),NX,NY) = false;
					 */
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2],0,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2],1,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2],2,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0]-1,v[1],v[2],0,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1]-1,v[2],1,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2]-1,2,NX,NY,NZ,3);

					write3DVector<bool>(delphi->idebmap,false,v[0],v[1],v[2],NX,NY,NZ);

				}
				if (accurateTriangulation && !isAvailableScalarField)
				{
					#if !defined(USE_COMPRESSED_GRIDS)
					if (!optimizeGrids)
					{
						verticesInsidenessMap[v[2]][v[1]][v[0]] = false;
						verticesInsidenessMap[v[2]][v[1]][v[0]+1] = false;
						verticesInsidenessMap[v[2]][v[1]+1][v[0]] = false;
						verticesInsidenessMap[v[2]][v[1]+1][v[0]+1] = false;

						verticesInsidenessMap[v[2]+1][v[1]][v[0]] = false;
						verticesInsidenessMap[v[2]+1][v[1]][v[0]+1] = false;
						verticesInsidenessMap[v[2]+1][v[1]+1][v[0]] = false;
						verticesInsidenessMap[v[2]+1][v[1]+1][v[0]+1] = false;
						/*
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2],NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2],NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2],NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2],NX,NY,NZ);
						 *
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2]+1,NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2]+1,NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2]+1,NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2]+1,NX,NY,NZ);
						 */
					}
					else
						#endif
					{
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1],v[2],NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1],v[2],NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1]+1,v[2],NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2],NX,NY,NZ);

						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1],v[2]+1,NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1],v[2]+1,NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1]+1,v[2]+1,NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2]+1,NX,NY,NZ);
					}
				}

				// update the scalar filed consistently with cavity shape detection
				if (accurateTriangulation && isAvailableScalarField)
				{
					scalarField[v[2]][v[1]][v[0]] = +INFINITY;
					scalarField[v[2]][v[1]][v[0]+1] = +INFINITY;
					scalarField[v[2]][v[1]+1][v[0]] = +INFINITY;
					scalarField[v[2]][v[1]+1][v[0]+1] = +INFINITY;

					scalarField[v[2]+1][v[1]][v[0]] = +INFINITY;
					scalarField[v[2]+1][v[1]][v[0]+1] = +INFINITY;
					scalarField[v[2]+1][v[1]+1][v[0]] = +INFINITY;
					scalarField[v[2]+1][v[1]+1][v[0]+1] = +INFINITY;
				}
			}

			// if this point must be inside and status is allowed to be modified then set it inside
			if (modStatus && read3DVector<int>(status,v[0],v[1],v[2],NX,NY,NZ) == STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK)
				write3DVector<int>(status,STATUS_POINT_INSIDE,v[0],v[1],v[2],NX,NY,NZ);
			// else restore the orginal value
			else if (read3DVector<int>(status,v[0],v[1],v[2],NX,NY,NZ) == STATUS_CAVITY_POINT_SHAPE_OK)
				write3DVector<int>(status,cavityId,v[0],v[1],v[2],NX,NY,NZ);

		} // end for each cavity point

		if (numRemoved && currentCavity->size() != numRemoved)
		{
			// thus they can get filtered again
			delphi->cavitiesFlag[i] = false;
			// cout << "Partially filled";
		}
		else if (numRemoved && currentCavity->size() == numRemoved)
		{
			delphi->cavitiesFlag[i] = true;
			// cout << "Completely filled";
		}
		//else
		//	cout << "Non filled";
	} // end for each cavity
}


void Surface::filterCavitiesWithBilevelStatusMap(bool modStatus)
{
	// cout << " Probe radius to filter cavities is " << probe_radius << " [A]";
	int cavityId;
	double dist2ref = probe_radius*probe_radius;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int64_t coarse_nx, coarse_ny;
	int64_t coarse_i, coarse_j, coarse_k;
	int64_t fine_i, fine_j, fine_k;

	coarse_nx = getCoarseN(NX);
	coarse_ny = getCoarseN(NY);

	int **bilevel_status = delphi->bilevel_status;

	int i = 0;


	// Connolly like filter
	for (vector<vector<int*>*>::iterator it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++, i++)
	{
		if (delphi->cavitiesFlag[i])
			continue;

		vector<int*> *currentCavity = (*it);
		cavityId = i + STATUS_FIRST_CAV;

		// mark each cavity point as temporary STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK code
		for (vector<int*>::iterator it2 = currentCavity->begin(); it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;

			writeBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK,v[0],v[1],v[2],NX,NY,NZ);
		}

		// for each cavity point 'march' around it in order to understand if the
		// water molecule can fit in around the current cavity point
		for (vector<int*>::iterator it2 = currentCavity->begin(); it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;

			double downx = delphi->x[v[0]] - probe_radius;
			double downy = delphi->y[v[1]] - probe_radius;
			double downz = delphi->z[v[2]] - probe_radius;

			double upx = delphi->x[v[0]] + probe_radius;
			double upy = delphi->y[v[1]] + probe_radius;
			double upz = delphi->z[v[2]] + probe_radius;

			// Determine which are the grid cells that
			// are occupied by the bounding box of the object
			int64_t ix_start = (int64_t)rintp((downx-delphi->xmin)*delphi->scale);
			int64_t iy_start = (int64_t)rintp((downy-delphi->ymin)*delphi->scale);
			int64_t iz_start = (int64_t)rintp((downz-delphi->zmin)*delphi->scale);

			int64_t ix_end = (int64_t)rintp((upx-delphi->xmin)*delphi->scale);
			int64_t iy_end = (int64_t)rintp((upy-delphi->ymin)*delphi->scale);
			int64_t iz_end = (int64_t)rintp((upz-delphi->zmin)*delphi->scale);

			if (ix_start < 0) ix_start = 0;
			if (iy_start < 0) iy_start = 0;
			if (iz_start < 0) iz_start = 0;

			if (ix_end < 0) ix_end = 0;
			if (iy_end < 0) iy_end = 0;
			if (iz_end < 0) iz_end = 0;

			if (ix_end >= NX) ix_end = NX-1;
			if (iy_end >= NY) iy_end = NY-1;
			if (iz_end >= NZ) iz_end = NZ-1;

			if (ix_start >= NX) ix_start = NX-1;
			if (iy_start >= NY) iy_start = NY-1;
			if (iz_start >= NZ) iz_start = NZ-1;

			bool go_out = false;

			for (int64_t iz = iz_start; iz<=iz_end; iz++)
			{
				coarse_k = getCoarseID(NZ, iz);
				fine_k = getFineID(coarse_k, iz);

				for (int64_t iy = iy_start; iy<=iy_end; iy++)
				{
					coarse_j = getCoarseID(NY, iy);
					fine_j = getFineID(coarse_j, iy);

					for (int64_t ix = ix_start; ix<=ix_end; ix++)
					{
						coarse_i = getCoarseID(NX, ix);
						fine_i = getFineID(coarse_i, ix);

						double a = (delphi->x[ix] - delphi->x[v[0]]);
						double b = (delphi->y[iy] - delphi->y[v[1]]);
						double c = (delphi->z[iz] - delphi->z[v[2]]);

						double dist2 = a*a + b*b + c*c;
						if (dist2 > dist2ref)
							continue;

						int val;

						int64_t coarse_index = coarse_k*(coarse_ny*coarse_nx) + coarse_j*coarse_nx + coarse_i;

						if (bilevel_status[coarse_index] == NULL)
							val = STATUS_POINT_TEMPORARY_OUT;
						else
						{
							int64_t fine_index = getUnrolledFineID(fine_i, fine_j, fine_k);
							val = bilevel_status[coarse_index][fine_index];
						}
						if (val != STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK &&
							val != STATUS_CAVITY_POINT_SHAPE_OK &&
							val != -cavityId)
						{
							go_out = true;
							break;
						}
					}
					if (go_out)
						break;
				}
				if (go_out)
					break;
			}

			// Record that a probe can stay completely inside
			// assign all points inside accordingly
			if (!go_out)
			{
				for (int64_t iz = iz_start; iz <= iz_end; iz++)
				{
					coarse_k = getCoarseID(NZ, iz);
					fine_k = getFineID(coarse_k, iz);

					for (int64_t iy = iy_start; iy <= iy_end; iy++)
					{
						coarse_j = getCoarseID(NY, iy);
						fine_j = getFineID(coarse_j, iy);

						for (int64_t ix = ix_start; ix <= ix_end; ix++)
						{
							double a = (delphi->x[ix]-delphi->x[v[0]]);
							double b = (delphi->y[iy]-delphi->y[v[1]]);
							double c = (delphi->z[iz]-delphi->z[v[2]]);

							double dist2 = a*a + b*b + c*c;
							if (dist2 > dist2ref)
								continue;

							// don't overwrite support cavity points
							coarse_i = getCoarseID(NX, ix);
							fine_i = getFineID(coarse_i, ix);

							int64_t coarse_index = coarse_k*(coarse_ny*coarse_nx) + coarse_j*coarse_nx + coarse_i;
							int64_t fine_index = getUnrolledFineID(fine_i, fine_j, fine_k);

							int val = STATUS_POINT_TEMPORARY_OUT;

							if (bilevel_status[coarse_index] != NULL)
								val = bilevel_status[coarse_index][fine_index];

							if (val != -cavityId)
							{
								if (bilevel_status[coarse_index] == NULL)
								{
									bilevel_status[coarse_index] = allocateBilevelMinigridCells<int>(STATUS_POINT_TEMPORARY_OUT);
								}
								bilevel_status[coarse_index][fine_index] = STATUS_CAVITY_POINT_SHAPE_OK;
							}
						}
					}
				}
				writeBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,-cavityId,v[0],v[1],v[2],NX,NY,NZ);
			}
		}

		// number of grid points removed
		int numRemoved = 0;

		// the parts of the cavity that remained unmarked are due to the
		// non possibility to fit on it the probe
		// thus these points are remarked as inside.
		for (vector<int*>::iterator it2 = currentCavity->begin(); it2 != currentCavity->end(); it2++)
		{
			int *v = *it2;

			int64_t coarse_i = getCoarseID(NX, v[0]);
			int64_t coarse_j = getCoarseID(NY, v[1]);
			int64_t coarse_k = getCoarseID(NZ, v[2]);
			int64_t fine_i = getFineID(coarse_i, v[0]);
			int64_t fine_j = getFineID(coarse_j, v[1]);
			int64_t fine_k = getFineID(coarse_k, v[2]);
			int64_t coarse_index = coarse_k*(coarse_ny*coarse_nx) + coarse_j*coarse_nx + coarse_i;
			int64_t fine_index = getUnrolledFineID(fine_i, fine_j, fine_k);

			int val;

			if (bilevel_status[coarse_index] == NULL)
				val = STATUS_POINT_TEMPORARY_OUT;
			else
				val = bilevel_status[coarse_index][fine_index];

			// here the probe does not fit thus put inside. this point is no more in cavity
			if (val == STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK)
			{
				numRemoved++;

				if (delphi->buildEpsMap)
				{
					/*
					 *				delphi->EPSMAP((v[0]),(v[1]),(v[2]),0,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]),(v[1]),(v[2]),1,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]),(v[1]),(v[2]),2,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]-1),(v[1]),(v[2]),0,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]),(v[1]-1),(v[2]),1,NX,NY,NZ) = inside;
					 *				delphi->EPSMAP((v[0]),(v[1]),(v[2]-1),2,NX,NY,NZ) = inside;
					 *
					 *				delphi->IDEBMAP((v[0]),(v[1]),(v[2]),NX,NY) = false;
					 */
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2],0,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2],1,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2],2,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0]-1,v[1],v[2],0,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1]-1,v[2],1,NX,NY,NZ,3);
					write4DVector<int>(delphi->epsmap,inside,v[0],v[1],v[2]-1,2,NX,NY,NZ,3);

					write3DVector<bool>(delphi->idebmap,false,v[0],v[1],v[2],NX,NY,NZ);

				}
				if (accurateTriangulation && !isAvailableScalarField)
				{
					#if !defined(USE_COMPRESSED_GRIDS)
					if (!optimizeGrids)
					{
						verticesInsidenessMap[v[2]][v[1]][v[0]] = false;
						verticesInsidenessMap[v[2]][v[1]][v[0]+1] = false;
						verticesInsidenessMap[v[2]][v[1]+1][v[0]] = false;
						verticesInsidenessMap[v[2]][v[1]+1][v[0]+1] = false;

						verticesInsidenessMap[v[2]+1][v[1]][v[0]] = false;
						verticesInsidenessMap[v[2]+1][v[1]][v[0]+1] = false;
						verticesInsidenessMap[v[2]+1][v[1]+1][v[0]] = false;
						verticesInsidenessMap[v[2]+1][v[1]+1][v[0]+1] = false;
						/*
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2],NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2],NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2],NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2],NX,NY,NZ);
						 *
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2]+1,NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2]+1,NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2]+1,NX,NY,NZ);
						 *					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2]+1,NX,NY,NZ);
						 */
					}
					else
						#endif
					{
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1],v[2],NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1],v[2],NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1]+1,v[2],NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2],NX,NY,NZ);

						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1],v[2]+1,NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1],v[2]+1,NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0],v[1]+1,v[2]+1,NX,NY,NZ);
						write32xCompressedGrid(compressed_verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2]+1,NX,NY,NZ);
					}
				}

				// update the scalar filed consistently with cavity shape detection
				if (accurateTriangulation && isAvailableScalarField)
				{
					scalarField[v[2]][v[1]][v[0]] = +INFINITY;
					scalarField[v[2]][v[1]][v[0]+1] = +INFINITY;
					scalarField[v[2]][v[1]+1][v[0]] = +INFINITY;
					scalarField[v[2]][v[1]+1][v[0]+1] = +INFINITY;

					scalarField[v[2]+1][v[1]][v[0]] = +INFINITY;
					scalarField[v[2]+1][v[1]][v[0]+1] = +INFINITY;
					scalarField[v[2]+1][v[1]+1][v[0]] = +INFINITY;
					scalarField[v[2]+1][v[1]+1][v[0]+1] = +INFINITY;
				}
			}

			int val_to_write = STATUS_POINT_TEMPORARY_OUT;

			// if this point must be inside and status is allowed to be modified then set it inside
			if (modStatus && val == STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK)
				val_to_write = STATUS_POINT_INSIDE;
			// else restore the orginal value
			else if (val == STATUS_CAVITY_POINT_SHAPE_OK)
				val_to_write = cavityId;

			if (val_to_write != STATUS_POINT_TEMPORARY_OUT)
			{
				if (bilevel_status[coarse_index] == NULL)
				{
					bilevel_status[coarse_index] = allocateBilevelMinigridCells<int>(STATUS_POINT_TEMPORARY_OUT);
				}
				bilevel_status[coarse_index][fine_index] = val_to_write;
			}
		} // end for each cavity point

		if (numRemoved && currentCavity->size() != numRemoved)
		{
			// thus they can get filtered again
			delphi->cavitiesFlag[i] = false;
			// cout << "Partially filled";
		}
		else if (numRemoved && currentCavity->size() == numRemoved)
		{
			delphi->cavitiesFlag[i] = true;
			// cout << "Completely filled";
		}
		//else
		//	cout << "Non filled";
	} // end for each cavity
}


/** queue based implementation*/
void Surface::floodFill(int ix, int iy, int iz, int idold, int idnew)
{
	queue<pair<int,pair<int,int>>> myqueue;
	pair<int,pair<int,int>> triplet;
	int cix,ciy,ciz;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int *status = delphi->status;
	int **bilevel_status = delphi->bilevel_status;


	myqueue.push( pair<int,pair<int,int>> (ix,pair<int,int>(iy,iz)) );

	while (myqueue.size() != 0)
	{
		// cout << endl << myqueue.size();
		triplet = myqueue.front();
		myqueue.pop();

		cix = triplet.first;
		ciy = triplet.second.first;
		ciz = triplet.second.second;

		// cout << endl << myqueue.size();
		// cout << endl <<cix << " " << ciy << " " << ciz;

		if (!optimizeGrids && read3DVector<int>(status,cix,ciy,ciz,NX,NY,NZ) != idold)
			continue;
		else if (optimizeGrids && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix,ciy,ciz,NX,NY,NZ) != idold)
			continue;
		else
		{
			if (!optimizeGrids)
				write3DVector<int>(status,idnew,cix,ciy,ciz,NX,NY,NZ);
			else
				writeBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,idnew,cix,ciy,ciz,NX,NY,NZ);

			// if this is a cavity then update the cavity list
			if (idnew >= 4)
			{
				// current cavity vector
				vector<int*> *vec;

				// first time this cavity is encountered
				if (delphi->cavitiesVec->size() < idnew - 4 + 1)
				{
					//cout << endl << INFO << "Detected cavity " << idnew-4;
					vec = new vector<int*>();
					if (vec == NULL)
					{
						cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
						cout.flush();
						exit(-1);
					}
					delphi->cavitiesVec->push_back(vec);
				}
				else
				{
					vec = delphi->cavitiesVec->at(idnew - 4);
				}
				int *v = allocateVector<int>(3);
				if (v == NULL)
				{
					cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
					cout.flush();
					exit(-1);
				}
				v[0] = cix;
				v[1] = ciy;
				v[2] = ciz;
				// printf("\n %d %d %d (%d,%d)",cix,ciy,ciz,idnew,idold);
				vec->push_back(v);

			}
		}

		if ((!optimizeGrids && (cix+1)<NX && read3DVector<int>(status,cix+1,ciy,ciz,NX,NY,NZ) == idold) ||
			( optimizeGrids && (cix+1)<NX && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix+1,ciy,ciz,NX,NY,NZ) == idold))
		{
			triplet.first = cix+1;
			triplet.second.first = ciy;
			triplet.second.second = ciz;
			myqueue.push(triplet);
		}
		if ((!optimizeGrids && (ciy+1)<NY && read3DVector<int>(status,cix,ciy+1,ciz,NX,NY,NZ) == idold) ||
			( optimizeGrids && (ciy+1)<NY && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix,ciy+1,ciz,NX,NY,NZ) == idold))
		{
			triplet.first = cix;
			triplet.second.first = ciy+1;
			triplet.second.second = ciz;
			myqueue.push(triplet);
		}
		// if (ciz+1<delphi->nz && (delphi->status[ciz+1][ciy][cix]) == idold)
		// if (ciz+1<NZ && (STATUSMAP(cix,ciy,(ciz+1),NX,NY)) == idold)
		if ((!optimizeGrids && (ciz+1)<NZ && read3DVector<int>(status,cix,ciy,ciz+1,NX,NY,NZ) == idold) ||
			( optimizeGrids && (ciz+1)<NZ && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix,ciy,ciz+1,NX,NY,NZ) == idold))
		{
			triplet.first = cix;
			triplet.second.first = ciy;
			triplet.second.second = ciz+1;
			myqueue.push(triplet);
		}
		// if (cix-1>=0 && (delphi->status[ciz][ciy][cix-1]) == idold)
		// if (cix-1>=0 && (STATUSMAP((cix-1),(ciy),(ciz),NX,NY)) == idold)
		if ((!optimizeGrids && (cix-1)>=0 && read3DVector<int>(status,cix-1,ciy,ciz,NX,NY,NZ) == idold) ||
			( optimizeGrids && (cix-1)>=0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix-1,ciy,ciz,NX,NY,NZ) == idold))
		{
			triplet.first = cix-1;
			triplet.second.first = ciy;
			triplet.second.second = ciz;
			myqueue.push(triplet);
		}
		// if (ciy-1>=0 && (delphi->status[ciz][ciy-1][cix]) == idold)
		// if (ciy-1>=0 && (STATUSMAP(cix,(ciy-1),ciz,NX,NY)) == idold)
		if ((!optimizeGrids && (ciy-1)>=0 && read3DVector<int>(status,cix,ciy-1,ciz,NX,NY,NZ) == idold) ||
			( optimizeGrids && (ciy-1)>=0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix,ciy-1,ciz,NX,NY,NZ) == idold))
		{
			triplet.first = cix;
			triplet.second.first = ciy-1;
			triplet.second.second = ciz;
			myqueue.push(triplet);
		}
		// if (ciz-1>=0 && (delphi->status[ciz-1][ciy][cix]) == idold)
		// if (ciz-1>=0 && (STATUSMAP(cix,ciy,(ciz-1),NX,NY)) == idold)
		if ((!optimizeGrids && (ciz-1)>=0 && read3DVector<int>(status,cix,ciy,ciz-1,NX,NY,NZ) == idold) ||
			( optimizeGrids && (ciz-1)>=0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix,ciy,ciz-1,NX,NY,NZ) == idold))
		{
			triplet.first = cix;
			triplet.second.first = ciy;
			triplet.second.second = ciz-1;
			myqueue.push(triplet);
		}
	}
	return;
}


void Surface::floodFill4(int ix, int iy, int iz, int idold, int idnew)
{
	int num_threads = conf.numThreads;

	int maxMoves = 10;
	vector<pair<pair<int,int>,int>> seeds_down;
	vector<pair<pair<int,int>,int>> seeds_up;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int *status = delphi->status;
	int **bilevel_status = delphi->bilevel_status;


	seeds_down.push_back( pair<pair<int,int>,int> (pair<int,int>(ix,iy),iz) );

	int ix_or = ix;
	int iy_or = iy;
	int iz_or = iz;

	cout << endl << INFO << "Z-percolation...";
	cout.flush();

	/////////////////////////////////
	// z-percolation to init threads
	/////////////////////////////////

	// lower percolation
	while (1)
	{
		// if (iz-1 >= 0 && STATUSMAP(ix,iy,iz-1,NX,NY) == idold)
		if ((!optimizeGrids && iz-1 >= 0 && read3DVector<int>(status,ix,iy,iz-1,NX,NY,NZ) == idold) ||
			( optimizeGrids && iz-1 >= 0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz-1,NX,NY,NZ) == idold))
		{
			seeds_down.push_back( pair<pair<int,int>,int> (pair<int,int> (ix,iy),iz-1) );
			iz--;
			continue;
		}

		bool found = false;

		// +x test
		for (int i=0; i<maxMoves; i++)
		{
			// feasibility of the move
			// if (!(ix+i < NX) || !(STATUSMAP(ix+i,iy,iz,NX,NY) == idold))
			if ((!optimizeGrids && (!(ix+i < NX) || !(read3DVector<int>(status,ix+1,iy,iz,NX,NY,NZ) == idold))) ||
				( optimizeGrids && (!(ix+i < NX) || !(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix+1,iy,iz,NX,NY,NZ) == idold))))
			{
				found = false;
				break;
			}
			// if (ix+i < NX && iz-1 >= 0 && STATUSMAP(ix+i,iy,iz-1,NX,NY) == idold)
			if ((!optimizeGrids && ix+i < NX && iz-1 >= 0 && read3DVector<int>(status,ix+1,iy,iz-1,NX,NY,NZ) == idold) ||
				( optimizeGrids && ix+i < NX && iz-1 >= 0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix+1,iy,iz-1,NX,NY,NZ) == idold))
			{
				seeds_down.push_back( pair<pair<int,int>,int> (pair<int,int> (ix+i,iy),iz-1) );
				found = true;
				ix += i;
				iz--;
				break;
			}
		}

		if (found)
			continue;

		// -x test
		for (int i=0; i<maxMoves; i++)
		{
			// feasibility of the move
			// if (!(ix-i >= 0) || !(STATUSMAP(ix-i,iy,iz,NX,NY) == idold))
			if ((!optimizeGrids && (!(ix-i >= 0) || !(read3DVector<int>(status,ix-i,iy,iz,NX,NY,NZ) == idold))) ||
				( optimizeGrids && (!(ix-i >= 0) || !(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix-i,iy,iz,NX,NY,NZ) == idold))))
			{
				found = false;
				break;
			}
			// if (ix-i >= 0 && iz-1 >= 0 && STATUSMAP(ix-i,iy,iz-1,NX,NY) == idold)
			if ((!optimizeGrids && ix-i >= 0 && iz-1 >= 0 && read3DVector<int>(status,ix-1,iy,iz-1,NX,NY,NZ) == idold) ||
				( optimizeGrids && ix-i >= 0 && iz-1 >= 0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix-1,iy,iz-1,NX,NY,NZ) == idold))
			{
				seeds_down.push_back( pair<pair<int,int>,int> (pair<int,int> (ix-i,iy),iz-1) );
				found = true;
				ix -= i;
				iz--;
				break;
			}
		}

		// +y test
		for (int i=0; i<maxMoves; i++)
		{
			// feasibility of the move
			// if (!(iy+i < NY) || !(STATUSMAP(ix,iy+i,iz,NX,NY) == idold))
			if ((!optimizeGrids && (!(iy+i < NY) || !(read3DVector<int>(status,ix,iy+i,iz,NX,NY,NZ) == idold))) ||
				( optimizeGrids && (!(iy+i < NY) || !(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy+i,iz,NX,NY,NZ) == idold))))
			{
				found = false;
				break;
			}
			// if (iy+i < NY && iz-1 >= 0 && STATUSMAP(ix,iy+i,iz-1,NX,NY) == idold)
			if ((!optimizeGrids && iy+i < NY && iz-1 >= 0 && read3DVector<int>(status,ix,iy+i,iz-1,NX,NY,NZ) == idold) ||
				( optimizeGrids && iy+i < NY && iz-1 >= 0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy+i,iz-1,NX,NY,NZ) == idold))
			{
				seeds_down.push_back( pair<pair<int,int>,int> (pair<int,int> (ix,iy+i),iz-1) );
				found = true;
				iy += i;
				iz--;
				break;
			}
		}

		if (found)
			continue;

		// -y test
		for (int i=0; i<maxMoves; i++)
		{
			// feasibility of the move
			// if (!(iy-i >= 0) || !(STATUSMAP(ix,iy-i,iz,NX,NY) == idold))
			if ((!optimizeGrids && (!(iy-i >= 0) || !(read3DVector<int>(status,ix,iy-i,iz,NX,NY,NZ) == idold))) ||
				( optimizeGrids && (!(iy-i >= 0) || !(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy-i,iz,NX,NY,NZ) == idold))))
			{
				found = false;
				break;
			}
			// if (iy-i >= 0 && iz-1 >= 0 && STATUSMAP(ix,iy-i,iz-1,NX,NY) == idold)
			if ((!optimizeGrids && iy-i >= 0 && iz-1 >= 0 && read3DVector<int>(status,ix,iy-i,iz-1,NX,NY,NZ) == idold) ||
				( optimizeGrids && iy-i >= 0 && iz-1 >= 0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy-i,iz-1,NX,NY,NZ) == idold))
			{
				seeds_down.push_back( pair<pair<int,int>,int> (pair<int,int> (ix,iy-i),iz-1) );
				found = true;
				iy -= i;
				iz--;
				break;
			}
		}

		if (!found)
			break;
	}

	ix = ix_or;
	iy = iy_or;
	iz = iz_or;

	// upper percolation
	while (1)
	{
		// if (iz+1 < NZ && STATUSMAP(ix,iy,iz+1,NX,NY) == idold)
		if ((!optimizeGrids && iz+1 < NZ && read3DVector<int>(status,ix,iy,iz+1,NX,NY,NZ) == idold) ||
			( optimizeGrids && iz+1 < NZ && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz+1,NX,NY,NZ) == idold))
		{
			seeds_up.push_back( pair<pair<int,int>,int> (pair<int,int> (ix,iy),iz+1) );
			iz++;
			continue;
		}

		bool found = false;

		// +x test
		for (int i=0; i<maxMoves; i++)
		{
			// feasibility of the move
			// if (!(ix+i < NX) || !(STATUSMAP(ix+i,iy,iz,NX,NY) == idold))
			if ((!optimizeGrids && (!(ix+i < NX) || !(read3DVector<int>(status,ix+i,iy,iz,NX,NY,NZ) == idold))) ||
				( optimizeGrids && (!(ix+i < NX) || !(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix+i,iy,iz,NX,NY,NZ) == idold))))
			{
				found = false;
				break;
			}
			// if (ix+i < NX && iz+1 < NZ && STATUSMAP(ix+i,iy,iz+1,NX,NY) == idold)
			if ((!optimizeGrids && ix+i < NX && iz+1 < NZ && read3DVector<int>(status,ix+i,iy,iz+1,NX,NY,NZ) == idold) ||
				( optimizeGrids && ix+i < NX && iz+1 < NZ && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix+i,iy,iz+1,NX,NY,NZ) == idold))
			{
				seeds_up.push_back( pair<pair<int,int>,int> (pair<int,int> (ix+i,iy),iz+1) );
				found = true;
				ix += i;
				iz++;
				break;
			}
		}

		if (found)
			continue;

		// -x test
		for (int i=0; i<maxMoves; i++)
		{
			// feasibility of the move
			// if (!(ix-i >= 0) || !(STATUSMAP(ix-i,iy,iz,NX,NY) == idold))
			if ((!optimizeGrids && (!(ix-i >= 0) || !(read3DVector<int>(status,ix-i,iy,iz,NX,NY,NZ) == idold))) ||
				( optimizeGrids && (!(ix-i >= 0) || !(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix-i,iy,iz,NX,NY,NZ) == idold))))
			{
				found = false;
				break;
			}
			// if (ix-i >= 0 && iz+1 < NZ && STATUSMAP(ix-i,iy,iz+1,NX,NY) == idold)
			if ((!optimizeGrids && ix-i >= 0 && iz+1 < NZ && read3DVector<int>(status,ix-i,iy,iz+1,NX,NY,NZ) == idold) ||
				( optimizeGrids && ix-i >= 0 && iz+1 < NZ && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix-i,iy,iz+1,NX,NY,NZ) == idold))
			{
				seeds_up.push_back( pair<pair<int,int>,int> (pair<int,int> (ix-i,iy),iz+1) );
				found = true;
				ix -= i;
				iz++;
				break;
			}
		}

		// +y test
		for (int i=0; i<maxMoves; i++)
		{
			// feasibility of the move
			// if (!(iy+i < NY) || !(STATUSMAP(ix,iy+i,iz,NX,NY) == idold))
			if ((!optimizeGrids && (!(iy+i < NY) || !(read3DVector<int>(status,ix,iy+i,iz,NX,NY,NZ) == idold))) ||
				( optimizeGrids && (!(iy+i < NY) || !(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy+i,iz,NX,NY,NZ) == idold))))
			{
				found = false;
				break;
			}
			// if (iy+i < NY && iz+1 < NZ && STATUSMAP(ix,iy+i,iz+1,NX,NY) == idold)
			if ((!optimizeGrids && iy+i < NY && iz+1 < NZ && read3DVector<int>(status,ix,iy+i,iz+1,NX,NY,NZ) == idold) ||
				( optimizeGrids && iy+i < NY && iz+1 < NZ && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy+i,iz+1,NX,NY,NZ) == idold))
			{
				seeds_up.push_back( pair<pair<int,int>,int> (pair<int,int> (ix,iy+i),iz+1) );
				found = true;
				iy += i;
				iz++;
				break;
			}
		}

		if (found)
			continue;

		// -y test
		for (int i=0; i<maxMoves; i++)
		{
			// feasibility of the move
			// if (!(iy-i >= 0) || !(STATUSMAP(ix,iy-i,iz,NX,NY) == idold))
			if ((!optimizeGrids && (!(iy-i >= 0) || !(read3DVector<int>(status,ix,iy-i,iz,NX,NY,NZ) == idold))) ||
				( optimizeGrids && (!(iy-i >= 0) || !(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy-i,iz,NX,NY,NZ) == idold))))
			{
				found = false;
				break;
			}
			// if (iy-i >= 0 && iz+1 < NZ && STATUSMAP(ix,iy-i,iz+1,NX,NY) == idold)
			if ((!optimizeGrids && iy-i >= 0 && iz+1 < NZ && read3DVector<int>(status,ix,iy-i,iz+1,NX,NY,NZ) == idold) ||
				( optimizeGrids && iy-i >= 0 && iz+1 < NZ && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy-i,iz+1,NX,NY,NZ) == idold))
			{
				seeds_up.push_back( pair<pair<int,int>,int> (pair<int,int> (ix,iy-i),iz+1) );
				found = true;
				iy -= i;
				iz++;
				break;
			}
		}

		if (!found)
			break;
	}

	pair<pair<int,int>,int> *zz = allocateVector<pair<pair<int,int>,int>> (seeds_down.size()+seeds_up.size());
	int index = 0;

	int dim = (int)seeds_down.size()-1;
	for (int i=dim; i>=0; i--)
		zz[index++] = seeds_down[i];

	dim = (int)seeds_up.size();
	for (int i=0; i<dim; i++)
		zz[index++] = seeds_up[i];

	// cout << endl << "ok!";

	int z_max = seeds_up[seeds_up.size()-1].second;
	int z_min = seeds_down[seeds_down.size()-1].second;

	// start thread partitioning
	// split load evenly between threads
	int num_z = z_max-z_min+1;

	cout << "ok!";
	cout.flush();

	/*
	int old_num_threads = num_threads;

	if (num_threads >= num_z>>1)
	{
		num_threads = min(4, (num_z>>1) - 1);
		cout << endl << INFO << "Setting " << num_threads << " threads for floodfill";
	}
	*/

	#ifdef ENABLE_BOOST_THREADS
	boost::thread_group thdGroup;
	#endif

	queue<pair<pair<int,int>,int>> **upper_queues;
	queue<pair<pair<int,int>,int>> **lower_queues;

	upper_queues = allocateVector<queue<pair<pair<int,int>,int>>*> (num_threads-1);
	lower_queues = allocateVector<queue<pair<pair<int,int>,int>>*> (num_threads-1);

	int chunk = num_z/num_threads;
	int rem = num_z%num_threads;

	// cout << endl << "c " << chunk << " r " << rem;
	// setup split
	int start = 0;
	int stop  = 0;

	pair<int,int> old_new;
	pair<pair<int,int>,int> ind;

	old_new.first = idold;
	old_new.second = idnew;

	int upper_queue_index = 0;
	int lower_queue_index = 0;

	for (int j=0; j<num_threads; j++)
	{
		if (j == 0)
		{
			start = 0;
			stop = chunk;
		}
		else
		{
			start = stop;
			stop = start+chunk;
		}

		if (j < rem)
			stop++;

		pair<int,int> limits(start+z_min,stop+z_min);
		ind = zz[start];

		cout << endl << "--------------------------------------";
		cout << endl << start << " " << stop;
		cout << endl << start+z_min << " " << stop+z_min;
		cout << endl << "Start at " << ind.first.first << " " << ind.first.second << " " << ind.second;

		pair<queue<pair<pair<int,int>,int>>*, queue<pair<pair<int,int>,int>>*> queues;

		if (j == 0)
		{
			queues.second = NULL;
			upper_queues[upper_queue_index] = new queue<pair<pair<int,int>,int>>();
			queues.first = upper_queues[upper_queue_index];
			upper_queue_index++;
		}
		else if (j == num_threads - 1)
		{
			queues.first = NULL;
			lower_queues[lower_queue_index] = new queue<pair<pair<int,int>,int>>();
			queues.second = lower_queues[lower_queue_index];
			lower_queue_index++;
		}
		else
		{
			lower_queues[lower_queue_index] = new queue<pair<pair<int,int>,int>>();
			upper_queues[upper_queue_index] = new queue<pair<pair<int,int>,int>>();

			queues.first  = upper_queues[upper_queue_index];
			queues.second = lower_queues[lower_queue_index];

			upper_queue_index++;
			lower_queue_index++;
		}

		// cout << endl << "queues " << queues.first << " " << queues.second;

		#ifdef ENABLE_BOOST_THREADS
		thdGroup.create_thread(boost::bind(&Surface::floodFill3,this,ind,limits,old_new,queues,(queue<pair<pair<int,int>,int>>*)NULL));
		#else
		floodFill3(ind,limits,old_new,queues,(queue<pair<pair<int,int>,int>>*)NULL);
		#endif
	}

	#ifdef ENABLE_BOOST_THREADS
	// join; final part of the computation of the volume within the surface
	thdGroup.join_all();
	#endif

	pair<int,int> limits(-1,-1);

	//cout << endl;
	// run on the remaining queues
	for (int i=0; i<num_threads-1; i++)
		if (upper_queues[i] != NULL && upper_queues[i]->size() != 0)
		{
			//cout << "*";
			ind = upper_queues[i]->front();
			floodFill3(ind,limits,old_new,
					   pair<queue<pair<pair<int,int>,int>>*, queue<pair<pair<int,int>,int>>*> (NULL,NULL),
					   upper_queues[i]);
		}

		// run on the remaining queues
		for (int i=0; i<num_threads-1; i++)
			if (lower_queues[i] != NULL && lower_queues[i]->size() != 0)
			{
				//cout << "*";
				ind = lower_queues[i]->front();
				floodFill3(ind,limits,old_new,
						   pair<queue<pair<pair<int,int>,int>>*, queue<pair<pair<int,int>,int>>*> (NULL,NULL),
						   lower_queues[i]);
			}

			deleteVector<pair<pair<int,int>,int>>(zz);

		// delete queues
		for (int i=0; i<num_threads-1; i++)
			if (lower_queues[i] != NULL)
				delete lower_queues[i];

	for (int i=0; i<num_threads-1; i++)
		if (upper_queues[i] != NULL)
			delete upper_queues[i];

	delete lower_queues;
	delete upper_queues;

	// num_threads = old_num_threads;
}


// partial flood fill routine to be run in parallel
void Surface::floodFill3(	pair<pair<int,int>,int > ind,
							pair<int,int> z_limits,
							pair<int,int> old_new,
							pair<queue<pair<pair<int,int>,int>>*,queue<pair<pair<int,int>,int>>*> queues,
							queue<pair<pair<int,int>,int>> *in_queue)
{
	int idold = old_new.first;
	int idnew = old_new.second;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int *status = delphi->status;
	int **bilevel_status = delphi->bilevel_status;

	int ix = ind.first.first;
	int iy = ind.first.second;
	int iz = ind.second;

	// already visited, stops
	// if (STATUSMAP(ix,iy,iz,NX,NY) == idnew)
	if (!optimizeGrids)
	{
		if (read3DVector<int>(status,ix,iy,iz,NX,NY,NZ) == idnew)
			return;
	}
	else
	{
		if (readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz,NX,NY,NZ) == idnew)
			return;
	}

	pair<int,int> duplo;
	pair<pair<int,int>,int>  triplet;

	int iz_min = z_limits.first;
	int iz_max = z_limits.second;

	queue<pair<pair<int,int>,int>> *z_queue;

	if (in_queue == NULL)
		z_queue = new queue<pair<pair<int,int>,int>>();
	else
		z_queue = in_queue;

	queue<pair<pair<int,int>,int>> *out_up_xy = queues.first;
	queue<pair<pair<int,int>,int>> *out_down_xy = queues.second;

	if (in_queue == NULL)
		z_queue->push( pair<pair<int,int>,int> (pair<int,int>(ix,iy),iz) );

	while (z_queue->size() != 0)
	{
		triplet = z_queue->front();
		z_queue->pop();

		int iz = triplet.second;
		int ix = triplet.first.first;
		int iy = triplet.first.second;

		// if (STATUSMAP(ix,iy,iz,NX,NY) == idnew)
		if (!optimizeGrids)
		{
			if (read3DVector<int>(status,ix,iy,iz,NX,NY,NZ) == idnew)
				continue;
		}
		else
		{
			if (readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz,NX,NY,NZ) == idnew)
				continue;
		}

		//cout << endl << iz;
		stack<pair<int,int>> xy_stack;
		xy_stack.push(pair<int,int >(ix,iy));
		int y1;

		bool spanLeft=0, spanRight=0;

		while (xy_stack.size() != 0)
		{
			duplo = xy_stack.top();
			xy_stack.pop();

			ix = duplo.first;
			iy = duplo.second;

			spanLeft = spanRight = 0;

			y1 = iy;
			if (!optimizeGrids)
			{
				// while (y1 >= 0 && STATUSMAP(ix,y1,iz,NX,NY) == idold)
				while (y1 >= 0 && read3DVector<int>(status,ix,y1,iz,NX,NY,NZ) == idold)
					y1--;
			}
			else
			{
				while (y1 >= 0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,y1,iz,NX,NY,NZ) == idold)
					y1--;
			}
			y1++;

			// while (y1 < NY && STATUSMAP(ix,y1,iz,NX,NY) == idold)
			while (y1 < NY)
			{
				if (!optimizeGrids)
				{
					if (read3DVector<int>(status,ix,y1,iz,NX,NY,NZ) != idold)
						break;
				}
				else
				{
					if (readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,y1,iz,NX,NY,NZ) != idold)
						break;
				}

				// class mutex to protect writing
				{
					#ifdef ENABLE_BOOST_THREADS
					boost::mutex::scoped_lock scopedLock(mutex);
					#endif

					// STATUSMAP(ix,y1,iz,NX,NY) = idnew;
					if (!optimizeGrids)
						write3DVector<int>(status,idnew,ix,y1,iz,NX,NY,NZ);
					else
						writeBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,idnew,ix,y1,iz,NX,NY,NZ);

					// if this is a cavity then update the cavity list
					if (idnew >= 4)
					{
						// current cavity vector
						vector<int*> *vec;

						// first time this cavity is encountered
						if (delphi->cavitiesVec->size() <= idnew - 4)
						{
							//cout << endl << INFO << "Detected cavity " << idnew-4;
							vec = new vector<int*>();
							if (vec == NULL)
							{
								cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
								cout.flush();
								exit(-1);
							}
							delphi->cavitiesVec->push_back(vec);
						}
						else
						{
							vec = delphi->cavitiesVec->at(idnew - 4);
						}

						int *v = allocateVector<int>(3);
						if (v == NULL)
						{
							cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
							cout.flush();
							exit(-1);
						}
						v[0] = ix;
						v[1] = y1;
						v[2] = iz;
						// printf("\n %d %d %d (%d,%d)",cix,ciy,ciz,idnew,idold);
						vec->push_back(v);
					}
				}

				// can go up or not?
				// if (iz < NZ-1 && STATUSMAP(ix,y1,iz+1,NX,NY) == idold)
				if ((!optimizeGrids && iz < NZ-1 && read3DVector<int>(status,ix,y1,iz+1,NX,NY,NZ) == idold) ||
					( optimizeGrids && iz < NZ-1 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,y1,iz+1,NX,NY,NZ) == idold))
				{
					if (iz+1 == -1 || iz+1 < iz_max)
						z_queue->push( pair<pair<int,int>,int> (pair<int,int>(ix,y1),iz+1) );
					// postpone
					else
						out_up_xy->push( pair<pair<int,int>,int> (pair<int,int>(ix,y1),iz+1) );
				}
				// can go down or not?
				// if (iz > 0 && STATUSMAP(ix,y1,iz-1,NX,NY) == idold)
				if ((!optimizeGrids && iz > 0 && read3DVector<int>(status,ix,y1,iz-1,NX,NY,NZ) == idold) ||
					( optimizeGrids && iz > 0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,y1,iz-1,NX,NY,NZ) == idold))
				{
					if (iz-1 == -1 || iz-1 >= iz_min)
						z_queue->push( pair<pair<int,int>,int> (pair<int,int> (ix,y1),iz-1) );
					// postpone
					else
						out_down_xy->push( pair<pair<int,int>,int> (pair<int,int>(ix,y1),iz-1) );
				}
				// if(!spanLeft && ix > 0 && STATUSMAP(ix-1,y1,iz,NX,NY) == idold)
				if ((!optimizeGrids && !spanLeft && ix > 0 && read3DVector<int>(status,ix-1,y1,iz,NX,NY,NZ) == idold) ||
					( optimizeGrids && !spanLeft && ix > 0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix-1,y1,iz,NX,NY,NZ) == idold))
				{
					xy_stack.push( pair<int,int> (ix-1,y1) );
					spanLeft = 1;
				}
				// else if(spanLeft && ix > 0 && STATUSMAP(ix-1,y1,iz,NX,NY) != idold)
				else
				{
					if ((!optimizeGrids && spanLeft && ix > 0 && read3DVector<int>(status,ix-1,y1,iz,NX,NY,NZ) != idold) ||
						( optimizeGrids && spanLeft && ix > 0 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix-1,y1,iz,NX,NY,NZ) != idold))
					{
						spanLeft = 0;
					}
				}
				// if(!spanRight && ix < NX-1 && STATUSMAP(ix+1,y1,iz,NX,NY) == idold)
				if ((!optimizeGrids && !spanRight && ix < NX-1 && read3DVector<int>(status,ix+1,y1,iz,NX,NY,NZ) == idold) ||
					( optimizeGrids && !spanRight && ix < NX-1 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix+1,y1,iz,NX,NY,NZ) == idold))
				{
					xy_stack.push( pair<int,int> (ix+1,y1) );
					spanRight = 1;
				}
				// else if(spanRight && ix < NX-1 && STATUSMAP(ix+1,y1,iz,NX,NY) != idold)
				else
				{
					if ((!optimizeGrids && spanRight && ix < NX-1 && read3DVector<int>(status,ix+1,y1,iz,NX,NY,NZ) != idold) ||
						( optimizeGrids && spanRight && ix < NX-1 && readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix+1,y1,iz,NX,NY,NZ) != idold))
					{
						spanRight = 0;
					}
				}
				y1++;
			}
		}
	}

	if (in_queue == NULL)
		delete z_queue;
}


// slightly more cache friendly along y than along x
void Surface::floodFill2(int ix, int iy, int iz, int idold, int idnew)
{
	queue<pair<pair<int,int>,int>> z_queue;
	pair<int,int> duplo;
	pair<pair<int,int>,int> triplet;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int *status = delphi->status;


	z_queue.push( pair<pair<int,int>,int> (pair<int,int>(ix,iy),iz) );

	while (z_queue.size() != 0)
	{
		triplet = z_queue.front();
		z_queue.pop();

		int iz = triplet.second;
		int ix = triplet.first.first;
		int iy = triplet.first.second;

		// if (STATUSMAP(ix,iy,iz,NX,NY) == idnew)
		if (read3DVector<int>(status,ix,iy,iz,NX,NY,NZ) == idnew)
			continue;

		stack<pair<int,int>> xy_stack;
		xy_stack.push( pair<int,int> (ix,iy) );
		int y1;

		bool spanLeft = 0, spanRight = 0;

		while (xy_stack.size() != 0)
		{
			duplo = xy_stack.top();
			xy_stack.pop();

			ix = duplo.first;
			iy = duplo.second;

			spanLeft = spanRight = 0;

			y1 = iy;
			// while(y1 >= 0 && STATUSMAP(ix,y1,iz,NX,NY) == idold)
			while (y1 >= 0 && read3DVector<int>(status,ix,y1,iz,NX,NY,NZ) == idold)
				y1--;
			y1++;

			// while (y1 < NY && STATUSMAP(ix,y1,iz,NX,NY) == idold)
			while (y1 < NY && read3DVector<int>(status,ix,y1,iz,NX,NY,NZ) == idold)
			{
				{
					// STATUSMAP(ix,y1,iz,NX,NY) = idnew;
					write3DVector<int>(status,idnew,ix,y1,iz,NX,NY,NZ);

					// if this is a cavity then update the cavity list
					if (idnew >= 4)
					{
						// current cavity vector
						vector<int*> *vec;

						// first time this cavity is encountered
						if (delphi->cavitiesVec->size() < idnew - 4 + 1)
						{
							// cout << endl << INFO << "Detected cavity " << idnew-4;
							vec = new vector<int*>();
							if (vec == NULL)
							{
								cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
								cout.flush();
								exit(-1);
							}
							delphi->cavitiesVec->push_back(vec);
						}
						else
						{
							vec = delphi->cavitiesVec->at(idnew - 4);
						}

						int *v = allocateVector<int>(3);
						if (v == NULL)
						{
							cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
							cout.flush();
							exit(-1);
						}
						v[0] = ix;
						v[1] = y1;
						v[2] = iz;
						// printf("\n %d %d %d (%d,%d)",cix,ciy,ciz,idnew,idold);
						vec->push_back(v);
					}
				}
				// can go up or not?
				// if (iz < NZ-1 && STATUSMAP(ix,y1,iz+1,NX,NY) == idold)
				if (iz < NZ-1 && read3DVector<int>(status,ix,y1,iz+1,NX,NY,NZ) == idold)
				{
					z_queue.push( pair<pair<int,int>,int> (pair<int,int> (ix,y1),iz+1) );
				}
				// can go down or not?
				// if (iz > 0 && STATUSMAP(ix,y1,iz-1,NX,NY) == idold)
				if (iz > 0 && read3DVector<int>(status,ix,y1,iz-1,NX,NY,NZ) == idold)
				{
					z_queue.push( pair<pair<int,int>,int> (pair<int,int> (ix,y1),iz-1) );
				}
				// if (!spanLeft && ix > 0 && STATUSMAP(ix-1,y1,iz,NX,NY) == idold)
				if (!spanLeft && ix > 0 && read3DVector<int>(status,ix-1,y1,iz,NX,NY,NZ) == idold)
				{
					xy_stack.push( pair<int,int> (ix-1,y1) );
					spanLeft = 1;
				}
				// else if (spanLeft && ix > 0 && STATUSMAP(ix-1,y1,iz,NX,NY) != idold)
				else if (spanLeft && ix > 0 && read3DVector<int>(status,ix-1,y1,iz,NX,NY,NZ) != idold)
				{
					spanLeft = 0;
				}
				// if (!spanRight && ix < NX-1 && STATUSMAP(ix+1,y1,iz,NX,NY) == idold)
				if (!spanRight && ix < NX-1 && read3DVector<int>(status,ix+1,y1,iz,NX,NY,NZ) == idold)
				{
					xy_stack.push( pair<int,int> (ix+1,y1) );
					spanRight = 1;
				}
				// else if(spanRight && ix < NX-1 && STATUSMAP(ix+1,y1,iz,NX,NY) != idold)
				else if (spanRight && ix < NX-1 && read3DVector<int>(status,ix+1,y1,iz,NX,NY,NZ) != idold)
				{
					spanRight = 0;
				}
				y1++;
			}
		}
	}
}


// slightly more cache friendly along y than along x
void Surface::floodFillWithBilevelStatusMap2(int i, int j, int k, int idold, int idnew)
{
	queue<pair<pair<int,int>,int>> z_queue;
	pair<int,int> duplo;
	pair<pair<int,int>,int> triplet;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;
	int64_t coarse_nx = getCoarseN(NX);
	int64_t coarse_ny = getCoarseN(NY);

	int **bilevel_status = delphi->bilevel_status;


	z_queue.push( pair<pair<int,int>,int> (pair<int,int>(i,j),k) );

	while (z_queue.size() != 0)
	{
		triplet = z_queue.front();
		z_queue.pop();

		int i = triplet.first.first;
		int j = triplet.first.second;
		int k = triplet.second;

		int64_t coarse_i = getCoarseID(NX, i);
		int64_t coarse_j = getCoarseID(NY, j);
		int64_t coarse_k = getCoarseID(NZ, k);
		int64_t fine_j = getFineID(coarse_j, j);
		int64_t fine_k = getFineID(coarse_k, k);
		coarse_j *= coarse_nx;
		coarse_k *= coarse_ny*coarse_nx;
		int64_t coarse_ijk = coarse_k + coarse_j + coarse_i;

		if (bilevel_status[coarse_ijk] != NULL)
		{
			int fine_i = getFineID(coarse_i, i);
			int fine_index = getUnrolledFineID(fine_i, fine_j, fine_k);
			if (bilevel_status[coarse_ijk][fine_index] == idnew)
				continue;
		}
		// computation of k+1, k-1 of coarse cells of bilevel status map
		int64_t coarse_k1  = getCoarseID(NZ, k+1);
		int64_t coarse_k_1 = getCoarseID(NZ, k-1);
		int64_t fine_k1  = getFineID(coarse_k1, k+1);
		int64_t fine_k_1 = getFineID(coarse_k_1, k-1);
		coarse_k1  *= coarse_ny*coarse_nx;
		coarse_k_1 *= coarse_ny*coarse_nx;

		stack<pair<int,int>> xy_stack;
		xy_stack.push( pair<int,int> (i,j) );
		int i1;

		bool spanLeft = 0, spanRight = 0;

		while (xy_stack.size() != 0)
		{
			duplo = xy_stack.top();
			xy_stack.pop();

			i = duplo.first;
			j = duplo.second;

			int64_t coarse_j   = getCoarseID(NY, j);
			int64_t coarse_j1  = getCoarseID(NY, j+1);
			int64_t coarse_j_1 = getCoarseID(NY, j-1);
			int64_t fine_j   = getFineID(coarse_j, j);
			int64_t fine_j1  = getFineID(coarse_j1, j+1);
			int64_t fine_j_1 = getFineID(coarse_j_1, j-1);
			coarse_j   *= coarse_nx;
			coarse_j1  *= coarse_nx;
			coarse_j_1 *= coarse_nx;
			int64_t coarse_jk = coarse_k + coarse_j;

			spanLeft = spanRight = 0;

			i1 = i;
			while (i1 >= 0) // read3DVector<int>(status,i1,j,k,NX,NY,NZ) == idold
			{
				int64_t coarse_i1 = getCoarseID(NX, i1);
				int64_t coarse_i1jk = coarse_jk + coarse_i1;
				if (bilevel_status[coarse_i1jk] != NULL)
				{
					int64_t fine_i1 = getFineID(coarse_i1, i1);
					int64_t fine_index = getUnrolledFineID(fine_i1, fine_j, fine_k);
					if (bilevel_status[coarse_i1jk][fine_index] == idold)
						i1--;
					else
						break;
				}
				else
					if (idold == STATUS_POINT_TEMPORARY_OUT)
						i1--;
				else
					break;
			}
			i1++;

			while (i1 < NX)
			{
				int64_t coarse_i1 = getCoarseID(NX, i1);
				int64_t coarse_i1jk = coarse_jk + coarse_i1;

				if (bilevel_status[coarse_i1jk] == NULL)
				{
					if (idold == STATUS_POINT_TEMPORARY_OUT)
					{
						bilevel_status[coarse_i1jk] = allocateBilevelMinigridCells<int>(STATUS_POINT_TEMPORARY_OUT);
					}
					else
						break;
				}

				int64_t fine_i1 = getFineID(coarse_i1, i1);
				int64_t fine_index = getUnrolledFineID(fine_i1, fine_j, fine_k);

				if (bilevel_status[coarse_i1jk][fine_index] != idold)
					break;

				bilevel_status[coarse_i1jk][fine_index] = idnew;

				// if this is a cavity then update the cavity list
				if (idnew >= 4)
				{
					// current cavity vector
					vector<int*> *vec;

					// first time this cavity is encountered
					if (delphi->cavitiesVec->size() < idnew - 4 + 1)
					{
						// cout << endl << INFO << "Detected cavity " << idnew-4;
						vec = new vector<int*>();
						if (vec == NULL)
						{
							cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
							cout.flush();
							exit(-1);
						}
						delphi->cavitiesVec->push_back(vec);
					}
					else
					{
						vec = delphi->cavitiesVec->at(idnew - 4);
					}

					int *v = allocateVector<int>(3);
					if (v == NULL)
					{
						cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
						cout.flush();
						exit(-1);
					}
					v[0] = i1;
					v[1] = j;
					v[2] = k;
					vec->push_back(v);
				}

				// can go up or not?
				// if (k < NZ-1 && STATUSMAP(i1,j,k+1,NX,NY) == idold)
				if (k < NZ-1) // read3DVector<int>(status,i1,j,k+1,NX,NY,NZ) == idold
				{
					int64_t coarse_i1jk1 = coarse_k1 + coarse_j + coarse_i1;
					if (bilevel_status[coarse_i1jk1] != NULL)
					{
						int64_t fine_index = getUnrolledFineID(fine_i1, fine_j, fine_k1);
						if (bilevel_status[coarse_i1jk1][fine_index] == idold)
							z_queue.push( pair<pair<int,int>,int> (pair<int,int> (i1,j),k+1) );
					}
					else
						if (idold == STATUS_POINT_TEMPORARY_OUT)
							z_queue.push( pair<pair<int,int>,int> (pair<int,int> (i1,j),k+1) );

				}

				// can go down or not?
				// if (k > 0 && STATUSMAP(i1,j,k-1,NX,NY) == idold)
				if (k > 0) // read3DVector<int>(status,i1,j,k-1,NX,NY,NZ) == idold)
				{
					int64_t coarse_i1jk_1 = coarse_k_1 + coarse_j + coarse_i1;
					if (bilevel_status[coarse_i1jk_1] != NULL)
					{
						int64_t fine_index = getUnrolledFineID(fine_i1, fine_j, fine_k_1);
						if (bilevel_status[coarse_i1jk_1][fine_index] == idold)
							z_queue.push( pair<pair<int,int>,int> (pair<int,int> (i1,j),k-1) );
					}
					else
						if (idold == STATUS_POINT_TEMPORARY_OUT)
							z_queue.push( pair<pair<int,int>,int> (pair<int,int> (i1,j),k-1) );
				}
				// if (!spanLeft && j > 0 && STATUSMAP(i1,j-1,k,NX,NY) == idold)
				if (!spanLeft && j > 0) // read3DVector<int>(status,i1,j-1,k,NX,NY,NZ) == idold)
				{
					int64_t coarse_i1j_1k = coarse_k + coarse_j_1 + coarse_i1;
					if (bilevel_status[coarse_i1j_1k] != NULL)
					{
						int64_t fine_index = getUnrolledFineID(fine_i1, fine_j_1, fine_k);
						if (bilevel_status[coarse_i1j_1k][fine_index] == idold)
						{
							xy_stack.push( pair<int,int> (i1,j-1) );
							spanLeft = 1;
						}
					}
					else if (idold == STATUS_POINT_TEMPORARY_OUT)
					{
						xy_stack.push( pair<int,int> (i1,j-1) );
						spanLeft = 1;
					}
				}
				// else if (spanLeft && j > 0 && STATUSMAP(i1,j-1,k,NX,NY) != idold)
				else if (spanLeft && j > 0) // read3DVector<int>(status,i1,j-1,k,NX,NY,NZ) != idold)
				{
					int64_t coarse_i1j_1k = coarse_k + coarse_j_1 + coarse_i1;
					if (bilevel_status[coarse_i1j_1k] != NULL)
					{
						int64_t fine_index = getUnrolledFineID(fine_i1, fine_j_1, fine_k);
						if (bilevel_status[coarse_i1j_1k][fine_index] != idold)
						{
							spanLeft = 0;
						}
					}
					else
					{
						if (idold != STATUS_POINT_TEMPORARY_OUT)
							spanLeft = 0;
					}
				}
				// if (!spanRight && j < NY-1 && STATUSMAP(i1,j+1,k,NX,NY) == idold)
				if (!spanRight && j < NY-1) // && read3DVector<int>(status,i1,j+1,k,NX,NY,NZ) == idold)
				{
					int64_t coarse_i1j1k = coarse_k + coarse_j1 + coarse_i1;
					if (bilevel_status[coarse_i1j1k] != NULL)
					{
						int64_t fine_index = getUnrolledFineID(fine_i1, fine_j1, fine_k);
						if (bilevel_status[coarse_i1j1k][fine_index] == idold)
						{
							xy_stack.push( pair<int,int> (i1,j+1) );
							spanRight = 1;
						}
					}
					else if (idold == STATUS_POINT_TEMPORARY_OUT)
					{
						xy_stack.push( pair<int,int> (i1,j+1) );
						spanRight = 1;
					}
				}
				// else if (spanRight && j < NY-1 && STATUSMAP(i1,j+1,k,NX,NY) != idold)
				else if (spanRight && j < NY-1) // read3DVector<int>(status,i1,j+1,k,NX,NY,NZ) != idold)
				{
					int64_t coarse_i1j1k = coarse_k + coarse_j1 + coarse_i1;
					if (bilevel_status[coarse_i1j1k] != NULL)
					{
						int64_t fine_j1 = getFineID(getCoarseID(NY, j+1), j+1);
						int64_t fine_index = getUnrolledFineID(fine_i1, fine_j1, fine_k);
						if (bilevel_status[coarse_i1j1k][fine_index] != idold)
						{
							spanRight = 0;
						}
					}
					else
					{
						if (idold != STATUS_POINT_TEMPORARY_OUT)
							spanRight = 0;
					}
				}
				i1++;
			}
		}
	}
}


void Surface::buildSternLayer()
{
	// for each atom build a bounding cube
	// compute in/out info for each atom
	// for (int i=0; i<delphi->numAtoms; i++)
	for (int i=0; i<delphi->atoms.size(); i++)
	{
		double *sphere_center = delphi->atoms[i].pos;
		double radius = delphi->atoms[i].radius + sternLayer;
		double ref = radius*radius;

		// compute the object's bounding box
		double downx = sphere_center[0] - radius;
		double downy = sphere_center[1] - radius;
		double downz = sphere_center[2] - radius;

		double upx = sphere_center[0] + radius;
		double upy = sphere_center[1] + radius;
		double upz = sphere_center[2] + radius;

		// Determine which are the grid cells that are
		// occupied by the bounding box of the object
		int64_t ix_start = (int64_t)rintp((downx-delphi->xmin)*delphi->scale);
		int64_t iy_start = (int64_t)rintp((downy-delphi->ymin)*delphi->scale);
		int64_t iz_start = (int64_t)rintp((downz-delphi->zmin)*delphi->scale);

		int64_t ix_end = (int64_t)rintp((upx-delphi->xmin)*delphi->scale);
		int64_t iy_end = (int64_t)rintp((upy-delphi->ymin)*delphi->scale);
		int64_t iz_end = (int64_t)rintp((upz-delphi->zmin)*delphi->scale);

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

		if (ix_end >= delphi->nx)
			ix_end = delphi->nx - 1;
		if (iy_end >= delphi->ny)
			iy_end = delphi->ny - 1;
		if (iz_end >= delphi->nz)
			iz_end = delphi->nz - 1;

		if (ix_start >= delphi->nx)
			ix_start = delphi->nx - 1;
		if (iy_start >= delphi->ny)
			iy_start = delphi->ny - 1;
		if (iz_start >= delphi->nz)
			iz_start = delphi->nz - 1;

		for (int64_t iz=iz_start; iz<=iz_end; iz++)
			for (int64_t iy=iy_start; iy<=iy_end; iy++)
				for (int64_t ix=ix_start; ix<=ix_end; ix++)
				{
					//bool val = delphi->IDEBMAP(ix,iy,iz,delphi->nx,delphi->ny);
					bool val = read3DVector<bool>(delphi->idebmap,ix,iy,iz,delphi->nx,delphi->ny,delphi->nz);

					// outside is val = true
					if (val)
					{
						double v[3],dist;
						// if distance wrt sphere center is less then radius, mark as in
						v[0] = delphi->x[ix];
						v[1] = delphi->y[iy];
						v[2] = delphi->z[iz];
						DIST2(dist,v,sphere_center);
						// if inside the sphere (radius+layer) then mark as inside
						if (dist < ref)
						{
							//delphi->IDEBMAP(ix,iy,iz,delphi->nx,delphi->ny) = false;
							write3DVector<bool>(delphi->idebmap,false,ix,iy,iz,delphi->nx,delphi->ny,delphi->nz);
						}
					}
				}
	}
}


bool Surface::getProjection(double p[3], double *proj1, double *proj2,
							double *proj3, double *normal1, double *normal2, double *normal3)
{
	cout << endl << WARN << "Projection not supported!";
	return false;
}


void Surface::getRayIntersection(double p1[3], double p2[3], vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, bool computeNormals, int thread_id)
{
	cout << endl << WARN << "Ray-Patch intersection not supported!";
}


void Surface::projector(int start, int end)
{	
	for (int i=start; i<end; i++)
	{
		double gridPoint[3];
		gridPoint[0] = delphi->x[delphi->ibgp[3*i]];
		gridPoint[1] = delphi->y[delphi->ibgp[3*i+1]];
		gridPoint[2] = delphi->z[delphi->ibgp[3*i+2]];

		// project on the surface
		if (bgp_type[i] == EXTERNAL_BGP)
		{		
			getProjection(gridPoint,
						  &(delphi->scspos[3*i]),
						  &(delphi->scspos[3*i+1]),
						  &(delphi->scspos[3*i+2]),
						  &(delphi->scsnor[3*i]),
						  &(delphi->scsnor[3*i+1]),
						  &(delphi->scsnor[3*i+2]));
		}
		// apply inner projection routine.
		// get the nearest two atoms and project on the weighted voronoi plane
		// this routine is only valid for molecules where one has a notion of atom
		else
		{
			double minDist = INFINITY, secondDist = INFINITY;
			int first = -1, second = -1;

			// move from grid point of delphi grid to auxiliary atoms grid
			int64_t ix = (int64_t)rintp((gridPoint[0]-gxmin)*gscale);
			int64_t iy = (int64_t)rintp((gridPoint[1]-gymin)*gscale);
			int64_t iz = (int64_t)rintp((gridPoint[2]-gzmin)*gscale);

			// get the nearest atoms to set the dielectric constant
			// the dielectric constant is mapped according to the additively weighted voronoi diagram
			// that is the signed distance from the point p is ||p-c||^2-r^2 where c is the center
			// of the atom and r is the radius. The minimum signed distance wins.
			for (int k=0; k<SHIFT_MAP; k++)
			{
				int64_t cx = ix + shift_map[k][0];
				int64_t cy = iy + shift_map[k][1];
				int64_t cz = iz + shift_map[k][2];

				if (cx >= ggrid || cy >= ggrid || cz >= ggrid || cx < 0 || cy < 0 || cz < 0)
					continue;

				int num_atoms = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ].size();

				for (int j=0; j<num_atoms; j++)
				{
					int atom_index = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ][j];

					double signed_dist = 0;
					DIST2(signed_dist,delphi->atoms[atom_index].pos,gridPoint)

					double radius2 = delphi->atoms[atom_index].radius2;
					signed_dist -= radius2;

					if (signed_dist < minDist)
					{
						secondDist = minDist;
						minDist = signed_dist;
						second = first;
						first = atom_index;
					}
					else
					{
						if (signed_dist < secondDist)
						{
							secondDist = signed_dist;
							second = atom_index;
						}
					}
				}
			}
			
			// exit(-1);
			// compute the voronoi plane using power distance
			double w[4];
			double *c1,*c2;

			c1 = delphi->atoms[first].pos;
			c2 = delphi->atoms[second].pos;

			SUB(w,c2,c1)
			w[0] = 2*w[0];
			w[1] = 2*w[1];
			w[2] = 2*w[2];
			w[3] = (DOT(c1,c1)) - (DOT(c2,c2)) - delphi->atoms[first].radius*delphi->atoms[first].radius + delphi->atoms[second].radius*delphi->atoms[second].radius;

			// project on the plane
			double dist, proj[3];
			point2plane(gridPoint, w, &dist,proj);
			
			// fprintf(fp,"%f %f %f\n", gridPoint[0],gridPoint[1],gridPoint[2]);
			// printf("\n%f %f %f", proj[0],proj[1],proj[2]);
			
			delphi->scspos[3*i  ] = proj[0];
			delphi->scspos[3*i+1] = proj[1];
			delphi->scspos[3*i+2] = proj[2];

			delphi->scsnor[3*i  ] = gridPoint[0] - proj[0];
			delphi->scsnor[3*i+1] = gridPoint[1] - proj[1];
			delphi->scsnor[3*i+2] = gridPoint[2] - proj[2];
		}
	}
}


void Surface::initPatchBasedRayTracing (int64_t nxyz[3], int panels[2])
{
	if (getNumPatches() == 0)
	{
		cout << endl << WARN << "Cannot get this surface type with patch-based ray-tracing!";
		return;
	}

	int64_t NX = nxyz[0];
	int64_t NY = nxyz[1];
	int64_t NZ = nxyz[2];
	int64_t N_MAX = MAX(NX, MAX(NY, NZ));

	int numPanels = panels[0] + panels[1];

	int numPatches = getNumPatches();

	#if !defined(SINGLE_PASS_RT)

	#if !defined(MULTITHREADING)
	num_pixel_intersections = allocateVector<int>(N_MAX*N_MAX* numPanels);

	memset(num_pixel_intersections, 0, N_MAX*N_MAX* numPanels *sizeof(int));
	#else
	num_patch_intersections = allocateVector<int>(numPatches*numPanels);

	memset(num_patch_intersections, 0, numPatches* numPanels *sizeof(int));

	first_patch_intersection_index = allocateVector<int>(numPatches*numPanels);
	#endif

	#else // SINGLE_PASS_RT

	temp_intersections_buffer = new vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> [ numPatches*numPanels ];

	intersection_pixel_id = new vector<int64_t> [ numPatches*numPanels ];

	#endif // SINGLE_PASS_RT

	#if defined(RAY_VS_CELL_TESTS_CULLING) || defined(RAY_VS_PATCH_TESTS_CULLING)
	screen_buffer = new bool [ conf.numThreads*N_MAX*N_MAX ];
	memset(screen_buffer, 0, conf.numThreads*N_MAX*N_MAX *sizeof(bool));

	tagged_pixel_ids = new int [ conf.numThreads*N_MAX*N_MAX ];
	#endif
}


#if !defined(SINGLE_PASS_RT)

void Surface::allocateIntermediatePatchBuffers (int64_t nxyz[3], int panels[2], int potentialIntersections)
{
	if (getNumPatches() == 0)
	{
		return;
	}

	int64_t NX = nxyz[0];
	int64_t NY = nxyz[1];
	int64_t NZ = nxyz[2];

	int numPanels = panels[0] + panels[1];

	int numPatches = getNumPatches();


	#if !defined(MULTITHREADING)
	int64_t N_MAX = MAX(NX, MAX(NY, NZ));

	pixel_intersections = allocateVector<pair<VERTEX_TYPE,VERTEX_TYPE*>*> (N_MAX*N_MAX* numPanels);

	for (int64_t id=0; id < N_MAX*N_MAX* numPanels; id++)
	{
		pixel_intersections[id] = NULL;

		int numInt = num_pixel_intersections[id];

		if (numInt > 0)
			pixel_intersections[id] = allocateVector<pair<VERTEX_TYPE,VERTEX_TYPE*>> (numInt);
		num_pixel_intersections[id] = 0;
	}
	#else // MULTITHREADING
	temp_intersections_buffer = allocateVector<pair<VERTEX_TYPE,VERTEX_TYPE*>> (potentialIntersections);

	intersection_pixel_id = allocateVector<int64_t>(potentialIntersections);

	int count = 0;

	for (int it=0; it<numPatches; it++)
	{
		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				int id = it*numPanels + (panels[0]*int_phase + panel);

				first_patch_intersection_index[id] = count;
				count += num_patch_intersections[id];
				// resetting of patch intersections (they are refilled during the below patch vs rays intersections)
				num_patch_intersections[id] = 0;
			}
		}
	}
	#endif // MULTITHREADING
}

#endif // SINGLE_PASS_RT


#if !defined(SINGLE_PASS_RT)

void Surface::countPixelIntersections (int64_t nxyz[3], int panels[2])
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

	int numPatches = getNumPatches();


	#if defined(RAY_VS_CELL_TESTS_CULLING) || defined(RAY_VS_PATCH_TESTS_CULLING)
	delete tagged_pixel_ids;
	tagged_pixel_ids = NULL;

	delete screen_buffer;
	screen_buffer = NULL;
	#endif

	#ifdef MULTITHREADING
	num_pixel_intersections = allocateVector<int>(N_MAX*N_MAX* numPanels);

	memset(num_pixel_intersections, 0, N_MAX*N_MAX* numPanels *sizeof(int));

	pixel_intersections = allocateVector<pair<VERTEX_TYPE,VERTEX_TYPE*>*> (N_MAX*N_MAX* numPanels);

	for (int it=0; it<numPatches; it++)
	{
		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);
				int first_id = first_patch_intersection_index[patch_id];

				for (int i=0; i<num_patch_intersections[patch_id]; i++)
				{
					int64_t panel_pixel_id = intersection_pixel_id[ first_id + i ];
					++num_pixel_intersections[panel_pixel_id];
				}
			}
		}
	}
	for (int64_t id=0; id < N_MAX*N_MAX* numPanels; id++)
	{
		pixel_intersections[id] = NULL;

		int numInt = num_pixel_intersections[id];

		if (numInt > 0)
			pixel_intersections[id] = allocateVector<pair<VERTEX_TYPE,VERTEX_TYPE*>> (numInt);
		num_pixel_intersections[id] = 0;
	}
	#endif // MULTITHREADING
}

#else // SINGLE_PASS_RT

void Surface::countPixelIntersections (int64_t nxyz[3], int panels[2])
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

	int numPatches = getNumPatches();


	#if defined(RAY_VS_CELL_TESTS_CULLING) || defined(RAY_VS_PATCH_TESTS_CULLING)
	delete tagged_pixel_ids;
	tagged_pixel_ids = NULL;

	delete screen_buffer;
	screen_buffer = NULL;
	#endif

	num_pixel_intersections = allocateVector<int>(N_MAX*N_MAX* numPanels);

	memset(num_pixel_intersections, 0, N_MAX*N_MAX* numPanels *sizeof(int));

	pixel_intersections = allocateVector<pair<VERTEX_TYPE,VERTEX_TYPE*>*> (N_MAX*N_MAX* numPanels);

	for (int it=0; it<numPatches; it++)
	{
		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);

				for (int i=0; i<intersection_pixel_id[patch_id].size(); i++)
				{
					int64_t panel_pixel_id = intersection_pixel_id[ patch_id ][i];
					++num_pixel_intersections[panel_pixel_id];
				}
			}
		}
	}
	for (int64_t id=0; id < N_MAX*N_MAX* numPanels; id++)
	{
		pixel_intersections[id] = NULL;

		int numInt = num_pixel_intersections[id];

		if (numInt > 0)
			pixel_intersections[id] = allocateVector<pair<VERTEX_TYPE,VERTEX_TYPE*>> (numInt);
		num_pixel_intersections[id] = 0;
	}
}

#endif // SINGLE_PASS_RT


#if !defined(SINGLE_PASS_RT)

void Surface::copyPatchIntersections (int nxyz[3], int panels[2])
{
	if (getNumPatches() == 0)
	{
		return;
	}

	int numPanels = panels[0] + panels[1];


	#ifdef MULTITHREADING
	for (int it = 0; it<getNumPatches(); it++)
	{
		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);
				int first_id = first_patch_intersection_index[patch_id];

				for (int i=0; i<num_patch_intersections[patch_id]; i++)
				{
					int64_t panel_pixel_id = intersection_pixel_id[ first_id + i ];

					pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p1 = &temp_intersections_buffer[ first_id + i ];
					pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p2 = &pixel_intersections[panel_pixel_id][ num_pixel_intersections[panel_pixel_id]++ ];

					int_p2->first = int_p1->first;
					int_p2->second = int_p1->second;
				}
			}
		}
	}
	#endif // MULTITHREADING
}

#else // SINGLE_PASS_RT

void Surface::copyPatchIntersections (int64_t nxyz[3], int panels[2])
{
	if (getNumPatches() == 0)
	{
		return;
	}

	int numPanels = panels[0] + panels[1];


	for (int it = 0; it<getNumPatches(); it++)
	{
		for (int int_phase = 0; int_phase < 2; int_phase++)
		{
			for (int panel=0; panel < panels[ int_phase ]; panel++)
			{
				int patch_id = it*numPanels + (panels[0]*int_phase + panel);

				for (int i=0; i<intersection_pixel_id[patch_id].size(); i++)
				{
					int64_t panel_pixel_id = intersection_pixel_id[ patch_id ][i];

					// copy from patch-centric buffer to pixel-centric one
					pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p1 = &temp_intersections_buffer[ patch_id ][i];
					pair<VERTEX_TYPE,VERTEX_TYPE*> *int_p2 = &pixel_intersections[panel_pixel_id][ num_pixel_intersections[panel_pixel_id]++ ];

					int_p2->first = int_p1->first;
					int_p2->second = int_p1->second;
				}
			}
		}
	}
}

#endif // SINGLE_PASS_RT


void Surface::reorderPatchIntersections (int pixel_start, int pixel_end)
{
	if (getNumPatches() == 0)
	{
		return;
	}

	for (int id = pixel_start; id < pixel_end; id += conf.numThreads)
	{
		if (num_pixel_intersections[id] > 0)
		{
			sort(pixel_intersections[id],
				 pixel_intersections[id] + num_pixel_intersections[id], compKeepIndex);
		}
	}
}


// This routine exploits the buffer previously filled with intersections' data via
// the patch-based ray-tracing routines getPatchIntersectionData(...)
// SD PB_NEW also save the grid packets now
void Surface::setVerticesAndGridsWithIntersectionData (int thread_id, int nb, int start, int end, int iters_block, int jump, packet pack, packet gridPack)
{
	vector<pair<int,int>> intersection_indices;

	intersection_indices.reserve(2000);

	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *v_int = pack.first;
	vector<coordVec> *v_norm = pack.second;
	#else
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *v_int = pack.first;
	#else
	vector<compressedCoordNormPacket> *v_int = pack.first;
	#endif
	#endif

	// SD PB_NEW
	bool enableSaveGridIntersections = false;
	
	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *v_int_grid = nullptr;
	vector<coordVec> *v_norm_grid = nullptr;
	#else
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *v_int_grid = nullptr;
	#else
	vector<compressedCoordNormPacket> *v_int_grid = nullptr;
	#endif
	#endif
	if  (gridPack.first != nullptr)
	{
		enableSaveGridIntersections = true;
		
		v_int_grid = gridPack.first;
		#if !defined(COORD_NORM_PACKING)
		v_norm_grid = gridPack.second;
		#endif
	}
	// END PB_NEW
	
	double pa[3],pb[3];
	double dir;
	
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	int N_MAX = MAX(NX, MAX(NY, NZ));
	

	int panels[] = {3,3};
	if (!collectFaceIntersections && !delphi->buildStatus && !delphi->buildEpsMap)
		panels[0] = 0;
	else if (!collectFaceIntersections && delphi->buildStatus && !delphi->buildEpsMap)
		panels[0] = 1;
	if (!accurateTriangulation || isAvailableScalarField)
		panels[1] = 0;
	int numPanels = panels[0] + panels[1];

	// if cavities and espilon map are not necessary then skip this step at all
	if (!collectFaceIntersections && !delphi->buildStatus && !delphi->buildEpsMap)
	{
		// skip
	}
	// if only cavity detection is required, only the first panel is ray cast
	else if (!collectFaceIntersections && (delphi->buildStatus && panel != 0) && !delphi->buildEpsMap)
	{
		// skip panel 1,2
	}
	// if both cavity detection and epsilon map are required, 1 or 3 panels must be analyzed
	else if (collectFaceIntersections || (delphi->buildStatus && panel == 0) || delphi->buildEpsMap)
	{
		double last_vol_integral = 0.0;

		if (thread_id == 0)
			panelVolumeFlag[panel][0] = 1;

		if (panel == 0) {
			pa[0] = delphi->x[0];
			pb[0] = delphi->x[NX-1];
		}
		else if (panel == 1) {
			pa[2] = delphi->z[0];
			pb[2] = delphi->z[NZ-1];
		}
		else
		{
			pa[1] = delphi->y[0];
			pb[1] = delphi->y[NY-1];
		}
		for (int nn = start; nn < end; nn += jump)
		{
			for (int n = nn; n < min(end, nn + iters_block); n++)
			{
				if (panel == 0)
				{
					// rays along axis X -> smallest striding and maximal caching efficiency
					pa[2] = delphi->z[n];
					pb[2] = pa[2];
				}
				else if (panel == 1)
				{
					// rays along axis Z -> largest striding and minimal caching efficiency
					pa[1] = delphi->y[n];
					pb[1] = pa[1];
				}
				else
				{
					pa[2] = delphi->z[n];
					pb[2] = pa[2];
				}

				for (int m=0; m<nb; m++)
				{
					++threadTotalRays[thread_id];

					int id = (n*N_MAX+m)*numPanels + panel;
					int vsize = num_pixel_intersections[id];

					if (vsize == 0)
						continue;

					intersection_indices.clear();

					// check if repetitions allow the ray to enter and exist. If the
					// ray only enters then this is an error
					double lastValid = INFINITY;
					bool closure = false;

					#if !defined(NEW_INTERSECTION_FILTERING)
					for (int int_id = 0; int_id < vsize; int_id++)
					{
						double entry_point = pixel_intersections[id][int_id].first;

						int entry_int_id = int_id;
						// get the entering point
						if (fabs(entry_point - lastValid) >= EPS_INT)
						{
							closure = false;
							// search the exiting point of the ray
							while (1)
							{
								++int_id;
								if (int_id == vsize)
									break;

								if (fabs(pixel_intersections[id][int_id].first - entry_point) >= EPS_INT)
								{
									// ok the ray is exiting
									intersection_indices.push_back( pair<int,int>(entry_int_id,int_id) );
									lastValid = pixel_intersections[id][int_id].first;
									closure = true;
									break;
								}
							}
						}
					}

					// single or multiple tangent intersections
					if (lastValid == INFINITY)
					{
						continue;
					}

					#else // NEW_INTERSECTION_FILTERING

					// if vsize is even ..., otherwise closure remains false;
					if (((vsize>>1)<<1) == vsize)
					{
						for (int i=0; i<vsize; i+=2)
						{
							// Any two non-excessively close intersections are stored and considered later
							if (fabs(pixel_intersections[id][i+1].first - pixel_intersections[id][i].first) >= EPS_INT)
							{
								intersection_indices.push_back( pair<int,int>(i,i+1) );
							}
						}
						closure = true;
					}

					if (intersection_indices.size() == 0)
					{
						continue;
					}
					#endif // NEW_INTERSECTION_FILTERING

					// check finished.
					if (!closure)
					{
						threadPanelVolume[panel][thread_id] += last_vol_integral;

						#if defined(REPORT_FAILED_RAYS)
						if (panel == 0)
						{
							// rays along axis X -> smallest striding and maximal caching efficiency
							pa[1] = delphi->y[m];
							pb[1] = pa[1];
						}
						else if (panel == 1)
						{
							// rays along axis Z -> largest striding and minimal caching efficiency
							pa[0] = delphi->x[m];
							pb[0] = pa[0];
						}
						else
						{
							pa[0] = delphi->x[m];
							pb[0] = pa[0];
						}
						#endif
						// approximate the current ray with the previous one
						if (m > 0 && panel == 0)
						{
							// eps, ideb, status
							if (delphi->buildStatus && delphi->buildEpsMap)
							{
								for (int ix=0; ix<NX; ix++)
								{
									const int eee = read4DVector<int>(delphi->epsmap,ix,m-1,n,0,NX,NY,NZ,3);
									write4DVector<int>(delphi->epsmap,eee,ix,m,n,0,NX,NY,NZ,3);
								
									if (!optimizeGrids)
									{
										const int sss = read3DVector<int>(delphi->status,ix,m-1,n,NX,NY,NZ);
										write3DVector<int>(delphi->status,sss,ix,m,n,NX,NY,NZ);
									}
									else
									{
										const int sss = readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,m-1,n,NX,NY,NZ);
										writeBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,sss,ix,m,n,NX,NY,NZ);
									}
									const bool bbb = read3DVector<bool>(delphi->idebmap,ix,m-1,n,NX,NY,NZ);
									write3DVector<bool>(delphi->idebmap,bbb,ix,m,n,NX,NY,NZ);
								}
							}
							// eps, idebmap
							else if (delphi->buildEpsMap)
							{
								for (int ix=0; ix<NX; ix++)
								{
									const int eee = read4DVector<int>(delphi->epsmap,ix,m-1,n,0,NX,NY,NZ,3);
									write4DVector<int>(delphi->epsmap,eee,ix,m,n,0,NX,NY,NZ,3);

									const bool bbb = read3DVector<bool>(delphi->idebmap,ix,m-1,n,NX,NY,NZ);
									write3DVector<bool>(delphi->idebmap,bbb,ix,m,n,NX,NY,NZ);
								}
							}
							// status
							else if (delphi->buildStatus)
							{
								for (int ix=0; ix<NX; ix++)
								{
									if (!optimizeGrids)
									{
										const int sss = read3DVector<int>(delphi->status,ix,m-1,n,NX,NY,NZ);
										write3DVector<int>(delphi->status,sss,ix,m,n,NX,NY,NZ);
									}
									else
									{
										const int sss = readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,m-1,n,NX,NY,NZ);
										writeBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,sss,ix,m,n,NX,NY,NZ);
									}
								}
							}
						}
						else if (m > 0 && panel == 1 && delphi->buildEpsMap)
						{
							for (int iz=0; iz<NZ; iz++)
							{
								const int eee = read4DVector<int>(delphi->epsmap,m-1,n,iz,2,NX,NY,NZ,3);
								write4DVector<int>(delphi->epsmap,eee,m,n,iz,2,NX,NY,NZ,3);
							}
						}
						else if (m > 0 && panel == 2 && delphi->buildEpsMap)
						{
							for (int iy=0; iy<NY; iy++)
							{
								const int eee = read4DVector<int>(delphi->epsmap,m-1,iy,n,1,NX,NY,NZ,3);
								write4DVector<int>(delphi->epsmap,eee,m,iy,n,1,NX,NY,NZ,3);
							}
						}

						++threadFailedRays[thread_id];

						continue;
					}
				
					// if there are intersections mark the epsmap cells that are inside
					// by counting the number of intersections. This is done for each panel
					// that is in the x,y,z directions respectively.
					// Eg. Tracing the ray one obtains inside/out information
					// 000000000|-1-1-1-1-1-1-1-1-1-1-1-1|000000000000|-1-1-1-1-1|000000
					if (panel == 0)
					{
						// rays along axis X -> smallest striding and maximal caching efficiency
						pa[1] = delphi->y[m];
						pb[1] = pa[1];
					}
					else if (panel == 1)
					{
						// rays along axis Z -> largest striding and minimal caching efficiency
						pa[0] = delphi->x[m];
						pb[0] = pa[0];
					}
					else
					{
						pa[0] = delphi->x[m];
						pb[0] = pa[0];
					}
					double lim1,lim2;

					last_vol_integral = 0.0;

					if (panel == 0)
					{
						dir = pb[0] - pa[0];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[0] + dir*pixel_intersections[id][ int_id1 ].first;
							lim2 = pa[0] + dir*pixel_intersections[id][ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							int i1 = (int)rintp((lim1-delphi->xmin)*delphi->scale);
							int i2 = (int)rintp((lim2-delphi->xmin)*delphi->scale);

							i1 = (i1 < 0)   ? 0    : i1;
							i1 = (i1 >= NX) ? NX-1 : i1;
							i2 = (i2 < 0)   ? 0    : i2;
							i2 = (i2 >= NX) ? NX-1 : i2;

							// second inters. after the first one .. it should not occur
							if (i2 < i1)
							{
								continue;
							}
							// SD PB_NEW collect also grid intersection data
							if (enableSaveGridIntersections)
							{
								verticesBuffers[thread_id].push_back(lim1);
								verticesBuffers[thread_id].push_back(pa[1]);
								verticesBuffers[thread_id].push_back(pa[2]);
								VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
								verticesBuffers[thread_id].push_back(lim2);
								verticesBuffers[thread_id].push_back(pa[1]);
								verticesBuffers[thread_id].push_back(pa[2]);
								VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

								#if !defined(COORD_NORM_PACKING)
								{
									coordVec cv(i1,m,n, intersec1, X_DIR);
									v_int_grid->push_back(cv);
								}
								{
									coordVec cv(i2,m,n, intersec2, X_DIR);
									v_int_grid->push_back(cv);
								}
								if (computeNormals && providesAnalyticalNormals)
								{
									VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);
									if (normal != NULL)
									{
										coordVec cv(i1,m,n, normal, X_DIR);
										v_norm_grid->push_back(cv);
									}
									normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
									if (normal != NULL)
									{
										coordVec cv(i2,m,n, normal, X_DIR);
										v_norm_grid->push_back(cv);
									}
								}
								#else // COORD_NORM_PACKING
								VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);
								VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(i1,m,n, intersec1, normal1, X_DIR);
									#else
									compressedCoordNormPacket cnv(i1,m,n, intersec1, normal1, X_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(i2,m,n, intersec2, normal2, X_DIR);
									#else
									compressedCoordNormPacket cnv(i2,m,n, intersec2, normal2, X_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								#endif // COORD_NORM_PACKING
							}
							// END PB_NEW

							// eps, status, idebmap
							if (delphi->buildEpsMap && delphi->buildStatus)
							{
								// for (int ix=0; ix<NX; ix++)
								for (int ix=i1; ix<=i2; ix++)
								{
									// epsmap
									if (delphi->x[ix] <= lim2 - delphi->hside) // delphi->x[ix] >= lim1-delphi->hside is always true if ix>=i1
									{
										write4DVector<int>(delphi->epsmap,inside,ix,m,n,0,NX,NY,NZ,3);
									}
									// status and idebmap
									// in the loop arguments, "ix<=i2" because, below, it could be true that delphi->x[i2] <= lim2
									if (delphi->x[ix] <= lim2 && delphi->x[ix] >= lim1)
									{
										if (!optimizeGrids)
											write3DVector<int>(delphi->status,STATUS_POINT_INSIDE,ix,m,n,NX,NY,NZ);
										else
											writeBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,STATUS_POINT_INSIDE,ix,m,n,NX,NY,NZ);
										write3DVector<bool>(delphi->idebmap,false,ix,m,n,NX,NY,NZ);
									}
								}
							}
							// epsmap, idebmap
							else if (delphi->buildEpsMap)
							{
								// for (int ix=0; ix<NZ; ix++)
								for (int ix=i1; ix<=i2; ix++)
								{
									if (delphi->x[ix] <= lim2 - delphi->hside) // delphi->x[ix] >= lim1-delphi->hside is always true if ix>=i1
									{
										write4DVector<int>(delphi->epsmap,inside,ix,m,n,0,NX,NY,NZ,3);
									}
									if (delphi->x[ix] <= lim2 && delphi->x[ix] >= lim1)
									{
										write3DVector<bool>(delphi->idebmap,false,ix,m,n,NX,NY,NZ);
									}
								}
							}
							// only status
							else if (delphi->buildStatus)
							{
								// for (int ix=0; ix<NX; ix++)
								for (int ix=i1; ix<=i2; ix++)
								{
									if (delphi->x[ix] <= lim2 && delphi->x[ix] >= lim1)
									{
										if (!optimizeGrids)
											write3DVector<int>(delphi->status,STATUS_POINT_INSIDE,ix,m,n,NX,NY,NZ);
										else
											writeBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,STATUS_POINT_INSIDE,ix,m,n,NX,NY,NZ);
									}
								}
							}
						}
					}
					else if (panel == 1)
					{
						dir = pb[2] - pa[2];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[2] + dir*pixel_intersections[id][ int_id1 ].first;
							lim2 = pa[2] + dir*pixel_intersections[id][ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							int k1 = (int)rintp((lim1-delphi->zmin)*delphi->scale);
							int k2 = (int)rintp((lim2-delphi->zmin)*delphi->scale);

							k1 = (k1 < 0)   ? 0    : k1;
							k1 = (k1 >= NZ) ? NZ-1 : k1;
							k2 = (k2 < 0)   ? 0    : k2;
							k2 = (k2 >= NZ) ? NZ-1 : k2;

							// second inters. after the first one .. it should not occur
							if (k2 < k1)
							{
								continue;
							}
							// SD PB_NEW collect also grid intersection data
							if (enableSaveGridIntersections)
							{
								verticesBuffers[thread_id].push_back(pa[0]);
								verticesBuffers[thread_id].push_back(pa[1]);
								verticesBuffers[thread_id].push_back(lim1);
								VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
								verticesBuffers[thread_id].push_back(pa[0]);
								verticesBuffers[thread_id].push_back(pa[1]);
								verticesBuffers[thread_id].push_back(lim2);
								VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

								#if !defined(COORD_NORM_PACKING)
								{
									coordVec cv(m,n,k1, intersec1, Z_DIR);
									v_int_grid->push_back(cv);
								}
								{
									coordVec cv(m,n,k2, intersec2, Z_DIR);
									v_int_grid->push_back(cv);
								}
								if (computeNormals && providesAnalyticalNormals)
								{
									VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);

									if (normal != NULL)
									{
										coordVec cv(m,n,k1, normal, Z_DIR);
										v_norm_grid->push_back(cv);
									}
									normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
									if (normal != NULL)
									{
										coordVec cv(m,n,k2, normal, Z_DIR);
										v_norm_grid->push_back(cv);
									}
								}
								#else // COORD_NORM_PACKING
								VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);
								VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(m,n,k1, intersec1, normal1, Z_DIR);
									#else
									compressedCoordNormPacket cnv(m,n,k1, intersec1, normal1, Z_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(m,n,k2, intersec2, normal2, Z_DIR);
									#else
									compressedCoordNormPacket cnv(m,n,k2, intersec2, normal2, Z_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								#endif // COORD_NORM_PACKING
							}
							if (delphi->buildEpsMap)
							{
								// for (int iz=0; iz<NZ; iz++)
								for (int iz=k1; iz<k2; iz++)
								{
									// "iz<i2" (and not "iz<=i2"): it cannot occur that delphi->z[ i2 ] <= lim2 - delphi->hside
									if (delphi->z[iz] <= lim2 - delphi->hside) // delphi->z[iz] >= lim1 - delphi->hside is always true if iz>=i1
									{
										write4DVector<int>(delphi->epsmap,inside,m,n,iz,2,NX,NY,NZ,3);
									}
								}
							}
						}
					}
					else
					{
						dir = pb[1] - pa[1];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[1] + dir*pixel_intersections[id][ int_id1 ].first;
							lim2 = pa[1] + dir*pixel_intersections[id][ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							int j1 = (int)rintp((lim1-delphi->ymin)*delphi->scale);
							int j2 = (int)rintp((lim2-delphi->ymin)*delphi->scale);

							j1 = (j1 < 0)   ? 0    : j1;
							j1 = (j1 >= NY) ? NY-1 : j1;
							j2 = (j2 < 0)   ? 0    : j2;
							j2 = (j2 >= NY) ? NY-1 : j2;

							// second inters. after the first one .. it should not occur
							if (j2 < j1)
							{
								continue;
							}
							// SD PB_NEW collect also grid intersection data
							if (enableSaveGridIntersections)
							{
								verticesBuffers[thread_id].push_back(pa[0]);
								verticesBuffers[thread_id].push_back(lim1);
								verticesBuffers[thread_id].push_back(pa[2]);
								VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
								verticesBuffers[thread_id].push_back(pa[0]);
								verticesBuffers[thread_id].push_back(lim2);
								verticesBuffers[thread_id].push_back(pa[2]);
								VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

								#if !defined(COORD_NORM_PACKING)
								{
									coordVec cv(m,j1,n, intersec1, Y_DIR);
									v_int_grid->push_back(cv);
								}
								{
									coordVec cv(m,j2,n, intersec2, Y_DIR);
									v_int_grid->push_back(cv);
								}
								if (computeNormals && providesAnalyticalNormals)
								{
									VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);

									if (normal != NULL)
									{
										coordVec cv(m,j1,n, normal, Y_DIR);
										v_norm_grid->push_back(cv);
									}
									normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
									if (normal != NULL)
									{
										coordVec cv(m,j2,n, normal, Y_DIR);
										v_norm_grid->push_back(cv);
									}
								}
								#else // COORD_NORM_PACKING
								VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);
								VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(m,j1,n, intersec1, normal1, Y_DIR);
									#else
									compressedCoordNormPacket cnv(m,j1,n, intersec1, normal1, Y_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(m,j2,n, intersec2, normal2, Y_DIR);
									#else
									compressedCoordNormPacket cnv(m,j2,n, intersec2, normal2, Y_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								#endif // COORD_NORM_PACKING
							}
							// END PB_NEW

							if (delphi->buildEpsMap)
							{
								// for (int iy=0; iy<NY; iy++)
								for (int iy=j1; iy<j2; iy++)
								{
									// "iy<j2" (and not "iy<=j2"): it cannot occur that delphi->y[ j2 ] <= lim2 - delphi->hside
									if (delphi->y[iy] <= lim2 - delphi->hside) // delphi->y[iy] >= lim1 - delphi->hside is always true if iy>=j1
									{
										write4DVector<int>(delphi->epsmap,inside,m,iy,n,1,NX,NY,NZ,3);
									}
								}
							}
						}
					}
					threadPanelVolume[panel][thread_id] += last_vol_integral;
				}
			}
		}
	}

	// set<string> orphans;
	// if a triangulation is needed and we don't already have a scalar field,
	// more rays need to be casted
	if (accurateTriangulation && !isAvailableScalarField)
	{
		double last_vol_integral = 0.0;

		if (thread_id == 0)
			panelVolumeFlag[panel][1] = 1;

		double delta = delta_accurate_triangulation - delphi->hside;

		if (panel == 0) {
			pa[0] = delphi->x[0] + delta;
			pb[0] = delphi->x[NX-1] + delta;
		}
		else if (panel == 1) {
			pa[2] = delphi->z[0] + delta;
			pb[2] = delphi->z[NZ-1] + delta;

		}
		else {
			pa[1] = delphi->y[0] + delta;
			pb[1] = delphi->y[NY-1] + delta;
		}
		for (int nn = start; nn < end; nn += jump)
		{
			for (int n = nn; n < min(end, nn + iters_block); n++)
			{
				if (panel == 0) {
					// rays along axis X -> smallest striding and maximal caching efficiency
					pa[2] = delphi->z[n] + delta;
					pb[2] = pa[2];
				}
				else if (panel == 1) {
					// rays along axis Z -> largest striding and minimal caching efficiency
					pa[1] = delphi->y[n] + delta;
					pb[1] = pa[1];
				}
				else {
					pa[2] = delphi->z[n] + delta;
					pb[2] = pa[2];
				}

				for (int m=0; m<nb; m++)
				{
					++threadTotalRays[thread_id];

					int id = (n*N_MAX+m)*numPanels + (panels[0]+panel);
					int vsize = num_pixel_intersections[id];

					if (vsize == 0)
						continue;

					intersection_indices.clear();

					// check if repetitions allow the ray to enter and exist. If the
					// ray only enters then this is an error
					double lastValid = INFINITY;
					bool closure = false;

					#if !defined(NEW_INTERSECTION_FILTERING)
					for (int int_id = 0; int_id < vsize; int_id++)
					{
						double entry_point = pixel_intersections[id][int_id].first;

						int entry_int_id = int_id;

						if (fabs(entry_point - lastValid) >= EPS_INT)
						{
							closure = false;
							// search the exiting point of the ray
							while (1)
							{
								++int_id;
								if (int_id == vsize)
									break;

								if (fabs(pixel_intersections[id][int_id].first - entry_point) >= EPS_INT)
								{
									// ok the ray is exiting
									intersection_indices.push_back( pair<int,int>(entry_int_id,int_id) );
									lastValid = pixel_intersections[id][int_id].first;
									closure = true;
									break;
								}
							}
						}
					}

					// single or multiple tangent intersections
					if (lastValid == INFINITY)
					{
						continue;
					}

					#else // NEW_INTERSECTION_FILTERING

					// if vsize is even ..., otherwise closure remains false;
					if (((vsize>>1)<<1) == vsize)
					{
						for (int i=0; i<vsize; i+=2)
						{
							// Any two non-excessively close intersections are stored and considered later
							if (fabs(pixel_intersections[id][i+1].first - pixel_intersections[id][i].first) >= EPS_INT)
							{
								intersection_indices.push_back( pair<int,int>(i,i+1) );
							}
						}
						closure = true;
					}

					if (intersection_indices.size() == 0)
					{
						continue;
					}
					#endif // NEW_INTERSECTION_FILTERING

					// check finished.
					if (!closure)
					{
						threadPanelVolume[panel][thread_id] += last_vol_integral;

						// approximate the current ray with the previous one
						// analytical intersections cannot be recovered, only in/out position is recorded
						// The following marching cubes will be semi analytical. Semi means that were the analytical intersection
						// is not present the usual marching cubes rule will be used

						if (m > 0 && panel == 0)
						{
							for (int i=0; i<NX; i++)
							{
								#if !defined(USE_COMPRESSED_GRIDS)
								if (!optimizeGrids)
								{
									verticesInsidenessMap[n][m][i] = verticesInsidenessMap[n][m-1][i];
									// bool val = read3DVector<bool>(verticesInsidenessMap,i,m-1,n,NX,NY,NZ);
									// write3DVector<bool>(verticesInsidenessMap,val,i,m,n,NX,NY,NZ);
								}
								else
								#endif
								{
									bool val = read32xCompressedGrid(compressed_verticesInsidenessMap,i,m-1,n,NX,NY,NZ);
									write32xCompressedGrid(compressed_verticesInsidenessMap,val,i,m,n,NX,NY,NZ);
								}
							}
						}
						else if (m > 0 && panel == 1)
						{
							for (int k=0; k<NZ; k++)
							{
								#if !defined(USE_COMPRESSED_GRIDS)
								if (!optimizeGrids)
								{
									verticesInsidenessMap[k][n][m] = verticesInsidenessMap[k][n][m-1];
									// bool val = read3DVector<bool>(verticesInsidenessMap,m-1,n,k,NX,NY,NZ);
									// write3DVector<bool>(verticesInsidenessMap,val,m,n,k,NX,NY,NZ);
								}
								else
								#endif
								{
									bool val = read32xCompressedGrid(compressed_verticesInsidenessMap,m-1,n,k,NX,NY,NZ);
									write32xCompressedGrid(compressed_verticesInsidenessMap,val,m,n,k,NX,NY,NZ);
								}
							}
						}
						else if (m > 0 && panel == 2)
						{
							for (int j=0; j<NY; j++)
							{
								#if !defined(USE_COMPRESSED_GRIDS)
								if (!optimizeGrids)
								{
									verticesInsidenessMap[n][j][m] = verticesInsidenessMap[n][j][m-1];
									// bool val = read3DVector<bool>(verticesInsidenessMap,m-1,j,n,NX,NY,NZ);
									// write3DVector<bool>(verticesInsidenessMap,val,m,j,n,NX,NY,NZ);
								}
								else
								#endif
								{
									bool val = read32xCompressedGrid(compressed_verticesInsidenessMap,m-1,j,n,NX,NY,NZ);
									write32xCompressedGrid(compressed_verticesInsidenessMap,val,m,j,n,NX,NY,NZ);
								}
							}
						}

						++threadFailedRays[thread_id];

						continue;
					}

					// if there are intersections mark the epsmap cells that are inside
					// by counting the number of intersections. This is done for each panel
					// that is in the x,y,z directions respectively.
					// Eg. Tracing the ray one obtains inside/out information
					// 000000000|-1-1-1-1-1-1-1-1-1-1-1-1|000000000000|-1-1-1-1-1|000000
					if (panel == 0) {
						pa[1] = delphi->y[m] + delta;
						pb[1] = pa[1];
					}
					else if (panel == 1) {
						// rays along axis Z -> largest striding and minimal caching efficiency
						pa[0] = delphi->x[m] + delta;
						pb[0] = pa[0];
					}
					else
					{
						pa[0] = delphi->x[m] + delta;
						pb[0] = pa[0];
					}
					double lim1,lim2;

					last_vol_integral = 0.0;

					if (panel == 0)
					{
						dir = pb[0] - pa[0];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[0] + dir*pixel_intersections[id][ int_id1 ].first;
							lim2 = pa[0] + dir*pixel_intersections[id][ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							// get the cube which the intersections belong
							int xa = (int)rintp((lim1-delphi->xmin)*delphi->scale);
							int xb = (int)rintp((lim2-delphi->xmin)*delphi->scale);

							// high frequency intersection
							if (xb <= xa)
							{
								continue;
							}

							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								for (int i=xa + 1; i<=xb; i++)
									verticesInsidenessMap[n][m][i] = false;
							}
							else
							#endif
							{
								for (int i=xa + 1; i<=xb; i++)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,m,n,NX,NY,NZ);
							}

							verticesBuffers[thread_id].push_back(lim1);
							verticesBuffers[thread_id].push_back(pa[1]);
							verticesBuffers[thread_id].push_back(pa[2]);
							VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
							verticesBuffers[thread_id].push_back(lim2);
							verticesBuffers[thread_id].push_back(pa[1]);
							verticesBuffers[thread_id].push_back(pa[2]);
							VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

							#if !defined(COORD_NORM_PACKING)
							{
								coordVec cv(xa,m,n, intersec1, X_DIR);
								v_int->push_back(cv);
							}
							{
								coordVec cv(xb,m,n, intersec2, X_DIR);
								v_int->push_back(cv);
							}
							if (computeNormals && providesAnalyticalNormals)
							{
								VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);
								if (normal != NULL)
								{
									coordVec cv(xa,m,n, normal, X_DIR);
									v_norm->push_back(cv);
								}
								normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
								if (normal != NULL)
								{
									coordVec cv(xb,m,n, normal, X_DIR);
									v_norm->push_back(cv);
								}
							}
							#else // COORD_NORM_PACKING
							VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);
							VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(xa,m,n, intersec1, normal1, X_DIR);
								#else
								compressedCoordNormPacket cnv(xa,m,n, intersec1, normal1, X_DIR);
								#endif
								v_int->push_back(cnv);
							}
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(xb,m,n, intersec2, normal2, X_DIR);
								#else
								compressedCoordNormPacket cnv(xb,m,n, intersec2, normal2, X_DIR);
								#endif
								v_int->push_back(cnv);

							}
							#endif // COORD_NORM_PACKING
						}
					}
					else if (panel == 1)
					{
						dir = pb[2] - pa[2];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[2] + dir*pixel_intersections[id][ int_id1 ].first;
							lim2 = pa[2] + dir*pixel_intersections[id][ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							int za = (int)rintp((lim1-delphi->zmin)*delphi->scale);
							int zb = (int)rintp((lim2-delphi->zmin)*delphi->scale);

							// high frequency intersection
							if (zb <= za)
							{
								continue;
							}

							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								for (int k=za + 1; k<=zb; k++)
									verticesInsidenessMap[k][n][m] = false;
							}
							else
							#endif
							{
								for (int k=za + 1; k<=zb; k++)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,m,n,k,NX,NY,NZ);
							}

							verticesBuffers[thread_id].push_back(pa[0]);
							verticesBuffers[thread_id].push_back(pa[1]);
							verticesBuffers[thread_id].push_back(lim1);
							VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
							verticesBuffers[thread_id].push_back(pa[0]);
							verticesBuffers[thread_id].push_back(pa[1]);
							verticesBuffers[thread_id].push_back(lim2);
							VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

							#if !defined(COORD_NORM_PACKING)
							{
								coordVec cv(m,n,za, intersec1, Z_DIR);
								v_int->push_back(cv);
							}
							{
								coordVec cv(m,n,zb, intersec2, Z_DIR);
								v_int->push_back(cv);
							}
							if (computeNormals && providesAnalyticalNormals)
							{
								VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);
								if (normal != NULL)
								{
									coordVec cv(m,n,za, normal, Z_DIR);
									v_norm->push_back(cv);
								}
								normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
								if (normal != NULL)
								{
									coordVec cv(m,n,zb, normal, Z_DIR);
									v_norm->push_back(cv);
								}
							}
							#else // COORD_NORM_PACKING
							VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);
							VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(m,n,za, intersec1, normal1, Z_DIR);
								#else
								compressedCoordNormPacket cnv(m,n,za, intersec1, normal1, Z_DIR);
								#endif
								v_int->push_back(cnv);
							}
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(m,n,zb, intersec2, normal2, Z_DIR);
								#else
								compressedCoordNormPacket cnv(m,n,zb, intersec2, normal2, Z_DIR);
								#endif
								v_int->push_back(cnv);
							}
							#endif // COORD_NORM_PACKING
						}
					}
					else
					{
						dir = pb[1] - pa[1];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[1] + dir*pixel_intersections[id][ int_id1 ].first;
							lim2 = pa[1] + dir*pixel_intersections[id][ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							int ya = (int)rintp((lim1-delphi->ymin)*delphi->scale);
							int yb = (int)rintp((lim2-delphi->ymin)*delphi->scale);

							// high frequency intersection
							if (yb <= ya)
							{
								continue;
							}

							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								for (int j=ya + 1; j<=yb; j++)
									verticesInsidenessMap[n][j][m] = false;
							}
							else
							#endif
							{
								for (int j=ya + 1; j<=yb; j++)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,m,j,n,NX,NY,NZ);
							}

							verticesBuffers[thread_id].push_back(pa[0]);
							verticesBuffers[thread_id].push_back(lim1);
							verticesBuffers[thread_id].push_back(pa[2]);
							VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
							verticesBuffers[thread_id].push_back(pa[0]);
							verticesBuffers[thread_id].push_back(lim2);
							verticesBuffers[thread_id].push_back(pa[2]);
							VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

							#if !defined(COORD_NORM_PACKING)
							{
								coordVec cv(m,ya,n, intersec1, Y_DIR);
								v_int->push_back(cv);
							}
							{
								coordVec cv(m,yb,n, intersec2, Y_DIR);
								v_int->push_back(cv);
							}
							if (computeNormals && providesAnalyticalNormals)
							{
								VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);

								if (normal != NULL)
								{
									coordVec cv(m,ya,n, normal, Y_DIR);
									v_norm->push_back(cv);
								}
								normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
								if (normal != NULL)
								{
									coordVec cv(m,yb,n, normal, Y_DIR);
									v_norm->push_back(cv);
								}
							}
							#else // COORD_NORM_PACKING
							VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id1].second);
							VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, pixel_intersections[id][int_id2].second);
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(m,ya,n, intersec1, normal1, Y_DIR);
								#else
								compressedCoordNormPacket cnv(m,ya,n, intersec1, normal1, Y_DIR);
								#endif
								v_int->push_back(cnv);
							}
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(m,yb,n, intersec2, normal2, Y_DIR);
								#else
								compressedCoordNormPacket cnv(m,yb,n, intersec2, normal2, Y_DIR);
								#endif
								v_int->push_back(cnv);
							}
							#endif // COORD_NORM_PACKING
						}
					}
					threadPanelVolume[panel][thread_id] += last_vol_integral;
				}
			}
		}
	}
}


// Conventional ray-based ray-tracing routine; it computes the intersections and thanks to them it stores vertices,
// normals and grid data (e.g, the status map), if required
// SD PB_NEW also save the grid packets now
void Surface::intersectWithRayBasedAlgorithm (int thread_id, int nb, int start, int end, int iters_block, int jump, packet pack, packet gridPack)
{
	vector<pair<VERTEX_TYPE,VERTEX_TYPE*>> intersections;
	vector<pair<int,int>> intersection_indices;

	intersections.reserve(2000);
	intersection_indices.reserve(2000);

	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *v_int = pack.first;
	vector<coordVec> *v_norm = pack.second;
	#else
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *v_int = pack.first;
	#else
	vector<compressedCoordNormPacket> *v_int = pack.first;
	#endif
	#endif

	// SD PB_NEW
	bool enableSaveGridIntersections = false;

	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *v_int_grid = nullptr;
	vector<coordVec> *v_norm_grid = nullptr;
	#else
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *v_int_grid = nullptr;
	#else
	vector<compressedCoordNormPacket> *v_int_grid = nullptr;
	#endif
	#endif
	if  (gridPack.first != nullptr)
	{
		enableSaveGridIntersections = true;

		v_int_grid = gridPack.first;
		#if !defined(COORD_NORM_PACKING)
		v_norm_grid = gridPack.second;
		#endif
	}
	// END PB_NEW

	double pa[3],pb[3];
	double dir;

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	// if cavities and espilon map are not necessary then skip this step at all
	if (!collectFaceIntersections && !delphi->buildStatus && !delphi->buildEpsMap)
	{
		// skip
	}
	// if only cavity detection is required, only the first panel is ray cast
	else if (!collectFaceIntersections && (delphi->buildStatus && panel != 0) && !delphi->buildEpsMap)
	{
		// skip panel 1,2
	}
	// if both cavity detection and epsilon map are required, 3 panels must be analyzed
	else if (collectFaceIntersections || (delphi->buildStatus && panel == 0) || delphi->buildEpsMap)
	{
		double last_vol_integral = 0.0;

		if (thread_id == 0)
			panelVolumeFlag[panel][0] = 1;

		if (panel == 0) {
			pa[0] = delphi->x[0];
			pb[0] = delphi->x[NX-1];
		}
		else if (panel == 1) {
			pa[2] = delphi->z[0];
			pb[2] = delphi->z[NZ-1];
		}
		else
		{
			pa[1] = delphi->y[0];
			pb[1] = delphi->y[NY-1];
		}
		for (int nn = start; nn < end; nn += jump)
		{
			for (int n = nn; n < min(end, nn + iters_block); n++)
			{
				if (panel == 0)
				{
					// rays along axis X -> smallest striding and maximal caching efficiency
					pa[2] = delphi->z[n];
					pb[2] = pa[2];
				}
				else if (panel == 1)
				{
					// rays along axis Z -> largest striding and minimal caching efficiency
					pa[1] = delphi->y[n];
					pb[1] = pa[1];
				}
				else
				{
					pa[2] = delphi->z[n];
					pb[2] = pa[2];
				}

				for (int m=0; m<nb; m++)
				{
					++threadTotalRays[thread_id];

					if (panel == 0)
					{
						// rays along axis X -> smallest striding and maximal caching efficiency
						pa[1] = delphi->y[m];
						pb[1] = pa[1];
					}
					else if (panel == 1)
					{
						// rays along axis Z -> largest striding and minimal caching efficiency
						pa[0] = delphi->x[m];
						pb[0] = pa[0];
					}
					else
					{
						pa[0] = delphi->x[m];
						pb[0] = pa[0];
					}
					intersections.clear();
					getRayIntersection(pa,pb,intersections,computeNormals,thread_id);

					if (intersections.size() == 0)
						continue;

					intersection_indices.clear();

					// check if repetitions allow the ray to enter and exist. If the
					// ray only enters then this is an error
					vector<pair<VERTEX_TYPE,VERTEX_TYPE*>>::iterator itt = intersections.begin();

					int vsize = (int)intersections.size();
					double lastValid = INFINITY;
					bool closure = false;

					#if !defined(NEW_INTERSECTION_FILTERING)

					for (int int_id = 0; itt != intersections.end(); itt++, int_id++)
					{
						double entry_point = itt->first;

						int entry_int_id = int_id;
						// get the entering point
						if (fabs(entry_point - lastValid) >= EPS_INT)
						{
							// ok the ray is entering
							closure = false;
							// search the exiting point of the ray
							while (1)
							{
								++int_id;
								++itt;
								if (int_id == vsize)
									break;

								if (fabs(itt->first - entry_point) >= EPS_INT)
								{
									// ok the ray is exiting
									intersection_indices.push_back( pair<int,int>(entry_int_id,int_id) );
									lastValid = itt->first;
									closure = true;
									break;
								}
							}
						}
					}

					// single or multiple tangent intersections
					if (lastValid == INFINITY)
					{
						continue;
					}

					#else // NEW_INTERSECTION_FILTERING

					// if vsize is even ..., otherwise closure remains false;
					if (((vsize>>1)<<1) == vsize)
					{
						int entry_int_id = 0;

						for (; itt != intersections.end(); itt++)
						{
							double entry_point = itt->first;
							++itt;
							double exit_point = itt->first;
							// Any two non-excessively close intersections are stored and considered later
							if (fabs(exit_point - entry_point) >= EPS_INT)
							{
								intersection_indices.push_back( pair<int,int>(entry_int_id,entry_int_id+1) );
							}
							entry_int_id += 2;
						}
						closure = true;
					}

					if (intersection_indices.size() == 0)
					{
						continue;
					}
					#endif // NEW_INTERSECTION_FILTERING

					// check finished.
					if (!closure)
					{
						threadPanelVolume[panel][thread_id] += last_vol_integral;

						#if defined(REPORT_FAILED_RAYS)
						printRayIntersection(pa,pb);
						#endif
						// approximate the current ray with the previous one
						if (m > 0 && panel == 0)
						{
							// eps, ideb, status
							if (delphi->buildStatus && delphi->buildEpsMap)
							{
								for (int ix=0; ix<NX; ix++)
								{
									const int eee = read4DVector<int>(delphi->epsmap,ix,m-1,n,0,NX,NY,NZ,3);
									write4DVector<int>(delphi->epsmap,eee,ix,m,n,0,NX,NY,NZ,3);

									if (!optimizeGrids)
									{
										const int sss = read3DVector<int>(delphi->status,ix,m-1,n,NX,NY,NZ);
										write3DVector<int>(delphi->status,sss,ix,m,n,NX,NY,NZ);
									}
									else
									{
										const int sss = readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,m-1,n,NX,NY,NZ);
										writeBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,sss,ix,m,n,NX,NY,NZ);
									}
									const bool bbb = read3DVector<bool>(delphi->idebmap,ix,m-1,n,NX,NY,NZ);
									write3DVector<bool>(delphi->idebmap,bbb,ix,m,n,NX,NY,NZ);
								}
							}
							// eps, idebmap
							else if (delphi->buildEpsMap)
							{
								for (int ix=0; ix<NX; ix++)
								{
									const int eee = read4DVector<int>(delphi->epsmap,ix,m-1,n,0,NX,NY,NZ,3);
									write4DVector<int>(delphi->epsmap,eee,ix,m,n,0,NX,NY,NZ,3);

									const bool bbb = read3DVector<bool>(delphi->idebmap,ix,m-1,n,NX,NY,NZ);
									write3DVector<bool>(delphi->idebmap,bbb,ix,m,n,NX,NY,NZ);
								}
							}
							// status
							else if (delphi->buildStatus)
							{
								for (int ix=0; ix<NX; ix++)
								{
									if (!optimizeGrids)
									{
										const int sss = read3DVector<int>(delphi->status,ix,m-1,n,NX,NY,NZ);
										write3DVector<int>(delphi->status,sss,ix,m,n,NX,NY,NZ);
									}
									else
									{
										const int sss = readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,m-1,n,NX,NY,NZ);
										writeBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,sss,ix,m,n,NX,NY,NZ);
									}
								}
							}
						}
						else if (m > 0 && panel == 1 && delphi->buildEpsMap)
						{
							for (int iz=0; iz<NZ; iz++)
							{
								const int eee = read4DVector<int>(delphi->epsmap,m-1,n,iz,2,NX,NY,NZ,3);
								write4DVector<int>(delphi->epsmap,eee,m,n,iz,2,NX,NY,NZ,3);
							}
						}
						else if (m > 0 && panel == 2 && delphi->buildEpsMap)
						{
							for (int iy=0; iy<NY; iy++)
							{
								const int eee = read4DVector<int>(delphi->epsmap,m-1,iy,n,1,NX,NY,NZ,3);
								write4DVector<int>(delphi->epsmap,eee,m,iy,n,1,NX,NY,NZ,3);
							}
						}

						++threadFailedRays[thread_id];

						continue;
					}

					// if there are intersections mark the epsmap cells that are inside
					// by counting the number of intersections. This is done for each panel
					// that is in the x,y,z directions respectively.
					// Eg. Tracing the ray one obtains inside/out information
					// 000000000|-1-1-1-1-1-1-1-1-1-1-1-1|000000000000|-1-1-1-1-1|000000
					double lim1,lim2;
					double lastValidLim2 = INFINITY;

					last_vol_integral = 0.0;

					if (panel == 0)
					{
						dir = pb[0] - pa[0];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[0] + dir*intersections[ int_id1 ].first;
							lim2 = pa[0] + dir*intersections[ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							// get the cube which the intersections belong
							int i1 = (int)rintp((lim1-delphi->xmin)*delphi->scale);
							int i2 = (int)rintp((lim2-delphi->xmin)*delphi->scale);

							i1 = (i1 < 0)   ? 0    : i1;
							i1 = (i1 >= NX) ? NX-1 : i1;
							i2 = (i2 < 0)   ? 0    : i2;
							i2 = (i2 >= NX) ? NX-1 : i2;

							// second inters. after the first one .. it should not occur
							if (i2 < i1)
							{
								continue;
							}
							// SD PB_NEW collect also grid intersection data
							if (enableSaveGridIntersections)
							{
								verticesBuffers[thread_id].push_back(lim1);
								verticesBuffers[thread_id].push_back(pa[1]);
								verticesBuffers[thread_id].push_back(pa[2]);
								VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
								verticesBuffers[thread_id].push_back(lim2);
								verticesBuffers[thread_id].push_back(pa[1]);
								verticesBuffers[thread_id].push_back(pa[2]);
								VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

								#if !defined(COORD_NORM_PACKING)
								{
									coordVec cv(i1,m,n, intersec1, X_DIR);
									v_int_grid->push_back(cv);
								}
								{
									coordVec cv(i2,m,n, intersec2, X_DIR);
									v_int_grid->push_back(cv);
								}
								if (computeNormals && providesAnalyticalNormals)
								{
									VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);

									if (normal != NULL)
									{
										coordVec cv(i1,m,n, normal, X_DIR);
										v_norm_grid->push_back(cv);
									}
									normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
									if (normal != NULL)
									{
										coordVec cv(i2,m,n, normal, X_DIR);
										v_norm_grid->push_back(cv);
									}
								}
								#else // COORD_NORM_PACKING
								VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);
								VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(i1,m,n, intersec1, normal1, X_DIR);
									#else
									compressedCoordNormPacket cnv(i1,m,n, intersec1, normal1, X_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(i2,m,n, intersec2, normal2, X_DIR);
									#else
									compressedCoordNormPacket cnv(i2,m,n, intersec2, normal2, X_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								#endif // COORD_NORM_PACKING
							}
							// END PB_NEW

							// eps, status, idebmap
							if (delphi->buildEpsMap && delphi->buildStatus)
							{
								for (int ix=i1; ix<=i2; ix++)
								{
									// epsmap
									if (delphi->x[ix] <= lim2 - delphi->hside) // delphi->x[ix] >= lim1 - delphi->hside is always true if ix >= i1
									{
										write4DVector<int>(delphi->epsmap,inside,ix,m,n,0,NX,NY,NZ,3);
									}
									// status and idebmap
									// in the loop arguments, "ix<=i2" because, below, it could be true that delphi->x[i2] <= lim2
									if (delphi->x[ix] <= lim2 && delphi->x[ix] >= lim1)
									{
										if (!optimizeGrids)
											write3DVector<int>(delphi->status,STATUS_POINT_INSIDE,ix,m,n,NX,NY,NZ);
										else
											writeBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,STATUS_POINT_INSIDE,ix,m,n,NX,NY,NZ);
										write3DVector<bool>(delphi->idebmap,false,ix,m,n,NX,NY,NZ);
									}
								}
							}
							// epsmap, idebmap
							else if (delphi->buildEpsMap)
							{
								for (int ix=i1; ix<=i2; ix++)
								{
									if (delphi->x[ix] <= lim2 - delphi->hside) // delphi->x[ix] >= lim1-delphi->hside is always true if ix>=i1
									{
										write4DVector<int>(delphi->epsmap,inside,ix,m,n,0,NX,NY,NZ,3);
									}
									if (delphi->x[ix] <= lim2 && delphi->x[ix] >= lim1)
									{
										write3DVector<bool>(delphi->idebmap,false,ix,m,n,NX,NY,NZ);
									}
								}
							}
							// only status
							else if (delphi->buildStatus)
							{
								for (int ix=i1; ix<=i2; ix++)
								{
									if (delphi->x[ix] <= lim2 && delphi->x[ix] >= lim1)
									{
										if (!optimizeGrids)
											write3DVector<int>(delphi->status,STATUS_POINT_INSIDE,ix,m,n,NX,NY,NZ);
										else
											writeBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,STATUS_POINT_INSIDE,ix,m,n,NX,NY,NZ);
									}
								}
							}
						}
					}
					else if (panel == 1)
					{
						dir = pb[2] - pa[2];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[2] + dir*intersections[ int_id1 ].first;
							lim2 = pa[2] + dir*intersections[ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							// get the cube which the intersections belong
							int k1 = (int)rintp((lim1-delphi->zmin)*delphi->scale);
							int k2 = (int)rintp((lim2-delphi->zmin)*delphi->scale);

							k1 = (k1 < 0)   ? 0    : k1;
							k1 = (k1 >= NZ) ? NZ-1 : k1;
							k2 = (k2 < 0)   ? 0    : k2;
							k2 = (k2 >= NZ) ? NZ-1 : k2;

							// second inters. after the first one .. it should not occur
							if (k2 < k1)
							{
								continue;
							}
							// SD PB_NEW collect also grid intersection data
							if (enableSaveGridIntersections)
							{
								verticesBuffers[thread_id].push_back(pa[0]);
								verticesBuffers[thread_id].push_back(pa[1]);
								verticesBuffers[thread_id].push_back(lim1);
								VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
								verticesBuffers[thread_id].push_back(pa[0]);
								verticesBuffers[thread_id].push_back(pa[1]);
								verticesBuffers[thread_id].push_back(lim2);
								VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

								#if !defined(COORD_NORM_PACKING)
								{
									coordVec cv(m,n,k1, intersec1, Z_DIR);
									v_int_grid->push_back(cv);
								}
								{
									coordVec cv(m,n,k2, intersec2, Z_DIR);
									v_int_grid->push_back(cv);
								}
								if (computeNormals && providesAnalyticalNormals)
								{
									VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);

									if (normal != NULL)
									{
										coordVec cv(m,n,k1, normal, Z_DIR);
										v_norm_grid->push_back(cv);
									}
									normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
									if (normal != NULL)
									{
										coordVec cv(m,n,k2, normal, Z_DIR);
										v_norm_grid->push_back(cv);
									}
								}
								#else // COORD_NORM_PACKING
								VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);
								VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(m,n,k1, intersec1, normal1, Z_DIR);
									#else
									compressedCoordNormPacket cnv(m,n,k1, intersec1, normal1, Z_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(m,n,k2, intersec2, normal2, Z_DIR);
									#else
									compressedCoordNormPacket cnv(m,n,k2, intersec2, normal2, Z_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								#endif // COORD_NORM_PACKING
							}
							// END PB_NEW

							if (delphi->buildEpsMap)
							{
								for (int iz=k1; iz<k2; iz++)
								{
									// "iz<i2" (and not "iz<=i2"): it cannot occur that delphi->z[ i2 ] <= lim2-delphi->hside
									if (delphi->z[iz] <= lim2 - delphi->hside) // delphi->z[iz] >= lim1-delphi->hside is always true if iz>=i1
									{
										write4DVector<int>(delphi->epsmap,inside,m,n,iz,2,NX,NY,NZ,3);
									}
								}
							}
						}
					}
					else
					{
						dir = pb[1] - pa[1];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[1] + dir*intersections[ int_id1 ].first;
							lim2 = pa[1] + dir*intersections[ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							// get the cube which the intersections belong
							int j1 = (int)rintp((lim1-delphi->ymin)*delphi->scale);
							int j2 = (int)rintp((lim2-delphi->ymin)*delphi->scale);

							j1 = (j1 < 0)   ? 0    : j1;
							j1 = (j1 >= NY) ? NY-1 : j1;
							j2 = (j2 < 0)   ? 0    : j2;
							j2 = (j2 >= NY) ? NY-1 : j2;

							// second inters. after the first one .. it should not occur
							if (j2 < j1)
							{
								continue;
							}
							// SD PB_NEW collect also grid intersection data
							if (enableSaveGridIntersections)
							{
								verticesBuffers[thread_id].push_back(pa[0]);
								verticesBuffers[thread_id].push_back(lim1);
								verticesBuffers[thread_id].push_back(pa[2]);
								VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
								verticesBuffers[thread_id].push_back(pa[0]);
								verticesBuffers[thread_id].push_back(lim2);
								verticesBuffers[thread_id].push_back(pa[2]);
								VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

								#if !defined(COORD_NORM_PACKING)
								{
									coordVec cv(m,j1,n, intersec1, Y_DIR);
									v_int_grid->push_back(cv);
								}
								{
									coordVec cv(m,j2,n, intersec2, Y_DIR);
									v_int_grid->push_back(cv);
								}
								if (computeNormals && providesAnalyticalNormals)
								{
									VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);

									if (normal != NULL)
									{
										coordVec cv(m,j1,n, normal, Y_DIR);
										v_norm_grid->push_back(cv);
									}
									normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
									if (normal != NULL)
									{
										coordVec cv(m,j2,n, normal, Y_DIR);
										v_norm_grid->push_back(cv);
									}
								}
								#else // COORD_NORM_PACKING
								VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);
								VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(m,j1,n, intersec1, normal1, Y_DIR);
									#else
									compressedCoordNormPacket cnv(m,j1,n, intersec1, normal1, Y_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								{
									#if !defined(COMPRESS_INTERSECTION_COORDS)
									coordNormPacket cnv(m,j2,n, intersec2, normal2, Y_DIR);
									#else
									compressedCoordNormPacket cnv(m,j2,n, intersec2, normal2, Y_DIR);
									#endif
									v_int_grid->push_back(cnv);
								}
								#endif // COORD_NORM_PACKING
							}
							// END PB_NEW

							if (delphi->buildEpsMap)
							{
								// for (int iy=0; iy<NY; iy++)
								for (int iy=j1; iy<j2; iy++)
								{
									// "iy<j2" (and not "iy<=j2"): it cannot occur that delphi->y[ j2 ] <= lim2-delphi->hside
									if (delphi->y[iy] <= lim2 - delphi->hside) // delphi->y[iy] >= lim1 - delphi->hside is always true if iy>=j1
									{
										write4DVector<int>(delphi->epsmap,inside,m,iy,n,1,NX,NY,NZ,3);
									}
								}
							}
						}
					}
					threadPanelVolume[panel][thread_id] += last_vol_integral;
				}
			}
		}
	}

	// set<string> orphans;
	// if a triangulation is needed and we don't already have a scalar field,
	// more rays need to be casted
	if (accurateTriangulation && !isAvailableScalarField)
	{
		double last_vol_integral = 0.0;

		if (thread_id == 0)
			panelVolumeFlag[panel][1] = 1;

		double delta = delta_accurate_triangulation - delphi->hside;

		if (panel == 0) {
			pa[0] = delphi->x[0] + delta;
			pb[0] = delphi->x[NX-1] + delta;
		}
		else if (panel == 1) {
			pa[2] = delphi->z[0] + delta;
			pb[2] = delphi->z[NZ-1] + delta;
		}
		else {
			pa[1] = delphi->y[0] + delta;
			pb[1] = delphi->y[NY-1] + delta;
		}
		for (int nn = start; nn < end; nn += jump)
		{
			for (int n = nn; n < min(end, nn + iters_block); n++)
			{
				if (panel == 0) {
					// rays along axis X -> smallest striding and maximal caching efficiency
					pa[2] = delphi->z[n] + delta;
					pb[2] = pa[2];
				}
				else if (panel == 1) {
					// rays along axis Z -> largest striding and minimal caching efficiency
					pa[1] = delphi->y[n] + delta;
					pb[1] = pa[1];
				}
				else {
					pa[2] = delphi->z[n] + delta;
					pb[2] = pa[2];
				}

				for (int m=0; m<nb; m++)
				{
					++threadTotalRays[thread_id];

					if (panel == 0) {
						pa[1] = delphi->y[m] + delta;
						pb[1] = pa[1];
					}
					else if (panel == 1) {
						// rays along axis Z -> largest striding and minimal caching efficiency
						pa[0] = delphi->x[m] + delta;
						pb[0] = pa[0];
					}
					else
					{
						pa[0] = delphi->x[m] + delta;
						pb[0] = pa[0];
					}
					intersections.clear();
					getRayIntersection(pa,pb,intersections,computeNormals,thread_id);

					if (intersections.size() == 0)
						continue;

					intersection_indices.clear();

					// check if repetitions allow the ray to enter and exist. If the
					// ray only enters then this is an error
					vector<pair<VERTEX_TYPE,VERTEX_TYPE*>>::iterator itt = intersections.begin();

					int vsize = (int)intersections.size();
					double lastValid = INFINITY;
					bool closure = false;

					#if !defined(NEW_INTERSECTION_FILTERING)
					for (int int_id = 0; itt != intersections.end(); itt++, int_id++)
					{
						double entry_point = itt->first;

						int entry_int_id = int_id;
						// get the entering point
						if (fabs(entry_point - lastValid) >= EPS_INT)
						{
							// ok the ray is entering
							closure = false;
							// search the exiting point of the ray
							while (1)
							{
								++int_id;
								++itt;
								if (int_id == vsize)
									break;

								if (fabs(itt->first - entry_point) >= EPS_INT)
								{
									// ok the ray is exiting
									intersection_indices.push_back( pair<int,int>(entry_int_id,int_id) );
									lastValid = itt->first;
									closure = true;
									break;
								}
							}
						}
					}

					// single or multiple tangent intersections
					if (lastValid == INFINITY)
					{
						continue;
					}

					#else // NEW_INTERSECTION_FILTERING

					// if vsize is even ..., otherwise closure remains false;
					if (((vsize>>1)<<1) == vsize)
					{
						int entry_int_id = 0;

						for (; itt != intersections.end(); itt++)
						{
							double entry_point = itt->first;
							++itt;
							double exit_point = itt->first;
							// Any two non-excessively close intersections are stored and considered later
							if (fabs(exit_point - entry_point) >= EPS_INT)
							{
								intersection_indices.push_back( pair<int,int>(entry_int_id,entry_int_id+1) );
							}
							entry_int_id += 2;
						}
						closure = true;
					}

					if (intersection_indices.size() == 0)
					{
						continue;
					}
					#endif // NEW_INTERSECTION_FILTERING

					// check finished.
					if (!closure)
					{
						threadPanelVolume[panel][thread_id] += last_vol_integral;

						#if defined(REPORT_FAILED_RAYS)
						printRayIntersection(pa,pb);
						#endif
						// approximate the current ray with the previous one
						// analytical intersections cannot be recovered, only in/out position is recorded
						// The following marching cubes will be semi analytical. Semi means that were the analytical intersection
						// is not present the usual marching cubes rule will be used

						if (m > 0 && panel == 0)
						{
							for (int i=0; i<NX; i++)
							{
								#if !defined(USE_COMPRESSED_GRIDS)
								if (!optimizeGrids)
								{
									verticesInsidenessMap[n][m][i] = verticesInsidenessMap[n][m-1][i];
									// bool val = read3DVector<bool>(verticesInsidenessMap,i,m-1,n,NX,NY,NZ);
									// write3DVector<bool>(verticesInsidenessMap,val,i,m,n,NX,NY,NZ);
								}
								else
								#endif
								{
									bool val = read32xCompressedGrid(compressed_verticesInsidenessMap,i,m-1,n,NX,NY,NZ);
									write32xCompressedGrid(compressed_verticesInsidenessMap,val,i,m,n,NX,NY,NZ);
								}
							}
						}
						else if (m > 0 && panel == 1)
						{
							for (int k=0; k<NZ; k++)
							{
								#if !defined(USE_COMPRESSED_GRIDS)
								if (!optimizeGrids)
								{
									verticesInsidenessMap[k][n][m] = verticesInsidenessMap[k][n][m-1];
									// bool val = read3DVector<bool>(verticesInsidenessMap,m-1,n,k,NX,NY,NZ);
									// write3DVector<bool>(verticesInsidenessMap,val,m,n,k,NX,NY,NZ);
								}
								else
								#endif
								{
									bool val = read32xCompressedGrid(compressed_verticesInsidenessMap,m-1,n,k,NX,NY,NZ);
									write32xCompressedGrid(compressed_verticesInsidenessMap,val,m,n,k,NX,NY,NZ);
								}
							}
						}
						else if (m > 0 && panel == 2)
						{
							for (int j=0; j<NY; j++)
							{
								#if !defined(USE_COMPRESSED_GRIDS)
								if (!optimizeGrids)
								{
									verticesInsidenessMap[n][j][m] = verticesInsidenessMap[n][j][m-1];
									// bool val = read3DVector<bool>(verticesInsidenessMap,m-1,j,n,NX,NY,NZ);
									// write3DVector<bool>(verticesInsidenessMap,val,m,j,n,NX,NY,NZ);
								}
								else
								#endif
								{
									bool val = read32xCompressedGrid(compressed_verticesInsidenessMap,m-1,j,n,NX,NY,NZ);
									write32xCompressedGrid(compressed_verticesInsidenessMap,val,m,j,n,NX,NY,NZ);
								}
							}
						}

						++threadFailedRays[thread_id];

						continue;
					}

					// if there are intersections mark the epsmap cells that are inside
					// by counting the number of intersections. This is done for each panel
					// that is in the x,y,z directions respectively.
					// Eg. Tracing the ray one obtains inside/out information
					// 000000000|-1-1-1-1-1-1-1-1-1-1-1-1|000000000000|-1-1-1-1-1|000000
					double lim1,lim2;
					double lastValidLim2 = INFINITY;

					last_vol_integral = 0.0;

					if (panel == 0)
					{
						dir = pb[0] - pa[0];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[0] + dir*intersections[ int_id1 ].first;
							lim2 = pa[0] + dir*intersections[ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							// get the cube which the intersections belong
							int xa = (int)rintp((lim1-delphi->xmin)*delphi->scale);
							int xb = (int)rintp((lim2-delphi->xmin)*delphi->scale);

							// high frequency intersection
							if (xb <= xa)
							{
								continue;
							}

							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								for (int i=xa + 1; i<=xb; i++)
									verticesInsidenessMap[n][m][i] = false;
							}
							else
							#endif
							{
								for (int i=xa + 1; i<=xb; i++)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,m,n,NX,NY,NZ);
							}

							verticesBuffers[thread_id].push_back(lim1);
							verticesBuffers[thread_id].push_back(pa[1]);
							verticesBuffers[thread_id].push_back(pa[2]);
							VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
							verticesBuffers[thread_id].push_back(lim2);
							verticesBuffers[thread_id].push_back(pa[1]);
							verticesBuffers[thread_id].push_back(pa[2]);
							VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

							#if !defined(COORD_NORM_PACKING)
							{
								coordVec cv(xa,m,n, intersec1, X_DIR);
								v_int->push_back(cv);
							}
							{
								coordVec cv(xb,m,n, intersec2, X_DIR);
								v_int->push_back(cv);
							}
							if (computeNormals && providesAnalyticalNormals)
							{
								VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);

								if (normal != NULL)
								{
									coordVec cv(xa,m,n, normal, X_DIR);
									v_norm->push_back(cv);
								}
								normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
								if (normal != NULL)
								{
									coordVec cv(xb,m,n, normal, X_DIR);
									v_norm->push_back(cv);
								}
							}
							#else // COORD_NORM_PACKING
							VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);
							VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(xa,m,n, intersec1, normal1, X_DIR);
								#else
								compressedCoordNormPacket cnv(xa,m,n, intersec1, normal1, X_DIR);
								#endif
								v_int->push_back(cnv);
							}
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(xb,m,n, intersec2, normal2, X_DIR);
								#else
								compressedCoordNormPacket cnv(xb,m,n, intersec2, normal2, X_DIR);
								#endif
								v_int->push_back(cnv);
							}
							#endif // COORD_NORM_PACKING
						}
					}
					else if (panel == 1)
					{
						dir = pb[2] - pa[2];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[2] + dir*intersections[ int_id1 ].first;
							lim2 = pa[2] + dir*intersections[ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							int za = (int)rintp((lim1-delphi->zmin)*delphi->scale);
							int zb = (int)rintp((lim2-delphi->zmin)*delphi->scale);

							// high frequency intersection
							if (zb <= za)
							{
								continue;
							}

							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								for (int k=za + 1; k<=zb; k++)
									verticesInsidenessMap[k][n][m] = false;
							}
							else
							#endif
							{
								for (int k=za + 1; k<=zb; k++)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,m,n,k,NX,NY,NZ);
							}

							verticesBuffers[thread_id].push_back(pa[0]);
							verticesBuffers[thread_id].push_back(pa[1]);
							verticesBuffers[thread_id].push_back(lim1);
							VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
							verticesBuffers[thread_id].push_back(pa[0]);
							verticesBuffers[thread_id].push_back(pa[1]);
							verticesBuffers[thread_id].push_back(lim2);
							VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

							#if !defined(COORD_NORM_PACKING)
							{
								coordVec cv(m,n,za, intersec1, Z_DIR);
								v_int->push_back(cv);
							}
							{
								coordVec cv(m,n,zb, intersec2, Z_DIR);
								v_int->push_back(cv);
							}
							if (computeNormals && providesAnalyticalNormals)
							{
								VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);

								if (normal != NULL)
								{
									coordVec cv(m,n,za, normal, Z_DIR);
									v_norm->push_back(cv);
								}
								normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
								if (normal != NULL)
								{
									coordVec cv(m,n,zb, normal, Z_DIR);
									v_norm->push_back(cv);
								}
							}
							#else // COORD_NORM_PACKING
							VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);
							VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(m,n,za, intersec1, normal1, Z_DIR);
								#else
								compressedCoordNormPacket cnv(m,n,za, intersec1, normal1, Z_DIR);
								#endif
								v_int->push_back(cnv);
							}
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(m,n,zb, intersec2, normal2, Z_DIR);
								#else
								compressedCoordNormPacket cnv(m,n,zb, intersec2, normal2, Z_DIR);
								#endif
								v_int->push_back(cnv);
							}
							#endif // COORD_NORM_PACKING
						}
					}
					else
					{
						dir = pb[1] - pa[1];

						for (int intersec_pair = 0; intersec_pair < intersection_indices.size(); intersec_pair++)
						{
							int int_id1 = intersection_indices[ intersec_pair ].first;
							int int_id2 = intersection_indices[ intersec_pair ].second;

							lim1 = pa[1] + dir*intersections[ int_id1 ].first;
							lim2 = pa[1] + dir*intersections[ int_id2 ].first;
							last_vol_integral += lim2-lim1;

							int ya = (int)rintp((lim1-delphi->ymin)*delphi->scale);
							int yb = (int)rintp((lim2-delphi->ymin)*delphi->scale);

							// high frequency intersection
							if (yb <= ya)
							{
								continue;
							}

							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								for (int j=ya + 1; j<=yb; j++)
									verticesInsidenessMap[n][j][m] = false;
							}
							else
							#endif
							{
								for (int j=ya + 1; j<=yb; j++)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,m,j,n,NX,NY,NZ);
							}

							verticesBuffers[thread_id].push_back(pa[0]);
							verticesBuffers[thread_id].push_back(lim1);
							verticesBuffers[thread_id].push_back(pa[2]);
							VERTEX_TYPE *intersec1 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];
							verticesBuffers[thread_id].push_back(pa[0]);
							verticesBuffers[thread_id].push_back(lim2);
							verticesBuffers[thread_id].push_back(pa[2]);
							VERTEX_TYPE *intersec2 = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

							#if !defined(COORD_NORM_PACKING)
							{
								coordVec cv(m,ya,n, intersec1, Y_DIR);
								v_int->push_back(cv);
							}
							{
								coordVec cv(m,yb,n, intersec2, Y_DIR);
								v_int->push_back(cv);
							}
							if (computeNormals && providesAnalyticalNormals)
							{
								VERTEX_TYPE *normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);

								if (normal != NULL)
								{
									coordVec cv(m,ya,n, normal, Y_DIR);
									v_norm->push_back(cv);
								}
								normal = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
								if (normal != NULL)
								{
									coordVec cv(m,yb,n, normal, Y_DIR);
									v_norm->push_back(cv);
								}
							}
							#else // COORD_NORM_PACKING
							VERTEX_TYPE *normal1 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id1].second);
							VERTEX_TYPE *normal2 = sanitizeNormalPtrForThread(normalsBuffers, thread_id, computeNormals, providesAnalyticalNormals, intersections[int_id2].second);
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(m,ya,n, intersec1, normal1, Y_DIR);
								#else
								compressedCoordNormPacket cnv(m,ya,n, intersec1, normal1, Y_DIR);
								#endif
								v_int->push_back(cnv);
							}
							{
								#if !defined(COMPRESS_INTERSECTION_COORDS)
								coordNormPacket cnv(m,yb,n, intersec2, normal2, Y_DIR);
								#else
								compressedCoordNormPacket cnv(m,yb,n, intersec2, normal2, Y_DIR);
								#endif
								v_int->push_back(cnv);
							}
							#endif // COORD_NORM_PACKING
						}
					}
					threadPanelVolume[panel][thread_id] += last_vol_integral;
				}
			}
		}
	}
}


void Surface::assembleVerticesList (packet pack, int *vertex_index)
{
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *v1 = pack.first;
	vector<coordVec> *v2 = pack.second;

	vector<coordVec>::iterator it1;
	vector<coordVec>::iterator it2;
	#else
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *v1 = pack.first;

	vector<coordNormPacket>::iterator it1;
	#else
	vector<compressedCoordNormPacket> *v1 = pack.first;

	vector<compressedCoordNormPacket>::iterator it1;
	#endif
	#endif

	VERTEX_TYPE *intersec, *normal;

	int ix_, iy_, iz_, dir;

	#if !defined(COORD_NORM_PACKING)
	if (computeNormals && providesAnalyticalNormals)
		it2 = v2->begin();
	#endif

	// unload intersections/normals buffer
	for (it1 = v1->begin(); it1 != v1->end(); it1++)
	{
		#if !defined(COORD_NORM_PACKING)
		const coordVec &cv = *it1;
		intersec = cv.vec;
		ix_ = cv.ix;
		iy_ = cv.iy;
		iz_ = cv.iz;
		dir = cv.dir;

		if (computeNormals && providesAnalyticalNormals)
		{
			const coordVec &cv2 = *(it2++);
			normal = cv2.vec;
		}
		#else
		#if !defined(COMPRESS_INTERSECTION_COORDS)
		coordNormPacket &cnv = *it1;
		ix_ = cnv.ix;
		iy_ = cnv.iy;
		iz_ = cnv.iz;
		dir = cnv.dir;
		intersec = cnv.vec;
		normal = cnv.nor;
		#else
		compressedCoordNormPacket &cnv = *it1;
		cnv.getCompressedCoords(&ix_, &iy_, &iz_, &dir);
		intersec = cnv.vec;
		normal = cnv.nor;
		#endif
		#endif

		// check if it is a dangling vertex. This can happen
		// if the vertex is sampled near a cube vertex
		// in this case a local inconsistency can happen.
		// this can be easily removed by removing the vertex
		// and letting win the new vertex.
		// That is the last traced panel wins.
		// If you change the order of ray-tracing panels, the dangling
		// vertices would be other ones; but still this procedure would
		// eliminate them.
		// Note that a dangling vertex is almost identical in position to
		// its non-dangling twin, thus the degree of approximations
		// is absolutely negligible

		bool remove = 0;

		if (dir == X_DIR)
		{
			// if (read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_,NX,NY,NZ) == read3DVector<bool>(verticesInsidenessMap,ix_+1,iy_,iz_,NX,NY,NZ))
			#if !defined(USE_COMPRESSED_GRIDS)
			if (!optimizeGrids)
			{
				if (verticesInsidenessMap[iz_][iy_][ix_] == verticesInsidenessMap[iz_][iy_][ix_+1])
					remove = 1;
			}
			else
				#endif
			{
				if (read32xCompressedGrid(compressed_verticesInsidenessMap,ix_  ,iy_,iz_,NX,NY,NZ) ==
					read32xCompressedGrid(compressed_verticesInsidenessMap,ix_+1,iy_,iz_,NX,NY,NZ))
					remove = 1;
			}
		}
		else if (dir == Y_DIR)
		{
			// if (read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_,NX,NY,NZ) == read3DVector<bool>(verticesInsidenessMap,ix_,iy_+1,iz_,NX,NY,NZ))
			#if !defined(USE_COMPRESSED_GRIDS)
			if (!optimizeGrids)
			{
				if (verticesInsidenessMap[iz_][iy_][ix_] == verticesInsidenessMap[iz_][iy_+1][ix_])
					remove = 1;
			}
			else
				#endif
			{
				if (read32xCompressedGrid(compressed_verticesInsidenessMap,ix_,iy_  ,iz_,NX,NY,NZ) ==
					read32xCompressedGrid(compressed_verticesInsidenessMap,ix_,iy_+1,iz_,NX,NY,NZ))
					remove = 1;
			}
		}
		else
		{
			// if (read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_,NX,NY,NZ) == read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_+1,NX,NY,NZ))
			#if !defined(USE_COMPRESSED_GRIDS)
			if (!optimizeGrids)
			{
				if (verticesInsidenessMap[iz_][iy_][ix_] == verticesInsidenessMap[iz_+1][iy_][ix_])
					remove = 1;
			}
			else
				#endif
			{
				if (read32xCompressedGrid(compressed_verticesInsidenessMap,ix_,iy_,iz_  ,NX,NY,NZ) ==
					read32xCompressedGrid(compressed_verticesInsidenessMap,ix_,iy_,iz_+1,NX,NY,NZ))
					remove = 1;
			}
		}
		if (remove)
		{
			continue;
		}

		int current_index = -1;

		// ok the vertex is not dangling
		// add to octrees/bilevel grids and remove duplicate vertices on the same edge, if any
		// TODO now duplicated vertices should be absent, check.
		#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)

		// here, octrees are employed for storing and checking vertices' indices

		if (dir == X_DIR)
			current_index = intersectionsMatrixAlongX->at(ix_,iy_,iz_);
		else if (dir == Y_DIR)
			current_index = intersectionsMatrixAlongY->at(ix_,iy_,iz_);
		else if (dir == Z_DIR)
			current_index = intersectionsMatrixAlongZ->at(ix_,iy_,iz_);
		else
		{
			cout << endl << ERR << "Non existing direction during assembling!";
			exit(-1);
		}

		if (current_index == -1)
		{
			if (dir == X_DIR)
				intersectionsMatrixAlongX->set(ix_,iy_,iz_, *vertex_index);
			else if (dir == Y_DIR)
				intersectionsMatrixAlongY->set(ix_,iy_,iz_, *vertex_index);
			else
				intersectionsMatrixAlongZ->set(ix_,iy_,iz_, *vertex_index);

			#if !defined(AVOID_NORMALS_MATRICES)
			if (computeNormals && providesAnalyticalNormals)
			{
				if (dir == X_DIR)
					normalsMatrixAlongX->set(ix_,iy_,iz_, *vertex_index);
				else if (dir == Y_DIR)
					normalsMatrixAlongY->set(ix_,iy_,iz_, *vertex_index);
				else
					normalsMatrixAlongZ->set(ix_,iy_,iz_, *vertex_index);
			}
			#endif

			++*vertex_index;

			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			vertList.push_back(intersec);
			#else
			vertList.push_back(intersec[0]);
			vertList.push_back(intersec[1]);
			vertList.push_back(intersec[2]);
			#endif

			if (computeNormals && providesAnalyticalNormals)
			{
				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				localNormals.push_back(normal);
				#else
				localNormals.push_back(normal[0]);
				localNormals.push_back(normal[1]);
				localNormals.push_back(normal[2]);
				#endif
			}
		}

		#else // OPTIMIZE_INTERSECTIONS_MANAGEMENT

		// here, bilevel grids are employed for storing and checking vertices' indices

		if (dir == X_DIR)
			current_index = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix_,iy_,iz_,NX,NY,NZ);
		else if (dir == Y_DIR)
			current_index = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix_,iy_,iz_,NX,NY,NZ);
		else if (dir == Z_DIR)
			current_index = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix_,iy_,iz_,NX,NY,NZ);
		else
		{
			cout << endl << ERR << "Non existing direction during assembling!";
			exit(-1);
		}

		if (current_index == -1)
		{
			if (dir == X_DIR)
				writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,*vertex_index,ix_,iy_,iz_,NX,NY,NZ);
			else if (dir == Y_DIR)
				writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,*vertex_index,ix_,iy_,iz_,NX,NY,NZ);
			else
				writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,*vertex_index,ix_,iy_,iz_,NX,NY,NZ);

			++*vertex_index;

			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			vertList.push_back(intersec);
			#else
			vertList.push_back(intersec[0]);
			vertList.push_back(intersec[1]);
			vertList.push_back(intersec[2]);
			#endif

			if (computeNormals && providesAnalyticalNormals)
			{
				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				normalsList.push_back(normal);
				#else
				normalsList.push_back(normal[0]);
				normalsList.push_back(normal[1]);
				normalsList.push_back(normal[2]);
				#endif
			}
		}

		#endif // OPTIMIZE_INTERSECTIONS_MANAGEMENT
	}
}


void Surface::assembleVerticesList (packet pack, vector<VERTEX_TYPE*> *localVert, vector<VERTEX_TYPE*> *localNormals,
									int *localIndex)
{
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *v1 = pack.first;
	vector<coordVec> *v2 = pack.second;

	vector<coordVec>::iterator it1;
	vector<coordVec>::iterator it2;
	#else
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *v1 = pack.first;

	vector<coordNormPacket>::iterator it1;
	#else
	vector<compressedCoordNormPacket> *v1 = pack.first;

	vector<compressedCoordNormPacket>::iterator it1;
	#endif
	#endif

	VERTEX_TYPE *intersec, *normal;

	int ix_, iy_, iz_, dir;

	#if !defined(COORD_NORM_PACKING)
	if (computeNormals && providesAnalyticalNormals)
		it2 = v2->begin();
	#endif

	// unload intersections/normals buffer
	for (it1 = v1->begin(); it1 != v1->end(); it1++)
	{
		#if !defined(COORD_NORM_PACKING)
		const coordVec &cv = *it1;
		intersec = cv.vec;
		ix_ = cv.ix;
		iy_ = cv.iy;
		iz_ = cv.iz;
		dir = cv.dir;

		if (computeNormals && providesAnalyticalNormals)
		{
			const coordVec &cv2 = *(it2++);
			normal = cv2.vec;
		}
		#else
		#if !defined(COMPRESS_INTERSECTION_COORDS)
		coordNormPacket &cnv = *it1;
		ix_ = cnv.ix;
		iy_ = cnv.iy;
		iz_ = cnv.iz;
		dir = cnv.dir;
		intersec = cnv.vec;
		normal = cnv.nor;
		#else
		compressedCoordNormPacket &cnv = *it1;
		cnv.getCompressedCoords(&ix_, &iy_, &iz_, &dir);
		intersec = cnv.vec;
		normal = cnv.nor;
		#endif
		#endif

		// check if it is a dangling vertex. This can happen
		// if the vertex is sampled near a cube vertex
		// in this case a local inconsistency can happen.
		// this can be easily removed by removing the vertex
		// and letting win the new vertex.
		// That is the last traced panel wins.
		// If you change the order of ray-tracing panels, the dangling
		// vertices would be other ones; but still this procedure would
		// eliminate them.
		// Note that a dangling vertex is almost identical in position to
		// its non-dangling twin, thus the degree of approximations
		// is absolutely negligible

		bool remove = 0;

		if (dir == X_DIR)
		{
			// if (read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_,NX,NY,NZ) == read3DVector<bool>(verticesInsidenessMap,ix_+1,iy_,iz_,NX,NY,NZ))
			#if !defined(USE_COMPRESSED_GRIDS)
			if (!optimizeGrids)
			{
				if (verticesInsidenessMap[iz_][iy_][ix_] == verticesInsidenessMap[iz_][iy_][ix_+1])
					remove = 1;
			}
			else
			#endif
			{
				if (read32xCompressedGrid(compressed_verticesInsidenessMap,ix_  ,iy_,iz_,NX,NY,NZ) ==
					read32xCompressedGrid(compressed_verticesInsidenessMap,ix_+1,iy_,iz_,NX,NY,NZ))
					remove = 1;
			}
		}
		else if (dir == Y_DIR)
		{
			// if (read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_,NX,NY,NZ) == read3DVector<bool>(verticesInsidenessMap,ix_,iy_+1,iz_,NX,NY,NZ))
			#if !defined(USE_COMPRESSED_GRIDS)
			if (!optimizeGrids)
			{
				if (verticesInsidenessMap[iz_][iy_][ix_] == verticesInsidenessMap[iz_][iy_+1][ix_])
					remove = 1;
			}
			else
			#endif
			{
				if (read32xCompressedGrid(compressed_verticesInsidenessMap,ix_,iy_  ,iz_,NX,NY,NZ) ==
					read32xCompressedGrid(compressed_verticesInsidenessMap,ix_,iy_+1,iz_,NX,NY,NZ))
					remove = 1;
			}
		}
		else
		{
			// if (read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_,NX,NY,NZ) == read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_+1,NX,NY,NZ))
			#if !defined(USE_COMPRESSED_GRIDS)
			if (!optimizeGrids)
			{
				if (verticesInsidenessMap[iz_][iy_][ix_] == verticesInsidenessMap[iz_+1][iy_][ix_])
					remove = 1;
			}
			else
			#endif
			{
				if (read32xCompressedGrid(compressed_verticesInsidenessMap,ix_,iy_,iz_  ,NX,NY,NZ) ==
					read32xCompressedGrid(compressed_verticesInsidenessMap,ix_,iy_,iz_+1,NX,NY,NZ))
					remove = 1;
			}
		}
		if (remove)
		{
			continue;
		}

		int current_index = -1;

		// ok the vertex is not dangling
		// add to octrees/bilevel grids and remove duplicate vertices on the same edge, if any
		// TODO now duplicated vertices should be absent, check.
		#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)

		// here, octrees are employed for storing and checking vertices' indices

		if (dir == X_DIR)
			current_index = intersectionsMatrixAlongX->at(ix_,iy_,iz_);
		else if (dir == Y_DIR)
			current_index = intersectionsMatrixAlongY->at(ix_,iy_,iz_);
		else if (dir == Z_DIR)
			current_index = intersectionsMatrixAlongZ->at(ix_,iy_,iz_);
		else
		{
			cout << endl << ERR << "Non existing direction during assembling!";
			exit(-1);
		}

		if (current_index == -1)
		{
			if (dir == X_DIR)
				intersectionsMatrixAlongX->set(ix_,iy_,iz_, *localIndex);
			else if (dir == Y_DIR)
				intersectionsMatrixAlongY->set(ix_,iy_,iz_, *localIndex);
			else
				intersectionsMatrixAlongZ->set(ix_,iy_,iz_, *localIndex);

			#if !defined(AVOID_NORMALS_MATRICES)
			if (computeNormals && providesAnalyticalNormals)
			{
				if (dir == X_DIR)
					normalsMatrixAlongX->set(ix_,iy_,iz_, *localIndex);
				else if (dir == Y_DIR)
					normalsMatrixAlongY->set(ix_,iy_,iz_, *localIndex);
				else
					normalsMatrixAlongZ->set(ix_,iy_,iz_, *localIndex);
			}
			#endif

			++*localIndex;

			localVert->push_back(intersec);

			if (computeNormals && providesAnalyticalNormals)
			{
				localNormals->push_back(normal);
			}
		}

		#else // OPTIMIZE_INTERSECTIONS_MANAGEMENT

		// here, bilevel grids are employed for storing and checking vertices' indices

		if (dir == X_DIR)
			current_index = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix_,iy_,iz_,NX,NY,NZ);
		else if (dir == Y_DIR)
			current_index = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix_,iy_,iz_,NX,NY,NZ);
		else if (dir == Z_DIR)
			current_index = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix_,iy_,iz_,NX,NY,NZ);
		else
		{
			cout << endl << ERR << "Non existing direction during assembling!";
			exit(-1);
		}

		if (current_index == -1)
		{
			if (dir == X_DIR)
				writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,*localIndex,ix_,iy_,iz_,NX,NY,NZ);
			else if (dir == Y_DIR)
				writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,*localIndex,ix_,iy_,iz_,NX,NY,NZ);
			else
				writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,*localIndex,ix_,iy_,iz_,NX,NY,NZ);

			++*localIndex;

			localVert->push_back(intersec);

			if (computeNormals && providesAnalyticalNormals)
			{
				localNormals->push_back(normal);
			}
		}

		#endif // OPTIMIZE_INTERSECTIONS_MANAGEMENT
	}
}


void Surface::convertLocalGridIndicesToGlobalIndices (packet pack, int indexOffset)
{
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	#if !defined(COORD_NORM_PACKING)
	vector<coordVec> *v1 = pack.first;

	vector<coordVec>::iterator it1;
	#else
	#if !defined(COMPRESS_INTERSECTION_COORDS)
	vector<coordNormPacket> *v1 = pack.first;

	vector<coordNormPacket>::iterator it1;
	#else
	vector<compressedCoordNormPacket> *v1 = pack.first;

	vector<compressedCoordNormPacket>::iterator it1;
	#endif
	#endif

	VERTEX_TYPE *intersec;

	int ix_, iy_, iz_, dir;
	int local_index;

	for (it1 = v1->begin(); it1 != v1->end(); it1++)
	{
		#if !defined(COORD_NORM_PACKING)
		const coordVec &cv = *it1;
		intersec = cv.vec;
		#else
		#if !defined(COMPRESS_INTERSECTION_COORDS)
		coordNormPacket &cnv = *it1;
		intersec = cnv.vec;
		#else
		compressedCoordNormPacket &cnv = *it1;
		intersec = cnv.vec;
		#endif
		#endif

		if (intersec == NULL)
			continue;

		#if !defined(COORD_NORM_PACKING)
		ix_ = cv.ix;
		iy_ = cv.iy;
		iz_ = cv.iz;
		dir = cv.dir;
		#else
		#if !defined(COMPRESS_INTERSECTION_COORDS)
		ix_ = cnv.ix;
		iy_ = cnv.iy;
		iz_ = cnv.iz;
		dir = cnv.dir;
		#else
		cnv.getCompressedCoords(&ix_, &iy_, &iz_, &dir);
		#endif
		#endif

		#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)

		if (dir == X_DIR)
		{
			local_index = intersectionsMatrixAlongX->at(ix_,iy_,iz_);
			intersectionsMatrixAlongX->set(ix_,iy_,iz_, local_index + indexOffset);
		}
		else if (dir == Y_DIR)
		{
			local_index = intersectionsMatrixAlongY->at(ix_,iy_,iz_);
			intersectionsMatrixAlongY->set(ix_,iy_,iz_, local_index + indexOffset);
		}
		else
		{
			local_index = intersectionsMatrixAlongZ->at(ix_,iy_,iz_);
			intersectionsMatrixAlongZ->set(ix_,iy_,iz_, local_index + indexOffset);
		}

		#else // OPTIMIZE_INTERSECTIONS_MANAGEMENT

		if (dir == X_DIR)
		{
			local_index = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix_,iy_,iz_,NX,NY,NZ);
			writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1, local_index + indexOffset, ix_,iy_,iz_,NX,NY,NZ);
		}
		else if (dir == Y_DIR)
		{
			local_index = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix_,iy_,iz_,NX,NY,NZ);
			writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1, local_index + indexOffset, ix_,iy_,iz_,NX,NY,NZ);
		}
		else
		{
			local_index = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix_,iy_,iz_,NX,NY,NZ);
			writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1, local_index + indexOffset, ix_,iy_,iz_,NX,NY,NZ);
		}

		#endif // OPTIMIZE_INTERSECTIONS_MANAGEMENT
	}
}


#if !defined(OPTIMIZE_VERTICES_ADDITION)

inline void Surface::getVertices(double isolevel, int start_z, int end_z, int jump,
								 vector<coordVec> *localVert,
								 vector<VERTEX_TYPE*> *localNormals)
{
	// Marching cubes vertex indices and edges convention; vertInd is indicated in the cube
	//
	//		   v4_________e4_________v5
	//			/|  				/|
	//		e7 / |                 / |
	//		  /  |             e5 /  |
	//		 /   | e8            /	 | e9
	//	  v7/____|_____e6_______/v6	 |
	//		|	 |              |	 |
	//	 	|  v0|______e0______|____|v1
	//	e11 |	/        		|   /
	//		|  /			e10	|  / 
	//		| /	e3				| / e1
	//		|/					|/
	//	  v3/_________e2________/v2
	//            

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	
	double **positions = allocateMatrix2D<double>(8,3);
	double *gx = delphi->x;
	double *gy = delphi->y;
	double *gz = delphi->z;
	double d_hside = delphi->hside;


	int thread_id = start_z - 1;

	// start_z and end_z have to be inputted >0 and NZ-1 so that the boundary of the grid
	// is skipped: there will never be triangles here if the grid is correctly built
	for (int k = start_z; k < end_z; k += jump)
	{
		// skip the boundary of the grid
		for (int j = 1; j < NY-1; j++)
			for (int i = 1; i < NX-1; i++)
			{
				double votes[8];
				double table[6];
				double xg, yg, zg;

				// if an accurate description is available each inside/outside flag
				// and edge surface intersections were analytically computed by the ray-tracing routine
				// this allows to perform a marching cubes where each vertex of the mesh is granted to
				// to belong to the surface. 
				// This is an analytical marching cubes
				if (accurateTriangulation && !isAvailableScalarField)
				{
					#if !defined(USE_COMPRESSED_GRIDS)
					if (!optimizeGrids)
					{
						if (verticesInsidenessMap[k][j][i])
							votes[0] = +1;
						else
							votes[0] = -1;

						if (verticesInsidenessMap[k][j][i+1])
							votes[1] = +1;
						else
							votes[1] = -1;

						if (verticesInsidenessMap[k][j+1][i+1])
							votes[2] = +1;
						else
							votes[2] = -1;

						if (verticesInsidenessMap[k][j+1][i])
							votes[3] = +1;
						else
							votes[3] = -1;

						if (verticesInsidenessMap[k+1][j][i])
							votes[4] = +1;
						else
							votes[4] = -1;

						if (verticesInsidenessMap[k+1][j][i+1])
							votes[5] = +1;
						else
							votes[5] = -1;

						if (verticesInsidenessMap[k+1][j+1][i+1])
							votes[6] = +1;
						else
							votes[6] = -1;

						if (verticesInsidenessMap[k+1][j+1][i])
							votes[7] = +1;
						else
							votes[7] = -1;
					}
					else
					#endif
					{
						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k,NX,NY,NZ))
							votes[0] = +1;
						else
							votes[0] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k,NX,NY,NZ))
							votes[1] = +1;
						else
							votes[1] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k,NX,NY,NZ))
							votes[2] = +1;
						else
							votes[2] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k,NX,NY,NZ))
							votes[3] = +1;
						else
							votes[3] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k+1,NX,NY,NZ))
							votes[4] = +1;
						else
							votes[4] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k+1,NX,NY,NZ))
							votes[5] = +1;
						else
							votes[5] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k+1,NX,NY,NZ))
							votes[6] = +1;
						else
							votes[6] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k+1,NX,NY,NZ))
							votes[7] = +1;
						else
							votes[7] = -1;
					}
				}
				// Classical marching cube interpolating scalar field values
				else if (accurateTriangulation && isAvailableScalarField)
				{
					votes[0] = scalarField[k][j][i];
					votes[1] = scalarField[k][j][i+1];
					votes[2] = scalarField[k][j+1][i+1];
					votes[3] = scalarField[k][j+1][i];
					votes[4] = scalarField[k+1][j][i];
					votes[5] = scalarField[k+1][j][i+1];
					votes[6] = scalarField[k+1][j+1][i+1];
					votes[7] = scalarField[k+1][j+1][i];

					xg = gx[i];
					yg = gy[j];
					zg = gz[k];
					table[0] = xg - d_hside;
					table[1] = xg + d_hside;
					table[2] = yg - d_hside;
					table[3] = yg + d_hside;
					table[4] = zg - d_hside;
					table[5] = zg + d_hside;

					for (int vertInd=0; vertInd<8; vertInd++)
					{
						positions[vertInd][0] = table[mulTable[vertInd][0]];
						positions[vertInd][1] = table[mulTable[vertInd][1]];
						positions[vertInd][2] = table[mulTable[vertInd][2]];
					}
				}
				else
				{
					for (int vertInd=0; vertInd<8; vertInd++)
						votes[vertInd] = getInsidness(i,j,k,vertInd);
				}
				
				int cubeindex = classifyCube(votes, isolevel);
		
				if (cubeindex != -1)
					// activeCubes[k][j][i] = true;
					#if !defined(USE_COMPRESSED_GRIDS)
					if (!optimizeGrids)
						write3DVector<bool>(activeCubes,true,i,j,k,NX,NY,NZ);
					else
					#endif
					{
						write32xCompressedGrid(compressed_activeCubes,true,i,j,k,NX,NY,NZ);
					}
				else
					continue;

				int vert_offset;
				int ix=i, iy=j, iz=k;

				// 5
				if (edgeTable[cubeindex] & 32)
				{
					vert_offset = -1;
					#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
					vert_offset = intersectionsMatrixAlongY->at(ix+1,iy,iz+1);
					#else
					vert_offset = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix+1,iy,iz+1,NX,NY,NZ);
					#endif

					if (vert_offset == -1)
					{
						verticesBuffers[thread_id].push_back(0.);
						verticesBuffers[thread_id].push_back(0.);
						verticesBuffers[thread_id].push_back(0.);

						VERTEX_TYPE *p = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

						if (accurateTriangulation && isAvailableScalarField)
							vertexInterp(isolevel,positions[5],positions[6],votes[5],votes[6],p);
						else
						{
							p[0] = delphi->x[ix+1]-delphi->hside;
							p[1] = delphi->y[iy];
							p[2] = delphi->z[iz+1]-delphi->hside;
						}
						coordVec v(ix+1,iy,iz+1,p,Y_DIR);
						localVert->push_back(v);
						
						if (computeNormals && providesAnalyticalNormals)
						{
							// the normal will be approximated later with the nil key values
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);

							VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

							localNormals->push_back(normal);
						}
					}
				}

				// 6
				if (edgeTable[cubeindex] & 64)
				{
					vert_offset = -1;
					#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
					vert_offset = intersectionsMatrixAlongX->at(ix,iy+1,iz+1);
					#else
					vert_offset = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy+1,iz+1,NX,NY,NZ);
					#endif

					if (vert_offset == -1)
					{
						verticesBuffers[thread_id].push_back(0.);
						verticesBuffers[thread_id].push_back(0.);
						verticesBuffers[thread_id].push_back(0.);

						VERTEX_TYPE *p = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

						if (accurateTriangulation && isAvailableScalarField)
							vertexInterp(isolevel,positions[6],positions[7],votes[6],votes[7],p);
						else
						{
							p[0] = delphi->x[ix];
							p[1] = delphi->y[iy+1]-delphi->hside;
							p[2] = delphi->z[iz+1]-delphi->hside;
						}
						coordVec v(ix,iy+1,iz+1,p,X_DIR);
						localVert->push_back(v);

						if (computeNormals && providesAnalyticalNormals)
						{
							// the normal will be approximated later with the nil key values
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);

							VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

							localNormals->push_back(normal);
						}
					}
				}
	
				// 10
				if (edgeTable[cubeindex] & 1024)
				{
					vert_offset = -1;
					#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
					vert_offset = intersectionsMatrixAlongZ->at(ix+1,iy+1,iz);
					#else
					vert_offset = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix+1,iy+1,iz,NX,NY,NZ);
					#endif

					if (vert_offset == -1)
					{
						verticesBuffers[thread_id].push_back(0.);
						verticesBuffers[thread_id].push_back(0.);
						verticesBuffers[thread_id].push_back(0.);

						VERTEX_TYPE *p = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

						if (accurateTriangulation && isAvailableScalarField)
							vertexInterp(isolevel,positions[2],positions[6],votes[2],votes[6],p);
						else
						{
							p[0] = delphi->x[ix+1]-delphi->hside;
							p[1] = delphi->y[iy+1]-delphi->hside;
							p[2] = delphi->z[iz];
						}

						coordVec v(ix+1,iy+1,iz,p,Z_DIR);
						localVert->push_back(v);

						if (computeNormals && providesAnalyticalNormals)
						{
							// the normal will be approximated later with the nil key values
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);
							normalsBuffers[thread_id].push_back(0.);

							VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

							localNormals->push_back(normal);
						}
					}
				}
			} // end y
	} // end k
	deleteMatrix2D<double>(8,3,positions);
}

#else // OPTIMIZE_VERTICES_ADDITION

inline void Surface::getVertices(double isolevel, int start_z, int end_z, int jump,
								 vector<coordVec> *localVert,
								 vector<VERTEX_TYPE*> *localNormals)
{
	// Marching cubes vertex indices and edges convention; vertInd is indicated in the cube
	//
	//		   v4_________e4_________v5
	//			/|  				/|
	//		e7 / |                 / |
	//		  /  |             e5 /  |
	//		 /   | e8            /	 | e9
	//	  v7/____|_____e6_______/v6	 |
	//		|	 |              |	 |
	//	 	|  v0|______e0______|____|v1
	//	e11 |	/        		|   /
	//		|  /			e10	|  /
	//		| /	e3				| / e1
	//		|/					|/
	//	  v3/_________e2________/v2
	//

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	double **positions = allocateMatrix2D<double>(8,3);
	double *gx = delphi->x;
	double *gy = delphi->y;
	double *gz = delphi->z;
	double d_hside = delphi->hside;


	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int numVertices = (int)vertList.size();
	#else
	int numVertices = (int)(vertList.size() / 3.);
	#endif

	for (int n = 0; n < numVertices n++)
	{
		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		int i = (int)rintp((vertList[n][0] - delphi->xmin) * delphi->scale);
		int j = (int)rintp((vertList[n][1] - delphi->ymin) * delphi->scale);
		int k = (int)rintp((vertList[n][2] - delphi->zmin) * delphi->scale);
		#else
		int i = (int)rintp((vertList[ n*3+0 ] - delphi->xmin) * delphi->scale);
		int j = (int)rintp((vertList[ n*3+1 ] - delphi->ymin) * delphi->scale);
		int k = (int)rintp((vertList[ n*3+2 ] - delphi->zmin) * delphi->scale);
		#endif

		#if !defined(USE_COMPRESSED_GRIDS)
		if (!optimizeGrids)
		{
			if (read3DVector<bool>(activeCubes,i,j,k,NX,NY,NZ))
				continue;
			write3DVector<bool>(activeCubes,true,i,j,k,NX,NY,NZ);
		}
		else
		#endif
		{
			if (read32xCompressedGrid(compressed_activeCubes,i,j,k,NX,NY,NZ))
				continue;
			write32xCompressedGrid(compressed_activeCubes,true,i,j,k,NX,NY,NZ);
		}

		double votes[8];
		double table[6];
		double xg, yg, zg;

		// if an accurate description is available each inside/outside flag
		// and edge surface intersections were analytically computed by the ray-tracing routine
		// this allows to perform a marching cubes where each vertex of the mesh is granted to
		// to belong to the surface.
		// This is an analytical marching cubes
		if (accurateTriangulation && !isAvailableScalarField)
		{
			#if !defined(USE_COMPRESSED_GRIDS)
			if (!optimizeGrids)
			{
				if (verticesInsidenessMap[k][j][i])
					votes[0] = +1;
				else
					votes[0] = -1;

				if (verticesInsidenessMap[k][j][i+1])
					votes[1] = +1;
				else
					votes[1] = -1;

				if (verticesInsidenessMap[k][j+1][i+1])
					votes[2] = +1;
				else
					votes[2] = -1;

				if (verticesInsidenessMap[k][j+1][i])
					votes[3] = +1;
				else
					votes[3] = -1;

				if (verticesInsidenessMap[k+1][j][i])
					votes[4] = +1;
				else
					votes[4] = -1;

				if (verticesInsidenessMap[k+1][j][i+1])
					votes[5] = +1;
				else
					votes[5] = -1;

				if (verticesInsidenessMap[k+1][j+1][i+1])
					votes[6] = +1;
				else
					votes[6] = -1;

				if (verticesIread3DVector<bool>nsidenessMap[k+1][j+1][i])
					votes[7] = +1;
				else
					votes[7] = -1;
			}
			else
			#endif
			{
				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k,NX,NY,NZ))
					votes[0] = +1;
				else
					votes[0] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k,NX,NY,NZ))
					votes[1] = +1;
				else
					votes[1] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k,NX,NY,NZ))
					votes[2] = +1;
				else
					votes[2] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k,NX,NY,NZ))
					votes[3] = +1;
				else
					votes[3] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k+1,NX,NY,NZ))
					votes[4] = +1;
				else
					votes[4] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k+1,NX,NY,NZ))
					votes[5] = +1;
				else
					votes[5] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k+1,NX,NY,NZ))
					votes[6] = +1;
				else
					votes[6] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k+1,NX,NY,NZ))
					votes[7] = +1;
				else
					votes[7] = -1;
			}
		}
		else if (!accurateTriangulation)
		{
			for (int vertInd=0; vertInd<8; vertInd++)
				votes[vertInd] = getInsidness(i,j,k,vertInd);
		}

		int cubeindex = classifyCube(votes, isolevel);

		if (cubeindex == -1)
			continue;

		int vert_offset;
		int ix=i, iy=j, iz=k;

		// 5
		if (edgeTable[cubeindex] & 32)
		{
			vert_offset = -1;
			#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
			vert_offset = intersectionsMatrixAlongY->at(ix+1,iy,iz+1);
			#else
			vert_offset = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix+1,iy,iz+1,NX,NY,NZ);
			#endif

			if (vert_offset == -1)
			{
				verticesBuffers[thread_id].push_back(0.);
				verticesBuffers[thread_id].push_back(0.);
				verticesBuffers[thread_id].push_back(0.);

				VERTEX_TYPE *p = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

				if (accurateTriangulation && isAvailableScalarField)
					vertexInterp(isolevel,positions[5],positions[6],votes[5],votes[6],p);
				else
				{
					p[0] = delphi->x[ix+1]-delphi->hside;
					p[1] = delphi->y[iy];
					p[2] = delphi->z[iz+1]-delphi->hside;
				}
				coordVec v(ix+1,iy,iz+1,p,Y_DIR);
				localVert->push_back(v);

				if (computeNormals && providesAnalyticalNormals)
				{
					// the normal will be approximated later with the nil key values
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);

					VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

					localNormals->push_back(normal);
				}
			}
		}

		// 6
		if (edgeTable[cubeindex] & 64)
		{
			vert_offset = -1;
			#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
			vert_offset = intersectionsMatrixAlongX->at(ix,iy+1,iz+1);
			#else
			vert_offset = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy+1,iz+1,NX,NY,NZ);
			#endif

			if (vert_offset == -1)
			{
				verticesBuffers[thread_id].push_back(0.);
				verticesBuffers[thread_id].push_back(0.);
				verticesBuffers[thread_id].push_back(0.);

				VERTEX_TYPE *p = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

				if (accurateTriangulation && isAvailableScalarField)
					vertexInterp(isolevel,positions[6],positions[7],votes[6],votes[7],p);
				else
				{
					p[0] = delphi->x[ix];
					p[1] = delphi->y[iy+1]-delphi->hside;
					p[2] = delphi->z[iz+1]-delphi->hside;
				}
				coordVec v(ix,iy+1,iz+1,p,X_DIR);
				localVert->push_back(v);

				if (computeNormals && providesAnalyticalNormals)
				{
					// the normal will be approximated later with the nil key values
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);

					VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

					localNormals->push_back(normal);
				}
			}
		}

		// 10
		if (edgeTable[cubeindex] & 1024)
		{
			vert_offset = -1;
			#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
			vert_offset = intersectionsMatrixAlongZ->at(ix+1,iy+1,iz);
			#else
			vert_offset = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix+1,iy+1,iz,NX,NY,NZ);
			#endif

			if (vert_offset == -1)
			{
				verticesBuffers[thread_id].push_back(0.);
				verticesBuffers[thread_id].push_back(0.);
				verticesBuffers[thread_id].push_back(0.);

				VERTEX_TYPE *p = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

				if (accurateTriangulation && isAvailableScalarField)
					vertexInterp(isolevel,positions[2],positions[6],votes[2],votes[6],p);
				else
				{
					p[0] = delphi->x[ix+1]-delphi->hside;
					p[1] = delphi->y[iy+1]-delphi->hside;
					p[2] = delphi->z[iz];
				}

				coordVec v(ix+1,iy+1,iz,p,Z_DIR);
				localVert->push_back(v);

				if (computeNormals && providesAnalyticalNormals)
				{
					// the normal will be approximated later with the nil key values
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);
					normalsBuffers[thread_id].push_back(0.);

					VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

					localNormals->push_back(normal);
				}
			}
		}
	}
	deleteMatrix2D<double>(8,3,positions);
}

#endif // OPTIMIZE_VERTICES_ADDITION


/*
// Cache-aware but slower version

inline void Surface::getVertices(double isolevel, int start_z, int end_z, int jump,
								 vector<coordVec> *localVert, vector<double*> *localNormals)
{
	// Marching cubes vertex indices and edges convention; vertInd is indicated in the cube
	//
	//		   v4_________e4_________v5
	//			/|  				/|
	//		e7 / |                 / |
	//		  /  |             e5 /  |
	//		 /   | e8            /	 | e9
	//	  v7/____|_____e6_______/v6	 |
	//		|	 |              |	 |
	//	 	|  v0|______e0______|____|v1
	//	e11 |	/        		|   /
	//		|  /			e10	|  /
	//		| /	e3				| / e1
	//		|/					|/
	//	  v3/_________e2________/v2
	//

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	double **positions = allocateMatrix2D<double>(8,3);
	double *gx = delphi->x;
	double *gy = delphi->y;
	double *gz = delphi->z;
	// double d_hside = delphi->hside;

	unsigned int *stack_verticesInsidenessMap = allocate32xCompressedGrid(true,NX,NY,5);

	// each macro-cube is active if all the spanned 4^3 cubes are active;
	// hence, the grid size is reduced by 4x4x4 (">>2" divides by 4)
	int comp_NX = (NX+4)>>2;
	int comp_NY = (NY+4)>>2;
	int comp_NZ = (NZ+4)>>2;
	// the following is allocated in triangulateSurface()
	// unsigned int *compressed_activeMacroCubes = allocate32xCompressedGrid(0,comp_NX,comp_NY,comp_NZ);

	// start_z and end_z have to be inputted >0 and NZ-1 so that the boundary of the grid
	// is skipped: there will never be triangles here if the grid is correctly built
	for (int kk = start_z; kk < end_z; kk += jump)
	{
		if (kk == 0)
			continue;

		int iters_block_z = 4;
		if (kk == 1) iters_block_z = 3;

		// presence of "1": halo layer has unitary thickness
		for (int k=kk; k < 1 + min(end_z, kk+iters_block_z); k++)
		{
			for (int j = 1; j < NY; j++)
				for (int i = 1; i < NX; i++)
				{
					bool val = read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k,NX,NY,NZ);
					write32xCompressedGrid(stack_verticesInsidenessMap,val,i,j, k - kk, NX,NY,4+1);
				}
		}
		int iters_block_y = 3;

		// skip the boundary of the grid
		for (int jj = 1; jj < NY-1; jj += iters_block_y)
		{
			if (jj >= 4) iters_block_y = 4;
			int iters_block_x = 3;

			// skip the boundary of the grid
			for (int ii = 1; ii < NX-1; ii += iters_block_x)
			{
				if (ii >= 4) iters_block_x = 4;

				bool any_active_cube = false;

				for (int k=kk; k < min(end_z, kk+iters_block_z); k++)
				{
					// skip the boundary of the grid
					for (int j = jj; j < min(NY-1, jj+iters_block_y); j++)
					{
						for (int i = ii; i < min(NX-1, ii+iters_block_x); i++)
						{
							double votes[8];
							double table[6];
							double xg, yg, zg;

							// if an accurate description is available each inside/outside flag
							// and edge surface intersections were analytically computed by the ray-tracing routine
							// this allows to perform a marching cubes where each vertex of the mesh is granted to
							// to belong to the surface.
							// This is an analytical marching cubes
							if (accurateTriangulation && !isAvailableScalarField)
							{
								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k][j][i])
										votes[0] = +1;
									else
										votes[0] = -1;
								}
								else
								{
									if (read32xCompressedGrid(stack_verticesInsidenessMap,i,j,k-kk,NX,NY,4+1))
										votes[0] = +1;
									else
										votes[0] = -1;
								}

								// bool val = read3DVector<bool>(verticesInsidenessMap,i,j,k,NX,NY,NZ);
								// votes[0] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k][j][i+1])
										votes[1] = +1;
									else
										votes[1] = -1;
								}
								else
								{
									if (read32xCompressedGrid(stack_verticesInsidenessMap,i+1,j,k-kk,NX,NY,4+1))
										votes[1] = +1;
									else
										votes[1] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i+1,j,k,NX,NY,NZ);
								// votes[1] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k][j+1][i+1])
										votes[2] = +1;
									else
										votes[2] = -1;
								}
								else
								{
									if (read32xCompressedGrid(stack_verticesInsidenessMap,i+1,j+1,k-kk,NX,NY,4+1))
										votes[2] = +1;
									else
										votes[2] = -1;'
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i+1,j+1,k,NX,NY,NZ);
								// votes[2] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k][j+1][i])
										votes[3] = +1;
									else
										votes[3] = -1;
								}
								else
								{
									if (read32xCompressedGrid(stack_verticesInsidenessMap,i,j+1,k-kk,NX,NY,4+1))
										votes[3] = +1;
									else
										votes[3] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i,j+1,k,NX,NY,NZ);
								// votes[3] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k+1][j][i])
										votes[4] = +1;
									else
										votes[4] = -1;
								}
								else
								{
									if (read32xCompressedGrid(stack_verticesInsidenessMap,i,j,k+1-kk,NX,NY,4+1))
										votes[4] = +1;
									else
										votes[4] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i,j,k+1,NX,NY,NZ);
								// votes[4] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k+1][j][i+1])
										votes[5] = +1;
									else
										votes[5] = -1;
								}
								else
								{
									if (read32xCompressedGrid(stack_verticesInsidenessMap,i+1,j,k+1-kk,NX,NY,4+1))
										votes[5] = +1;
									else
										votes[5] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i+1,j,k+1,NX,NY,NZ);
								// votes[5] = ((val) ? (+1) : (-1));


								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k+1][j+1][i+1])
										votes[6] = +1;
									else
										votes[6] = -1;
								}
								else
								{
									if (read32xCompressedGrid(stack_verticesInsidenessMap,i+1,j+1,k+1-kk,NX,NY,4+1))
										votes[6] = +1;
									else
										votes[6] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i+1,j+1,k+1,NX,NY,NZ);
								// votes[6] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k+1][j+1][i])
										votes[7] = +1;
									else
										votes[7] = -1;
								}
								else
								{
									if (read32xCompressedGrid(stack_verticesInsidenessMap,i,j+1,k+1-kk,NX,NY,4+1))
										votes[7] = +1;
									else
										votes[7] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i,j+1,k+1,NX,NY,NZ);
								// votes[7] = ((val) ? (+1) : (-1));
							}
							// Classical marching cube interpolating scalar field values
							else if (accurateTriangulation && isAvailableScalarField)
							{
								votes[0] = scalarField[k][j][i];
								votes[1] = scalarField[k][j][i+1];
								votes[2] = scalarField[k][j+1][i+1];
								votes[3] = scalarField[k][j+1][i];
								votes[4] = scalarField[k+1][j][i];
								votes[5] = scalarField[k+1][j][i+1];
								votes[6] = scalarField[k+1][j+1][i+1];
								votes[7] = scalarField[k+1][j+1][i];

								xg = gx[i];
								yg = gy[j];
								zg = gz[k];
								table[0] = xg - d_hside;
								table[1] = xg + d_hside;
								table[2] = yg - d_hside;
								table[3] = yg + d_hside;
								table[4] = zg - d_hside;
								table[5] = zg + d_hside;
							}

							if (!accurateTriangulation)
								for (int vertInd=0; vertInd<8; vertInd++)
									votes[vertInd] = getInsidness(i,j,k,vertInd);

							if (accurateTriangulation && isAvailableScalarField)
							{
								for (int vertInd=0; vertInd<8; vertInd++)
								{
									positions[vertInd][0] = table[mulTable[vertInd][0]];
									positions[vertInd][1] = table[mulTable[vertInd][1]];
									positions[vertInd][2] = table[mulTable[vertInd][2]];
								}
							}

							int cubeindex = classifyCube(votes, isolevel);

							if (cubeindex != -1)
							{
								// activeCubes[k][j][i] = true;
								if (!optimizeGrids)
									write3DVector<bool>(activeCubes,true,i,j,k,NX,NY,NZ);
								else
									write32xCompressedGrid(compressed_activeCubes,true,i,j,k,NX,NY,NZ);
								any_active_cube = true;
							}
							else
								continue;

							int vert_offset;
							int ix=i, iy=j, iz=k;

							// 5
							if (edgeTable[cubeindex] & 32)
							{
								vert_offset = -1;
								vert_offset = intersectionsMatrixAlongY->at(ix+1,iy,iz+1);

								if (vert_offset == -1)
								{
									verticesBuffers[thread_id].push_back(0.);
									verticesBuffers[thread_id].push_back(0.);
									verticesBuffers[thread_id].push_back(0.);

									VERTEX_TYPE *p = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

									if (accurateTriangulation && isAvailableScalarField)
										vertexInterp(isolevel,positions[5],positions[6],votes[5],votes[6],p);
									else
									{
										p[0] = delphi->x[ix+1]-delphi->hside;
										p[1] = delphi->y[iy];
										p[2] = delphi->z[iz+1]-delphi->hside;
									}
									coordVec v(ix+1,iy,iz+1,p,Y_DIR);
									localVert->push_back(v);

									if (computeNormals && providesAnalyticalNormals)
									{
										// the normal will be approximated later with the nil key values
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);

										VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

										localNormals->push_back(normal);
									}
								}
							}

							// 6
							if (edgeTable[cubeindex] & 64)
							{
								vert_offset = -1;
								vert_offset = intersectionsMatrixAlongX->at(ix,iy+1,iz+1);

								if (vert_offset == -1)
								{
									verticesBuffers[thread_id].push_back(0.);
									verticesBuffers[thread_id].push_back(0.);
									verticesBuffers[thread_id].push_back(0.);

									VERTEX_TYPE *p = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

									if (accurateTriangulation && isAvailableScalarField)
										vertexInterp(isolevel,positions[6],positions[7],votes[6],votes[7],p);
									else
									{
										p[0] = delphi->x[ix];
										p[1] = delphi->y[iy+1]-delphi->hside;
										p[2] = delphi->z[iz+1]-delphi->hside;
									}
									coordVec v(ix,iy+1,iz+1,p,X_DIR);
									localVert->push_back(v);

									if (computeNormals && providesAnalyticalNormals)
									{
										// the normal will be approximated later with the nil key values
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);

										VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

										localNormals->push_back(normal);
									}
								}
							}

							// 10
							if (edgeTable[cubeindex] & 1024)
							{
								vert_offset = -1;
								vert_offset = intersectionsMatrixAlongZ->at(ix+1,iy+1,iz);

								if (vert_offset == -1)
								{
									verticesBuffers[thread_id].push_back(0.);
									verticesBuffers[thread_id].push_back(0.);
									verticesBuffers[thread_id].push_back(0.);

									VERTEX_TYPE *p = &verticesBuffers[thread_id][ verticesBuffers[thread_id].size()-3 ];

									if (accurateTriangulation && isAvailableScalarField)
										vertexInterp(isolevel,positions[2],positions[6],votes[2],votes[6],p);
									else
									{
										p[0] = delphi->x[ix+1]-delphi->hside;
										p[1] = delphi->y[iy+1]-delphi->hside;
										p[2] = delphi->z[iz];
									}

									coordVec v(ix+1,iy+1,iz,p,Z_DIR);
									localVert->push_back(v);

									if (computeNormals && providesAnalyticalNormals)
									{
										// the normal will be approximated later with the nil key values
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);
										normalsBuffers[thread_id].push_back(0.);

										VERTEX_TYPE *normal = &normalsBuffers[thread_id][ normalsBuffers[thread_id].size()-3 ];

										localNormals->push_back(normal);
									}
								}
							}
						}
					}
				}
				if (any_active_cube)
				{
					write32xCompressedGrid(compressed_activeMacroCubes,1,ii>>2,jj>>2,kk>>2, comp_NX,comp_NY,comp_NZ);
				}
			}
		}
	}
	deleteMatrix2D<double>(8,3,positions);

	deleteVector<unsigned int>(stack_verticesInsidenessMap);
}
*/


#if !defined(OPTIMIZE_VERTICES_ADDITION)

inline double Surface::triangulationKernel(double isolevel, bool revert, int start_z, int end_z, int jump,
										   vector<int> *localTriList, VERTEX_TYPE *localArea)
{
	// VERTEX_TYPE **positions = allocateMatrix2D<double>(8,3);
	int **triangles = allocateMatrix2D<int>(5,3);
	VERTEX_TYPE **vertPointers = allocateVector<VERTEX_TYPE*>(3);

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	double *gx = delphi->x;
	double *gy = delphi->y;
	double *gz = delphi->z;
	// double d_hside = delphi->hside;

	double area = 0;


	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int numVertices = (int)vertList.size();
	#else
	int numVertices = (int)(vertList.size() / 3.);
	#endif

	// jump is the number of threads used here
	localTriList->reserve( 3 * max(10, (int)(2*numVertices / (double)jump)) );

	// start_z and end_z have to be inputted >0 and NZ-1 so that the boundary of the grid
	// is skipped: there will never be triangles here if the grid is correctly built
	for (int k = start_z; k < end_z; k += jump)
	{
		// skip the boundary of the grid
		for (int j = 1; j < NY-1; j++)
			for (int i = 1; i < NX-1; i++)
			{
				// early voxel culling
				// if (!activeCubes[k][j][i])
				#if !defined(USE_COMPRESSED_GRIDS)
				if (!optimizeGrids)
				{
					if (!read3DVector<bool>(activeCubes,i,j,k,NX,NY,NZ))
						continue;
				}
				else
				#endif
				{
					if (!read32xCompressedGrid(compressed_activeCubes,i,j,k,NX,NY,NZ))
						continue;
				}

				double votes[8];

				// if an accurate description is available each inside/outside flag
				// and edge surface intersections were analytically computed by the ray-tracing routine
				// this allows to perform a marching cubes where each vertex of the mesh is granted to
				// to belong to the surface. 
				// This is an analytical marching cubes
				if (accurateTriangulation && !isAvailableScalarField)
				{
					#if !defined(USE_COMPRESSED_GRIDS)
					if (!optimizeGrids)
					{
						if (verticesInsidenessMap[k][j][i])
							votes[0] = +1;
						else
							votes[0] = -1;

						if (verticesInsidenessMap[k][j][i+1])
							votes[1] = +1;
						else
							votes[1] = -1;

						if (verticesInsidenessMap[k][j+1][i+1])
							votes[2] = +1;
						else
							votes[2] = -1;

						if (verticesInsidenessMap[k][j+1][i])
							votes[3] = +1;
						else
							votes[3] = -1;

						if (verticesInsidenessMap[k+1][j][i])
							votes[4] = +1;
						else
							votes[4] = -1;

						if (verticesInsidenessMap[k+1][j][i+1])
							votes[5] = +1;
						else
							votes[5] = -1;

						if (verticesInsidenessMap[k+1][j+1][i+1])
							votes[6] = +1;
						else
							votes[6] = -1;

						if (verticesInsidenessMap[k+1][j+1][i])
							votes[7] = +1;
						else
							votes[7] = -1;
					}
					else
					#endif
					{
						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k,NX,NY,NZ))
							votes[0] = +1;
						else
							votes[0] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k,NX,NY,NZ))
							votes[1] = +1;
						else
							votes[1] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k,NX,NY,NZ))
							votes[2] = +1;
						else
							votes[2] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k,NX,NY,NZ))
							votes[3] = +1;
						else
							votes[3] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k+1,NX,NY,NZ))
							votes[4] = +1;
						else
							votes[4] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k+1,NX,NY,NZ))
							votes[5] = +1;
						else
							votes[5] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k+1,NX,NY,NZ))
							votes[6] = +1;
						else
							votes[6] = -1;

						if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k+1,NX,NY,NZ))
							votes[7] = +1;
						else
							votes[7] = -1;
					}
				}
				// Classical marching cube interpolating scalar field values
				else if (accurateTriangulation && isAvailableScalarField)
				{
					votes[0] = scalarField[k][j][i];
					votes[1] = scalarField[k][j][i+1];
					votes[2] = scalarField[k][j+1][i+1];
					votes[3] = scalarField[k][j+1][i];
					votes[4] = scalarField[k+1][j][i];
					votes[5] = scalarField[k+1][j][i+1];
					votes[6] = scalarField[k+1][j+1][i+1];
					votes[7] = scalarField[k+1][j+1][i];
				}
				int numTriangles;

				// double xg = gx[i];
				// double yg = gy[j];
				// double zg = gz[k];

				// double table[6];
				// table[0] = xg - d_hside;
				// table[1] = xg + d_hside;
				// table[2] = yg - d_hside;
				// table[3] = yg + d_hside;
				// table[4] = zg - d_hside;
				// table[5] = zg + d_hside;
				
				if (!accurateTriangulation)
					for (int vertInd=0; vertInd<8; vertInd++)
					{
						votes[vertInd] = getInsidness(i,j,k,vertInd);
						// positions[vertInd][0] = table[mulTable[vertInd][0]];
						// positions[vertInd][1] = table[mulTable[vertInd][1]];
						// positions[vertInd][2] = table[mulTable[vertInd][2]];
					}
				// else
				// 	for (int vertInd=0; vertInd<8; vertInd++)
				// 	{
				// 		positions[vertInd][0] = table[mulTable[vertInd][0]];
				// 		positions[vertInd][1] = table[mulTable[vertInd][1]];
				// 		positions[vertInd][2] = table[mulTable[vertInd][2]];
				// 	}
				
				// now that the votes and positions are obtained 
				// triangulate the cell by marching cube rule.
				numTriangles = getTriangles(votes,isolevel,triangles,i,j,k,NX,NY,NZ);
				
				// printf("\n Num tri %d",numTriangles);
				// save triangles
				for (int kk=0; kk<numTriangles; kk++)
				{
					int tri_temp[3];

					for (int ll=0; ll<3; ll++)
					{
						tri_temp[ll] = triangles[kk][ll];
						#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
						vertPointers[ll] = vertList[tri_temp[ll]];
						#else
						vertPointers[ll] = &vertList[ tri_temp[ll]*3 ];
						#endif
					}
					VERTEX_TYPE a=0, b=0, c=0;
					DIST(a,vertPointers[0],vertPointers[1]);
					DIST(b,vertPointers[0],vertPointers[2]);
					DIST(c,vertPointers[1],vertPointers[2]);
					VERTEX_TYPE ttt = (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c);

					if (ttt > 0)
						area += 0.25*sqrt(ttt);

					if (!isAvailableScalarField)
					{
						localTriList->push_back(tri_temp[2]);
						localTriList->push_back(tri_temp[1]);
						localTriList->push_back(tri_temp[0]);
					}
					else
					{
						localTriList->push_back(tri_temp[0]);
						localTriList->push_back(tri_temp[1]);
						localTriList->push_back(tri_temp[2]);
					}
				}
			}
	}
	// end of MC kernel
	// cout << endl <<"I am finished " << start_z << " " << end_z;
	// cout.flush();

	// deleteMatrix2D<VERTEX_TYPE>(8,3,positions);
	deleteMatrix2D<int>(5,3,triangles);
	deleteVector<VERTEX_TYPE*>(vertPointers);

	*localArea = area;

	return area;
}

#else // OPTIMIZE_VERTICES_ADDITION

inline double Surface::triangulationKernel(double isolevel, bool revert, int start_z, int end_z, int jump,
										   vector<int> *localTriList, VERTEX_TYPE *localArea)
{
	// VERTEX_TYPE **positions = allocateMatrix2D<VERTEX_TYPE>(8,3);
	int **triangles = allocateMatrix2D<int>(5,3);
	VERTEX_TYPE **vertPointers = allocateVector<VERTEX_TYPE*>(3);

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	double *gx = delphi->x;
	double *gy = delphi->y;
	double *gz = delphi->z;
	// double d_hside = delphi->hside;

	double area = 0;


	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int numVertices = (int)vertList.size();
	#else
	int numVertices = (int)(vertList.size() / 3.);
	#endif

	// jump is the number of threads used here
	localTriList->reserve( 3 * max(10, (int)(2*numVertices / (double)jump)) );

	#if !defined(USE_COMPRESSED_GRIDS)
	if (!optimizeGrids)
	{
		memset(activeCubes, 0, NX*NY*NZ*sizeof(bool));
	}
	else
	#endif
	{
		int64_t compressed_nx = NX >> 5L;
		if ((compressed_nx << 5L) < NX) ++compressed_nx;
		int64_t tot = compressed_nx*NY*NZ;

		memset(compressed_activeCubes, 0, tot*sizeof(unsigned int));
	}

	for (int n = 0; n < numVertices; n++)
	{
		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		int i = (int)rintp((vertList[n][0] - delphi->xmin) * delphi->scale);
		int j = (int)rintp((vertList[n][1] - delphi->ymin) * delphi->scale);
		int k = (int)rintp((vertList[n][2] - delphi->zmin) * delphi->scale);
		#else
		int i = (int)rintp((vertList[ n*3+0 ] - delphi->xmin) * delphi->scale);
		int j = (int)rintp((vertList[ n*3+1 ] - delphi->ymin) * delphi->scale);
		int k = (int)rintp((vertList[ n*3+2 ] - delphi->zmin) * delphi->scale);
		#endif

		#if !defined(USE_COMPRESSED_GRIDS)
		if (!optimizeGrids)
		{
			if (read3DVector<bool>(activeCubes,i,j,k,NX,NY,NZ))
				continue;
			write3DVector<bool>(activeCubes,true,i,j,k,NX,NY,NZ);
		}
		else
		#endif
		{
			if (read32xCompressedGrid(compressed_activeCubes,i,j,k,NX,NY,NZ))
				continue;
			write32xCompressedGrid(compressed_activeCubes,true,i,j,k,NX,NY,NZ);
		}

		double votes[8];

		// if an accurate description is available each inside/outside flag
		// and edge surface intersections were analytically computed by the ray-tracing routine
		// this allows to perform a marching cubes where each vertex of the mesh is granted to
		// to belong to the surface.
		// This is an analytical marching cubes
		if (accurateTriangulation && !isAvailableScalarField)
		{
			#if !defined(USE_COMPRESSED_GRIDS)
			if (!optimizeGrids)
			{
				if (verticesInsidenessMap[k][j][i])
					votes[0] = +1;
				else
					votes[0] = -1;

				if (verticesInsidenessMap[k][j][i+1])
					votes[1] = +1;
				else
					votes[1] = -1;

				if (verticesInsidenessMap[k][j+1][i+1])
					votes[2] = +1;
				else
					votes[2] = -1;

				if (verticesInsidenessMap[k][j+1][i])
					votes[3] = +1;
				else
					votes[3] = -1;

				if (verticesInsidenessMap[k+1][j][i])
					votes[4] = +1;
				else
					votes[4] = -1;

				if (verticesInsidenessMap[k+1][j][i+1])
					votes[5] = +1;
				else
					votes[5] = -1;

				if (verticesInsidenessMap[k+1][j+1][i+1])
					votes[6] = +1;
				else
					votes[6] = -1;

				if (verticesInsidenessMap[k+1][j+1][i])
					votes[7] = +1;
				else
					votes[7] = -1;
			}
			else
			#endif
			{
				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k,NX,NY,NZ))
					votes[0] = +1;
				else
					votes[0] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k,NX,NY,NZ))
					votes[1] = +1;
				else
					votes[1] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k,NX,NY,NZ))
					votes[2] = +1;
				else
					votes[2] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k,NX,NY,NZ))
					votes[3] = +1;
				else
					votes[3] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k+1,NX,NY,NZ))
					votes[4] = +1;
				else
					votes[4] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k+1,NX,NY,NZ))
					votes[5] = +1;
				else
					votes[5] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k+1,NX,NY,NZ))
					votes[6] = +1;
				else
					votes[6] = -1;

				if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k+1,NX,NY,NZ))
					votes[7] = +1;
				else
					votes[7] = -1;
			}
		}
		int numTriangles;

		double xg = gx[i];
		double yg = gy[j];
		double zg = gz[k];

		// double table[6];
		// table[0] = xg - d_hside;
		// table[1] = xg + d_hside;
		// table[2] = yg - d_hside;
		// table[3] = yg + d_hside;
		// table[4] = zg - d_hside;
		// table[5] = zg + d_hside;

		if (!accurateTriangulation)
			for (int vertInd=0; vertInd<8; vertInd++)
			{
				votes[vertInd] = getInsidness(i,j,k,vertInd);
				// positions[vertInd][0] = table[mulTable[vertInd][0]];
				// positions[vertInd][1] = table[mulTable[vertInd][1]];
				// positions[vertInd][2] = table[mulTable[vertInd][2]];
			}
		// else
		// 	for (int vertInd=0; vertInd<8; vertInd++)
		// 	{
		// 		positions[vertInd][0] = table[mulTable[vertInd][0]];
		// 		positions[vertInd][1] = table[mulTable[vertInd][1]];
		// 		positions[vertInd][2] = table[mulTable[vertInd][2]];
		// 	}

		// now that the votes and positions are obtained
		// triangulate the cell by marching cube rule.
		numTriangles = getTriangles(votes,isolevel,triangles,i,j,k,NX,NY,NZ);

		// printf("\n Num tri %d",numTriangles);
		// save triangles
		for (int kk=0; kk<numTriangles; kk++)
		{
			int tri_temp[3];

			for (int ll=0; ll<3; ll++)
			{
				tri_temp[ll] = triangles[kk][ll];
				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertPointers[ll] = vertList[tri_temp[ll]];
				#else
				vertPointers[ll] = &vertList[ tri_temp[ll]*3 ];
				#endif
			}
			VERTEX_TYPE a=0, b=0, c=0;
			DIST(a,vertPointers[0],vertPointers[1]);
			DIST(b,vertPointers[0],vertPointers[2]);
			DIST(c,vertPointers[1],vertPointers[2]);
			VERTEX_TYPE ttt = (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c);

			if (ttt > 0)
				area += 0.25*sqrt(ttt);
			if (!isAvailableScalarField)
			{
				localTriList->push_back(tri_temp[2]);
				localTriList->push_back(tri_temp[1]);
				localTriList->push_back(tri_temp[0]);
			}
			else
			{
				localTriList->push_back(tri_temp[0]);
				localTriList->push_back(tri_temp[1]);
				localTriList->push_back(tri_temp[2]);
			}
		}
	}
	// end of MC kernel
	// cout << endl <<"I am finished " << start_z << " " << end_z;
	// cout.flush();

	// deleteMatrix2D<double>(8,3,positions);
	deleteMatrix2D<int>(5,3,triangles);
	deleteVector<VERTEX_TYPE*>(vertPointers);

	*localArea = area;

	return area;
}

#endif // OPTIMIZE_VERTICES_ADDITION


/*
// Cache-aware version but slightly slower than the above version

inline double Surface::triangulationKernel(double isolevel, bool revert, int start_z, int end_z, int jump,
										   vector<int> *localTriList, double *localArea)
{
	// VERTEX_TYPE **positions = allocateMatrix2D<VERTEX_TYPE>(8,3);
	int **triangles = allocateMatrix2D<int>(5,3);
	VERTEX_TYPE **vertPointers = allocateVector<VERTEX_TYPE*>(3);

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int numVertices = (int)vertList.size();
	#else
	int numVertices = (int)(vertList.size() / 3.);
	#endif

	// jump is the number of threads used here
	localTriList->reserve( 3 * max(10, (int)(2*numVertices / (double)jump)) );

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	(*localArea) = 0;

	double *gx = delphi->x;
	double *gy = delphi->y;
	double *gz = delphi->z;
	double d_hside = delphi->hside;

	int comp_NX = (NX+4)>>2;
	int comp_NY = (NY+4)>>2;
	int comp_NZ = (NZ+4)>>2;

	for (int kk = start_z; kk < end_z; kk += jump)
	{
		if (kk == 0)
			continue;

		int iters_block_z = 4;
		// if (kk == 1) iters_block_z = 3;

		int iters_block_y = 3;

		// skip the boundary of the grid
		for (int jj = 1; jj < NY-1; jj += iters_block_y)
		{
			if (jj >= 4) iters_block_y = 4;
			int iters_block_x = 3;

			// skip the boundary of the grid
			for (int ii = 1; ii < NX-1; ii += iters_block_x)
			{
				if (ii >= 4) iters_block_x = 4;

				bool any_active_cube = read32xCompressedGrid(compressed_activeMacroCubes,ii>>2,jj>>2,kk>>2, comp_NX,comp_NY,comp_NZ);

				if (!any_active_cube)
					continue;

				// skip the boundary of the grid
				for (int k=kk; k < min(end_z, kk+iters_block_z); k++)
				{
					for (int j = jj; j < min(NY-1, jj+iters_block_y); j++)
					{
						for (int i = ii; i < min(NX-1, ii+iters_block_x); i++)
						{
							// early voxel culling
							// if (!activeCubes[k][j][i])
							if (!optimizeGrids)
							{
								if (!read3DVector<bool>(activeCubes,i,j,k,NX,NY,NZ))
									continue;
							}
							else
								if (!read32xCompressedGrid(compressed_activeCubes,i,j,k,NX,NY,NZ))
									continue;
							}
							double votes[8];

							// if an accurate description is available each inside/outside flag
							// and edge surface intersections were analytically computed by the ray-tracing routine
							// this allows to perform a marching cubes where each vertex of the mesh is granted to
							// to belong to the surface.
							// This is an analytical marching cubes
							if (accurateTriangulation && !isAvailableScalarField)
							{
								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k][j][i])
										votes[0] = +1;
									else
										votes[0] = -1;
								}
								else
								{
									if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k,NX,NY,NZ))
										votes[0] = +1;
									else
										votes[0] = -1;
								}

								// bool val = read3DVector<bool>(verticesInsidenessMap,i,j,k,NX,NY,NZ);
								// votes[0] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k][j][i+1])
										votes[1] = +1;
									else
										votes[1] = -1;
								}
								else
								{
									if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k,NX,NY,NZ))
										votes[1] = +1;
									else
										votes[1] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i+1,j,k,NX,NY,NZ);
								// votes[1] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k][j+1][i+1])
										votes[2] = +1;
									else
										votes[2] = -1;
								}
								else
								{
									if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k,NX,NY,NZ))
										votes[2] = +1;
									else
										votes[2] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i+1,j+1,k,NX,NY,NZ);
								// votes[2] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k][j+1][i])
										votes[3] = +1;
									else
										votes[3] = -1;
								}
								else
								{
									if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k,NX,NY,NZ))
										votes[3] = +1;
									else
										votes[3] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i,j+1,k,NX,NY,NZ);
								// votes[3] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k+1][j][i])
										votes[4] = +1;
									else
										votes[4] = -1;
								}
								else
								{
									if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j,k+1,NX,NY,NZ))
										votes[4] = +1;
									else
										votes[4] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i,j,k+1,NX,NY,NZ);
								// votes[4] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k+1][j][i+1])
										votes[5] = +1;
									else
										votes[5] = -1;
								}
								else
								{
									if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j,k+1,NX,NY,NZ))
										votes[5] = +1;
									else
										votes[5] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i+1,j,k+1,NX,NY,NZ);
								// votes[5] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k+1][j+1][i+1])
										votes[6] = +1;
									else
										votes[6] = -1;
								}
								else
								{
									if (read32xCompressedGrid(compressed_verticesInsidenessMap,i+1,j+1,k+1,NX,NY,NZ))
										votes[6] = +1;
									else
										votes[6] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i+1,j+1,k+1,NX,NY,NZ);
								// votes[6] = ((val) ? (+1) : (-1));

								if (!optimizeGrids)
								{
									if (verticesInsidenessMap[k+1][j+1][i])
										votes[7] = +1;
									else
										votes[7] = -1;
								}
								else
								{
									if (read32xCompressedGrid(compressed_verticesInsidenessMap,i,j+1,k+1,NX,NY,NZ))
										votes[7] = +1;
									else
										votes[7] = -1;
								}

								// val = read3DVector<bool>(verticesInsidenessMap,i,j+1,k+1,NX,NY,NZ);
								// votes[7] = ((val) ? (+1) : (-1));
							}
							// Classical marching cube interpolating scalar field values
							else if (accurateTriangulation && isAvailableScalarField)
							{
								votes[0] = scalarField[k][j][i];
								votes[1] = scalarField[k][j][i+1];
								votes[2] = scalarField[k][j+1][i+1];
								votes[3] = scalarField[k][j+1][i];
								votes[4] = scalarField[k+1][j][i];
								votes[5] = scalarField[k+1][j][i+1];
								votes[6] = scalarField[k+1][j+1][i+1];
								votes[7] = scalarField[k+1][j+1][i];
							}
							int numTriangles;

							double xg = gx[i];
							double yg = gy[j];
							double zg = gz[k];

							// double table[6];
							// table[0] = xg - d_hside;
							// table[1] = xg + d_hside;
							// table[2] = yg - d_hside;
							// table[3] = yg + d_hside;
							// table[4] = zg - d_hside;
							// table[5] = zg + d_hside;

							if (!accurateTriangulation)
								for (int vertInd=0; vertInd<8; vertInd++)
								{
									votes[vertInd] = getInsidness(i,j,k,vertInd);
									// positions[vertInd][0] = table[mulTable[vertInd][0]];
									// positions[vertInd][1] = table[mulTable[vertInd][1]];
									// positions[vertInd][2] = table[mulTable[vertInd][2]];
								}
							// else
							// 	for (int vertInd=0; vertInd<8; vertInd++)
							// 	{
							// 		positions[vertInd][0] = table[mulTable[vertInd][0]];
							// 		positions[vertInd][1] = table[mulTable[vertInd][1]];
							// 		positions[vertInd][2] = table[mulTable[vertInd][2]];
							// 	}

							// now that the votes and positions are obtained
							// triangulate the cell by marching cube rule.
							numTriangles = getTriangles(votes,isolevel,triangles,i,j,k,NX,NY,NZ);

							// printf("\n Num tri %d",numTriangles);
							// save triangles
							for (int m=0; m<numTriangles; m++)
							{
								int tri_temp[3];

								for (int ll=0; ll<3; ll++)
								{
									tri_temp[ll] = triangles[m][ll];
									#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
									vertPointers[ll] = vertList[tri_temp[ll]];
									#else
									vertPointers[ll] = &vertList[ tri_temp[ll]*3 ];
									#endif
								}

								double a=0, b=0, c=0;
								DIST(a,vertPointers[0],vertPointers[1]);
								DIST(b,vertPointers[0],vertPointers[2]);
								DIST(c,vertPointers[1],vertPointers[2]);
								double ttt = (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c);

								if (ttt > 0)
									(*localArea) += 0.25*sqrt(ttt);

								if (!isAvailableScalarField)
								{
									localTriList->push_back(tri_temp[2]);
									localTriList->push_back(tri_temp[1]);
									localTriList->push_back(tri_temp[0]);
								}
								else
								{
									localTriList->push_back(tri_temp[0]);
									localTriList->push_back(tri_temp[1]);
									localTriList->push_back(tri_temp[2]);
								}
							}
						}
					}
				}
			}
		}
	}
	// end of MC kernel
	// cout << endl <<"I am finished " << start_z << " " << end_z;
	// cout.flush();

	// deleteMatrix2D<VERTEX_TYPE>(8,3,positions);
	deleteMatrix2D<int>(5,3,triangles);
	deleteVector<VERTEX_TYPE*>(vertPointers);

	return (*localArea);
}
*/


double Surface::triangulateSurface(bool outputMesh, bool buildAtomsMapHere,
								   double isolevel, const char *fileName, bool revert)
{
	int num_threads = conf.numThreads;

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	auto chrono_start = chrono::high_resolution_clock::now();


	#if !defined(USE_COMPRESSED_GRIDS)
	if (!optimizeGrids)
	{
		if (verticesInsidenessMap == NULL && accurateTriangulation && !isAvailableScalarField)
		{
			cout << endl << WARN << "Cannot triangulate without inside/out info for grid points";
			return 0;
		}
	}
	else
	#endif
	{
		if (compressed_verticesInsidenessMap == NULL && accurateTriangulation && !isAvailableScalarField)
		{
			cout << endl << WARN << "Cannot triangulate without inside/out info for grid points";
			return 0;
		}
	}

	#ifdef ENABLE_BOOST_THREADS
	boost::thread_group thdGroup;
	#endif
	
	// in parallel complete missing vertices and mark active cubes to make faster the second pass
	// in case of classical MC, all vertices are computed on this pass
	vector<vector<int>*> localTri;
	vector<vector<coordVec>*> localVert;
	vector<vector<VERTEX_TYPE*>*> localNormals;

	VERTEX_TYPE *area = allocateVector<VERTEX_TYPE>(num_threads);

	localTri.reserve(num_threads);
	localVert.reserve(num_threads);
	localNormals.reserve(num_threads);

	for (int i=0; i<num_threads; i++)
	{
		localTri.push_back(new vector<int>());
		localVert.push_back(new vector<coordVec>());
		localNormals.push_back(new vector<VERTEX_TYPE*>());

		localVert[i]->reserve(1000);
		localNormals[i]->reserve(1000);
	}

	// activeCubes = allocateMatrix3D<bool>(NX,NY,NZ);
	#if !defined(USE_COMPRESSED_GRIDS)
	if (!optimizeGrids)
	{
		if (activeCubes != NULL)
			deleteVector<bool>(activeCubes);

		activeCubes = allocateVector<bool>(NX*NY*NZ);

		for (int64_t i=0; i < NX*NY*NZ; i++)
			activeCubes[i] = 0;
	}
	else
	#endif
	{
		if (compressed_activeCubes != NULL)
			deleteVector<unsigned int>(compressed_activeCubes);

		compressed_activeCubes = allocate32xCompressedGrid(0U,NX,NY,NZ);
	}

	///////////// starts parallel vertices retrieval ///////////////////////
	// this is time consuming only if the scalar field is available
	// differently, in case of analytical intersections, this step will be very fast
	// and only used to repair variations done by the cavity filling or by rays that missed the target

	cout << endl << INFO << "Generating MC vertices...";
	cout.flush();


	#if !defined(OPTIMIZE_VERTICES_ADDITION)
	// load balanced and cache friendly thread dispatch
	for (int j=0; j<num_threads; j++)
	{
		// voxels with Z coordinates equal to 0 and NZ-1 are skipped
		#ifdef ENABLE_BOOST_THREADS
		thdGroup.create_thread(boost::bind(&Surface::getVertices,this,isolevel,j+1,NZ-1,num_threads,localVert[j],localNormals[j]));
		#else
		getVertices(isolevel,1,NZ-1,1,localVert[0],localNormals[0]);
		#endif
	}
	#else
	getVertices(isolevel,1,NZ-1,1,localVert[0],localNormals[0]);
	#endif

	/*
	int comp_NX = (NX+4)>>2;
	int comp_NY = (NY+4)>>2;
	int comp_NZ = (NZ+4)>>2;
	compressed_activeMacroCubes = allocate32xCompressedGrid(0,comp_NX,comp_NY,comp_NZ);

	// conflict-free cross-thread writes are possible if each thread proceeds launching XY slabs
	// whose thickness is multiple of 4
	int fine_grid_size = 4;
	int jump = num_threads * fine_grid_size;

	for (int j=0; j<num_threads; j++)
	{
		// voxels with Z coordinates equal to 0 and NZ-1 are skipped

		int start = j*fine_grid_size;

		// voxels with Z coordinates equal to 0 and NZ-1 are skipped
		#ifdef ENABLE_BOOST_THREADS
		thdGroup.create_thread(boost::bind(&Surface::getVertices,this,isolevel,start,NZ-1,jump,localVert[j],localNormals[j]));
		#else
		getVertices(isolevel,start,NZ-1,jump,localVert[0],localNormals[0]);
		#endif
	}
	*/
	///////////////////////////////////////////////////////////////////////////////
		
	#ifdef ENABLE_BOOST_THREADS
	// join; final part of the computation of the volume within the surface
	thdGroup.join_all();
	#endif

	int addedVertices = 0;
	// int dup = 0;

	for (int j=0; j<num_threads; j++)
	{
		vector<coordVec> *lv = localVert[j];
		vector<coordVec>::iterator it2;

		addedVertices += (int)localVert[j]->size();

		vector<VERTEX_TYPE*> *ln;
		vector<VERTEX_TYPE*>::iterator it3;

		if (computeNormals && providesAnalyticalNormals)
		{
			ln = localNormals[j];
			it3 = ln->begin();
		}

		// set<string> brutal;

		for (it2 = lv->begin(); it2 != lv->end(); it2++)
		{
			coordVec &v = *it2;

			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			int ind = (int)vertList.size();
			vertList.push_back(v.vec);
			#else
			int ind = (int)(vertList.size() / 3.);
			vertList.push_back(v.vec[0]);
			vertList.push_back(v.vec[1]);
			vertList.push_back(v.vec[2]);
			#endif

			#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
			if (v.dir == X_DIR)
				intersectionsMatrixAlongX->set(v.ix,v.iy,v.iz,ind);
			else if (v.dir == Y_DIR)
				intersectionsMatrixAlongY->set(v.ix,v.iy,v.iz,ind);
			else
				intersectionsMatrixAlongZ->set(v.ix,v.iy,v.iz,ind);
			#else
			if (v.dir == X_DIR)
				writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ind,v.ix,v.iy,v.iz,NX,NY,NZ);
			else if (v.dir == Y_DIR)
				writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ind,v.ix,v.iy,v.iz,NX,NY,NZ);
			else
				writeBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ind,v.ix,v.iy,v.iz,NX,NY,NZ);
			#endif

			// char buff[100];
			// sprintf(buff, "%d_%d_%d_%d", v.ix,v.iy,v.iz,v.dir);
			// set<string>::iterator it = brutal.find(string(buff));
			// if (it != brutal.end())
			// 	dup++;
			// else
			// 	brutal.insert(string(buff));

			if (computeNormals && providesAnalyticalNormals)
			{
				VERTEX_TYPE *vn = *it3;

				#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT) && !defined(AVOID_NORMALS_MATRICES)
				if (v.dir == X_DIR)
					normalsMatrixAlongX->set(v.ix,v.iy,v.iz,ind);
				else if (v.dir == Y_DIR)
					normalsMatrixAlongY->set(v.ix,v.iy,v.iz,ind);
				else
					normalsMatrixAlongZ->set(v.ix,v.iy,v.iz,ind);
				#endif

				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				normalsList.push_back(vn);
				#else
				normalsList.push_back(vn[0]);
				normalsList.push_back(vn[1]);
				normalsList.push_back(vn[2]);
				#endif

				it3++;
			}
		}

		#if defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		verticesBuffers[j].clear();
		#endif
		localVert[j]->clear();
		delete localVert[j];

		if (computeNormals && providesAnalyticalNormals)
		{
			#if defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			normalsBuffers[j].clear();
			#endif
			localNormals[j]->clear();
			delete localNormals[j];
		}
	}

	// cout << endl << "duplicated "<<dup;
	cout << "ok!";
	cout.flush();

	cout << endl << INFO << "MC added "<< addedVertices << " non analytical vertices";
	/*
	int orphans2 = 0;
	int orphans1 = 0;

	// search for orphan cubes
	for (int iz=0; iz<NZ; iz++)
		for (int iy=0; iy<NY; iy++)
			for (int ix=0; ix<NX; ix++)
			{
				int count = 0;

				#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)
				if (intersectionsMatrixAlongX->at(ix,iy,iz) != -1) count++;
				if (intersectionsMatrixAlongY->at(ix+1,iy,iz) != -1) count++;
				if (intersectionsMatrixAlongX->at(ix,iy+1,iz) != -1) count++;
				if (intersectionsMatrixAlongY->at(ix,iy,iz) != -1) count++;
				if (intersectionsMatrixAlongX->at(ix,iy,iz+1) != -1) count++;
				if (intersectionsMatrixAlongY->at(ix+1,iy,iz+1) != -1) count++;
				if (intersectionsMatrixAlongX->at(ix,iy+1,iz+1) != -1) count++;
				if (intersectionsMatrixAlongY->at(ix,iy,iz+1) != -1) count++;
				if (intersectionsMatrixAlongZ->at(ix,iy,iz) != -1) count++;
				if (intersectionsMatrixAlongZ->at(ix+1,iy,iz) != -1) count++;
				if (intersectionsMatrixAlongZ->at(ix+1,iy+1,iz) != -1) count++;
				if (intersectionsMatrixAlongZ->at(ix,iy+1,iz) != -1) count++;
				#else
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy,iz,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix+1,iy,iz,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy+1,iz,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix,iy,i,NX,NY,NZz) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy,iz+1,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix+1,iy,iz+1,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy+1,iz+1,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix,iy,iz+1,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix,iy,iz,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix+1,iy,iz,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix+1,iy+1,iz,NX,NY,NZ) != -1) count++;
				if (readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix,iy+1,iz,NX,NY,NZ) != -1) count++;
				#endif

				if (count == 2)
					orphans2++;

				if (count == 1)
					orphans1++;
			}

	cout << endl << "orphans1 " << orphans1;
	cout << endl << "orphans2 " << orphans2;
	*/

	////////////////////////////// generate triangles /////////////////////////
	// all vertices are computed, stored and uniquely indexed, now get triangles

	#if !defined(OPTIMIZE_VERTICES_ADDITION)
	for (int j=0; j<num_threads; j++)
	{
		// only octrees version;
		// voxels with Z coordinates equal to 0 and NZ-1 are skipped
		#ifdef ENABLE_BOOST_THREADS
		thdGroup.create_thread(boost::bind(&Surface::triangulationKernel,this,isolevel,revert,j+1,NZ-1,num_threads,localTri[j],&area[j]));
		#else
		triangulationKernel(isolevel,revert,1,NZ-1,1,localTri[0],&area[0]);
		#endif
	}
	#else
	triangulationKernel(isolevel,revert,1,NZ-1,1,localTri[0],&area[0]);
	#endif

	/*
	// The following was designed so as to use the cache-aware version of triangulationKernel(...);
	// see the versions of triangulationKernel(...)
	fine_grid_size = 4;
	jump = num_threads * fine_grid_size;

	for (int j=0; j<num_threads; j++)
	{
		// voxels with Z coordinates equal to 0 and NZ-1 are skipped

		int start = j*fine_grid_size;

		// voxels with Z coordinates equal to 0 and NZ-1 are skipped
		#ifdef ENABLE_BOOST_THREADS
		thdGroup.create_thread(boost::bind(&Surface::triangulationKernel,this,isolevel,revert,start,NZ-1,jump,localTri[j],&area[j]));
		#else
		triangulationKernel(isolevel,revert,start,NZ-1,jump,localTri[0],&area[0]);
		#endif
	}
	// deleteVector<unsigned int>(compressed_activeMacroCubes);
	*/

	#ifdef ENABLE_BOOST_THREADS
	// join; final part of the computation of the volume within the surface
	thdGroup.join_all();
	#endif


	// deleteMatrix3D<bool>(delphi->nx,delphi->ny,delphi->nz,activeCubes);
	#if !defined(USE_COMPRESSED_GRIDS)
	if (!optimizeGrids)
		deleteVector<bool>(activeCubes);
	else
	#endif
	{
		if (compressed_activeCubes != NULL)
			deleteVector<unsigned int>(compressed_activeCubes);
	}

	double surf_area = 0;

	int total_num_triangle_data = 0;

	for (int j=0; j<num_threads; j++)
	{
		total_num_triangle_data += localTri[j]->size();
	}
	triList.reserve(total_num_triangle_data);

	#ifdef ENABLE_BOOST_THREADS

	// Reduction copies
	for (int j=0; j<num_threads; j++)
	{
		surf_area += area[j];

		vector<int> *lt = localTri[j];

		for (vector<int>::iterator it = lt->begin(); it != lt->end(); it++)
		{
			triList.push_back(*it);
		}
		delete localTri[j];
	}
	#else
	surf_area = area[0];

	vector<int> *lt = localTri[0];

	// Copies
	for (vector<int>::iterator it = lt->begin(); it != lt->end(); it++)
	{
		triList.push_back(*it);
	}
	delete localTri[0];
	#endif

	localTri.clear();
	localVert.clear();
	localNormals.clear();

	deleteVector<VERTEX_TYPE>(area);

	auto chrono_end = chrono::high_resolution_clock::now();

	chrono::duration<double> MC_time = chrono_end - chrono_start;
	cout << endl << INFO << "MC time is ";
	printf ("%.4e [s]", MC_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	if (!conf.parallelPocketLoop)
	{
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after MC is " << current_mem_in_MB << " MB";
	}
	#endif
	cout << endl << INFO << "Total, grid conformant, surface area is " << setprecision(10) << surf_area << " [A^2]";


	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int numVertices = (int)vertList.size();
	#else
	int numVertices = (int)(vertList.size() / 3.);
	#endif
	int numTriangles = (int)(triList.size() / 3.);

	// DEBUG
	/*
	int dang = 0;
	for (int jj=0; jj<numVertices; jj++)
	{
		bool found = false;
		for (int tt=0; tt<numTriangles; tt++)
		{
			if (triList[ tt*3+0 ] == jj || triList[ tt*3+1 ] == jj || triList[ tt*3+2 ] == jj)
			{
				found = true;
				break;
			}
		}

		if (!found)
			dang++;
	}

	cout << endl << "dangling " << dang;
	*/

	cout << endl << INFO << "Number of vertices " << numVertices << " number of triangles " << numTriangles;

	// check atoms flag
	if (buildAtomsMapHere && vertexAtomsMapFlag)
	{
		auto chrono_start = chrono::high_resolution_clock::now();

		if (vertexAtomsMap != NULL)
			deleteVector<int>(vertexAtomsMap);

		vertexAtomsMap = allocateVector<int>(numVertices);
		buildAtomsMap();

		cout << endl << INFO << "Connecting vertices to atoms...";
		cout.flush();

		for (int i=0; i<numVertices; i++)
		{
			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			#if !defined(FLOAT_VERTICES)
			vdwAccessible(vertList[i], vertexAtomsMap[i]);
			#else
			double v[3];
			v[0] = vertList[i][0];
			v[1] = vertList[i][1];
			v[2] = vertList[i][2];
			vdwAccessible(v, vertexAtomsMap[i]);
			#endif
			#else
			#if !defined(FLOAT_VERTICES)
			vdwAccessible(&vertList[i*3], vertexAtomsMap[i]);
			#else
			double v[3];
			v[0] = vertList[ i*3+0 ];
			v[1] = vertList[ i*3+1 ];
			v[2] = vertList[ i*3+2 ];
			vdwAccessible(v, vertexAtomsMap[i]);
			#endif
			#endif

			if (vertexAtomsMap[i] == -1)
			{
				cout << endl << WARN << "Cannot detect nearest atom for vertex " << i;
			}
		}
		cout << "ok!";
		cout.flush();

		chrono_end = chrono::high_resolution_clock::now();
		chrono::duration<double> other_time = chrono_end - chrono_start;
		cout << endl << INFO << "Vertex atom mapping time is ";
		printf ("%.4e [s]", other_time.count());

		#if !defined(AVOID_MEM_CHECKS)
		if (!conf.parallelPocketLoop)
		{
			double current_mem_in_MB, peak_mem_in_MB;
			getMemSpace (current_mem_in_MB, peak_mem_in_MB);
			cout << endl << INFO << "Memory required after vertex atom mapping is " << current_mem_in_MB << " MB";
		}
		#endif
	}

	// Some more deallocations before the below triList.reserve()
	deleteIntersectionsMatrices();
	#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT) && !defined(AVOID_NORMALS_MATRICES)
	deleteNormalsMatrices();
	#endif

	chrono_start = chrono::high_resolution_clock::now();

	vector<int> appNormals;

	appNormals.reserve(1000);

	// check if it has to approximate normals or they are already available
	if (computeNormals && providesAnalyticalNormals)
	{
		// search for normals to be approximated
		if (normalsList.size() != 0)
		{
			int ind = 0;

			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			for (vector<VERTEX_TYPE*>::iterator it = normalsList.begin(); it != normalsList.end(); it++)
			{
				VERTEX_TYPE *p = *it;
				if (p[0] == 0 && p[1] == 0 && p[2] == 0)
					appNormals.push_back(ind);
				++ind;
			}
			#else
			for (vector<VERTEX_TYPE>::iterator it = normalsList.begin(); it != normalsList.end(); it += 3)
			{
				VERTEX_TYPE *p = &(*it);
				if (p[0] == 0 && p[1] == 0 && p[2] == 0)
					appNormals.push_back(ind);
				++ind;
			}
			#endif
		}
	}

	if (computeNormals && providesAnalyticalNormals)
	{
		// cout << endl << appNormals.size();
		if (appNormals.size() != 0)
		{
			cout << endl << INFO << "Some analytical normals will be approximated...";
			approximateNormals(appNormals, true);
		}
	}
	else if (computeNormals && !providesAnalyticalNormals)
	{
		cout << endl << INFO << "Approximating vertices normals by triangulation...";
		approximateNormals(appNormals, false);
	}

	chrono_end = chrono::high_resolution_clock::now();
	chrono::duration<double> other_time = chrono_end - chrono_start;
	cout << endl << INFO << "Normals' approximations time is ";
	printf ("%.4e [s]", other_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	if (!conf.parallelPocketLoop)
	{
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after normals' approximations is " << current_mem_in_MB << " MB";
	}
	#endif

	#if !defined(AVOID_SAVING_MESH)
	if (outputMesh)
	{
		chrono_start = chrono::high_resolution_clock::now();

		int format = deduceFormat();
		bool f = saveMesh(format, revert, fileName, vertList, triList, normalsList);

		if (!f)
		{
			cout << endl << ERR << "Errors in saving the mesh!";
		}
		else
		{
			auto chrono_end = chrono::high_resolution_clock::now();
			chrono::duration<double> saveMesh_time = chrono_end - chrono_start;
			cout << endl << INFO << "Outputting mesh time (in triangulateSurface()) ";
			printf ("%.4e [s]", saveMesh_time.count());
		}

		/*
		// Code usable for outputting mesh data in binary format

		chrono_start = chrono::high_resolution_clock::now();

		f = saveMeshBinary(format, revert, fileName, vertList, triList, normalsList);

		if (!f)
		{
			cout << endl << ERR << "Errors in saving the mesh in binary format!";
		}
		else
		{
			auto chrono_end = chrono::high_resolution_clock::now();
			chrono::duration<double> saveMesh_time = chrono_end - chrono_start;
			cout << endl << INFO << "Outputting mesh time in binary format (in triangulateSurface()) ";
			printf ("%.4e [s]", saveMesh_time.count());
		}
		*/
	}
	#endif

	if (buildAtomsMapHere && vertexAtomsMapFlag)
		disposeAtomsMap();

	return surf_area;
}


void Surface::updateVertexTrianglesLists (int **vertexTrianglesList, double **planes, unsigned int vertex_flag[], bool doOnlyList)
{
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int nv = (int)vertList.size();
	#else
	int nv = (int)(vertList.size() / 3.);
	#endif
	int nt = (int)(triList.size() / 3.);

	if (!doOnlyList)
	{
		double w[4], *p1, *p2, *p3;

		for (int i=0; i<nt; i++)
		{
			int v1 = triList[ i*3+0 ];
			int v2 = triList[ i*3+1 ];
			int v3 = triList[ i*3+2 ];

			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			p1 = vertList[v1];
			p2 = vertList[v2];
			p3 = vertList[v3];
			#else
			p1 = &vertList[v1*3];
			p2 = &vertList[v2*3];
			p3 = &vertList[v3*3];
			#endif

			vertexTrianglesList[v1][ ++vertexTrianglesList[v1][0] ] = i;
			vertexTrianglesList[v2][ ++vertexTrianglesList[v2][0] ] = i;
			vertexTrianglesList[v3][ ++vertexTrianglesList[v3][0] ] = i;

			// compute plane normal from triangles
			plane3points(p1,p2,p3,w,false);
			planes[i][0] = w[0];
			planes[i][1] = w[1];
			planes[i][2] = w[2];
		}
	}
	else
	{
		for (int i=0; i<nt; i++)
		{
			int v1 = triList[ i*3+0 ];
			int v2 = triList[ i*3+1 ];
			int v3 = triList[ i*3+2 ];

			// read the proper bit
			if (read32xCompressedVector (vertex_flag, v1, nv))
				vertexTrianglesList[v1][ ++vertexTrianglesList[v1][0] ] = i;

			if (read32xCompressedVector (vertex_flag, v2, nv))
				vertexTrianglesList[v2][ ++vertexTrianglesList[v2][0] ] = i;

			if (read32xCompressedVector (vertex_flag, v3, nv))
				vertexTrianglesList[v3][ ++vertexTrianglesList[v3][0] ] = i;
		}
	}
}


void Surface::approximateNormals(vector<int> &appNormals, bool doOnlyList)
{
	#define MAX_NEIGHBOURS 20

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int nv = (int)vertList.size();
	#else
	int nv = (int)(vertList.size() / 3.);
	#endif
	int nt = (int)(triList.size() / 3.);

	if (!doOnlyList)
	{
		// the most optimised version resides in the other IF body
		int **vertexTrianglesList = allocateVector<int*>(nv);

		for (int i=0; i<nv; i++)
		{
			vertexTrianglesList[i] = allocateVector<int>(MAX_NEIGHBOURS + 1);

			// use the first element to update/read length
			vertexTrianglesList[i][0] = 0;
		}

		double **planes = allocateMatrix2D<double>(nt,3);

		unsigned int *vertex_flag;

		updateVertexTrianglesLists (vertexTrianglesList, planes, vertex_flag, doOnlyList);

		for (int iv=0; iv<nv; iv++)
		{
			VERTEX_TYPE norm;

			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			normalsBuffers[0].push_back(0.);
			normalsBuffers[0].push_back(0.);
			normalsBuffers[0].push_back(0.);

			VERTEX_TYPE *meanNormal = &normalsBuffers[0][ normalsBuffers[0].size()-3 ];

			normalsList.push_back(meanNormal);
			#else
			VERTEX_TYPE meanNormal[3];
			#endif

			for (int j=0; j < vertexTrianglesList[iv][0]; j++)
			{
				int triangleID = vertexTrianglesList[iv][ j+1 ];
				ADD(meanNormal,meanNormal,planes[triangleID]);
			}
			NORMALIZE(meanNormal,norm);

			#if defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			normalsList.push_back(meanNormal[0]);
			normalsList.push_back(meanNormal[1]);
			normalsList.push_back(meanNormal[2]);
			#endif
		}

		for (int i = 0; i < nv; i++)
			deleteVector<int>(vertexTrianglesList[i]);
		deleteVector<int*>(vertexTrianglesList);

		deleteMatrix2D<double>(nt,3,planes);
	}
	else
	{
		if (appNormals.size() == 0)
			return;

		// Most optimised version with vertex flagging and memory-saving techniques:
		// Many vertexTrianglesList[*] are not allocated and vertex_flag[] is compressed
		int **vertexTrianglesList = allocateVector<int*>(nv);

		// compressed mask buffer in which each unsigned int stores 32 flags
		unsigned int *vertex_flag = allocate32xCompressedVector(0, nv);

		for (vector<int>::iterator it = appNormals.begin(); it != appNormals.end(); it++)
		{
			unsigned int iv = (*it);

			// switch on the proper bit
			write32xCompressedVector(vertex_flag, 1, iv, nv);

			vertexTrianglesList[iv] = allocateVector<int>(MAX_NEIGHBOURS + 1);

			// use the first element to update/read length
			vertexTrianglesList[iv][0] = 0;
		}

		double **planes;

		updateVertexTrianglesLists(vertexTrianglesList, planes, vertex_flag, doOnlyList);


		double w[4], *p1, *p2, *p3;

		for (vector<int>::iterator it = appNormals.begin(); it != appNormals.end(); it++)
		{
			unsigned int iv = (*it);

			VERTEX_TYPE norm;

			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			VERTEX_TYPE *meanNormal = normalsList[iv];
			#else
			VERTEX_TYPE meanNormal[3];
			#endif
			meanNormal[0] = 0;
			meanNormal[1] = 0;
			meanNormal[2] = 0;

			for (int j=0; j < vertexTrianglesList[iv][0]; j++)
			{
				int triangleID = vertexTrianglesList[iv][ j+1 ];

				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				p1 = vertList[ triList[triangleID*3+0] ];
				p2 = vertList[ triList[triangleID*3+1] ];
				p3 = vertList[ triList[triangleID*3+2] ];
				#else
				p1 = &vertList[ triList[triangleID*3+0]*3 ];
				p2 = &vertList[ triList[triangleID*3+1]*3 ];
				p3 = &vertList[ triList[triangleID*3+2]*3 ];
				#endif

				// compute plane normal from triangles;
				// some normals are recalculated, but we avoid pre-calculating
				// and storing them in an nt*(3 doubles)-sized buffer
				plane3points(p1,p2,p3,w,false);
				ADD(meanNormal,meanNormal,w);
			}
			NORMALIZE(meanNormal,norm);

			#if defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			normalsList[ iv*3+0 ] = meanNormal[0];
			normalsList[ iv*3+1 ] = meanNormal[1];
			normalsList[ iv*3+2 ] = meanNormal[2];
			#endif

			deleteVector<int>(vertexTrianglesList[iv]);
		}
		deleteVector<int*>(vertexTrianglesList);

		deleteVector<unsigned int>(vertex_flag);
	}
}


#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
bool Surface::savePLYMesh(int format, bool revert, const char *fileName,
						  vector<VERTEX_TYPE*> &vertList, vector<int> &triList, vector<VERTEX_TYPE*> &normalsList)
#else
bool Surface::savePLYMesh(int format, bool revert, const char *fileName,
						  vector<VERTEX_TYPE> &vertList, vector<int> &triList, vector<VERTEX_TYPE> &normalsList)
#endif
{
#ifdef PLY_ENABLED
	(void)format;

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	const int numVertices = (int)vertList.size();
	#else
	const int numVertices = (int)(vertList.size() / 3);
	#endif
	const int numTriangles = (int)(triList.size() / 3);

	char fullName[100];
	sprintf(fullName, "%s%s.ply", rootFile.c_str(), fileName);

	std::filebuf fb;
	fb.open(fullName, std::ios::out | std::ios::binary);
	std::ostream outstream(&fb);
	if (outstream.fail())
	{
		cout << endl << WARN << "Cannot write file " << fileName;
		return false;
	}
	cout << endl << INFO << "Writing triangulated surface in PLY binary file format in " << fileName << "...";

	struct vertex3
	{
		VERTEX_TYPE x, y, z;
	};
	struct int3
	{
		int x, y, z;
	};

	std::vector<vertex3> vertStructVec;
	std::vector<int3> triStructVec;
	std::vector<vertex3> normalsStructVec;

	vertStructVec.reserve(numVertices);
	triStructVec.reserve(numTriangles);

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	for (int i=0; i<numVertices; i++)
		vertStructVec.push_back({vertList[i][0], vertList[i][1], vertList[i][2]});
	#else
	for (int i=0; i<numVertices; i++)
		vertStructVec.push_back({vertList[i*3+0], vertList[i*3+1], vertList[i*3+2]});
	#endif

	for (int i=0; i<numTriangles; i++)
	{
		if (!revert)
			triStructVec.push_back({triList[i*3+0], triList[i*3+1], triList[i*3+2]});
		else
			triStructVec.push_back({triList[i*3+2], triList[i*3+1], triList[i*3+0]});
	}

	tinyply::Type coordType = (sizeof(VERTEX_TYPE) == sizeof(double)) ? tinyply::Type::FLOAT64 : tinyply::Type::FLOAT32;
	tinyply::PlyFile mesh_ply;
	mesh_ply.add_properties_to_element("vertex", {"x", "y", "z"},
									   coordType, numVertices,
									   reinterpret_cast<uint8_t *>(vertStructVec.data()),
									   tinyply::Type::INVALID, 0);
	mesh_ply.add_properties_to_element("face", {"vertex_indices"},
									   tinyply::Type::INT32, numTriangles,
									   reinterpret_cast<uint8_t *>(triStructVec.data()),
									   tinyply::Type::UINT8, 3);

	bool hasNormals = false;
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	hasNormals = ((int)normalsList.size() == numVertices);
	#else
	hasNormals = ((int)normalsList.size() == numVertices*3);
	#endif

	if (!normalsList.empty() && !hasNormals)
	{
		cout << endl << WARN << "Skipping normals in PLY output due to inconsistent normal buffer size";
	}

	if (hasNormals)
	{
		normalsStructVec.reserve(numVertices);
		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		for (int i=0; i<numVertices; i++)
			normalsStructVec.push_back({normalsList[i][0], normalsList[i][1], normalsList[i][2]});
		#else
		for (int i=0; i<numVertices; i++)
			normalsStructVec.push_back({normalsList[i*3+0], normalsList[i*3+1], normalsList[i*3+2]});
		#endif

		mesh_ply.add_properties_to_element("vertex", {"nx", "ny", "nz"},
										   coordType, numVertices,
										   reinterpret_cast<uint8_t *>(normalsStructVec.data()),
										   tinyply::Type::INVALID, 0);
	}

	mesh_ply.write(outstream, true);
	return true;
#else
	(void)format;
	(void)revert;
	(void)fileName;
	(void)vertList;
	(void)triList;
	(void)normalsList;
	cout << endl << ERR << "PLY output requested but NanoShaper was built without tinyply support";
	return false;
#endif
}


#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
bool Surface::saveMesh(int format, bool revert, const char *fileName,
					   vector<VERTEX_TYPE*> &vertList, vector<int> &triList, vector<VERTEX_TYPE*> &normalsList)
#else
bool Surface::saveMesh(int format, bool revert, const char *fileName,
					   vector<VERTEX_TYPE> &vertList, vector<int> &triList, vector<VERTEX_TYPE> &normalsList)
#endif
{
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int numVertices = (int)vertList.size();
	#else
	int numVertices = (int)(vertList.size() / 3.);
	#endif
	int numTriangles = (int)(triList.size() / 3.);

	char fullName[100];

	if (format == PLY)
	{
		if (!savePLYMesh(format, revert, fileName, vertList, triList, normalsList))
			return false;
	}
	else if (format == OFF || format == OFF_A || format == OFF_N || format == OFF_N_A)
	{
		// save all in OFF format
		FILE *fp;
		sprintf(fullName, "%s%s.off", rootFile.c_str(), fileName);
		fp = fopen(fullName, "w");

		if (fp == NULL)
		{
			cout << endl << WARN << "Cannot write file " << fileName;
			return false;
		}
		if (format == OFF_A)
		{
			if (vertexAtomsMap == NULL)
			{
				cout << endl << ERR << "Cannot save in OFF+A format if nearest atom info is not available";
				fclose(fp);
				return false;
			}
			fprintf(fp, "OFF+A\n");
			cout << endl << INFO << "Writing triangulated surface in OFF+A file format in " << fileName << "...";
		}
		else if (format == OFF_N)
		{			
			if (normalsList.size() == 0)
			{
				cout << endl << ERR << "Cannot save in OFF+N format if normals are not available";
				fclose(fp);
				return false;
			}
			fprintf(fp, "OFF+N\n");
			cout << endl << INFO << "Writing triangulated surface in OFF+N file format in " << fileName << "...";
		}
		else if (format == OFF_N_A)
		{
			if (normalsList.size() == 0)
			{
				cout << endl << ERR << "Cannot save in OFF+N+A format if normals are not available";
				fclose(fp);
				return false;
			}
			if (vertexAtomsMap == NULL)
			{
				cout << endl << ERR << "Cannot save in OFF+N+A format if nearest atom info is not available";
				fclose(fp);
				return false;
			}
			fprintf(fp, "OFF+N+A\n");
			cout << endl << INFO << "Writing triangulated surface in OFF+N+A file format in " << fileName << "...";
		}
		else
		{
			fprintf(fp, "OFF\n");
			cout << endl << INFO << "Writing triangulated surface in OFF file format in " << fileName << "...";
		}
		cout.flush();

		time_t pt;
		time(&pt);

		fprintf(fp, "# File created by %s version %s date %s\n", PROGNAME, VERSION, ctime(&pt));
		fprintf(fp, "%d %d 0\n", numVertices, numTriangles);

		for (int i=0; i<numVertices; i++)
		{
			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			if (format == OFF_A)
				fprintf(fp, "%.3f %.3f %.3f %d\n", vertList[i][0],vertList[i][1],vertList[i][2],vertexAtomsMap[i]);
			else if (format == OFF_N)
				fprintf(fp,"%.3f %.3f %.3f %.3f %.3f %.3f\n", vertList[i][0],vertList[i][1],vertList[i][2],normalsList[i][0],normalsList[i][1],normalsList[i][2]);
			else if (format == OFF_N_A)
				fprintf(fp, "%.3f %.3f %.3f %.3f %.3f %.3f %d\n", vertList[i][0],vertList[i][1],vertList[i][2],normalsList[i][0],normalsList[i][1],normalsList[i][2],vertexAtomsMap[i]);
			else
				fprintf(fp, "%.3f %.3f %.3f\n", vertList[i][0],vertList[i][1],vertList[i][2]);
			#else
			if (format == OFF_A)
				fprintf(fp, "%.3f %.3f %.3f %d\n", vertList[ i*3+0 ],vertList[ i*3+1 ],vertList[ i*3+2 ],vertexAtomsMap[i]);
			else if (format == OFF_N)
				fprintf(fp,"%.3f %.3f %.3f %.3f %.3f %.3f\n", vertList[ i*3+0 ],vertList[ i*3+1 ],vertList[ i*3+2 ],normalsList[ i*3+0 ],normalsList[ i*3+1 ],normalsList[ i*3+2 ]);
			else if (format == OFF_N_A)
				fprintf(fp, "%.3f %.3f %.3f %.3f %.3f %.3f %d\n", vertList[ i*3+0 ],vertList[ i*3+1 ],vertList[ i*3+2 ],normalsList[ i*3+0 ],normalsList[ i*3+1 ],normalsList[ i*3+2 ],vertexAtomsMap[i]);
			else
				fprintf(fp, "%.3f %.3f %.3f\n", vertList[ i*3+0 ],vertList[ i*3+1 ],vertList[ i*3+2 ]);
			#endif
		}
		for (int i=0; i<numTriangles; i++)
			if (!revert)
				fprintf(fp, "3 %d %d %d\n", triList[ i*3+0 ],triList[ i*3+1 ],triList[ i*3+2 ]);
			else
				fprintf(fp, "3 %d %d %d\n", triList[ i*3+2 ],triList[ i*3+1 ],triList[ i*3+0 ]);

		fclose(fp);
	}
	else if (format == MSMS || format == MSMS_NO_A)
	{
		if (normalsList.size() == 0)
		{
			cout << endl << ERR << "Cannot save in MSMS format if normals are not available";
			return false;
		}
		if (format == MSMS)
		{
			cout << endl << INFO << "Saving in MSMS format, no patch info...";
			if (vertexAtomsMap == NULL)
			{
				cout << endl << ERR << "Cannot save in MSMS format if nearest atom info is not available";
				return false;
			}
		}
		else if (format == MSMS_NO_A)
			cout << endl << INFO << "Saving in MSMS format, no patch info, no nearest atom...";

		FILE *fp1, *fp2;
		sprintf(fullName, "%s%s.face", rootFile.c_str(),fileName);
		fp1 = fopen(fullName, "w");
		sprintf(fullName, "%s%s.vert", rootFile.c_str(),fileName);
		fp2 = fopen(fullName, "w");

		if (fp1 == NULL || fp2 == NULL)
		{
			cout << endl << ERR << "Error in writing MSMS files";
			return false;
		}
		time_t pt;
		time(&pt);    
		
		fprintf(fp1, "# File created by %s version %s date %s", PROGNAME,VERSION,ctime(&pt));
		fprintf(fp1, "#faces\n");
		fprintf(fp1, "%d\n", numTriangles);

		for (int i=0; i<numTriangles; i++)
			if (!revert)
				fprintf(fp1, "%6d %6d %6d %2d %6d\n", triList[ i*3+0 ]+1,triList[ i*3+1 ]+1,triList[ i*3+2 ]+1,1,1);
			else
				fprintf(fp1, "%6d %6d %6d %2d %6d\n", triList[ i*3+2 ]+1,triList[ i*3+1 ]+1,triList[ i*3+0 ]+1,1,1);

		fclose(fp1);

		fprintf(fp2, "# File created by %s version %s date %s", PROGNAME,VERSION,ctime(&pt));
		fprintf(fp2, "#vertex\n");
		fprintf(fp2, "%d\n", numVertices);
	
		for (int i=0; i<numVertices; i++)
		{
			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			if (format == MSMS)
				fprintf(fp2, "%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d \n",
						vertList[i][0],vertList[i][1],vertList[i][2],normalsList[i][0],normalsList[i][1],normalsList[i][2],0,vertexAtomsMap[i]+1,0);
			else
				fprintf(fp2, "%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d \n",
						vertList[i][0],vertList[i][1],vertList[i][2],normalsList[i][0],normalsList[i][1],normalsList[i][2],0,1,0);
			#else
			if (format == MSMS)
				fprintf(fp2, "%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d \n",
						vertList[ i*3+0 ],vertList[ i*3+1 ],vertList[ i*3+2 ],normalsList[ i*3+0 ],normalsList[ i*3+1 ],normalsList[ i*3+2 ],0,vertexAtomsMap[i]+1,0);
				else
					fprintf(fp2, "%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d \n",
							vertList[ i*3+0 ],vertList[ i*3+1 ],vertList[ i*3+2 ],normalsList[ i*3+0 ],normalsList[ i*3+1 ],normalsList[ i*3+2 ],0,1,0);
			#endif
		}
		fclose(fp2);
	}

	sprintf(fullName, "%s%striangleAreas.txt", rootFile.c_str(), fileName);
	FILE *fp = fopen(fullName, "w");
	// compute per triangle area, skip degenerate ones
	for (int i=0; i<numTriangles; i++)
	{
		VERTEX_TYPE a=0, b=0, c=0;
		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		DIST(a,(vertList[triList[ i*3+0 ]]),(vertList[triList[ i*3+1 ]]))
		DIST(b,(vertList[triList[ i*3+0 ]]),(vertList[triList[ i*3+2 ]]))
		DIST(c,(vertList[triList[ i*3+1 ]]),(vertList[triList[ i*3+2 ]]))
		#else
		DIST(a,(&vertList[triList[ i*3+0 ]*3]),(&vertList[triList[ i*3+1 ]*3]))
		DIST(b,(&vertList[triList[ i*3+0 ]*3]),(&vertList[triList[ i*3+2 ]*3]))
		DIST(c,(&vertList[triList[ i*3+1 ]*3]),(&vertList[triList[ i*3+2 ]*3]))
		#endif
		VERTEX_TYPE ttt = (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c);

		if (ttt > 0)
			fprintf(fp, "%f\n", 0.25*sqrt(ttt));
		else
			fprintf(fp, "0.\n");
	}
	fclose(fp);
	
	return true;
}


#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
bool Surface::saveMeshBinary(int format, bool revert, const char *fileName,
							 vector<VERTEX_TYPE*> &vertList, vector<int> &triList, vector<VERTEX_TYPE*> &normalsList)
#else
bool Surface::saveMeshBinary(int format, bool revert, const char *fileName,
							 vector<VERTEX_TYPE> &vertList, vector<int> &triList, vector<VERTEX_TYPE> &normalsList)
#endif
{
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int numVertices = (int)vertList.size();
	#else
	int numVertices = (int)(vertList.size() / 3.);
	#endif
	int numTriangles = (int)(triList.size() / 3.);

	char fullName[100];

	if (format == PLY)
	{
		cout << endl << ERR << "PLY format is only supported through saveMesh(), not saveMeshBinary()";
		return false;
	}

	if (format == OFF || format == OFF_A || format == OFF_N || format == OFF_N_A)
	{
		// save all in OFF format
		FILE *fp;
		sprintf(fullName, "%s%s.off.bin", rootFile.c_str(), fileName);
		fp = fopen(fullName, "wb");

		if (fp == NULL)
		{
			cout << endl << WARN << "Cannot write file " << fileName;
			return false;
		}
		if (format == OFF_A)
		{
			if (vertexAtomsMap == NULL)
			{
				cout << endl << ERR << "Cannot save in OFF+A format if nearest atom info is not available";
				fclose(fp);
				return false;
			}
			fprintf(fp, "OFF+A\n");
			cout << endl << INFO << "Writing triangulated surface in OFF+A binary file format in " << fileName << "...";
		}
		else if (format == OFF_N)
		{
			if (normalsList.size() == 0)
			{
				cout << endl << ERR << "Cannot save in OFF+N binary format if normals are not available";
				fclose(fp);
				return false;
			}
			fprintf(fp, "OFF+N\n");
			cout << endl << INFO << "Writing triangulated surface in OFF+N binary file format in " << fileName << "...";
		}
		else if (format == OFF_N_A)
		{
			if (normalsList.size() == 0)
			{
				cout << endl << ERR << "Cannot save in OFF+N+A binary format if normals are not available";
				fclose(fp);
				return false;
			}

			if (vertexAtomsMap == NULL)
			{
				cout << endl << ERR << "Cannot save in OFF+N+A binary format if nearest atom info is not available";
				fclose(fp);
				return false;
			}
			fprintf(fp, "OFF+N+A\n");
			cout << endl << INFO << "Writing triangulated surface in OFF+N+A binary file format in " << fileName << "...";
		}
		else
		{
			fprintf(fp, "OFF\n");
			cout << endl << INFO << "Writing triangulated surface in OFF binary file format in " << fileName << "...";
		}
		cout.flush();

		time_t pt;
		time(&pt);

		fprintf(fp, "# File created by %s version %s date %s\n", PROGNAME, VERSION, ctime(&pt));
		fprintf(fp, "%d %d 0\n", numVertices, numTriangles);

		/*
		// Fastest version, but more memory costly due to the copy of data into arrays with other formats
		struct vertexData1
		{
			float vtx1, vtx2, vtx3;

			int atomID;
		};
		struct vertexData2
		{
			float vtx1, vtx2, vtx3;
			float nor1, nor2, nor3;
		};
		struct vertexData3
		{
			float vtx1, vtx2, vtx3;
			float nor1, nor2, nor3;

			int atomID;
		};
		struct vertexData4
		{
			float vtx1, vtx2, vtx3;
		};

		if (format == OFF_A)
		{
			vertexData1 *vertex_data = allocateVector<vertexData1>(numVertices);

			for (int i=0; i<numVertices; i++)
			{
				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertex_data[i].vtx1 = vertList[i][0];
				vertex_data[i].vtx2 = vertList[i][1];
				vertex_data[i].vtx3 = vertList[i][2];
				#else
				vertex_data[i].vtx1 = vertList[ i*3+0 ];
				vertex_data[i].vtx2 = vertList[ i*3+1 ];
				vertex_data[i].vtx3 = vertList[ i*3+2 ];
				#endif

				vertex_data[i].atomID = vertexAtomsMap[i];
			}
			fwrite(vertex_data, sizeof(struct vertexData1), numVertices, fp);
			deleteVector<vertexData1>(vertex_data);
		}
		else if (format == OFF_N)
		{
			vertexData2 *vertex_data = allocateVector<vertexData2>(numVertices);

			for (int i=0; i<numVertices; i++)
			{
				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertex_data[i].vtx1 = vertList[i][0];
				vertex_data[i].vtx2 = vertList[i][1];
				vertex_data[i].vtx3 = vertList[i][2];
				vertex_data[i].nor1 = normalsList[i][0];
				vertex_data[i].nor2 = normalsList[i][1];
				vertex_data[i].nor3 = normalsList[i][2];
				#else
				vertex_data[i].vtx1 = vertList[ i*3+0 ];
				vertex_data[i].vtx2 = vertList[ i*3+1 ];
				vertex_data[i].vtx3 = vertList[ i*3+2 ];
				vertex_data[i].nor1 = normalsList[ i*3+0 ];
				vertex_data[i].nor2 = normalsList[ i*3+1 ];
				vertex_data[i].nor3 = normalsList[ i*3+2 ];
				#endif
			}
			fwrite(vertex_data, sizeof(struct vertexData2), numVertices, fp);
			deleteVector<vertexData2>(vertex_data);
		}
		else if (format == OFF_N_A)
		{
			vertexData3 *vertex_data = allocateVector<vertexData3>(numVertices);

			for (int i=0; i<numVertices; i++)
			{
				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertex_data[i].vtx1 = vertList[i][0];
				vertex_data[i].vtx2 = vertList[i][1];
				vertex_data[i].vtx3 = vertList[i][2];
				vertex_data[i].nor1 = normalsList[i][0];
				vertex_data[i].nor2 = normalsList[i][1];
				vertex_data[i].nor3 = normalsList[i][2];
				#else
				vertex_data[i].vtx1 = vertList[ i*3+0 ];
				vertex_data[i].vtx2 = vertList[ i*3+1 ];
				vertex_data[i].vtx3 = vertList[ i*3+2 ];
				vertex_data[i].nor1 = normalsList[ i*3+0 ];
				vertex_data[i].nor2 = normalsList[ i*3+1 ];
				vertex_data[i].nor3 = normalsList[ i*3+2 ];
				#endif

				vertex_data[i].atomID = vertexAtomsMap[i];
			}
			fwrite(vertex_data, sizeof(struct vertexData3), numVertices, fp);
			deleteVector<vertexData3>(vertex_data);
		}
		else
		{
			vertexData4 *vertex_data = allocateVector<vertexData4>(numVertices);

			for (int i=0; i<numVertices; i++)
			{
				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertex_data[i].vtx1 = vertList[i][0];
				vertex_data[i].vtx2 = vertList[i][1];
				vertex_data[i].vtx3 = vertList[i][2];
				#else
				vertex_data[i].vtx1 = vertList[ i*3+0 ];
				vertex_data[i].vtx2 = vertList[ i*3+1 ];
				vertex_data[i].vtx3 = vertList[ i*3+2 ];
				#endif
			}
			fwrite(vertex_data, sizeof(struct vertexData4), numVertices, fp);
			deleteVector<vertexData4>(vertex_data);
		}

		struct triangleData {
			int id1, id2, id3;
		};
		triangleData *triangle_data = allocateVector<triangleData>(numTriangles);

		if (!revert)
		{
			for (int i=0; i<numTriangles; i++)
			{
				triangle_data[i].id1 = triList[ i*3+0 ];
				triangle_data[i].id2 = triList[ i*3+1 ];
				triangle_data[i].id3 = triList[ i*3+2 ];
			}
		}
		else
		{
			for (int i=0; i<numTriangles; i++)
			{
				triangle_data[i].id1 = triList[ i*3+2 ];
				triangle_data[i].id2 = triList[ i*3+1 ];
				triangle_data[i].id3 = triList[ i*3+0 ];
			}
		}
		fwrite(triangle_data, sizeof(struct triangleData), numTriangles, fp);
		deleteVector<triangleData>(triangle_data);
		*/

		if (format == OFF_A)
		{
			for (int i=0; i<numVertices; i++)
			{
				float vertex_data[3];

				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertex_data[0] = vertList[i][0];
				vertex_data[1] = vertList[i][1];
				vertex_data[2] = vertList[i][2];
				#else
				vertex_data[0] = vertList[ i*3+0 ];
				vertex_data[1] = vertList[ i*3+1 ];
				vertex_data[2] = vertList[ i*3+2 ];
				#endif
				fwrite(vertex_data, sizeof(float), 3, fp);

				fwrite(&vertexAtomsMap[i], sizeof(int), 1, fp);
			}
		}
		else if (format == OFF_N)
		{
			for (int i=0; i<numVertices; i++)
			{
				float vertex_data[6];

				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertex_data[0] = vertList[i][0];
				vertex_data[1] = vertList[i][1];
				vertex_data[2] = vertList[i][2];
				vertex_data[3] = normalsList[i][0];
				vertex_data[4] = normalsList[i][1];
				vertex_data[5] = normalsList[i][2];
				#else
				vertex_data[0] = vertList[ i*3+0 ];
				vertex_data[1] = vertList[ i*3+1 ];
				vertex_data[2] = vertList[ i*3+2 ];
				vertex_data[3] = normalsList[ i*3+0 ];
				vertex_data[4] = normalsList[ i*3+1 ];
				vertex_data[5] = normalsList[ i*3+2 ];
				#endif
				fwrite(vertex_data, sizeof(float), 6, fp);
			}
		}
		else if (format == OFF_N_A)
		{
			for (int i=0; i<numVertices; i++)
			{
				float vertex_data[6];

				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertex_data[0] = vertList[i][0];
				vertex_data[1] = vertList[i][1];
				vertex_data[2] = vertList[i][2];
				vertex_data[3] = normalsList[i][0];
				vertex_data[4] = normalsList[i][1];
				vertex_data[5] = normalsList[i][2];
				#else
				vertex_data[0] = vertList[ i*3+0 ];
				vertex_data[1] = vertList[ i*3+1 ];
				vertex_data[2] = vertList[ i*3+2 ];
				vertex_data[3] = normalsList[ i*3+0 ];
				vertex_data[4] = normalsList[ i*3+1 ];
				vertex_data[5] = normalsList[ i*3+2 ];
				#endif
				fwrite(vertex_data, sizeof(float), 6, fp);

				fwrite(&vertexAtomsMap[i], sizeof(int), 1, fp);
			}
		}
		else
		{
			for (int i=0; i<numVertices; i++)
			{
				float vertex_data[3];

				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				vertex_data[0] = vertList[i][0];
				vertex_data[1] = vertList[i][1];
				vertex_data[2] = vertList[i][2];
				#else
				vertex_data[0] = vertList[ i*3+0 ];
				vertex_data[1] = vertList[ i*3+1 ];
				vertex_data[2] = vertList[ i*3+2 ];
				#endif
				fwrite(vertex_data, sizeof(float), 3, fp);
			}
		}
		if (!revert)
		{
			for (int i=0; i<numTriangles; i++)
			{
				fwrite(&triList[i*3], sizeof(int), 3, fp);
			}
		}
		else
		{
			for (int i=0; i<numTriangles; i++)
			{
				fwrite(&triList[ i*3+2 ], sizeof(int), 1, fp);
				fwrite(&triList[ i*3+1 ], sizeof(int), 1, fp);
				fwrite(&triList[ i*3+0 ], sizeof(int), 1, fp);
			}
		}

		fclose(fp);
	}

	sprintf(fullName, "%s%striangleAreas.bin", rootFile.c_str(), fileName);
	FILE *fp = fopen(fullName, "wb");

	/*
	// Fastest version, but more memory costly due to the copy of data into arrays with other formats

	float *triangle_area = allocateVector<float>(numTriangles);
	// compute per triangle area, skip degenerate ones

	for (int i=0; i<numTriangles; i++)
	{
		VERTEX_TYPE a=0, b=0, c=0;
		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		DIST(a,(vertList[triList[ i*3+0 ]]),(vertList[triList[ i*3+1 ]]))
		DIST(b,(vertList[triList[ i*3+0 ]]),(vertList[triList[ i*3+2 ]]))
		DIST(c,(vertList[triList[ i*3+1 ]]),(vertList[triList[ i*3+2 ]]))
		#else
		DIST(a,(&vertList[triList[ i*3+0 ]*3]),(&vertList[triList[ i*3+1 ]*3]))
		DIST(b,(&vertList[triList[ i*3+0 ]*3]),(&vertList[triList[ i*3+2 ]*3]))
		DIST(c,(&vertList[triList[ i*3+1 ]*3]),(&vertList[triList[ i*3+2 ]*3]))
		#endif
		VERTEX_TYPE ttt = (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c);

		if (ttt > 0)
			triangle_area[i] = 0.25*sqrt(ttt);
		else
			triangle_area[i] = 0.0;
	}
	fwrite(triangle_area, sizeof(float), numTriangles, fp);
	deleteVector<float>(triangle_area);
	*/
	for (int i=0; i<numTriangles; i++)
	{
		VERTEX_TYPE a=0, b=0, c=0;
		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		DIST(a,(vertList[triList[ i*3+0 ]]),(vertList[triList[ i*3+1 ]]))
		DIST(b,(vertList[triList[ i*3+0 ]]),(vertList[triList[ i*3+2 ]]))
		DIST(c,(vertList[triList[ i*3+1 ]]),(vertList[triList[ i*3+2 ]]))
		#else
		DIST(a,(&vertList[triList[ i*3+0 ]*3]),(&vertList[triList[ i*3+1 ]*3]))
		DIST(b,(&vertList[triList[ i*3+0 ]*3]),(&vertList[triList[ i*3+2 ]*3]))
		DIST(c,(&vertList[triList[ i*3+1 ]*3]),(&vertList[triList[ i*3+2 ]*3]))
		#endif
		VERTEX_TYPE ttt = (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c);

		float area;

		if (ttt > 0)
			area = 0.25*sqrt(ttt);
		else
			area = 0.0;

		fwrite(&area, sizeof(float), 1, fp);
	}
	fclose(fp);

	return true;
}


bool Surface::difference(Surface *surf)
{
	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	// S1 is the fat probe
	// S2 is the regular probe

	// int inside1 = inside;

	// grid consistency checks
	if (NX != surf->delphi->nx ||
		NY != surf->delphi->ny ||
		NZ != surf->delphi->nz )
	{
		return false;
	}

	if (delphi->buildEpsMap != surf->delphi->buildEpsMap)
	{
		return false;
	}

	if (delphi->buildStatus != surf->delphi->buildStatus)
	{
		return false;
	}

	// apply difference rule to epsmap, status map, idebmap, insideness map

	// skip epsmap usually not used
	/*
	// epsmap
	if (delphi->buildEpsMap)
	{
		for (int k=0; k<NZ; k++)
			for (int j=0; j<NY; j++)
				for (int i=0; i<NX; i++)
				{
					for (int l=0; l<3; l++)
					{
						int stat1 = delphi->EPSMAP(i,j,k,l,NX,NY,NZ);
						int stat2 = surf->delphi->EPSMAP(i,j,k,l,NX,NY,NZ);
						// if S1 is not out and and S2 is out then that's the pocket
						// the new pocket is marked as outside
						// This is done such that the new objects can be detected as cavities
						// by the cavity detector in a further step such that they can be
						// volume filtered
						if ((stat1 != 0 && stat2 == 0))
							delphi->EPSMAP(i,j,k,l,NX,NY,NZ) = 0;
						else
							delphi->EPSMAP(i,j,k,l,NX,NY,NZ) = inside1;

					}
				}
	}
	*/
	// status, idebmap and verticesInsidenessMap
	if (delphi->buildStatus)
	{
		// free memory
		if (delphi->cavitiesVec != NULL)
		{
			vector<vector<int*>*>::iterator it;
			for (it=delphi->cavitiesVec->begin(); it!=delphi->cavitiesVec->end(); it++)
			{
				vector<int*> *inner = (*it);
				vector<int*>::iterator it2;
				for (it2 = inner->begin(); it2 != inner->end(); it2++)
					deleteVector<int>((*it2));
				delete inner;
			}
			delete delphi->cavitiesVec;
		}

		delphi->cavitiesVec = new vector <vector<int*>*>();

		// add the unique 'cavity'.
		// At this stage the 'cavity' is simply the result of the difference map
		// such that it will be filtered by the digital Connolly filter
		vector<int*> *vv = new vector<int*>();

		delphi->cavitiesVec->push_back(vv);

		// clean memory
		delphi->cavitiesSize.clear();
		delphi->cavitiesFlag.clear();

		// keep it
		delphi->cavitiesFlag.push_back(false);

		int countCubes = 0;

		for (int k=0; k<NZ; k++)
			for (int j=0; j<NY; j++)
				for (int i=0; i<NX; i++)
				{
					// TODO idebmap

					// int stat1 = delphi->status[k][j][i];
					// int stat1 = delphi->STATUSMAP(i,j,k,NX,NY);
					// int stat2 = surf->delphi->status[k][j][i];
					// int stat2 = surf->delphi->STATUSMAP(i,j,k,NX,NY);

					int stat1, stat2;

					stat1 = read3DVector<int>(delphi->status,i,j,k,NX,NY,NZ);
					stat2 = read3DVector<int>(surf->delphi->status,i,j,k,NX,NY,NZ);

					// if S1 is not out and and S2 is out then that's the pocket
					if ((stat1 != STATUS_POINT_TEMPORARY_OUT && stat1 != STATUS_POINT_OUT) &&
						(stat2 == STATUS_POINT_TEMPORARY_OUT || stat2 == STATUS_POINT_OUT))
					{
						// mark as outside temporary, such that cavity detection can work on it
						// delphi->STATUSMAP(i,j,k,NX,NY) = STATUS_POINT_TEMPORARY_OUT;
						// mark as first cavity such that it can be directly filtered
						// delphi->STATUSMAP(i,j,k,NX,NY) = STATUS_FIRST_CAV;
						write3DVector<int>(delphi->status,STATUS_FIRST_CAV,i,j,k,NX,NY,NZ);

						countCubes++;

						int *v = allocateVector<int>(3);
						v[0] = i;
						v[1] = j;
						v[2] = k;
						// store the point as a cavity
						vv->push_back(v);

						// switch surrounding points

						if (accurateTriangulation)
						{
							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								verticesInsidenessMap[k][j][i] = true;
								verticesInsidenessMap[k][j][i+1] = true;
								verticesInsidenessMap[k][j+1][i] = true;
								verticesInsidenessMap[k][j+1][i+1] = true;
								verticesInsidenessMap[k+1][j][i] = true;
								verticesInsidenessMap[k+1][j][i+1] = true;
								verticesInsidenessMap[k+1][j+1][i] = true;
								verticesInsidenessMap[k+1][j+1][i+1] = true;
								/*
								write3DVector<bool>(verticesInsidenessMap,true,i,j,k,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i+1,j,k,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i,j+1,k,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i+1,j+1,k,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i,j,k+1,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i+1,j,k+1,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i,j+1,k+1,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i+1,j+1,k+1,NX,NY,NZ);
								*/
							}
							else
							#endif
							{
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i,j,k,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i+1,j,k,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i,j+1,k,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i+1,j+1,k,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i,j,k+1,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i+1,j,k+1,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i,j+1,k+1,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i+1,j+1,k+1,NX,NY,NZ);
							}
						}
					}
					else
					{
						// delphi->status[k][j][i] = STATUS_POINT_INSIDE;
						// delphi->STATUSMAP(i,j,k,NX,NY) = STATUS_POINT_INSIDE;
						write3DVector<int>(delphi->status,STATUS_POINT_INSIDE,i,j,k,NX,NY,NZ);

						if (accurateTriangulation)
						{
							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								verticesInsidenessMap[k][j][i] = false;
								// write3DVector<bool>(verticesInsidenessMap,false,i,j,k,NX,NY,NZ);

								if ((i+1)<NX)
									verticesInsidenessMap[k][j][i+1] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i+1,j,k,NX,NY,NZ);

								if ((j+1)<NY)
									verticesInsidenessMap[k][j+1][i] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i,j+1,k,NX,NY,NZ);

								if ((i+1)<NX && (j+1)<NY)
									verticesInsidenessMap[k][j+1][i+1] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i+1,j+1,k,NX,NY,NZ);

								if ((k+1)<NZ)
									verticesInsidenessMap[k+1][j][i] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i,j,k+1,NX,NY,NZ);

								if ((i+1)<NX && (k+1)<NZ)
									verticesInsidenessMap[k+1][j][i+1] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i+1,j,k+1,NX,NY,NZ);

								if ((j+1)<NY && (k+1)<NZ)
									verticesInsidenessMap[k+1][j+1][i] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i,j+1,k+1,NX,NY,NZ);

								if ((i+1)<NX && (j+1)<NY && (k+1)<NZ)
									verticesInsidenessMap[k+1][j+1][i+1] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i+1,j+1,k+1,NX,NY,NZ);
							}
							else
							#endif
							{
								write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,j,k,NX,NY,NZ);

								if ((i+1)<NX)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i+1,j,k,NX,NY,NZ);

								if ((j+1)<NY)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,j+1,k,NX,NY,NZ);

								if ((i+1)<NX && (j+1)<NY)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i+1,j+1,k,NX,NY,NZ);

								if ((k+1)<NZ)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,j,k+1,NX,NY,NZ);

								if ((i+1)<NX && (k+1)<NZ)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i+1,j,k+1,NX,NY,NZ);

								if ((j+1)<NY && (k+1)<NZ)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,j+1,k+1,NX,NY,NZ);

								if ((i+1)<NX && (j+1)<NY && (k+1)<NZ)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i+1,j+1,k+1,NX,NY,NZ);
							}
						}
					}
				}
		delphi->cavitiesSize.push_back((delphi->side)*(delphi->side)*(delphi->side)*countCubes);
	}

	// TODO aggregate analytical triangulation points
	// TODO aggregate boundary grid points
	return true;
}


bool Surface::differenceWithBilevelStatusMap(Surface *surf)
{
	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	int64_t coarse_nx, coarse_ny, coarse_nz;
	int64_t coarse_i, coarse_j, coarse_k;
	int64_t fine_i, fine_j, fine_k;

	coarse_nx = getCoarseN(NX);
	coarse_ny = getCoarseN(NY);
	coarse_nz = getCoarseN(NZ);

	// S1 is the fat probe
	// S2 is the regular probe

	// int inside1 = inside;

	// grid consistency checks
	if (NX != surf->delphi->nx ||
		NY != surf->delphi->ny ||
		NZ != surf->delphi->nz )
	{
		return false;
	}

	if (delphi->buildEpsMap != surf->delphi->buildEpsMap)
	{
		return false;
	}

	if (delphi->buildStatus != surf->delphi->buildStatus)
	{
		return false;
	}

	// apply difference rule to epsmap, status map, idebmap, insideness map

	// skip epsmap usually not used

	// epsmap
	//if (delphi->buildEpsMap)
	//{
	//	for (int k=0; k<NZ; k++)
	//		for (int j=0; j<NY; j++)
	//			for (int i=0; i<NX; i++)
	//			{
	//				for (int l=0; l<3; l++)
	//				{
	//					int stat1 = delphi->EPSMAP(i,j,k,l,NX,NY,NZ);
	//					int stat2 = surf->delphi->EPSMAP(i,j,k,l,NX,NY,NZ);
	//					// if S1 is not out and and S2 is out then that's the pocket
	//					// the new pocket is marked as outside
	//					// This is done such that the new objects can be detected as cavities
	//					// by the cavity detector in a further step such that they can be
	//					// volume filtered
	//					if ((stat1 != 0 && stat2 == 0))
	//						delphi->EPSMAP(i,j,k,l,NX,NY,NZ) = 0;
	//					else
	//						delphi->EPSMAP(i,j,k,l,NX,NY,NZ) = inside1;
	//				}
	//			}
	//}

	// status, idebmap and verticesInsidenessMap
	if (delphi->buildStatus)
	{
		// free memory
		if (delphi->cavitiesVec != NULL)
		{
			vector<vector<int*>*>::iterator it;
			for (it=delphi->cavitiesVec->begin(); it!=delphi->cavitiesVec->end(); it++)
			{
				vector<int*> *inner = (*it);
				vector<int*>::iterator it2;
				for (it2 = inner->begin(); it2 != inner->end(); it2++)
					deleteVector<int>((*it2));
				delete inner;
			}
			delete delphi->cavitiesVec;
		}

		delphi->cavitiesVec = new vector <vector<int*>*>();

		// add the unique 'cavity'.
		// At this stage the 'cavity' is simply the result of the difference map
		// such that it will be filtered by the digital Connolly filter
		vector<int*> *vv = new vector<int*>();

		delphi->cavitiesVec->push_back(vv);

		// clean memory
		delphi->cavitiesSize.clear();
		delphi->cavitiesFlag.clear();

		// keep it
		delphi->cavitiesFlag.push_back(false);

		int countCubes = 0;

		for (int64_t k=0; k<NZ; k++)
		{
			coarse_k = getCoarseID(NZ, k);
			fine_k = getFineID(coarse_k, k);

			for (int64_t j=0; j<NY; j++)
			{
				coarse_j = getCoarseID(NY, j);
				fine_j = getFineID(coarse_j, j);

				for (int64_t i=0; i<NX; i++)
				{
					coarse_i = getCoarseID(NX, i);
					fine_i = getFineID(coarse_i, i);
					// TODO idebmap

					// int stat1 = delphi->status[k][j][i];
					// int stat1 = delphi->STATUSMAP(i,j,k,NX,NY);
					// int stat2 = surf->delphi->status[k][j][i];
					// int stat2 = surf->delphi->STATUSMAP(i,j,k,NX,NY);

					int stat1, stat2;

					int64_t coarse_index = coarse_k*(coarse_ny*coarse_nx) + coarse_j*coarse_nx + coarse_i;
					int64_t fine_index = getUnrolledFineID(fine_i, fine_j, fine_k);

					if (delphi->bilevel_status[coarse_index] == NULL)
						stat1 = STATUS_POINT_TEMPORARY_OUT;
					else
						stat1 = delphi->bilevel_status[coarse_index][fine_index];

					if (surf->delphi->bilevel_status[coarse_index] == NULL)
						stat2 = STATUS_POINT_TEMPORARY_OUT;
					else
						stat2 = surf->delphi->bilevel_status[coarse_index][fine_index];

					// if S1 is not out and and S2 is out then that's the pocket
					if ((stat1 != STATUS_POINT_TEMPORARY_OUT && stat1 != STATUS_POINT_OUT) &&
						(stat2 == STATUS_POINT_TEMPORARY_OUT || stat2 == STATUS_POINT_OUT))
					{
						// mark as outside temporary, such that cavity detection can work on it
						// delphi->STATUSMAP(i,j,k,NX,NY) = STATUS_POINT_TEMPORARY_OUT;
						// mark as first cavity such that it can be directly filtered
						// delphi->STATUSMAP(i,j,k,NX,NY) = STATUS_FIRST_CAV;
						if (delphi->bilevel_status[coarse_index] == NULL)
						{
							delphi->bilevel_status[coarse_index] = allocateBilevelMinigridCells<int>(STATUS_POINT_TEMPORARY_OUT);
						}
						delphi->bilevel_status[coarse_index][fine_index] = STATUS_FIRST_CAV;

						countCubes++;

						int *v = allocateVector<int>(3);
						v[0] = i;
						v[1] = j;
						v[2] = k;
						// store the point as a cavity
						vv->push_back(v);

						// switch surrounding points

						if (accurateTriangulation)
						{
							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								verticesInsidenessMap[k][j][i] = true;
								verticesInsidenessMap[k][j][i+1] = true;
								verticesInsidenessMap[k][j+1][i] = true;
								verticesInsidenessMap[k][j+1][i+1] = true;
								verticesInsidenessMap[k+1][j][i] = true;
								verticesInsidenessMap[k+1][j][i+1] = true;
								verticesInsidenessMap[k+1][j+1][i] = true;
								verticesInsidenessMap[k+1][j+1][i+1] = true;
								/*
								write3DVector<bool>(verticesInsidenessMap,true,i,j,k,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i+1,j,k,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i,j+1,k,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i+1,j+1,k,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i,j,k+1,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i+1,j,k+1,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i,j+1,k+1,NX,NY,NZ);
								write3DVector<bool>(verticesInsidenessMap,true,i+1,j+1,k+1,NX,NY,NZ);
								*/
							}
							else
							#endif
							{
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i,j,k,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i+1,j,k,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i,j+1,k,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i+1,j+1,k,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i,j,k+1,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i+1,j,k+1,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i,j+1,k+1,NX,NY,NZ);
								write32xCompressedGrid(compressed_verticesInsidenessMap,true,i+1,j+1,k+1,NX,NY,NZ);
							}
						}
					}
					else
					{
						// delphi->status[k][j][i] = STATUS_POINT_INSIDE;
						// delphi->STATUSMAP(i,j,k,NX,NY) = STATUS_POINT_INSIDE;
						if (delphi->bilevel_status[coarse_index] == NULL)
						{
							delphi->bilevel_status[coarse_index] = allocateBilevelMinigridCells<int>(STATUS_POINT_TEMPORARY_OUT);
						}
						delphi->bilevel_status[coarse_index][fine_index] = STATUS_POINT_INSIDE;

						if (accurateTriangulation)
						{
							#if !defined(USE_COMPRESSED_GRIDS)
							if (!optimizeGrids)
							{
								verticesInsidenessMap[k][j][i] = false;
								// write3DVector<bool>(verticesInsidenessMap,false,i,j,k,NX,NY,NZ);

								if ((i+1)<NX)
									verticesInsidenessMap[k][j][i+1] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i+1,j,k,NX,NY,NZ);

								if ((j+1)<NY)
									verticesInsidenessMap[k][j+1][i] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i,j+1,k,NX,NY,NZ);

								if ((i+1)<NX && (j+1)<NY)
									verticesInsidenessMap[k][j+1][i+1] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i+1,j+1,k,NX,NY,NZ);

								if ((k+1)<NZ)
									verticesInsidenessMap[k+1][j][i] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i,j,k+1,NX,NY,NZ);

								if ((i+1)<NX && (k+1)<NZ)
									verticesInsidenessMap[k+1][j][i+1] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i+1,j,k+1,NX,NY,NZ);

								if ((j+1)<NY && (k+1)<NZ)
									verticesInsidenessMap[k+1][j+1][i] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i,j+1,k+1,NX,NY,NZ);

								if ((i+1)<NX && (j+1)<NY && (k+1)<NZ)
									verticesInsidenessMap[k+1][j+1][i+1] = false;
									// write3DVector<bool>(verticesInsidenessMap,false,i+1,j+1,k+1,NX,NY,NZ);
							}
							else
							#endif
							{
								write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,j,k,NX,NY,NZ);

								if ((i+1)<NX)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i+1,j,k,NX,NY,NZ);

								if ((j+1)<NY)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,j+1,k,NX,NY,NZ);

								if ((i+1)<NX && (j+1)<NY)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i+1,j+1,k,NX,NY,NZ);

								if ((k+1)<NZ)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,j,k+1,NX,NY,NZ);

								if ((i+1)<NX && (k+1)<NZ)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i+1,j,k+1,NX,NY,NZ);

								if ((j+1)<NY && (k+1)<NZ)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i,j+1,k+1,NX,NY,NZ);

								if ((i+1)<NX && (j+1)<NY && (k+1)<NZ)
									write32xCompressedGrid(compressed_verticesInsidenessMap,false,i+1,j+1,k+1,NX,NY,NZ);
							}
						}
					}
				}
			}
		}
		delphi->cavitiesSize.push_back((delphi->side)*(delphi->side)*(delphi->side)*countCubes);
	}

	// TODO aggregate analytical triangulation points
	// TODO aggregate boundary grid points
	return true;
}


void Surface::tri2Balls()
{
	if (vertList.size() == 0)
	{
		cout << endl << WARN << "Triangulation not available, cannot build set of balls approximation!";
		cout << endl << REMARK << "Please enable triangulation by 'Triangulation = true' ";
		return;
	}
	
	if (!optimizeGrids)
	{
		if (delphi->status == NULL)
		{
			cout << endl << WARN << "Status map not built, cannot build set of balls approximation!";
			cout << endl << REMARK << "Please enable status map by 'Build_status_map = true'";
			return;
		}
	}
	else
	{
		if (delphi->bilevel_status == NULL)
		{
			cout << endl << WARN << "bilevel_status map not built, cannot build set of balls approximation!";
			cout << endl << REMARK << "Please enable bilevel status map by 'Build_status_map = true'";
			return;
		}
	}

	int64_t NX = delphi->nx;
	int64_t NY = delphi->ny;
	int64_t NZ = delphi->nz;

	#ifndef ENABLE_CGAL
	cout << endl << WARN << "CGAL is required by tri2Balls function!";
	cout << endl << REMARK << "Recompile enabling CGAL support";
	return;
	#else

	#if !defined(NO_CGAL_PATCHING)
	vector<_Weighted_point> l;
	#else
	// nopatch
	vector<std::pair<_Weighted_point,int>> l;
	#endif

	double max_x=-1e6, min_x=1e6, max_y=-1e6, min_y=1e6, max_z=-1e6, min_z=1e6;

	double x, y, z;

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int numVertices = (int)vertList.size();
	#else
	int numVertices = (int)(vertList.size() / 3.);
	#endif

	l.reserve(numVertices+6);

	// cout << endl << INFO << "Converting triangulation to set of balls...";
			
	for (unsigned int i=0; i<numVertices; i++)
	{
		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		x = vertList[i][0];
		y = vertList[i][1];
		z = vertList[i][2];
		#else
		x = vertList[ i*3+0 ];
		y = vertList[ i*3+1 ];
		z = vertList[ i*3+2 ];
		#endif

		#if !defined(NO_CGAL_PATCHING)
		l.push_back(_Weighted_point(_Point3(x,y,z), 0., i));
		#else
		// nopatch
		l.push_back(std::make_pair(_Weighted_point(_Point3(x,y,z), 0.), i));
		#endif

		max_x = max(max_x,x);
		max_y = max(max_y,y);
		max_z = max(max_z,z);

		min_x = min(min_x,x);
		min_y = min(min_y,y);
		min_z = min(min_z,z);
	}
	
	// Regular Triangulation object 
	_Rt rT;

	double mid_x = (max_x+min_x)*0.5;
	double mid_y = (max_y+min_y)*0.5;
	double mid_z = (max_z+min_z)*0.5;

	min_x -= fabs(mid_x-min_x)*2.;
	min_y -= fabs(mid_y-min_y)*2.;
	min_z -= fabs(mid_z-min_z)*2.;

	max_x += fabs(mid_x-max_x)*2.;
	max_y += fabs(mid_y-max_y)*2.;
	max_z += fabs(mid_z-max_z)*2.;

	// add bounding box using -1 weights to let easy detection of these virtual points
	#if !defined(NO_CGAL_PATCHING)
	l.push_back(_Weighted_point(_Point3(min_x,mid_y,mid_z),-1,numVertices));
	l.push_back(_Weighted_point(_Point3(max_x,mid_y,mid_z),-1,numVertices+1));
	l.push_back(_Weighted_point(_Point3(mid_x,min_y,mid_z),-1,numVertices+2));
	l.push_back(_Weighted_point(_Point3(mid_x,max_y,mid_z),-1,numVertices+3));
	l.push_back(_Weighted_point(_Point3(mid_x,mid_y,min_z),-1,numVertices+4));
	l.push_back(_Weighted_point(_Point3(mid_x,mid_y,max_z),-1,numVertices+5));
	#else
    // nopatch
	l.push_back(std::make_pair(_Weighted_point(_Point3(min_x,mid_y,mid_z),-1),numVertices));
	l.push_back(std::make_pair(_Weighted_point(_Point3(max_x,mid_y,mid_z),-1),numVertices+1));
	l.push_back(std::make_pair(_Weighted_point(_Point3(mid_x,min_y,mid_z),-1),numVertices+2));
	l.push_back(std::make_pair(_Weighted_point(_Point3(mid_x,max_y,mid_z),-1),numVertices+3));
	l.push_back(std::make_pair(_Weighted_point(_Point3(mid_x,mid_y,min_z),-1),numVertices+4));
	l.push_back(std::make_pair(_Weighted_point(_Point3(mid_x,mid_y,max_z),-1),numVertices+5));
	#endif

	cout << endl << INFO << "Computing triangulation...";
	time_t rt_start, rt_end;
    time(&rt_start);

	#ifdef CGAL_LINKED_WITH_TBB
	// Construct the locking data-structure, using the bounding-box of the points
	_Rt::Lock_data_structure locking_ds(CGAL::Bbox_3(min_x,min_y, min_z, max_x, max_y, max_z), 50);
	rT.set_lock_data_structure(&locking_ds);
	rT.insert (l.begin(), l.end());
	#else
	rT.insert (l.begin(), l.end());
	#endif

	time(&rt_end);
    double insertion_time = difftime(rt_end, rt_start);

	assert( rT.is_valid() );
	assert( rT.dimension() == 3 );

	// cout << "ok!";
	// cout << endl << INFO << "Computing voronoi points...";

	for (_Finite_Cells_Iterator fcit = rT.finite_cells_begin(); fcit!=rT.finite_cells_end(); fcit++)
	{
		const _Point3 &p = rT.geom_traits().construct_weighted_circumcenter_3_object()(fcit->vertex(0)->point(),fcit->vertex(1)->point(),fcit->vertex(2)->point(),fcit->vertex(3)->point());
		VorPoint *&mc = (VorPoint*&)fcit->info();
		mc = new VorPoint();
		mc->vor[0] = p.x();
		mc->vor[1] = p.y();
		mc->vor[2] = p.z();
	}

	// cout << "ok!";
	// cout << endl << INFO << "Collecting polar balls...";
	
	vector<_Cell_handle> cells;
	cells.reserve(1000);

	vector<double*> polarBalls;
	vector<double> polarDist;

	for(_Finite_Vertex_Iterator fvit = rT.finite_vertices_begin(); fvit != rT.finite_vertices_end(); fvit++)
	{	  
		const _Weight& w0 = fvit->point().weight();
		
		// skip bounding box points
		if (w0 == -1)
			continue;
	  
		if (rT.is_infinite(fvit->cell()))
			continue;

		// current point id around which we are moving
		#if !defined(NO_CGAL_PATCHING)
		const int &currentId = fvit->point().index();
		#else
		// nopatch
		const int &currentId = fvit->info();
		#endif
		
		cells.clear();
		rT.incident_cells(fvit,std::back_inserter(cells));

		bool infiniteCell = false;

		// start check if all tetraedra are feasible
		vector<_Cell_handle>::iterator it;
		for (it=cells.begin(); it!=cells.end(); it++)
		{
			if (rT.is_infinite((*it)))
			{
				infiniteCell = true;
				break;
			}
		}
	 
		if (infiniteCell)
			continue;

		double ref[3], winner[3], maxDist=0, dist;
		VorPoint *winner_vp;

		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		ref[0] = vertList[currentId][0];
		ref[1] = vertList[currentId][1];
		ref[2] = vertList[currentId][2];
		#else
		ref[0] = vertList[ currentId*3+0 ];
		ref[1] = vertList[ currentId*3+1 ];
		ref[2] = vertList[ currentId*3+2 ];
		#endif
		
		// compute the internal farthest ball (the inner polar ball)
		for (it=cells.begin(); it!=cells.end(); it++)
		{
			const VorPoint *vp = (*it)->info();
			DIST2(dist,vp->vor,ref);

			// check it is internal by using the status map
			int64_t ix = (int64_t)rintp((vp->vor[0]-delphi->xmin)*delphi->scale);
			int64_t iy = (int64_t)rintp((vp->vor[1]-delphi->ymin)*delphi->scale);
			int64_t iz = (int64_t)rintp((vp->vor[2]-delphi->zmin)*delphi->scale);

			if (ix<0 || iy<0 || iz<0 || ix>=delphi->nx || iy>=delphi->ny || iz>=delphi->nz)
				continue;

			if (!optimizeGrids)
			{
				/*
				if ((delphi->STATUSMAP(ix,iy,iz,NX,NY) == STATUS_POINT_INSIDE) &&
					((ix+1)<NX && delphi->STATUSMAP((ix+1),(iy),(iz),NX,NY) == STATUS_POINT_INSIDE) &&
					((iy+1)<NY && delphi->STATUSMAP(ix,(iy+1),iz,NX,NY) == STATUS_POINT_INSIDE) &&
					((iz+1)<NZ && delphi->STATUSMAP(ix,iy,(iz+1),NX,NY) == STATUS_POINT_INSIDE) &&
					((ix-1)>=0 && delphi->STATUSMAP((ix-1),iy,iz,NX,NY) == STATUS_POINT_INSIDE) &&
					((iy-1)>=0 && delphi->STATUSMAP(ix,(iy-1),iz,NX,NY) == STATUS_POINT_INSIDE) &&
					((iz-1)>=0 && delphi->STATUSMAP(ix,iy,(iz-1),NX,NY) == STATUS_POINT_INSIDE))
				*/
				if ((read3DVector<int>(delphi->status,ix,iy,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((ix+1)<NX && read3DVector<int>(delphi->status,ix+1,iy,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((iy+1)<NY && read3DVector<int>(delphi->status,ix,iy+1,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((iz+1)<NZ && read3DVector<int>(delphi->status,ix,iy,iz+1,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((ix-1)>=0 && read3DVector<int>(delphi->status,ix-1,iy,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((iy-1)>=0 && read3DVector<int>(delphi->status,ix,iy-1,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((iz-1)>=0 && read3DVector<int>(delphi->status,ix,iy,iz-1,NX,NY,NZ) == STATUS_POINT_INSIDE))
				{
					if (dist > maxDist)
					{
						maxDist = dist;
						ASSIGN(winner, vp->vor);
						winner_vp = (VorPoint*)vp;
					}
				}
			}
			else
			{
				if ((readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((ix+1)<NX && readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix+1,iy,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((iy+1)<NY && readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy+1,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((iz+1)<NZ && readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz+1,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((ix-1)>=0 && readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix-1,iy,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((iy-1)>=0 && readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy-1,iz,NX,NY,NZ) == STATUS_POINT_INSIDE) &&
					((iz-1)>=0 && readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz-1,NX,NY,NZ) == STATUS_POINT_INSIDE))
				{
					if (dist > maxDist)
					{
						maxDist = dist;
						ASSIGN(winner, vp->vor);
						winner_vp = (VorPoint*)vp;
					}
				}
			}
		}

		if (maxDist == 0)
			continue;

		if (winner_vp->visited == true)
			continue;
		else
			winner_vp->visited = true;

		double *p = allocateVector<double>(3);
		polarBalls.push_back(p);
		ASSIGN(p,winner);
		polarDist.push_back(sqrt(maxDist));
	}
  
	// cout << "ok!";

	char fullName[1024];
	sprintf(fullName, "%spolar.txt", rootFile.c_str());
	FILE *fp = fopen(fullName, "w");
	
	sprintf(fullName, "%spolar.pqr", rootFile.c_str());
	FILE *fp2 = fopen(fullName,"w");
	for (unsigned int i=0; i<polarBalls.size(); i++)
	{
		fprintf(fp, "%f %f %f %f \n", polarBalls[i][0], polarBalls[i][1], polarBalls[i][2], polarDist[i]);
		fprintf(fp2, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f%8.4f%8.4f\n", i%100000, "C","NSR", 'A', 1, polarBalls[i][0], polarBalls[i][1], polarBalls[i][2], 0., polarDist[i]);
		deleteVector<double>(polarBalls[i]);
	}

	fclose(fp);
	fclose(fp2);

	#endif // ENABLE_CGAL
}


/*
// This function is replaced by the below optimised version

void Surface::smoothSurface(bool outputMesh, bool buildAtomsMapHere, char *fn, bool revert)
{
	auto chrono_start = chrono::high_resolution_clock::now();

	#define MAX_NEIGHBOURS 20

	vector<VERTEX_TYPE*> tempVertices;
	vector<VERTEX_TYPE*> tempNormals;

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int nv = (int)vertList.size();
	#else
	int nv = (int)(vertList.size() / 3.);
	#endif
	int nt = (int)(triList.size() / 3.);

	tempVertices.resize(nv);
	
	if (computeNormals)
		tempNormals.resize(nv);

	for (int i=0; i<nv; i++)
		tempVertices[i] = allocateVector<VERTEX_TYPE>(3);

	if (computeNormals)
		for (int i=0; i<nv; i++)
			tempNormals[i] = allocateVector<VERTEX_TYPE>(3);

	int **vertdeg = allocateVector<int*>(nv);

	for (int i=0; i<nv; i++)
	{
		vertdeg[i] = allocateVector<int>(MAX_NEIGHBOURS + 1);
		vertdeg[i][0] = 0;
	}

	bool flagvert;

	for (int i=0; i<nt; i++)
	{
		// first vertex
		flagvert = true;

		for (int j=0; j < vertdeg[ triList[i*3+0] ][0]; j++)
		{
			if (triList[i*3+1] == vertdeg[ triList[i*3+0] ][j+1])
			{
				flagvert = false;
				break;
			}
		}
		if (flagvert)
		{
			vertdeg[ triList[i*3+0] ][ ++vertdeg[triList[i*3+0]][0] ] = triList[i*3+1];
		}
		flagvert = true;

		for (int j=0; j < vertdeg[ triList[i*3+0] ][0]; j++)
		{
			if (triList[i*3+2] == vertdeg[ triList[i*3+0] ][j+1])
			{
				flagvert = false;
				break;
			}
		}
		if (flagvert)
		{
			vertdeg[ triList[i*3+0] ][ ++vertdeg[triList[i*3+0]][0] ] = triList[i*3+2];
		}
		// second vertex
		flagvert = true;

		for (int j=0; j < vertdeg[ triList[i*3+1] ][0]; j++)
		{
			if (triList[i*3+0] == vertdeg[ triList[i*3+1] ][j+1])
			{
				flagvert = false;
				break;
			}
		}
		if (flagvert)
		{
			vertdeg[ triList[i*3+1] ][ ++vertdeg[triList[i*3+1]][0] ] = triList[i*3+0];
		}
		flagvert = true;

		for (int j=0; j < vertdeg[ triList[i*3+1] ][0]; j++)
		{
			if (triList[i*3+2] == vertdeg[ triList[i*3+1] ][j+1])
			{
				flagvert = false;
				break;
			}
		}
		if (flagvert)
		{
			vertdeg[ triList[i*3+1] ][ ++vertdeg[triList[i*3+1]][0] ] = triList[i*3+2];
		}
		// third vertex
		flagvert = true;

		for (int j=0; j < vertdeg[ triList[i*3+2] ][0]; j++)
		{
			if (triList[i*3+0] == vertdeg[ triList[i*3+2] ][j+1])
			{
				flagvert = false;
				break;
			}
		}
		if (flagvert)
		{
			vertdeg[ triList[i*3+2] ][ ++vertdeg[triList[i*3+2]][0] ] = triList[i*3+0];
		}
		flagvert = true;

		for (int j=0; j < vertdeg[ triList[i*3+2] ][0]; j++)
		{
			if (triList[i*3+1] == vertdeg[ triList[i*3+2] ][j+1])
			{
				flagvert = false;
				break;
			}
		}
		if (flagvert)
		{
			vertdeg[ triList[i*3+2] ][ ++vertdeg[triList[i*3+2]][0] ] = triList[i*3+1];
		}
	}

	for (int i=0; i<nv; i++)
	{
		tempVertices[i][0] = 0;
		tempVertices[i][1] = 0;
		tempVertices[i][2] = 0;

		int div = vertdeg[i][0];

		// forward
		// primary neighbours
		for (int j=0; j<vertdeg[i][0]; j++)
		{
			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			tempVertices[i][0] += vertList[ vertdeg[i][j+1] ][0];
			tempVertices[i][1] += vertList[ vertdeg[i][j+1] ][1];
			tempVertices[i][2] += vertList[ vertdeg[i][j+1] ][2];
			#else
			tempVertices[i][0] += vertList[ vertdeg[i][j+1]*3+0 ];
			tempVertices[i][1] += vertList[ vertdeg[i][j+1]*3+1 ];
			tempVertices[i][2] += vertList[ vertdeg[i][j+1]*3+2 ];
			#endif
		}
		tempVertices[i][0] /= div;
		tempVertices[i][1] /= div;
		tempVertices[i][2] /= div;
		
		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		tempVertices[i][0] = 0.5*(vertList[i][0] + tempVertices[i][0]);
		tempVertices[i][1] = 0.5*(vertList[i][1] + tempVertices[i][1]);
		tempVertices[i][2] = 0.5*(vertList[i][2] + tempVertices[i][2]);
		#else
		tempVertices[i][0] = 0.5*(vertList[ i*3+0 ] + tempVertices[i][0]);
		tempVertices[i][1] = 0.5*(vertList[ i*3+1 ] + tempVertices[i][1]);
		tempVertices[i][2] = 0.5*(vertList[ i*3+2 ] + tempVertices[i][2]);
		#endif

		if (computeNormals)
		{
			tempNormals[i][0] = 0;
			tempNormals[i][1] = 0;
			tempNormals[i][2] = 0;
			// forward
			// primary neighbours
			for (int j=0; j<vertdeg[i][0]; j++)
			{
				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				tempNormals[i][0] += normalsList[ vertdeg[i][j+1] ][0];
				tempNormals[i][1] += normalsList[ vertdeg[i][j+1] ][1];
				tempNormals[i][2] += normalsList[ vertdeg[i][j+1] ][2];
				#else
				tempNormals[i][0] += normalsList[ vertdeg[i][j+1]*3+0 ];
				tempNormals[i][1] += normalsList[ vertdeg[i][j+1]*3+1 ];
				tempNormals[i][2] += normalsList[ vertdeg[i][j+1]*3+2 ];
				#endif
			}
			tempNormals[i][0] /= div;
			tempNormals[i][1] /= div;
			tempNormals[i][2] /= div;

			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			tempNormals[i][0] = 0.5*(normalsList[i][0] + tempNormals[i][0]);
			tempNormals[i][1] = 0.5*(normalsList[i][1] + tempNormals[i][1]);
			tempNormals[i][2] = 0.5*(normalsList[i][2] + tempNormals[i][2]);
			#else
			tempNormals[i][0] = 0.5*(normalsList[ i*3+0 ] + tempNormals[i][0]);
			tempNormals[i][1] = 0.5*(normalsList[ i*3+1 ] + tempNormals[i][1]);
			tempNormals[i][2] = 0.5*(normalsList[ i*3+2 ] + tempNormals[i][2]);
			#endif

			VERTEX_TYPE tt;
			NORMALIZE(tempNormals[i],tt);
		}
	}

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	// save all
	for (int i=0; i<nv; i++)
	{
		vertList[i][0] = tempVertices[i][0];
		vertList[i][1] = tempVertices[i][1];
		vertList[i][2] = tempVertices[i][2];
	}
	if (computeNormals)
	{
		for (int i=0; i<nv; i++)
		{
			normalsList[i][0] = tempNormals[i][0];
			normalsList[i][1] = tempNormals[i][1];
			normalsList[i][2] = tempNormals[i][2];
		}
	}
	#else
	// save all
	for (int i=0; i<nv; i++)
	{
		vertList[ i*3+0 ] = tempVertices[i][0];
		vertList[ i*3+1 ] = tempVertices[i][1];
		vertList[ i*3+2 ] = tempVertices[i][2];
	}
	if (computeNormals)
	{
		for (int i=0; i<nv; i++)
		{
			normalsList[ i*3+0 ] = tempNormals[i][0];
			normalsList[ i*3+1 ] = tempNormals[i][1];
			normalsList[ i*3+2 ] = tempNormals[i][2];
		}
	}
	#endif

	// delete all
	for (int i=0; i<nv; i++)
		deleteVector<VERTEX_TYPE>(tempVertices[i]);

	if (computeNormals)
		for (int i=0; i<nv; i++)
			deleteVector<VERTEX_TYPE>(tempNormals[i]);

	for (int i=0; i<nv; i++)
		deleteVector<int>(vertdeg[i]);
	deleteVector<int*>(vertdeg);

	auto chrono_end = chrono::high_resolution_clock::now();
	chrono::duration<double> smoothing_time = chrono_end - chrono_start;
	cout << endl << INFO << "Smoothing time ";
	printf ("%.4e [s]", smoothing_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	if (!conf.parallelPocketLoop)
	{
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after smoothing is " << current_mem_in_MB << " MB";
	}
	#endif

	// check atoms flag
	if (buildAtomsMapHere && vertexAtomsMapFlag)
	{
		auto chrono_start = chrono::high_resolution_clock::now();

		if (vertexAtomsMap != NULL)
			deleteVector<int>(vertexAtomsMap);

		vertexAtomsMap = allocateVector<int>(nv);
		buildAtomsMap();

		cout << endl << INFO << "Connecting vertices to atoms...";
		cout.flush();

		for (int i=0; i<nv; i++)
		{
			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			vdwAccessible(vertList[i], vertexAtomsMap[i]);
			#else
			vdwAccessible(&vertList[i*3], vertexAtomsMap[i]);
			#endif

			if (vertexAtomsMap[i] == -1)
			{
				cout << endl << WARN << "Cannot detect nearest atom for vertex " << i;
			}
		}
		cout << "ok!";
		cout.flush();

		auto chrono_end = chrono::high_resolution_clock::now();
		chrono::duration<double> other_time = chrono_end - chrono_start;
		cout << endl << INFO << "Vertex atom mapping time is ";
		printf ("%.4e [s]", other_time.count());

		#if !defined(AVOID_MEM_CHECKS)
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after vertex atom mapping is " << current_mem_in_MB << " MB";
		#endif
	}

	if (outputMesh)
	{
		auto chrono_start = chrono::high_resolution_clock::now();
		int format = deduceFormat();
		bool f = saveMesh(format,revert,fn,vertList,triList,normalsList);

		if (!f)
			cout << endl << ERR << "Problems in saving the mesh!";
		else
		{
			auto chrono_end = chrono::high_resolution_clock::now();
			chrono::duration<double> saveMesh_time = chrono_end - chrono_start;
			cout << endl << INFO << "Outputting mesh time (in smoothSurface()) ";
			printf ("%.4e [s]", saveMesh_time.count());
		}
	}
}
*/



// This version can be called by multiple threads and circumevents scattered stores which would normally
// require atomics with private index-based writes

void Surface::getTempVerticesAndNormals (VERTEX_TYPE tempVertices[], VERTEX_TYPE tempNormals[],
										 int numNeighbours[], int thread_id, int num_threads)
{
	int nt = (int)(triList.size() / 3.);


	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	for (int i=0; i<nt; i++)
	{
		int v1 = triList[ i*3+0 ];
		int v2 = triList[ i*3+1 ];
		int v3 = triList[ i*3+2 ];

		if (v1 % num_threads == thread_id)
		{
			tempVertices[ v1*3+0 ] += vertList[v2][0] + vertList[v3][0];
			tempVertices[ v1*3+1 ] += vertList[v2][1] + vertList[v3][1];
			tempVertices[ v1*3+2 ] += vertList[v2][2] + vertList[v3][2];

			if (computeNormals)
			{
				tempNormals[ v1*3+0 ] += normalsList[v2][0] + normalsList[v3][0];
				tempNormals[ v1*3+1 ] += normalsList[v2][1] + normalsList[v3][1];
				tempNormals[ v1*3+2 ] += normalsList[v2][2] + normalsList[v3][2];
			}
			// Note: every edge is shared by 2 triangles which then contribute twice to each
			// vertex, but also the vertex and normal data sums are overestimated by 2X
			numNeighbours[ v1 ] += 2;
		}
		if (v2 % num_threads == thread_id)
		{
			tempVertices[ v2*3+0 ] += vertList[v1][0] + vertList[v3][0];
			tempVertices[ v2*3+1 ] += vertList[v1][1] + vertList[v3][1];
			tempVertices[ v2*3+2 ] += vertList[v1][2] + vertList[v3][2];

			if (computeNormals)
			{
				tempNormals[ v2*3+0 ] += normalsList[v1][0] + normalsList[v3][0];
				tempNormals[ v2*3+1 ] += normalsList[v1][1] + normalsList[v3][1];
				tempNormals[ v2*3+2 ] += normalsList[v1][2] + normalsList[v3][2];
			}
			numNeighbours[ v2 ] += 2;
		}
		if (v3 % num_threads == thread_id)
		{
			tempVertices[ v3*3+0 ] += vertList[v1][0] + vertList[v2][0];
			tempVertices[ v3*3+1 ] += vertList[v1][1] + vertList[v2][1];
			tempVertices[ v3*3+2 ] += vertList[v1][2] + vertList[v2][2];

			if (computeNormals)
			{
				tempNormals[ v3*3+0 ] += normalsList[v1][0] + normalsList[v2][0];
				tempNormals[ v3*3+1 ] += normalsList[v1][1] + normalsList[v2][1];
				tempNormals[ v3*3+2 ] += normalsList[v1][2] + normalsList[v2][2];
			}
			numNeighbours[ v3 ] += 2;
		}
	}
	#else // USE_OPTIMIZED_VERTICES_BUFFERING
	for (int i=0; i<nt; i++)
	{
		int v1 = triList[ i*3+0 ];
		int v2 = triList[ i*3+1 ];
		int v3 = triList[ i*3+2 ];

		if (v1 % num_threads == thread_id)
		{
			tempVertices[ v1*3+0 ] += vertList[ v2*3+0 ] + vertList[ v3*3+0 ];
			tempVertices[ v1*3+1 ] += vertList[ v2*3+1 ] + vertList[ v3*3+1 ];
			tempVertices[ v1*3+2 ] += vertList[ v2*3+2 ] + vertList[ v3*3+2 ];

			if (computeNormals)
			{
				tempNormals[ v1*3+0 ] += normalsList[ v2*3+0 ] + normalsList[ v3*3+0 ];
				tempNormals[ v1*3+1 ] += normalsList[ v2*3+1 ] + normalsList[ v3*3+1 ];
				tempNormals[ v1*3+2 ] += normalsList[ v2*3+2 ] + normalsList[ v3*3+2 ];
			}
			// Note: every edge is shared by 2 triangles which then contribute twice to each
			// vertex, but also the vertex and normal data sums are overestimated by 2X
			numNeighbours[ v1 ] += 2;
		}
		if (v2 % num_threads == thread_id)
		{
			tempVertices[ v2*3+0 ] += vertList[ v1*3+0 ] + vertList[ v3*3+0 ];
			tempVertices[ v2*3+1 ] += vertList[ v1*3+1 ] + vertList[ v3*3+1 ];
			tempVertices[ v2*3+2 ] += vertList[ v1*3+2 ] + vertList[ v3*3+2 ];

			if (computeNormals)
			{
				tempNormals[ v2*3+0 ] += normalsList[ v1*3+0 ] + normalsList[ v3*3+0 ];
				tempNormals[ v2*3+1 ] += normalsList[ v1*3+1 ] + normalsList[ v3*3+1 ];
				tempNormals[ v2*3+2 ] += normalsList[ v1*3+2 ] + normalsList[ v3*3+2 ];
			}
			numNeighbours[ v2 ] += 2;
		}
		if (v3 % num_threads == thread_id)
		{
			tempVertices[ v3*3+0 ] += vertList[ v1*3+0 ] + vertList[ v2*3+0 ];
			tempVertices[ v3*3+1 ] += vertList[ v1*3+1 ] + vertList[ v2*3+1 ];
			tempVertices[ v3*3+2 ] += vertList[ v1*3+2 ] + vertList[ v2*3+2 ];

			if (computeNormals)
			{
				tempNormals[ v3*3+0 ] += normalsList[ v1*3+0 ] + normalsList[ v2*3+0 ];
				tempNormals[ v3*3+1 ] += normalsList[ v1*3+1 ] + normalsList[ v2*3+1 ];
				tempNormals[ v3*3+2 ] += normalsList[ v1*3+2 ] + normalsList[ v2*3+2 ];
			}
			numNeighbours[ v3 ] += 2;
		}
	}
	#endif // USE_OPTIMIZED_VERTICES_BUFFERING
}


void Surface::getTempVerticesAndNormals (VERTEX_TYPE tempVertices[], VERTEX_TYPE tempNormals[],
										 int numNeighbours[])
{
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	for (int i=0; i < (int)(triList.size() / 3.); i++)
	{
		int v1 = triList[ i*3+0 ];
		int v2 = triList[ i*3+1 ];
		int v3 = triList[ i*3+2 ];

		tempVertices[ v1*3+0 ] += vertList[v2][0] + vertList[v3][0];
		tempVertices[ v1*3+1 ] += vertList[v2][1] + vertList[v3][1];
		tempVertices[ v1*3+2 ] += vertList[v2][2] + vertList[v3][2];
		tempVertices[ v2*3+0 ] += vertList[v1][0] + vertList[v3][0];
		tempVertices[ v2*3+1 ] += vertList[v1][1] + vertList[v3][1];
		tempVertices[ v2*3+2 ] += vertList[v1][2] + vertList[v3][2];
		tempVertices[ v3*3+0 ] += vertList[v1][0] + vertList[v2][0];
		tempVertices[ v3*3+1 ] += vertList[v1][1] + vertList[v2][1];
		tempVertices[ v3*3+2 ] += vertList[v1][2] + vertList[v2][2];

		if (computeNormals)
		{
			tempNormals[ v1*3+0 ] += normalsList[v2][0] + normalsList[v3][0];
			tempNormals[ v1*3+1 ] += normalsList[v2][1] + normalsList[v3][1];
			tempNormals[ v1*3+2 ] += normalsList[v2][2] + normalsList[v3][2];
			tempNormals[ v2*3+0 ] += normalsList[v1][0] + normalsList[v3][0];
			tempNormals[ v2*3+1 ] += normalsList[v1][1] + normalsList[v3][1];
			tempNormals[ v2*3+2 ] += normalsList[v1][2] + normalsList[v3][2];
			tempNormals[ v3*3+0 ] += normalsList[v1][0] + normalsList[v2][0];
			tempNormals[ v3*3+1 ] += normalsList[v1][1] + normalsList[v2][1];
			tempNormals[ v3*3+2 ] += normalsList[v1][2] + normalsList[v2][2];
		}
		// Note: every edge is shared by 2 triangles which then contribute twice to each
		// vertex, but also the vertex and normal data sums are overestimated by 2X
		numNeighbours[ v1 ] += 2;
		numNeighbours[ v2 ] += 2;
		numNeighbours[ v3 ] += 2;
	}
	#else // USE_OPTIMIZED_VERTICES_BUFFERING
	for (int i=0; i < (int)(triList.size() / 3.); i++)
	{
		int v1 = triList[ i*3+0 ];
		int v2 = triList[ i*3+1 ];
		int v3 = triList[ i*3+2 ];

		tempVertices[ v1*3+0 ] += vertList[ v2*3+0 ] + vertList[ v3*3+0 ];
		tempVertices[ v1*3+1 ] += vertList[ v2*3+1 ] + vertList[ v3*3+1 ];
		tempVertices[ v1*3+2 ] += vertList[ v2*3+2 ] + vertList[ v3*3+2 ];
		tempVertices[ v2*3+0 ] += vertList[ v1*3+0 ] + vertList[ v3*3+0 ];
		tempVertices[ v2*3+1 ] += vertList[ v1*3+1 ] + vertList[ v3*3+1 ];
		tempVertices[ v2*3+2 ] += vertList[ v1*3+2 ] + vertList[ v3*3+2 ];
		tempVertices[ v3*3+0 ] += vertList[ v1*3+0 ] + vertList[ v2*3+0 ];
		tempVertices[ v3*3+1 ] += vertList[ v1*3+1 ] + vertList[ v2*3+1 ];
		tempVertices[ v3*3+2 ] += vertList[ v1*3+2 ] + vertList[ v2*3+2 ];

		if (computeNormals)
		{
			tempNormals[ v1*3+0 ] += normalsList[ v2*3+0 ] + normalsList[ v3*3+0 ];
			tempNormals[ v1*3+1 ] += normalsList[ v2*3+1 ] + normalsList[ v3*3+1 ];
			tempNormals[ v1*3+2 ] += normalsList[ v2*3+2 ] + normalsList[ v3*3+2 ];
			tempNormals[ v2*3+0 ] += normalsList[ v1*3+0 ] + normalsList[ v3*3+0 ];
			tempNormals[ v2*3+1 ] += normalsList[ v1*3+1 ] + normalsList[ v3*3+1 ];
			tempNormals[ v2*3+2 ] += normalsList[ v1*3+2 ] + normalsList[ v3*3+2 ];
			tempNormals[ v3*3+0 ] += normalsList[ v1*3+0 ] + normalsList[ v2*3+0 ];
			tempNormals[ v3*3+1 ] += normalsList[ v1*3+1 ] + normalsList[ v2*3+1 ];
			tempNormals[ v3*3+2 ] += normalsList[ v1*3+2 ] + normalsList[ v2*3+2 ];
		}
		// Note: every edge is shared by 2 triangles which then contribute twice to each
		// vertex, but also the vertex and normal data sums are overestimated by 2X
		numNeighbours[ v1 ] += 2;
		numNeighbours[ v2 ] += 2;
		numNeighbours[ v3 ] += 2;
	}
	#endif // USE_OPTIMIZED_VERTICES_BUFFERING
}


void Surface::getNewVerticesAndNormals (VERTEX_TYPE tempVertices[], VERTEX_TYPE tempNormals[],
										int numNeighbours[], int thread_id, int num_threads)
{
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int nv = (int)vertList.size();
	#else
	int nv = (int)(vertList.size() / 3.);
	#endif

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	for (int i=thread_id; i < nv; i += num_threads)
	{
		VERTEX_TYPE scale = 1.0 / (VERTEX_TYPE)numNeighbours[i];

		vertList[i][0] = 0.5*(vertList[i][0] + tempVertices[ i*3+0 ]*scale);
		vertList[i][1] = 0.5*(vertList[i][1] + tempVertices[ i*3+1 ]*scale);
		vertList[i][2] = 0.5*(vertList[i][2] + tempVertices[ i*3+2 ]*scale);

		if (computeNormals)
		{
			// a 0.5X front multiplication is not needed because of the below normalisation
			tempNormals[ i*3+0 ] = normalsList[i][0] + tempNormals[ i*3+0 ]*scale;
			tempNormals[ i*3+1 ] = normalsList[i][1] + tempNormals[ i*3+1 ]*scale;
			tempNormals[ i*3+2 ] = normalsList[i][2] + tempNormals[ i*3+2 ]*scale;

			VERTEX_TYPE tt = 1.0 / sqrt(
				tempNormals[i*3+0] * tempNormals[i*3+0] +
				tempNormals[i*3+1] * tempNormals[i*3+1] +
				tempNormals[i*3+2] * tempNormals[i*3+2]);

			normalsList[i][0] = tempNormals[ i*3+0 ] * tt;
			normalsList[i][1] = tempNormals[ i*3+1 ] * tt;
			normalsList[i][2] = tempNormals[ i*3+2 ] * tt;
		}
	}
	#else
	for (int i=thread_id; i < nv; i += num_threads)
	{
		VERTEX_TYPE scale = 1.0 / (VERTEX_TYPE)numNeighbours[i];

		vertList[ i*3+0 ] = 0.5*(vertList[ i*3+0 ] + tempVertices[ i*3+0 ]*scale);
		vertList[ i*3+1 ] = 0.5*(vertList[ i*3+1 ] + tempVertices[ i*3+1 ]*scale);
		vertList[ i*3+2 ] = 0.5*(vertList[ i*3+2 ] + tempVertices[ i*3+2 ]*scale);

		if (computeNormals)
		{
			// a 0.5X front multiplication is not needed because of the below normalisation
			tempNormals[ i*3+0 ] = normalsList[ i*3+0 ] + tempNormals[ i*3+0 ]*scale;
			tempNormals[ i*3+1 ] = normalsList[ i*3+1 ] + tempNormals[ i*3+1 ]*scale;
			tempNormals[ i*3+2 ] = normalsList[ i*3+2 ] + tempNormals[ i*3+2 ]*scale;

			VERTEX_TYPE tt = 1.0 / sqrt(
				tempNormals[i*3+0] * tempNormals[i*3+0] +
				tempNormals[i*3+1] * tempNormals[i*3+1] +
				tempNormals[i*3+2] * tempNormals[i*3+2]);

			normalsList[ i*3+0 ] = tempNormals[ i*3+0 ] * tt;
			normalsList[ i*3+1 ] = tempNormals[ i*3+1 ] * tt;
			normalsList[ i*3+2 ] = tempNormals[ i*3+2 ] * tt;
		}
	}
	#endif
}

/*
// This version may be drastically faster using coherent memory buffering

void Surface::smoothSurface(bool outputMesh, bool buildAtomsMapHere,
							const char *fileName, bool revert)
{
	auto chrono_start = chrono::high_resolution_clock::now();

	VERTEX_TYPE *tempVertices;
	VERTEX_TYPE *tempNormals;

	int *numNeighbours;

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int nv = (int)vertList.size();
	#else
	int nv = (int)(vertList.size() / 3.);
	#endif
	int nt = (int)(triList.size() / 3.);

	tempVertices = allocateVector<VERTEX_TYPE>(nv*3);
	memset(tempVertices, 0, nv*3*sizeof(VERTEX_TYPE));

	if (computeNormals)
	{
		tempNormals = allocateVector<VERTEX_TYPE>(nv*3);
		memset(tempNormals, 0, nv*3*sizeof(VERTEX_TYPE));
	}
	numNeighbours = allocateVector<int>(nv);
	memset(numNeighbours, 0, nv*sizeof(int));


	int num_threads = conf.numThreads;

	#ifdef ENABLE_BOOST_THREADS
	boost::thread_group thdGroup;
	#endif

	// Try to use the below multithreaded code
	getTempVerticesAndNormals (tempVertices, tempNormals, numNeighbours);

	// for (int j=0; j < num_threads; j++)
	// {
	// 	#ifdef ENABLE_BOOST_THREADSoptimizeGrids
	// 	thdGroup.create_thread(boost::bind(&Surface::getTempVerticesAndNormals, this, tempVertices, tempNormals, numNeighbours, j, num_threads));
	// 	#else
	// 	getTempVerticesAndNormals (tempVertices, tempNormals, numNeighbours, j, num_threads);
	// 	#endif
	// }
	// #ifdef ENABLE_BOOST_THREADS
	// thdGroup.join_all();
	// #endif

	for (int j=0; j < num_threads; j++)
	{
		#ifdef ENABLE_BOOST_THREADSoptimizeGrids
		thdGroup.create_thread(boost::bind(&Surface::getNewVerticesAndNormals, this, tempVertices, tempNormals, numNeighbours, j, num_threads));
		#else
		getNewVerticesAndNormals (tempVertices, tempNormals, numNeighbours, j, num_threads);
		#endif
	}
	#ifdef ENABLE_BOOST_THREADS
	thdGroup.join_all();
	#endif

	// delete all
	deleteVector<VERTEX_TYPE>(tempVertices);

	if (computeNormals)
		deleteVector<VERTEX_TYPE>(tempNormals);

	deleteVector<int>(numNeighbours);

	auto chrono_end = chrono::high_resolution_clock::now();
	chrono::duration<double> smoothing_time = chrono_end - chrono_start;
	cout << endl << INFO << "Smoothing time ";
	printf ("%.4e [s]", smoothing_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	double current_mem_in_MB, peak_mem_in_MB;
	getMemSpace (current_mem_in_MB, peak_mem_in_MB);
	cout << endl << INFO << "Memory required after smoothing is " << current_mem_in_MB << " MB";
	#endif

	// check atoms flag
	if (buildAtomsMapHere && vertexAtomsMapFlag)
	{
		auto chrono_start = chrono::high_resolution_clock::now();

		if (vertexAtomsMap != NULL)
			deleteVector<int>(vertexAtomsMap);

		vertexAtomsMap = allocateVector<int>(nv);
		buildAtomsMap();

		cout << endl << INFO << "Connecting vertices to atoms...";
		cout.flush();

		for (int i=0; i<nv; i++)
		{
			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			#if !defined(FLOAT_VERTICES)
			vdwAccessible(vertList[i], vertexAtomsMap[i]);
			#else
			double v[3];
			v[0] = vertList[i][0];
			v[1] = vertList[i][1];
			v[2] = vertList[i][2];
			vdwAccessible(v, vertexAtomsMap[i]);
			#endif
			#else
			#if !defined(FLOAT_VERTICES)
			vdwAccessible(&vertList[i*3], vertexAtomsMap[i]);
			#else
			double v[3];
			v[0] = vertList[ i*3+0 ];
			v[1] = vertList[ i*3+1 ];
			v[2] = vertList[ i*3+2 ];
			vdwAccessible(v, vertexAtomsMap[i]);
			#endif
			#endif

			if (vertexAtomsMap[i] == -1)
			{
				cout << endl << WARN << "Cannot detect nearest atom for vertex " << i;
			}
		}
		cout << "ok!";
		cout.flush();

		auto chrono_end = chrono::high_resolution_clock::now();
		chrono::duration<double> other_time = chrono_end - chrono_start;
		cout << endl << INFO << "Vertex atom mapping time is ";
		printf ("%.4e [s]", other_time.count());

		#if !defined(AVOID_MEM_CHECKS)
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after vertex atom mapping is " << current_mem_in_MB << " MB";
		#endif
	}

	#if !defined(AVOID_SAVING_MESH)
	if (outputMesh)
	{
		auto chrono_start = chrono::high_resolution_clock::now();
		int format = deduceFormat();
		bool f = saveMesh(format, revert, fileName, vertList, triList, normalsList);

		if (!f)
		{
			cout << endl << ERR << "Errors in saving the mesh!";
		}
		else
		{
			auto chrono_end = chrono::high_resolution_clock::now();
			chrono::duration<double> saveMesh_time = chrono_end - chrono_start;
			cout << endl << INFO << "Outputting mesh time (in smoothSurface()) ";
			printf ("%.4e [s]", saveMesh_time.count());
		}

		// Code usable for outputting mesh data in binary format

		// chrono_start = chrono::high_resolution_clock::now();
		// f = saveMeshBinary(format, revert, fileName, vertList, triList, normalsList);

		// if (!f)
		// {
		// 	cout << endl << ERR << "Errors in saving the mesh in binary format!";
		// }
		// else
		// {
		// 	auto chrono_end = chrono::high_resolution_clock::now();
		// 	chrono::duration<double> saveMesh_time = chrono_end - chrono_start;
		// 	cout << endl << INFO << "Outputting mesh time in binary format (in smoothSurface()) ";
		// 	printf ("%.4e [s]", saveMesh_time.count());
		// }
	}
	#endif
}
*/


// Try to use the above version if noted to be profitable

void Surface::smoothSurface(bool outputMesh, bool buildAtomsMapHere,
							const char *fileName, bool revert)
{
	auto chrono_start = chrono::high_resolution_clock::now();

	VERTEX_TYPE *tempVertices;
	VERTEX_TYPE *tempNormals;

	int *num_neighbours;

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int nv = (int)vertList.size();
	#else
	int nv = (int)(vertList.size() / 3.);
	#endif
	int nt = (int)(triList.size() / 3.);

	tempVertices = allocateVector<VERTEX_TYPE>(nv*3);
	memset(tempVertices, 0, nv*3*sizeof(VERTEX_TYPE));

	if (computeNormals)
	{
		tempNormals = allocateVector<VERTEX_TYPE>(nv*3);
		memset(tempNormals, 0, nv*3*sizeof(VERTEX_TYPE));
	}
	num_neighbours = allocateVector<int>(nv);
	memset(num_neighbours, 0, nv*sizeof(int));

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	for (int i=0; i<nt; i++)
	{
		int v1 = triList[ i*3+0 ];
		int v2 = triList[ i*3+1 ];
		int v3 = triList[ i*3+2 ];

		tempVertices[ v1*3+0 ] += vertList[v2][0] + vertList[v3][0];
		tempVertices[ v1*3+1 ] += vertList[v2][1] + vertList[v3][1];
		tempVertices[ v1*3+2 ] += vertList[v2][2] + vertList[v3][2];
		tempVertices[ v2*3+0 ] += vertList[v1][0] + vertList[v3][0];
		tempVertices[ v2*3+1 ] += vertList[v1][1] + vertList[v3][1];
		tempVertices[ v2*3+2 ] += vertList[v1][2] + vertList[v3][2];
		tempVertices[ v3*3+0 ] += vertList[v1][0] + vertList[v2][0];
		tempVertices[ v3*3+1 ] += vertList[v1][1] + vertList[v2][1];
		tempVertices[ v3*3+2 ] += vertList[v1][2] + vertList[v2][2];

		if (computeNormals)
		{
			tempNormals[ v1*3+0 ] += normalsList[v2][0] + normalsList[v3][0];
			tempNormals[ v1*3+1 ] += normalsList[v2][1] + normalsList[v3][1];
			tempNormals[ v1*3+2 ] += normalsList[v2][2] + normalsList[v3][2];
			tempNormals[ v2*3+0 ] += normalsList[v1][0] + normalsList[v3][0];
			tempNormals[ v2*3+1 ] += normalsList[v1][1] + normalsList[v3][1];
			tempNormals[ v2*3+2 ] += normalsList[v1][2] + normalsList[v3][2];
			tempNormals[ v3*3+0 ] += normalsList[v1][0] + normalsList[v2][0];
			tempNormals[ v3*3+1 ] += normalsList[v1][1] + normalsList[v2][1];
			tempNormals[ v3*3+2 ] += normalsList[v1][2] + normalsList[v2][2];
		}
		// Note: every edge is shared by 2 triangles which then contribute twice to each
		// vertex, but also the vertex and normal data sums are overestimated by 2X
		num_neighbours[ v1 ] += 2;
		num_neighbours[ v2 ] += 2;
		num_neighbours[ v3 ] += 2;
	}
	for (int i=0; i<nv; i++)
	{
		VERTEX_TYPE scale = 1.0 / (VERTEX_TYPE)num_neighbours[i];

		vertList[i][0] = 0.5*(vertList[i][0] + tempVertices[ i*3+0 ]*scale);
		vertList[i][1] = 0.5*(vertList[i][1] + tempVertices[ i*3+1 ]*scale);
		vertList[i][2] = 0.5*(vertList[i][2] + tempVertices[ i*3+2 ]*scale);

		if (computeNormals)
		{
			// a 0.5X front multiplication is not needed because of the below normalisation
			tempNormals[ i*3+0 ] = normalsList[i][0] + tempNormals[ i*3+0 ]*scale;
			tempNormals[ i*3+1 ] = normalsList[i][1] + tempNormals[ i*3+1 ]*scale;
			tempNormals[ i*3+2 ] = normalsList[i][2] + tempNormals[ i*3+2 ]*scale;

			VERTEX_TYPE tt = 1.0 / sqrt(
				tempNormals[i*3+0] * tempNormals[i*3+0] +
				tempNormals[i*3+1] * tempNormals[i*3+1] +
				tempNormals[i*3+2] * tempNormals[i*3+2]);

			normalsList[i][0] = tempNormals[ i*3+0 ] * tt;
			normalsList[i][1] = tempNormals[ i*3+1 ] * tt;
			normalsList[i][2] = tempNormals[ i*3+2 ] * tt;
		}
	}
	#else // USE_OPTIMIZED_VERTICES_BUFFERING
	for (int i=0; i<nt; i++)
	{
		int v1 = triList[ i*3+0 ];
		int v2 = triList[ i*3+1 ];
		int v3 = triList[ i*3+2 ];

		tempVertices[ v1*3+0 ] += vertList[ v2*3+0 ] + vertList[ v3*3+0 ];
		tempVertices[ v1*3+1 ] += vertList[ v2*3+1 ] + vertList[ v3*3+1 ];
		tempVertices[ v1*3+2 ] += vertList[ v2*3+2 ] + vertList[ v3*3+2 ];
		tempVertices[ v2*3+0 ] += vertList[ v1*3+0 ] + vertList[ v3*3+0 ];
		tempVertices[ v2*3+1 ] += vertList[ v1*3+1 ] + vertList[ v3*3+1 ];
		tempVertices[ v2*3+2 ] += vertList[ v1*3+2 ] + vertList[ v3*3+2 ];
		tempVertices[ v3*3+0 ] += vertList[ v1*3+0 ] + vertList[ v2*3+0 ];
		tempVertices[ v3*3+1 ] += vertList[ v1*3+1 ] + vertList[ v2*3+1 ];
		tempVertices[ v3*3+2 ] += vertList[ v1*3+2 ] + vertList[ v2*3+2 ];

		if (computeNormals)
		{
			tempNormals[ v1*3+0 ] += normalsList[ v2*3+0 ] + normalsList[ v3*3+0 ];
			tempNormals[ v1*3+1 ] += normalsList[ v2*3+1 ] + normalsList[ v3*3+1 ];
			tempNormals[ v1*3+2 ] += normalsList[ v2*3+2 ] + normalsList[ v3*3+2 ];
			tempNormals[ v2*3+0 ] += normalsList[ v1*3+0 ] + normalsList[ v3*3+0 ];
			tempNormals[ v2*3+1 ] += normalsList[ v1*3+1 ] + normalsList[ v3*3+1 ];
			tempNormals[ v2*3+2 ] += normalsList[ v1*3+2 ] + normalsList[ v3*3+2 ];
			tempNormals[ v3*3+0 ] += normalsList[ v1*3+0 ] + normalsList[ v2*3+0 ];
			tempNormals[ v3*3+1 ] += normalsList[ v1*3+1 ] + normalsList[ v2*3+1 ];
			tempNormals[ v3*3+2 ] += normalsList[ v1*3+2 ] + normalsList[ v2*3+2 ];
		}
		// Note: every edge is shared by 2 triangles which then contribute twice to each
		// vertex, but also the vertex and normal data sums are overestimated by 2X
		num_neighbours[ v1 ] += 2;
		num_neighbours[ v2 ] += 2;
		num_neighbours[ v3 ] += 2;
	}
	for (int i=0; i<nv; i++)
	{
		VERTEX_TYPE scale = 1.0 / (VERTEX_TYPE)num_neighbours[i];

		vertList[ i*3+0 ] = 0.5*(vertList[ i*3+0 ] + tempVertices[ i*3+0 ]*scale);
		vertList[ i*3+1 ] = 0.5*(vertList[ i*3+1 ] + tempVertices[ i*3+1 ]*scale);
		vertList[ i*3+2 ] = 0.5*(vertList[ i*3+2 ] + tempVertices[ i*3+2 ]*scale);

		if (computeNormals)
		{
			// a 0.5X front multiplication is not needed because of the below normalisation
			tempNormals[ i*3+0 ] = normalsList[ i*3+0 ] + tempNormals[ i*3+0 ]*scale;
			tempNormals[ i*3+1 ] = normalsList[ i*3+1 ] + tempNormals[ i*3+1 ]*scale;
			tempNormals[ i*3+2 ] = normalsList[ i*3+2 ] + tempNormals[ i*3+2 ]*scale;

			VERTEX_TYPE tt = 1.0 / sqrt(
				tempNormals[i*3+0] * tempNormals[i*3+0] +
				tempNormals[i*3+1] * tempNormals[i*3+1] +
				tempNormals[i*3+2] * tempNormals[i*3+2]);

			normalsList[ i*3+0 ] = tempNormals[ i*3+0 ] * tt;
			normalsList[ i*3+1 ] = tempNormals[ i*3+1 ] * tt;
			normalsList[ i*3+2 ] = tempNormals[ i*3+2 ] * tt;
		}
	}
	#endif // USE_OPTIMIZED_VERTICES_BUFFERING

	// delete all
	deleteVector<VERTEX_TYPE>(tempVertices);
	if (computeNormals)
		deleteVector<VERTEX_TYPE>(tempNormals);

	deleteVector<int>(num_neighbours);

	auto chrono_end = chrono::high_resolution_clock::now();
	chrono::duration<double> smoothing_time = chrono_end - chrono_start;
	cout << endl << INFO << "Smoothing time ";
	printf ("%.4e [s]", smoothing_time.count());

	#if !defined(AVOID_MEM_CHECKS)
	double current_mem_in_MB, peak_mem_in_MB;
	getMemSpace (current_mem_in_MB, peak_mem_in_MB);
	cout << endl << INFO << "Memory required after smoothing is " << current_mem_in_MB << " MB";
	#endif

	// check atoms flag
	if (buildAtomsMapHere && vertexAtomsMapFlag)
	{
		auto chrono_start = chrono::high_resolution_clock::now();

		if (vertexAtomsMap != NULL)
			deleteVector<int>(vertexAtomsMap);

		vertexAtomsMap = allocateVector<int>(nv);
		buildAtomsMap();

		cout << endl << INFO << "Connecting vertices to atoms...";
		cout.flush();

		for (int i=0; i<nv; i++)
		{
				#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
				#if !defined(FLOAT_VERTICES)
				vdwAccessible(vertList[i], vertexAtomsMap[i]);
				#else
				double v[3];
				v[0] = vertList[i][0];
				v[1] = vertList[i][1];
				v[2] = vertList[i][2];
				vdwAccessible(v, vertexAtomsMap[i]);
				#endif
				#else
				#if !defined(FLOAT_VERTICES)
				vdwAccessible(&vertList[i*3], vertexAtomsMap[i]);
				#else
				double v[3];
				v[0] = vertList[ i*3+0 ];
				v[1] = vertList[ i*3+1 ];
				v[2] = vertList[ i*3+2 ];
				vdwAccessible(v, vertexAtomsMap[i]);
				#endif
				#endif

			if (vertexAtomsMap[i] == -1)
			{
				cout << endl << WARN << "Cannot detect nearest atom for vertex " << i;
			}
		}
		cout << "ok!";
		cout.flush();

		auto chrono_end = chrono::high_resolution_clock::now();
		chrono::duration<double> other_time = chrono_end - chrono_start;
		cout << endl << INFO << "Vertex atom mapping time is ";
		printf ("%.4e [s]", other_time.count());

		#if !defined(AVOID_MEM_CHECKS)
		double current_mem_in_MB, peak_mem_in_MB;
		getMemSpace (current_mem_in_MB, peak_mem_in_MB);
		cout << endl << INFO << "Memory required after vertex atom mapping is " << current_mem_in_MB << " MB";
		#endif
	}

	#if !defined(AVOID_SAVING_MESH)
	if (outputMesh)
	{
		auto chrono_start = chrono::high_resolution_clock::now();
		int format = deduceFormat();
		bool f = saveMesh(format, revert, fileName, vertList, triList, normalsList);

		if (!f)
		{
			cout << endl << ERR << "Errors in saving the mesh!";
		}
		else
		{
			auto chrono_end = chrono::high_resolution_clock::now();
			chrono::duration<double> saveMesh_time = chrono_end - chrono_start;
			cout << endl << INFO << "Outputting mesh time (in smoothSurface()) ";
			printf ("%.4e [s]", saveMesh_time.count());
		}

		// Code usable for outputting mesh data in binary format

		// chrono_start = chrono::high_resolution_clock::now();
		// f = saveMeshBinary(format, revert, fileName, vertList, triList, normalsList);

		// if (!f)
		// {
		// 	cout << endl << ERR << "Errors in saving the mesh in binary format!";
		// }
		// else
		// {
		// 	auto chrono_end = chrono::high_resolution_clock::now();
		// 	chrono::duration<double> saveMesh_time = chrono_end - chrono_start;
		// 	cout << endl << INFO << "Outputting mesh time in binary format (in smoothSurface()) ";
		// 	printf ("%.4e [s]", saveMesh_time.count());
		// }
	}
	#endif
}


int Surface::deduceFormat()
{
	int format = -1;

	if (savePLY)
	{
		format = PLY;
	}
	else if (!saveMSMS)
	{
		if (vertexAtomsMapFlag && computeNormals)
			format = OFF_N_A;
		else
		{
			if (vertexAtomsMapFlag)
				format = OFF_A;
			else if (computeNormals)
				format = OFF_N;
			else
				format = OFF;
		}
	}
	else
	{
		if (vertexAtomsMapFlag)
			format = MSMS;
		else
			format = MSMS_NO_A;
	}
	return format;
}


char Surface::getInsidness(int i, int j, int k, int vertInd)
{
	int insideness[8];

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	int *status = delphi->status;
	int **bilevel_status = delphi->bilevel_status;


	if (!optimizeGrids && vertInd == 3)
	{
		//insideness[0] = delphi->status[k][j][i];
		//insideness[0] = STATUSMAP(i,j,k,NX,NY);
		insideness[0] = read3DVector<int>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i-1>=0)
			//insideness[1] = delphi->status[k][j][i-1];
			//insideness[1] = STATUSMAP((i-1),j,k,NX,NY);
			insideness[1] = read3DVector<int>(status,i-1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i-1>=0 && k-1>=0)
			//insideness[2] = delphi->status[k-1][j][i-1];
			//insideness[2] = STATUSMAP((i-1),j,(k-1),NX,NY);
			insideness[2] = read3DVector<int>(status,i-1,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (k-1>=0)
			//insideness[3] = delphi->status[k-1][j][i];
			//insideness[3] = STATUSMAP(i,j,(k-1),NX,NY);
			insideness[3] = read3DVector<int>(status,i,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (j+1<NY)
			//insideness[4] = delphi->status[k][j+1][i];
			//insideness[4] = STATUSMAP(i,(j+1),k,NX,NY);
			insideness[4] = read3DVector<int>(status,i,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (j+1<NY && i-1>=0)
			//insideness[5] = delphi->status[k][j+1][i-1];
			//insideness[5] = STATUSMAP((i-1),(j+1),k,NX,NY);
			insideness[5] = read3DVector<int>(status,i-1,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (j+1<NY && i-1>=0 && k-1>=0)
			//insideness[6] = delphi->status[k-1][j+1][i-1];
			//insideness[6] = STATUSMAP((i-1),(j+1),(k-1),NX,NY);
			insideness[6] = read3DVector<int>(status,i-1,j+1,k-1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (j+1<NY && k-1>=0)
			//insideness[7] = delphi->status[k-1][j+1][i];
			//insideness[7] = STATUSMAP(i,(j+1),(k-1),NX,NY);
			insideness[7] = read3DVector<int>(status,i,j+1,k-1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (!optimizeGrids && vertInd == 2)
	{
		//insideness[0] = STATUSMAP(i,j,k,NX,NY);
		insideness[0] = read3DVector<int>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i+1<NX)
			//insideness[1] = STATUSMAP((i+1),(j),(k),NX,NY);
			insideness[1] = read3DVector<int>(status,i+1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i+1<NX && j+1<NY)
			//insideness[2] = STATUSMAP((i+1),(j+1),(k),NX,NY);
			insideness[2] = read3DVector<int>(status,i+1,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (j+1<NY)
			//insideness[3] = STATUSMAP((i),(j+1),(k),NX,NY);
			insideness[3] = read3DVector<int>(status,i,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (k-1>=0)
			//insideness[4] = STATUSMAP((i),(j),(k-1),NX,NY);
			insideness[4] = read3DVector<int>(status,i,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (k-1>=0 && i+1<NX)
			//insideness[5] = STATUSMAP((i+1),(j),(k-1),NX,NY);
			insideness[5] = read3DVector<int>(status,i+1,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (k-1>=0 && i+1<NX && j+1<NY)
			//insideness[6] = STATUSMAP((i+1),(j+1),(k-1),NX,NY);
			insideness[6] = read3DVector<int>(status,i+1,j+1,k-1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (k-1>=0 && j+1<NY)
			//insideness[7] = STATUSMAP((i),(j+1),(k-1),NX,NY);
			insideness[7] = read3DVector<int>(status,i,j+1,k-1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (!optimizeGrids && vertInd == 1)
	{
		//insideness[0] = STATUSMAP(i,j,k,NX,NY);
		insideness[0] = read3DVector<int>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i+1<NX)
			//insideness[1] = STATUSMAP((i+1),(j),(k),NX,NY);
			insideness[1] = read3DVector<int>(status,i+1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;

		if (insideness[1] > 0)
			insideness[1] = 1;

		// if (i+1<NX && j-1>=0)
			//insideness[2] = STATUSMAP((i+1),(j-1),(k),NX,NY);
			insideness[2] = read3DVector<int>(status,i+1,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (j-1>=0)
			//insideness[3] = STATUSMAP((i),(j-1),(k),NX,NY);
			insideness[3] = read3DVector<int>(status,i,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (k-1>=0)
			//insideness[4] = STATUSMAP((i),(j),(k-1),NX,NY);
			insideness[4] = read3DVector<int>(status,i,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (k-1>=0 && i+1<NX)
			//insideness[5] = STATUSMAP((i+1),(j),(k-1),NX,NY);
			insideness[5] = read3DVector<int>(status,i+1,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (k-1>=0 && i+1<NX && j-1>=0)
			//insideness[6] = STATUSMAP((i+1),(j-1),(k-1),NX,NY);
			insideness[6] = read3DVector<int>(status,i+1,j-1,k-1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (k-1>=0 && j-1>=0)
			//insideness[7] = STATUSMAP((i),(j-1),(k-1),NX,NY);
			insideness[7] = read3DVector<int>(status,i,j-1,k-1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (!optimizeGrids && vertInd == 0)
	{
		//insideness[0] = STATUSMAP(i,j,k,NX,NY);
		insideness[0] = read3DVector<int>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i-1>=0)
			//insideness[1] = STATUSMAP((i-1),j,k,NX,NY);
			insideness[1] = read3DVector<int>(status,i-1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i-1>=0 && k-1>=0)
			//insideness[2] = STATUSMAP((i-1),(j),(k-1),NX,NY);
			insideness[2] = read3DVector<int>(status,i-1,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (k-1>=0)
			//insideness[3] = STATUSMAP((i),(j),(k-1),NX,NY);
			insideness[3] = read3DVector<int>(status,i,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (j-1>=0)
			//insideness[4] = STATUSMAP((i),(j-1),(k),NX,NY);
			insideness[4] = read3DVector<int>(status,i,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (j-1>=0 && i-1>=0)
			//insideness[5] = STATUSMAP((i-1),(j-1),(k),NX,NY);
			insideness[5] = read3DVector<int>(status,i-1,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (j-1>=0 && i-1>=0 && k-1>=0)
			//insideness[6] = STATUSMAP((i-1),(j-1),(k-1),NX,NY);
			insideness[6] = read3DVector<int>(status,i-1,j-1,k-1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (j-1>=0 && k-1>=0)
			//insideness[7] = STATUSMAP((i),(j-1),(k-1),NX,NY);
			insideness[7] = read3DVector<int>(status,i,j-1,k-1,NX,NY,NZ);
		// else
		// 	i	#else // OPTIMIZE_GRIDSnsideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (!optimizeGrids && vertInd == 7)
	{
		//insideness[0] = STATUSMAP(i,j,k,NX,NY);
		insideness[0] = read3DVector<int>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i-1>=0)
			//insideness[1] = STATUSMAP((i-1),j,k,NX,NY);
			insideness[1] = read3DVector<int>(status,i-1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i-1>=0 && k+1<NZ)
			//insideness[2] = STATUSMAP((i-1),j,(k+1),NX,NY);
			insideness[2] = read3DVector<int>(status,i-1,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (k+1<NZ)
			//insideness[3] = STATUSMAP((i),(j),(k+1),NX,NY);
			insideness[3] = read3DVector<int>(status,i,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (j+1<NY)
			//insideness[4] = STATUSMAP((i),(j+1),(k),NX,NY);
			insideness[4] = read3DVector<int>(status,i,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (j+1<NY && i-1>=0)
			//insideness[5] = STATUSMAP((i-1),(j+1),(k),NX,NY);
			insideness[5] = read3DVector<int>(status,i-1,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (j+1<NY && i-1>=0 && k+1<NZ)
			//insideness[6] = STATUSMAP((i-1),(j+1),(k+1),NX,NY);
			insideness[6] = read3DVector<int>(status,i-1,j+1,k+1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (j+1<NY && k+1<NZ)
			//insideness[7] = STATUSMAP((i),(j+1),(k+1),NX,NY);
			insideness[7] = read3DVector<int>(status,i,j+1,k+1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (!optimizeGrids && vertInd == 6)
	{
		//insideness[0] = STATUSMAP(i,j,k,NX,NY);
		insideness[0] = read3DVector<int>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i+1<NX)
			//insideness[1] = STATUSMAP((i+1),(j),(k),NX,NY);
			insideness[1] = read3DVector<int>(status,i+1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i+1<NX && j+1<NY)
			//insideness[2] = STATUSMAP((i+1),(j+1),(k),NX,NY);
			insideness[2] = read3DVector<int>(status,i+1,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (j+1<NY)
			//insideness[3] = STATUSMAP((i),(j+1),(k),NX,NY);
			insideness[3] = read3DVector<int>(status,i,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (k+1<NZ)
			//insideness[4] = STATUSMAP((i),(j),(k+1),NX,NY);
			insideness[4] = read3DVector<int>(status,i,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (k+1<NZ && i+1<NX)
			//insideness[5] = STATUSMAP((i+1),(j),(k+1),NX,NY);
			insideness[5] = read3DVector<int>(status,i+1,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (k+1<NZ && i+1<NX && j+1<NY)
			//insideness[6] = STATUSMAP((i+1),(j+1),(k+1),NX,NY);
			insideness[6] = read3DVector<int>(status,i+1,j+1,k+1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (k+1<NZ && j+1<NY)
			//insideness[7] = STATUSMAP((i),(j+1),(k+1),NX,NY);
			insideness[7] = read3DVector<int>(status,i,j+1,k+1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (!optimizeGrids && vertInd == 5)
	{
		//insideness[0] = STATUSMAP(i,j,k,NX,NY);
		insideness[0] = read3DVector<int>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i+1<NX)
			//insideness[1] = STATUSMAP((i+1),(j),(k),NX,NY);
			insideness[1] = read3DVector<int>(status,i+1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i+1<NX && j-1>=0)
			//insideness[2] = STATUSMAP((i+1),(j-1),(k),NX,NY);
			insideness[2] = read3DVector<int>(status,i+1,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (j-1>=0)
			//insideness[3] = STATUSMAP((i),(j-1),(k),NX,NY);
			insideness[3] = read3DVector<int>(status,i,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (k+1<NZ)
			//insideness[4] = STATUSMAP((i),(j),(k+1),NX,NY);
			insideness[4] = read3DVector<int>(status,i,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (k+1<NZ && i+1<NX)
			//insideness[5] = STATUSMAP((i+1),(j),(k+1),NX,NY);
			insideness[5] = read3DVector<int>(status,i+1,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (k+1<NZ && i+1<NX && j-1>=0)
			//insideness[6] = STATUSMAP((i+1),(j-1),(k+1),NX,NY);
			insideness[6] = read3DVector<int>(status,i+1,j-1,k+1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (k+1<NZ && j-1>=0)
			//insideness[7] = STATUSMAP((i),(j-1),(k+1),NX,NY);
			insideness[7] = read3DVector<int>(status,i,j-1,k+1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (!optimizeGrids && vertInd == 4)
	{
		//insideness[0] = STATUSMAP(i,j,k,NX,NY);
		insideness[0] = read3DVector<int>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i-1>=0)
			//insideness[1] = STATUSMAP((i-1),j,k,NX,NY);
			insideness[1] = read3DVector<int>(status,i-1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i-1>=0 && k+1<NZ)
			//insideness[2] = STATUSMAP((i-1),j,(k+1),NX,NY);
			insideness[2] = read3DVector<int>(status,i-1,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (k+1<NZ)
			//insideness[3] = STATUSMAP((i),(j),(k+1),NX,NY);
			insideness[3] = read3DVector<int>(status,i,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (j-1>=0)
			//insideness[4] = STATUSMAP((i),(j-1),(k),NX,NY);
			insideness[4] = read3DVector<int>(status,i,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (j-1>=0 && i-1>=0)
			//insideness[5] = STATUSMAP((i-1),(j-1),(k),NX,NY);
			insideness[5] = read3DVector<int>(status,i-1,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (j-1>=0 && i-1>=0 && k+1<NZ)
			//insideness[6] = STATUSMAP((i-1),(j-1),(k+1),NX,NY);
			insideness[6] = read3DVector<int>(status,i-1,j-1,k+1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (j-1>=0 && k+1<NZ)
			//insideness[7] = STATUSMAP((i),(j-1),(k+1),NX,NY);
			insideness[7] = read3DVector<int>(status,i,j-1,k+1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (optimizeGrids && vertInd == 3)
	{
		insideness[0] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i-1>=0)
			insideness[1] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i-1>=0 && k-1>=0)
			insideness[2] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (k-1>=0)
			insideness[3] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (j+1<NY)
			insideness[4] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (j+1<NY && i-1>=0)
			insideness[5] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (j+1<NY && i-1>=0 && k-1>=0)
			insideness[6] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j+1,k-1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (j+1<NY && k-1>=0)
			insideness[7] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j+1,k-1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (optimizeGrids && vertInd == 2)
	{
		insideness[0] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i+1<NX)
			insideness[1] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i+1<NX && j+1<NY)
			insideness[2] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (j+1<NY)
			insideness[3] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (k-1>=0)
			insideness[4] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (k-1>=0 && i+1<NX)
			insideness[5] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (k-1>=0 && i+1<NX && j+1<NY)
			insideness[6] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j+1,k-1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (k-1>=0 && j+1<NY)
			insideness[7] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j+1,k-1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (optimizeGrids && vertInd == 1)
	{
		insideness[0] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i+1<NX)
			insideness[1] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i+1<NX && j-1>=0)
			insideness[2] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (j-1>=0)
			insideness[3] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (k-1>=0)
			insideness[4] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (k-1>=0 && i+1<NX)
			insideness[5] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (k-1>=0 && i+1<NX && j-1>=0)
			insideness[6] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j-1,k-1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (k-1>=0 && j-1>=0)
			insideness[7] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j-1,k-1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (optimizeGrids && vertInd == 0)
	{
		insideness[0] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i-1>=0)
			insideness[1] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i-1>=0 && k-1>=0)
			insideness[2] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (k-1>=0)
			insideness[3] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k-1,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (j-1>=0)
			insideness[4] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (j-1>=0 && i-1>=0)
			insideness[5] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (j-1>=0 && i-1>=0 && k-1>=0)
			insideness[6] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j-1,k-1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (j-1>=0 && k-1>=0)
			insideness[7] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j-1,k-1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (optimizeGrids && vertInd == 7)
	{
		insideness[0] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i-1>=0)
			insideness[1] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i-1>=0 && k+1<NZ)
			insideness[2] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (k+1<NZ)
			insideness[3] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (j+1<NY)
			insideness[4] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (j+1<NY && i-1>=0)
			insideness[5] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (j+1<NY && i-1>=0 && k+1<NZ)
			insideness[6] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j+1,k+1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (j+1<NY && k+1<NZ)
			insideness[7] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j+1,k+1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (optimizeGrids && vertInd == 6)
	{
		insideness[0] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i+1<NX)
			insideness[1] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i+1<NX && j+1<NY)
			insideness[2] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (j+1<NY)
			insideness[3] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j+1,k,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (k+1<NZ)
			insideness[4] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (k+1<NZ && i+1<NX)
			insideness[5] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (k+1<NZ && i+1<NX && j+1<NY)
			insideness[6] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j+1,k+1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (k+1<NZ && j+1<NY)
			insideness[7] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j+1,k+1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (optimizeGrids && vertInd == 5)
	{
		insideness[0] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i+1<NX)
			insideness[1] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i+1<NX && j-1>=0)
			insideness[2] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (j-1>=0)
			insideness[3] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (k+1<NZ)
			insideness[4] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (k+1<NZ && i+1<NX)
			insideness[5] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (k+1<NZ && i+1<NX && j-1>=0)
			insideness[6] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i+1,j-1,k+1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (k+1<NZ && j-1>=0)
			insideness[7] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j-1,k+1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else if (optimizeGrids && vertInd == 4)
	{
		insideness[0] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0] > 0)
			insideness[0] = 1;

		// if (i-1>=0)
			insideness[1] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j,k,NX,NY,NZ);
		// else
		// 	insideness[1] = 1;
		
		if (insideness[1] > 0)
			insideness[1] = 1;
	
		// if (i-1>=0 && k+1<NZ)
			insideness[2] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[2] = 1;

		if (insideness[2] > 0)
			insideness[2] = 1;

		// if (k+1<NZ)
			insideness[3] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j,k+1,NX,NY,NZ);
		// else
		// 	insideness[3] = 1;

		if (insideness[3] > 0)
			insideness[3] = 1;

		// if (j-1>=0)
			insideness[4] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[4] = 1;

		if (insideness[4] > 0)
			insideness[4] = 1;

		// if (j-1>=0 && i-1>=0)
			insideness[5] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j-1,k,NX,NY,NZ);
		// else
		// 	insideness[5] = 1;

		if (insideness[5] > 0)
			insideness[5] = 1;

		// if (j-1>=0 && i-1>=0 && k+1<NZ)
			insideness[6] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i-1,j-1,k+1,NX,NY,NZ);
		// else
		// 	insideness[6] = 1;

		if (insideness[6] > 0)
			insideness[6] = 1;

		// if (j-1>=0 && k+1<NZ)
			insideness[7] = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,i,j-1,k+1,NX,NY,NZ);
		// else
		// 	insideness[7] = 1;

		if (insideness[7] > 0)
			insideness[7] = 1;
	}
	else
	{
		cout << endl << ERR << "Unknown vertex index!";
		return false;
	}
	
	int acc = 0;
	// decide insidness
	for (int i=0; i<8; i++)
	{
		acc += (int)insideness[i];
	}
	// positively biased decision
	if (acc > 0)
		return +1;
	else 
		return -1;
}


void Surface::vertexInterp(double isolevel, double *p1, double *p2, double valp1, double valp2, VERTEX_TYPE *p)
{
   double mu;
   mu = (isolevel - valp1) / (valp2 - valp1);
   double temp[3];
   SUB(temp,p2,p1);
   ADD_MUL(p,p1,temp,mu);
}


inline int Surface::classifyCube(double *vertexValues, double isolevel)
{
	int cubeindex = 0;

	if (vertexValues[0] < isolevel) cubeindex |= 1;
	if (vertexValues[1] < isolevel) cubeindex |= 2;
	if (vertexValues[2] < isolevel) cubeindex |= 4;
	if (vertexValues[3] < isolevel) cubeindex |= 8;
	if (vertexValues[4] < isolevel) cubeindex |= 16;
	if (vertexValues[5] < isolevel) cubeindex |= 32;
	if (vertexValues[6] < isolevel) cubeindex |= 64;
	if (vertexValues[7] < isolevel) cubeindex |= 128;

	// Cube is entirely in/out of the surface 
	if (edgeTable[cubeindex] == 0)
		return -1;
	
	return cubeindex;
}

/**
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the matrix triangles
   will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
   of totally below the isolevel. Can act using interpolation or querying
   the analytically computed intersection point. In this last case it is guaranted
   that every point of the mesh belongs to the surface; ix,iy,iz are passed in order to query
   the intersection matrices
*/
// inline int Surface::getTriangles(double *vertexValues, double isolevel, int **triangles, int ix, int iy, int iz,
// 								 int NX, int NY, int NZ, int *xedge, int *yedge, int *zedge,int *xedge_down, int *yedge_down)
inline int Surface::getTriangles(double *vertexValues, double isolevel, int **triangles,
								 int ix, int iy, int iz, int NX, int NY, int NZ)

{
	// Marching cubes vertex indices and edges convention; vertInd is indicated in the cube
	//
	//		   v4_________e4_________v5
	//			/|  				/|
	//		e7 / |                 / |
	//		  /  |             e5 /  |
	//		 /   | e8            /	 | e9
	//	  v7/____|_____e6_______/v6	 |
	//		|	 |              |	 |
	//	 	|  v0|______e0______|____|v1
	//	e11 |	/        		|   /
	//		|  /			e10	|  /
	//		| /	e3				| / e1
	//		|/					|/
	//	  v3/_________e2________/v2
	//            
	int ntriang;
	int cubeindex;
	
	int vertlist_indexes[12] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
		
	//   Determine the index into the edge table which
	//   tells us which vertices are inside of the surface   
	cubeindex = classifyCube(vertexValues, isolevel);
	
	// Cube is entirely in/out of the surface 
	if (cubeindex == -1)
		return 0;
   
	// int vert_offset, i=ix, j=iy;
	int i=ix; // , j=iy;

	// the vertices are already analytically computed, they are only recovered.

	#if !defined(OPTIMIZE_INTERSECTIONS_MANAGEMENT)

	// here octrees are employed for retrieving vertices indices

	if (edgeTable[cubeindex] & 1)
	{
		/*
		vert_offset = XEDGE_DOWN(i,j-1,NX);
		vertlist_indexes[0] = vert_offset;
		if (vert_offset == -1)
		*/
		vertlist_indexes[0] = intersectionsMatrixAlongX->at(ix,iy,iz);
	}

	if (edgeTable[cubeindex] & 2)
	{
		/*
		vert_offset = YEDGE_DOWN(i,j,NX);
		vertlist_indexes[1] = vert_offset;
		if (vert_offset == -1)
		*/
		vertlist_indexes[1] = intersectionsMatrixAlongY->at(ix+1,iy,iz);
	}

	if (edgeTable[cubeindex] & 4)
	{
		/*
		vert_offset = XEDGE_DOWN(i,j,NX);
		vertlist_indexes[2] = vert_offset;
		if (vert_offset == -1)
		*/
		vertlist_indexes[2] = intersectionsMatrixAlongX->at(ix,iy+1,iz);
	}

	if (edgeTable[cubeindex] & 8)
	{
		/*
		vert_offset = YEDGE_DOWN(i-1,j,NX);
		vertlist_indexes[3] = vert_offset;
		if (vert_offset == -1)
		*/
		vertlist_indexes[3] = intersectionsMatrixAlongY->at(ix,iy,iz);
	}
	
	if (edgeTable[cubeindex] & 16)
	{
		/*
		vert_offset = XEDGE(i,j-1,NX);
		vertlist_indexes[4] = vert_offset;
		if (vert_offset == -1)
		*/
		vertlist_indexes[4] = intersectionsMatrixAlongX->at(ix,iy,iz+1);
	}

	if (edgeTable[cubeindex] & 32)
	{
		vertlist_indexes[5] = intersectionsMatrixAlongY->at(ix+1,iy,iz+1);
		//YEDGE(i,j,NX) = vertlist_indexes[5];
	}

	if (edgeTable[cubeindex] & 64)
	{
		vertlist_indexes[6] = intersectionsMatrixAlongX->at(ix,iy+1,iz+1);
		//XEDGE(i,j,NX) = vertlist_indexes[6];
	}

	if (edgeTable[cubeindex] & 128)
	{
		/*
		vert_offset = YEDGE(i-1,j,NX);
		vertlist_indexes[7] = vert_offset;
		if (vert_offset == -1)
		*/
		vertlist_indexes[7] = intersectionsMatrixAlongY->at(ix,iy,iz+1);
	}

	if (edgeTable[cubeindex] & 256)
	{
		/*
		vert_offset = ZEDGE(i-1,j-1,NX);
		vertlist_indexes[8] = vert_offset;
		if (vert_offset == -1)
		*/
		vertlist_indexes[8] = intersectionsMatrixAlongZ->at(ix,iy,iz);
	}

	if (edgeTable[cubeindex] & 512)
	{
		/*
		vert_offset = ZEDGE(i,j-1,NX);
		vertlist_indexes[9] = vert_offset;
		if (vert_offset == -1)
		*/
		vertlist_indexes[9] = intersectionsMatrixAlongZ->at(ix+1,iy,iz);
	}
	
	if (edgeTable[cubeindex] & 1024)
	{
		vertlist_indexes[10] = intersectionsMatrixAlongZ->at(ix+1,iy+1,iz);
		//ZEDGE(i,j,NX)=vertlist_indexes[10];
	}

	if (edgeTable[cubeindex] & 2048)
	{ 
		/*
		vert_offset = ZEDGE(i-1,j,NX);
		vertlist_indexes[11] = vert_offset;
		if (vert_offset == -1)
		*/
		vertlist_indexes[11] = intersectionsMatrixAlongZ->at(ix,iy+1,iz);
	}

	#else // OPTIMIZE_INTERSECTIONS_MANAGEMENT

	// here bilevel grids are employed for retrieving vertices indices

	if (edgeTable[cubeindex] & 1)
		vertlist_indexes[0] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy,iz,NX,NY,NZ);

	if (edgeTable[cubeindex] & 2)
		vertlist_indexes[1] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix+1,iy,iz,NX,NY,NZ);

	if (edgeTable[cubeindex] & 4)
		vertlist_indexes[2] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy+1,iz,NX,NY,NZ);

	if (edgeTable[cubeindex] & 8)
		vertlist_indexes[3] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix,iy,iz,NX,NY,NZ);

	if (edgeTable[cubeindex] & 16)
		vertlist_indexes[4] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy,iz+1,NX,NY,NZ);

	if (edgeTable[cubeindex] & 32)
		vertlist_indexes[5] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix+1,iy,iz+1,NX,NY,NZ);

	if (edgeTable[cubeindex] & 64)
		vertlist_indexes[6] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongX,-1,ix,iy+1,iz+1,NX,NY,NZ);

	if (edgeTable[cubeindex] & 128)
		vertlist_indexes[7] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongY,-1,ix,iy,iz+1,NX,NY,NZ);

	if (edgeTable[cubeindex] & 256)
		vertlist_indexes[8] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix,iy,iz,NX,NY,NZ);

	if (edgeTable[cubeindex] & 512)
		vertlist_indexes[9] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix+1,iy,iz,NX,NY,NZ);

	if (edgeTable[cubeindex] & 1024)
		vertlist_indexes[10] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix+1,iy+1,iz,NX,NY,NZ);

	if (edgeTable[cubeindex] & 2048)
		vertlist_indexes[11] = readBilevelGrid<int>(bilevel_intersectionsMatrixAlongZ,-1,ix,iy+1,iz,NX,NY,NZ);

	#endif // OPTIMIZE_INTERSECTIONS_MANAGEMENT

    // Create the triangles
    ntriang = 0;

    for (i=0; triTable[cubeindex][i] != -1; i += 3)
    {
		if (vertlist_indexes[triTable[cubeindex][i  ]] == -1)
		{
			cout << endl << WARN << "Mesh with hole!";
			continue;
		}
		if (vertlist_indexes[triTable[cubeindex][i+1]] == -1)
		{
			cout << endl << WARN << "Mesh with hole!";
			continue;
		}
		if (vertlist_indexes[triTable[cubeindex][i+2]] == -1)
		{
			cout << endl << WARN << "Mesh with hole!";
			continue;
		}

		// save all the assigned indexes
		triangles[ntriang][0] = vertlist_indexes[triTable[cubeindex][i  ]];
		triangles[ntriang][1] = vertlist_indexes[triTable[cubeindex][i+1]];
		triangles[ntriang][2] = vertlist_indexes[triTable[cubeindex][i+2]];
		ntriang++;
	}
	return ntriang;
}


void Surface::backupStatus()
{
	if (!optimizeGrids)
	{
		if (delphi == NULL)
		{
			cout << endl << WARN << "Cannot backup status if grid is not allocated";
			return;
		}
		int64_t tot = delphi->nx*delphi->ny*delphi->nz;

		if (delphi->tempStatus != NULL)
			deleteVector<int>(delphi->tempStatus);

		delphi->tempStatus = allocateVector<int>(tot);
		
		for (int64_t i=0; i<tot; i++)
			delphi->tempStatus[i] = delphi->status[i];
	}
	else
	{
		if (delphi == NULL)
		{
			cout << endl << WARN << "Cannot backup bilevel_status if grid is not allocated";
			return;
		}
		if (delphi->bilevel_tempStatus != NULL)
		{
			deleteBilevelGridCells<int>(delphi->bilevel_tempStatus, delphi->nx, delphi->ny, delphi->nz);
			deleteVector<int *>(delphi->bilevel_tempStatus);
		}
		// allocating a coarse grid, delphi->bilevel_tempStatus, of nx/4 * ny/4 * nz/4 coarse cells
		delphi->bilevel_tempStatus = copyBilevelGridCells(delphi->bilevel_status, delphi->nx, delphi->ny, delphi->nz);
	}
}


void Surface::removeBackupStatus()
{
	if (!optimizeGrids)
	{
		if (delphi->tempStatus != NULL)
			deleteVector<int>(delphi->tempStatus);
		delphi->tempStatus = NULL;
	}
	else
	{
		if (delphi->bilevel_tempStatus != NULL)
		{
			deleteBilevelGridCells<int>(delphi->bilevel_tempStatus, delphi->nx, delphi->ny, delphi->nz);
			deleteVector<int *>(delphi->bilevel_tempStatus);
		}
	}
}


int Surface::linkCavities(int *st1, int *st2)
{
	// this is the status of the slim surface
	int *tempStatus = st1;
	int *status = delphi->status;

	if (tempStatus == NULL)
	{
		cout << endl << WARN << "I cannot reconstruct links if I have not reference status";
		return 0;
	}

	// analyze each cavity before and after fill
	vector<vector<int*>*>::iterator it;

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	// cavity under link analysis
	int cavityId = STATUS_FIRST_CAV;

	vector <double>::iterator sizeIt;
	vector <bool>::iterator flagIt;

	// buildAtomsMap();

	for (it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++, cavityId++)
	{
		if (delphi->cavitiesFlag[ cavityId - STATUS_FIRST_CAV ] == true)
			continue;

		vector<int*> *vec1 = (*it);

		if (vec1->size() == 0)
			continue;

		// paired cavity under analysis
		int checkCavId = cavityId+1;

		// check if, for the all the other current cavities, exist a point
		// which had the same id before filtering. If yes a link is estabilished.
		for (checkCavId = cavityId+1;  (checkCavId-STATUS_FIRST_CAV) < (int)delphi->cavitiesVec->size(); checkCavId++)
		{
			if (delphi->cavitiesFlag[ checkCavId - STATUS_FIRST_CAV ] == true)
				continue;

			vector<int*> *vec2 = delphi->cavitiesVec->at( checkCavId - STATUS_FIRST_CAV );

			// this may happen if this cavity have been already aggregated
			if (vec2->size() == 0)
				continue;

			// check the distance between the two nearest bgps of the two cavities.
			// a bgp in this case is defined as a grid point in which at least
			// one grid point is out

			// these are respectively the two current grid points which are compared
			int *v1, *v2;

			vector<int*>::iterator it_c1;
			vector<int*>::iterator it_c2;
			double minDist2 = INFINITY;
			unsigned int winner1, winner2;

			for (unsigned int c1 = 0; c1<vec1->size(); c1++)
			{
				v1 = vec1->at(c1);
				// int ref = STATUSMAP(v1[0],v1[1],v1[2],NX,NY);
				int ref = read3DVector<int>(status,v1[0],v1[1],v1[2],NX,NY,NZ);

				// if is bgp get distance
				/*
				if ((STATUSMAP(v1[0]+1,v1[1],v1[2],NX,NY) != ref) ||
					(STATUSMAP(v1[0]-1,v1[1],v1[2],NX,NY) != ref) ||
					(STATUSMAP(v1[0],v1[1]+1,v1[2],NX,NY) != ref) ||
					(STATUSMAP(v1[0],v1[1]-1,v1[2],NX,NY) != ref) ||
					(STATUSMAP(v1[0],v1[1],v1[2]+1,NX,NY) != ref) ||
					(STATUSMAP(v1[0],v1[1],v1[2]-1,NX,NY) != ref))
				*/
				if ((read3DVector<int>(status,v1[0]+1,v1[1],v1[2],NX,NY,NZ) != ref) ||
					(read3DVector<int>(status,v1[0]-1,v1[1],v1[2],NX,NY,NZ) != ref) ||
					(read3DVector<int>(status,v1[0],v1[1]+1,v1[2],NX,NY,NZ) != ref) ||
					(read3DVector<int>(status,v1[0],v1[1]-1,v1[2],NX,NY,NZ) != ref) ||
					(read3DVector<int>(status,v1[0],v1[1],v1[2]+1,NX,NY,NZ) != ref) ||
					(read3DVector<int>(status,v1[0],v1[1],v1[2]-1,NX,NY,NZ) != ref))
				{
					// do check
				}
				else
					continue;

				for (unsigned int c2 = 0; c2<vec2->size(); c2++)
				{
					v2 = vec2->at(c2);
					// ref = STATUSMAP(v2[0],v2[1],v2[2],NX,NY);
					ref = read3DVector<int>(status,v2[0],v2[1],v2[2],NX,NY,NZ);

					// if is bgp get distance
					/*
					if ((STATUSMAP(v2[0]+1,v2[1],v2[2],NX,NY) != ref) ||
						(STATUSMAP(v2[0]-1,v2[1],v2[2],NX,NY) != ref) ||
						(STATUSMAP(v2[0],v2[1]+1,v2[2],NX,NY) != ref) ||
						(STATUSMAP(v2[0],v2[1]-1,v2[2],NX,NY) != ref) ||
						(STATUSMAP(v2[0],v2[1],v2[2]+1,NX,NY) != ref) ||
						(STATUSMAP(v2[0],v2[1],v2[2]-1,NX,NY) != ref))
					*/
					if ((read3DVector<int>(status,v2[0]+1,v2[1],v2[2],NX,NY,NZ) != ref) ||
						(read3DVector<int>(status,v2[0]-1,v2[1],v2[2],NX,NY,NZ) != ref) ||
						(read3DVector<int>(status,v2[0],v2[1]+1,v2[2],NX,NY,NZ) != ref) ||
						(read3DVector<int>(status,v2[0],v2[1]-1,v2[2],NX,NY,NZ) != ref) ||
						(read3DVector<int>(status,v2[0],v2[1],v2[2]+1,NX,NY,NZ) != ref) ||
						(read3DVector<int>(status,v2[0],v2[1],v2[2]-1,NX,NY,NZ) != ref))
					{
						// get distance
						double a = (delphi->x[v1[0]]-delphi->x[v2[0]])*(delphi->x[v1[0]]-delphi->x[v2[0]]);
						double b = (delphi->y[v1[1]]-delphi->y[v2[1]])*(delphi->y[v1[1]]-delphi->y[v2[1]]);
						double c = (delphi->z[v1[2]]-delphi->z[v2[2]])*(delphi->z[v1[2]]-delphi->z[v2[2]]);
						double dist2 = a+b+c;

						if (dist2 < minDist2)
						{
							winner1 = c1;
							winner2 = c2;
							minDist2 = dist2;
						}
					}
					else
						continue;
				}
			}

			if (minDist2 == INFINITY)
			{
				cout << endl << ERR << "During linkage, two non null cavities/pockets have no minimum distance bgps";
				exit(-1);
			}

			bool merge = false;
			// printf("\n Min distance between %d %d is %f", cavityId-4, checkCavId-4, sqrt(minDist2));

			// bound on the max path
			double refDist = 6;

			if (sqrt(minDist2) > refDist)
				merge = false;

			// if we are under a threshold distance a full check is needed
			// check if one can move freely inside the small probe surface between the two nearest bgps
			// while at the mean time never going out the 1.4 surface.
			else
			{
				int maxMoves = (int)(refDist*delphi->scale + 0.5);
				int moves = maxMoves;
				v1 = vec1->at(winner1);
				v2 = vec2->at(winner2);

				bool debug = false;
				//if (cavityId-4 == 6 && checkCavId-4 == 17)
				//	debug = true;

				// printf("\n %d %d %d -> %d %d %d",v1[0],v1[1],v1[2],v2[0],v2[1],v2[2]);
				merge = innerFloodFill(v1,v2,st1,st2,moves,debug);
				double geodetic = moves*delphi->side;

				if (!merge)
				{
				//	printf("\n Non merging, num moves %d",moves);
				}
				if (merge && geodetic > 1.5*sqrt(minDist2))
				{
				// printf("\n Geodetic is %f and max dist is %f rejecting link",geodetic,1.5*sqrt(minDist2));
					merge = false;
				}

				//if (merge)
				//	printf("\n Merged in %d moves with geodetic %f",moves,geodetic);
			}

			// these two new cavities need to be linked
			// because they were linked by small probe surface accessibility
			if (merge)
			{
				// merge all the grid points
				vector <int*>::iterator it_vec2;

				for (it_vec2 = vec2->begin(); it_vec2 != vec2->end(); it_vec2++)
				{
					int *v = *it_vec2;
					vec1->push_back(v);
				}

				// pointers are not deleted
				vec2->clear();

				// sum volumes of the aggregated cavities
				delphi->cavitiesSize[ cavityId - STATUS_FIRST_CAV ] += delphi->cavitiesSize[ checkCavId - STATUS_FIRST_CAV ];
				// nullify the volume to flag it
				delphi->cavitiesSize[ checkCavId - STATUS_FIRST_CAV ] = 0;
				// assure it is active the merged cavity
				delphi->cavitiesFlag[ cavityId - STATUS_FIRST_CAV ] = false;
			}

		}
	}

	// disposeAtomsMap();

	sizeIt = delphi->cavitiesSize.begin();
	flagIt = delphi->cavitiesFlag.begin();
	it = delphi->cavitiesVec->begin();

	int numRemoved = 0;

	// remove all the null cavities
	while (1)
	{
		if (it == delphi->cavitiesVec->end())
			break;

		vector<int*> *vec1 = (*it);

		// the merged cavities are removed
		if (vec1->size() == 0)
		{
			numRemoved++;
			// delete vector of pointers to points
			delete vec1;
			vec1 = NULL;

			// delete entry
			delphi->cavitiesVec->erase(it);

			// double check it is the right guy
			if (*sizeIt != 0)
			{
				cout << endl << ERR << "Inconsistent volume aggregation";
				exit(-1);
			}

			// delete size (already aggregated)
			delphi->cavitiesSize.erase(sizeIt);
			// delete flag
			delphi->cavitiesFlag.erase(flagIt);

			// restart iterators
			it = delphi->cavitiesVec->begin();
			sizeIt = delphi->cavitiesSize.begin();
			flagIt = delphi->cavitiesFlag.begin();
		}
		else
		{
			it++;
			sizeIt++;
			flagIt++;
		}
	}

	int new_id = STATUS_FIRST_CAV;

	// renumber cavities accordingly
	for (it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++)
	{
		vector<int*>::iterator inner;
		vector<int*> *cav = (*it);

		for (inner = cav->begin(); inner != cav->end(); inner++)
		{
			int *v = *inner;
			// conserve the support cavity coding

			/*
			if (STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) > 0)
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = new_id;
			else
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = -new_id;
			*/
			if (read3DVector<int>(status,v[0],v[1],v[2],NX,NY,NZ) > 0)
				write3DVector<int>(status,new_id,v[0],v[1],v[2],NX,NY,NZ);
			else
				write3DVector<int>(status,-new_id,v[0],v[1],v[2],NX,NY,NZ);
		}
		new_id++;
	}
	return numRemoved;
}


// This version is identical to the above one except for the use of bilevel grids for memory space reduction
int Surface::linkCavities(int **st1, int **st2)
{
	int **bilevel_tempStatus = st1;
	int **bilevel_status = delphi->bilevel_status;
	
	if (bilevel_tempStatus == NULL)
	{
		cout << endl << WARN << "I cannot reconstruct links if I have not reference bilevel_status";
		return 0;
	}

	// analyze each cavity before and after fill
	vector<vector<int*>*>::iterator it;
	
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
		
	// cavity under link analysis
	int cavityId = STATUS_FIRST_CAV;

	vector <double>::iterator sizeIt;
	vector <bool>::iterator flagIt;

	// buildAtomsMap();

	for (it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++, cavityId++)
	{
		if (delphi->cavitiesFlag[ cavityId - STATUS_FIRST_CAV ] == true)
			continue;

		vector<int*> *vec1 = (*it);

		if (vec1->size() == 0)
			continue;

		// paired cavity under analysis
		int checkCavId = cavityId+1;

		// check if, for the all the other current cavities, exist a point
		// which had the same id before filtering. If yes a link is estabilished.
		for (checkCavId = cavityId+1; (checkCavId-STATUS_FIRST_CAV) < (int)delphi->cavitiesVec->size(); checkCavId++)
		{				
			if (delphi->cavitiesFlag[ checkCavId - STATUS_FIRST_CAV ] == true)
				continue;
				
			vector<int*> *vec2 = delphi->cavitiesVec->at( checkCavId - STATUS_FIRST_CAV );

			// this may happen if this cavity have been already aggregated
			if (vec2->size() == 0)
				continue;
				
			// check the distance between the two nearest bgps of the two cavities.
			// a bgp in this case is defined as a grid point in which at least 
			// one grid point is out

			// these are respectively the two current grid points which are compared
			int *v1, *v2;

			vector<int*>::iterator it_c1;
			vector<int*>::iterator it_c2;
			double minDist2 = INFINITY;
			unsigned int winner1, winner2;

			for (unsigned int c1 = 0; c1<vec1->size(); c1++)
			{
				v1 = vec1->at(c1);
				
				int ref = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v1[0],v1[1],v1[2],NX,NY,NZ);
				
				if ((readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v1[0]+1,v1[1],v1[2],NX,NY,NZ) != ref) ||
					(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v1[0]-1,v1[1],v1[2],NX,NY,NZ) != ref) ||
					(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v1[0],v1[1]+1,v1[2],NX,NY,NZ) != ref) ||
					(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v1[0],v1[1]-1,v1[2],NX,NY,NZ) != ref) ||
					(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v1[0],v1[1],v1[2]+1,NX,NY,NZ) != ref) ||
					(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v1[0],v1[1],v1[2]-1,NX,NY,NZ) != ref))
				{
					// do check
				}
				else
					continue;

				for (unsigned int c2 = 0; c2<vec2->size(); c2++)
				{
					v2 = vec2->at(c2);
					ref = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v2[0],v2[1],v2[2],NX,NY,NZ);
					
					if ((readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v2[0]+1,v2[1],v2[2],NX,NY,NZ) != ref) ||
						(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v2[0]-1,v2[1],v2[2],NX,NY,NZ) != ref) ||
						(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v2[0],v2[1]+1,v2[2],NX,NY,NZ) != ref) ||
						(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v2[0],v2[1]-1,v2[2],NX,NY,NZ) != ref) ||
						(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v2[0],v2[1],v2[2]+1,NX,NY,NZ) != ref) ||
						(readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v2[0],v2[1],v2[2]-1,NX,NY,NZ) != ref))
					{
						// get distance
						double a = (delphi->x[v1[0]]-delphi->x[v2[0]])*(delphi->x[v1[0]]-delphi->x[v2[0]]);
						double b = (delphi->y[v1[1]]-delphi->y[v2[1]])*(delphi->y[v1[1]]-delphi->y[v2[1]]);
						double c = (delphi->z[v1[2]]-delphi->z[v2[2]])*(delphi->z[v1[2]]-delphi->z[v2[2]]);
						double dist2 = a+b+c;

						if (dist2 < minDist2)
						{
							winner1 = c1;
							winner2 = c2;
							minDist2 = dist2;
						}
					}
					else
						continue;
				}
			}

			if (minDist2 == INFINITY)
			{
				cout << endl << ERR << "During linkage, two non null cavities/pockets have no minimum distance bgps";
				exit(-1);
			}

			bool merge = false;
			// printf("\n Min distance between %d %d is %f", cavityId-4, checkCavId-4, sqrt(minDist2));

			// bound on the max path
			double refDist = 6;

			if (sqrt(minDist2) > refDist)
				merge = false;

			// if we are under a threshold distance a full check is needed
			// check if one can move freely inside the small probe surface between the two nearest bgps
			// while at the mean time never going out the 1.4 surface.
			else
			{
				int maxMoves = (int)(refDist*delphi->scale + 0.5);
				int moves = maxMoves;
				v1 = vec1->at(winner1);
				v2 = vec2->at(winner2);

				bool debug = false;
				//if (cavityId-4 == 6 && checkCavId-4 == 17)
				//	debug = true;

				// printf("\n %d %d %d -> %d %d %d",v1[0],v1[1],v1[2],v2[0],v2[1],v2[2]);
				merge = innerFloodFill(v1,v2,st1,st2,moves,debug);
				double geodetic = moves*delphi->side;
				
				if (!merge)
				{
				//	printf("\n Non merging, num moves %d",moves);
				}
				if (merge && geodetic > 1.5*sqrt(minDist2))
				{
				// printf("\n Geodetic is %f and max dist is %f rejecting link",geodetic,1.5*sqrt(minDist2));
					merge = false;
				}

				//if (merge)
				//	printf("\n Merged in %d moves with geodetic %f",moves,geodetic);
			}
			
			// these two new cavities need to be linked
			// because they were linked by small probe surface accessibility
			if (merge)
			{
				// merge all the grid points
				vector <int*>::iterator it_vec2;

				for (it_vec2 = vec2->begin(); it_vec2 != vec2->end(); it_vec2++)
				{
					int *v = *it_vec2;
					vec1->push_back(v);
				}
					
				// pointers are not deleted
				vec2->clear();

				// sum volumes of the aggregated cavities
				delphi->cavitiesSize[ cavityId - STATUS_FIRST_CAV ] += delphi->cavitiesSize[ checkCavId - STATUS_FIRST_CAV ];
				// nullify the volume to flag it
				delphi->cavitiesSize[ checkCavId - STATUS_FIRST_CAV ] = 0;
				// assure it is active the merged cavity
				delphi->cavitiesFlag[ cavityId - STATUS_FIRST_CAV ] = false;
			}

		}
	}

	// disposeAtomsMap();

	sizeIt = delphi->cavitiesSize.begin();
	flagIt = delphi->cavitiesFlag.begin();
	it = delphi->cavitiesVec->begin();

	int numRemoved = 0;

	// remove all the null cavities
	while (1)
	{
		if (it == delphi->cavitiesVec->end())
			break;

		vector<int*> *vec1 = (*it);

		// the merged cavities are removed
		if (vec1->size() == 0)
		{
			numRemoved++;
			// delete vector of pointers to points
			delete vec1;
			vec1 = NULL;
			
			// delete entry
			delphi->cavitiesVec->erase(it);
			
			// double check it is the right guy
			if (*sizeIt != 0)
			{
				cout << endl << ERR << "Inconsistent volume aggregation";
				exit(-1);
			}

			// delete size (already aggregated)
			delphi->cavitiesSize.erase(sizeIt);
			// delete flag
			delphi->cavitiesFlag.erase(flagIt);

			// restart iterators
			it = delphi->cavitiesVec->begin();
			sizeIt = delphi->cavitiesSize.begin();
			flagIt = delphi->cavitiesFlag.begin();
		}
		else
		{
			it++;
			sizeIt++;
			flagIt++;
		}
	}

	int new_id = STATUS_FIRST_CAV;

	// renumber cavities accordingly 
	for (it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++)
	{
		vector<int*>::iterator inner;
		vector<int*> *cav = (*it);

		for (inner = cav->begin(); inner != cav->end(); inner++)
		{
			int *v = *inner;
			// conserve the support cavity coding
			
			if (readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,v[0],v[1],v[2],NX,NY,NZ) > 0)
				writeBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,new_id,v[0],v[1],v[2],NX,NY,NZ);
			else
				writeBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,-new_id,v[0],v[1],v[2],NX,NY,NZ);
		}
		new_id++;
	}
	return numRemoved;
}


/* old method
int Surface::linkCavities()
{
	int *tempStatus = delphi->tempStatus;
	int *status = delphi->status;

	if (tempStatus == NULL)
	{
		cout << endl << WARN << "I cannot reconstruct links if I have not reference status";
		return 0;
	}
		
	// analyze each cavity before and after fill
	vector<vector<int*>*>::iterator it;
	int i=0;
	
	int NX = delphi->nx;
	int NY = delphi->ny;	
	int NZ = delphi->nz;
		
	// cavity under link analysis
	int cavityId = STATUS_FIRST_CAV;

	vector<double>::iterator sizeIt;
	vector<bool>::iterator flagIt;

	for (it=delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++, cavityId++)
	{
		if (delphi->cavitiesFlag[cavityId-STATUS_FIRST_CAV] == true)
			continue;

		vector<int*> *vec1 = (*it);

		if (vec1->size() == 0)
			continue;

		// get the first point in cavity
		int *v = *(vec1->begin());
			
		// cavity id in the old map
		int oldCavityId = TEMP_STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY);

		// paired cavity under analysis
		int checkCavId = cavityId+1;

		// check if, for the all the other current cavities, exist a point
		// which had the same id before filtering. If yes a link is estabilished.
		for (checkCavId = cavityId+1;  (checkCavId-STATUS_FIRST_CAV) < (int)delphi->cavitiesVec->size(); checkCavId++)
		{
			if (delphi->cavitiesFlag[ checkCavId - STATUS_FIRST_CAV ] == true)
				continue;
				
			vector<int*> *vec2 = delphi->cavitiesVec->at(checkCavId - STATUS_FIRST_CAV);

			// this may happen if this cavity have been agregated
			if (vec2->size() == 0)
				continue;
				
			// take a random point and check the previous map
			int *v = *(vec2->begin());
			int oldCavityId2 = TEMP_STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY);

			// these two new cavities need to be linked
			// because they were linked before filtering
			if (oldCavityId == oldCavityId2)
			{
				// printf("\n%sLinking cavities (%d,%d)", INFO, checkCavId-STATUS_FIRST_CAV, cavityId-STATUS_FIRST_CAV);
				// printf("\n old value %d",oldCavityId);
				// getchar();
					
				// merge the grid points
				vector<int*>::iterator it_vec2;

				for (it_vec2 = vec2->begin(); it_vec2 != vec2->end(); it_vec2++)
				{
					int *v = *it_vec2;
					vec1->push_back(v);
				}
					
				// pointers are not deleted
				vec2->clear();

				// sum volumes of the agregate cavities
				delphi->cavitiesSize[ cavityId-STATUS_FIRST_CAV ] += delphi->cavitiesSize[ checkCavId - STATUS_FIRST_CAV ];
				// nullify the volume to flag it
				delphi->cavitiesSize[ checkCavId - STATUS_FIRST_CAV ] = 0;
				// assure it is active the merged cavity
				delphi->cavitiesFlag[ cavityId-STATUS_FIRST_CAV ] = false;
			}

		}
	}

	sizeIt = delphi->cavitiesSize.begin();
	flagIt = delphi->cavitiesFlag.begin();
	it = delphi->cavitiesVec->begin();

	int numRemoved = 0;

	// remove all the null cavities
	while (1)
	{
		if (it == delphi->cavitiesVec->end())
			break;

		vector<int*> *vec1 = (*it);

		// the merged cavities are removed
		if (vec1->size() == 0)
		{
			numRemoved++;
			// delete vector of pointers to points
			delete vec1;
			vec1 = NULL;
			
			// delete entry
			delphi->cavitiesVec->erase(it);
			
			// double check it is the right guy
			if (*sizeIt != 0)
			{
				cout << endl << ERR << "Inconsistent volume agregation";
				exit(-1);
			}

			// delete size (already aggregated)
			delphi->cavitiesSize.erase(sizeIt);
			// delete flag
			delphi->cavitiesFlag.erase(flagIt);

			// restart iterators
			it = delphi->cavitiesVec->begin();
			sizeIt = delphi->cavitiesSize.begin();
			flagIt = delphi->cavitiesFlag.begin();
		}
		else
		{
			it++;
			sizeIt++;
			flagIt++;
		}
	}

	int new_id = STATUS_FIRST_CAV;
		
	// renumber cavities accordingly 
	for (it = delphi->cavitiesVec->begin(); it != delphi->cavitiesVec->end(); it++)
	{
		vector<int*>::iterator inner;
		vector<int*> *cav = (*it);

		for (inner = cav->begin(); inner != cav->end(); inner++)
		{
			int *v = *inner;
			// conserve the support cavity coding
			if (STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) > 0)
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = new_id;
			else
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = -new_id;
		}
		new_id++;
	}
	return numRemoved;
}
*/


Surface &Surface::operator -= (Surface &surf2)
{
	auto chrono_start = chrono::high_resolution_clock::now();
	
	if (!optimizeGrids)
		difference (&surf2);
	else
		differenceWithBilevelStatusMap (&surf2);

	auto chrono_end = chrono::high_resolution_clock::now();
	
	chrono::duration<double> diff1_time = chrono_end - chrono_start;
	cout << endl << INFO << "Diff. Step 1 ";
	printf ("%.4e [s]", diff1_time.count());
	cout.flush();
	

	///////////////////// connolly filter /////////////////////////////////
	setProbeRadius(1.4);
	// setProbeRadius(0.3);
	
	auto chrono2_start = chrono::high_resolution_clock::now();

	if (!optimizeGrids)
	{
		filterCavities(true);
		cav2out();
	}
	else
	{
		filterCavitiesWithBilevelStatusMap(true);
		cav2outWithBilevelStatusMap();
	}

	auto chrono2_end = chrono::high_resolution_clock::now();

	chrono::duration<double> diff2_time = chrono2_end - chrono2_start;
	cout << endl << INFO << "Diff. Step 2 ";
	printf ("%.4e [s]", diff2_time.count());
	cout.flush();
	

	auto chrono3_start = chrono::high_resolution_clock::now();
	
	///////////////////// get filtered cavities ///////////////////////////
	int ff;
	if (!optimizeGrids)
		ff = getCavities(STATUS_POINT_OUT);
	else
		ff = getCavitiesWithBilevelStatusMap(STATUS_POINT_OUT);
	if (!ff)
		// cout << endl << WARN <<  "Zero difference after filtration!";
		return *this;

	auto chrono3_end = chrono::high_resolution_clock::now();

	chrono::duration<double> diff3_time = chrono3_end - chrono3_start;
	cout << endl << INFO << "Diff. Step 3 ";
	printf ("%.4e [s]", diff3_time.count());
	cout.flush();


	return *this;  
}


bool Surface::innerFloodFill(int *start, int *target, int *tempStatus, int *tempStatus2, int &maxMoves, bool debug)
{
	int *status = delphi->status;

	// ix,iy,iz and integer distance between the origin and that point
	queue<pair<pair<int,int>,pair<int,int>>> myqueue;
	pair<pair<int,int>,pair<int,int>> quartet;

	int cix,ciy,ciz,dist;

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	int ix = start[0];
	int iy = start[1];
	int iz = start[2];

	#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
	int octree_side_size = 2;
	bool found_label = 0;

	while (!found_label)
	{
		if (octree_side_size < NX ||
			octree_side_size < NY ||
			octree_side_size < NZ)
		{
			octree_side_size *= 2;
		}
		else
		{
			found_label = 1;
		}
	}
	Octree<bool> visited(octree_side_size, false);
	#else
	unsigned int **compressed_bilevel_visited = allocate32xCompressedBilevelGridCells(NX, NY, NZ);
	#endif

	myqueue.push( pair<pair<int,int>,pair<int,int>> (pair<int,int>(ix,iy), pair<int,int>(iz,0)));

	while (myqueue.size() != 0)
	{
		quartet = myqueue.front();
		myqueue.pop();

		cix = quartet.first.first;
		ciy = quartet.first.second;
		ciz = quartet.second.first;
		dist = quartet.second.second;

		if (dist+1 > maxMoves)
		{
			/*
			FILE *fp;
			if (debug) fp = fopen("debugFlood.txt","w");
			printf("\n Target not found, too many moves");
			// restore old status
			for (unsigned int ll=0; ll<visitedPoints.size(); ll++)
			{
				pair<int,pair<int,int>> triplet = visitedPoints[ll];
				// STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY) = -STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY);
				if (debug) fprintf(fp,"\n %f %f %f", X[triplet.first],Y[triplet.second.first], Z[triplet.second.second]);
			}
			if (debug) fclose(fp);
			*/
			#if defined(OPTIMIZE_INNER_FLOOD_FILLING)
			delete32xCompressedBilevelGridCells(compressed_bilevel_visited, NX, NY, NZ);
			deleteVector<unsigned int *>(compressed_bilevel_visited);
			#endif

			maxMoves = dist;
			return false;
		}

		// target obtained
		if ((cix == target[0]) && (ciy == target[1]) && (ciz == target[2]))
		{
			/*
			printf("\n Target acquired!");
			FILE *fp;
			if (debug) fp = fopen("debugFlood.txt","w");

			// restore old status
			for (unsigned int ll=0; ll < visitedPoints.size(); ll++)
			{
				pair<int,pair<int,int>> triplet = visitedPoints[ll];
				//STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY) = -STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY);
				if (debug) fprintf(fp,"\n %f %f %f", X[triplet.first],Y[triplet.second.first],Z[triplet.second.second]);
			}
			if (debug) fclose(fp);
			*/
			#if defined(OPTIMIZE_INNER_FLOOD_FILLING)
			delete32xCompressedBilevelGridCells(compressed_bilevel_visited, NX, NY, NZ);
			deleteVector<unsigned int *>(compressed_bilevel_visited);
			#endif

			maxMoves = dist;
			return true;
		}

		// already visited?
		#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
		if (visited.at(cix,ciy,ciz) == true)
			continue;
		else
			visited.set(cix,ciy,ciz, true);
		#else
		if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy,ciz, NX, NY, NZ) == true)
			continue;
		else
			write32xCompressedBilevelGrid(compressed_bilevel_visited, false, true, cix,ciy,ciz, NX, NY, NZ);
		#endif

		// can still move, go on
		int newDist = dist+1;
		bool notVisited,inner,isAccessible;

		// double p[3];

		if (cix+1 < NX)
		{
			// p[0] = X[cix+1];
			// p[1] = Y[ciy];
			// p[2] = Z[ciz];

			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix+1,ciy,ciz) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix+1,ciy,ciz, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				// if (TEMP_STATUSMAP2(cix+1,ciy,ciz,NX,NY) == STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix+1,ciy,ciz,NX,NY) == STATUS_POINT_OUT)
				const int ts = read3DVector<int>(tempStatus2,cix+1,ciy,ciz,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				// inner = TEMP_STATUSMAP(cix+1,ciy,ciz,NX,NY) == STATUS_POINT_INSIDE || STATUSMAP(cix+1,ciy,ciz,NX,NY) >= STATUS_FIRST_CAV;

				const int tss = read3DVector<int>(status,cix+1,ciy,ciz,NX,NY,NZ);
				inner = ((read3DVector<int>(tempStatus,cix+1,ciy,ciz,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix + 1;
					quartet.first.second = ciy;
					quartet.second.first = ciz;

					// if in cavity/pocket does not count the move, the flow is still trying
					// to percolate through the linking surf but still the linking surf path connecting the two cavities
					// has not been found
					// if (STATUSMAP(cix+1,ciy,ciz,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;

					myqueue.push(quartet);
				}
			}
		}

		if (ciy+1 < NY)
		{
			// p[0] = X[cix];
			// p[1] = Y[ciy+1];
			// p[2] = Z[ciz];

			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix,ciy+1,ciz) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy+1,ciz, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				// if (TEMP_STATUSMAP2(cix,ciy+1,ciz,NX,NY) == STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix,ciy+1,ciz,NX,NY) == STATUS_POINT_OUT)
				const int ts = read3DVector<int>(tempStatus2,cix,ciy+1,ciz,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				// inner = TEMP_STATUSMAP(cix,ciy+1,ciz,NX,NY) == STATUS_POINT_INSIDE || STATUSMAP(cix,ciy+1,ciz,NX,NY) >= STATUS_FIRST_CAV;

				const int tss = read3DVector<int>(status,cix,ciy+1,ciz,NX,NY,NZ);
				inner = ((read3DVector<int>(tempStatus,cix,ciy+1,ciz,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy+1;
					quartet.second.first = ciz;

					// if (STATUSMAP(cix,ciy+1,ciz,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);
				}
			}
		}

		if (ciz+1 < NZ)
		{
			// p[0] = X[cix];
			// p[1] = Y[ciy];
			// p[2] = Z[ciz+1];

			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix,ciy,ciz+1) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy,ciz+1, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				// if (TEMP_STATUSMAP2(cix,ciy,ciz+1,NX,NY) == STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix,ciy,ciz+1,NX,NY) == STATUS_POINT_OUT)
				const int ts = read3DVector<int>(tempStatus2,cix,ciy,ciz+1,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				// inner = TEMP_STATUSMAP(cix,ciy,ciz+1,NX,NY) == STATUS_POINT_INSIDE || STATUSMAP(cix,ciy,ciz+1,NX,NY) >= STATUS_FIRST_CAV;

				const int tss = read3DVector<int>(status,cix,ciy,ciz+1,NX,NY,NZ);
				inner = ((read3DVector<int>(tempStatus,cix,ciy,ciz+1,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy;
					quartet.second.first = ciz+1;

					// if (STATUSMAP(cix,ciy,ciz+1,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);
				}
			}
		}


		if (cix-1 >= 0)
		{
			// p[0] = X[cix-1];
			// p[1] = Y[ciy];
			// p[2] = Z[ciz];

			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix-1,ciy,ciz) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix-1,ciy,ciz, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				// if (TEMP_STATUSMAP2(cix-1,ciy,ciz,NX,NY) == STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix-1,ciy,ciz,NX,NY) == STATUS_POINT_OUT)
				const int ts = read3DVector<int>(tempStatus2,cix-1,ciy,ciz,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				// inner = TEMP_STATUSMAP(cix-1,ciy,ciz,NX,NY) == STATUS_POINT_INSIDE || STATUSMAP(cix-1,ciy,ciz,NX,NY) >= STATUS_FIRST_CAV;
				const int tss = read3DVector<int>(status,cix-1,ciy,ciz,NX,NY,NZ);
				inner = ((read3DVector<int>(tempStatus,cix-1,ciy,ciz,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix-1;
					quartet.first.second = ciy;
					quartet.second.first = ciz;

					// if (STATUSMAP(cix-1,ciy,ciz,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);
				}
			}
		}

		if (ciy-1 >= 0)
		{
			// p[0] = X[cix];
			// p[1] = Y[ciy-1];
			// p[2] = Z[ciz];

			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix,ciy-1,ciz) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy-1,ciz, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				// if (TEMP_STATUSMAP2(cix,ciy-1,ciz,NX,NY) == STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix,ciy-1,ciz,NX,NY) == STATUS_POINT_OUT)
				const int ts = read3DVector<int>(tempStatus2,cix,ciy-1,ciz,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				// inner = TEMP_STATUSMAP(cix,ciy-1,ciz,NX,NY) == STATUS_POINT_INSIDE || STATUSMAP(cix,ciy-1,ciz,NX,NY) >= STATUS_FIRST_CAV;

				const int tss = read3DVector<int>(status,cix,ciy-1,ciz,NX,NY,NZ);
				inner = ((read3DVector<int>(tempStatus,cix,ciy-1,ciz,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy-1;
					quartet.second.first = ciz;

					// if (STATUSMAP(cix,ciy-1,ciz,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);
				}
			}
		}

		if (ciz-1 >= 0)
		{
			// p[0] = X[cix];
			// p[1] = Y[ciy];
			// p[2] = Z[ciz-1];

			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix,ciy,ciz-1) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy,ciz-1, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				// if (TEMP_STATUSMAP2(cix,ciy,ciz-1,NX,NY) == STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix,ciy,ciz-1,NX,NY) == STATUS_POINT_OUT)
				const int ts = read3DVector<int>(tempStatus2,cix,ciy,ciz-1,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				// inner = TEMP_STATUSMAP(cix,ciy,ciz-1,NX,NY) == STATUS_POINT_INSIDE || STATUSMAP(cix,ciy,ciz-1,NX,NY) >= STATUS_FIRST_CAV;

				const int tss = read3DVector<int>(status,cix,ciy,ciz-1,NX,NY,NZ);
				inner = ((read3DVector<int>(tempStatus,cix,ciy,ciz-1,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy;
					quartet.second.first = ciz-1;

					// if (STATUSMAP(cix,ciy,ciz-1,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;

					myqueue.push(quartet);
				}
			}
		}
	}

	/* FILE *fp;
	if (debug) fp = fopen("debugFlood.txt","w");

	// restore old status
	for (unsigned int ll=0; ll<visitedPoints.size(); ll++)
	{
		pair<int,pair<int,int>>triplet = visitedPoints[ll];
		// STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY) = -STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY);
		if (debug) fprintf(fp,"\n %f %f %f", X[triplet.first],Y[triplet.second.first],Z[triplet.second.second]);
	}
	if (debug) fclose(fp);
	printf("\n Target not found!");
	*/

	#if defined(OPTIMIZE_INNER_FLOOD_FILLING)
	delete32xCompressedBilevelGridCells(compressed_bilevel_visited, NX, NY, NZ);
	deleteVector<unsigned int *>(compressed_bilevel_visited);
	#endif

	maxMoves = dist;
	return false;
}


bool Surface::innerFloodFill(int *start, int *target, int **bilevel_tempStatus, int **bilevel_tempStatus2, int &maxMoves, bool debug)
{
	int **bilevel_status = delphi->bilevel_status;

	// ix,iy,iz and integer distance between the origin and that point
	queue<pair<pair<int,int>,pair<int,int>>> myqueue;
	pair<pair<int,int>,pair<int,int>> quartet;

	int cix, ciy, ciz, dist;
	
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	int ix = start[0];
	int iy = start[1];
	int iz = start[2];

	#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
	int octree_side_size = 2;
	bool found_label = 0;

	while (!found_label)
	{
		if (octree_side_size < NX ||
			octree_side_size < NY ||
			octree_side_size < NZ)
		{
			octree_side_size *= 2;
		}
		else
		{
			found_label = 1;
		}
	}
	Octree<bool> visited(octree_side_size, false);
	#else
	unsigned int **compressed_bilevel_visited = allocate32xCompressedBilevelGridCells(NX, NY, NZ);
	#endif

	myqueue.push( pair<pair<int,int>,pair<int,int>> (pair<int,int>(ix,iy),pair<int,int>(iz,0)) );

	while (myqueue.size() != 0)
	{
		quartet = myqueue.front();
		myqueue.pop();
		
		cix = quartet.first.first;
		ciy = quartet.first.second;
		ciz = quartet.second.first;
		dist = quartet.second.second;

		if (dist+1 > maxMoves)
		{
			#if defined(OPTIMIZE_INNER_FLOOD_FILLING)
			delete32xCompressedBilevelGridCells(compressed_bilevel_visited, NX, NY, NZ);
			deleteVector<unsigned int *>(compressed_bilevel_visited);
			#endif

			maxMoves = dist;
			return false;
		}

		// target obtained
		if ((cix == target[0]) && (ciy == target[1]) && (ciz == target[2]))
		{
			#if defined(OPTIMIZE_INNER_FLOOD_FILLING)
			delete32xCompressedBilevelGridCells(compressed_bilevel_visited, NX, NY, NZ);
			deleteVector<unsigned int *>(compressed_bilevel_visited);
			#endif

			maxMoves = dist;
			return true;
		}

		// already visited?
		#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
		if (visited.at(cix,ciy,ciz) == true)
			continue;
		else
			visited.set(cix,ciy,ciz, true);
		#else
		if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy,ciz, NX, NY, NZ) == true)
			continue;
		else
			write32xCompressedBilevelGrid(compressed_bilevel_visited, false, true, cix,ciy,ciz, NX, NY, NZ);
		#endif

		// can still move, go on
		int newDist = dist+1;
		bool notVisited,inner,isAccessible;

		// double p[3];

		if (cix+1 < NX)
		{
			// p[0] = X[cix+1];
			// p[1] = Y[ciy];
			// p[2] = Z[ciz];

			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix+1,ciy,ciz) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix+1,ciy,ciz, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				const int ts = readBilevelGrid<int>(bilevel_tempStatus2,STATUS_POINT_TEMPORARY_OUT,cix+1,ciy,ciz,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				const int tss = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix+1,ciy,ciz,NX,NY,NZ);
				inner = ((readBilevelGrid<int>(bilevel_tempStatus,STATUS_POINT_TEMPORARY_OUT,cix+1,ciy,ciz,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);
				
				if (isAccessible && inner)
				{
					quartet.first.first = cix+1;
					quartet.first.second = ciy;
					quartet.second.first = ciz;

					// if in cavity/pocket does not count the move, the flow is still trying
					// to percolate through the linking surf but still the linking surf path connecting the two cavities
					// has not been found
					// if (STATUSMAP(cix+1,ciy,ciz,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
						
					myqueue.push(quartet);
				}
			}
		}

		if (ciy+1 < NY)
		{
			// p[0] = X[cix];
			// p[1] = Y[ciy+1];
			// p[2] = Z[ciz];

			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix,ciy+1,ciz) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy+1,ciz, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				const int ts = readBilevelGrid<int>(bilevel_tempStatus2,STATUS_POINT_TEMPORARY_OUT,cix,ciy+1,ciz,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				const int tss = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix,ciy+1,ciz,NX,NY,NZ);
				inner = ((readBilevelGrid<int>(bilevel_tempStatus,STATUS_POINT_TEMPORARY_OUT,cix,ciy+1,ciz,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy+1;
					quartet.second.first = ciz;

					// if (STATUSMAP(cix,ciy+1,ciz,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);
				}
			}
		}

		if (ciz+1 < NZ)
		{
			// p[0] = X[cix];
			// p[1] = Y[ciy];
			// p[2] = Z[ciz+1];
			
			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix,ciy,ciz+1) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy,ciz+1, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				const int ts = readBilevelGrid<int>(bilevel_tempStatus2,STATUS_POINT_TEMPORARY_OUT,cix,ciy,ciz+1,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				const int tss = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix,ciy,ciz+1,NX,NY,NZ);
				inner = ((readBilevelGrid<int>(bilevel_tempStatus,STATUS_POINT_TEMPORARY_OUT,cix,ciy,ciz+1,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy;
					quartet.second.first = ciz+1;

					// if (STATUSMAP(cix,ciy,ciz+1,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);
				}
			}
		}
		

		if (cix-1 >= 0)
		{
			// p[0] = X[cix-1];
			// p[1] = Y[ciy];
			// p[2] = Z[ciz];
			
			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix-1,ciy,ciz) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix-1,ciy,ciz, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				const int ts = readBilevelGrid<int>(bilevel_tempStatus2,STATUS_POINT_TEMPORARY_OUT,cix-1,ciy,ciz,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				const int tss = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix-1,ciy,ciz,NX,NY,NZ);
				inner = ((readBilevelGrid<int>(bilevel_tempStatus,STATUS_POINT_TEMPORARY_OUT,cix-1,ciy,ciz,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix-1;
					quartet.first.second = ciy;
					quartet.second.first = ciz;

					// if (STATUSMAP(cix-1,ciy,ciz,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);
				}
			}
		}

		if (ciy-1 >= 0)
		{
			// p[0] = X[cix];
			// p[1] = Y[ciy-1];
			// p[2] = Z[ciz];
		
			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix,ciy-1,ciz) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy-1,ciz, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				const int ts = readBilevelGrid<int>(bilevel_tempStatus2,STATUS_POINT_TEMPORARY_OUT,cix,ciy-1,ciz,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;
								
				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				const int tss = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix,ciy-1,ciz,NX,NY,NZ);
				inner = ((readBilevelGrid<int>(bilevel_tempStatus,STATUS_POINT_TEMPORARY_OUT,cix,ciy-1,ciz,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy-1;
					quartet.second.first = ciz;

					// if (STATUSMAP(cix,ciy-1,ciz,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);
				}
			}
		}

		if (ciz-1 >= 0)
		{
			// p[0] = X[cix];
			// p[1] = Y[ciy];
			// p[2] = Z[ciz-1];

			#if !defined(OPTIMIZE_INNER_FLOOD_FILLING)
			if (visited.at(cix,ciy,ciz-1) == false)
				notVisited = true;
			else
				notVisited = false;
			#else
			if (read32xCompressedBilevelGrid(compressed_bilevel_visited, false, cix,ciy,ciz-1, NX, NY, NZ) == false)
				notVisited = true;
			else
				notVisited = false;
			#endif

			if (notVisited)
			{
				const int ts = readBilevelGrid<int>(bilevel_tempStatus2,STATUS_POINT_TEMPORARY_OUT,cix,ciy,ciz-1,NX,NY,NZ);

				if (ts == STATUS_POINT_TEMPORARY_OUT || ts == STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				const int tss = readBilevelGrid<int>(bilevel_status,STATUS_POINT_TEMPORARY_OUT,cix,ciy,ciz-1,NX,NY,NZ);
				inner = ((readBilevelGrid<int>(bilevel_tempStatus,STATUS_POINT_TEMPORARY_OUT,cix,ciy,ciz-1,NX,NY,NZ) == STATUS_POINT_INSIDE) || tss >= STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy;
					quartet.second.first = ciz-1;

					// if (STATUSMAP(cix,ciy,ciz-1,NX,NY) >= STATUS_FIRST_CAV)
					if (tss >= STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;

					myqueue.push(quartet);
				}
			}
		}
	}

	/*
	FILE *fp;
	if (debug) fp = fopen("debugFlood.txt","w");

	// restore old status
	for (unsigned int ll=0; ll<visitedPoints.size(); ll++)
	{
		pair<int,pair<int,int>>triplet = visitedPoints[ll];
		// cSTATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY) = -STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY);
		if (debug) fprintf(fp,"\n %f %f %f", X[triplet.first],Y[triplet.second.first],Z[triplet.second.second]);
	}
	if (debug) fclose(fp);
	printf("\n Target not found!");
	*/

	#if defined(OPTIMIZE_INNER_FLOOD_FILLING)
	delete32xCompressedBilevelGridCells(compressed_bilevel_visited, NX, NY, NZ);
	deleteVector<unsigned int *>(compressed_bilevel_visited);
	#endif

	maxMoves = dist;
	return false;
}


bool Surface::isCompletelyOut(VERTEX_TYPE *pos)
{
	if (!optimizeGrids)
	{
		if (delphi->status == NULL)
		{
			cout << endl << ERR << "Cannot compute if a point is completely out without the status map";
			exit(-1);
		}
	}
	else
	{
		if (delphi->bilevel_status == NULL)
		{
			cout << endl << ERR << "Cannot compute if a point is completely out without the bilevel status map";
			exit(-1);
		}
	}
	int64_t nx = delphi->nx;
	int64_t ny = delphi->ny;
	int64_t nz = delphi->nz;

	int64_t ix = (int64_t)rintp((pos[0]-delphi->xmin)*delphi->scale);
	int64_t iy = (int64_t)rintp((pos[1]-delphi->ymin)*delphi->scale);
	int64_t iz = (int64_t)rintp((pos[2]-delphi->zmin)*delphi->scale);

	int test;
		
	for (int k=0; k<SHIFT_MAP; k++)
	{
		int64_t cx = ix+shift_map[k][0];
		int64_t cy = iy+shift_map[k][1];
		int64_t cz = iz+shift_map[k][2];

		if (cx >= nx || cy >= ny || cz >= nz || cx<0 || cy<0 || cz<0)
			continue;

		if (!optimizeGrids)
			test = read3DVector<int>(delphi->status,cx,cy,cz,nx,ny,nz);
		else
			test = readBilevelGrid<int>(delphi->bilevel_status,STATUS_POINT_TEMPORARY_OUT,cx,cy,cz,nx,ny,nz);

		// ok
		if (test == STATUS_POINT_TEMPORARY_OUT || test == STATUS_POINT_OUT)
			continue;
		// at least one is not completely out
		else
			return false;
	}
	return true;
}


void Surface::triangulationPointsAreCompletelyOut(Surface *surf, vector<bool> &results)
{
	const double normalProbeOffset = 0.10;
	const double normalNormEps2 = 1e-20;

	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int nv = (int)vertList.size();
	const bool hasNormals = ((int)normalsList.size() == nv);

	for (unsigned int i=0; i<nv; i++)
	{
		VERTEX_TYPE *v = vertList[i];
		bool isOut = surf->isCompletelyOut(v);

		if (isOut && hasNormals)
		{
			VERTEX_TYPE *n = normalsList[i];
			double nn2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
			if (nn2 > normalNormEps2)
			{
				double invN = 1.0 / sqrt(nn2);
				VERTEX_TYPE shifted[3];
				shifted[0] = v[0] + (VERTEX_TYPE)(normalProbeOffset * n[0] * invN);
				shifted[1] = v[1] + (VERTEX_TYPE)(normalProbeOffset * n[1] * invN);
				shifted[2] = v[2] + (VERTEX_TYPE)(normalProbeOffset * n[2] * invN);
				isOut = surf->isCompletelyOut(shifted);
			}
		}

		results.push_back(isOut);
	}
	#else
	int nv = (int)(vertList.size() / 3.);
	const bool hasNormals = ((int)normalsList.size() == nv*3);

	for (unsigned int i=0; i<nv; i++)
	{
		VERTEX_TYPE *v = &vertList[i*3];
		bool isOut = surf->isCompletelyOut(v);

		if (isOut && hasNormals)
		{
			VERTEX_TYPE *n = &normalsList[i*3];
			double nn2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
			if (nn2 > normalNormEps2)
			{
				double invN = 1.0 / sqrt(nn2);
				VERTEX_TYPE shifted[3];
				shifted[0] = v[0] + (VERTEX_TYPE)(normalProbeOffset * n[0] * invN);
				shifted[1] = v[1] + (VERTEX_TYPE)(normalProbeOffset * n[1] * invN);
				shifted[2] = v[2] + (VERTEX_TYPE)(normalProbeOffset * n[2] * invN);
				isOut = surf->isCompletelyOut(shifted);
			}
		}

		results.push_back(isOut);
	}
	#endif
}


void Surface::saveSelectedPoints(vector<bool> &selection, const char *file, const char *file2)
{
	char fullName[1024];
	sprintf(fullName, "%s%s", rootFile.c_str(),file);
	FILE *fp = fopen(fullName, "w");
	FILE *fp2 = NULL;

	if (fp == NULL)
	{
		cout << endl << WARN << "Cannot write file " << fullName;
		return;
	}

	bool legacyDualOutput = (file2 != NULL && file2[0] != '\0');
	if (legacyDualOutput)
	{
		sprintf(fullName, "%s%s", rootFile.c_str(),file2);
		if (normalsList.size() != 0)
			fp2 = fopen(fullName,"w");
		if (normalsList.size() != 0 && fp2 == NULL)
			cout << endl << WARN << "Cannot write file " << fullName;
	}

	int count = 0;
	for (unsigned int i=0; i<selection.size(); i++)
	{
		if (selection[i])
			count++;
	}
	// In legacy mode keep the XYZ file VMD-friendly.
	// In single-output mode keep the same header and append normals per point when available.
	fprintf(fp, "%d\n", count);
	fprintf(fp, "pocket_entrance_as_a_set_of_dummy_atoms_by_NanoShaper\n");

	for (unsigned int i=0; i<selection.size(); i++)
	{
		if (selection[i])
		{
			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			VERTEX_TYPE *v = vertList[i];
			VERTEX_TYPE *n = NULL;

			if (normalsList.size() != 0)
				n = normalsList[i];
			#else
			VERTEX_TYPE *v = &vertList[i*3];
			VERTEX_TYPE *n = NULL;

			if (normalsList.size() != 0)
				n = &normalsList[i*3];
			#endif

			if (legacyDualOutput || n == NULL)
				fprintf(fp, "C %f %f %f\n", v[0],v[1],v[2]);
			else
				fprintf(fp, "C %f %f %f %f %f %f\n", v[0],v[1],v[2],n[0],n[1],n[2]);

			if (fp2 != NULL && n != NULL)
				fprintf(fp2, "%f %f %f %f %f %f\n", v[0],v[1],v[2],n[0],n[1],n[2]);
		}
	}
	fclose(fp);

	if (fp2 != NULL)
		fclose(fp2);
}


bool Surface::vdwAccessible(double *pos, int &winner)
{
	double minDist = INFINITY;
	winner = -1;
	
	int64_t ix = (int64_t)rintp((pos[0]-gxmin)*gscale);
	int64_t iy = (int64_t)rintp((pos[1]-gymin)*gscale);
	int64_t iz = (int64_t)rintp((pos[2]-gzmin)*gscale);

	// get the nearest atom and says if it is in or out
	// the signed distance from the point p is ||p-c||^2-r^2 where c is the center
	// of the atom and r is the radius. The minimum signed distance wins and if this
	// negative we are inside (false) and if it is positive we are outside (true)
	for (int k=0; k<SHIFT_MAP; k++)
	{
		int64_t cx = ix+shift_map[k][0];
		int64_t cy = iy+shift_map[k][1];
		int64_t cz = iz+shift_map[k][2];

		// multidielectric map is square 
		if (cx >= ggrid || cy >= ggrid || cz >= ggrid || cx<0 || cy<0 || cz<0)
			continue;

		int num_atoms = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ].size();

		for (int j=0; j<num_atoms; j++)
		{
			int atom_index = gridMultiMap[ cz*(ggrid*ggrid) + cy*ggrid + cx ][j];

			double signed_dist = 0;
			// DIST2(signed_dist,delphi->atoms[atom_index]->pos,pos)
			// double rad = delphi->atoms[atom_index]->radius;
			DIST2(signed_dist,delphi->atoms[atom_index].pos,pos)
			double rad = delphi->atoms[atom_index].radius;
			signed_dist -= rad*rad;
			
			if (signed_dist < minDist)
			{
				minDist = signed_dist;
				winner = atom_index;
			}
		}	
	}
	
	if (winner == -1)
	{
		return true;
	}
	if (minDist < 0)
		return false;
	else
		return true;
}


double Surface::saveTriSubSet(char *triSubset, vector<bool> &results, bool revert)
{
	double area = 0;
	
	// previous vertex index, new vertex index
	// -1 if the vertex is removed
	map<int,int> indicesMap;

	// reduced list of triangles, vertices and normals
	vector<int> triListRed;
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	vector<VERTEX_TYPE*> vertListRed;
	vector<VERTEX_TYPE*> normalsListRed;
	#else
	vector<VERTEX_TYPE> vertListRed;
	vector<VERTEX_TYPE> normalsListRed;
	#endif


	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	int nv = (int)vertList.size();
	#else
	int nv = (int)(vertList.size() / 3.);
	#endif

	if (triList.size() == 0)
	{
		cout << ERR << "Cannot save a reduced set of triangles if triangulation is absent";
		return 0.;
	}

	if (results.size() != nv)
	{
		cout << ERR << "Cannot save a reduced set of triangles if number of vertices is unconsistent between current mesh and vertices flags";
		return 0.;
	}

	for (unsigned int i=0; i<results.size(); i++)
	{
		if (results[i])
		{
			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			indicesMap.insert( pair<int,int> (i, (int)vertListRed.size()) );

			vertListRed.push_back(vertList[i]);
			#else
			indicesMap.insert( pair<int,int> (i, (int)(vertListRed.size() / 3.)) );
			vertListRed.push_back(vertList[ i*3+0 ]);
			vertListRed.push_back(vertList[ i*3+1 ]);
			vertListRed.push_back(vertList[ i*3+2 ]);
			#endif

			#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
			if (normalsList.size() != 0)
				normalsListRed.push_back(normalsList[i]);
			#else
			if (normalsList.size() != 0)
			{
				normalsListRed.push_back(normalsList[ i*3+0 ]);
				normalsListRed.push_back(normalsList[ i*3+1 ]);
				normalsListRed.push_back(normalsList[ i*3+2 ]);
			}
			#endif
		}
		else
			indicesMap.insert(pair<int,int>(i,-1));
	}

	int nt = (int)(triList.size() / 3.);

	for (unsigned int i=0; i<nt; i++)
	{
		map<int,int>::iterator it;
		
		it = indicesMap.find(triList[ i*3 ]);
		int i1 = it->second;
		if (i1 == -1)
			continue;
		
		it = indicesMap.find(triList[ i*3+1 ]);
		int i2 = it->second;
		if (i2 == -1)
			continue;

		it = indicesMap.find(triList[ i*3+2 ]);
		int i3 = it->second;
		if (i3 == -1)
			continue;

		VERTEX_TYPE a=0, b=0, c=0;

		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		VERTEX_TYPE *v0 = vertListRed[i1];
		VERTEX_TYPE *v1 = vertListRed[i2];
		VERTEX_TYPE *v2 = vertListRed[i3];
		#else
		VERTEX_TYPE *v0 = &vertListRed[ i1*3 ];
		VERTEX_TYPE *v1 = &vertListRed[ i2*3 ];
		VERTEX_TYPE *v2 = &vertListRed[ i3*3 ];
		#endif

		DIST(a,v0,v1);
		DIST(b,v0,v2);
		DIST(c,v1,v2);

		double ttt = (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c);

		if (ttt > 0)
			area += 0.25*sqrt(ttt);

		triListRed.push_back(i1);
		triListRed.push_back(i2);
		triListRed.push_back(i3);
	}

	auto chrono_start = chrono::high_resolution_clock::now();

	int format = deduceFormat();
	saveMesh(format, revert, triSubset, vertListRed, triListRed, normalsListRed);

	auto chrono_end = chrono::high_resolution_clock::now();

	chrono::duration<double> saveMesh_time = chrono_end - chrono_start;

	cout << endl << INFO << "Outputting mesh time (in saveTriSubSet()) ";
	printf ("%.4e [s]", saveMesh_time.count());

	return area;
}


int Surface::getTriangulatedNumVertices(void) const
{
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	return (int)vertList.size();
	#else
	return (int)(vertList.size() / 3.);
	#endif
}


int Surface::getTriangulatedNumTriangles(void) const
{
	return (int)(triList.size() / 3.);
}


void Surface::getTriangulatedVertex(int idx, double out[3]) const
{
	#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
	out[0] = vertList[idx][0];
	out[1] = vertList[idx][1];
	out[2] = vertList[idx][2];
	#else
	out[0] = vertList[idx*3+0];
	out[1] = vertList[idx*3+1];
	out[2] = vertList[idx*3+2];
	#endif
}


void Surface::getTriangulatedTriangle(int triIdx, int &v0, int &v1, int &v2) const
{
	v0 = triList[triIdx*3+0];
	v1 = triList[triIdx*3+1];
	v2 = triList[triIdx*3+2];
}


bool Surface::nearestAtomForPoint(double *pos, int &winner)
{
	return vdwAccessible(pos, winner);
}


int Surface::getNumLoadedAtoms(void) const
{
	return (int)delphi->atoms.size();
}

int Surface::getLoadedAtomOutputSerial(int atomIdx) const
{
	if (atomIdx < 0 || atomIdx >= (int)delphi->atoms.size())
		return atomIdx + 1;

	const int serial = delphi->atoms[atomIdx].ai.getSerial();
	if (serial > 0)
		return serial;

	return atomIdx + 1;
}


void Surface::prepareAtomsMap(void)
{
	buildAtomsMap();
}


void Surface::clearAtomsMap(void)
{
	disposeAtomsMap();
}


///////////////////////////////////////////// Visualisation functions //////////////////////////////////////////////

#if defined(USE_VIS_TOOLS)

#ifndef NO_OPENGL
void Surface::openWindow (Vis *vis)
{
	glutInitDisplayMode (GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition (0, 0);
	glutInitWindowSize (vis->screen.pixels[0], vis->screen.pixels[1]);

	glutCreateWindow ("NanoShaper");

	// glDisable (GL_DEPTH_TEST);
	glEnable (GL_DEPTH_TEST);
	glDisable (GL_BLEND);
	glShadeModel (GL_FLAT);
	glDisable (GL_DITHER);
	glDisable (GL_LIGHTING);

	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}
#endif


void Surface::Projection (Vis *vis)
{

	double temp;


	vis->screen.dim[0] = vis->perspective.ortho[0] / vis->perspective.zoom;
	vis->screen.dim[1] = vis->perspective.ortho[1] / vis->perspective.zoom;

	vis->viewpoint.sin_longitude = sinf(vis->perspective.longitude * DEG_TO_RAD);
	vis->viewpoint.cos_longitude = cosf(vis->perspective.longitude * DEG_TO_RAD);

	vis->viewpoint.sin_latitude = sinf(vis->perspective.latitude * DEG_TO_RAD);
	vis->viewpoint.cos_latitude = cosf(vis->perspective.latitude * DEG_TO_RAD);

	temp = vis->perspective.radius * vis->viewpoint.cos_latitude;

	vis->viewpoint.pos[0] = temp * vis->viewpoint.sin_longitude;
	vis->viewpoint.pos[1] = vis->perspective.radius * vis->viewpoint.sin_latitude;
	vis->viewpoint.pos[2] = temp * vis->viewpoint.cos_longitude;

	for (int l = 0; l < 3; l++)
		vis->viewpoint.pos[l] += vis->scene_center[l];

	vis->viewpoint.screen_dist = vis->perspective.radius / 2.0;

	temp = vis->viewpoint.screen_dist / vis->perspective.radius;

	for (int l = 0; l < 3; l++)
		vis->screen.ctr[l] = vis->viewpoint.pos[l] + temp*(vis->scene_center[l] - vis->viewpoint.pos[l]);

	Rotate (vis->screen.dim[0], 0.0, 0.0,
			vis->viewpoint.sin_longitude, vis->viewpoint.cos_longitude,
			vis->viewpoint.sin_latitude, vis->viewpoint.cos_latitude,
			&vis->screen.dir[0][0], &vis->screen.dir[0][1], &vis->screen.dir[0][2]);

	Rotate (0.0, vis->screen.dim[1], 0.0,
			vis->viewpoint.sin_longitude, vis->viewpoint.cos_longitude,
			vis->viewpoint.sin_latitude, vis->viewpoint.cos_latitude,
			&vis->screen.dir[1][0], &vis->screen.dir[1][1], &vis->screen.dir[1][2]);

	for (int l = 0; l < 3; l++)
    {
		vis->screen.vtx[l] = temp * (0.0 - vis->viewpoint.pos[l]) - vis->screen.dir[0][l] - vis->screen.dir[1][l];

		vis->screen.dir[0][l] *= (2.0 / (double)vis->screen.pixels[0]);
		vis->screen.dir[1][l] *= (2.0 / (double)vis->screen.pixels[1]);
    }

	#ifndef NO_OPENGL
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity ();

	// gluOrtho2D (0, vis->screen.pixels[0], 0, vis->screen.pixels[1]);
	glFrustum (-vis->screen.dim[0], vis->screen.dim[0],
			   -vis->screen.dim[1], vis->screen.dim[1],
			   vis->viewpoint.screen_dist, 1.0e+30F);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	if ((int)((fabs(vis->perspective.latitude) + 90.0) / 180.0) % 2 == 0)
	{
		gluLookAt (vis->viewpoint.pos[0], vis->viewpoint.pos[1], vis->viewpoint.pos[2],
				   vis->screen.ctr[0], vis->screen.ctr[1], vis->screen.ctr[2], 0.0, 1.0, 0.0);
	}
	else
	{
		gluLookAt (vis->viewpoint.pos[0], vis->viewpoint.pos[1], vis->viewpoint.pos[2],
				   vis->screen.ctr[0], vis->screen.ctr[1], vis->screen.ctr[2], 0.0, -1.0, 0.0);
	}

	glClearColor (vis->screen.col[0], vis->screen.col[1], vis->screen.col[2], 0.0);
	#endif
}


void Surface::Project (double p1[], double p2[], Vis *vis)
{
	double x1[4], x2[4];
	double temp;


	for (int l = 0; l < 4; l++)
		x1[l] = p1[l] - vis->viewpoint.pos[l];

	temp = vis->viewpoint.cos_longitude * x1[2] + vis->viewpoint.sin_longitude * x1[0];

	x2[0] = vis->viewpoint.cos_longitude * x1[0] - vis->viewpoint.sin_longitude * x1[2];
	x2[1] = vis->viewpoint.cos_latitude  * x1[1] - vis->viewpoint.sin_latitude * temp;
	x2[2] = vis->viewpoint.cos_latitude  * temp  + vis->viewpoint.sin_latitude * x1[1];

	p2[2] = -x2[2];
	temp = vis->viewpoint.screen_dist / p2[2];

	p2[0] = temp * x2[0];
	p2[1] = temp * x2[1];
}


#ifndef NO_OPENGL

void Surface::visualizeString (double r, double g, double b, int x, int y, char *string, void *font)
{
	glColor3f (r, g, b);
	glWindowPos2i (x, y);

	for (int i = 0; i < (int)strlen(string); i++)
	{
		glutBitmapCharacter (font, string[i]);
	}
	glEnd ();
}


void Surface::visualizeTriangulation (void)
{
	glColor3d (0., 0., 0.);

	int numTriangles = (int)(triList.size() / 3.);

	for (int n = 0; n < numTriangles; n++)
	{
		glBegin (GL_LINE_LOOP);

		#if !defined(USE_OPTIMIZED_VERTICES_BUFFERING)
		glVertex3dv (vertList[ triList[n*3+0] ]);
		glVertex3dv (vertList[ triList[n*3+1] ]);
		glVertex3dv (vertList[ triList[n*3+2] ]);
		#else
		glVertex3dv (&vertList[ triList[n*3+0]*3 ]);
		glVertex3dv (&vertList[ triList[n*3+1]*3 ]);
		glVertex3dv (&vertList[ triList[n*3+2]*3 ]);
		#endif

		glEnd ();
	}
}


void Surface::initVisualization (int argc, char *argv[], Vis *vis)
{
	vis->screen.pixels[0] = 512;
	vis->screen.pixels[1] = 512;

	vis->screen.col[0] = 1.0;
	vis->screen.col[1] = 1.0;
	vis->screen.col[2] = 1.0;

	vis->perspective.longitude = 45.0;
	vis->perspective.latitude  = 45.0;

	vis->perspective.zoom = 1.0;

	double box_size_x = delphi->xmax - delphi->xmin;
	double box_size_y = delphi->ymax - delphi->ymin;
	double box_size_z = delphi->zmax - delphi->zmin;
	double max_box_size = fmax(box_size_x, fmax(box_size_y, box_size_z));

	vis->perspective.ortho[0] = 0.5 * max_box_size;
	vis->perspective.ortho[1] = 0.5 * max_box_size;

	vis->perspective.radius = 4.0 * max_box_size;

	vis->scene_center[0] = 0.5 * (delphi->xmin + delphi->xmax);
	vis->scene_center[1] = 0.5 * (delphi->ymin + delphi->ymax);
	vis->scene_center[2] = 0.5 * (delphi->zmin + delphi->zmax);


	glutInit (&argc, argv);

	openWindow (vis);

	Projection (vis);

	initMenu (vis);
}

#endif // NO_OPENGL


void Surface::Display (void)
{
	#ifndef NO_OPENGL

	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	visualizeTriangulation ();

	glutSwapBuffers ();

	#endif // NO_OPENGL


	#ifndef NO_OPENGL

	++vis.editing.pause;

	#endif // NO_OPENGL
}


#ifndef NO_OPENGL

void Surface::rotateViewpoint (double dx, double dy, Vis *vis)
{
	vis->perspective.longitude += dx / DEG_TO_RAD;
	vis->perspective.latitude  += dy / DEG_TO_RAD;

	Projection (vis);
}

#endif


void Surface::endVisualization (Vis *vis)
{
	;
}


#ifndef NO_OPENGL


void Surface::initMenu (Vis *vis)
{
	vis->editing.sub_menu[1] = glutCreateMenu(changeView);
	glutAddMenuEntry ("Rotate", ROTATE);
	glutAddMenuEntry ("Zoom", ZOOM);
	glutAddMenuEntry ("Distance", DISTANCE);

	vis->editing.menu = glutCreateMenu (processMenuEvents);
	glutAddSubMenu ("View", vis->editing.sub_menu[1]);
	glutAddMenuEntry ("Quit", QUIT);

	glutAttachMenu (GLUT_RIGHT_BUTTON);
}


void processMenuEvents (int option)
{
	vis.menu.option = option;

	if (option == QUIT)
	{
		glutDestroyMenu (vis.editing.menu);

		exit(0);
	}
	vis.editing.pause = 0;
}


void changeView (int option)
{
	vis.editing.option = option;
	vis.editing.pause = 0;
}


void Surface::motionFunction (int x, int y)
{
	double mouse_dx, mouse_dy;


	y = vis.screen.pixels[1] - y - 1;

	mouse_dx = (double)(x - vis.mouse.x[0]) / vis.screen.pixels[0];
	mouse_dy = (double)(y - vis.mouse.x[1]) / vis.screen.pixels[1];

	if (vis.editing.option == ROTATE)
	{
		rotateViewpoint (-1.0e+4 * mouse_dx, -1.0e+4 * mouse_dy, &vis);
	}
	else if (vis.editing.option == ZOOM)
	{
		vis.perspective.zoom *= (1.0 + mouse_dy);

		Projection (&vis);
	}
	else if (vis.editing.option == DISTANCE)
	{
		double box_size_x = delphi->xmax - delphi->xmin;
		double box_size_y = delphi->ymax - delphi->ymin;
		double box_size_z = delphi->zmax - delphi->zmin;
		double max_box_size = fmax(box_size_x, fmax(box_size_y, box_size_z));

		vis.perspective.radius *= (1.0 + mouse_dy);
		vis.perspective.radius = fmax(vis.perspective.radius, 2.0 * max_box_size);

		Projection (&vis);
	}

	vis.mouse.x[0] = x;
	vis.mouse.x[1] = y;

	vis.editing.pause = 0;
}


void Surface::mouseFunction (int button, int state, int x, int y)
{
	if (button != GLUT_LEFT_BUTTON || state != GLUT_DOWN)
	{
		return;
	}

	y = vis.screen.pixels[1] - y - 1;

	vis.mouse.x[0] = x;
	vis.mouse.x[1] = y;

	vis.editing.pause = 0;
}


void Surface::keybordFunction (unsigned char key, int x, int y)
{
   ;
}


void Surface::Reshape (GLsizei w, GLsizei h)
{
	vis.perspective.ortho[0] *= (double)w / (double)vis.screen.pixels[0];
	vis.perspective.ortho[1] *= (double)h / (double)vis.screen.pixels[1];

	vis.screen.pixels[0] = w;
	vis.screen.pixels[1] = h;

	glViewport (0, 0, w, h);

	Projection (&vis);

	vis.editing.pause = 0;
}


#endif // NO_OPENGL


#endif // USE_VIS_TOOLS
