
//---------------------------------------------------------
/**    @file		Surface.cpp
*     @brief	Surface.cpp is the Surface CLASS
*															*/
//---------------------------------------------------------

#include "Surface.h"

void Surface::init()
{
	delphi = NULL;
	panel = 0;	
	totalSurfaceArea=0;
	totalVolume=0;	
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
	last_rows_ind = -1;
	last_cols_ind = -1;

	intersectionsMatrixAlongX = NULL;	
	intersectionsMatrixAlongY = NULL;
	intersectionsMatrixAlongZ = NULL;

	normalsMatrixAlongX = NULL;	
	normalsMatrixAlongY = NULL;
	normalsMatrixAlongZ = NULL;

	bgp_type = NULL;
	gridMultiMap = NULL;
	verticesInsidenessMap = NULL;
	sternLayer = -1;
	isRCbased = true;
	randDisplacement = RAND_DISPLACEMENT;
	gridLoad = NULL;
	useLoadBalancing = true;
	totalLoad = 0;
	vertexAtomsMap=NULL;
	vertexAtomsMapFlag=false;
	computeNormals = false;
	saveMSMS = false;
	providesAnalyticalNormals = false;
	activeCubes = NULL;
	MAX_ATOMS_MULTI_GRID = 100;
}

void Surface::init(ConfigFile* cf)
{
	bool projBGP = cf->read<bool>( "Project_boundary_grid_points", false );	
	bool accTri = cf->read<bool>( "Accurate_Triangulation", false );
	bool checkDuplicatedVertices = cf->read<bool>( "Check_duplicated_vertices", true );
	bool wellShaped = cf->read<bool>( "Keep_Water_Shaped_Cavities", false );
	double probeRadius = cf->read<double>( "Probe_Radius", 1.4 );
	bool lb = cf->read<bool>( "Load_Balancing",true);
	bool vaFlag = cf->read<bool>( "Vertex_Atom_Info",false);		
	bool computeNormals = cf->read<bool>( "Compute_Vertex_Normals",false);		
	bool saveMSMS = cf->read<bool>( "Save_Mesh_MSMS_Format",false);		
	double sternLayer = cf->read<double>( "Stern_layer", -1. );
	MAX_ATOMS_MULTI_GRID = cf->read<int>( "Max_Atoms_Multi_Grid", 100 );

	setProjBGP(projBGP);
	setTriangulationFlag(accTri);
	setCheckDuplicatedVertices(checkDuplicatedVertices);
	setKeepWellShapedCavities(wellShaped);
	setProbeRadius(probeRadius);
	setLoadBalancing(lb);
	setVertexAtomsMap(vaFlag);
	setComputeNormals(computeNormals);
	setSaveMSMS(saveMSMS);		
	
	// if >0 enable stern layer, else disabled by default
	if (sternLayer>0)
		setSternLayer(sternLayer);
}

Surface::Surface()
{
	init();
}


Surface::Surface(ConfigFile* cf)
{
	init();
	init(cf);
}

void Surface::clear()
{
	if (triList.size()>0)
	{
		vector<int*>::iterator it;
		for (it = triList.begin();it!=triList.end();it++)
			deleteVector<int>((*it));
	}

	if (vertList.size()>0)
	{
		vector<double*>::iterator it;
		for (it = vertList.begin();it!=vertList.end();it++)
			deleteVector<double>((*it));
	}

	if (normalsList.size()>0)
	{
		vector<double*>::iterator it;
		for (it = normalsList.begin();it!=normalsList.end();it++)
			deleteVector<double>((*it));
	}

	if (intersectionsMatrixAlongX!=NULL)
	{
		for (int i=0;i<delphi->nx;i++)
			for (int j=0;j<delphi->ny;j++)
				for (int k=0;k<delphi->nz;k++)
					if (intersectionsMatrixAlongX->at(i,j,k)!=-1)
						intersectionsMatrixAlongX->erase(i,j,k);
					
		delete intersectionsMatrixAlongX;
	}

	if (intersectionsMatrixAlongY!=NULL)
	{
		for (int i=0;i<delphi->nx;i++)
			for (int j=0;j<delphi->ny;j++)
				for (int k=0;k<delphi->nz;k++)
					if (intersectionsMatrixAlongY->at(i,j,k)!=-1)
						intersectionsMatrixAlongY->erase(i,j,k);					
		delete intersectionsMatrixAlongY;
	}

	if (intersectionsMatrixAlongZ!=NULL)
	{
		for (int i=0;i<delphi->nx;i++)
			for (int j=0;j<delphi->ny;j++)
				for (int k=0;k<delphi->nz;k++)
					if (intersectionsMatrixAlongZ->at(i,j,k)!=-1)
						intersectionsMatrixAlongZ->erase(i,j,k);
					
		delete intersectionsMatrixAlongZ;
	}

	if (normalsMatrixAlongX!=NULL)
	{
		for (int i=0;i<delphi->nx;i++)
			for (int j=0;j<delphi->ny;j++)
				for (int k=0;k<delphi->nz;k++)
					if (normalsMatrixAlongX->at(i,j,k)!=-1)
						normalsMatrixAlongX->erase(i,j,k);
					
		delete normalsMatrixAlongX;
	}

	if (normalsMatrixAlongY!=NULL)
	{
		for (int i=0;i<delphi->nx;i++)
			for (int j=0;j<delphi->ny;j++)
				for (int k=0;k<delphi->nz;k++)
					if (normalsMatrixAlongY->at(i,j,k)!=-1)
						normalsMatrixAlongY->erase(i,j,k);
					
		delete normalsMatrixAlongY;
	}

	if (normalsMatrixAlongZ!=NULL)
	{
		for (int i=0;i<delphi->nx;i++)
			for (int j=0;j<delphi->ny;j++)
				for (int k=0;k<delphi->nz;k++)
					if (normalsMatrixAlongZ->at(i,j,k)!=-1)
						normalsMatrixAlongZ->erase(i,j,k);
					
		delete normalsMatrixAlongZ;
	}
	
	if (verticesInsidenessMap!=NULL)
		deleteMatrix3D<bool>(delphi->nx,delphi->ny,verticesInsidenessMap);
		//deleteVector<bool>(verticesInsidenessMap);

	if (scalarField!=NULL)
		deleteMatrix3D<double>(delphi->nx,delphi->ny,scalarField);

	if (gridLoad!=NULL)
		deleteVector<int>(gridLoad);

	if (vertexAtomsMap!=NULL)
		deleteVector<int>(vertexAtomsMap);
}

Surface::~Surface()
{	
	clear();
}

bool Surface::build()
{
	cout << endl << WARN << "Build surface not supported!";
	return false;
}


bool Surface::save(char* fileName)
{
	cout << endl << WARN << "Save surface not supported!";
	return false;
}

bool Surface::load(char* fileName)
{
	cout << endl << WARN << "Load surface not supported";
	return false;
}

void Surface::printSummary()
{
	cout << endl << WARN << "Print summary not supported!";
}


void Surface::allocIntersectionsMatrices()
{
	if (intersectionsMatrixAlongX!=NULL)
		delete intersectionsMatrixAlongX;
				
	intersectionsMatrixAlongX = new Octree<int>(1024,-1);

	if (intersectionsMatrixAlongY!=NULL)
		delete intersectionsMatrixAlongY;

	intersectionsMatrixAlongY = new Octree<int>(1024,-1);

	if (intersectionsMatrixAlongZ!=NULL)
		delete intersectionsMatrixAlongZ;
	
	intersectionsMatrixAlongZ = new Octree<int>(1024,-1);
}


void Surface::allocNormalsMatrices()
{
	if (normalsMatrixAlongX!=NULL)
		delete normalsMatrixAlongX;
					
	normalsMatrixAlongX = new Octree<int>(1024,-1);

	if (normalsMatrixAlongY!=NULL)
		delete normalsMatrixAlongY;

	normalsMatrixAlongY = new Octree<int>(1024,-1);

	if (normalsMatrixAlongZ!=NULL)
		delete normalsMatrixAlongZ;
	
	normalsMatrixAlongZ = new Octree<int>(1024,-1);
}


/** build epsmap and compute boundary grid points by intersection
and projection routines. The minimal number of projections is performed.
The ray tracing part is parallelized with boost threading if enabled*/
bool Surface::getSurf(bool fillCav,double vol)
{
	double volPanel[3]={0,0,0};	
	double duration = 0;
	double** volPanelTHD;
	int* numIntersections;

	// per thread independent data structures. this avoid to manage
	// a mutex on the octree, thus it is faster and simpler		
	vector<coordVec>* buffersIntersections;
	vector<coordVec>* buffersNormals;

	Timer chrono_ray;

	fillCavitiesFlag = fillCav;

	if (delphi == NULL)
	{
		cout << endl << WARN << "Cannot get surface without DelPhi environment!";
		return false;
	}

	int num_thd=1;
	int old_num_thd = num_thd;

	// if Ray casting based perform parallel ray casting
	// if not assumes that grids colouring took place into the build phase
	if (isRCbased)
	{
		#ifdef ENABLE_BOOST_THREADS

		cout << endl << INFO << "Use load balancing: " << useLoadBalancing;

		if (num_cores<=0)
		{
			cout << endl << INFO << "Detected " << (int)boost::thread::hardware_concurrency() << " logical cores";
			num_thd = MIN(MIN(MIN(delphi->nx,delphi->ny),delphi->nz), (int)boost::thread::hardware_concurrency());
			cout << endl << INFO << "Setting " << num_thd << " threads";
		}
		else
		{
			cout << endl << INFO << "User selected num threads " << num_cores;
			num_thd = num_cores;
		}

		old_num_thd = num_thd;

		#endif
		
		// this is always needed. Vertices are stored here in any case
		allocIntersectionsMatrices();

		// this is used only if requested and are analytically computed
		// if not analytically compute octree will not be used
		if (computeNormals && providesAnalyticalNormals)
			allocNormalsMatrices();
		
		if (accurateTriangulation)
		{	
			if (verticesInsidenessMap!=NULL)
				deleteMatrix3D<bool>(last_nx,last_ny,verticesInsidenessMap);
				//deleteVector<bool>(verticesInsidenessMap);
			
			verticesInsidenessMap = allocateMatrix3D<bool>(delphi->nx,delphi->ny,delphi->nz); 
			//verticesInsidenessMap = allocateVector<bool>(delphi->nx*delphi->ny*delphi->nz);
			
			/*
			size_t tot = delphi->nx*delphi->ny*delphi->nz;
			for (size_t i=0;i<tot;i++)
				verticesInsidenessMap[i]=true;				
			*/

			for (int i=0;i<delphi->nx;i++)
				for (int j=0;j<delphi->ny;j++)
					for (int k=0;k<delphi->nz;k++)
						verticesInsidenessMap[i][j][k]=true;				
		}
		else
			verticesInsidenessMap = NULL;

		last_nx = delphi->nx;
		last_ny = delphi->ny;
		last_nz = delphi->nz;
		
		int na,nb;
				
		if (!delphi->getMultiDiel())
			cout << endl << INFO << "Inside id value is " << inside;
		
		#ifdef DEBUG_SURFACE
			// force 1 thread for debug
			num_thd = 1;
		#endif

		volPanelTHD = allocateMatrix2D<double>(num_thd,3);
		numIntersections = allocateVector<int>(num_thd);
		buffersIntersections = new vector<coordVec>[num_thd];		
		buffersNormals = new vector<coordVec>[num_thd];

		for (int i=0;i<num_thd;i++)
		{
			buffersIntersections[i].reserve(1000);
			buffersNormals[i].reserve(1000);
		}	
			
		int numint = 0;
		
		for (int l=0;l<num_thd;l++)
			numIntersections[l]=0;

		chrono_ray.start();		
		
		// Phase 1, ray trace from each coordinate plane if necessary
		for (panel=0;panel<3;panel++)
		{		
			if (!accurateTriangulation || isAvailableScalarField)
			{
				// if cavities and espilon map is not necessary then skip this step at all
				if (!delphi->buildStatus && !delphi->buildEpsMap)
				{
					// skip
					break;
				}

				// if only cavity detection is required, only the first panel is ray cast
				else if (delphi->buildStatus && !delphi->buildEpsMap && panel!=0) 
				{
					// skip panel 1,2
					continue;
				}
			}

			preProcessPanel();
			cout.flush();
			cout << endl << INFO << "Ray-tracing panel " << panel << "...";
			volPanel[panel]=0;
		
			// left y,z panel
			if (panel==0)
			{
				na = delphi->ny;
				nb = delphi->nz;
			}
			// lower y,x panel
			else if (panel==1)
			{
				na = delphi->ny;
				nb = delphi->nx;
			}
			// back z,x panel
			else
			{
				na = delphi->nz;
				nb = delphi->nx;		
			}
				
			// trivial split
			int chunk = na/num_thd;
			int rem = na%num_thd;

			#ifdef ENABLE_BOOST_THREADS
				boost::thread_group thdGroup; 
			#endif
			
			// setup split
			int j=0;
			int start = 0;
			int stop  = 0;
			
			if (!useLoadBalancing)
			{
				for (j=0;j<num_thd;j++)
				{
					if (j==0)
					{
						start = 0;
						stop = chunk;
					}
					else
					{
						start = stop;
						stop = start+chunk;
					}
										
					if (j<rem)
						stop++;

					packet pack;
					pack.first = &buffersIntersections[j];
					pack.second = &buffersNormals[j];
						
					#ifdef ENABLE_BOOST_THREADS
						thdGroup.create_thread(boost::bind(&Surface::intersector,this,volPanelTHD[j],nb,start,stop,1,&(numIntersections[j]),pack));
					#else
						intersector(volPanelTHD[j],nb,start,stop,1,&(numIntersections[j]),pack);
					#endif					
				}
			}
			else
			{
				int jump = num_thd;
				for (int i=0;i<num_thd;i++)
				{				
					packet pack;
					pack.first = &buffersIntersections[i];
					pack.second = &buffersNormals[i];

					#ifdef ENABLE_BOOST_THREADS
						thdGroup.create_thread(boost::bind(&Surface::intersector,this,volPanelTHD[i],nb,i,na,jump,&(numIntersections[i]),pack));
					#else
						intersector(volPanelTHD[i],nb,i,na,jump,&(numIntersections[i]),pack);
					#endif				
				}
			}
			
			// end setup

			#ifdef ENABLE_BOOST_THREADS
				// join
				thdGroup.join_all();
				// reduce
				for (int j=0;j<num_thd;j++)
				{	
					volPanel[panel]+=volPanelTHD[j][panel];
					numint+=numIntersections[j];
				}				
			#else			
				volPanel[panel]=volPanelTHD[0][panel];
				numint = numIntersections[0];
			#endif
			
			cout << "ok!";		
		}       

		double ray_time = chrono_ray.stop();
		cout << endl << INFO << "Ray-tracing computation time.. " << ray_time << " [s]"; 	
		
		// assuming squared grid for this stat
		if (accurateTriangulation && !isAvailableScalarField)
			printf("\n%sApproximated %d rays (%.5f %%)",INFO,numint,numint/(6*(float)delphi->nx*delphi->ny)*100);
		else
			printf("\n%sApproximated %d rays (%.5f %%)",INFO,numint,numint/(3*(float)delphi->nx*delphi->ny)*100);
	}

	Timer chrono_cav;
	chrono_cav.start();

	// before bgp identification and vertices storage, cavities are filled if requested
	if (fillCav)
	{
		cout << endl << INFO << "Performing cavity detection and conditional filling...";
		cout.flush();
		int cav = getCavities();
		cout << "ok!";
		cout << endl << INFO << "Detected " << cav << " cavitiy[ies]";		
		cout.flush();
		fillCavities(vol);
		if (wellShaped)
		{
			cout << endl << INFO << "Performing cavities shape filtering..";
			cout.flush();
			filterCavities();
		}
		cout << endl << INFO << "Recovering cavities atoms....";
		cout.flush();
		getCavitiesAtoms();
		cout << "ok!";
	}

	duration = chrono_cav.stop();
	cout << endl << INFO << "Cavity detection time is " << duration <<  " [s]";

	vector<int*> bgp;
	vector<int> bgp_type_temp;

	// after cavity detection assemble Octree, such that vertices that are switched off, can be removed now
	if (isRCbased)
	{
		cout << endl << INFO << "Assembling octrees..";
		
		// clean if required
		
		vector<int*>::iterator triIt;
		for (triIt = triList.begin();triIt!=triList.end();triIt++)
			free(*triIt);		
		triList.clear();

		vector<double*>::iterator ite;
		for (ite = vertList.begin();ite!=vertList.end();ite++)		
			free(*ite);		
		vertList.clear();

		for (ite = normalsList.begin();ite!=normalsList.end();ite++)		
			free(*ite);		
		normalsList.clear();

		// end clean

		int vert_index = 0;
		int norm_index = 0;
			
		// assembling final octrees and vertices/normals matrix
		for (int l=0;l<num_thd;l++)
		{			
			vector<coordVec>* v1 = &buffersIntersections[l];
			vector<coordVec>* v2 = &buffersNormals[l];

			int ix_,iy_,iz_,dir;
			double *intersec,*normal;
		
			vector<coordVec>::iterator it1;		
			vector<coordVec>::iterator it2;		
			
			if (computeNormals && providesAnalyticalNormals)
				it2 = v2->begin();

			// unload intersections/normals buffer
			for (it1=v1->begin();it1!=v1->end();it1++)
			{
				const coordVec& cv = (*it1);
				ix_ = cv.ix;
				iz_ = cv.iz;
				iy_ = cv.iy;
				dir = cv.dir;
				intersec = cv.vec;

				if (computeNormals && providesAnalyticalNormals)
				{
					const coordVec& cv2 = (*it2);
					normal = cv2.vec;
				}

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
				// it is absolutely negligible

				if (dir==X_DIR)
				{
					if (verticesInsidenessMap[ix_][iy_][iz_]==verticesInsidenessMap[ix_+1][iy_][iz_])
					//if (read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_,delphi->nx,delphi->ny,delphi->nz)==read3DVector<bool>(verticesInsidenessMap,ix_+1,iy_,iz_,delphi->nx,delphi->ny,delphi->nz))
					{
						free(intersec);
						if (computeNormals && providesAnalyticalNormals)
						{
							free(normal);
							it2++;
						}
						continue;
					}
				}
				else if (dir==Y_DIR)				
				{
					if (verticesInsidenessMap[ix_][iy_][iz_]==verticesInsidenessMap[ix_][iy_+1][iz_])
					//if (read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_,delphi->nx,delphi->ny,delphi->nz)==read3DVector<bool>(verticesInsidenessMap,ix_,iy_+1,iz_,delphi->nx,delphi->ny,delphi->nz))
					{
						free(intersec);
						if (computeNormals && providesAnalyticalNormals)
						{
							free(normal);
							it2++;
						}
						continue;
					}				
				}
				else
				{
					if (verticesInsidenessMap[ix_][iy_][iz_]==verticesInsidenessMap[ix_][iy_][iz_+1])
					//if (read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_,delphi->nx,delphi->ny,delphi->nz)==read3DVector<bool>(verticesInsidenessMap,ix_,iy_,iz_+1,delphi->nx,delphi->ny,delphi->nz))
					{
						free(intersec);
						if (computeNormals && providesAnalyticalNormals)
						{
							free(normal);
							it2++;
						}
						continue;
					}								
				}
								
				// ok the vertex is not dangling
				// add to octrees and remove duplicate vertices on the same edge, if any
				// TODO now duplicated vertices should be absent, check.
				if (dir==X_DIR)
				{					
					if (intersectionsMatrixAlongX->at(ix_,iy_,iz_)!=-1)
					{
						free(vertList[vertList.size()-1]);
						vertList.erase(vertList.end()-1);						
					}
					else
					{						
						intersectionsMatrixAlongX->set(ix_,iy_,iz_,vert_index);
						vert_index++;
					}

					if (computeNormals && providesAnalyticalNormals)
					{
						if (normalsMatrixAlongX->at(ix_,iy_,iz_)!=-1)
						{
							free(normalsList[normalsList.size()-1]);
							normalsList.erase(normalsList.end()-1);						
						}
						else
						{						
							normalsMatrixAlongX->set(ix_,iy_,iz_,norm_index);
							norm_index++;
						}
					}
				}	
				else if (dir==Y_DIR)				
				{
					if (intersectionsMatrixAlongY->at(ix_,iy_,iz_)!=-1)
					{
						free(vertList[vertList.size()-1]);
						vertList.erase(vertList.end()-1);						
					}
					else
					{
						intersectionsMatrixAlongY->set(ix_,iy_,iz_,vert_index);
						vert_index++;
					}

					if (computeNormals && providesAnalyticalNormals)
					{
						if (normalsMatrixAlongY->at(ix_,iy_,iz_)!=-1)
						{
							free(normalsList[normalsList.size()-1]);
							normalsList.erase(normalsList.end()-1);						
						}
						else
						{
							normalsMatrixAlongY->set(ix_,iy_,iz_,norm_index);
							norm_index++;
						}
					}

				}
				else if (dir==Z_DIR)
				{
					if (intersectionsMatrixAlongZ->at(ix_,iy_,iz_)!=-1)
					{
						free(vertList[vertList.size()-1]);
						vertList.erase(vertList.end()-1);
					}
					else
					{
						intersectionsMatrixAlongZ->set(ix_,iy_,iz_,vert_index);
						vert_index++;
					}

					if (computeNormals && providesAnalyticalNormals)
					{
						if (normalsMatrixAlongZ->at(ix_,iy_,iz_)!=-1)
						{
							free(normalsList[normalsList.size()-1]);
							normalsList.erase(normalsList.end()-1);
						}
						else
						{
							normalsMatrixAlongZ->set(ix_,iy_,iz_,norm_index);
							norm_index++;
						}
					}

				}
				else
				{
					cout << endl << ERR << "Non existing direction during octree assembling!";
					exit(-1);
				}

				vertList.push_back(intersec);
				if (computeNormals && providesAnalyticalNormals)
				{
					normalsList.push_back(normal);
					it2++;
				}
			}
		}

		cout << "ok!";
		
		num_thd = old_num_thd;

		delete[] buffersIntersections;
		delete[] buffersNormals;
		deleteMatrix2D<double>(num_thd,volPanelTHD);	
		deleteVector<int>(numIntersections);
	}

	// if only buildEpsMap is active, this correction due to cavity detection
	// is not necessary
	
	if (delphi->buildEpsMap && delphi->buildStatus)
	{
		cout << endl << INFO << "Writing idebmap...";
		int NX = delphi->nx;
		int NY = delphi->ny;
		int NZ = delphi->nz;
		//short*** status = delphi->status;
		short* status = delphi->status;
		bool* idebmap = delphi->idebmap;
		// write idebmap
		for (int k=0;k<NZ;k++)
			for (int j=0;j<NY;j++)				
					for (int i=0;i<NX;i++)
					{
						// inside
						//if (status[i][j][k]==STATUS_POINT_INSIDE) IDEBMAP(i,j,k,nx,ny)=false;
						//if ((STATUSMAP(i,j,k,NX,NY))==STATUS_POINT_INSIDE) 
						if (read3DVector<short>(status,i,j,k,NX,NY,NZ)==STATUS_POINT_INSIDE)
						{
							//IDEBMAP(i,j,k,NX,NY)=false;
							write3DVector<bool>(idebmap,false,i,j,k,NX,NY,NZ);
						}
						// outside
						else 
						{
							//IDEBMAP(i,j,k,NX,NY)=true;
							write3DVector<bool>(idebmap,true,i,j,k,NX,NY,NZ);
						}
					}
		cout << "ok!";
	}

	// Stern Layer
	// build it if enabled and if a epsilon map flag was enabled
	// TODO if (sternLayer>0 && (delphi->buildEpsMap or getDelphiBinding()))
	if (sternLayer>0 && delphi->buildEpsMap)
	{
		cout << endl << INFO << "Computing Stern Layer...";	
		buildSternLayer();
		cout << "ok!";
	}

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	// rest bgp count
	delphi->nbgp = 0;

	// custom clean up
	postRayCasting();

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
			cout << endl << INFO << "Applying multiple dielectric correction..";
			buildAtomsMap();
			applyMultidielectric();
			cout << "ok!";
		}
	}

	int internal_bgps =0;
	int external_bgps =0;

	// Phase 2. boundary grid points are located and projected if requested
	if (projBGP)
	{
		preBoundaryProjection();
		cout << endl << INFO << "Detecting boundary grid points...";
		for (int iz=0;iz<NZ;iz++)	
		{
			for (int iy=0;iy<NY;iy++)
			{
				for (int ix=0;ix<NX;ix++)
				{				
					//int kost = delphi->EPSMAP(ix,iy,iz,0,NX,NY,NZ);
					const int konst = read4DVector<int>(delphi->epsmap,ix,iy,iz,0,NX,NY,NZ,3);

					if ((ix-1<0) || (iy-1<0) || (iz-1<0))
						continue;

					/*if (konst != delphi->EPSMAP(ix,iy,iz,1,NX,NY,NZ) ||
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
						int* temp = allocateVector<int>(3);
						temp[0] = ix;
						temp[1] = iy;
						temp[2] = iz;
						bgp.push_back(temp);

						// if only one point is external than the bgp is external (surface bgp)
						// else it is internal bgp (different dielectrics)							
						/*if (delphi->EPSMAP(ix,iy,iz,0,NX,NY,NZ) == 0 ||
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
			delphi->scsarea = allocateVector<double>(delphi->nbgp);
			delphi->ibgp	= allocateVector<int>(3*delphi->nbgp);
		}

		if (delphi->getDelphiBinding() && bgp.size()>=((unsigned)delphi->maxbgp))
		{
			cout << endl << ERR << "Number of bgp is " << bgp.size() << " and the maximum allowed is " << delphi->maxbgp;
			cout << endl << ERR << "Please increase ibmx in DelPhi and recompile";
			exit(-1);
		}

		if (bgp_type!=NULL)
			deleteVector<int>(bgp_type);

		bgp_type = allocateVector<int>(delphi->nbgp);

		// freeze vector
		vector<int*>::iterator it;
		int i=0;
		for (it=bgp.begin();it!=bgp.end();it++)
		{
			//delphi->ibgp[i] = *it;
			int* vv = *it;
			delphi->ibgp[3*i]=vv[0];
			delphi->ibgp[3*i+1]=vv[1];
			delphi->ibgp[3*i+2]=vv[2];
			deleteVector<int>(vv);
			i++;
		}

		for (unsigned int l=0;l<bgp_type_temp.size();l++)
			bgp_type[l]=bgp_type_temp[l];

		cout << "ok!";
		cout << endl << INFO << "Detected " << delphi->nbgp << " boundary grid points (bgp)";
		cout << endl << INFO << "Detected " << external_bgps << " external bgps, and " << internal_bgps << " internal bgps";
		cout << endl << INFO << "Scaling bgps...";

		#ifdef ENABLE_BOOST_THREADS
			boost::thread_group thdGroup; 
		#endif

		#ifdef ENABLE_BOOST_THREADS
			if (num_cores<=0)
			{
				// the overhead of thread creation is quite low, so it is
				// convenient to dispatch as thread as possible in order both
				// to saturate all cores and to minimize their unbalancing 
				num_thd = MIN(64,MAX(1,delphi->nbgp));
				//cout << endl << INFO << "Automatically dispatching " << num_thd << " threads";
			}
			else
			{
				//cout << endl << INFO << "User selected num threads " << num_cores;
				num_thd = num_cores;
			}
		#else
			num_thd = 1;
		#endif
		
		int chunk = delphi->nbgp/num_thd;
		int rem = delphi->nbgp%num_thd;

		// setup split
		int j=0,start,stop;
		for (j=0;j<num_thd;j++)
		{
			if (j==0)
			{
				start = 0;
				stop = chunk;
			}
			else
			{
				start = stop;
				stop = start+chunk;
			}

			if (j<rem)
				stop++;

			//printf("\n %d %d",start,stop);

			#ifdef ENABLE_BOOST_THREADS
				thdGroup.create_thread(boost::bind(&Surface::projector,this,start,stop));
			#else
				projector(start,stop);
			#endif
		}				
		// end setup
		
		#ifdef ENABLE_BOOST_THREADS
			thdGroup.join_all();
		#endif

		cout << "ok!";		

		deleteVector<int>(bgp_type);
		bgp_type = NULL;

		int* atsurf = delphi->atsurf;

		// setting atsurf = NULL it is a way to avoid this computation
		if (delphi->getDelphiBinding() && atsurf!=NULL)
		{
			cout << endl << INFO << "Linking boundary grid points to nearest atom";
			// if multi-diel is enabled the atoms map is already available
			// if not we have to build it
			if (!delphi->getMultiDiel())
				buildAtomsMap();
						
			double* l_scspos = delphi->scspos;
						
			for (int i=0;i<delphi->nbgp;i++)
			{				
				int nearestAtom=-1;
				double* v = &(l_scspos[3*i]);
				vdwAccessible(v,nearestAtom);
				if (nearestAtom==-1) 
					cout << endl << WARN << "Cannot detect nearest atom for bgp index "<<i;	
				atsurf[i]=nearestAtom;
			}
		}		
	}                   		
	
	// if multi-dielectric or stern layer was computed then dispose atoms map
	// in any case if delphi is bound and bgps have been computed we have to drop it
	if (delphi->getMultiDiel() || (delphi->getDelphiBinding() && delphi->atsurf!=NULL && projBGP))
		disposeAtomsMap();

	// assumes that the volume is computed in a ray casting based way
	if (isRCbased)
	{
		//printf("\n%f %f %f\n",volPanel[0],volPanel[1],volPanel[2]);
		volPanel[0]*=delphi->A;		
		volPanel[1]*=delphi->A;
		volPanel[2]*=delphi->A;
		
		// compute the total volume only considering one plane
		if (delphi->buildStatus && !delphi->buildEpsMap && panel!=0) 
			totalVolume = volPanel[0];
		// compute the total volume as the average of the 3 volumes
		else
			totalVolume=(volPanel[0]+volPanel[1]+volPanel[2])/3;
	}

	return true;
}

void Surface::buildAtomsMap()
{
	// Build a 3D accelaration grid for the atoms	
	double rmax = 0;
	// get the biggest atom 
	for (int i=0;i<delphi->numAtoms;i++)
		rmax = MAX((delphi->atoms[i]->radius),rmax);
	
	gside = 2*rmax;	
	gscale = 1./gside;
	ggrid = 1;
	
	while(1)
	{
		ggrid = (unsigned int)floor(gscale*1.3*delphi->rmaxdim);
		
		if (ggrid<=100)
		{		
			if (ggrid<=0)
			{
				gside=gside/2.;
				gscale = 1./gside;
			}
			else
				break;
		}
		else
		{
			gside=gside*2;
			gscale = 1./gside;
		}
	}	

	gxmin = delphi->baricenter[0]-(ggrid-1)/(2*gscale);
	gymin = delphi->baricenter[1]-(ggrid-1)/(2*gscale);
	gzmin = delphi->baricenter[2]-(ggrid-1)/(2*gscale);

	if (gridMultiMap!=NULL)
		deleteVector<int>(gridMultiMap);

	gridMultiMap = allocateVector<int>(MAX_ATOMS_MULTI_GRID*ggrid*ggrid*ggrid);
	for (unsigned int k=0;k<ggrid;k++)	
		for (unsigned int j=0;j<ggrid;j++)
			for (unsigned int i=0;i<ggrid;i++)
				// init the number of atoms mapped in that place 
				//GRID_MULTI_MAP(i,j,k,0,ggrid,ggrid,ggrid)=0;
				write4DVector<int>(gridMultiMap,0,i,j,k,0,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

	// fill 3D acceleration grid with atoms
	for (int i=0;i<delphi->numAtoms;i++)
	{
		int ix = (int)rintp((delphi->atoms[i]->pos[0]-gxmin)/gside);
		int iy = (int)rintp((delphi->atoms[i]->pos[1]-gymin)/gside);
		int iz = (int)rintp((delphi->atoms[i]->pos[2]-gzmin)/gside);			
		//int max_ind = (size_t)GRID_MULTI_MAP(ix,iy,iz,0,ggrid,ggrid,ggrid);
		int max_ind = (size_t)read4DVector<int>(gridMultiMap,ix,iy,iz,0,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

		max_ind++;
		if (max_ind>=MAX_ATOMS_MULTI_GRID)
		{
			cout << endl << ERR << "Increase Max_Atoms_Multi_Grid. Current value is " << MAX_ATOMS_MULTI_GRID;
			exit(-1);
		}

		//GRID_MULTI_MAP(ix,iy,iz,0,ggrid,ggrid,ggrid)=max_ind;
		write4DVector<int>(gridMultiMap,max_ind,ix,iy,iz,0,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

		//printf("\n %d %d %d atom %d",ix,iy,iz,i);
		//GRID_MULTI_MAP(ix,iy,iz,max_ind,ggrid,ggrid,ggrid) = i;
		write4DVector<int>(gridMultiMap,i,ix,iy,iz,max_ind,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);
	}
}

void Surface::disposeAtomsMap()
{
	// remove acceleration grid
	deleteVector<int>(gridMultiMap);
	gridMultiMap = NULL;
}

void Surface::applyMultidielectric()
{
	// for each internal point in the eps map gets the nearest atom
	// and fix the epsilon map accordingly.

	if (gridMultiMap==NULL)
	{
		cout << endl << WARN << "Cannot apply multi-dielectric correction without atoms map";
		return;
	}
		
	// for each internal point get nearest atom and change eps value accordingly
	for (int k=0;k<delphi->nz;k++)
		for (int j=0;j<delphi->ny;j++)
			for (int i=0;i<delphi->nx;i++)			
			{				
				// for each cube face fix the epsmap value
				//int value = delphi->EPSMAP(i,j,k,0,delphi->nx,delphi->ny,delphi->nz);
				int value = read4DVector<int>(delphi->epsmap,i,j,k,0,delphi->nx,delphi->ny,delphi->nz,3);

				if (value==inside)
					swap2multi(gxmin,gymin,gzmin,gside,ggrid,gridMultiMap,i,j,k,0);

				//value = delphi->EPSMAP(i,j,k,1,(delphi->nx),(delphi->ny),(delphi->nz));
				value = read4DVector<int>(delphi->epsmap,i,j,k,1,delphi->nx,delphi->ny,delphi->nz,3);

				if (value==inside)
					swap2multi(gxmin,gymin,gzmin,gside,ggrid,gridMultiMap,i,j,k,1);
				
				//value = delphi->EPSMAP(i,j,k,2,(delphi->nx),(delphi->ny),(delphi->nz));
				value = read4DVector<int>(delphi->epsmap,i,j,k,2,delphi->nx,delphi->ny,delphi->nz,3);
				
				if (value==inside)
					swap2multi(gxmin,gymin,gzmin,gside,ggrid,gridMultiMap,i,j,k,2);
				
			}			
}


void Surface::swap2multi(double gxmin,double gymin,double gzmin,double gside,unsigned int ggrid,int* gridMultiMap,int i,int j,int k,int l)
{
	// generate point position
	double pos[3];

	pos[0] = delphi->x[i];
	pos[1] = delphi->y[j];
	pos[2] = delphi->z[k];

	int ori = i, orj=j,ork=k,orl=l;
	// position the point on the corrseponding cube side. Pos now is the position of the midpoint
	pos[l]+=delphi->hside;

	// move from midpoint of delphi grid to auxiliary atoms grid
	int ix = (int)rintp((pos[0]-gxmin)/gside);
	int iy = (int)rintp((pos[1]-gymin)/gside);
	int iz = (int)rintp((pos[2]-gzmin)/gside);

	double minDist=INFINITY;
	int winner = -1;
	
	// get the nearest atom to set the dielectric constant
	// the dielectric constant is mapped according to the additively weighted voronoi diagram
	// that is the signed distance from the point p is ||p-c||^2-r^2 where c is the center
	// of the atom and r is the radius. The minimum signed distance wins.
	for (k=0;k<SHIFT_MAP;k++)
	{
		unsigned int cx = ix+shift_map[k][0];
		unsigned int cy = iy+shift_map[k][1];
		unsigned int cz = iz+shift_map[k][2];

		if (cx>(ggrid-1) || cy>(ggrid-1) || cz>(ggrid-1) || cx<0 || cy<0 || cz<0)
			continue;

		//int max_ind = GRID_MULTI_MAP(cx,cy,cz,0,ggrid,ggrid,ggrid);
		int max_ind = read4DVector<int>(gridMultiMap,cx,cy,cz,0,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

		for (int j=0;j<max_ind;j++)
		{
			//int atom_index = GRID_MULTI_MAP(cx,cy,cz,j+1,ggrid,ggrid,ggrid);
			int atom_index = read4DVector<int>(gridMultiMap,cx,cy,cz,j+1,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

			double signed_dist = 0;
			DIST2(signed_dist,delphi->atoms[atom_index]->pos,pos)
			double rad = delphi->atoms[atom_index]->radius;
			signed_dist-=(rad*rad);
			if (signed_dist<minDist)
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
	//delphi->EPSMAP(ori,orj,ork,orl,delphi->nx,delphi->ny,delphi->nz) = delphi->atoms[winner]->dielectric;
	write4DVector<int>(delphi->epsmap,delphi->atoms[winner]->dielectric,ori,orj,ork,orl,delphi->nx,delphi->ny,delphi->nz,3);
	return;	
}

double Surface::getVolume()
{
	return totalVolume;
}

double Surface::getArea()
{
	return totalSurfaceArea;
}

int Surface::getCavities(int idStart)
{
	// search cavities in status vector
	bool cavities=true;
	int id = idStart;
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	short* status = delphi->status;

	if (status==NULL)
	{
		cout << endl << ERR << "Cannot do cavity detection without a status map";
		cout << endl << REMARK << "Please set Build_status_map = true";
		cout << endl;
		exit(-1);	
	}

	// free memory 
	if (delphi->cavitiesVec !=NULL)
	{
		vector< vector < int* >* >::iterator it;
		for (it=delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++)
		{
			vector < int* >* inner = (*it);
			vector < int* >::iterator it2;
			for (it2=inner->begin();it2!=inner->end();it2++)
				free((*it2));
			delete inner;
		}
		delete delphi->cavitiesVec;
	}
	
	// clean memory
	delphi->cavitiesSize.clear();	
	delphi->cavitiesFlag.clear();

	// allocate empty vector
	delphi->cavitiesVec = new vector< vector < int* >* >();
	delphi->cavitiesVec->reserve(50);
	
	int times = 0;
	int startZ = 0;
	int endZ = NZ-1;
	int stepZ = 1;
	int sig = 1;
	
	while(1)
	{
		int i,j,k;
		cavities = false;
		
		if (times%2==0)
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
			stepZ=-1;
			sig=-1;
		}
		
		times++;
		
		// get starting point
		for (k=startZ;(sig*k)<=endZ;k=k+stepZ)
		{
			for (j=0;j<NY;j++)
			{
				for (i=0;i<NX;i++)
				{
					//if (STATUSMAP(i,j,k,NX,NY)==STATUS_POINT_TEMPORARY_OUT)
					if (read3DVector<short>(status,i,j,k,NX,NY,NZ)==STATUS_POINT_TEMPORARY_OUT)
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
		// if idStart = STATUS_POINT_TEMPORARY_OUT the first pass is from STATUS_POINT_TEMPORARY_OUT -> STATUS_POINT_OUT this is not cavity, this is only real outside
		// from that moment on, if there are still cavities they are all marked with STATUS_POINT_TEMPORARY_OUT
		// and one by one they are detected and associated to a cavity index
		//floodFill(i,j,k,STATUS_POINT_TEMPORARY_OUT,id);	
		//floodFill4(i,j,k,STATUS_POINT_TEMPORARY_OUT,id);	
		floodFill2(i,j,k,STATUS_POINT_TEMPORARY_OUT,id);			
	}

	int numCavities = id-STATUS_POINT_OUT;	
	delphi->cavitiesSize.resize(numCavities);
	delphi->cavitiesFlag.resize(numCavities);
	for (int i=0;i<numCavities;i++)
		delphi->cavitiesFlag[i]=false;

	return numCavities;
}

void Surface::getCavitiesAtoms()
{	
	delphi->initCav2Atoms();
	vector< set<int>* >& cav2atoms = delphi->cav2atoms;
	vector < vector <int* >* >::iterator it;	
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	short* status = delphi->status;

	if (delphi->buildStatus==false)
	{
		cout << endl << WARN << "Cannot get cavity atoms without status map";
		return;
	}

	if (delphi->cavitiesVec->size()==0)			
		return;
	
	buildAtomsMap();
	int i=0;

	//cout << endl << cav2atoms.size() << " " << delphi->cavitiesVec->size();

	for (it=delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++)
	{	
		vector <int* >::iterator it2;
		vector<int*>* vec = (*it);
		
		// Get bgps for each cavity. And obtain the nearest atom
		// in order to avoid using epsmap a slightly changed notion of bgp is employed
		for (it2=vec->begin();it2!=vec->end();it2++)
		{
			int* v = (*it2);
			unsigned int ix=v[0],iy=v[1],iz=v[2];

			//cout << endl << ix << " " << iy << " " << iz;
			
			//int kost = STATUSMAP(ix,iy,iz,NX,NY);
			const int kost = read3DVector<short>(status,ix,iy,iz,NX,NY,NZ);

			//(cout << endl << ix << " " << iy << " " << iz << " " << NX).flush();
			//(cout << endl << "prima if").flush();

			// is it pseudo-bgp? if yes continue else skip.
			/*if (kost != STATUSMAP((ix+1),iy,iz,NX,NY) ||
				kost != STATUSMAP((ix-1),iy,iz,NX,NY) ||
				kost != STATUSMAP(ix,(iy+1),iz,NX,NY) ||
				kost != STATUSMAP(ix,(iy-1),iz,NX,NY) ||
				kost != STATUSMAP(ix,iy,(iz+1),NX,NY) ||
				kost != STATUSMAP(ix,iy,(iz-1),NX,NY))*/

			if (kost != read3DVector<short>(status,ix+1,iy,iz,NX,NY,NZ) ||
				kost != read3DVector<short>(status,ix-1,iy,iz,NX,NY,NZ) ||
				kost != read3DVector<short>(status,ix,iy+1,iz,NX,NY,NZ) ||
				kost != read3DVector<short>(status,ix,iy-1,iz,NX,NY,NZ) ||
				kost != read3DVector<short>(status,ix,iy,iz+1,NX,NY,NZ) ||
				kost != read3DVector<short>(status,ix,iy,iz-1,NX,NY,NZ))
			{
				//(cout << endl << "dopo if").flush();
			}
			else
			{
				//(cout << endl << "dopo if").flush();
				continue;
			}

			double gridPoint[3];

			gridPoint[0] = delphi->x[ix];
			gridPoint[1] = delphi->y[iy];
			gridPoint[2] = delphi->z[iz];

			// move from grid point of delphi grid to auxiliary atoms grid
			ix = (int)rintp((gridPoint[0]-gxmin)/gside);
			iy = (int)rintp((gridPoint[1]-gymin)/gside);
			iz = (int)rintp((gridPoint[2]-gzmin)/gside);
		
			double minDist = INFINITY;
			unsigned int first=-1;
			
			// get the nearest atom
			for (int k=0;k<SHIFT_MAP;k++)
			{
				unsigned int cx = ix+shift_map[k][0];
				unsigned int cy = iy+shift_map[k][1];
				unsigned int cz = iz+shift_map[k][2];

				if (cx>(ggrid-1) || cy>(ggrid-1) || cz>(ggrid-1) || cx<0 || cy<0 || cz<0)
					continue;

				//int max_ind = GRID_MULTI_MAP(cx,cy,cz,0,ggrid,ggrid,ggrid);
				int max_ind = read4DVector<int>(gridMultiMap,cx,cy,cz,0,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

				for (int j=0;j<max_ind;j++)
				{					
					//int atom_index = GRID_MULTI_MAP(cx,cy,cz,j+1,ggrid,ggrid,ggrid);					
					int atom_index = read4DVector<int>(gridMultiMap,cx,cy,cz,j+1,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

					double signed_dist = 0;					
					DIST2(signed_dist,delphi->atoms[atom_index]->pos,gridPoint)										
					double radius2 = delphi->atoms[atom_index]->radius2;			
					signed_dist-=radius2;
					
					if (signed_dist<minDist)
					{
						minDist = signed_dist;					
						first = atom_index;
					}
				}	
			}
			if (first==-1)
				cout << endl << WARN << "No nearest atom in cavity/pocket!";
						
			cav2atoms[i]->insert(first);			
		}
		i++;
	}
	disposeAtomsMap();
}
void Surface::fillCavities(double vol,bool silent)
{
	if (vol<0)
	{
		cout << endl << WARN << "Cannot fill with a negative volume. Setting " << DEFAULT_VOLUME;
		vol = DEFAULT_VOLUME;
	}
	// analyze each cavity and fill if required
	vector < vector <int* >* >::iterator it;
	int i=0;
	double cubeVol = delphi->side*delphi->side*delphi->side;

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	if (!silent)
	{
		cout << endl << INFO << "Threshold volume is " << vol;
		cout << endl << INFO << "Tot num cavities is " << delphi->cavitiesVec->size();
	}
	for (it=delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++)
	{
		if (!silent)
			cout << endl << INFO << "Cavity " << i ;
		
		double cavVol = (*it)->size()*cubeVol;
		if (!silent)
			printf("\tvol is %.4lf [A^3] \t",cavVol);

		if (internals!=NULL)
		{
			*internals << endl << "cav " << cavVol;
		}

		delphi->cavitiesSize[i]=cavVol;
		delphi->cavitiesFlag[i]=false;

		if (cavVol<=vol)
		{
			delphi->cavitiesFlag[i]=true;
			vector <int* >::iterator it2;
			vector<int*>* vec = (*it);
			// filling eps map
			for (it2=vec->begin();it2!=vec->end();it2++)
			{
				int* v = (*it2);
				// that grid point is filled

				// apply correction only if epsmap is used
				if (delphi->buildEpsMap)
				{
					/*
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),0,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),1,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),2,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]-1),(v[1]),(v[2]),0,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]-1),(v[2]),1,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]-1),2,NX,NY,NZ)=inside;
					*/
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
					
					verticesInsidenessMap[v[0]][v[1]][v[2]] = false;
					verticesInsidenessMap[v[0]+1][v[1]][v[2]] = false;
					verticesInsidenessMap[v[0]][v[1]+1][v[2]] = false;
					verticesInsidenessMap[v[0]+1][v[1]+1][v[2]] = false;

					verticesInsidenessMap[v[0]][v[1]][v[2]+1] = false;
					verticesInsidenessMap[v[0]+1][v[1]][v[2]+1] = false;
					verticesInsidenessMap[v[0]][v[1]+1][v[2]+1] = false;
					verticesInsidenessMap[v[0]+1][v[1]+1][v[2]+1] = false;
					/*
					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2],delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2],delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2],delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2],delphi->nx,delphi->ny,delphi->nz);

					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2]+1,delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2]+1,delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2]+1,delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2]+1,delphi->nx,delphi->ny,delphi->nz);
					*/
				}
				// update the scalar filed consistently with cavity detection
				if (accurateTriangulation && isAvailableScalarField)
				{
					scalarField[v[0]][v[1]][v[2]] = +INFINITY;
					scalarField[v[0]+1][v[1]][v[2]] = +INFINITY;
					scalarField[v[0]][v[1]+1][v[2]] = +INFINITY;
					scalarField[v[0]+1][v[1]+1][v[2]] = +INFINITY;

					scalarField[v[0]][v[1]][v[2]+1] = +INFINITY;
					scalarField[v[0]+1][v[1]][v[2]+1] = +INFINITY;
					scalarField[v[0]][v[1]+1][v[2]+1] = +INFINITY;
					scalarField[v[0]+1][v[1]+1][v[2]+1] = +INFINITY;
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

// conservative filter
/*
void Surface::filterCavities()
{	
	int ngrid = (int)rintp(2.*probe_radius/delphi->side);
	if ((ngrid%2)==0)
		ngrid++;

	int k = (ngrid-1)/2;
	//cout << " Probe radius is " << probe_radius << " [A], grid points " << k;
	int cavityId;

	double dist2ref = probe_radius*probe_radius;
	vector < vector <int* >* >::iterator it;
	int i=0;

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	// analyze each cavity
	for (it=delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++,i++)
	{
		cout << endl << INFO << "Cavity " << i ;
		if (delphi->cavitiesFlag[i])
		{
			cout << " already filled";
			continue;
		}
		else
			cout << " Under analysis...";

		vector<int*>* currentCavity = (*it);
		vector <int* >::iterator it2;		

		// mark each cavity point as temporary STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK code
		for (it2=currentCavity->begin();it2!=currentCavity->end();it2++)
		{
			int* v = *it2;
			cavityId = delphi->status[v[0]][v[1]][v[2]];
			delphi->status[v[0]][v[1]][v[2]] = STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK;
		}

		// for each cavity point 'march' around it in order to understand if the 
		// water molecule can fit in the cavity
		
		// a cavity point is filled if a probe sphere cannot stay into the cavity.

		bool keepCavity = false;

		for (it2=currentCavity->begin();it2!=currentCavity->end();it2++)
		{
			int* v = *it2;			
			int start_ix = v[0]-k;
			int start_iy = v[1]-k;
			int start_iz = v[2]-k;

			int end_ix = v[0]+k;
			int end_iy = v[1]+k;
			int end_iz = v[2]+k;
		
			bool oneProbeFits = true;

			for (int ix = start_ix;ix<end_ix;ix++)
			{
				for (int iy = start_iy;iy<end_iy;iy++)
				{
					for (int iz = start_iz;iz<end_iz;iz++)
					{												
						double a = (delphi->x[ix]-delphi->x[v[0]])*(delphi->x[ix]-delphi->x[v[0]]);
						double b = (delphi->y[iy]-delphi->y[v[1]])*(delphi->y[iy]-delphi->y[v[1]]);
						double c = (delphi->z[iz]-delphi->z[v[2]])*(delphi->z[iz]-delphi->z[v[2]]);

						double dist2 = a+b+c;
						if (dist2>dist2ref)
							continue;

						// there is no space we went outside the cavity, go on to the next grid point
						if (delphi->status[ix][iy][iz]!=STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK)
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

		for (it2=currentCavity->begin();it2!=currentCavity->end();it2++)
		{	
			int* v = *it2;								
			// clean cavity
			if (!keepCavity)
			{		

				if (delphi->buildEpsMap)
				{
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),0,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),1,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),2,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]-1),(v[1]),(v[2]),0,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]-1),(v[2]),1,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]-1),2,NX,NY,NZ)=inside;
					
					delphi->IDEBMAP((v[0]),(v[1]),(v[2]),NX,NY) = false;
				}

				if (accurateTriangulation && !isAvailableScalarField)
				{
					verticesInsidenessMap[v[0]][v[1]][v[2]] = false;
					verticesInsidenessMap[v[0]+1][v[1]][v[2]] = false;
					verticesInsidenessMap[v[0]][v[1]+1][v[2]] = false;
					verticesInsidenessMap[v[0]+1][v[1]+1][v[2]] = false;

					verticesInsidenessMap[v[0]][v[1]][v[2]+1] = false;
					verticesInsidenessMap[v[0]+1][v[1]][v[2]+1] = false;
					verticesInsidenessMap[v[0]][v[1]+1][v[2]+1] = false;
					verticesInsidenessMap[v[0]+1][v[1]+1][v[2]+1] = false;
				}
				
				// update the scalar filed consistently with cavity shape detection
				if (accurateTriangulation && isAvailableScalarField)
				{
					scalarField[v[0]][v[1]][v[2]] = +INFINITY;
					scalarField[v[0]+1][v[1]][v[2]] = +INFINITY;
					scalarField[v[0]][v[1]+1][v[2]] = +INFINITY;
					scalarField[v[0]+1][v[1]+1][v[2]] = +INFINITY;

					scalarField[v[0]][v[1]][v[2]+1] = +INFINITY;
					scalarField[v[0]+1][v[1]][v[2]+1] = +INFINITY;
					scalarField[v[0]][v[1]+1][v[2]+1] = +INFINITY;
					scalarField[v[0]+1][v[1]+1][v[2]+1] = +INFINITY;
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
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	short* status = delphi->status;
	int i=0;

	// for each cavity if it is not filled revert to OUT
	// if it is filled revert to IN
	vector < vector <int* >* >::iterator it;

	if (delphi->cavitiesVec==NULL)
	{
		cout << endl << "WARN" << "Cannot convert cavities to out flags if cavities are absent!";
		return;
	}

	for (i=0,it=delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++,i++)
	{
		vector<int*>* currentCavity = (*it);
		vector<int*>::iterator it2;

		for (it2 = currentCavity->begin();it2!=currentCavity->end();it2++)
		{
			int* v = *it2;
			
			// if a confirmed cavity (it can be STATUS_POINT_INISDE if filtered) then switch to out
			// if it is a point that is a support of the cavity, switch it too
			/*
			if (STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY)>=STATUS_FIRST_CAV || STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY)<=STATUS_FIRST_SUPPORT_CAV)
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = STATUS_POINT_TEMPORARY_OUT;
			*/

			short val = read3DVector<short>(status,v[0],v[1],v[2],NX,NY,NZ);
			if (val >=STATUS_FIRST_CAV || val<=STATUS_FIRST_SUPPORT_CAV)
				write3DVector<short>(status,STATUS_POINT_TEMPORARY_OUT,v[0],v[1],v[2],NX,NY,NZ);
		}
	}
}

void Surface::filterCavities(bool modStatus)
{	
	//cout << " Probe radius to filter cavities is " << probe_radius << " [A]";
	int cavityId;	
	double dist2ref = probe_radius*probe_radius;

	vector < vector <int* >* >::iterator it;
	int i=0;
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	short* status = delphi->status;

	// Connolly like filter
	for (i=0,it=delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++,i++)
	{	
		if (delphi->cavitiesFlag[i])
			continue;

		vector<int*>* currentCavity = (*it);
		vector <int* >::iterator it2;

		cavityId = i+STATUS_FIRST_CAV;

		// mark each cavity point as temporary STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK code
		for (it2=currentCavity->begin();it2!=currentCavity->end();it2++)
		{
			int* v = *it2;
			//STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK;
			write3DVector<short>(status,STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK,v[0],v[1],v[2],NX,NY,NZ);
		}

		// for each cavity point 'march' around it in order to understand if the 
		// water molecule can fit in around the current cavity point
		for (it2=currentCavity->begin();it2!=currentCavity->end();it2++)
		{
			int* v = *it2;
			
			double downx = delphi->x[v[0]]-probe_radius;
			double downy = delphi->y[v[1]]+probe_radius;
			double downz = delphi->z[v[2]]-probe_radius;

			double upx = delphi->x[v[0]]+probe_radius;
			double upy = delphi->y[v[1]]-probe_radius;
			double upz = delphi->z[v[2]]+probe_radius;
	          
			// resolve which are the grid cubes that
			// are cut by the bounding box object. These
			// cubes see the cell bounding box    
			int ix_start = (int)rintp((downx-delphi->xmin)/delphi->side);
			int iy_start = (int)rintp((downy-delphi->ymin)/delphi->side);
			int iz_start = (int)rintp((downz-delphi->zmin)/delphi->side);
	     	     
			int ix_end = (int)rintp((upx-delphi->xmin)/delphi->side);
			int iy_end = (int)rintp((upy-delphi->ymin)/delphi->side);
			int iz_end = (int)rintp((upz-delphi->zmin)/delphi->side);

			if (ix_start<0) ix_start=0;
			if (iy_start<0) iy_start=0;
			if (iz_start<0) iz_start=0;

			if (ix_end<0) ix_end=0;
			if (iy_end<0) iy_end=0;
			if (iz_end<0) iz_end=0;

			if (ix_end>=(int)NX) ix_end=NX-1;
			if (iy_end>=(int)NY) iy_end=NY-1;	          
			if (iz_end>=(int)NZ) iz_end=NZ-1;	          

			if (ix_start>=(int)NX) ix_start=NX-1;
			if (iy_start>=(int)NY) iy_start=NY-1;	          
			if (iz_start>=(int)NZ) iz_start=NZ-1;	          
			
			int tot = 0;
			int count = 0;
			bool go_out = false;

			for (int ix = ix_start;ix<=ix_end;ix++)
			{
				for (int iy = iy_start;iy>=iy_end;iy--)
				{
					for (int iz = iz_start;iz<=iz_end;iz++)
					{	
						/*
						double a = (delphi->x[ix]-delphi->x[v[0]])*(delphi->x[ix]-delphi->x[v[0]]);
						double b = (delphi->y[iy]-delphi->y[v[1]])*(delphi->y[iy]-delphi->y[v[1]]);
						double c = (delphi->z[iz]-delphi->z[v[2]])*(delphi->z[iz]-delphi->z[v[2]]);
						*/

						double a = (delphi->x[ix]-delphi->x[v[0]]);
						a = a*a;
						double b = (delphi->y[iy]-delphi->y[v[1]]);
						b = b*b;
						double c = (delphi->z[iz]-delphi->z[v[2]]);
						c = c*c;

						double dist2 = a+b+c;
						if (dist2>dist2ref)
							continue;
									
						short val = read3DVector<short>(status,ix,iy,iz,NX,NY,NZ);
						// this is inside the probe and it is an out point, stop immediately, go to next cavity point						
						/*if (STATUSMAP(ix,iy,iz,NX,NY)!=STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK && 
							STATUSMAP(ix,iy,iz,NX,NY)!=STATUS_CAVITY_POINT_SHAPE_OK && 
							STATUSMAP(ix,iy,iz,NX,NY)!= -cavityId)*/

						if (val!=STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK && val!=STATUS_CAVITY_POINT_SHAPE_OK && val!= -cavityId)
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
				for (int ix = ix_start;ix<=ix_end;ix++)
					for (int iy = iy_start;iy>=iy_end;iy--)
						for (int iz = iz_start;iz<=iz_end;iz++)
						{	
							/*
							double a = (delphi->x[ix]-delphi->x[v[0]])*(delphi->x[ix]-delphi->x[v[0]]);
							double b = (delphi->y[iy]-delphi->y[v[1]])*(delphi->y[iy]-delphi->y[v[1]]);
							double c = (delphi->z[iz]-delphi->z[v[2]])*(delphi->z[iz]-delphi->z[v[2]]);
							*/
							double a = (delphi->x[ix]-delphi->x[v[0]]);
							a = a*a;
							double b = (delphi->y[iy]-delphi->y[v[1]]);
							b = b*b;
							double c = (delphi->z[iz]-delphi->z[v[2]]);
							c = c*c;

							double dist2 = a+b+c;
							if (dist2>dist2ref)
								continue;
							else
							{
								// don't overwrite support cavity points
								//if (STATUSMAP(ix,iy,iz,NX,NY)!= -cavityId)
								//	STATUSMAP(ix,iy,iz,NX,NY)=STATUS_CAVITY_POINT_SHAPE_OK;

								if (read3DVector<short>(status,ix,iy,iz,NX,NY,NZ)!=-cavityId)
									write3DVector<short>(status,STATUS_CAVITY_POINT_SHAPE_OK,ix,iy,iz,NX,NY,NZ);
							}
						}				
				// the point is marked as a support of the cavity
				//STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY)= -cavityId;
				write3DVector<short>(status,-cavityId,v[0],v[1],v[2],NX,NY,NZ);
			}			
		}

		// number of grid points removed
		int numRemoved = 0;

		// the parts of the cavity that remained unmarked are due to the
		// non possibility to fit on it the probe
		// thus these points are remarked as inside.
		for (it2=currentCavity->begin();it2!=currentCavity->end();it2++)
		{
			int* v = *it2;
			// here the probe does not fit thus put inside. this point is no more in cavity		
			//if (STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY)==STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK)
			if (read3DVector<short>(status,v[0],v[1],v[2],NX,NY,NZ)==STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK)
			{
				numRemoved ++;
			
				if (delphi->buildEpsMap)
				{
					/*
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),0,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),1,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]),2,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]-1),(v[1]),(v[2]),0,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]-1),(v[2]),1,NX,NY,NZ)=inside;
					delphi->EPSMAP((v[0]),(v[1]),(v[2]-1),2,NX,NY,NZ)=inside;
					
					delphi->IDEBMAP((v[0]),(v[1]),(v[2]),NX,NY) = false;
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
					
					verticesInsidenessMap[v[0]][v[1]][v[2]] = false;
					verticesInsidenessMap[v[0]+1][v[1]][v[2]] = false;
					verticesInsidenessMap[v[0]][v[1]+1][v[2]] = false;
					verticesInsidenessMap[v[0]+1][v[1]+1][v[2]] = false;

					verticesInsidenessMap[v[0]][v[1]][v[2]+1] = false;
					verticesInsidenessMap[v[0]+1][v[1]][v[2]+1] = false;
					verticesInsidenessMap[v[0]][v[1]+1][v[2]+1] = false;
					verticesInsidenessMap[v[0]+1][v[1]+1][v[2]+1] = false;
					
					/*
					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2],delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2],delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2],delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2],delphi->nx,delphi->ny,delphi->nz);

					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1],v[2]+1,delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1],v[2]+1,delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0],v[1]+1,v[2]+1,delphi->nx,delphi->ny,delphi->nz);
					write3DVector<bool>(verticesInsidenessMap,false,v[0]+1,v[1]+1,v[2]+1,delphi->nx,delphi->ny,delphi->nz);
					*/
				}
				
				// update the scalar filed consistently with cavity shape detection
				if (accurateTriangulation && isAvailableScalarField)
				{
					scalarField[v[0]][v[1]][v[2]] = +INFINITY;
					scalarField[v[0]+1][v[1]][v[2]] = +INFINITY;
					scalarField[v[0]][v[1]+1][v[2]] = +INFINITY;
					scalarField[v[0]+1][v[1]+1][v[2]] = +INFINITY;

					scalarField[v[0]][v[1]][v[2]+1] = +INFINITY;
					scalarField[v[0]+1][v[1]][v[2]+1] = +INFINITY;
					scalarField[v[0]][v[1]+1][v[2]+1] = +INFINITY;
					scalarField[v[0]+1][v[1]+1][v[2]+1] = +INFINITY;
				}
			}
			
			// if this point must be inside and status is allowed to be modified then set it inside
			/*
			if (STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY)==STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK && modStatus)
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY)=STATUS_POINT_INSIDE;*/

			if (read3DVector<short>(status,v[0],v[1],v[2],NX,NY,NZ)==STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK && modStatus)
				write3DVector<short>(status,STATUS_POINT_INSIDE,v[0],v[1],v[2],NX,NY,NZ);
			
			// else restore the orginal value
			/*
			else if (STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY)==STATUS_CAVITY_POINT_SHAPE_OK)			
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = cavityId;
			*/
			else if (read3DVector<short>(status,v[0],v[1],v[2],NX,NY,NZ)==STATUS_CAVITY_POINT_SHAPE_OK)			
				write3DVector<short>(status,cavityId,v[0],v[1],v[2],NX,NY,NZ);			
			
		} // end for each cavity point

		if (numRemoved && currentCavity->size()!=numRemoved)
		{
			// thus they can be get filtered again
			delphi->cavitiesFlag[i] = false;
			//cout << "Partially filled";
		}
		else if(numRemoved && currentCavity->size()==numRemoved)
		{
			delphi->cavitiesFlag[i] = true;
			//cout << "Completely filled";
		}
		//else
		//	cout << "Non filled";

	} // end for each cavity	
}


/** queue based implementation*/
void Surface::floodFill(int ix,int iy,int iz,int idold,int idnew)
{
	queue<pair<int,pair<int,int> > > myqueue;
	pair<int,pair<int,int> > triplet;
	int cix,ciy,ciz;
	short* status = delphi->status;
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	myqueue.push(pair<int,pair<int,int > >(ix,pair<int,int>(iy,iz)));

	while(myqueue.size()!=0)
	{	
		//cout << endl << myqueue.size();
		triplet = myqueue.front();
		myqueue.pop();
		
		cix = triplet.first;
		ciy = triplet.second.first;
		ciz = triplet.second.second;

		//cout << endl << myqueue.size();
		//cout << endl <<cix << " " << ciy << " " << ciz;

		//if (delphi->status[cix][ciy][ciz]==idold)
		//if (STATUSMAP(cix,ciy,ciz,NX,NY)==idold)
		if (read3DVector<short>(status,cix,ciy,ciz,NX,NY,NZ)==idold)
		{
			//delphi->status[cix][ciy][ciz]=idnew;
			//STATUSMAP(cix,ciy,ciz,NX,NY)=idnew;
			write3DVector<short>(status,idnew,cix,ciy,ciz,NX,NY,NZ);

			// if this is a cavity then update the cavity list
			if (idnew>=4)
			{
				// current cavity vector
				vector < int* >* vec;

				// first time this cavity is encountered
				if (delphi->cavitiesVec->size()<((unsigned)(idnew-4+1)))
				{
					//cout << endl << INFO << "Detected cavity " << idnew-4;
					vec = new vector<int*>();
					if (vec==NULL)
					{
						cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
						cout.flush();						
						exit(-1);
					}
					delphi->cavitiesVec->push_back(vec);
				}
				else
				{
					vec = delphi->cavitiesVec->at(idnew-4);
				}				
				
				int* v = allocateVector<int>(3);
				if (v==NULL)
				{
					cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
					cout.flush();
					exit(-1);
				}
				v[0] = cix;
				v[1] = ciy;
				v[2] = ciz;
				//printf("\n %d %d %d (%d,%d)",cix,ciy,ciz,idnew,idold);
				vec->push_back(v);
				
			}
		}
		else
			continue;
		
		//if (cix+1<delphi->nx && (delphi->status[cix+1][ciy][ciz])==idold)
		//if (cix+1<NX && (STATUSMAP((cix+1),(ciy),(ciz),NX,NY))==idold)
		if ((cix+1)<NX && read3DVector<short>(status,cix+1,ciy,ciz,NX,NY,NZ)==idold)
		{
			triplet.first = cix+1;
			triplet.second.first = ciy;
			triplet.second.second = ciz;
			myqueue.push(triplet);
		}
		//if (ciy+1<delphi->ny && (delphi->status[cix][ciy+1][ciz])==idold)
		//if (ciy+1<NY && (STATUSMAP(cix,(ciy+1),(ciz),NX,NY))==idold)
		if ((ciy+1)<NY && read3DVector<short>(status,cix,ciy+1,ciz,NX,NY,NZ)==idold)
		{
			triplet.first = cix;
			triplet.second.first = ciy+1;
			triplet.second.second = ciz;
			myqueue.push(triplet);		
		}		
		//if (ciz+1<delphi->nz && (delphi->status[cix][ciy][ciz+1])==idold)
		//if (ciz+1<NZ && (STATUSMAP(cix,ciy,(ciz+1),NX,NY))==idold)
		if ((ciz+1)<NZ && read3DVector<short>(status,cix,ciy,ciz+1,NX,NY,NZ)==idold)
		{
			triplet.first = cix;
			triplet.second.first = ciy;
			triplet.second.second = ciz+1;
			myqueue.push(triplet);				
		}

		//if (cix-1>=0 && (delphi->status[cix-1][ciy][ciz])==idold)
		//if (cix-1>=0 && (STATUSMAP((cix-1),(ciy),(ciz),NX,NY))==idold)
		if ((cix-1)>=0 && read3DVector<short>(status,cix-1,ciy,ciz,NX,NY,NZ)==idold)
		{
			triplet.first = cix-1;
			triplet.second.first = ciy;
			triplet.second.second = ciz;
			myqueue.push(triplet);
		}

		//if (ciy-1>=0 && (delphi->status[cix][ciy-1][ciz])==idold)
		//if (ciy-1>=0 && (STATUSMAP(cix,(ciy-1),ciz,NX,NY))==idold)
		if ((ciy-1)>=0 && read3DVector<short>(status,cix,ciy-1,ciz,NX,NY,NZ)==idold)
		{
			triplet.first = cix;
			triplet.second.first = ciy-1;
			triplet.second.second = ciz;
			myqueue.push(triplet);
		}

		//if (ciz-1>=0 && (delphi->status[cix][ciy][ciz-1])==idold)
		//if (ciz-1>=0 && (STATUSMAP(cix,ciy,(ciz-1),NX,NY))==idold)
		if ((ciz-1)>=0 && read3DVector<short>(status,cix,ciy,ciz-1,NX,NY,NZ)==idold)
		{
			triplet.first = cix;
			triplet.second.first = ciy;
			triplet.second.second = ciz-1;
			myqueue.push(triplet);
		}
	}
	return;
}

void Surface::floodFill4(int ix,int iy,int iz, int idold, int idnew)
{	
	int numMoves = 0;
	int maxMoves = 10;
	vector<pair< pair<int,int>, int> > seeds_down;
	vector<pair< pair<int,int>, int> > seeds_up;
	short* status = delphi->status;
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	
	seeds_down.push_back(pair< pair<int,int>, int>(pair<int,int>(ix,iy),iz));
	
	int ix_or = ix;
	int iy_or = iy;
	int iz_or = iz;

	cout << endl << INFO << "Z-percolation...";
	cout.flush();

	/////////////////////////////////
	// z-percolation to init threads
	/////////////////////////////////
	
	// lower percolation
	while(1)
	{
		//if ((iz-1)>=0 && STATUSMAP(ix,iy,(iz-1),NX,NY) == idold)
		if ((iz-1)>=0 && read3DVector<short>(status,ix,iy,iz-1,NX,NY,NZ)==idold)
		{
			seeds_down.push_back(pair< pair<int,int>, int>(pair<int,int>(ix,iy),iz-1));
			iz--;
			continue;
		}

		bool found = false;
		
		// +x test
		for (int i=0;i<maxMoves;i++)
		{
			// feasibility of the move
			//if (!((ix+i)<NX) || !(STATUSMAP((ix+i),iy,iz,NX,NY) == idold))
			if (!((ix+i)<NX) || !(read3DVector<short>(status,ix+1,iy,iz,NX,NY,NZ)==idold))
			{
				found = false;
				break;
			}
			
			//if ((ix+i)<NX && (iz-1)>=0 && STATUSMAP((ix+i),iy,(iz-1),NX,NY) == idold)
			if ((ix+i)<NX && (iz-1)>=0 && read3DVector<short>(status,ix+1,iy,iz-1,NX,NY,NZ)==idold)
			{
				seeds_down.push_back(pair< pair<int,int>, int>(pair<int,int>(ix+i,iy),iz-1));
				found = true;
				ix+=i;
				iz--;
				break;
			}
		}

		if (found)
			continue;

		// -x test
		for (int i=0;i<maxMoves;i++)
		{
			// feasibility of the move
			//if (!((ix-i)>=0) || !(STATUSMAP((ix-i),iy,iz,NX,NY) == idold))
			if (!((ix-i)>=0) || !(read3DVector<short>(status,ix-i,iy,iz,NX,NY,NZ)==idold))
			{
				found = false;
				break;
			}
			
			//if ((ix-i)>=0 && (iz-1)>=0 && STATUSMAP((ix-i),iy,(iz-1),NX,NY) == idold)
			if ((ix-i)>=0 && (iz-1)>=0 && read3DVector<short>(status,ix-1,iy,iz-1,NX,NY,NZ)==idold)
			{
				seeds_down.push_back(pair< pair<int,int>, int>(pair<int,int>(ix-i,iy),iz-1));
				found = true;
				ix-=i;
				iz--;
				break;
			}
		}


		// +y test
		for (int i=0;i<maxMoves;i++)
		{
			// feasibility of the move
			//if (!((iy+i)<NY) || !(STATUSMAP(ix,(iy+i),iz,NX,NY) == idold))
			if (!((iy+i)<NY) || !(read3DVector<short>(status,ix,iy+i,iz,NX,NY,NZ)==idold))
			{
				found = false;
				break;
			}
			
			//if ((iy+i)<NY && (iz-1)>=0 && STATUSMAP(ix,(iy+i),(iz-1),NX,NY) == idold)
			if ((iy+i)<NY && (iz-1)>=0 && read3DVector<short>(status,ix,iy+i,iz-1,NX,NY,NZ)==idold)
			{
				seeds_down.push_back(pair< pair<int,int>, int>(pair<int,int>(ix,iy+i),iz-1));
				found = true;
				iy+=i;
				iz--;
				break;
			}
		}

		if (found)
			continue;

		// -y test
		for (int i=0;i<maxMoves;i++)
		{
			// feasibility of the move
			//if (!((iy-i)>=0) || !(STATUSMAP(ix,(iy-i),iz,NX,NY) == idold))
			if (!((iy-i)>=0) || !(read3DVector<short>(status,ix,iy-i,iz,NX,NY,NZ)==idold))
			{
				found = false;
				break;
			}
			
			//if ((iy-i)>=0 && (iz-1)>=0 && STATUSMAP(ix,(iy-i),(iz-1),NX,NY) == idold)
			if ((iy-i)>=0 && (iz-1)>=0 && read3DVector<short>(status,ix,iy-i,iz-1,NX,NY,NZ)==idold)
			{
				seeds_down.push_back(pair< pair<int,int>, int>(pair<int,int>(ix,iy-i),iz-1));
				found = true;
				iy-=i;
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
	while(1)
	{
		//if ((iz+1)<NZ && STATUSMAP(ix,iy,(iz+1),NX,NY) == idold)
		if ((iz+1)<NZ && read3DVector<short>(status,ix,iy,iz+1,NX,NY,NZ)==idold)
		{
			seeds_up.push_back(pair< pair<int,int>, int>(pair<int,int>(ix,iy),iz+1));
			iz++;
			continue;
		}

		bool found = false;
		
		// +x test
		for (int i=0;i<maxMoves;i++)
		{
			// feasibility of the move
			//if (!((ix+i)<NX) || !(STATUSMAP((ix+i),iy,iz,NX,NY) == idold))
			if (!((ix+i)<NX) || !(read3DVector<short>(status,ix+i,iy,iz,NX,NY,NZ) == idold))
			{
				found = false;
				break;
			}
			
			//if ((ix+i)<NX && (iz+1)<NZ && STATUSMAP((ix+i),iy,(iz+1),NX,NY) == idold)
			if ((ix+i)<NX && (iz+1)<NZ && read3DVector<short>(status,ix+i,iy,iz+1,NX,NY,NZ) == idold)
			{
				seeds_up.push_back(pair< pair<int,int>, int>(pair<int,int>(ix+i,iy),iz+1));
				found = true;
				ix+=i;
				iz++;
				break;
			}
		}

		if (found)
			continue;

		// -x test
		for (int i=0;i<maxMoves;i++)
		{
			// feasibility of the move
			//if (!((ix-i)>=0) || !(STATUSMAP((ix-i),iy,iz,NX,NY) == idold))
			if (!((ix-i)>=0) || !(read3DVector<short>(status,ix-i,iy,iz,NX,NY,NZ) == idold))
			{
				found = false;
				break;
			}
			
			//if ((ix-i)>=0 && (iz+1)<NZ && STATUSMAP((ix-i),iy,(iz+1),NX,NY) == idold)
			if ((ix-i)>=0 && (iz+1)<NZ && read3DVector<short>(status,ix-i,iy,iz+1,NX,NY,NZ) == idold)
			{
				seeds_up.push_back(pair< pair<int,int>, int>(pair<int,int>(ix-i,iy),iz+1));
				found = true;
				ix-=i;
				iz++;
				break;
			}
		}


		// +y test
		for (int i=0;i<maxMoves;i++)
		{
			// feasibility of the move
			//if (!((iy+i)<NY) || !(STATUSMAP(ix,(iy+i),iz,NX,NY) == idold))
			if (!((iy+i)<NY) || !(read3DVector<short>(status,ix,iy+i,iz,NX,NY,NZ) == idold))
			{
				found = false;
				break;
			}
			
			//if ((iy+i)<NY && (iz+1)<NZ && STATUSMAP(ix,(iy+i),(iz+1),NX,NY) == idold)
			if ((iy+i)<NY && (iz+1)<NZ && read3DVector<short>(status,ix,iy+i,iz+1,NX,NY,NZ) == idold)
			{
				seeds_up.push_back(pair< pair<int,int>, int>(pair<int,int>(ix,iy+i),iz+1));
				found = true;
				iy+=i;
				iz++;
				break;
			}
		}

		if (found)
			continue;

		// -y test
		for (int i=0;i<maxMoves;i++)
		{
			// feasibility of the move
			//if (!((iy-i)>=0) || !(STATUSMAP(ix,(iy-i),iz,NX,NY) == idold))
			if (!((iy-i)>=0) || !(read3DVector<short>(status,ix,iy-i,iz,NX,NY,NZ) == idold))
			{
				found = false;
				break;
			}
			
			//if ((iy-i)>=0 && (iz+1)<NZ && STATUSMAP(ix,(iy-i),(iz+1),NX,NY) == idold)
			if ((iy-i)>=0 && (iz+1)<NZ && read3DVector<short>(status,ix,iy-i,iz+1,NX,NY,NZ) == idold)
			{
				seeds_up.push_back(pair< pair<int,int>, int>(pair<int,int>(ix,iy-i),iz+1));
				found = true;
				iy-=i;
				iz++;
				break;
			}
		}

		if (!found)
			break;
	}

	pair< pair<int,int>,int >* zz = allocateVector<pair< pair<int,int>,int > >(seeds_down.size()+seeds_up.size());
	int index = 0;

	int dim = (int)seeds_down.size()-1;
	for (int i=dim;i>=0;i--)	
		zz[index++]=seeds_down[i];

	dim = (int)seeds_up.size();
	for (int i=0;i<dim;i++)
		zz[index++]=seeds_up[i];
	
	//cout << endl << "ok!";
	int z_max,z_min;

	z_max = seeds_up[seeds_up.size()-1].second ;
	z_min = seeds_down[seeds_down.size()-1].second ;

	int num_thd = num_cores;
	
	// start thread partitioning
	// split load evenly between threads
	int num_z = z_max-z_min+1;

	cout << "ok!";
	cout.flush();	

	if (num_cores<=0 || num_cores>=num_z/2)
	{
		num_thd = 4;
		cout << endl << INFO << "Setting " << num_thd << " threads for floodfill";		
	}

	#ifdef ENABLE_BOOST_THREADS
		boost::thread_group thdGroup; 
	#else
		num_thd = 1;
	#endif

	queue< pair<pair<int,int>,int> >** upper_queues;
	queue< pair<pair<int,int>,int> >** lower_queues;

	upper_queues = allocateVector< queue< pair<pair<int,int>,int> >* > (num_thd-1);
	lower_queues = allocateVector< queue< pair<pair<int,int>,int> >* > (num_thd-1);

	int chunk = num_z/num_thd;
	int rem = num_z%num_thd;

	//cout << endl << "c " << chunk << " r " << rem;				
	// setup split
	int j=0;
	int start = 0;
	int stop  = 0;

	pair<int,int> old_new;
	pair<pair<int,int>,int> ind;
	
	old_new.first = idold;
	old_new.second = idnew;

	int upper_queue_index = 0;
	int lower_queue_index = 0;


	for (j=0;j<num_thd;j++)
	{
		if (j==0)
		{
			start = 0;
			stop = chunk;
		}
		else
		{
			start = stop;
			stop = start+chunk;
		}
									
		if (j<rem)
			stop++;
		
		pair<int,int> limits(start+z_min,stop+z_min);
		ind = zz[start];

		cout << endl << "--------------------------------------";
		cout << endl << start << " " << stop;
		cout << endl << start+z_min << " " << stop+z_min;		
		cout << endl << "Start at " << ind.first.first << " " << ind.first.second << " " << ind.second; 

		pair< queue< pair<pair<int,int>,int> > *, queue< pair<pair<int,int>,int > > * > queues;
		
		if (j==0)
		{
			queues.second = NULL;			
			upper_queues[upper_queue_index] = new queue< pair<pair<int,int>,int> >();
			queues.first = upper_queues[upper_queue_index];
			upper_queue_index++;
		}
		else if (j==num_thd-1)
		{
			queues.first = NULL;
			lower_queues[lower_queue_index] =  new queue< pair<pair<int,int>,int> >();
			queues.second = lower_queues[lower_queue_index];
			lower_queue_index++;			
		}
		else
		{
			lower_queues[lower_queue_index] =  new queue< pair<pair<int,int>,int> >();
			upper_queues[upper_queue_index] = new queue< pair<pair<int,int>,int> >();

			queues.first  = upper_queues[upper_queue_index];
			queues.second = lower_queues[lower_queue_index];

			upper_queue_index++;
			lower_queue_index++;			
		}

		//cout << endl << "queues " << queues.first << " " << queues.second;
		
		#ifdef ENABLE_BOOST_THREADS
			thdGroup.create_thread(boost::bind(&Surface::floodFill3,this,ind,limits,old_new,queues,(queue<pair<pair<int,int>,int> >*)NULL));
		#else
			floodFill3(ind,limits,old_new,queues,(queue<pair<pair<int,int>,int> >*)NULL);
		#endif	
	}
	
	#ifdef ENABLE_BOOST_THREADS
		thdGroup.join_all();
	#endif

	pair<int,int> limits(-1,-1);

	//cout << endl ;
	// run on the remaining queues
	for (int i=0;i<num_thd-1;i++)
		if (upper_queues[i]!=NULL && upper_queues[i]->size()!=0)
		{
			//cout << "*";
			ind = upper_queues[i]->front();
			floodFill3(ind,limits,old_new,
				pair< queue< pair<pair<int,int>,int> > *, queue< pair<pair<int,int>,int > > * >(NULL,NULL),
				upper_queues[i]);
		}

	// run on the remaining queues
	for (int i=0;i<num_thd-1;i++)	
		if (lower_queues[i]!=NULL && lower_queues[i]->size()!=0)
		{
			//cout << "*";
			ind = lower_queues[i]->front();
			floodFill3(ind,limits,old_new,
				pair< queue< pair<pair<int,int>,int> > *, queue< pair<pair<int,int>,int > > * >(NULL,NULL),
				lower_queues[i]);
		}

	// delete queues
	for (int i=0;i<num_thd-1;i++)	
		if (lower_queues[i]!=NULL)
			delete lower_queues[i];

	for (int i=0;i<num_thd-1;i++)	
		if (upper_queues[i]!=NULL)
			delete upper_queues[i];

	delete zz;
	delete lower_queues;
	delete upper_queues;

}
// partial flood fill routine to be run in parallel
void Surface::floodFill3(	pair<pair<int,int>,int > ind, 
							pair<int,int> z_limits, 
							pair<int,int> old_new, 
							pair< queue< pair<pair<int,int>,int> > *, queue< pair<pair<int,int>,int > > * > queues,
							queue< pair<pair<int,int>,int> > * in_queue)
{	

	int idold = old_new.first;
	int idnew = old_new.second;
	short* status = delphi->status;
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	int ix = ind.first.first;
	int iy = ind.first.second;
	int iz = ind.second;

	// already visited, stops
	//if (STATUSMAP(ix,iy,iz,NX,NY) == idnew)
	if (read3DVector<short>(status,ix,iy,iz,NX,NY,NZ)==idnew)
		return;

	pair<int,int> duplo;
	pair < pair<int,int>, int>  triplet;

	int iz_min = z_limits.first;
	int iz_max = z_limits.second; 

	queue<pair<pair<int,int>,int> >* z_queue;

	if (in_queue==NULL)
		z_queue = new queue<pair<pair<int,int>,int> >();
	else
		z_queue = in_queue;

	queue< pair< pair<int,int>,int > >* out_up_xy = queues.first;
	queue< pair< pair<int,int>,int > >* out_down_xy = queues.second;

	if (in_queue==NULL)
		z_queue->push(pair<pair<int,int>,int>(pair<int,int>(ix,iy),iz));

	while(z_queue->size()!=0)
	{
		triplet = z_queue->front();
		z_queue->pop();
				
		int iz = triplet.second;
		int ix = triplet.first.first;
		int iy = triplet.first.second;

		//if (STATUSMAP(ix,iy,iz,NX,NY) == idnew)	
		if (read3DVector<short>(status,ix,iy,iz,NX,NY,NZ)==idnew)
			continue;

		//cout << endl << iz;
		stack<pair<int,int> >  xy_stack;
		xy_stack.push(pair<int,int >(ix,iy));
		int y1;

		bool spanLeft=0,spanRight=0;
		
		while(xy_stack.size()!=0)
		{	
			duplo = xy_stack.top();
			xy_stack.pop();
			
			ix = duplo.first;
			iy = duplo.second;

			spanLeft = spanRight = 0;

			y1 = iy;
			//while(y1 >= 0 && STATUSMAP(ix,y1,iz,NX,NY) == idold) 
			while(y1 >= 0 && read3DVector<short>(status,ix,y1,iz,NX,NY,NZ) == idold) 
				y1--;
			y1++;
						
			//while(y1 < NY && STATUSMAP(ix,y1,iz,NX,NY) == idold)
			while(y1 < NY && read3DVector<short>(status,ix,y1,iz,NX,NY,NZ) == idold)
			{
				
				// class mutex to protect writing
				{  
					#ifdef ENABLE_BOOST_THREADS
						boost::mutex::scoped_lock scopedLock(mutex);
					#endif									

					//STATUSMAP(ix,y1,iz,NX,NY) = idnew;			
					write3DVector<short>(status,idnew,ix,y1,iz,NX,NY,NZ);

					// if this is a cavity then update the cavity list
					if (idnew>=4)
					{
						// current cavity vector
						vector < int* >* vec;

						// first time this cavity is encountered
						if (delphi->cavitiesVec->size()<((unsigned)(idnew-4+1)))
						{
							//cout << endl << INFO << "Detected cavity " << idnew-4;
							vec = new vector<int*>();
							if (vec==NULL)
							{
								cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
								cout.flush();						
								exit(-1);
							}
							delphi->cavitiesVec->push_back(vec);
						}
						else
						{
							vec = delphi->cavitiesVec->at(idnew-4);
						}				
						
						int* v = allocateVector<int>(3);
						if (v==NULL)
						{
							cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
							cout.flush();
							exit(-1);
						}
						v[0] = ix;
						v[1] = y1;
						v[2] = iz;
						//printf("\n %d %d %d (%d,%d)",cix,ciy,ciz,idnew,idold);
						vec->push_back(v);									
					}				
				}
				
				// can go up or not?
				//if (iz<(NZ-1) && STATUSMAP(ix,y1,iz+1,NX,NY) == idold)
				if (iz<(NZ-1) && read3DVector<short>(status,ix,y1,iz+1,NX,NY,NZ) == idold)
				{	
					if ((iz+1)==-1 || (iz+1)<iz_max)
						z_queue->push(pair<pair<int,int>,int >(pair<int,int>(ix,y1),iz+1));
					// postpone
					else
						out_up_xy->push(pair<pair<int,int>,int >(pair<int,int>(ix,y1),iz+1));
				}
				// can go down or not?
				//if (iz>0 && STATUSMAP(ix,y1,iz-1,NX,NY) == idold)
				if (iz>0 && read3DVector<short>(status,ix,y1,iz-1,NX,NY,NZ) == idold)
				{
					if ((iz-1)==-1 || (iz-1)>=iz_min)
						z_queue->push(pair<pair<int,int>,int >(pair<int,int>(ix,y1),iz-1));
					// postpone
					else
						out_down_xy->push(pair<pair<int,int>,int >(pair<int,int>(ix,y1),iz-1));
				}

				//if(!spanLeft && ix > 0 && STATUSMAP(ix-1,y1,iz,NX,NY) == idold)
				if(!spanLeft && ix > 0 && read3DVector<short>(status,ix-1,y1,iz,NX,NY,NZ) == idold)
				{
					xy_stack.push(pair<int,int>(ix-1,y1));
					spanLeft = 1;
				}
				//else if(spanLeft && ix > 0 && STATUSMAP(ix-1,y1,iz,NX,NY) != idold)
				else if(spanLeft && ix > 0 && read3DVector<short>(status,ix-1,y1,iz,NX,NY,NZ) != idold)
				{
					spanLeft = 0;
				}

				//if(!spanRight && ix < NX-1 && STATUSMAP(ix+1,y1,iz,NX,NY) == idold) 
				if(!spanRight && ix < NX-1 && read3DVector<short>(status,ix+1,y1,iz,NX,NY,NZ) == idold) 
				{
					xy_stack.push(pair<int,int>(ix+1,y1));
					spanRight = 1;
				}
				//else if(spanRight && ix < NX-1 && STATUSMAP(ix+1,y1,iz,NX,NY) != idold)
				else if(spanRight && ix < NX-1 && read3DVector<short>(status,ix+1,y1,iz,NX,NY,NZ) != idold)
				{
	                spanRight = 0;
				} 
				y1++;
			}
		}
	}
	
	if (in_queue==NULL)
		delete z_queue;
}

// slightly more cache friendly along y than on x
void Surface::floodFill2(int ix,int iy,int iz,int idold,int idnew)
{	
	queue< pair< pair<int,int>,int> >  z_queue;
	pair<int,int>  duplo;
	pair < pair<int,int>, int>  triplet;
	short* status = delphi->status;
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	z_queue.push(pair<pair<int,int>,int>(pair<int,int>(ix,iy),iz));

	while(z_queue.size()!=0)
	{
		triplet = z_queue.front();
		z_queue.pop();
				
		int iz = triplet.second;
		int ix = triplet.first.first;
		int iy = triplet.first.second;

		//if (STATUSMAP(ix,iy,iz,NX,NY) == idnew)	
		if (read3DVector<short>(status,ix,iy,iz,NX,NY,NZ)==idnew)
			continue;

		stack<pair<int,int> >  xy_stack;
		xy_stack.push(pair<int,int >(ix,iy));
		int y1;

		bool spanLeft=0,spanRight=0;
		
		while(xy_stack.size()!=0)
		{	
			duplo = xy_stack.top();
			xy_stack.pop();
			
			ix = duplo.first;
			iy = duplo.second;

			spanLeft = spanRight = 0;

			y1 = iy;
			//while(y1 >= 0 && STATUSMAP(ix,y1,iz,NX,NY) == idold) 
			while(y1 >= 0 && read3DVector<short>(status,ix,y1,iz,NX,NY,NZ) == idold) 
				y1--;
			y1++;
						
			//while(y1 < NY && STATUSMAP(ix,y1,iz,NX,NY) == idold)
			while(y1 < NY && read3DVector<short>(status,ix,y1,iz,NX,NY,NZ) == idold)
			{
				{  
					//STATUSMAP(ix,y1,iz,NX,NY) = idnew;
					write3DVector<short>(status,idnew,ix,y1,iz,NX,NY,NZ);
					
					// if this is a cavity then update the cavity list
					if (idnew>=4)
					{
						// current cavity vector
						vector < int* >* vec;

						// first time this cavity is encountered
						if (delphi->cavitiesVec->size()<((unsigned)(idnew-4+1)))
						{
							//cout << endl << INFO << "Detected cavity " << idnew-4;
							vec = new vector<int*>();
							if (vec==NULL)
							{
								cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
								cout.flush();						
								exit(-1);
							}
							delphi->cavitiesVec->push_back(vec);
						}
						else
						{
							vec = delphi->cavitiesVec->at(idnew-4);
						}				
						
						int* v = allocateVector<int>(3);
						if (v==NULL)
						{
							cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
							cout.flush();
							exit(-1);
						}
						v[0] = ix;
						v[1] = y1;
						v[2] = iz;
						//printf("\n %d %d %d (%d,%d)",cix,ciy,ciz,idnew,idold);
						vec->push_back(v);					
					}
				}
				// can go up or not?
				//if (iz<NZ-1 && STATUSMAP(ix,y1,iz+1,NX,NY) == idold)
				if (iz<(NZ-1) && read3DVector<short>(status,ix,y1,iz+1,NX,NY,NZ) == idold)
				{						
					z_queue.push(pair<pair<int,int>,int >(pair<int,int>(ix,y1),iz+1));					
				}
				// can go down or not?
				//if (iz>0 && STATUSMAP(ix,y1,iz-1,NX,NY) == idold)
				if (iz>0 && read3DVector<short>(status,ix,y1,iz-1,NX,NY,NZ) == idold)
				{
					z_queue.push(pair<pair<int,int>,int >(pair<int,int>(ix,y1),iz-1));
				}

				//if(!spanLeft && ix > 0 && STATUSMAP(ix-1,y1,iz,NX,NY) == idold)
				if(!spanLeft && ix > 0 && read3DVector<short>(status,ix-1,y1,iz,NX,NY,NZ) == idold)
				{
					xy_stack.push(pair<int,int>(ix-1,y1));
					spanLeft = 1;
				}
				//else if(spanLeft && ix > 0 && STATUSMAP(ix-1,y1,iz,NX,NY) != idold)
				else if(spanLeft && ix > 0 && read3DVector<short>(status,ix-1,y1,iz,NX,NY,NZ) != idold)
				{
					spanLeft = 0;
				}
				//if(!spanRight && ix < NX-1 && STATUSMAP(ix+1,y1,iz,NX,NY) == idold) 
				if(!spanRight && ix < NX-1 && read3DVector<short>(status,ix+1,y1,iz,NX,NY,NZ) == idold) 
				{
					xy_stack.push(pair<int,int>(ix+1,y1));
					spanRight = 1;
				}
				//else if(spanRight && ix < NX-1 && STATUSMAP(ix+1,y1,iz,NX,NY) != idold)
				else if(spanRight && ix < NX-1 && read3DVector<short>(status,ix+1,y1,iz,NX,NY,NZ) != idold)
				{
	                spanRight = 0;
				} 
				y1++;
			}
		}
	}
}

/*
void Surface::floodFill2(int ix,int iy,int iz,int idold,int idnew)
{	
	queue< pair< pair<int,int>,int> >  z_queue;
	pair<int,int>  duplo;
	pair < pair<int,int>, int>  triplet;
	short* status = delphi->status;
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	z_queue.push(pair<pair<int,int>,int>(pair<int,int>(ix,iy),iz));
	
	while(z_queue.size()!=0)
	{
		triplet = z_queue.front();
		z_queue.pop();
		
		stack<pair<int,int> >  xy_stack;
		int iz = triplet.second;
		int ix = triplet.first.first;
		int iy = triplet.first.second;

		//if (STATUSMAP(ix,iy,iz,NX,NY) == idnew)	
		if (read3DVector<short>(status,ix,iy,iz,NX,NY,NZ)==idnew)
			continue;

		//cout << endl << "Slice " << iz << " queue size" << z_queue.size();
		xy_stack.push(pair<int,int >(ix,iy));
		int x1;

		bool spanLeft=0,spanRight=0;
		
		while(xy_stack.size()!=0)
		{	
			duplo = xy_stack.top();
			xy_stack.pop();
			
			ix = duplo.first;
			iy = duplo.second;

			spanLeft = spanRight = 0;

			x1 = ix;
			//while(x1 >= 0 && STATUSMAP(x1,iy,iz,NX,NY) == idold) 
			while(x1 >= 0 && read3DVector<short>(status,x1,iy,iz,NX,NY,NZ) == idold) 
				x1--;
			x1++;
						
			//while(x1 < NX && STATUSMAP(x1,iy,iz,NX,NY) == idold)
			while(x1 < NX && read3DVector<short>(status,x1,iy,iz,NX,NY,NZ) == idold)
			{
				//STATUSMAP(x1,iy,iz,NX,NY) = idnew;
				write3DVector<short>(status,idnew,x1,iy,iz,NX,NY,NZ);
				
				// if this is a cavity then update the cavity list
				if (idnew>=4)
				{
					// current cavity vector
					vector < int* >* vec;

					// first time this cavity is encountered
					if (delphi->cavitiesVec->size()<((unsigned)(idnew-4+1)))
					{
						//cout << endl << INFO << "Detected cavity " << idnew-4;
						vec = new vector<int*>();
						if (vec==NULL)
						{
							cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
							cout.flush();						
							exit(-1);
						}
						delphi->cavitiesVec->push_back(vec);
					}
					else
					{
						vec = delphi->cavitiesVec->at(idnew-4);
					}				
					
					int* v = allocateVector<int>(3);
					if (v==NULL)
					{
						cout << endl << ERR << "Not enough memory to complete cavity detection, stopping";
						cout.flush();
						exit(-1);
					}
					v[0] = x1;
					v[1] = iy;
					v[2] = iz;
					//printf("\n %d %d %d (%d,%d)",cix,ciy,ciz,idnew,idold);
					vec->push_back(v);					
				}
			
				// can go up or not?
				//if (iz<NZ-1 && STATUSMAP(x1,iy,iz+1,NX,NY) == idold)
				if (iz<NZ-1 && read3DVector<short>(status,x1,iy,iz+1,NX,NY,NZ) == idold)
					z_queue.push(pair<pair<int,int>,int >(pair<int,int>(x1,iy),iz+1));					

				// can go down or not?
				//if (iz>0 && STATUSMAP(x1,iy,iz-1,NX,NY) == idold)
				if (iz>0 && read3DVector<short>(status,x1,iy,iz-1,NX,NY,NZ) == idold)
					z_queue.push(pair<pair<int,int>,int >(pair<int,int>(x1,iy),iz-1));

				//if(!spanLeft && iy > 0 && STATUSMAP(x1,iy-1,iz,NX,NY) == idold)
				if(!spanLeft && iy > 0 && read3DVector<short>(status,x1,iy-1,iz,NX,NY,NZ) == idold)
				{
					xy_stack.push(pair<int,int>(x1,iy-1));
					spanLeft = 1;
				}
				//else if(spanLeft && iy > 0 && STATUSMAP(x1,iy-1,iz,NX,NY) != idold)
				else if(spanLeft && iy > 0 && read3DVector<short>(status,x1,iy-1,iz,NX,NY,NZ) != idold)
				{
					spanLeft = 0;
				}
				//if(!spanRight && iy < NY-1 && STATUSMAP(x1,iy+1,iz,NX,NY) == idold) 
				if(!spanRight && iy < NY-1 && read3DVector<short>(status,x1,iy+1,iz,NX,NY,NZ) == idold) 
				{
					xy_stack.push(pair<int,int>(x1,iy+1));
					spanRight = 1;
				}
				//else if(spanRight && iy < NY-1 && STATUSMAP(x1,iy+1,iz,NX,NY) != idold)
				else if(spanRight && iy < NY-1 && read3DVector<short>(status,x1,iy+1,iz,NX,NY,NZ) != idold)
				{
	                spanRight = 0;
				} 
				x1++;
			}
		}
	}
}
*/
void Surface::buildSternLayer()
{
	// for each atom build a bounding cube
	// compute in/out info for each atom
	for (int i=0;i<delphi->numAtoms;i++)
	{
		// compute bounding box object
		double downx = INFINITY;
		double downy = -INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = INFINITY;
		double upz = -INFINITY;
		
		double* sphere_center = delphi->atoms[i]->pos;
		double radius = delphi->atoms[i]->radius+sternLayer;
		double ref = radius*radius;

		downx = sphere_center[0]-radius;
		downy = sphere_center[1]+radius;
		downz = sphere_center[2]-radius;

		upx = sphere_center[0]+radius;
		upy = sphere_center[1]-radius;
		upz = sphere_center[2]+radius;

		// resolve which are the grid cubes that
		// are cut by the bounding box object. These
		// cubes see the cell bounding box    
		int ix_start = (int)rintp((downx-delphi->xmin)/delphi->side);
		int iy_start = (int)rintp((downy-delphi->ymin)/delphi->side);
		int iz_start = (int)rintp((downz-delphi->zmin)/delphi->side);
	     	     
		int ix_end = (int)rintp((upx-delphi->xmin)/delphi->side);
		int iy_end = (int)rintp((upy-delphi->ymin)/delphi->side);
		int iz_end = (int)rintp((upz-delphi->zmin)/delphi->side);

		if (ix_start<0)
		 ix_start=0;
		if (iy_start<0)
		 iy_start=0;
		if (iz_start<0)
		 iz_start=0;

		if (ix_end<0)
		 ix_end=0;
		if (iy_end<0)
		 iy_end=0;
		if (iz_end<0)
		 iz_end=0;

		if (ix_end>=(int)delphi->nx)
		 ix_end=delphi->nx-1;
		if (iy_end>=(int)delphi->ny)
		 iy_end=delphi->ny-1;	          
		if (iz_end>=(int)delphi->nz)
		 iz_end=delphi->nz-1;	          

		if (ix_start>=(int)delphi->nx)
		 ix_start=delphi->nx-1;
		if (iy_start>=(int)delphi->ny)
		 iy_start=delphi->ny-1;	          
		if (iz_start>=(int)delphi->nz)
		 iz_start=delphi->nz-1;	          

		for (int iz=iz_start;iz<=iz_end;iz++)
		 for (int iy=iy_start;iy>=iy_end;iy--)
			 for (int ix=ix_start;ix<=ix_end;ix++)
			 {
				 // if already inside skip
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
					if (dist<ref)
					{
						//delphi->IDEBMAP(ix,iy,iz,delphi->nx,delphi->ny) = false;
						write3DVector<bool>(delphi->idebmap,false,ix,iy,iz,delphi->nx,delphi->ny,delphi->nz);
					}
				 }
			 }				 
	}
}
bool Surface::getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3)
{
	cout << endl << WARN << "Projection not supported!";
	return false;
}

void Surface::getRayIntersection(double p1[3],double p2[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals)
{
	cout << endl << WARN << "Ray-Patch intersection not supported!";
}


void Surface::projector(int start,int end)
{	
	for (int i=start;i<end;i++)
	{
		double gridPoint[3];
		gridPoint[0] = delphi->x[delphi->ibgp[3*i]];
		gridPoint[1] = delphi->y[delphi->ibgp[3*i+1]];
		gridPoint[2] = delphi->z[delphi->ibgp[3*i+2]];

		// project on the surface
		if (bgp_type[i]==EXTERNAL_BGP)
		{		
			getProjection(gridPoint,&(delphi->scspos[3*i]),
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
			// move from grid point of delphi grid to auxiliary atoms grid
			int ix = (int)rintp((gridPoint[0]-gxmin)/gside);
			int iy = (int)rintp((gridPoint[1]-gymin)/gside);
			int iz = (int)rintp((gridPoint[2]-gzmin)/gside);

			double minDist=INFINITY,secondDist=INFINITY;
			int first=-1,second=-1;

			// get the nearest atoms to set the dielectric constant
			// the dielectric constant is mapped according to the additively weighted voronoi diagram
			// that is the signed distance from the point p is ||p-c||^2-r^2 where c is the center
			// of the atom and r is the radius. The minimum signed distance wins.
			for (int k=0;k<SHIFT_MAP;k++)
			{
				unsigned int cx = ix+shift_map[k][0];
				unsigned int cy = iy+shift_map[k][1];
				unsigned int cz = iz+shift_map[k][2];

				if (cx>(ggrid-1) || cy>(ggrid-1) || cz>(ggrid-1) || cx<0 || cy<0 || cz<0)
					continue;

				//int max_ind = GRID_MULTI_MAP(cx,cy,cz,0,ggrid,ggrid,ggrid);
				int max_ind = read4DVector<int>(gridMultiMap,cx,cy,cz,0,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

				for (int j=0;j<max_ind;j++)
				{
					//int atom_index = GRID_MULTI_MAP(cx,cy,cz,j+1,ggrid,ggrid,ggrid);
					int atom_index = read4DVector<int>(gridMultiMap,cx,cy,cz,j+1,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);
					double signed_dist = 0;
					DIST2(signed_dist,delphi->atoms[atom_index]->pos,gridPoint)
					double radius2 = delphi->atoms[atom_index]->radius2;
					signed_dist-=radius2;

					//printf("\n%d %f %d %d",atom_index,signed_dist,first,second);

					if (signed_dist<minDist)
					{
						secondDist = minDist;
						minDist = signed_dist;
						second = first;
						first = atom_index;
					}
					else
					{
						if (signed_dist<secondDist)
						{
							secondDist = signed_dist;
							second = atom_index;
						}					
					}
				}	
			}
			
			//exit(-1);
			// compute the voronoi plane using power distance
			double w[4];
			double *c1,*c2;
			c1 = delphi->atoms[first]->pos;
			c2 = delphi->atoms[second]->pos;
			SUB(w,c2,c1)
			w[0] = 2*w[0];
			w[1] = 2*w[1];
			w[2] = 2*w[2];
			w[3] = (DOT(c1,c1))-(DOT(c2,c2))-delphi->atoms[first]->radius*delphi->atoms[first]->radius+delphi->atoms[second]->radius*delphi->atoms[second]->radius;

			//printf("\n %f",w[3]);
			// project on the plane
			double dist,proj[3];
			point2plane(gridPoint,w,&dist,proj);
			
			//fprintf(fp,"%f %f %f\n",gridPoint[0],gridPoint[1],gridPoint[2]);
			//printf("\n%f %f %f",proj[0],proj[1],proj[2]);
			
			delphi->scspos[3*i]=proj[0];
			delphi->scspos[3*i+1]=proj[1];
			delphi->scspos[3*i+2]=proj[2];

			//debug check
			/*
			delphi->scspos[3*i]=gridPoint[0];
			delphi->scspos[3*i+1]=gridPoint[1];
			delphi->scspos[3*i+2]=gridPoint[2];
			*/
			delphi->scsnor[3*i] = gridPoint[0]-proj[0];
			delphi->scsnor[3*i+1] = gridPoint[1]-proj[1];
			delphi->scsnor[3*i+2] = gridPoint[2]-proj[2];			
		}
	}
}


void Surface::intersector(double* volPanel,int nb,int start,int end,int jump,int* numIntersections,packet pack)
{
	double diff[3];
	vector<pair < double, double* > > intersections;	
	intersections.reserve(2000);

	vector<coordVec>* v_int = pack.first;
	vector<coordVec>* v_norm = pack.second;

	volPanel[panel]=0;
	double pa[3],pb[3];

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	short* status = delphi->status;
	
	// if cavities and espilon map is not necessary then skip this step at all
	if (!delphi->buildStatus && !delphi->buildEpsMap)
	{
		// skip
	}

	// if only cavity detection is required, only the first panel is ray cast
	else if (delphi->buildStatus && !delphi->buildEpsMap && panel!=0) 
	{
		// skip panel 1,2
	}

	// if both cavity detection and epsilon map are required, 3 panels must be analyzed
	else
	{
		for (int i=start;i<end;i+=jump)
		{				
			int trials = 0;
			for (int j=0;j<nb;j++)
			{
				//printf("\n %d %d %d",i,j,trials);
				// light source and ray along coordinate direction
				if (panel==0)
				{
					pa[0] = delphi->x[0];
					pa[1] = delphi->y[i];        
					pa[2] = delphi->z[j];
					pb[0] = delphi->x[delphi->nx-1];
					pb[1] = pa[1];
					pb[2] = pa[2];
				}
				else if (panel==1)
				{										
					pa[0] = delphi->x[j];
					pa[1] = delphi->y[i];        
					pa[2] = delphi->z[0];
					pb[0] = pa[0];
					pb[1] = pa[1];
					pb[2] = delphi->z[delphi->nz-1];
				}
				else
				{
					pa[0] = delphi->x[j];
					pa[2] = delphi->z[i];
					pa[1] = delphi->y[0];   
					pb[0] = pa[0];	
					pb[1] = delphi->y[delphi->ny-1];
					pb[2] = pa[2];
				}
												
				SUB(diff,pb,pa);
			
				intersections.clear();

				getRayIntersection(pa,pb,intersections,0,false);

				if (intersections.size()==0)
					continue;
					
				// check if repetions allow the ray to enter and exist. If the
				// ray only enters then this is an error
				vector<pair<double,double*> >::iterator itt = intersections.begin();
				unsigned int vsize = (unsigned int)intersections.size(),ind=0;
				bool closure = false;
				double lastValid = INFINITY;
				double a;
				
				for (;itt!=intersections.end();itt++)
				{
					a = (*itt).first;
					if (fabs(a-lastValid)>=EPS_INT)
					{
						// ok the ray is entering
						closure = false;
						// search the exiting point of the ray
						while(1)
						{							
							ind++;							
							if (ind==vsize)
								break;
							itt++;

							if (fabs(a-(*itt).first)>=EPS_INT)
							{
								// ok the ray is exiting 

								//printf(" (ok!)");
								closure = true;
								lastValid = (*itt).first;
								break;	
							}
						}			
					}
					ind++;
				}

				/* single or multiple tangent intersection.*/
				if (lastValid==INFINITY)
				{
					//cout << " * ";				
					continue;
				}
				
				// check finished.
				if (!closure)
				{
					// approximate the current ray with the previous one
						if (panel==0)
						{
							if (delphi->buildStatus && delphi->buildEpsMap)
							{
								for (int ix=0;ix<NX;ix++)
								{
									//delphi->EPSMAP((ix),(i),(j),0,NX,NY,NZ)= delphi->EPSMAP((ix),(i),(j-1),0,NX,NY,NZ);
									const int eee = read4DVector<int>(delphi->epsmap,ix,i,j-1,0,NX,NY,NZ,3);
									write4DVector<int>(delphi->epsmap,eee,ix,i,j,0,NX,NY,NZ,3);

									//delphi->status[ix][i][j] = delphi->status[ix][i][j-1];
									//STATUSMAP(ix,i,j,NX,NY) = STATUSMAP(ix,i,(j-1),NX,NY);
									const short sss = read3DVector<short>(delphi->status,ix,i,j-1,NX,NY,NZ);
									write3DVector<short>(delphi->status,sss,ix,i,j,NX,NY,NZ);
								}
							}
							else if (delphi->buildStatus)
							{
								for (int ix=0;ix<NX;ix++)
								{
									//delphi->status[ix][i][j] = delphi->status[ix][i][j-1];
									//STATUSMAP(ix,i,j,NX,NY) = STATUSMAP(ix,i,(j-1),NX,NY);
									const short sss = read3DVector<short>(delphi->status,ix,i,j-1,NX,NY,NZ);
									write3DVector<short>(delphi->status,sss,ix,i,j,NX,NY,NZ);
								}
							}
							// in this case idebmap is directly written because it is no more derived from status map
							else if (delphi->buildEpsMap)
							{
								for (int ix=0;ix<NX;ix++)						
								{
									//delphi->EPSMAP((ix),(i),(j),0,NX,NY,NZ)= delphi->EPSMAP((ix),(i),(j-1),0,NX,NY,NZ);
									const int eee = read4DVector<int>(delphi->epsmap,ix,i,j-1,0,NX,NY,NZ,3);
									write4DVector<int>(delphi->epsmap,eee,ix,i,j,0,NX,NY,NZ,3);

									//delphi->IDEBMAP(ix,i,j,NX,NY)=delphi->IDEBMAP(ix,i,j-1,NX,NY);
									const bool bbb = read3DVector<bool>(delphi->idebmap,ix,i,j-1,NX,NY,NZ);
									write3DVector<bool>(delphi->idebmap,bbb,ix,i,j,NX,NY,NZ);
								}
							}
						}
						else if (panel==1)
						{
							for (int iz=0;iz<NZ;iz++)
							{
								//delphi->EPSMAP((j),(i),(iz),2,NX,NY,NZ)=delphi->EPSMAP((j-1),(i),(iz),2,NX,NY,NZ);
								const int eee = read4DVector<int>(delphi->epsmap,j-1,i,iz,2,NX,NY,NZ,3);
								write4DVector<int>(delphi->epsmap,eee,j,i,iz,2,NX,NY,NZ,3);
							}
						}
						else if (panel==2)
						{
							for (int iy=0;iy<NY;iy++)
							{
								//delphi->EPSMAP((j),(iy),(i),1,NX,NY,NZ)=delphi->EPSMAP((j-1),(iy),(i),1,NX,NY,NZ);	
								const int eee = read4DVector<int>(delphi->epsmap,j-1,iy,i,1,NX,NY,NZ,3);
								write4DVector<int>(delphi->epsmap,eee,j,iy,i,1,NX,NY,NZ,3);
							}
						}
				
						(*numIntersections)++;
						continue;
					
				}
											
				// if there are intersections mark the epsmap cells that are inside
				// by counting the number of intersections. This is done for each panel
				// that is in the x,y,z directions respectively.
				// Eg. Tracing the ray one obtains inside/out information
				// 000000000|-1-1-1-1-1-1-1-1-1-1-1-1|000000000000|-1-1-1-1-1|000000			
				if (intersections.size()>1)
				{					
					if (panel==0)
					{
						double lim1,lim2;
						double lastValidLim2 = INFINITY;					

						for (unsigned int ig=0;ig<intersections.size();ig+=2)
						{										
							if ((ig+1)>=vsize)
								break;
							
							lim1 = pa[0]+diff[0]*intersections[ig].first;
							lim2 = pa[0]+diff[0]*intersections[ig+1].first;
						
							// skip repeated intersections 
							//if (lim2 == lim1 || lim1 == lastValidLim2)
							if (fabs(lim2-lim1)<EPS_INT || fabs(lim1-lastValidLim2)<EPS_INT)
							{																
								ig--;
								continue;
							}
							
							//printf("\n lim2 %f lim1 %f",lim2,lim1);
							volPanel[panel]+=lim2-lim1;												
							lastValidLim2 = lim2;
													
							if (delphi->buildEpsMap && delphi->buildStatus)
							{
								for (int ix=0;ix<NX;ix++)
								{			
									if (delphi->x[ix]+delphi->hside<=lim2 && delphi->x[ix]+delphi->hside>=lim1)
									{
										//delphi->EPSMAP((ix),(i),(j),0,NX,NY,NZ)=inside;
										write4DVector<int>(delphi->epsmap,inside,ix,i,j,0,NX,NY,NZ,3);
									}

									if (delphi->x[ix]<=lim2 && delphi->x[ix]>=lim1)
									{
										//delphi->status[ix][i][j]=STATUS_POINT_INSIDE;
										//STATUSMAP(ix,i,j,NX,NY)=STATUS_POINT_INSIDE;
										write3DVector<short>(delphi->status,STATUS_POINT_INSIDE,ix,i,j,NX,NY,NZ);
									}
								}
							}
							else if (delphi->buildEpsMap)
							{
								for (int ix=0;ix<NX;ix++)
								{	
									// also idebmap is directly written because in this case status map is not used
									if (delphi->x[ix]+delphi->hside<=lim2 && delphi->x[ix]+delphi->hside>=lim1)
									{
										//delphi->EPSMAP((ix),(i),(j),0,NX,NY,NZ)=inside;									
										write4DVector<int>(delphi->epsmap,inside,ix,i,j,0,NX,NY,NZ,3);
									}

									if (delphi->x[ix]<=lim2 && delphi->x[ix]>=lim1)
									{
										//delphi->IDEBMAP((ix),(i),(j),NX,NY)=false;
										write3DVector<bool>(delphi->idebmap,false,ix,i,j,NX,NY,NZ);
									}
								}
							}
							else if (delphi->buildStatus)
							{
								for (int ix=0;ix<NX;ix++)
									if (delphi->x[ix]<=lim2 && delphi->x[ix]>=lim1)
									{
										//delphi->status[ix][i][j]=STATUS_POINT_INSIDE;
										//STATUSMAP(ix,i,j,NX,NY)=STATUS_POINT_INSIDE;
										write3DVector<short>(delphi->status,STATUS_POINT_INSIDE,ix,i,j,NX,NY,NZ);
									}
							}
						}
					}
					else if (panel==1)
					{
						double lim1,lim2;
						double lastValidLim2 = INFINITY;
						
						for (unsigned int ig=0;ig<intersections.size();ig+=2)
						{															
							if ((ig+1)>=vsize)
								break;
							
							lim1 = pa[2]+diff[2]*intersections[ig].first;
							lim2 = pa[2]+diff[2]*intersections[ig+1].first;						
					
							// skip repeated intersections 
							//if (lim2 == lim1 || lim1 == lastValidLim2)
							if (fabs(lim2-lim1)<EPS_INT || fabs(lim1-lastValidLim2)<EPS_INT)
							{
								ig--;
								continue;
							}
							
							volPanel[panel]+=lim2-lim1;
							lastValidLim2 = lim2;
					
							for (int iz=0;iz<NZ;iz++)
							{						
								if (delphi->z[iz]+delphi->hside<=lim2 && delphi->z[iz]+delphi->hside>=lim1)
								{
									//delphi->EPSMAP((j),(i),(iz),2,NX,NY,NZ)=inside;
									write4DVector<int>(delphi->epsmap,inside,j,i,iz,2,NX,NY,NZ,3);
								}
							}
						}
					}
					else
					{
						double lim1,lim2;
						double lastValidLim2 = INFINITY;

						for (unsigned int ig=0;ig<intersections.size();ig+=2)
						{										
							if ((ig+1)>=vsize)
								break;						

							lim1 = pa[1]+diff[1]*intersections[ig].first;
							lim2 = pa[1]+diff[1]*intersections[ig+1].first;

							// skip repeated intersections 											
							//if (lim2 == lim1 || lim1 == lastValidLim2)
							if (fabs(lim2-lim1)<EPS_INT || fabs(lim1-lastValidLim2)<EPS_INT)
							{
								ig--;
								continue;
							}	
							
							volPanel[panel]+=lim2-lim1;									
							lastValidLim2 = lim2;
							
							for (int iy=0;iy<NY;iy++)
							{							
								if (delphi->y[iy]+delphi->hside<=lim2 && delphi->y[iy]+delphi->hside>=lim1)									
								{
									//delphi->EPSMAP((j),(iy),(i),1,NX,NY,NZ)=inside;
									write4DVector<int>(delphi->epsmap,inside,j,iy,i,1,NX,NY,NZ,3);
								}
							}
						}												
					}									
				}
			}
		}
	}
		
	//set<string> orphans;
	// if a triangulation is needed and we have not already a scalar field,
	// more rays need to be casted
	if (accurateTriangulation && !isAvailableScalarField)
	{
		double delta = delta_accurate_triangulation;		
		for (int i=start;i<end;i+=jump)
		{				
			int trials = 0;
			for (int j=0;j<nb;j++)
			{
				if (panel==0)
				{
					pa[0] = delphi->x[0]-delphi->hside+delta;
					pa[1] = delphi->y[i]-delphi->hside+delta;        
					pa[2] = delphi->z[j]-delphi->hside+delta;
					pb[0] = delphi->x[delphi->nx-1]-delphi->hside+delta;
					pb[1] = pa[1];
					pb[2] = pa[2];
				}
				else if (panel==1)
				{										
					pa[0] = delphi->x[j]-delphi->hside+delta;
					pa[1] = delphi->y[i]-delphi->hside+delta;        
					pa[2] = delphi->z[0]-delphi->hside+delta;
					pb[0] = pa[0];
					pb[1] = pa[1];
					pb[2] = delphi->z[delphi->nz-1]-delphi->hside+delta;
				}
				else
				{
					pa[0] = delphi->x[j]-delphi->hside+delta;
					pa[1] = delphi->y[0]-delphi->hside+delta;   
					pa[2] = delphi->z[i]-delphi->hside+delta;					
					pb[0] = pa[0];	
					pb[1] = delphi->y[delphi->ny-1]-delphi->hside+delta;
					pb[2] = pa[2];
				}
						
				SUB(diff,pb,pa);
				intersections.clear();
				getRayIntersection(pa,pb,intersections,0,computeNormals);

				if (intersections.size()==0)
					continue;

				vector<pair<double,double*> >::iterator itt = intersections.begin();
				unsigned int vsize = (unsigned int)intersections.size(),ind=0;
				bool closure = false;
				double lastValid = INFINITY;
				double a;
				
				for (;itt!=intersections.end();itt++)
				{
					a = (*itt).first;
					//printf("\n try to enter at -> %f",(*itt).first);
					// get the entering point
					//if (a!=lastValid)
					if (fabs(a-lastValid)>=EPS_INT)
					{
						//printf(" ok!");

						// ok the ray is entering
						closure = false;
						// search the exiting point of the ray
						while(1)
						{							
							ind++;							
							if (ind==vsize)
								break;
							itt++;
							
							//printf("\n try to exit at -> %f",(*itt).first);
							//if (a!=(*itt).first)
							if (fabs(a-(*itt).first)>=EPS_INT)
							{
								// ok the ray is exiting 

								//printf(" (ok!)");
								closure = true;
								lastValid = (*itt).first;
								break;	
							}
						}			
					}
					ind++;
				}

				/* single or multiple tangent intersection.*/
				if (lastValid==INFINITY)
				{
					//cout << " * ";				
					continue;
				}
								
				// check finished.
				if (!closure)
				{						
					// approximate the current ray with the previous one
					// analytical intersections cannot be recovered, only in/out position is recorded
					// The following marching cubes will be semi analytical. Semi means that were the analytical intersection
					// is not present the usual marching cubes rule will be used																	
						if (panel==0)
							for (int ix=0;ix<delphi->nx;ix++)						
							{
								verticesInsidenessMap[ix][i][j] = verticesInsidenessMap[ix][i][j-1];							
								//bool val = read3DVector<bool>(verticesInsidenessMap,ix,i,j-1,delphi->nx,delphi->ny,delphi->nz);
								//write3DVector<bool>(verticesInsidenessMap,val,ix,i,j,delphi->nx,delphi->ny,delphi->nz);
							}
						if (panel==1)
							for (int iz=0;iz<delphi->nz;iz++)							
							{
								verticesInsidenessMap[j][i][iz] = verticesInsidenessMap[j-1][i][iz];
								//bool val = read3DVector<bool>(verticesInsidenessMap,j-1,i,iz,delphi->nx,delphi->ny,delphi->nz);
								//write3DVector<bool>(verticesInsidenessMap,val,j,i,iz,delphi->nx,delphi->ny,delphi->nz);
							}
						if (panel==2)
							for (int iy=0;iy<delphi->ny;iy++)							
							{
								verticesInsidenessMap[j][iy][i] = verticesInsidenessMap[j-1][iy][i];
								//bool val = read3DVector<bool>(verticesInsidenessMap,j-1,iy,i,delphi->nx,delphi->ny,delphi->nz);
								//write3DVector<bool>(verticesInsidenessMap,val,j,iy,i,delphi->nx,delphi->ny,delphi->nz);
							}
						(*numIntersections)++;							
						continue;
				}

				// if there are intersections mark the epsmap cells that are inside
				// by counting the number of intersections. This is done for each panel
				// that is in the x,y,z directions respectively.
				// Eg. Tracing the ray one obtains inside/out information
				// 000000000|-1-1-1-1-1-1-1-1-1-1-1-1|000000000000|-1-1-1-1-1|000000			
				if (intersections.size()>1)
				{					
					if (panel==0)
					{
						bool state  = false;
						double lim1,lim2;
						double lastValidLim2 = INFINITY;					

						//printf("\n-----------");
						for (unsigned int ig=0;ig<intersections.size();ig+=2)
						{										
							if ((ig+1)>=vsize)
								break;
							
							lim1 = pa[0]+diff[0]*intersections[ig].first;
							lim2 = pa[0]+diff[0]*intersections[ig+1].first;
						
							//printf("\n lim2 %f lim1 %f",lim2,lim1);

							// skip repeated intersections 
							//if (lim2 == lim1 || lim1 == lastValidLim2)
							if (fabs(lim2-lim1)<EPS_INT || fabs(lim1-lastValidLim2)<EPS_INT)
							{																
								ig--;
								continue;
							}
							
							lastValidLim2 = lim2;

							// if no cavity detection and no epsmap is enabled the volume must be computed here
							if (!delphi->buildStatus && !delphi->buildEpsMap)
								volPanel[panel]+=lim2-lim1;									

							double *intersec1,*intersec2;
							intersec1 = allocateVector<double>(3);
							intersec2 = allocateVector<double>(3);
							intersec1[0] = lim1;
							intersec1[1] = pa[1];
							intersec1[2] = pa[2];

							intersec2[0] = lim2;
							intersec2[1] = pa[1];
							intersec2[2] = pa[2];			
														
							// get the cube which the intersections belong. 
							int ix_ = (int)rintp((intersec1[0]-delphi->xmin)/delphi->side);
							int xa = ix_+1;
							int iy_ = i;
							int iz_ = j;
							int xb = (int)rintp((intersec2[0]-delphi->xmin)/delphi->side);

							// high frequency intersection
							if (xb<xa)
							{
								deleteVector<double>(intersec1);
								deleteVector<double>(intersec2);
								continue;
							}

							for (int y=xa;y<=xb;y++)
								verticesInsidenessMap[y][i][j]=false;
								//write3DVector<bool>(verticesInsidenessMap,false,y,i,j,delphi->nx,delphi->ny,delphi->nz);

																					
							{
							coordVec cv(ix_,iy_,iz_,intersec1,X_DIR);
							v_int->push_back(cv);
							}
											
							if (computeNormals && providesAnalyticalNormals)
							{
								double* n = intersections[ig].second;								
								if (n!=NULL)
								{
									coordVec cv(ix_,iy_,iz_,n,X_DIR);
									v_norm->push_back(cv);
								}
							}
														
							ix_ = (int)rintp((intersec2[0]-delphi->xmin)/delphi->side);
							iy_ = i;
							iz_ = j;

							{
							coordVec cv(ix_,iy_,iz_,intersec2,X_DIR);
							v_int->push_back(cv);
							}

							if (computeNormals && providesAnalyticalNormals)
							{
								double* n = intersections[ig+1].second;
								if (n!=NULL)
								{
									coordVec cv(ix_,iy_,iz_,n,X_DIR);
									v_norm->push_back(cv);
								}
							}							
						}
						//printf("\n-----------");
					}
					else if (panel==1)
					{
						double lim1,lim2;
						double lastValidLim2 = INFINITY;
						
						for (unsigned int ig=0;ig<intersections.size();ig+=2)
						{															
							if ((ig+1)>=vsize)
								break;
							
							lim1 = pa[2]+diff[2]*intersections[ig].first;
							lim2 = pa[2]+diff[2]*intersections[ig+1].first;						
					
							// skip repeated intersections 
							//if (lim2 == lim1 || lim1 == lastValidLim2)
							if (fabs(lim2-lim1)<EPS_INT || fabs(lim1-lastValidLim2)<EPS_INT)
							{
								ig--;
								continue;
							}
							
							
							lastValidLim2 = lim2;

							// if no cavity detection and no epsmap is enabled the volume must be computed here
							if (!delphi->buildStatus && !delphi->buildEpsMap)
								volPanel[panel]+=lim2-lim1;									

							double *intersec1,*intersec2;
							intersec1 = allocateVector<double>(3);
							intersec2 = allocateVector<double>(3);
							intersec1[0] = pa[0] ;
							intersec1[1] = pa[1] ;
							intersec1[2] = lim1 ;

							intersec2[0] = pa[0] ;
							intersec2[1] = pa[1] ;
							intersec2[2] = lim2 ;								
							
							// get the cube which the intersections belong. The little 1e-7 is 
							// used in order to be sure to get the right cube. Consider that
							// we are at the edge so by performing a little move on the right
							// direction we are sure to identify the right cube
							int iz_ = (int)rintp((intersec1[2]-delphi->zmin)/delphi->side);
							int za = iz_+1;
							int ix_ = j;
							int iy_ = i;
							int zb = (int)rintp((intersec2[2]-delphi->zmin)/delphi->side);

							if (zb<za)
							{
								deleteVector<double>(intersec1);
								deleteVector<double>(intersec2);
								continue;
							}

							
							{
							coordVec cv(ix_,iy_,iz_,intersec1,Z_DIR);
							v_int->push_back(cv);
							}

							if (computeNormals && providesAnalyticalNormals)
							{
								double* n = intersections[ig].second;
								if (n!=NULL)
								{
									coordVec cv(ix_,iy_,iz_,n,Z_DIR);
									v_norm->push_back(cv);
								}
							}
							
							iz_ = (int)rintp((intersec2[2]-delphi->zmin)/delphi->side);
							ix_ = j;
							iy_ = i;
							
							for (int y=za;y<=zb;y++)
								verticesInsidenessMap[j][i][y]=false;
								//write3DVector<bool>(verticesInsidenessMap,false,j,i,y,delphi->nx,delphi->ny,delphi->nz);


							{
							coordVec cv(ix_,iy_,iz_,intersec2,Z_DIR);
							v_int->push_back(cv);
							}

							if (computeNormals && providesAnalyticalNormals)
							{
								double* n = intersections[ig+1].second;
								if (n!=NULL)
								{
									coordVec cv(ix_,iy_,iz_,n,Z_DIR);
									v_norm->push_back(cv);
								}
							}
													

						}
					}
					else
					{
						double lim1,lim2;
						double lastValidLim2 = INFINITY;

						for (unsigned int ig=0;ig<intersections.size();ig+=2)
						{										
							if ((ig+1)>=vsize)
								break;						

							lim1 = pa[1]+diff[1]*intersections[ig].first;
							lim2 = pa[1]+diff[1]*intersections[ig+1].first;

							// skip repeated intersections 											
							//if (lim2 == lim1 || lim1 == lastValidLim2)
							if (fabs(lim2-lim1)<EPS_INT || fabs(lim1-lastValidLim2)<EPS_INT)
							{
								ig--;
								continue;
							}	
							
							lastValidLim2 = lim2;

							// if no cavity detection and no epsmap is enabled the volume must be computed here
							if (!delphi->buildStatus && !delphi->buildEpsMap)
								volPanel[panel]+=lim2-lim1;									

							double *intersec1,*intersec2;
							intersec1 = allocateVector<double>(3);
							intersec2 = allocateVector<double>(3);
							intersec1[0] = pa[0] ;
							intersec1[1] = lim1 ;
							intersec1[2] = pa[2] ;

							intersec2[0] = pa[0] ;
							intersec2[1] = lim2 ;
							intersec2[2] = pa[2];								
							
							// get the cube which the intersections belong. The little 1e-7 is 
							// used in order to be sure to get the right cube. Consider that
							// we are at the edge so by performing a little move on the right
							// direction we are sure to identify the right cube
							int iy_ = (int)rintp((intersec1[1]-delphi->ymin)/delphi->side);
							int ix_ = j;
							int iz_ = i;
							int ya = iy_+1;
							int yb = (int)rintp((intersec2[1]-delphi->ymin)/delphi->side);

							if (yb<ya)
							{
								deleteVector<double>(intersec1);
								deleteVector<double>(intersec2);
								continue;
							}
		
							{
							coordVec cv(ix_,iy_,iz_,intersec1,Y_DIR);
							v_int->push_back(cv);
							}

							if (computeNormals && providesAnalyticalNormals)
							{
								double* n = intersections[ig].second;
								if (n!=NULL)
								{
									coordVec cv(ix_,iy_,iz_,n,Y_DIR);
									v_norm->push_back(cv);
								}
							}
							

							//ix_ = (int)rintp((intersec2[0]-delphi->xmin)/delphi->side);
							iy_ = (int)rintp((intersec2[1]-delphi->ymin)/delphi->side);

							
							//iz_ = (int)rintp((intersec2[2]-delphi->zmin)/delphi->side);

							ix_ = j;
							iz_ = i;
							

							for (int y=ya;y<=yb;y++)
								verticesInsidenessMap[j][y][i]=false;
								//write3DVector<bool>(verticesInsidenessMap,false,j,y,i,delphi->nx,delphi->ny,delphi->nz);

							{
							coordVec cv(ix_,iy_,iz_,intersec2,Y_DIR);
							v_int->push_back(cv);
							}

							if (computeNormals && providesAnalyticalNormals)
							{
								double* n = intersections[ig+1].second;
								if (n!=NULL)
								{
									coordVec cv(ix_,iy_,iz_,n,Y_DIR);
									v_norm->push_back(cv);
								}
							}
													
						}
					}
				}
			}
		}
	}

}  

inline void Surface::getVertices(double isolevel,int start_z,int end_z,int jump,vector<coordVec>* vertList,vector<double*>* normalsList)
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
	
	double** positions = allocateMatrix2D<double>(8,3);	
	double* gx = delphi->x;
	double* gy = delphi->y;
	double* gz = delphi->z;
	double d_hside = delphi->hside;
		
	for (int k=start_z;k<end_z;k+=jump)	
	{	
		for (int i=0;i<NX;i++)
			for (int j=0;j<NY;j++)			
			{

				//cout << endl << i << " " << j << " " << k;

				// skip the boundary of the grid, there will never be triangles here if
				// the grid is correctly built
				if (i==0 || i==(NX-1)|| j==0 || j==(NY-1) || k==0 || k==(NZ-1))				
					continue;

				double votes[8];
				double xg;
				double yg;
				double zg;
				double table[6];

				// if an accurate description is available each inside/outside flag
				// and edge surface intersections were analytically computed by the ray-tracing routine
				// this allows to perform a marching cubes where each vertex of the mesh is granted to
				// to belong to the surface. 
				// This is an analytical marching cubes				
				if (accurateTriangulation && !isAvailableScalarField)
				{					
					if (verticesInsidenessMap[i][j][k])					
						votes[0] = +1;
					else
						votes[0] = -1;
					
					//bool val = read3DVector<bool>(verticesInsidenessMap,i,j,k,NX,NY,NZ);
					//votes[0] = ((val) ? (+1) : (-1));
					
					
					if (verticesInsidenessMap[i+1][j][k])
						votes[1] = +1;
					else
						votes[1] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i+1,j,k,NX,NY,NZ);
					//votes[1] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i+1][j+1][k])
						votes[2] = +1;
					else
						votes[2] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i+1,j+1,k,NX,NY,NZ);
					//votes[2] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i][j+1][k])
						votes[3] = +1;
					else
						votes[3] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i,j+1,k,NX,NY,NZ);
					//votes[3] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i][j][k+1])
						votes[4] = +1;
					else
						votes[4] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i,j,k+1,NX,NY,NZ);
					//votes[4] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i+1][j][k+1])
						votes[5] = +1;
					else
						votes[5] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i+1,j,k+1,NX,NY,NZ);
					//votes[5] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i+1][j+1][k+1])
						votes[6] = +1;
					else
						votes[6] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i+1,j+1,k+1,NX,NY,NZ);
					//votes[6] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i][j+1][k+1])
						votes[7] = +1;
					else
						votes[7] = -1;
					

					//val = read3DVector<bool>(verticesInsidenessMap,i,j+1,k+1,NX,NY,NZ);
					//votes[7] = ((val) ? (+1) : (-1));
				}
				// Classical marching cube interpolating scalar field values
				else if (accurateTriangulation && isAvailableScalarField)
				{
					votes[0] = scalarField[i][j][k];
					votes[1] = scalarField[i+1][j][k];
					votes[2] = scalarField[i+1][j+1][k];
					votes[3] = scalarField[i][j+1][k];
					votes[4] = scalarField[i][j][k+1];
					votes[5] = scalarField[i+1][j][k+1];
					votes[6] = scalarField[i+1][j+1][k+1];
					votes[7] = scalarField[i][j+1][k+1];

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
					for (int vertInd =0;vertInd<8;vertInd++)
						votes[vertInd] = getInsidness(i,j,k,vertInd);					

				if (accurateTriangulation && isAvailableScalarField)
				{
					for (int vertInd =0;vertInd<8;vertInd++)
					{
						positions[vertInd][0]= table[mulTable[vertInd][0]];
						positions[vertInd][1]= table[mulTable[vertInd][1]];
						positions[vertInd][2]= table[mulTable[vertInd][2]];
					}
				}
				
				int cubeindex = classifyCube(votes,isolevel);
		
				if (cubeindex!=-1)
					//activeCubes[i][j][k]=true;				
					write3DVector<bool>(activeCubes,true,i,j,k,NX,NY,NZ);
				else
					continue;			
										
				int vert_offset;
				int ix=i,iy=j,iz=k;

				// 5
				if (edgeTable[cubeindex] & 32)
				{
					vert_offset = -1;
					vert_offset	= intersectionsMatrixAlongY->at(ix+1,iy,iz+1);
										
					if (vert_offset==-1)
					{												
						double* p = allocateVector<double>(3);

						if (accurateTriangulation && isAvailableScalarField)
							vertexInterp(isolevel,positions[5],positions[6],votes[5],votes[6],p);
						else
						{
							p[0] = delphi->x[ix+1]-delphi->hside;
							p[1] = delphi->y[iy];
							p[2] = delphi->z[iz+1]-delphi->hside;
						}

						coordVec v(ix+1,iy,iz+1,p,Y_DIR);
						vertList->push_back(v);
						
						// mark that this normal should be approximated later by the triangulation
						if (computeNormals && providesAnalyticalNormals)
						{
							double* normal = allocateVector<double>(3);
							normal[0]=0;
							normal[1]=0;
							normal[2]=0;
							normalsList->push_back(normal);
						}
					}								
				}

				// 6
				if (edgeTable[cubeindex] & 64)
				{
					vert_offset = -1;
					vert_offset	=	intersectionsMatrixAlongX->at(ix,iy+1,iz+1);					
					
					if (vert_offset==-1)
					{					
						double* p = allocateVector<double>(3);
						if (accurateTriangulation && isAvailableScalarField)
							vertexInterp(isolevel,positions[6],positions[7],votes[6],votes[7],p);						
						else
						{
							p[0] = delphi->x[ix];
							p[1] = delphi->y[iy+1]-delphi->hside;
							p[2] = delphi->z[iz+1]-delphi->hside;
						}
						coordVec v(ix,iy+1,iz+1,p,X_DIR);
						vertList->push_back(v);
												
						// mark that this normal should be approximated
						if (computeNormals && providesAnalyticalNormals)
						{
							double* normal = allocateVector<double>(3);
							normal[0]=0;
							normal[1]=0;
							normal[2]=0;
							// approximate the normal by 
							normalsList->push_back(normal);
						}
					}								
				}
	
				// 10
				if (edgeTable[cubeindex] & 1024)
				{				
					vert_offset = -1;
					vert_offset =	intersectionsMatrixAlongZ->at(ix+1,iy+1,iz);
										
					if (vert_offset==-1)
					{					
						double* p = allocateVector<double>(3);
						
						if (accurateTriangulation && isAvailableScalarField)
							vertexInterp(isolevel,positions[2],positions[6],votes[2],votes[6],p);						
						else
						{
							p[0] = delphi->x[ix+1]-delphi->hside;
							p[1] = delphi->y[iy+1]-delphi->hside;
							p[2] = delphi->z[iz];
						}

						coordVec v(ix+1,iy+1,iz,p,Z_DIR);
						vertList->push_back(v);

						// mark that this normal should be approximated
						if (computeNormals && providesAnalyticalNormals)
						{
							double* normal = allocateVector<double>(3);
							normal[0]=0;
							normal[1]=0;
							normal[2]=0;
							// approximate the normal by 
							normalsList->push_back(normal);
						}
					}								
				}
			} // end y
	} // end k 


	deleteMatrix2D<double>(8,positions);		
}

inline double Surface::triangulationKernel(double isolevel,bool revert,int start_z,int end_z,int jump,vector<int*>* localTriList,double* localArea)
{	
	double** positions = allocateMatrix2D<double>(8,3);
	int** triangles = allocateMatrix2D<int>(5,3);
	double** vertPointers = allocateVector<double*>(3);

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	(*localArea) = 0;

	double* gx = delphi->x;
	double* gy = delphi->y;
	double* gz = delphi->z;
	double d_hside = delphi->hside;

	// surprisingly this is the cache friendly set of loops

	for (int k=start_z;k<end_z;k+=jump)	
	{		
		for (int i=0;i<NX;i++)
			for (int j=0;j<NY;j++)			
			{

				// skip the boundary of the grid, there will never be triangles here if
				// the grid is correctly built
				if (i==0 || i==(NX-1)|| j==0 || j==(NY-1) || k==0 || k==(NZ-1))				
					continue;

				// early reject
				//if (!(activeCubes[i][j][k]))
				if (!(read3DVector<bool>(activeCubes,i,j,k,NX,NY,NZ)))
					continue;

				double votes[8];

				// if an accurate description is available each inside/outside flag
				// and edge surface intersections were analytically computed by the ray-tracing routine
				// this allows to perform a marching cubes where each vertex of the mesh is granted to
				// to belong to the surface. 
				// This is an analytical marching cubes				
				if (accurateTriangulation && !isAvailableScalarField)
				{
					if (verticesInsidenessMap[i][j][k])					
						votes[0] = +1;
					else
						votes[0] = -1;
					
					//bool val = read3DVector<bool>(verticesInsidenessMap,i,j,k,NX,NY,NZ);
					//votes[0] = ((val) ? (+1) : (-1));
					
					
					if (verticesInsidenessMap[i+1][j][k])
						votes[1] = +1;
					else
						votes[1] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i+1,j,k,NX,NY,NZ);
					//votes[1] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i+1][j+1][k])
						votes[2] = +1;
					else
						votes[2] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i+1,j+1,k,NX,NY,NZ);
					//votes[2] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i][j+1][k])
						votes[3] = +1;
					else
						votes[3] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i,j+1,k,NX,NY,NZ);
					//votes[3] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i][j][k+1])
						votes[4] = +1;
					else
						votes[4] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i,j,k+1,NX,NY,NZ);
					//votes[4] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i+1][j][k+1])
						votes[5] = +1;
					else
						votes[5] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i+1,j,k+1,NX,NY,NZ);
					//votes[5] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i+1][j+1][k+1])
						votes[6] = +1;
					else
						votes[6] = -1;
					
					//val = read3DVector<bool>(verticesInsidenessMap,i+1,j+1,k+1,NX,NY,NZ);
					//votes[6] = ((val) ? (+1) : (-1));

					
					if (verticesInsidenessMap[i][j+1][k+1])
						votes[7] = +1;
					else
						votes[7] = -1;
					

					//val = read3DVector<bool>(verticesInsidenessMap,i,j+1,k+1,NX,NY,NZ);
					//votes[7] = ((val) ? (+1) : (-1));
				}
				// Classical marching cube interpolating scalar field values
				else if (accurateTriangulation && isAvailableScalarField)
				{
					votes[0] = scalarField[i][j][k];
					votes[1] = scalarField[i+1][j][k];
					votes[2] = scalarField[i+1][j+1][k];
					votes[3] = scalarField[i][j+1][k];
					votes[4] = scalarField[i][j][k+1];
					votes[5] = scalarField[i+1][j][k+1];
					votes[6] = scalarField[i+1][j+1][k+1];
					votes[7] = scalarField[i][j+1][k+1];
				}
				
				int numTriangles;

				double xg = gx[i];
				double yg = gy[j];
				double zg = gz[k];

				double table[6];
				table[0] = xg - d_hside;
				table[1] = xg + d_hside;
				table[2] = yg - d_hside;
				table[3] = yg + d_hside;
				table[4] = zg - d_hside;
				table[5] = zg + d_hside;
				
				if (!accurateTriangulation)
					for (int vertInd =0;vertInd<8;vertInd++)
					{
						votes[vertInd] = getInsidness(i,j,k,vertInd);
						positions[vertInd][0]= table[mulTable[vertInd][0]];
						positions[vertInd][1]= table[mulTable[vertInd][1]];
						positions[vertInd][2]= table[mulTable[vertInd][2]];
					}
				else
					for (int vertInd =0;vertInd<8;vertInd++)
					{
						positions[vertInd][0]= table[mulTable[vertInd][0]];
						positions[vertInd][1]= table[mulTable[vertInd][1]];
						positions[vertInd][2]= table[mulTable[vertInd][2]];
					}
				
				// now that the votes and positions are obtained 
				// triangulate the cell by marching cube rule.											
				numTriangles = getTriangles(votes,positions,isolevel,triangles,i,j,k,NX,NY,NZ);				
				
				//printf("\n Num tri %d",numTriangles);
				// save triangles
				for (int kk = 0;kk<numTriangles;kk++)
				{
					int tri_temp[3];					

					for (int ll = 0;ll<3;ll++)
					{										
						tri_temp[ll] = triangles[kk][ll];
						vertPointers[ll]=vertList[tri_temp[ll]];
					}
							
					double a =0,b=0,c=0;
					DIST(a,vertPointers[0],vertPointers[1]);
					DIST(b,vertPointers[0],vertPointers[2]);
					DIST(c,vertPointers[1],vertPointers[2]);
					double ttt = (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c);
					if (ttt>0)						
						(*localArea) += 0.25*sqrt(ttt);

					int* tri = allocateVector<int>(3);
					if (!isAvailableScalarField)
					{
						tri[2]=tri_temp[0];
						tri[1]=tri_temp[1];
						tri[0]=tri_temp[2];						
					}
					else
					{
						tri[0]=tri_temp[0];
						tri[1]=tri_temp[1];
						tri[2]=tri_temp[2];						
					}
					localTriList->push_back(tri);
				}				
				
			}
			
	}
	// end of MC kernel
	//cout << endl <<"I am finished " << start_z << " " << end_z;
	//cout.flush();

	deleteMatrix2D<double>(8,positions);
	deleteMatrix2D<int>(5,triangles);
	deleteVector<double*>(vertPointers);

	return (*localArea);
}

double Surface::triangulateSurface(double isolevel,const char* fileName,bool revert)
{
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	Timer chrono;
	chrono.start();

	if (verticesInsidenessMap==NULL && accurateTriangulation && !isAvailableScalarField)
	{
		cout << endl << WARN << "Cannot triangulate without inside/out info for grid points";
		return 0;
	}

	int num_thd = 1;

	if (num_cores<=0)
	{
		#ifdef ENABLE_BOOST_THREADS
			num_thd = (int)boost::thread::hardware_concurrency();		
		#endif
		cout << endl << INFO << "Setting " << num_thd << " threads";
	}
	else
	{
		cout << endl << INFO << "User selected num threads " << num_cores;
		num_thd = num_cores;
	}

	#ifdef ENABLE_BOOST_THREADS
			boost::thread_group thdGroup; 
	#endif
	
	// in parallel complete missing vertices and mark active cubes to make faster the second pass
	// in case of classical MC, all vertices are computed on this pass
	double* area = NULL;
	vector < vector<int*>* > localTri;
	vector < vector<coordVec>* > localVert;
	vector < vector<double*>* > localNormals;

	area = allocateVector<double>(num_thd);

	localTri.reserve(num_thd);	
	localVert.reserve(num_thd);	
	localNormals.reserve(num_thd);

	for (int i=0;i<num_thd;i++)
	{
		localTri.push_back(new vector<int*>());
		localVert.push_back(new vector<coordVec>());
		localNormals.push_back(new vector<double*>());
	}

	//activeCubes = allocateMatrix3D<bool>(delphi->nx,delphi->ny,delphi->nz);
	if (activeCubes!=NULL)
		deleteVector<bool>(activeCubes);

	activeCubes = allocateVector<bool>(delphi->nx*delphi->ny*delphi->nz);
								
	///////////// starts parallel vertices retrieval ///////////////////////
	// this is time consuming only if the scalar field is available
	// differently, in case of analytical intersections, this step will be very fast
	// and only used to repair variations done by the cavity filling or by rays that missed the target

	cout << endl << INFO << "Generating MC vertices...";
	cout.flush();

	// load balanced and cache friendly thread dispatch
	for (int j=0;j<num_thd;j++)
	{
		#ifdef ENABLE_BOOST_THREADS
			thdGroup.create_thread(boost::bind(&Surface::getVertices,this,isolevel,j,NZ,num_thd,localVert[j],localNormals[j]));
		#else
			getVertices(isolevel,0,NZ,1,localVert[0],localNormals[0]);
		#endif						
	}
	///////////////////////////////////////////////////////////////////////////////
		
	#ifdef ENABLE_BOOST_THREADS
		// join
		thdGroup.join_all();
	#endif

	int addedVertices = 0;
//	int dup=0;
		for (int j=0;j<num_thd;j++)
		{
			vector<coordVec>* lv = localVert[j];
			vector<coordVec>::iterator it2;

			addedVertices+=(int)localVert[j]->size();

			vector<double*>* ln;
			vector<double*>::iterator it3;

			if (computeNormals && providesAnalyticalNormals)
			{
				ln = localNormals[j];
				it3 = ln->begin();
			}

//			set<string> brutal;			

			for (it2=lv->begin();it2!=lv->end();it2++)
			{
				coordVec& v = *it2;
				vertList.push_back(v.vec);
				int ind = (int)vertList.size()-1;
				if (v.dir==X_DIR)
					intersectionsMatrixAlongX->set(v.ix,v.iy,v.iz,ind);				
				else if (v.dir==Y_DIR)
					intersectionsMatrixAlongY->set(v.ix,v.iy,v.iz,ind);			
				else
					intersectionsMatrixAlongZ->set(v.ix,v.iy,v.iz,ind);

//				char buff[100];
//				sprintf(buff,"%d_%d_%d_%d",v.ix,v.iy,v.iz,v.dir);
//				set<string>::iterator it = brutal.find(string(buff));
//				if (it!=brutal.end())
//					dup++;
//				else
//					brutal.insert(string(buff));

				if (computeNormals && providesAnalyticalNormals)
				{
					double* vn = *it3;

					if (v.dir==X_DIR)
						normalsMatrixAlongX->set(v.ix,v.iy,v.iz,ind);
					else if (v.dir==Y_DIR)
						normalsMatrixAlongY->set(v.ix,v.iy,v.iz,ind);
					else
						normalsMatrixAlongZ->set(v.ix,v.iy,v.iz,ind);					
					
					normalsList.push_back(vn);
					it3++;
				}
			}
		}
	
	//cout << endl << "duplicated "<<dup;
	cout << "ok!";
	cout.flush();

	cout << endl << INFO << "MC added "<< addedVertices << " non analytical vertices";
	/*
	int orphans2 = 0;
	int orphans1 = 0;

	// search for orphan cubes
	for (int ix=0;ix<delphi->nx;ix++)
		for (int iy=0;iy<delphi->ny;iy++)
			for (int iz=0;iz<delphi->nz;iz++)
			{				
				int count = 0;

				int test = intersectionsMatrixAlongX->at(ix,iy,iz);						
				if (test!=-1) count++;

				test = intersectionsMatrixAlongY->at(ix+1,iy,iz);
				if (test!=-1) count++;

				test = intersectionsMatrixAlongX->at(ix,iy+1,iz);
				if (test!=-1) count++;

				test = intersectionsMatrixAlongY->at(ix,iy,iz);
				if (test!=-1) count++;
				
				test = intersectionsMatrixAlongX->at(ix,iy,iz+1);
				if (test!=-1) count++;
	
				test = intersectionsMatrixAlongY->at(ix+1,iy,iz+1);
				if (test!=-1) count++;

				test	=	intersectionsMatrixAlongX->at(ix,iy+1,iz+1);
				if (test!=-1) count++;

				test = intersectionsMatrixAlongY->at(ix,iy,iz+1);
				if (test!=-1) count++;

				test = intersectionsMatrixAlongZ->at(ix,iy,iz);
				if (test!=-1) count++;
	
				test = intersectionsMatrixAlongZ->at(ix+1,iy,iz);
				if (test!=-1) count++;
	
				test =	intersectionsMatrixAlongZ->at(ix+1,iy+1,iz);			
				if (test!=-1) count++;
		
				test = intersectionsMatrixAlongZ->at(ix,iy+1,iz);
				if (test!=-1) count++;

				if (count==2)
					orphans2++;

				if (count==1)
					orphans1++;

			}

	cout << endl << "orphans1 " << orphans1;
	cout << endl << "orphans2 " << orphans2;
	*/

	// trivial split
	int chunk = NZ/num_thd;
	int rem = NZ%num_thd;

	// setup split
	int start = 0;
	int stop  = 0;

	////////////////////////////// generate triangles /////////////////////////
	// all vertices are computed, stored and uniquely indexed, now get triangles
	for (int j=0;j<num_thd;j++)
	{
		// only octrees version
		#ifdef ENABLE_BOOST_THREADS		
			thdGroup.create_thread(boost::bind(&Surface::triangulationKernel,this,isolevel,revert,j,NZ,num_thd,localTri[j],&(area[j])));
		#else
			triangulationKernel(isolevel,revert,0,NZ,1,localTri[0],&(area[0]));						
		#endif					
	}
		
	totalSurfaceArea = 0;

	#ifdef ENABLE_BOOST_THREADS
		// join
		thdGroup.join_all();

		cout << endl << INFO << "Triangles done";
		cout.flush();

		// reduce
		for (int j=0;j<num_thd;j++)
		{
			totalSurfaceArea+=area[j];				
			vector<int*>* lt = localTri[j];						
			vector<int*>::iterator it;
			for (it=lt->begin();it!=lt->end();it++)			
				triList.push_back(*it);
		}
	#else			
		totalSurfaceArea = area[0];
		vector<int*>* lt = localTri[0];
		vector<int*>::iterator it;
		for (it=lt->begin();it!=lt->end();it++)			
			triList.push_back(*it);
	#endif

	//deleteMatrix3D<bool>(delphi->nx,delphi->ny,activeCubes);
	deleteVector<bool>(activeCubes);

	for (int i=0;i<num_thd;i++)
	{
		delete localTri[i];
		delete localVert[i];
		delete localNormals[i];
	}

	double duration = chrono.stop();
	cout << endl << INFO << "MC time is " << duration << " [s]";

	cout << endl << INFO << "Total, grid conformant, surface area is " << totalSurfaceArea << " [A^2]";
	
	int numVertexes = (int)vertList.size();
	int numTriangles = (int)triList.size();

	// DEBUG
	/*
	int dang =0;
	for (int jj=0;jj<numVertexes;jj++)
	{
		bool found = false;
		for (int tt=0;tt<numTriangles;tt++)
		{
			if (triList[tt][0]==jj || triList[tt][1]==jj || triList[tt][2]==jj)
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

	cout << endl << INFO << "Number of vertices " << numVertexes << " number of triangles " << numTriangles;

	if (vertexAtomsMapFlag && vertexAtomsMap!=NULL)
		deleteVector<int>(vertexAtomsMap);

	/// check atoms flag
	if (vertexAtomsMapFlag)
	{	
		vertexAtomsMap = allocateVector<int>(numVertexes);
		buildAtomsMap();
		cout << endl << INFO << "Connecting vertices to atoms..";
		cout.flush();
		for (int i=0;i<numVertexes;i++)			
		{			
			vdwAccessible(vertList[i],vertexAtomsMap[i]);
			if (vertexAtomsMap[i]==-1)
			{
				cout << endl << WARN << "Cannot detect nearest atom for vertex " << i;
			}
		}
		cout << "ok!";
		cout.flush();
	}

	vector<int> appNormals;
	// check if it has to approximate normals or they are already available
	if (computeNormals && providesAnalyticalNormals)
	{
		// serch for normals to be approximated
		if (normalsList.size()!=0)
		{
			vector<double*>::iterator it;
			bool found = false;
			int ind=0;
			for (it = normalsList.begin();it!=normalsList.end();it++)
			{
				double* p = *it;
				if (p[0]==0 && p[1]==0 && p[2]==0)
					appNormals.push_back(ind);				
				ind++;
			}
		}
	}

	if (computeNormals && providesAnalyticalNormals)
	{	
		//cout << endl << appNormals.size();
		if (appNormals.size()!=0)
		{
			cout << endl << INFO << "Some analytical normals will be approximated...";		
			approximateNormals(appNormals,true);
		}
	}	
	else if (computeNormals && !providesAnalyticalNormals)
	{
		cout << endl << INFO << "Approximating vertices normals by triangulation...";
		approximateNormals(appNormals,false);
	}

	int format = deduceFormat();
	bool f = saveMesh(format,revert,fileName,vertList,triList,normalsList);

	if (vertexAtomsMapFlag)
		disposeAtomsMap();

	if (!f)
	{
		cout << endl << ERR << "Errors in saving the mesh!";
	}

	return totalSurfaceArea;
}

void Surface::approximateNormals(vector<int>& appNormals,bool doOnlyList)
{
	int nv = (int)vertList.size();
	int nt = (int)triList.size();

	vector<int>** vertexTrianglesList = new vector<int>*[nv];

	for (int i=0;i<nv;i++)
		vertexTrianglesList[i] = new vector<int>();

	double** planes = allocateMatrix2D<double>(nt,3);
	 
	for (int i=0;i<nt;i++)
	{
		vertexTrianglesList[triList[i][0]]->push_back(i);
		vertexTrianglesList[triList[i][1]]->push_back(i);
		vertexTrianglesList[triList[i][2]]->push_back(i);
	}
     
	// compute triangle planes		
	double w[4],p1[3],p2[3],p3[3];
	
	// compute plane normals from triangles 
	for(int it=0;it<nt;it++)
	{
	  ASSIGN(p1,vertList[triList[it][0]]);
	  ASSIGN(p2,vertList[triList[it][1]]);
	  ASSIGN(p3,vertList[triList[it][2]]);	  
	  plane3points(p1,p2,p3,w,false);
	  planes[it][0]=w[0];
	  planes[it][1]=w[1];
	  planes[it][2]=w[2];
	}

	if (!doOnlyList)
	{				
		normalsList.reserve(nv);
				
		double* meanNormal;
		for (int iv=0;iv<nv;iv++)
		{
			meanNormal = allocateVector<double>(3);			
			normalsList.push_back(meanNormal);
			meanNormal[0]=0;
			meanNormal[1]=0;
			meanNormal[2]=0;
		    
			vector<int>::iterator it; 

			for (it = vertexTrianglesList[iv]->begin();it!=vertexTrianglesList[iv]->end();it++)
			{
				int triangleID = *it;
				ADD(meanNormal,meanNormal,planes[triangleID]);
			}

			double vlen = (double)vertexTrianglesList[iv]->size();
			
			meanNormal[0]=meanNormal[0]/vlen;
			meanNormal[1]=meanNormal[1]/vlen;
			meanNormal[2]=meanNormal[2]/vlen;
			
			double norm;
			NORMALIZE(meanNormal,norm);
		}
	}
	else
	{
		double* meanNormal;
		vector<int>::iterator it;
		for (it=appNormals.begin();it!=appNormals.end();it++)		
		{
			int iv = (*it);
			meanNormal = normalsList[iv];
			meanNormal[0]=0;
			meanNormal[1]=0;
			meanNormal[2]=0;
		    
			vector<int>::iterator it; 

			for (it = vertexTrianglesList[iv]->begin();it!=vertexTrianglesList[iv]->end();it++)
			{
				int triangleID = *it;
				ADD(meanNormal,meanNormal,planes[triangleID]);
			}

			double vlen = (double)vertexTrianglesList[iv]->size();
			
			meanNormal[0]=meanNormal[0]/vlen;
			meanNormal[1]=meanNormal[1]/vlen;
			meanNormal[2]=meanNormal[2]/vlen;
			
			double norm;
			NORMALIZE(meanNormal,norm);
		}	
	}

	for (int i=0;i<nv;i++)
		delete vertexTrianglesList[i];
	delete[] vertexTrianglesList;
		
	deleteMatrix2D<double>(nt,planes);
}
bool Surface::saveMesh(int format,bool revert,const char* fileName,vector<double*>& vertList,vector<int*>& triList, vector<double*>& normalsList)
{
	int numVertexes = (int)vertList.size();
	int numTriangles = (int)triList.size();
	char fullName[100];

	if (format==OFF || format==OFF_A || format==OFF_N || format==OFF_N_A)
	{
		// save all in OFF format
		FILE* fp;
		sprintf(fullName,"%s.off",fileName);	
		fp = fopen(fullName,"w");

		if (fp==NULL)
		{
			cout << endl << WARN << "Cannot write file " << fileName;
			return false;
		}

		if (format==OFF_A)
		{
			if (vertexAtomsMap==NULL)
			{
				cout << endl << ERR << "Cannot save in OFF+A format if nearest atom info is not available";
				fclose(fp);
				return false;
			}
			fprintf(fp,"OFF+A\n");
			cout << endl << INFO << "Writing triangulated surface in OFF+A file format in " << fileName << "...";
		}
		else if (format==OFF_N)
		{			
			if (normalsList.size()==0)
			{
				cout << endl << ERR << "Cannot save in OFF+N format if normals are not available";
				fclose(fp);
				return false;
			}
			fprintf(fp,"OFF+N\n");
			cout << endl << INFO << "Writing triangulated surface in OFF+N file format in " << fileName << "...";
		}
		else if (format==OFF_N_A)
		{		
			if (normalsList.size()==0)
			{
				cout << endl << ERR << "Cannot save in OFF+N+A format if normals are not available";
				fclose(fp);
				return false;
			}

			if (vertexAtomsMap==NULL)
			{
				cout << endl << ERR << "Cannot save in OFF+N+A format if nearest atom info is not available";
				fclose(fp);
				return false;
			}

			fprintf(fp,"OFF+N+A\n");
			cout << endl << INFO << "Writing triangulated surface in OFF+N+A file format in " << fileName << "...";
		}
		else
		{
			fprintf(fp,"OFF\n");
			cout << endl << INFO << "Writing triangulated surface in OFF file format in " << fileName << "...";
		}

		cout.flush();

		time_t pt;
		time(&pt);    
		fprintf(fp,"# File created by %s version %s date %s\n",PROGNAME,VERSION,ctime(&pt));	
		fprintf(fp,"%d %d 0\n",numVertexes,numTriangles);

		for (int i=0;i<numVertexes;i++)		
		{
			if (format==OFF_A)
				fprintf(fp,"%.3f %.3f %.3f %d\n",vertList[i][0],vertList[i][1],vertList[i][2],vertexAtomsMap[i]);
			else if (format==OFF_N)
				fprintf(fp,"%.3f %.3f %.3f %.3f %.3f %.3f\n",vertList[i][0],vertList[i][1],vertList[i][2],normalsList[i][0],normalsList[i][1],normalsList[i][2]);
			else if (format==OFF_N_A)
				fprintf(fp,"%.3f %.3f %.3f %.3f %.3f %.3f %d\n",vertList[i][0],vertList[i][1],vertList[i][2],normalsList[i][0],normalsList[i][1],normalsList[i][2],vertexAtomsMap[i]);
			else
				fprintf(fp,"%.3f %.3f %.3f\n",vertList[i][0],vertList[i][1],vertList[i][2]);
		}

		for (int i=0;i<numTriangles;i++)
			if (!revert)
				fprintf(fp,"3 %d %d %d\n",triList[i][0],triList[i][1],triList[i][2]);
			else
				fprintf(fp,"3 %d %d %d\n",triList[i][2],triList[i][1],triList[i][0]);

		fclose(fp);
	}
	else if (format==MSMS || format == MSMS_NO_A)
	{
		if (normalsList.size()==0)
		{
			cout << endl << ERR << "Cannot save in MSMS format if normals are not available";		
			return false;
		}

		if (format==MSMS)
		{
			cout << endl << INFO << "Saving in MSMS format, no patch info...";
			if (vertexAtomsMap==NULL)
			{
				cout << endl << ERR << "Cannot save in MSMS format if nearest atom info is not available";
				return false;
			}
		}
		else if (format == MSMS_NO_A)
			cout << endl << INFO << "Saving in MSMS format, no patch info, no nearest atom..";

		FILE *fp1,*fp2;
		sprintf(fullName,"%s.face",fileName);
		fp1 = fopen(fullName,"w"); 
		sprintf(fullName,"%s.vert",fileName);
		fp2 = fopen(fullName,"w");

		if (fp1==NULL || fp2==NULL)
		{
			cout << endl << ERR << "Error in writing MSMS files";
			return false;
		}

		time_t pt;
		time(&pt);    
		
		fprintf(fp1,"# File created by %s version %s date %s",PROGNAME,VERSION,ctime(&pt));
		fprintf(fp1,"#faces\n");
		fprintf(fp1,"%d\n",numTriangles);

		for (int i=0;i<numTriangles;i++)
			if (!revert)
				fprintf(fp1,"%6d %6d %6d %2d %6d\n",triList[i][0]+1,triList[i][1]+1,triList[i][2]+1,1,1);
			else
				fprintf(fp1,"%6d %6d %6d %2d %6d\n",triList[i][2]+1,triList[i][1]+1,triList[i][0]+1,1,1);

		fclose(fp1);

		fprintf(fp2,"# File created by %s version %s date %s",PROGNAME,VERSION,ctime(&pt));
		fprintf(fp2,"#vertex\n");
		fprintf(fp2,"%d\n",numVertexes);
	
		for (int i=0;i<numVertexes;i++)
			if (format==MSMS)
				fprintf(fp2,"%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d \n",vertList[i][0],vertList[i][1],vertList[i][2],normalsList[i][0],normalsList[i][1],normalsList[i][2],0,vertexAtomsMap[i]+1,0);
			else
				fprintf(fp2,"%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d \n",vertList[i][0],vertList[i][1],vertList[i][2],normalsList[i][0],normalsList[i][1],normalsList[i][2],0,1,0);
		fclose(fp2);
	}
	return true;
}

bool Surface::difference(Surface* surf)
{
	int NX,NY,NZ;
	NX = delphi->nx;
	NY = delphi->ny;
	NZ = delphi->nz;

	double* X = delphi->x;
	double* Y = delphi->y;
	double* Z = delphi->z;

	// S1 is the fat probe
	// S2 is the regular probe

	int inside1 = inside;
	int inside2 = surf->inside;

	// grid consistency checks
	if (NX != surf->delphi->nx ||
		NY != surf->delphi->ny ||
		NZ != surf->delphi->nz )
	{
		return false;
	}

	if (delphi->buildEpsMap!=surf->delphi->buildEpsMap)
	{
		return false;
	}

	if (delphi->buildStatus!=surf->delphi->buildStatus)
	{
		return false;
	}

	// apply difference rule to epsmap, status map, idebmap, insideness map

	// skip epsmap usually not used
/*
	// epsmap
	if (delphi->buildEpsMap)
	{				
		for (int i=0;i<NX;i++)
			for (int j=0;j<NY;j++)
				for (int k=0;k<NZ;k++)
				{
					for (int l=0;l<3;l++)
					{
						int stat1 = delphi->EPSMAP(i,j,k,l,NX,NY,NZ);
						int stat2 = surf->delphi->EPSMAP(i,j,k,l,NX,NY,NZ);
						// if S1 is not out and and S2 is out then that's the pocket
						// the new pocket is marked as outside
						// This is done such that the new objects can be detected as cavities
						// by the cavity detector in a further step such that they can be
						// volume filtered						
						if ((stat1!=0 && stat2==0))
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
		if (delphi->cavitiesVec !=NULL)
		{
			vector< vector < int* >* >::iterator it;
			for (it=delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++)
			{
				vector < int* >* inner = (*it);
				vector < int* >::iterator it2;
				for (it2=inner->begin();it2!=inner->end();it2++)
					free((*it2));
				delete inner;
			}
			delete delphi->cavitiesVec;
		}

		delphi->cavitiesVec = new vector < vector <int* > * >();

		// add the unique 'cavity'.
		// At this stage the 'cavity' is simply the result of the difference map
		// such that it will be filtered by the digital Connolly filter
		vector< int* >* vv = new vector<int*>();

		delphi->cavitiesVec->push_back(vv);

		// clean memory
		delphi->cavitiesSize.clear();	
		delphi->cavitiesFlag.clear();

		// keep it
		delphi->cavitiesFlag.push_back(false);

		int countCubes=0;

		for (int k=0;k<NZ;k++)		
			for (int j=0;j<NY;j++)
				for (int i=0;i<NX;i++)
				{
					// TODO idebmap

					//int stat1 = delphi->status[i][j][k];
					//int stat1 = delphi->STATUSMAP(i,j,k,NX,NY);
					int stat1 = read3DVector<short>(delphi->status,i,j,k,NX,NY,NZ);

					//int stat2 = surf->delphi->status[i][j][k];
					//int stat2 = surf->delphi->STATUSMAP(i,j,k,NX,NY);
					int stat2 = read3DVector<short>(surf->delphi->status,i,j,k,NX,NY,NZ);
										
					// if S1 is not out and and S2 is out then that's the pocket					
					if ((stat1!=STATUS_POINT_TEMPORARY_OUT && stat1!=STATUS_POINT_OUT) && (stat2==STATUS_POINT_TEMPORARY_OUT || stat2==STATUS_POINT_OUT))
					{

						// mark as outside temporary, such that cavity detection can work on it						
						//delphi->STATUSMAP(i,j,k,NX,NY)=STATUS_POINT_TEMPORARY_OUT;

						// mark as first cavity such that it can be directly filtered
						//delphi->STATUSMAP(i,j,k,NX,NY)=STATUS_FIRST_CAV;
						write3DVector<short>(delphi->status,STATUS_FIRST_CAV,i,j,k,NX,NY,NZ);
						//cout << endl << "marking";

						countCubes++;

						int* v = allocateVector<int>(3);
						v[0]=i;
						v[1]=j;
						v[2]=k;
						// store the point as a cavity
						vv->push_back(v);
						
						// switch surrounding points
						
						if (accurateTriangulation)
						{
							verticesInsidenessMap[i][j][k] = true;
							verticesInsidenessMap[i+1][j][k] = true;
							verticesInsidenessMap[i][j+1][k] = true;
							verticesInsidenessMap[i+1][j+1][k] = true;
							verticesInsidenessMap[i][j][k+1] = true;
							verticesInsidenessMap[i+1][j][k+1] = true;
							verticesInsidenessMap[i][j+1][k+1] = true;
							verticesInsidenessMap[i+1][j+1][k+1] = true;

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
					}
					else
					{
						//delphi->status[i][j][k]=STATUS_POINT_INSIDE;
						//delphi->STATUSMAP(i,j,k,NX,NY)=STATUS_POINT_INSIDE;
						write3DVector<short>(delphi->status,STATUS_POINT_INSIDE,i,j,k,NX,NY,NZ);

						if (accurateTriangulation)
						{
							verticesInsidenessMap[i][j][k] = false;
							//write3DVector<bool>(verticesInsidenessMap,false,i,j,k,NX,NY,NZ);
							
							if ((i+1)<NX) 
								verticesInsidenessMap[i+1][j][k] = false;
								//write3DVector<bool>(verticesInsidenessMap,false,i+1,j,k,NX,NY,NZ);

							if ((j+1)<NY)
								verticesInsidenessMap[i][j+1][k] = false;
								//write3DVector<bool>(verticesInsidenessMap,false,i,j+1,k,NX,NY,NZ);

							if ((i+1)<NX && (j+1)<NY)
								verticesInsidenessMap[i+1][j+1][k] = false;
								//write3DVector<bool>(verticesInsidenessMap,false,i+1,j+1,k,NX,NY,NZ);

							if ((k+1)<NZ)
								verticesInsidenessMap[i][j][k+1] = false;
								//write3DVector<bool>(verticesInsidenessMap,false,i,j,k+1,NX,NY,NZ);

							if ((i+1)<NX && (k+1)<NZ)
								verticesInsidenessMap[i+1][j][k+1] = false;
								//write3DVector<bool>(verticesInsidenessMap,false,i+1,j,k+1,NX,NY,NZ);

							if ((j+1)<NY && (k+1)<NZ)
								verticesInsidenessMap[i][j+1][k+1] = false;
								//write3DVector<bool>(verticesInsidenessMap,false,i,j+1,k+1,NX,NY,NZ);

							if ((i+1)<NX && (j+1)<NY && (k+1)<NZ)
								verticesInsidenessMap[i+1][j+1][k+1] = false;
								//write3DVector<bool>(verticesInsidenessMap,false,i+1,j+1,k+1,NX,NY,NZ);						
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
	if (vertList.size()==0)
	{
		cout << endl << WARN << "Triangulation not available, cannot build set of balls approximation!";
		cout << endl << REMARK << "Please enable triangulation by 'Triangulation = true' ";
		return;
	}
	
	if (delphi->status==NULL)
	{
		cout << endl << WARN << "Status map not computed, cannot build set of balls approximation!";
		cout << endl << REMARK << "Please enable status map by 'Build_status_map = true'";
		return;
	}

	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	#ifndef ENABLE_CGAL 	
	cout << endl << WARN << "CGAL is required by tri2Balls function!";
	cout << endl << REMARK << "Recompile enabling CGAL support";
	return;
	#else

	vector<_Weighted_point> l;
	float x,y,z;
	float max_x=-1e6,min_x=1e6,max_y=-1e6,min_y=1e6,max_z=-1e6,min_z=1e6;

	l.reserve(vertList.size()+6);

	//cout << endl << INFO << "Converting triangulation to set of balls..";
			
	for (unsigned int i=0;i<vertList.size();i++)
	{
		x = (float)vertList[i][0];
		y = (float)vertList[i][1];
		z = (float)vertList[i][2];
				
		l.push_back(_Weighted_point(_Point3(x,y,z),0.f,i)); 
		max_x = max(max_x,x);
		max_y = max(max_y,y);
		max_z = max(max_z,z);

		min_x = std::min(min_x,x);
		min_y = std::min(min_y,y);
		min_z = std::min(min_z,z);	
	}
	
	// Regular Triangulation object 
	_Rt rT;

	float mid_x = (max_x+min_x)/2.f;
	float mid_y = (max_y+min_y)/2.f;
	float mid_z = (max_z+min_z)/2.f;

	min_x -= fabs(mid_x-min_x)*2;
	min_y -= fabs(mid_y-min_y)*2;
	min_z -= fabs(mid_z-min_z)*2;

	max_x += fabs(mid_x-max_x)*2;
	max_y += fabs(mid_y-max_y)*2;
	max_z += fabs(mid_z-max_z)*2;

	// add bounding box using -1 weights to let easy detection of these virtual points
	l.push_back(_Weighted_point(_Point3(min_x,mid_y,mid_z),-1,(const int)vertList.size()));
	l.push_back(_Weighted_point(_Point3(max_x,mid_y,mid_z),-1,(const int)vertList.size()+1));
	l.push_back(_Weighted_point(_Point3(mid_x,min_y,mid_z),-1,(const int)vertList.size()+2));
	l.push_back(_Weighted_point(_Point3(mid_x,max_y,mid_z),-1,(const int)vertList.size()+3));
	l.push_back(_Weighted_point(_Point3(mid_x,mid_y,min_z),-1,(const int)vertList.size()+4));
	l.push_back(_Weighted_point(_Point3(mid_x,mid_y,max_z),-1,(const int)vertList.size()+5));

	//cout << endl << INFO << "Computing triangulation....";

	rT.insert (l.begin(), l.end());

	assert( rT.is_valid() );
	assert( rT.dimension() == 3 );

	//cout << "ok!";
	//cout << endl << INFO << "Computing voronoi points....";

	_Finite_Cells_Iterator fcit = rT.finite_cells_begin(); 

	for (;fcit!=rT.finite_cells_end();fcit++)    
	{
		const _Point3& p = rT.geom_traits().construct_weighted_circumcenter_3_object()(fcit->vertex(0)->point(),fcit->vertex(1)->point(),fcit->vertex(2)->point(),fcit->vertex(3)->point());
		VorPoint*& mc = (VorPoint*&)fcit->info();
		mc = new VorPoint();
		mc->vor[0] = p.x();
		mc->vor[1] = p.y();
		mc->vor[2] = p.z();
	}

	//cout << "ok!";
	//cout << endl << INFO << "Collecting polar balls....";
	
	vector<_Cell_handle> cells;
	cells.reserve(1000);

	vector<double*> polarBalls;
	vector<double> polarDist;

	_Finite_Vertex_Iterator fvit = rT.finite_vertices_begin();
	
	for(;fvit!=rT.finite_vertices_end();fvit++)
	{	  
		const _Weight& w0 = fvit->point().weight();
		
		// skip bounding box points
		if ((w0)==-1)
			continue;
	  
		int numPoints = 0;

		if (rT.is_infinite(fvit->cell()))
			continue;

		// current point id around which we are moving
		const int& currentId = fvit->point().index();
		
		cells.clear();
		rT.incident_cells(fvit,std::back_inserter(cells));

		bool infiniteCell = false;

		// start check if all tetraedra are feasible
		std::vector<_Cell_handle>::iterator it;
		for (it=cells.begin();it!=cells.end();it++)
		{
			if (rT.is_infinite((*it)))
			{
				infiniteCell = true;
				break;
			}
		}
	 
		if (infiniteCell)
			continue;

		double ref[3],winner[3],maxDist=0,dist;
		VorPoint* winner_vp;
		ref[0]=vertList[currentId][0];
		ref[1]=vertList[currentId][1];
		ref[2]=vertList[currentId][2];
		
		// compute the internal farthest ball (the inner polar ball)
		for (it=cells.begin();it!=cells.end();it++)
		{
			const VorPoint* vp = (*it)->info();
			DIST2(dist,vp->vor,ref);

			// check it is internal by using the status map
			int ix = (int)rintp((vp->vor[0]-delphi->xmin)/delphi->side);
			int iy = (int)rintp((vp->vor[1]-delphi->ymin)/delphi->side);
			int iz = (int)rintp((vp->vor[2]-delphi->zmin)/delphi->side);

			if (ix<0 || iy<0 || iz<0 || ix>=delphi->nx || iy>=delphi->ny || iz>=delphi->nz)
				continue;

			int  count = 0;

			/*
			if ((delphi->STATUSMAP(ix,iy,iz,NX,NY)==STATUS_POINT_INSIDE) &&
				((ix+1)<NX && delphi->STATUSMAP((ix+1),(iy),(iz),NX,NY)==STATUS_POINT_INSIDE) &&
				((iy+1)<NY && delphi->STATUSMAP(ix,(iy+1),iz,NX,NY)==STATUS_POINT_INSIDE) &&
				((iz+1)<NZ && delphi->STATUSMAP(ix,iy,(iz+1),NX,NY)==STATUS_POINT_INSIDE) &&
				((ix-1)>=0 && delphi->STATUSMAP((ix-1),iy,iz,NX,NY)==STATUS_POINT_INSIDE) &&
				((iy-1)>=0 && delphi->STATUSMAP(ix,(iy-1),iz,NX,NY)==STATUS_POINT_INSIDE) &&
				((iz-1)>=0 && delphi->STATUSMAP(ix,iy,(iz-1),NX,NY)==STATUS_POINT_INSIDE) )
			*/
			if ((read3DVector<short>(delphi->status,ix,iy,iz,NX,NY,NZ)==STATUS_POINT_INSIDE) &&
				((ix+1)<NX && read3DVector<short>(delphi->status,ix+1,iy,iz,NX,NY,NZ)==STATUS_POINT_INSIDE) &&
				((iy+1)<NY && read3DVector<short>(delphi->status,ix,iy+1,iz,NX,NY,NZ)==STATUS_POINT_INSIDE) &&
				((iz+1)<NZ && read3DVector<short>(delphi->status,ix,iy,iz+1,NX,NY,NZ)==STATUS_POINT_INSIDE) &&
				((ix-1)>=0 && read3DVector<short>(delphi->status,ix-1,iy,iz,NX,NY,NZ)==STATUS_POINT_INSIDE) &&
				((iy-1)>=0 && read3DVector<short>(delphi->status,ix,iy-1,iz,NX,NY,NZ)==STATUS_POINT_INSIDE) &&
				((iz-1)>=0 && read3DVector<short>(delphi->status,ix,iy,iz-1,NX,NY,NZ)==STATUS_POINT_INSIDE) )
			{
				if (dist>maxDist)
				{
					maxDist = dist;
					ASSIGN(winner,vp->vor);
					winner_vp = (VorPoint*)vp;
				}
			}
		}

		if (maxDist==0)
			continue;

		if (winner_vp->visited==true)
			continue;
		else
			winner_vp->visited = true;

		double* p = allocateVector<double>(3);
		polarBalls.push_back(p);
		ASSIGN(p,winner);
		polarDist.push_back(sqrt(maxDist));
	}
  
	//cout << "ok!";

	FILE* fp = fopen("polar.txt","w");
	FILE* fp2 = fopen("polar.pqr","w");
	for (unsigned int i=0;i<polarBalls.size();i++)	
	{
		fprintf(fp,"%f %f %f %f \n",polarBalls[i][0],polarBalls[i][1],polarBalls[i][2],polarDist[i]);	
		fprintf(fp2,"ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f%8.4f%8.4f\n",i%100000,"C","NSR",'A',1,polarBalls[i][0],polarBalls[i][1],polarBalls[i][2],0.,polarDist[i]);
		deleteVector<double>(polarBalls[i]);
	}

	fclose(fp);
	fclose(fp2);

	#endif
}

void Surface::smoothSurface(const char* fn,bool revert)
{
	#define MAX_NEIGHBOURS 20

	vector<double*> tempVertices;
	vector<double*> tempNormals;

	int nv = (int)vertList.size();
	int nt = (int)triList.size();

	tempVertices.resize(nv);
	
	if (computeNormals)
		tempNormals.resize(nv);

	for (int i=0;i<nv;i++)
		tempVertices[i]=allocateVector<double>(3);

	if (computeNormals)
		for (int i=0;i<nv;i++)
			tempNormals[i]=allocateVector<double>(3);

	int *vertdeg[MAX_NEIGHBOURS];	
	bool flagvert;
	for(int i=0;i<MAX_NEIGHBOURS;i++)
		vertdeg[i]=allocateVector<int>(nv);		

	for(int i=0;i<nv;i++)
		vertdeg[0][i]=0;

	for(int i=0;i<nt;i++)
	{
		// first vertex
		flagvert=true;
		for(int j=0;j<vertdeg[0][triList[i][0]];j++)
		{
			if(triList[i][1]==vertdeg[j+1][triList[i][0]])
			{
				flagvert=false;
				break;
			}
		}
		if(flagvert)
		{
			vertdeg[0][triList[i][0]]++;
			vertdeg[vertdeg[0][triList[i][0]]][triList[i][0]]=triList[i][1];
			
		}
		flagvert=true;
		for(int j=0;j<vertdeg[0][triList[i][0]];j++)
		{
			if(triList[i][2]==vertdeg[j+1][triList[i][0]])
			{
				flagvert=false;
				break;
			}
		}
		if(flagvert)
		{
			vertdeg[0][triList[i][0]]++;
			vertdeg[vertdeg[0][triList[i][0]]][triList[i][0]]=triList[i][2];
			
		}
		//second vertex
		flagvert=true;
		for(int j=0;j<vertdeg[0][triList[i][1]];j++)
		{
			if(triList[i][0]==vertdeg[j+1][triList[i][1]])
			{
				flagvert=false;
				break;
			}
		}
		if(flagvert)
		{
			vertdeg[0][triList[i][1]]++;
			vertdeg[vertdeg[0][triList[i][1]]][triList[i][1]]=triList[i][0];
			
		}
		flagvert=true;
		for(int j=0;j<vertdeg[0][triList[i][1]];j++)
		{
			if(triList[i][2]==vertdeg[j+1][triList[i][1]])
			{
				flagvert=false;
				break;
			}
		}
		if(flagvert)
		{
			vertdeg[0][triList[i][1]]++;
			vertdeg[vertdeg[0][triList[i][1]]][triList[i][1]]=triList[i][2];
			
		}
		//third vertex
		flagvert=true;
		for(int j=0;j<vertdeg[0][triList[i][2]];j++)
		{
			if(triList[i][0]==vertdeg[j+1][triList[i][2]])
			{
				flagvert=false;
				break;
			}
		}
		if(flagvert)
		{
			vertdeg[0][triList[i][2]]++;
			vertdeg[vertdeg[0][triList[i][2]]][triList[i][2]]=triList[i][0];
			
		}
		flagvert=true;
		for(int j=0;j<vertdeg[0][triList[i][2]];j++)
		{
			if(triList[i][1]==vertdeg[j+1][triList[i][2]])
			{
				flagvert=false;
				break;
			}
		}
		if(flagvert)
		{
			vertdeg[0][triList[i][2]]++;
			vertdeg[vertdeg[0][triList[i][2]]][triList[i][2]]=triList[i][1];
			
		}
	}
	
	for(int i=0;i<nv;i++)
	{
		tempVertices[i][0]=0;
		tempVertices[i][1]=0;
		tempVertices[i][2]=0;

		if (computeNormals)
		{
			tempNormals[i][0]=0;
			tempNormals[i][1]=0;
			tempNormals[i][2]=0;		
		}
			
		int div = 0;		

		// forward
		// primary neighbours	
		for(int j=0;j<vertdeg[0][i];j++)
		{
			tempVertices[i][0]+=vertList[vertdeg[j+1][i]][0];
			tempVertices[i][1]+=vertList[vertdeg[j+1][i]][1];
			tempVertices[i][2]+=vertList[vertdeg[j+1][i]][2];
								
			div++;
		}
		tempVertices[i][0]/=(div);
		tempVertices[i][1]/=(div);
		tempVertices[i][2]/=(div);			
		
		tempVertices[i][0]=0.5*(vertList[i][0]+tempVertices[i][0]);
		tempVertices[i][1]=0.5*(vertList[i][1]+tempVertices[i][1]);
		tempVertices[i][2]=0.5*(vertList[i][2]+tempVertices[i][2]);

		if (computeNormals)
		{
			div = 0;		
			// forward
			// primary neighbours	
			for(int j=0;j<vertdeg[0][i];j++)
			{
				tempNormals[i][0]+=normalsList[vertdeg[j+1][i]][0];
				tempNormals[i][1]+=normalsList[vertdeg[j+1][i]][1];
				tempNormals[i][2]+=normalsList[vertdeg[j+1][i]][2];									
				div++;
			}
			tempNormals[i][0]/=(div);
			tempNormals[i][1]/=(div);
			tempNormals[i][2]/=(div);			

			tempNormals[i][0]=0.5*(normalsList[i][0]+tempNormals[i][0]);
			tempNormals[i][1]=0.5*(normalsList[i][1]+tempNormals[i][1]);
			tempNormals[i][2]=0.5*(normalsList[i][2]+tempNormals[i][2]);

			double tt;
			NORMALIZE(tempNormals[i],tt);
		}
	}

	// save all
	for(int i=0;i<nv;i++)
	{
		vertList[i][0]=tempVertices[i][0];
		vertList[i][1]=tempVertices[i][1];
		vertList[i][2]=tempVertices[i][2];
	}	

	if (computeNormals)
		for(int i=0;i<nv;i++)
		{
			normalsList[i][0]=tempNormals[i][0];
			normalsList[i][1]=tempNormals[i][1];
			normalsList[i][2]=tempNormals[i][2];
		}	

	//delete all
	for (int i=0;i<nv;i++)
		deleteVector<double>(tempVertices[i]);

	if (computeNormals)
		for (int i=0;i<nv;i++)
			deleteVector<double>(tempNormals[i]);

	for (int i=0;i<MAX_NEIGHBOURS;i++)
		deleteVector<int>(vertdeg[i]);

	int format = deduceFormat();		
	bool f = saveMesh(format,revert,fn,vertList,triList,normalsList);
	if (!f)
		cout << endl << ERR << "Problems in saving the mesh!";

	return;
}

int Surface::deduceFormat()
{
	int format=-1;
	if (!saveMSMS)
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
	short insideness[8];
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;
	short* status = delphi->status;

	if (vertInd == 3)
	{
		//insideness[0]=delphi->status[i][j][k];
		//insideness[0]=STATUSMAP(i,j,k,NX,NY);
		insideness[0]=read3DVector<short>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0]>0)
			insideness[0]=+1;

		if (i-1>=0)
			//insideness[1]=delphi->status[i-1][j][k];
			//insideness[1]=STATUSMAP((i-1),j,k,NX,NY);
			insideness[1]=read3DVector<short>(status,i-1,j,k,NX,NY,NZ);
		else
			insideness[1]=+1;
		
		if (insideness[1]>0)
			insideness[1]=+1;
	
		if (i-1>=0 && k-1>=0)
			//insideness[2]=delphi->status[i-1][j][k-1];
			//insideness[2]=STATUSMAP((i-1),j,(k-1),NX,NY);
			insideness[2]=read3DVector<short>(status,i-1,j,k-1,NX,NY,NZ);
		else
			insideness[2]=+1;

		if (insideness[2]>0)
			insideness[2]=+1;


		if (k-1>=0)
			//insideness[3]=delphi->status[i][j][k-1];
			//insideness[3]=STATUSMAP(i,j,(k-1),NX,NY);
			insideness[3]=read3DVector<short>(status,i,j,k-1,NX,NY,NZ);
		else
			insideness[3]=+1;

		if (insideness[3]>0)
			insideness[3]=+1;

		if (j+1<NY)
			//insideness[4]=delphi->status[i][j+1][k];
			//insideness[4]=STATUSMAP(i,(j+1),k,NX,NY);
			insideness[4]=read3DVector<short>(status,i,j+1,k,NX,NY,NZ);
		else
			insideness[4]=+1;

		if (insideness[4]>0)
			insideness[4]=+1;

		if (j+1<NY && i-1>=0)
			//insideness[5]=delphi->status[i-1][j+1][k];
			//insideness[5]=STATUSMAP((i-1),(j+1),k,NX,NY);
			insideness[5]=read3DVector<short>(status,i-1,j+1,k,NX,NY,NZ);
		else
			insideness[5]=+1;

		if (insideness[5]>0)
			insideness[5]=+1;

		if (j+1<NY && i-1>=0 && k-1>=0)
			//insideness[6]=delphi->status[i-1][j+1][k-1];
			//insideness[6]=STATUSMAP((i-1),(j+1),(k-1),NX,NY);
			insideness[6]=read3DVector<short>(status,i-1,j+1,k-1,NX,NY,NZ);
		else
			insideness[6]=+1;

		if (insideness[6]>0)
			insideness[6]=+1;

		if (j+1<NY && k-1>=0)
			//insideness[7]=delphi->status[i][j+1][k-1];
			//insideness[7]=STATUSMAP(i,(j+1),(k-1),NX,NY);
			insideness[7]=read3DVector<short>(status,i,j+1,k-1,NX,NY,NZ);
		else
			insideness[7]=+1;

		if (insideness[7]>0)
			insideness[7]=+1;

	}
	else if (vertInd==2)
	{
		//insideness[0]=STATUSMAP(i,j,k,NX,NY);
		insideness[0]=read3DVector<short>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0]>0)
			insideness[0]=+1;

		if (i+1<NX)
			//insideness[1]=STATUSMAP((i+1),(j),(k),NX,NY);
			insideness[1]=read3DVector<short>(status,i+1,j,k,NX,NY,NZ);
		else
			insideness[1]=+1;
		
		if (insideness[1]>0)
			insideness[1]=+1;
	
		if (i+1<NX && j+1<NY)
			//insideness[2]=STATUSMAP((i+1),(j+1),(k),NX,NY);
			insideness[2]=read3DVector<short>(status,i+1,j+1,k,NX,NY,NZ);
		else
			insideness[2]=+1;

		if (insideness[2]>0)
			insideness[2]=+1;


		if (j+1<NY)
			//insideness[3]=STATUSMAP((i),(j+1),(k),NX,NY);
			insideness[3]=read3DVector<short>(status,i,j+1,k,NX,NY,NZ);
		else
			insideness[3]=+1;

		if (insideness[3]>0)
			insideness[3]=+1;

		if (k-1>=0)
			//insideness[4]=STATUSMAP((i),(j),(k-1),NX,NY);
			insideness[4]=read3DVector<short>(status,i,j,k-1,NX,NY,NZ);
		else
			insideness[4]=+1;

		if (insideness[4]>0)
			insideness[4]=+1;

		if (k-1>=0 && i+1<NX)
			//insideness[5]=STATUSMAP((i+1),(j),(k-1),NX,NY);
			insideness[5]=read3DVector<short>(status,i+1,j,k-1,NX,NY,NZ);
		else
			insideness[5]=+1;

		if (insideness[5]>0)
			insideness[5]=+1;

		if (k-1>=0 && i+1<NX && j+1<NY)
			//insideness[6]=STATUSMAP((i+1),(j+1),(k-1),NX,NY);
			insideness[6]=read3DVector<short>(status,i+1,j+1,k-1,NX,NY,NZ);
		else
			insideness[6]=+1;

		if (insideness[6]>0)
			insideness[6]=+1;

		if (k-1>=0 && j+1<NY)
			//insideness[7]=STATUSMAP((i),(j+1),(k-1),NX,NY);
			insideness[7]=read3DVector<short>(status,i,j+1,k-1,NX,NY,NZ);
		else
			insideness[7]=+1;

		if (insideness[7]>0)
			insideness[7]=+1;
	}
	else if (vertInd==1)
	{
		//insideness[0]=STATUSMAP(i,j,k,NX,NY);
		insideness[0]=read3DVector<short>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0]>0)
			insideness[0]=+1;

		if (i+1<NX)
			//insideness[1]=STATUSMAP((i+1),(j),(k),NX,NY);
			insideness[1]=read3DVector<short>(status,i+1,j,k,NX,NY,NZ);
		else
			insideness[1]=+1;
		
		if (insideness[1]>0)
			insideness[1]=+1;
	
		if (i+1<NX && j-1>=0)
			//insideness[2]=STATUSMAP((i+1),(j-1),(k),NX,NY);
			insideness[2]=read3DVector<short>(status,i+1,j-1,k,NX,NY,NZ);
		else
			insideness[2]=+1;

		if (insideness[2]>0)
			insideness[2]=+1;


		if (j-1>=0)
			//insideness[3]=STATUSMAP((i),(j-1),(k),NX,NY);
			insideness[3]=read3DVector<short>(status,i,j-1,k,NX,NY,NZ);
		else
			insideness[3]=+1;

		if (insideness[3]>0)
			insideness[3]=+1;

		if (k-1>=0)
			//insideness[4]=STATUSMAP((i),(j),(k-1),NX,NY);
			insideness[4]=read3DVector<short>(status,i,j,k-1,NX,NY,NZ);
		else
			insideness[4]=+1;

		if (insideness[4]>0)
			insideness[4]=+1;

		if (k-1>=0 && i+1<NX)
			//insideness[5]=STATUSMAP((i+1),(j),(k-1),NX,NY);
			insideness[5]=read3DVector<short>(status,i+1,j,k-1,NX,NY,NZ);
		else
			insideness[5]=+1;

		if (insideness[5]>0)
			insideness[5]=+1;

		if (k-1>=0 && i+1<NX && j-1>=0)
			//insideness[6]=STATUSMAP((i+1),(j-1),(k-1),NX,NY);
			insideness[6]=read3DVector<short>(status,i+1,j-1,k-1,NX,NY,NZ);
		else
			insideness[6]=+1;

		if (insideness[6]>0)
			insideness[6]=+1;

		if (k-1>=0 && j-1>=0)
			//insideness[7]=STATUSMAP((i),(j-1),(k-1),NX,NY);
			insideness[7]=read3DVector<short>(status,i,j-1,k-1,NX,NY,NZ);
		else
			insideness[7]=+1;

		if (insideness[7]>0)
			insideness[7]=+1;	
	}
	else if (vertInd==0)
	{
		//insideness[0]=STATUSMAP(i,j,k,NX,NY);
		insideness[0]=read3DVector<short>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0]>0)
			insideness[0]=+1;

		if (i-1>=0)
			//insideness[1]=STATUSMAP((i-1),j,k,NX,NY);
			insideness[1]=read3DVector<short>(status,i-1,j,k,NX,NY,NZ);
		else
			insideness[1]=+1;
		
		if (insideness[1]>0)
			insideness[1]=+1;
	
		if (i-1>=0 && k-1>=0)
			//insideness[2]=STATUSMAP((i-1),(j),(k-1),NX,NY);
			insideness[2]=read3DVector<short>(status,i-1,j,k-1,NX,NY,NZ);
		else
			insideness[2]=+1;

		if (insideness[2]>0)
			insideness[2]=+1;


		if (k-1>=0)
			//insideness[3]=STATUSMAP((i),(j),(k-1),NX,NY);
			insideness[3]=read3DVector<short>(status,i,j,k-1,NX,NY,NZ);
		else
			insideness[3]=+1;

		if (insideness[3]>0)
			insideness[3]=+1;

		if (j-1>=0)
			//insideness[4]=STATUSMAP((i),(j-1),(k),NX,NY);
			insideness[4]=read3DVector<short>(status,i,j-1,k,NX,NY,NZ);
		else
			insideness[4]=+1;

		if (insideness[4]>0)
			insideness[4]=+1;

		if (j-1>=0 && i-1>=0)
			//insideness[5]=STATUSMAP((i-1),(j-1),(k),NX,NY);
			insideness[5]=read3DVector<short>(status,i-1,j-1,k,NX,NY,NZ);
		else
			insideness[5]=+1;

		if (insideness[5]>0)
			insideness[5]=+1;

		if (j-1>=0 && i-1>=0 && k-1>=0)
			//insideness[6]=STATUSMAP((i-1),(j-1),(k-1),NX,NY);
			insideness[6]=read3DVector<short>(status,i-1,j-1,k-1,NX,NY,NZ);
		else
			insideness[6]=+1;

		if (insideness[6]>0)
			insideness[6]=+1;

		if (j-1>=0 && k-1>=0)
			//insideness[7]=STATUSMAP((i),(j-1),(k-1),NX,NY);
			insideness[7]=read3DVector<short>(status,i,j-1,k-1,NX,NY,NZ);
		else
			insideness[7]=+1;

		if (insideness[7]>0)
			insideness[7]=+1;
	
	}
	else if (vertInd==7)
	{
		//insideness[0]=STATUSMAP(i,j,k,NX,NY);
		insideness[0]=read3DVector<short>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0]>0)
			insideness[0]=+1;

		if (i-1>=0)
			//insideness[1]=STATUSMAP((i-1),j,k,NX,NY);
			insideness[1]=read3DVector<short>(status,i-1,j,k,NX,NY,NZ);
		else
			insideness[1]=+1;
		
		if (insideness[1]>0)
			insideness[1]=+1;
	
		if (i-1>=0 && k+1<NZ)
			//insideness[2]=STATUSMAP((i-1),j,(k+1),NX,NY);
			insideness[2]=read3DVector<short>(status,i-1,j,k+1,NX,NY,NZ);
		else
			insideness[2]=+1;

		if (insideness[2]>0)
			insideness[2]=+1;


		if (k+1<NZ)
			//insideness[3]=STATUSMAP((i),(j),(k+1),NX,NY);
			insideness[3]=read3DVector<short>(status,i,j,k+1,NX,NY,NZ);
		else
			insideness[3]=+1;

		if (insideness[3]>0)
			insideness[3]=+1;

		if (j+1<NY)
			//insideness[4]=STATUSMAP((i),(j+1),(k),NX,NY);
			insideness[4]=read3DVector<short>(status,i,j+1,k,NX,NY,NZ);
		else
			insideness[4]=+1;

		if (insideness[4]>0)
			insideness[4]=+1;

		if (j+1<NY && i-1>=0)
			//insideness[5]=STATUSMAP((i-1),(j+1),(k),NX,NY);
			insideness[5]=read3DVector<short>(status,i-1,j+1,k,NX,NY,NZ);
		else
			insideness[5]=+1;

		if (insideness[5]>0)
			insideness[5]=+1;

		if (j+1<NY && i-1>=0 && k+1<NZ)
			//insideness[6]=STATUSMAP((i-1),(j+1),(k+1),NX,NY);
			insideness[6]=read3DVector<short>(status,i-1,j+1,k+1,NX,NY,NZ);
		else
			insideness[6]=+1;

		if (insideness[6]>0)
			insideness[6]=+1;

		if (j+1<NY && k+1<NZ)
			//insideness[7]=STATUSMAP((i),(j+1),(k+1),NX,NY);
			insideness[7]=read3DVector<short>(status,i,j+1,k+1,NX,NY,NZ);
		else
			insideness[7]=+1;

		if (insideness[7]>0)
			insideness[7]=+1;

	}
	else if (vertInd==6)
	{
		//insideness[0]=STATUSMAP(i,j,k,NX,NY);
		insideness[0]=read3DVector<short>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0]>0)
			insideness[0]=+1;

		if (i+1<NX)
			//insideness[1]=STATUSMAP((i+1),(j),(k),NX,NY);
			insideness[1]=read3DVector<short>(status,i+1,j,k,NX,NY,NZ);
		else
			insideness[1]=+1;
		
		if (insideness[1]>0)
			insideness[1]=+1;
	
		if (i+1<NX && j+1<NY)
			//insideness[2]=STATUSMAP((i+1),(j+1),(k),NX,NY);
			insideness[2]=read3DVector<short>(status,i+1,j+1,k,NX,NY,NZ);
		else
			insideness[2]=+1;

		if (insideness[2]>0)
			insideness[2]=+1;


		if (j+1<NY)
			//insideness[3]=STATUSMAP((i),(j+1),(k),NX,NY);
			insideness[3]=read3DVector<short>(status,i,j+1,k,NX,NY,NZ);
		else
			insideness[3]=+1;

		if (insideness[3]>0)
			insideness[3]=+1;

		if (k+1<NZ)
			//insideness[4]=STATUSMAP((i),(j),(k+1),NX,NY);
			insideness[4]=read3DVector<short>(status,i,j,k+1,NX,NY,NZ);
		else
			insideness[4]=+1;

		if (insideness[4]>0)
			insideness[4]=+1;

		if (k+1<NZ && i+1<NX)
			//insideness[5]=STATUSMAP((i+1),(j),(k+1),NX,NY);
			insideness[5]=read3DVector<short>(status,i+1,j,k+1,NX,NY,NZ);
		else
			insideness[5]=+1;

		if (insideness[5]>0)
			insideness[5]=+1;

		if (k+1<NZ && i+1<NX && j+1<NY)
			//insideness[6]=STATUSMAP((i+1),(j+1),(k+1),NX,NY);
			insideness[6]=read3DVector<short>(status,i+1,j+1,k+1,NX,NY,NZ);
		else
			insideness[6]=+1;

		if (insideness[6]>0)
			insideness[6]=+1;

		if (k+1<NZ && j+1<NY)
			//insideness[7]=STATUSMAP((i),(j+1),(k+1),NX,NY);
			insideness[7]=read3DVector<short>(status,i,j+1,k+1,NX,NY,NZ);
		else
			insideness[7]=+1;

		if (insideness[7]>0)
			insideness[7]=+1;
	
	}
	else if (vertInd==5)
	{
		//insideness[0]=STATUSMAP(i,j,k,NX,NY);
		insideness[0]=read3DVector<short>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0]>0)
			insideness[0]=+1;

		if (i+1<NX)
			//insideness[1]=STATUSMAP((i+1),(j),(k),NX,NY);
			insideness[1]=read3DVector<short>(status,i+1,j,k,NX,NY,NZ);
		else
			insideness[1]=+1;
		
		if (insideness[1]>0)
			insideness[1]=+1;
	
		if (i+1<NX && j-1>=0)
			//insideness[2]=STATUSMAP((i+1),(j-1),(k),NX,NY);
			insideness[2]=read3DVector<short>(status,i+1,j-1,k,NX,NY,NZ);
		else
			insideness[2]=+1;

		if (insideness[2]>0)
			insideness[2]=+1;


		if (j-1>=0)
			//insideness[3]=STATUSMAP((i),(j-1),(k),NX,NY);
			insideness[3]=read3DVector<short>(status,i,j-1,k,NX,NY,NZ);
		else
			insideness[3]=+1;

		if (insideness[3]>0)
			insideness[3]=+1;

		if (k+1<NZ)
			//insideness[4]=STATUSMAP((i),(j),(k+1),NX,NY);
			insideness[4]=read3DVector<short>(status,i,j,k+1,NX,NY,NZ);
		else
			insideness[4]=+1;

		if (insideness[4]>0)
			insideness[4]=+1;

		if (k+1<NZ && i+1<NX)
			//insideness[5]=STATUSMAP((i+1),(j),(k+1),NX,NY);
			insideness[5]=read3DVector<short>(status,i+1,j,k+1,NX,NY,NZ);
		else
			insideness[5]=+1;

		if (insideness[5]>0)
			insideness[5]=+1;

		if (k+1<NZ && i+1<NX && j-1>=0)
			//insideness[6]=STATUSMAP((i+1),(j-1),(k+1),NX,NY);
			insideness[6]=read3DVector<short>(status,i+1,j-1,k+1,NX,NY,NZ);
		else
			insideness[6]=+1;

		if (insideness[6]>0)
			insideness[6]=+1;

		if (k+1<NZ && j-1>=0)
			//insideness[7]=STATUSMAP((i),(j-1),(k+1),NX,NY);
			insideness[7]=read3DVector<short>(status,i,j-1,k+1,NX,NY,NZ);
		else
			insideness[7]=+1;

		if (insideness[7]>0)
			insideness[7]=+1;		
	}
	else if (vertInd==4)
	{
		//insideness[0]=STATUSMAP(i,j,k,NX,NY);
		insideness[0]=read3DVector<short>(status,i,j,k,NX,NY,NZ);

		// saturate outvalue to +1 such that cavities or temporary out points are all the same
		if (insideness[0]>0)
			insideness[0]=+1;

		if (i-1>=0)
			//insideness[1]=STATUSMAP((i-1),j,k,NX,NY);
			insideness[1]=read3DVector<short>(status,i-1,j,k,NX,NY,NZ);
		else
			insideness[1]=+1;
		
		if (insideness[1]>0)
			insideness[1]=+1;
	
		if (i-1>=0 && k+1<NZ)
			//insideness[2]=STATUSMAP((i-1),j,(k+1),NX,NY);
			insideness[2]=read3DVector<short>(status,i-1,j,k+1,NX,NY,NZ);
		else
			insideness[2]=+1;

		if (insideness[2]>0)
			insideness[2]=+1;


		if (k+1<NZ)
			//insideness[3]=STATUSMAP((i),(j),(k+1),NX,NY);
			insideness[3]=read3DVector<short>(status,i,j,k+1,NX,NY,NZ);
		else
			insideness[3]=+1;

		if (insideness[3]>0)
			insideness[3]=+1;

		if (j-1>=0)
			//insideness[4]=STATUSMAP((i),(j-1),(k),NX,NY);
			insideness[4]=read3DVector<short>(status,i,j-1,k,NX,NY,NZ);
		else
			insideness[4]=+1;

		if (insideness[4]>0)
			insideness[4]=+1;

		if (j-1>=0 && i-1>=0)
			//insideness[5]=STATUSMAP((i-1),(j-1),(k),NX,NY);
			insideness[5]=read3DVector<short>(status,i-1,j-1,k,NX,NY,NZ);
		else
			insideness[5]=+1;

		if (insideness[5]>0)
			insideness[5]=+1;

		if (j-1>=0 && i-1>=0 && k+1<NZ)
			//insideness[6]=STATUSMAP((i-1),(j-1),(k+1),NX,NY);
			insideness[6]=read3DVector<short>(status,i-1,j-1,k+1,NX,NY,NZ);
		else
			insideness[6]=+1;

		if (insideness[6]>0)
			insideness[6]=+1;

		if (j-1>=0 && k+1<NZ)
			//insideness[7]=STATUSMAP((i),(j-1),(k+1),NX,NY);
			insideness[7]=read3DVector<short>(status,i,j-1,k+1,NX,NY,NZ);
		else
			insideness[7]=+1;

		if (insideness[7]>0)
			insideness[7]=+1;
	
	}
	else
	{
		cout << endl << ERR << "Unknown vertex index!";
		return false;
	}
	
	int acc = 0;
	// decide insidness
	for (int i=0;i<8;i++)
	{
		acc += (int)insideness[i];
	}
	// positively biased decision
	if (acc>0)
		return +1;
	else 
		return -1;
}


void Surface::vertexInterp(double isolevel, double* p1, double* p2, double valp1, double valp2,double* p) 
{
   double mu;
   mu = (isolevel - valp1) / (valp2 - valp1);
   double temp[3];
   SUB(temp,p2,p1);
   ADD_MUL(p,p1,temp,mu);
}


inline int Surface::classifyCube(double* vertexValues,double isolevel)
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
   the analytically computed intersection point. In this last case it is garanted
   that every point of the mesh belongs to the surface; ix,iy,iz are passed in order to query
   the intersection matrices
*/
//inline int Surface::getTriangles(double* vertexValues,double** vertexPos,double isolevel, int** triangles,int ix,int iy,
//						  int iz,int NX, int NY,int NZ,int* xedge,int* yedge,int* zedge,int* xedge_down,int* yedge_down)
inline int Surface::getTriangles(double* vertexValues,double** vertexPos,double isolevel, int** triangles,int ix,int iy,
						  int iz,int NX, int NY,int NZ)

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
	
	int vertlist_indexes[12]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
		
	//   Determine the index into the edge table which
	//   tells us which vertices are inside of the surface   
	cubeindex = classifyCube(vertexValues,isolevel);
	
	// Cube is entirely in/out of the surface 
	if (cubeindex == -1)
		return 0;
   
	//int vert_offset,i=ix,j=iy;
	int i=ix,j=iy;

	// the vertices are already analytically computed, they are only recovered. 
	// 0 buffered
	if (edgeTable[cubeindex] & 1)	
	{
		/*
		vert_offset = XEDGE_DOWN(i,j-1,NX);
		vertlist_indexes[0] = vert_offset;		
		if (vert_offset==-1)
		*/
			vertlist_indexes[0] = intersectionsMatrixAlongX->at(ix,iy,iz);						
	}

	// 1 buffered
	if (edgeTable[cubeindex] & 2)
	{
		/*
		vert_offset = YEDGE_DOWN(i,j,NX);
		vertlist_indexes[1] = vert_offset;
		if (vert_offset==-1)
		*/
			vertlist_indexes[1] = intersectionsMatrixAlongY->at(ix+1,iy,iz);
	}

	// 2 buffered
	if (edgeTable[cubeindex] & 4)
	{
		/*
		vert_offset = XEDGE_DOWN(i,j,NX);
		vertlist_indexes[2] = vert_offset;
		if (vert_offset==-1)
		*/
			vertlist_indexes[2] = intersectionsMatrixAlongX->at(ix,iy+1,iz);
	}

	// 3 buffered
	if (edgeTable[cubeindex] & 8)
	{
		/*
		vert_offset = YEDGE_DOWN(i-1,j,NX);
		vertlist_indexes[3] = vert_offset;
		if (vert_offset==-1)
		*/
			vertlist_indexes[3] = intersectionsMatrixAlongY->at(ix,iy,iz);
	}
	
	// 4 buffered
	if (edgeTable[cubeindex] & 16)
	{
		/*
		vert_offset = XEDGE(i,j-1,NX);
		vertlist_indexes[4] = vert_offset;
		if (vert_offset==-1)
		*/
			vertlist_indexes[4] = intersectionsMatrixAlongX->at(ix,iy,iz+1);
	}

	// 5
	if (edgeTable[cubeindex] & 32)
	{
		vertlist_indexes[5]	= intersectionsMatrixAlongY->at(ix+1,iy,iz+1);
		//YEDGE(i,j,NX) = vertlist_indexes[5];
	}

	// 6
	if (edgeTable[cubeindex] & 64)
	{
		vertlist_indexes[6]	=	intersectionsMatrixAlongX->at(ix,iy+1,iz+1);
		//XEDGE(i,j,NX) = vertlist_indexes[6];
	}

	// 7 buffered
	if (edgeTable[cubeindex] & 128)
	{
		/*
		vert_offset = YEDGE(i-1,j,NX);
		vertlist_indexes[7] = vert_offset;
		if (vert_offset==-1)
		*/
			vertlist_indexes[7] = intersectionsMatrixAlongY->at(ix,iy,iz+1);
	}

	// 8 buffered
	if (edgeTable[cubeindex] & 256)
	{
		/*
		vert_offset = ZEDGE(i-1,j-1,NX);			
		vertlist_indexes[8] = vert_offset;
		if (vert_offset==-1)
		*/
			vertlist_indexes[8] = intersectionsMatrixAlongZ->at(ix,iy,iz);
	}

	// 9 buffered
	if (edgeTable[cubeindex] & 512)
	{
		/*
		vert_offset = ZEDGE(i,j-1,NX);
		vertlist_indexes[9] = vert_offset;
		if (vert_offset==-1)
		*/
			vertlist_indexes[9] = intersectionsMatrixAlongZ->at(ix+1,iy,iz);
	}
	
	// 10
	if (edgeTable[cubeindex] & 1024)
	{
		vertlist_indexes[10] =	intersectionsMatrixAlongZ->at(ix+1,iy+1,iz);			
		//ZEDGE(i,j,NX)=vertlist_indexes[10];
	}

	// 11 buffered
	if (edgeTable[cubeindex] & 2048)
	{ 
		/*
		vert_offset = ZEDGE(i-1,j,NX);
		vertlist_indexes[11] = vert_offset;
		if (vert_offset==-1)
		*/
			vertlist_indexes[11] = intersectionsMatrixAlongZ->at(ix,iy+1,iz);
	}

   // Create the triangles
   ntriang = 0;
   for (i=0;triTable[cubeindex][i]!=-1;i+=3) 
   {
      if (vertlist_indexes[triTable[cubeindex][i  ]]==-1)
	  {		  
		  cout << endl << WARN << "Mesh with hole!";
		  continue;
	  }
	  if (vertlist_indexes[triTable[cubeindex][i+1  ]]==-1)
	  {
		  cout << endl << WARN << "Mesh with hole!";
		  continue;
	  }
	  if (vertlist_indexes[triTable[cubeindex][i+2  ]]==-1)
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
	if (delphi==NULL)
	{
		cout << endl << WARN << "Cannot backup status if grid is not allocated";
		return;
	}
	int tot = delphi->nx*delphi->ny*delphi->nz;

	if (delphi->tempStatus!=NULL)
		deleteVector<short>(delphi->tempStatus);

	delphi->tempStatus = allocateVector<short>(tot);
		
	for (int i=0;i<tot;i++)
		delphi->tempStatus[i]=delphi->status[i];
}

void Surface::removeBackupStatus()
{
	if (delphi->tempStatus!=NULL)
		deleteVector<short>(delphi->tempStatus);
	delphi->tempStatus=NULL;
}

int Surface::linkCavities(short* st1,short* st2)
{
	// this is the status of the slim surface
	short* tempStatus = st1;
	short* tempStatus2 = st2;
	short* status = delphi->status;

	if (tempStatus==NULL)
	{
		cout << endl << WARN << "I cannot reconstruct links if I have not reference status";
		return 0;
	}
		
	// analyze each cavity before and after fill
	vector < vector <int* >* >::iterator it;		
	int i=0;
	
	int NX = delphi->nx;
	int NY = delphi->ny;	
	int NZ = delphi->nz;
		
	// cavity under link analysis
	int cavityId = STATUS_FIRST_CAV;

	vector <double>::iterator sizeIt;// cavitiesSize;				
	vector <bool>::iterator flagIt; //cavitiesFlag;

	//buildAtomsMap();
				
	for (it=delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++,cavityId++)
	{
		if (delphi->cavitiesFlag[cavityId-STATUS_FIRST_CAV]==true)
				continue;

		vector<int*>* vec1 = (*it);

		if (vec1->size()==0)
			continue;

		// paired cavity under analysis
		int checkCavId = cavityId+1;

		// check if, for the all the other current cavities, exist a point
		// which had the same id before filtering. If yes a link is estabilished.
		for (checkCavId = cavityId+1;  (checkCavId-STATUS_FIRST_CAV)<(int)delphi->cavitiesVec->size(); checkCavId++)			
		{				
			if (delphi->cavitiesFlag[checkCavId-STATUS_FIRST_CAV]==true)
				continue;
				
			vector<int*>* vec2 = delphi->cavitiesVec->at(checkCavId-STATUS_FIRST_CAV);

			// this may happen if this cavity have been already aggregated
			if (vec2->size()==0)
				continue;
				
			// check the distance between the two nearest bgps of the two cavities.
			// a bgp in this case is defined as a grid point in which at least 
			// one grid point is out

			// these are respectively the two current grid points which are compared
			int *v1,*v2;

			vector<int*>::iterator it_c1;
			vector<int*>::iterator it_c2;
			double minDist2 = INFINITY;
			unsigned int winner1, winner2;

			for (unsigned int c1 = 0; c1<vec1->size(); c1++)
			{
				v1 = vec1->at(c1);
				
				//int ref = STATUSMAP(v1[0],v1[1],v1[2],NX,NY);
				int ref = read3DVector<short>(status,v1[0],v1[1],v1[2],NX,NY,NZ);

				// if is bgp get distance
				/*
				if ((STATUSMAP(v1[0]+1,v1[1],v1[2],NX,NY)!=ref) ||
					(STATUSMAP(v1[0]-1,v1[1],v1[2],NX,NY)!=ref) ||
					(STATUSMAP(v1[0],v1[1]+1,v1[2],NX,NY)!=ref) ||
					(STATUSMAP(v1[0],v1[1]-1,v1[2],NX,NY)!=ref) ||
					(STATUSMAP(v1[0],v1[1],v1[2]+1,NX,NY)!=ref) ||
					(STATUSMAP(v1[0],v1[1],v1[2]-1,NX,NY)!=ref))*/

				if ((read3DVector<short>(status,v1[0]+1,v1[1],v1[2],NX,NY,NZ)!=ref) ||
					(read3DVector<short>(status,v1[0]-1,v1[1],v1[2],NX,NY,NZ)!=ref) ||
					(read3DVector<short>(status,v1[0],v1[1]+1,v1[2],NX,NY,NZ)!=ref) ||
					(read3DVector<short>(status,v1[0],v1[1]-1,v1[2],NX,NY,NZ)!=ref) ||
					(read3DVector<short>(status,v1[0],v1[1],v1[2]+1,NX,NY,NZ)!=ref) ||
					(read3DVector<short>(status,v1[0],v1[1],v1[2]-1,NX,NY,NZ)!=ref))

				{
					// do check
				}
				else
					continue;

				for (unsigned int c2 = 0;c2<vec2->size();c2++)
				{
					v2 = vec2->at(c2);
					//ref = STATUSMAP(v2[0],v2[1],v2[2],NX,NY);
					ref = read3DVector<short>(status,v2[0],v2[1],v2[2],NX,NY,NZ);

					// if is bgp get distance
					/*
					if ((STATUSMAP(v2[0]+1,v2[1],v2[2],NX,NY)!=ref) ||
						(STATUSMAP(v2[0]-1,v2[1],v2[2],NX,NY)!=ref) ||
						(STATUSMAP(v2[0],v2[1]+1,v2[2],NX,NY)!=ref) ||
						(STATUSMAP(v2[0],v2[1]-1,v2[2],NX,NY)!=ref) ||
						(STATUSMAP(v2[0],v2[1],v2[2]+1,NX,NY)!=ref) ||
						(STATUSMAP(v2[0],v2[1],v2[2]-1,NX,NY)!=ref))
					*/
					if ((read3DVector<short>(status,v2[0]+1,v2[1],v2[2],NX,NY,NZ)!=ref) ||
						(read3DVector<short>(status,v2[0]-1,v2[1],v2[2],NX,NY,NZ)!=ref) ||
						(read3DVector<short>(status,v2[0],v2[1]+1,v2[2],NX,NY,NZ)!=ref) ||
						(read3DVector<short>(status,v2[0],v2[1]-1,v2[2],NX,NY,NZ)!=ref) ||
						(read3DVector<short>(status,v2[0],v2[1],v2[2]+1,NX,NY,NZ)!=ref) ||
						(read3DVector<short>(status,v2[0],v2[1],v2[2]-1,NX,NY,NZ)!=ref))
					{
						// get distance
						double a = (delphi->x[v1[0]]-delphi->x[v2[0]])*(delphi->x[v1[0]]-delphi->x[v2[0]]);
						double b = (delphi->y[v1[1]]-delphi->y[v2[1]])*(delphi->y[v1[1]]-delphi->y[v2[1]]);
						double c = (delphi->z[v1[2]]-delphi->z[v2[2]])*(delphi->z[v1[2]]-delphi->z[v2[2]]);

						double dist2 = a+b+c;
						if (dist2<minDist2)
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

			if (minDist2==INFINITY)
			{
				cout << endl << ERR << "During linkage, two non null cavities/pockets have no minimum distance bgps";
				exit(-1);
			}

			bool merge = false;
			//printf("\n Min distance between %d %d is %f",cavityId-4,checkCavId-4,sqrt(minDist2));

			// bound on the max path
			double refDist = 6;
			if (sqrt(minDist2)>refDist)
				merge = false;

			// if we are under a threshold distance a full check is needed
			// check if one can move freely inside the small probe surface between the two nearest bgps
			// while at the mean time never going out the 1.4 surface.
			else
			{
				int maxMoves = (int)(refDist*delphi->scale+0.5); 
				int moves = maxMoves;
				v1 = vec1->at(winner1);
				v2 = vec2->at(winner2);				

				bool debug = false;
				//if (cavityId-4==6 && checkCavId-4==17)
				//	debug = true;

				//printf("\n %d %d %d -> %d %d %d",v1[0],v1[1],v1[2],v2[0],v2[1],v2[2]);
				merge = innerFloodFill(v1,v2,st1,st2,moves,debug);
				double geodetic = moves*delphi->side;
				
				if (!merge)
				{
				//	printf("\n Non merging, num moves %d",moves);
				}
				if (merge && geodetic>1.5*sqrt(minDist2))
				{
				//	printf("\n Geodetic is %f and max dist is %f rejecting link",geodetic,1.5*sqrt(minDist2));
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
				vector < int* >::iterator it_vec2;
				for (it_vec2=vec2->begin();it_vec2!=vec2->end();it_vec2++)
				{
					int* v = *it_vec2;
					vec1->push_back(v);
				}
					
				// pointers are not deleted
				vec2->clear();			

				// sum volumes of the aggregated cavities
				delphi->cavitiesSize[cavityId-STATUS_FIRST_CAV]+=delphi->cavitiesSize[checkCavId-STATUS_FIRST_CAV];
				// nullify the volume to flag it
				delphi->cavitiesSize[checkCavId-STATUS_FIRST_CAV] = 0;
				// assure it is active the merged cavity
				delphi->cavitiesFlag[cavityId-STATUS_FIRST_CAV]=false;
			}

		}
	}

	//disposeAtomsMap();

	sizeIt = delphi->cavitiesSize.begin();
	flagIt = delphi->cavitiesFlag.begin();
	it = delphi->cavitiesVec->begin();

	int numRemoved = 0;

	// remove all the null cavities
	while(1)		
	{
		if (it==delphi->cavitiesVec->end())
			break;

		vector<int*>* vec1 = (*it);

		// the merged cavities are removed
		if (vec1->size()==0)
		{
			numRemoved++;
			// delete vector of pointers to points
			delete vec1;
			
			// delete entry
			delphi->cavitiesVec->erase(it);
			
			// double check it is the right guy
			if ((*sizeIt)!=0)
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
	for (it = delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++)
	{
		vector<int*>::iterator inner;
		vector<int*>* cav = (*it);			
		for (inner = cav->begin();inner!=cav->end();inner++)
		{
			int* v = *inner;
			// conserve the support cavity coding
			
			/*
			if (STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY)>0)
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = new_id;
			else
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = -new_id;
			*/
			
			if (read3DVector<short>(status,(v[0]),(v[1]),(v[2]),NX,NY,NZ)>0)
				write3DVector<short>(status,new_id,(v[0]),(v[1]),(v[2]),NX,NY,NZ);
			else
				write3DVector<short>(status,-new_id,(v[0]),(v[1]),(v[2]),NX,NY,NZ);

		}
		new_id++;
	}
	return numRemoved;
}

/* old method
int Surface::linkCavities()
{
	short* tempStatus = delphi->tempStatus;
	short* status = delphi->status;

	if (tempStatus==NULL)
	{
		cout << endl << WARN << "I cannot reconstruct links if I have not reference status";
		return 0;
	}
		
	// analyze each cavity before and after fill
	vector < vector <int* >* >::iterator it;		
	int i=0;
	
	int NX = delphi->nx;
	int NY = delphi->ny;	
	int NZ = delphi->nz;
		
	// cavity under link analysis
	int cavityId = STATUS_FIRST_CAV;

	vector <double>::iterator sizeIt;// cavitiesSize;				
	vector <bool>::iterator flagIt; //cavitiesFlag;
				
	for (it=delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++,cavityId++)
	{
		if (delphi->cavitiesFlag[cavityId-STATUS_FIRST_CAV]==true)
				continue;

		vector<int*>* vec1 = (*it);

		if (vec1->size()==0)
			continue;

		// get the first point in cavity 			
		int* v=*(vec1->begin());
			
		// cavity id in the old map
		int oldCavityId = TEMP_STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY);

		// paired cavity under analysis
		int checkCavId = cavityId+1;

		// check if, for the all the other current cavities, exist a point
		// which had the same id before filtering. If yes a link is estabilished.
		for (checkCavId = cavityId+1;  (checkCavId-STATUS_FIRST_CAV)<(int)delphi->cavitiesVec->size(); checkCavId++)			
		{				
			if (delphi->cavitiesFlag[checkCavId-STATUS_FIRST_CAV]==true)
				continue;
				
			vector<int*>* vec2 = delphi->cavitiesVec->at(checkCavId-STATUS_FIRST_CAV);

			// this may happen if this cavity have been agregated
			if (vec2->size()==0)
				continue;
				
			// take a random point and check the previous map
			int* v = *(vec2->begin());
			int oldCavityId2 = TEMP_STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY);

			// these two new cavities need to be linked
			// because they were linked before filtering
			if (oldCavityId==oldCavityId2)
			{
				//printf("\n%sLinking cavities (%d,%d)",INFO,checkCavId-STATUS_FIRST_CAV,cavityId-STATUS_FIRST_CAV);
				//printf("\n old value %d",oldCavityId);
				//getchar();
					
				// merge the grid points
				vector < int* >::iterator it_vec2;
				for (it_vec2=vec2->begin();it_vec2!=vec2->end();it_vec2++)
				{
					int* v = *it_vec2;
					vec1->push_back(v);
				}
					
				// pointers are not deleted
				vec2->clear();			

				// sum volumes of the agregate cavities
				delphi->cavitiesSize[cavityId-STATUS_FIRST_CAV]+=delphi->cavitiesSize[checkCavId-STATUS_FIRST_CAV];
				// nullify the volume to flag it
				delphi->cavitiesSize[checkCavId-STATUS_FIRST_CAV] = 0;
				// assure it is active the merged cavity
				delphi->cavitiesFlag[cavityId-STATUS_FIRST_CAV]=false;
			}

		}
	}

	sizeIt = delphi->cavitiesSize.begin();
	flagIt = delphi->cavitiesFlag.begin();
	it = delphi->cavitiesVec->begin();

	int numRemoved = 0;

	// remove all the null cavities
	while(1)		
	{
		if (it==delphi->cavitiesVec->end())
			break;

		vector<int*>* vec1 = (*it);

		// the merged cavities are removed
		if (vec1->size()==0)
		{
			numRemoved++;
			// delete vector of pointers to points
			delete vec1;
			
			// delete entry
			delphi->cavitiesVec->erase(it);
			
			// double check it is the right guy
			if ((*sizeIt)!=0)
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
	for (it = delphi->cavitiesVec->begin();it!=delphi->cavitiesVec->end();it++)
	{
		vector<int*>::iterator inner;
		vector<int*>* cav = (*it);			
		for (inner = cav->begin();inner!=cav->end();inner++)
		{
			int* v = *inner;
			// conserve the support cavity coding
			if (STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY)>0)
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = new_id;
			else
				STATUSMAP((v[0]),(v[1]),(v[2]),NX,NY) = -new_id;
		}
		new_id++;
	}
	return numRemoved;
}
*/

Surface& Surface::operator-=(Surface& surf2) 
{
	Timer chrono;
	chrono.start();
	
	double duration = 0;	
	difference(&surf2);
	
	duration = chrono.stop();
	cout << endl << INFO << "Diff. Step 1 " << duration << " [s]";
	
	///////////////////// connolly filter /////////////////////////////////
	setProbeRadius(1.4);			
	
	Timer chrono2;
	chrono2.start();

	filterCavities(true);			
	cav2out();

	duration = chrono2.stop();
	
	cout << endl << INFO << "Diff. Step 2 " << duration << " [s]";
	cout.flush();
	
	chrono2.start();	
	
	///////////////////// get filtered cavities ///////////////////////////
	int ff = getCavities(STATUS_POINT_OUT);
	if (!ff)
	{
		//cout << endl << WARN <<  "Zero difference after filtration!";
		return *this;
	}

	duration = chrono2.stop();
	cout << endl << INFO << "Diff. Step 3 " << duration << " [s]";
	

	cout.flush();
	return *this;  
}

bool Surface::innerFloodFill(int* start,int* target,short* tempStatus,short* tempStatus2,int& maxMoves,bool debug)
{	
	short* status = delphi->status;	
	Octree<int> visited(1024,0);

	// ix,iy,iz and integer distance between the origin and that point
	queue<pair< pair<int,int>,pair<int,int> > > myqueue;
	pair<pair<int,int>,pair<int,int> > quartet;
	int cix,ciy,ciz,dist;
	
	int NX = delphi->nx;
	int NY = delphi->ny;
	int NZ = delphi->nz;

	double* X = delphi->x;
	double* Y = delphi->y;
	double* Z = delphi->z;

	int ix = start[0];
	int iy = start[1];
	int iz = start[2];

	myqueue.push(pair<pair<int,int>,pair<int,int > >(pair<int,int>(ix,iy),pair<int,int>(iz,0)));

	while(myqueue.size()!=0)
	{		
		quartet = myqueue.front();
		myqueue.pop();
		
		cix = quartet.first.first;
		ciy = quartet.first.second;
		ciz = quartet.second.first;
		dist = quartet.second.second;

		if ((dist+1)>maxMoves)
		{
/*			FILE* fp;
			if (debug)	fp = fopen("debugFlood.txt","w");
			printf("\n Target not found, too many moves");
			// restore old status
			for (unsigned int ll=0;ll<visitedPoints.size();ll++)
			{
				pair< int, pair <int,int>> triplet = visitedPoints[ll];
				//STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY) = -STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY);
				if (debug) fprintf(fp,"\n %f %f %f",X[triplet.first],Y[triplet.second.first],Z[triplet.second.second]);
			}
			if (debug) fclose(fp);*/
			maxMoves = dist;
			return false;
		}

		// target obtained
		if ((cix == target[0]) && (ciy == target[1]) && (ciz == target[2]))
		{
/*			printf("\n Target acquired!");
			FILE* fp;
			if (debug)	fp = fopen("debugFlood.txt","w");

			// restore old status
			for (unsigned int ll=0;ll<visitedPoints.size();ll++)
			{
				pair< int, pair <int,int>> triplet = visitedPoints[ll];
				//STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY) = -STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY);
				if (debug) fprintf(fp,"\n %f %f %f",X[triplet.first],Y[triplet.second.first],Z[triplet.second.second]);
			}
			if (debug) fclose(fp);*/
			maxMoves = dist;
			return true;
		}

		// already visited?
		if (visited.at(cix,ciy,ciz)==1)
			continue;
		else
			visited.set(cix,ciy,ciz,1);

		// can still move, go on
		int newDist = dist+1;
		bool notVisited,inner,isAccessible;

		double p[3];

		if (cix+1<NX)
		{
			p[0] = X[cix+1];
			p[1] = Y[ciy];
			p[2] = Z[ciz];

			if (visited.at(cix+1,ciy,ciz)==0)
				notVisited = true;
			else
				notVisited = false;

			if (notVisited)
			{				
				//if (TEMP_STATUSMAP2(cix+1,ciy,ciz,NX,NY)==STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix+1,ciy,ciz,NX,NY)==STATUS_POINT_OUT)
				const short ts = read3DVector<short>(tempStatus2,cix+1,ciy,ciz,NX,NY,NZ);
				if (ts==STATUS_POINT_TEMPORARY_OUT || ts==STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				//inner = TEMP_STATUSMAP(cix+1,ciy,ciz,NX,NY)==STATUS_POINT_INSIDE || STATUSMAP(cix+1,ciy,ciz,NX,NY)>=STATUS_FIRST_CAV;
				
				const int tss = read3DVector<short>(status,cix+1,ciy,ciz,NX,NY,NZ);
				inner = ((read3DVector<short>(tempStatus,cix+1,ciy,ciz,NX,NY,NZ)==STATUS_POINT_INSIDE) || (tss)>=STATUS_FIRST_CAV);
				
				if (isAccessible && inner)
				{
					quartet.first.first = cix+1;
					quartet.first.second = ciy;
					quartet.second.first = ciz;

					// if in cavity/pocket does not count the move, the flow is still trying
					// to percolate through the linking surf but still the linking surf path connecting the two cavities
					// has been not found
					//if (STATUSMAP(cix+1,ciy,ciz,NX,NY)>=STATUS_FIRST_CAV)
					if (tss>=STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
						
					myqueue.push(quartet);					
				}
			}
		}

		if (ciy+1<NY)
		{
			p[0] = X[cix];
			p[1] = Y[ciy+1];
			p[2] = Z[ciz];

			if (visited.at(cix,ciy+1,ciz)==0)
				notVisited = true;
			else
				notVisited = false;

			if (notVisited)
			{				
				//if (TEMP_STATUSMAP2(cix,ciy+1,ciz,NX,NY)==STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix,ciy+1,ciz,NX,NY)==STATUS_POINT_OUT)
				const short ts = read3DVector<short>(tempStatus2,cix,ciy+1,ciz,NX,NY,NZ);
				if (ts==STATUS_POINT_TEMPORARY_OUT || ts==STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				//inner = TEMP_STATUSMAP(cix,ciy+1,ciz,NX,NY)==STATUS_POINT_INSIDE || STATUSMAP(cix,ciy+1,ciz,NX,NY)>=STATUS_FIRST_CAV;

				const int tss = read3DVector<short>(status,cix,ciy+1,ciz,NX,NY,NZ);
				inner = ((read3DVector<short>(tempStatus,cix,ciy+1,ciz,NX,NY,NZ)==STATUS_POINT_INSIDE) || (tss)>=STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy+1;
					quartet.second.first = ciz;

					//if (STATUSMAP(cix,ciy+1,ciz,NX,NY)>=STATUS_FIRST_CAV)
					if (tss>=STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);					
				}
			}
		}

		if (ciz+1<NZ)
		{
			p[0] = X[cix];
			p[1] = Y[ciy];
			p[2] = Z[ciz+1];
			
			if (visited.at(cix,ciy,ciz+1)==0)
				notVisited = true;
			else
				notVisited = false;

			if (notVisited)
			{
				//if (TEMP_STATUSMAP2(cix,ciy,ciz+1,NX,NY)==STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix,ciy,ciz+1,NX,NY)==STATUS_POINT_OUT)
				const short ts = read3DVector<short>(tempStatus2,cix,ciy,ciz+1,NX,NY,NZ);
				if (ts==STATUS_POINT_TEMPORARY_OUT || ts==STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				//inner = TEMP_STATUSMAP(cix,ciy,ciz+1,NX,NY)==STATUS_POINT_INSIDE || STATUSMAP(cix,ciy,ciz+1,NX,NY)>=STATUS_FIRST_CAV;

				const int tss = read3DVector<short>(status,cix,ciy,ciz+1,NX,NY,NZ);
				inner = ((read3DVector<short>(tempStatus,cix,ciy,ciz+1,NX,NY,NZ)==STATUS_POINT_INSIDE) || (tss)>=STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy;
					quartet.second.first = ciz+1;

					//if (STATUSMAP(cix,ciy,ciz+1,NX,NY)>=STATUS_FIRST_CAV)
					if (tss>=STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);					
				}
			}
		}
		

		if (cix-1>=0)
		{
			p[0] = X[cix-1];
			p[1] = Y[ciy];
			p[2] = Z[ciz];
			
			if (visited.at(cix-1,ciy,ciz)==0)
				notVisited = true;
			else
				notVisited = false;

			if (notVisited)
			{
				//if (TEMP_STATUSMAP2(cix-1,ciy,ciz,NX,NY)==STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix-1,ciy,ciz,NX,NY)==STATUS_POINT_OUT)
				const short ts = read3DVector<short>(tempStatus2,cix-1,ciy,ciz,NX,NY,NZ);
				if (ts==STATUS_POINT_TEMPORARY_OUT || ts==STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				//inner = TEMP_STATUSMAP(cix-1,ciy,ciz,NX,NY)==STATUS_POINT_INSIDE || STATUSMAP(cix-1,ciy,ciz,NX,NY)>=STATUS_FIRST_CAV;
				const int tss = read3DVector<short>(status,cix-1,ciy,ciz,NX,NY,NZ);
				inner = ((read3DVector<short>(tempStatus,cix-1,ciy,ciz,NX,NY,NZ)==STATUS_POINT_INSIDE) || (tss)>=STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix-1;
					quartet.first.second = ciy;
					quartet.second.first = ciz;

					//if (STATUSMAP(cix-1,ciy,ciz,NX,NY)>=STATUS_FIRST_CAV)
					if (tss>=STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);					
				}
			}
		}

		if (ciy-1>=0)
		{
			p[0] = X[cix];
			p[1] = Y[ciy-1];
			p[2] = Z[ciz];
		
			if (visited.at(cix,ciy-1,ciz)==0)
				notVisited = true;
			else
				notVisited = false;


			if (notVisited)
			{
				//if (TEMP_STATUSMAP2(cix,ciy-1,ciz,NX,NY)==STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix,ciy-1,ciz,NX,NY)==STATUS_POINT_OUT)
				const short ts = read3DVector<short>(tempStatus2,cix,ciy-1,ciz,NX,NY,NZ);
				if (ts==STATUS_POINT_TEMPORARY_OUT || ts==STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;
								
				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				//inner = TEMP_STATUSMAP(cix,ciy-1,ciz,NX,NY)==STATUS_POINT_INSIDE || STATUSMAP(cix,ciy-1,ciz,NX,NY)>=STATUS_FIRST_CAV;

				const int tss = read3DVector<short>(status,cix,ciy-1,ciz,NX,NY,NZ);
				inner = ((read3DVector<short>(tempStatus,cix,ciy-1,ciz,NX,NY,NZ)==STATUS_POINT_INSIDE) || (tss)>=STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy-1;
					quartet.second.first = ciz;

					//if (STATUSMAP(cix,ciy-1,ciz,NX,NY)>=STATUS_FIRST_CAV)
					if (tss>=STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;
					myqueue.push(quartet);					
				}
			}
		}

		if (ciz-1>=0)
		{
			p[0] = X[cix];
			p[1] = Y[ciy];
			p[2] = Z[ciz-1];

			//notVisited = (STATUSMAP(cix,ciy,ciz-1,NX,NY)!=-STATUS_POINT_INSIDE && !(STATUSMAP(cix,ciy,ciz-1,NX,NY)<=-STATUS_FIRST_CAV));

			if (visited.at(cix,ciy,ciz-1)==0)
				notVisited = true;
			else
				notVisited = false;

			if (notVisited)
			{
				//if (TEMP_STATUSMAP2(cix,ciy,ciz-1,NX,NY)==STATUS_POINT_TEMPORARY_OUT || TEMP_STATUSMAP2(cix,ciy,ciz-1,NX,NY)==STATUS_POINT_OUT)
				const short ts = read3DVector<short>(tempStatus2,cix,ciy,ciz-1,NX,NY,NZ);
				if (ts==STATUS_POINT_TEMPORARY_OUT || ts==STATUS_POINT_OUT)
					isAccessible = true;
				else
					isAccessible = false;

				// if inside 1.4 rp surf or in cavity of current surf then it is ok
				//inner = TEMP_STATUSMAP(cix,ciy,ciz-1,NX,NY)==STATUS_POINT_INSIDE || STATUSMAP(cix,ciy,ciz-1,NX,NY)>=STATUS_FIRST_CAV;

				const int tss = read3DVector<short>(status,cix,ciy,ciz-1,NX,NY,NZ);
				inner = ((read3DVector<short>(tempStatus,cix,ciy,ciz-1,NX,NY,NZ)==STATUS_POINT_INSIDE) || (tss)>=STATUS_FIRST_CAV);

				if (isAccessible && inner)
				{
					quartet.first.first = cix;
					quartet.first.second = ciy;
					quartet.second.first = ciz-1;

					//if (STATUSMAP(cix,ciy,ciz-1,NX,NY)>=STATUS_FIRST_CAV)
					if (tss>=STATUS_FIRST_CAV)
						quartet.second.second = dist;
					else
						quartet.second.second = newDist;

					myqueue.push(quartet);
				}
			}
		}
	}
/*
	FILE* fp;
	if (debug)	fp = fopen("debugFlood.txt","w");

	// restore old status
	for (unsigned int ll=0;ll<visitedPoints.size();ll++)
	{
		pair< int, pair <int,int>> triplet = visitedPoints[ll];
		//STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY) = -STATUSMAP((triplet.first),(triplet.second.first),(triplet.second.second),NX,NY);
		if (debug) fprintf(fp,"\n %f %f %f",X[triplet.first],Y[triplet.second.first],Z[triplet.second.second]);
	}
	if (debug) fclose(fp);
	printf("\n Target not found!");*/
	maxMoves = dist;	
	return false;
}

bool Surface::isCompletelyOut(double* pos)
{	
	if (delphi->status == NULL)
	{
		cout << endl << ERR << "Cannot compute if a point is completely out without the status map";
		exit(-1);
	}

	int nx = delphi->nx;
	int ny = delphi->ny;
	int nz = delphi->nz;

	int ix = (int)rintp((pos[0]-delphi->xmin)/delphi->side);
	int iy = (int)rintp((pos[1]-delphi->ymin)/delphi->side);
	int iz = (int)rintp((pos[2]-delphi->zmin)/delphi->side);

	short test;
		
	for (int k=0;k<SHIFT_MAP;k++)
	{
		int cx = ix+shift_map[k][0];
		int cy = iy+shift_map[k][1];
		int cz = iz+shift_map[k][2];

		if (cx>(nx-1) || cy>(ny-1) || cz>(nz-1) || cx<0 || cy<0 || cz<0)
			continue;

		test = read3DVector<short>(delphi->status,cx,cy,cz,nx,ny,nz);

		// ok
		if (test == STATUS_POINT_TEMPORARY_OUT || test == STATUS_POINT_OUT)
			continue;
		// at least one is not completely out
		else
			return false;
	}
	return true;
}

void Surface::triangulationPointsAreCompletelyOut(Surface* surf,vector<bool>& results)
{	
	for (unsigned int i=0;i<vertList.size();i++)
	{
		results.push_back(surf->isCompletelyOut(vertList[i]));
	}
}

void Surface::saveSelectedPoints(vector<bool>& selection,char* file,char* file2)
{
	FILE* fp = fopen(file,"w");
	FILE* fp2 = NULL;
	
	if (normalsList.size()!=0)
		fp2 = fopen(file2,"w");

	int count = 0;
	for (unsigned int i=0;i<selection.size();i++)
	{
		if (selection[i])
			count ++;
	}
	// write header only for VMD file (without normals)
	fprintf(fp,"%d\n",count);
	fprintf(fp,"pocket_entrance_as_a_set_of_dummy_atoms_by_NanoShaper\n");

	for (unsigned int i=0;i<selection.size();i++)
	{
		if (selection[i])
		{
			double* v = vertList[i];
			double* n = NULL;

			if (fp2!=NULL)
				n = normalsList[i];

			fprintf(fp,"C %f %f %f\n",v[0],v[1],v[2]);

			if (n!=NULL)
				fprintf(fp2,"%f %f %f %f %f %f\n",v[0],v[1],v[2],n[0],n[1],n[2]);			
		}
	}
	fclose(fp);

	if (fp2!=NULL)
		fclose(fp2);
}

bool Surface::vdwAccessible(double* pos,int& winner)
{
	double minDist=INFINITY;
	winner = -1;
	
	int ix = (int)rintp((pos[0]-gxmin)/gside);
	int iy = (int)rintp((pos[1]-gymin)/gside);
	int iz = (int)rintp((pos[2]-gzmin)/gside);
	// get the nearest atom and says if it is in or out	
	// the signed distance from the point p is ||p-c||^2-r^2 where c is the center
	// of the atom and r is the radius. The minimum signed distance wins and if this
	// negative we are inside (false) and if it is positive we are outside (true)
	for (int k=0;k<SHIFT_MAP;k++)
	{
		unsigned int cx = ix+shift_map[k][0];
		unsigned int cy = iy+shift_map[k][1];
		unsigned int cz = iz+shift_map[k][2];

		// multidielectric map is square 
		if (cx>(ggrid-1) || cy>(ggrid-1) || cz>(ggrid-1) || cx<0 || cy<0 || cz<0)
			continue;

		//int max_ind = GRID_MULTI_MAP(cx,cy,cz,0,ggrid,ggrid,ggrid);
		int max_ind = read4DVector<int>(gridMultiMap,cx,cy,cz,0,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

		for (int j=0;j<max_ind;j++)
		{
			//int atom_index = GRID_MULTI_MAP(cx,cy,cz,j+1,ggrid,ggrid,ggrid);
			int atom_index = read4DVector<int>(gridMultiMap,cx,cy,cz,j+1,ggrid,ggrid,ggrid,MAX_ATOMS_MULTI_GRID);

			double signed_dist = 0;
			DIST2(signed_dist,delphi->atoms[atom_index]->pos,pos)
			double rad = delphi->atoms[atom_index]->radius;
			signed_dist-=(rad*rad);
			if (signed_dist<minDist)
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

	if (minDist<0) return false;
	else return true;
}


double Surface::saveTriSubSet(char* triSubset,vector<bool>& results,bool revert)
{
	double area = 0;
	
	// previous vertex index, new vertex index
	// -1 if the vertex is removed
	map<int,int> indicesMap;

	// reduced list of triangles, vertices and normals
	vector<int*> triListRed;	
	vector<double*> vertListRed;
	vector<double*> normalsListRed;

	if (triList.size()==0)
	{
		cout << ERR << "Cannot save a reduced set of triangles if triangulation is absent";
		return 0.;
	}

	if (results.size()!=vertList.size())
	{
		cout << ERR << "Cannot save a reduced set of triangles if number of vertices is unconsistent between current mesh and vertices flags";
		return 0.;
	}

	for (unsigned int i=0;i<results.size();i++)
	{
		if (results[i])
		{
			vertListRed.push_back(vertList[i]);				
			indicesMap.insert(std::pair<int,int>(i,(int)vertListRed.size()-1));
			if (normalsList.size()!=0)
				normalsListRed.push_back(normalsList[i]);
		}
		else
			indicesMap.insert(std::pair<int,int>(i,-1));
	}

	for (unsigned int i=0;i<triList.size();i++)
	{
		int* tri = triList[i];
		map<int,int>::iterator it;
		
		it=indicesMap.find(tri[0]);
		int i1 = (*it).second;
		if (i1==-1)
			continue;
		
		it=indicesMap.find(tri[1]);
		int i2 = (*it).second;
		if (i2==-1)
			continue;

		it=indicesMap.find(tri[2]);
		int i3 = (*it).second;
		if (i3==-1)
			continue;

		int* newTri = allocateVector<int>(3);
		newTri[0]=i1;
		newTri[1]=i2;
		newTri[2]=i3;

		double a =0,b=0,c=0;

		double* v0 = vertListRed[newTri[0]];
		double* v1 = vertListRed[newTri[1]];
		double* v2 = vertListRed[newTri[2]];

		DIST(a,v0,v1);
		DIST(b,v0,v2);
		DIST(c,v1,v2);

		double ttt = (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c);

		if (ttt>0)						
			area += 0.25*sqrt(ttt);

		triListRed.push_back(newTri);
	}

	int format = deduceFormat();
	saveMesh(format,revert,triSubset,vertListRed,triListRed,normalsListRed);
	
	// clean up
	for (unsigned int i=0;i<triListRed.size();i++)
		deleteVector(triListRed[i]);	

	return area;
}