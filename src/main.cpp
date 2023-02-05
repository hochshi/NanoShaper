
#include "globals.h"

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif

#include "SurfaceFactory.h"
#include "Surface.h"
#include "tools.h"
#include "DelphiShared.h"
#include "ConfigFile.h"

// number of cores
int num_cores;

// have to use a pointer because Swig will try to make a copy and fstream does
// not have a copy semantic
fstream* errorStream;
fstream* internals;
Configuration conf;

ConfigFile* init(string argv);
void dispose(ConfigFile* cf);
void stopDebug();
void restartDebug();
void cite();
void normalMode(Surface* surf,DelPhiShared* dg);
void pocketMode(bool hasAtomInfo,ConfigFile* cf);

#ifndef DELPHI_BIND
/** @mainpage 
	@brief 
	NanoShaper is a framework to analyze an arbitrary surface; a particular attention is given to molecular surfaces.
	NanoShaper presents several features: 
	<ul>
	<li> Loading of a triangulated mesh in OFF, PLY or MSMS format 
	<li> Analytical build and triangulation of the Skin surface, the Blobby surface and the Connolly surface
	<li> Colour with in/out info of a 3D grid
	<li> Pockets and cavities detection.
	<li> Easy to plug a new surface: it is sufficient to expand Surface class, no other modifications are required,
	just writing the .h and .cpp. The surface registration is automatically managed by a template based mechanism.
	<li> It can be compiled as a Standalone, python lib or as a DelPhi plug-in.
	</ul>
	
	History:
	
	<ul>
	<li> Version 0.3.1:
	<ul>
	<li> Much faster Skin Surface build-up. Acceleration grid 2d also for the mesh ray casting.
	</ul>
	
	<li> Version 0.3.2:
	<ul>
	<li> bug fix on Connolly Surface.
	</ul>
	
	<li> Version 0.4:
	<ul>
	<li> Using octrees to reduce memory requirements for triangulation and using an improved multi-threading strategy. 
	<li> Atom based Multi-dielectric is supported (DelPhi).
	<li> Stern Layer (DelPhi).
	<li> An efficient, ray-casting based algorithm to convert a set of surface samples to a set of balls whose envelope approximatess the surface.
	<li> Atoms in the neighboroud of a detected cavity are returned. 
	</ul>

	<li> Version 0.5:
	<ul>
	<li> "-=" Filtered Operator to support difference map to detect pockets not detectable by cavity detector 
	<li> Operative_Mode keyword. Values = normal, difference, membrane 
	<li> Introduced experimental Fast Van Der Waals surface. Analytical for FD PDE, not analytical triangulation 
	<li> Introduced experimental Coulombic surface 
	<li> From DelPhi atinf and atsurf are imported. atsurf is filled  (DelPhi)
	<li> Fast projection flag to switch between fast (approximate) and slow (accurate) SkinSurface bgp projections 
	<li> Output for pocket detection in ProShape style added when used with DelPhi 
	<li> Improved threads load balancing
	<li> Every vertex has its own reference atoms (OFF+A format defined) 
	<li> Faster flood fill algorithm 
	<li> Analytical normals computation or approximation by triangulation 
	<li> OFF+N, OFF+N+A, MSMS formats 
	<li> Vertex normals are computed, analytically if possible.
	</ul>
	
	<li> Version 0.6:
	<ul>
	<li> Refactoring using metaprogramming to instantiate surfaces. 
		To introduce a new surface, just define that class and register it using the registration macro.
		Now surface constructors access directly the configuration file structure. Any new surface thus
		can define each own keywords.
	</ul>
	
	<li> Version 0.7
	<ul>
	<li> Parallel Marching Cubes: 2 steps MC, cache friendly and trivially parallel 
	<li> Dangling vertices now are removed a priori. Dangling vertices are those that were sampled, gave
		an in/out map that is in contrast with another ray for numerical reasons. If this happen these
		vertices, after ray tracing, are identified and removed before marching cubes. In this way
		the exact number of MC vertices is known a priori.
		Another pathological situation that may happen on previous versions is that two vertices are false
		vertices on the same edge (let's say an high frequency) now these vertices are removed.
		The previous version of NanoShaper had this added spurious vertices in memory; however they are not
		present in the mesh, so the previous meshes are still correct.
		Moreover now a slightly faster way to color the grid in ray tracing is defined. 
	<li> If a normal cannot be computed analytically due to any kind of problem, than this is approximated
		in a subsequent step such that any file with normals will have all the normals computed. Most
		of them will be analytical (usually less than 1 over 10000 are not analytical). 
	<li> High accuracy timing boost::chrono is optionally available 
	<li> Bug fix in blobby surface due to surface factory refactoring 
	<li> Introduced the ExampleSurface to let the user better understand how to add a surface
	<li> Using templates instead of macros for memory access/write, check bounds and added p=NULL semantic on deletion.
		from now on vectors of primitive types should be allocated, accessed, written, deleted with the given template functions
	<li> Refactoring of main method. 
	<li> Introduced pure virtuals init/clear methods to force the Developer to write clean C++ code for the Surfaces
	<li> Optimization on the floodfill. Ping pong to get as much as possible cache friendly.
	<li> Bug fix of ConnollySurface, now Probe_Radius is not ignored.
	</ul>

	<li> Version 0.7.4
	<ul>
	<li> Introduced Max_Atoms_Multi_Grid keyword to control atomsmap max size.
	<li> Introduced an algorithm to detected the entrance of a pocket. The entrance is saved as a set of points in entranceX.xyz
	and entranceX.xyzn where the first file can be visualized in vmd.
	<li> If the Connolly surface is computed a file name called exposed.xyz saves as Carbons in xyz format (VMD) the list of exposed
	atoms. Additionally a file named exposedIndices.txt is the list of 0-based indices of exposed atoms.
	<li> The area of pockets, excluding the mouth is returned. Triangulated files are named cav_triX_body.TRI where is triangulation
	extension and X is the index of the pocket
	<li> Bugfix for MSMS save for blobby

	</ul>

	<li> Version 0.7.5
	<ul>
	<li> Bug fix: Save_Status_map flag was not read. Fixed

	</ul>


	</ul>

	@author Sergio Decherchi
	@date 25-07-2014
	@version 0.7.5
	*/
int main(int argc,char* argv[])
{
	#ifdef  DBGMEM_CRT
//		 _crtBreakAlloc = 16222;
		_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );		
	#endif	

	int numargs = argc-1;
	char confFile[BUFLEN];

	if (!numargs)
		strcpy(confFile,"surfaceConfiguration.prm");
	else
		strcpy(confFile,argv[1]);

	if (numargs>1)
		cout << endl << WARN << "Ignoring additional non required parameters";

	// check configuration consistency, init error stream, get configuration
	ConfigFile* cf = init(string(confFile));
	
	// operative modes are the NanoShaper lib utilization modes.

	// just build the surface
	if (!conf.operativeMode.compare("normal"))
	{		
		// Set up DelPhi-like environment
		DelPhiShared* dg = new DelPhiShared(conf.scale,conf.perfill,conf.molFile,conf.buildEpsmaps,conf.buildStatus,conf.multi_diel);
		// Get surface
		Surface* surf = surfaceFactory().create(cf,dg);
		normalMode(surf,dg);		
		cout << endl << INFO << "Cleaning memory...";
		cout.flush();
		delete surf;		
		delete dg;
		cout << "ok!";
	}
	// detect pockets
	else if (!conf.operativeMode.compare("pockets"))
	{		
		pocketMode(false,cf);						
	}
	else
	{
		cout << endl << INFO << "Unknown operative mode";
		cout << endl;
		return -1;
	}

	cite();	
	dispose(cf);
	cout << endl << endl ;

	// Memory leak detection 
	#ifdef DBGMEM_CRT
		_CrtDumpMemoryLeaks();
	#endif

	return 0;
}

#endif

	#ifdef DELPHI_BIND
	
	/** This is the DelPhi Fortran-C binding. Returns the number of bgp */
	// check visual, borland c++, xlc, watcom, for windows
	#if  (defined _WIN32) || (defined __WIN32__) || (defined __TOS_WIN__) || (defined __WINDOWS__) 
			extern "C" __declspec( dllexport )		
	// linux,unix,apple,android 
	#elif (defined __linux) || (defined __unix) || (defined macintosh) || (defined Macintosh) || (defined __APPLE__ && defined __MACH__) || (defined __ANDROID__)
			extern "C" 		
	#endif
					
	void epsmakemodule(double xmin,double ymin,double zmin,
					  double xmax,double ymax,double zmax,
					  double c1,double c2,double c3,	double rmax,double perf,							 
					  int* i_epsmap,int igrid,double scale,double* i_scspos,
				      double* i_scsnor,bool* i_idebmap,int* i_ibgp,int inside,int* ibnum,int natom,
					  double* xn1,double* rad, int* d,int maxbgp,double probe,double exrad,char* atinf,int* atsurf)
	{
		
		// check configuration consistency, init error stream and get configuration
		ConfigFile* cf = init(string("surfaceConfiguration.prm"));

		// override settings coming from conf with that coming from DelPhi
		if (cf->keyExists("Grid_scale"))
			cf->remove("Grid_scale");
		cf->add<double>("Grid_scale",scale);

		if (cf->keyExists("Grid_perfil"))
			cf->remove("Grid_perfil");
		cf->add<double>("Grid_perfil",perf);

		conf.scale = scale;
		conf.perfill = perf;

		ConfigFile* cf2;		
		// get data from custom configuration file
		try
		{
			cf2 = new ConfigFile( "custom.prm" );
		}
		catch(...)
		{
			cout << endl << ERR << "Cannot read custom.prm";
			cout.flush();
			exit(-1);			
		}		
		string mode = cf2->read<string>( "Surface", "ses" );
		delete cf2;

		if (cf->keyExists("Surface"))
			cf->remove("Surface");

		cf->add<string>("Surface",mode);

		if (conf.operativeMode=="normal")
		{			
			cout << endl << INFO << "Binding with DelPhi..";		
			
			// Set up delphi environment
			DelPhiShared* dg = new DelPhiShared();
			dg->DelPhiBinding(	xmin,ymin,zmin,xmax,ymax,zmax,
								c1,c2,c3,rmax,perf,
								i_epsmap,igrid,scale,
								i_scspos,i_scsnor,i_idebmap,
								i_ibgp,maxbgp,conf.buildStatus,atsurf);

			// load atoms info if present
			dg->loadAtoms(natom,xn1,rad,NULL,NULL,atinf);

			// here only means populate epsilon map
			dg->buildEpsmap(true);
			
			cout << endl << INFO << "DelPhi grid is " << igrid;
			
			Surface* surf = surfaceFactory().create(cf,dg);
			surf->setProjBGP(true);
			surf->setInsideCode(inside);

			if (exrad>0)
			{
				cout << endl << INFO << "Setting Stern Layer -> " << exrad << " Angstrom";
				surf->setSternLayer(exrad);
			}

			// normal protocol
			normalMode(surf,dg);
						
			// avoid to destroy things that now belong to DelPhi
			dg->finalizeBinding(ibnum);
			
			cout << endl << INFO << "Cleaning memory...";	
			cout.flush();

			delete surf;
			delete dg;
		
			dispose(cf);

			cout << "ok!";
			cout << endl << endl ;
			return;
		}
	
	// In this mode use the difference on grid operator to infer pockets and cavities (if requested)
	else if (!conf.operativeMode.compare("pockets"))
	{	
		DelPhiShared::parseAtomInfo(atinf,conf.molFile,natom,xn1,rad);		
		
		pocketMode(true,cf);
				
		dispose(cf);
		
		cout << endl << INFO << "Brute force exiting to avoid DelPhi solver";
		cout << endl << INFO << "Closing " << PROGNAME << "\n";
		exit(-1);
	}
	else
	{
		cout << endl << INFO << "Unknown operative mode: " << conf.operativeMode;
		exit(-1);
	}

	cite();		
	dispose(cf);	
	cout << endl << endl;	
	}
	
	#endif

// init streams, check configuration file for errors and read variables
ConfigFile* init(string confFile)
{
	cout << endl << endl << INFO << "Starting " << PROGNAME << " " << VERSION;
			
	ConfigFile* cf;

	// get data from configuration file
	try
	{
		cf = new ConfigFile(confFile.c_str());
	}
	catch(...)
	{
		cout << endl << ERR << "Cannot read " << confFile;
		exit(-1);			
	}
							
	errorStream = new fstream("stderror.txt",fstream::out);

	conf.buildEpsmaps = cf->read<bool>( "Build_epsilon_maps", false );
	conf.saveEpsmaps = cf->read<bool>( "Save_eps_maps", false );
	
	
	conf.saveStatusMap = cf->read<bool>( "Save_Status_map", false );

	conf.buildStatus = cf->read<bool>( "Build_status_map", false );
	conf.saveIdebmap = cf->read<bool>( "Save_ideb_map", false );
	conf.fillCavities = cf->read<bool>( "Cavity_Detection_Filling", false );

	conf.projBGP = cf->read<bool>( "Project_boundary_grid_points", false );

	conf.accTri = cf->read<bool>( "Accurate_Triangulation", false );
	conf.tri = cf->read<bool>( "Triangulation", false );
			
	conf.operativeMode = cf->read<string>( "Operative_Mode","normal");

	bool dbg = cf->read<bool>( "Debug_Internals", false );
	
	if (dbg)
		internals = new fstream("internals.txt",fstream::out);
	else
		internals = NULL;

	if (cf!=NULL)
	{
		if (!conf.buildEpsmaps && conf.saveEpsmaps)
		{
			cout << endl << ERR << "Asked to save epsmap without builiding it";
			cout << endl << REMARK << "Please set Build_epsilon_maps = true";
			cout << endl;
			exit(-1);
		}

		if (!conf.buildEpsmaps && conf.projBGP)
		{
			cout << endl << ERR << "Cannot project boundary grid points without an epsilon map.";
			cout << endl << REMARK << "Please set Build_epsilon_maps = true";
			cout << endl;
			exit(-1);	
		}

		if (!conf.accTri && !conf.buildStatus && conf.tri)
		{
			// status map is needed to deduced in/out vertices			
			cout << endl << ERR << "If non analytical triangulation is enabled status map is needed.";
			cout << endl << REMARK << "Please set Build_status_map = true";
			cout << endl;
			exit(-1);
		}

		if (conf.fillCavities && !conf.buildStatus)
		{
			// status map is needed to search cavities
			cout << endl << ERR << "If cavity detection is enabled status map is needed.";
			cout << endl << REMARK << "Please set Build_status_map = true";
			cout << endl;
			exit(-1);
		}

		if (conf.saveIdebmap && !conf.buildEpsmaps)
		{
			cout << endl << ERR << "Idebmap is computed only if epsilon map is enabled";
			cout << endl << REMARK << "Please set Build_epsilon_maps = true";
			cout << endl;
			exit(-1);	
		}

		if (!conf.operativeMode.compare("pockets") && !conf.buildStatus)
		{
			cout << endl << WARN << "Cannot do pocket detection without status map";
			cout << endl << REMARK << "Please set Build_status_map = true";
			cout << endl;
			exit(-1);	
		}
	}

	conf.cavVol = cf->read<double>( "Conditional_Volume_Filling_Value", 11.4 );
	conf.numMol = cf->read<int>( "Num_Wat_Pocket",2);
	
	// grid (DelPhi) params
	conf.scale = cf->read<double>( "Grid_scale", 2.0 );
	conf.perfill = cf->read<double>( "Grid_perfil", 80.0 );
	conf.molFile = cf->read<string>( "XYZR_FileName","temp.xyzr");
	conf.multi_diel = cf->read<bool>( "Multi_Dielectric",false);

	// tri
	conf.smoothing = cf->read<bool>( "Smooth_Mesh", false );	
	conf.tri2balls = cf->read<bool>( "Tri2Balls", false );
	
	// save data
	conf.saveEpsmaps = cf->read<bool>( "Save_eps_maps", false );
	conf.saveBgps = cf->read<bool>( "Save_bgps", false );	
	conf.saveCavities = cf->read<bool>( "Save_Cavities", false );

	// globals	
	conf.sysName = cf->read<string>( "Sys_Name","mol");
	conf.numthd = cf->read<int>( "Number_thread", -1 );	
	conf.printAvailSurf = cf->read<bool>( "Print_Available_Surfaces", false );	
	conf.currentSeed = cf->read<int>( "Seed", 1 );

	// pocket detection
	conf.cavAndPockets = cf->read<bool>( "Pockets_And_Cavities",true);
	conf.linkPockets = cf->read<bool>( "Link_Pockets",false);	
	conf.pocketRadiusBig = cf->read<double>( "Pocket_Radius_Big",3.0);		
	conf.pocketRadiusSmall = cf->read<double>( "Pocket_Radius_Small",1.4);		
	conf.pocketRadiusLink = cf->read<double>( "Pocket_Radius_Link",1.0);	
			
	num_cores = conf.numthd;

	return cf;
}

void dispose(ConfigFile* cf)
{
	delete cf;		
	errorStream->close();
	delete errorStream;
	if (internals!=NULL)
	{
		internals->close();
		delete internals;
	}
}

void stopDebug()
{
	if (internals!=NULL)
	{
		internals->close();
		internals = NULL;
	}
}

void restartDebug()
{
	if (internals==NULL)
		internals = new	fstream("internals.txt",fstream::app);	
}


void cite()
{
	cout << endl;
	cout << endl << INFO << "If you use NanoShaper please cite these works:";
	cout << endl << CITE << "\tS. Decherchi, W. Rocchia, \"A general and Robust Ray-Casting-Based Algorithm for Triangulating Surfaces at the Nanoscale\"; PlosOne";
	cout << endl << CITE << "\tlink: http://www.plosone.org/article/metrics/info%3Adoi%2F10.1371%2Fjournal.pone.0059744";
	cout << endl;
}	


/** the set of operations in the usual mode of usage. This function is not responsible for
Surface or grid memory. The caller is the responsible.*/
void normalMode(Surface* surf,DelPhiShared* dg)
{
	if (conf.printAvailSurf)
		surfaceFactory().print();

	Timer chrono;
	chrono.start();

	char refName[100];
	strcpy(refName,conf.sysName.c_str());

	// Pre-process surface
	bool outsurf = surf->build();

	if (!outsurf)
	{	
		cout << endl << ERR << "Surface construction failed!";
		exit(-1);
	}

	// Build DelPhi stuff
	surf->getSurf(conf.fillCavities,conf.cavVol);		
		
	double duration = chrono.stop();

	cout << endl << INFO << "Surface computation time.. " << duration << " [s]"; 	
	cout << endl << INFO << "Estimated volume " << surf->getVolume() << " [A^3]"; 

	if (internals!=NULL)
		(*internals) << endl << "volume " << surf->getVolume();
			
	if (conf.tri)
	{
		cout << endl << INFO << "Triangulating Surface...";				
			
		Timer chrono;
		chrono.start();

		surf->triangulateSurface();		

		if (internals!=NULL)
		{
			(*internals) << endl << "area " << surf->getArea();
			(*internals) << endl << "nv " << surf->getNumVertices();
			(*internals) << endl << "nt " << surf->getNumTriangles();
		}

		if (conf.smoothing)
		{
			cout << endl << INFO << "Smoothing surface...";
			surf->smoothSurface();			
		}

		double duration = chrono.stop();
		cout << "ok!";
			
		cout << endl << INFO << "Total Triangulation time " << duration << " [s]";
	}

	if (conf.tri2balls)
	{
		cout << endl << INFO << "Converting triangulation to balls...";				
		surf->tri2Balls();
		cout << "ok!";
	}

	if (conf.saveEpsmaps)
	{
		cout << endl << INFO << "Saving epsmaps...";
		// Save epsmap
		dg->saveEpsMaps(refName);
		cout << "ok!";
	}

	if (conf.saveBgps)
	{
		cout << endl << INFO << "Saving bgpmap...";
		dg->saveBGP(refName);
		cout << "ok!";
	}

	if (conf.saveStatusMap)
	{
		cout << endl << INFO << "Saving statusmap and cavities...";
		dg->saveStatus(refName);			
		cout << "ok!";	
	}

	if (conf.saveCavities)
	{
		cout << endl << INFO << "Saving cavities...";
		dg->saveCavities(false);
		cout << "ok!";	
	}
		
	if (conf.saveIdebmap)
	{
		cout << endl << INFO << "Saving idebmap...";
		dg->saveIdebMap(refName);
		cout << "ok!";		
	}		
}

void pocketMode(bool hasAtomInfo,ConfigFile* cf)
{
	double duration = 0;

	bool localEpsMap = false;
	bool localStatusMap = true;
	bool localMulti = false;

	// only debug pockets/cavities
	if (conf.debug) stopDebug();

	Timer chrono;
	chrono.start();

	Surface *surf1,*surf2,*surf3;
	double areaSurf=0,volSurf=0;
						
	char mol[100]="temp.txt",refName[100]="mol";
	
	/////////////////////////////////////////////////////////////////////////////////////////
	// Set up Surface 1 (fat probe)
	cout << endl ;
	cout << endl << INFO << "Step 1 -> fat probe";	
	DelPhiShared* dg1 = new DelPhiShared(conf.scale,conf.perfill,conf.molFile,localEpsMap,localStatusMap,localMulti,hasAtomInfo);
	int natom = dg1->getNumAtoms();
	cf->remove(string("Surface"));
	cf->add<string>("Surface","ses");
	surf1 = surfaceFactory().create(cf,dg1);
	surf1->setInsideCode(5);	
	surf1->setProbeRadius(conf.pocketRadiusBig);
	surf1->setProjBGP(false);
	surf1->setKeepWellShapedCavities(false);

	// Pre-process surface
	bool outsurf = surf1->build();

	if (!outsurf)
	{	
		cout << endl << ERR << "Surface 1 construction failed!";
		exit(-1);
	}

	// fat connolly cancel each cavity (anyway they are smaller than they should be at the end)
	surf1->getSurf(true,INFINITY);

	//////////////////////////////////////////////////////////////////////////////////////////
	// the subsequent surfaces don't need to read atom info

	// Set up Surface 2 (regular probe 1.4)
	// Set up DelPhi grid  2 
	cout << endl ;
	cout << endl << INFO << "Step 2 -> small probe";	
	
	DelPhiShared* dg2 = new DelPhiShared(conf.scale,conf.perfill,conf.molFile,localEpsMap,localStatusMap,localMulti,false);
	surf2 = surfaceFactory().create(cf,dg2);
	surf2->setInsideCode(10);
	surf2->setProjBGP(false);	
	surf2->setProbeRadius(conf.pocketRadiusSmall);
	surf2->setKeepWellShapedCavities(false);

	/*
	//// to check percolation ///////////////////////
	DelPhiShared* dg2 = new DelPhiShared(conf.scale,conf.perfill,conf.mol,localStatusMap,localMulti,false);
	surf2 = surfaceFactory().create(cf,dg2);		
	surf2->inside = 10;			
	////////////////////////////////////////////////
	*/
	
	// Pre-process surface
	outsurf = surf2->build();

	if (!outsurf)
	{	
		cout << endl << ERR << "Surface 2 construction failed!";
		exit(-1);
	}

	// if cav and pockets together -> do not perform cavity detection (keep all)
	// if only pockets then perform cavity detection and remove all cavities
	// by default keep both cavities and pockets
	surf2->getSurf(!conf.cavAndPockets,INFINITY);		

	// if triangulation is enabled we triangulate the small probe surface, get surface and area
	if (conf.tri)
	{			
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
*/	//////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////		
	DelPhiShared* dg3 = NULL;
	if (conf.linkPockets)
	{
		// Set up Surface 3 (accessibility probe)			
		dg3 = new DelPhiShared(conf.scale,conf.perfill,conf.molFile,localEpsMap,localStatusMap,localMulti,false);
		surf3 = surfaceFactory().create(cf,dg3);
		surf2->setProjBGP(false);	
		surf2->setProbeRadius(conf.pocketRadiusSmall);
		surf2->setKeepWellShapedCavities(false);
		surf3->setInsideCode(15);

		// Pre-process surface
		outsurf = surf3->build();

		if (!outsurf)
		{	
			cout << endl << ERR << "Surface 3 construction failed!";
			exit(-1);
		}
		// keep original surface
		surf3->getSurf(false);
	}

	cout << endl ;
	cout << endl << INFO << "Step 3 -> differential map";	

	cout << endl << INFO << "Building pockets by difference map..." ;
	(*surf1)-=(*surf2);

	////////////////////// recover split cavities links ///////////////		
	if (conf.linkPockets)
	{
		int nr = 1;			
		double duration1 = 0;
			
		Timer chrono1;
		chrono1.start();

		cout << endl ;
		cout << endl << INFO << "Step 3.1 -> linking";	

		cout << endl << INFO << "Linking cavities/pockets...";
		cout.flush();
		while(nr!=0)
		{			
			// check cavities links, use the link status map as reference map
			// to check accessibility
			nr = surf1->linkCavities(dg2->status,dg3->status);		
			cout << endl << INFO << "Merged " << nr << " cavities";		
			cout.flush();
		}

		duration1 = chrono1.stop();
		cout << endl << INFO << "Diff. Step 4 " << duration1 << " [s]";
			
		delete surf3;
		delete dg3;
	}
	///////////////////////////////////////////////////////////////////

	cout << endl ;
	cout << endl << INFO << "Step 4 -> filtering, envelope building";	

	///////////////////// volume filter ///////////////////////////////
	surf1->fillCavities(11.4*conf.numMol,false);
	surf1->getCavitiesAtoms();		
	cout << endl << INFO << "Saving cavities info..";

	if (hasAtomInfo)
	{
		// save cavities in ProShape format
		dg1->saveCavities2(true,conf.sysName);
	}
	else
	{		
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
		
	if (conf.cavAndPockets)
	{
		cout << endl;
		// do cavity detection to recover which is cavity in the reference
		cout << endl << INFO << "Recovering cavity/pocket distinction..";
		cout.flush();
		surf2->getCavities();
		// mark which are pockets and which are cavities
		dg1->markPockets(dg2->status,isPocket);
	}

	cout << endl << INFO << "Build the envelope of each cavity/pocket";
	cout.flush();

	cf->remove("Surface");
	cf->add<string>("Surface","skin");

	// TODO make it an option
	bool saveEntranceInfo = true;

	for (int i=0;i<nc;i++)
	{
		char mol[100],tri_[100];
		sprintf(mol,"cav%d.txt",i);
		sprintf(tri_,"cav_tri%d",i);								
		DelPhiShared* dg_temp = new DelPhiShared(2.0,50,mol,false,false,false);
		// reset seed
		srand(conf.currentSeed);
		Surface* cs_temp = surfaceFactory().create(cf,dg_temp);			
		// relatively big such that non degenerate configurations are avoided
		// due to the high packing of atoms. this has a minimal impact
		// on the estimation of the pocket/surface volume
		cs_temp->setRandDisplacement(0.1);
		cs_temp->setTriangulationFlag(true);
		cs_temp->setCheckDuplicatedVertices(false);			
		cs_temp->build();
		cs_temp->getSurf(false);
		volumes.push_back(cs_temp->getVolume());
		if (conf.tri)
		{
			cout << endl << INFO << "Triangulating enveloped cavity/pockets " << i;				
			cs_temp->triangulateSurface(0.0,tri_);
			areas.push_back(cs_temp->getArea());			

			// NanoShaper identifies the set of grid points that are nearest 
			// to the given triangulation points and save those points that are completely
			// solvent exposed (i.e. external) in surf2 (the slim surface).
			// Those points represent a point-wise representation of the entrance of the pocket
			if (saveEntranceInfo)
			{
				if (!isPocket[i])
				{
					//cout << endl << ERR << "Cannot find the entrance of a cavity; it must be a pocket.";
					if (!conf.cavAndPockets)
					{
						cout << endl << REMARK << "You have to enable Cavities and pockets flag to do a distinction between a pocket and a cavity";
						exit(-1);
					}
				}
				else
				{
					// is a pocket, we can identify the entrance
					// true means that it is part of the entrance
					vector<bool> results;
					cs_temp->triangulationPointsAreCompletelyOut(surf2,results);

					// the points that are completely out are those that are represent the entrance
					// together with them we must save the normals
				
					// surf2 save those normals and and save those vertices selectively of this surface
					char buff[100],buff2[100];
					sprintf(buff,"entrance%d.xyz",i);
					sprintf(buff2,"entrance%d.xyzn",i);
					cs_temp->saveSelectedPoints(results,buff,buff2);

					char triSubset[100];				
					sprintf(triSubset,"cav_tri_body%d",i);								
					
					// save triangulation and estimate area removing entrance triangulation points
					vector<bool> results2;
					for (unsigned int i=0;i<results.size();i++)
						results2.push_back(!results[i]);					
					
					bool revert = false;
					double reducedArea = cs_temp->saveTriSubSet(triSubset,results2,revert);
					areas2.push_back(reducedArea);
				}
			}
		}				
		delete cs_temp;
		delete dg_temp;
	}		
		
	duration = chrono.stop();

	if (conf.debug)	restartDebug();

	if (hasAtomInfo)
	{
		char name2[100];
		sprintf(name2,"%s.info",conf.sysName.c_str());
		FILE* fp = fopen(name2,"w");
		fprintf(fp,"\n Protein             : ");
		fprintf(fp,"\n Probe radius        : %.2f",conf.pocketRadiusSmall);
		fprintf(fp,"\n Number of atoms     : %d",natom);
		fprintf(fp,"\n Total Surface Area  : %.4f",areaSurf);
		fprintf(fp,"\n Total Volume        : %.4f",volSurf);
		fprintf(fp,"\n\n");
		fprintf(fp,"\n Pockets :");
		fprintf(fp,"\n");
		fprintf(fp,"\n Id\t\tN_mth\tSurface\t\t\tVolume\n");
		
		double sumA=0,sumV=0;
		for (int i=0;i<nc;i++)
		{
			if (conf.tri)
			{
				if (conf.cavAndPockets)
				{
					if (isPocket[i]==true)
						fprintf(fp," %d\t\t>0\t\t%.4f\t\t%.4f\n",i+1,areas[i],volumes[i]);
					else
						fprintf(fp," %d\t\t0\t\t%.4f\t\t%.4f\n",i+1,areas[i],volumes[i]);
				}
				else
					fprintf(fp," %d\t\t>0\t\t%.4f\t\t%.4f\n",i+1,areas[i],volumes[i]);

				sumA+=areas[i];
				sumV+=volumes[i];
			}
			else
			{
				if (conf.cavAndPockets)
				{				
					if (isPocket[i]==true)
						fprintf(fp,"\t%d\t>0\t\t%.4f\t\t%.4f\n",i+1,0.0,volumes[i]);
					else
						fprintf(fp,"\t%d\t0\t\t%.4f\t\t%.4f\n",i+1,0.0,volumes[i]);
				}
				else
					fprintf(fp," %d\t\t>0\t\t%.4f\t\t%.4f\n",i+1,areas[i],volumes[i]);

				sumV+=volumes[i];
			}
		}
		fprintf(fp,"\nTot\t\t%.4f\t\t%.4f",sumA,sumV);
		fclose(fp);	
	}
	else
	{
		cout << endl << INFO;
		cout << endl << INFO << "------------------------------------";
		cout << endl << INFO << "     Pocket Detection Summary       ";
		cout << endl << INFO << "------------------------------------";
		cout << endl << INFO;
		cout << endl << INFO << "Detected a total of " << nc << " pockets/cavities having at least the volume of " <<conf.numMol << " water molecules" ;
			
		for (int i=0,ii=0;i<nc;i++)
		{
			if (conf.cavAndPockets)
			{
				if (conf.tri)
				{
					if (isPocket[i])
					{					
						cout << endl << INFO << "Pocket " << i << " vol " << volumes[i] << " area " << areas[i] << " body area "<< areas2[ii];
						if (internals!=NULL)
						{
							*internals << endl << "pocket_vol " << volumes[i];
							*internals << endl << "pocket_area " << areas[i] << " pocket_body_area " << areas2[i];
						}
						ii++;
					}
					else
					{
						cout << endl << INFO << "Cavity " << i << " vol " << volumes[i] << " area " << areas[i];
						if (internals!=NULL)
						{
							*internals << endl << "cav_vol " << volumes[i];
							*internals << endl << "cav_area " << areas[i];
						}
					}
				}
				else
				{
					if (isPocket[i])
					{
						cout << endl << INFO << "Pocket " << i << " vol " << volumes[i];
						if (internals!=NULL)					
							*internals << endl << "pocket_vol " << volumes[i];						
					}
					else
					{
						cout << endl << INFO << "Cavity " << i << " vol " << volumes[i];
						if (internals!=NULL)
							*internals << endl << "cav_vol " << volumes[i];												
					}
				}
			}
			else
			{
				if (conf.tri)
				{
					cout << endl << INFO << "Pocket " << i << " vol " << volumes[i] << " area " << areas[i] << " body area " << areas2[ii];
					if (internals!=NULL)
					{
						*internals << endl << "pocket_vol " << volumes[i];
						*internals << endl << "pocket_area " << areas[i];
					}
					ii++;
				}
				else
				{
					cout << endl << INFO << "Pocket " << i << " vol " << volumes[i];
					if (internals!=NULL)
						*internals << endl << "pocket_vol " << volumes[i];
				}
			}
		}
	}
								
	cout << endl << INFO << "Pocket detection time.. " << duration << " [s]"; 	

	cout << endl << INFO << "Cleaning memory...";
	cout.flush();

	delete surf1;
	delete dg1;		
	delete surf2;
	delete dg2;
	cout << "ok!";

}
