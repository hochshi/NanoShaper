#include "ExampleSurface.h"

// denstructor should always call a clear function
// this makes the code clearer and it allows to reset an object
// without really destroying it
ExampleSurface::~ExampleSurface()
{
	clear();
}

void ExampleSurface::clear()
{
}

// constructors should always use a helper function. 
// such that after a clear a reinit can be done without destroying the object
void ExampleSurface::init()
{
	radius = DEFAULT_RADIUS;
	// every constructor should notify the kind of surface
	surfType = MOLECULAR_SURFACE;
	// every constructor should notify to the mother class if it is
	// able to provide the analytical normals
	// if the class declares it is not able to give analytical normals
	// then the normal can still be approximated from the triangulation
	providesAnalyticalNormals = true;
}

ExampleSurface::ExampleSurface(DelPhiShared* ds):Surface()
{
	init();
	// set the environment, i.e. the grid
	delphi = ds;	
}

ExampleSurface::ExampleSurface():Surface()
{
	init();
}

// here you can parse all the custom key words that you add in the configuration file
void ExampleSurface::init(ConfigFile* cf)
{
	// we have the configuration file; so we can query the radius parameter (or our custom parameters)
	// the second argument is the default value if the radius keyword is not found
	radius = cf->read<double>( "Example_Surface_Parameter", DEFAULT_RADIUS);	
}

ExampleSurface::ExampleSurface(ConfigFile* cf,DelPhiShared* ds):Surface(cf)
{
	init();
	init(cf);
	// set environment
	delphi = ds;	
}


inline void ExampleSurface::getRayIntersection(double pa[3],double pb[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals)
{
	// the ray starts from pa, arrives at pb.
	// if computeNormals is active you should provide normals
	// thdID is the thread identifier. We will not use it.
	// intersections contains as first value the parameter t of the ray pa+t*(pb-pa) and the second is the normal vector
	double dir[3];	
	dir[0]=pb[0]-pa[0];
	dir[1]=pb[1]-pa[1];
	dir[2]=pb[2]-pa[2];
	
	// we know the maximal number of intersections in this case; if we don't know that we can make a guess
	intersections.reserve(2);

	// in real scenario you should preliminarly identify a superset of patches that could
	// intersect the ray and then test only them. In this trival case we have only one patch, the sphere

	// these are two t parameters of the ray. For one sphere we can have up to 2 intersections
	double t1,t2;

	bool hasIntersection = raySphere(pa,dir,center,radius,&t1,&t2);
	
	// we will compute the intersection point explicitly if the normal is required
	double intPoint1[3],intPoint2[3];

	// the variable 'panel' is controlled by the mother class. You should not touch it
	// it identifies the panel from which we are casting rays
	
	if (panel==PANEL_YZ)
	{							
		if (hasIntersection && computeNormals)
		{
			intPoint1[0]=pa[0]+dir[0]*t1;
			intPoint1[1]=pa[1];
			intPoint1[2]=pa[2];

			intPoint2[0]=pa[0]+dir[0]*t2;
			intPoint2[1]=pa[1];
			intPoint2[2]=pa[2];
		}
	}			
	else if(panel==PANEL_XY)
	{
		if (hasIntersection && computeNormals)
		{			
			intPoint1[0]=pa[0];
			intPoint1[1]=pa[1];
			intPoint1[2]=pa[2]+dir[2]*t1;

			intPoint2[0]=pa[0];
			intPoint2[1]=pa[1];
			intPoint2[2]=pa[2]+dir[2]*t2;
		}
	}	
	else if (panel==PANEL_XZ)
	{
		if (hasIntersection && computeNormals)
		{			
			intPoint1[0]=pa[0];
			intPoint1[1]=pa[1]+dir[1]*t1;
			intPoint1[2]=pa[2];

			intPoint2[0]=pa[0];
			intPoint2[1]=pa[1]+dir[1]*t2;
			intPoint2[2]=pa[2];
		}
	}
	else
	{
		cout << endl << ERR << "Panel does not exist!";
		exit(-1);
	}

	if (hasIntersection)
	{
		if (computeNormals)
		{
			double* normal1 = allocateVector<double>(3);
			double* normal2 = allocateVector<double>(3);

			getNormalToSphere(intPoint1,center,radius,normal1);
			getNormalToSphere(intPoint2,center,radius,normal2);
					
			intersections.push_back(pair<double,double*>(t1,normal1));
			intersections.push_back(pair<double,double*>(t2,normal2));
		}
		else
		{			
			intersections.push_back(pair<double,double*>(t1,(double*)NULL));
			intersections.push_back(pair<double,double*>(t2,(double*)NULL));
		}

		// we must sort the intersections such that Surface class can deduce in/out info
		// and do checksum
		sort(intersections.begin(), intersections.end(),compKeepIndex);	
	}
}					

bool ExampleSurface::getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3)
{
	// You have to implement this function if you plan to interface this surface to DelPhi.
	// Only in this case it is needed that you really imnplement this method
	// As long as triangulation/cavity detectionm, area and volume are need only the intersection
	// routine is needed.

	// in a real scenario you should identify the nearest patch and then project to it;
	// you don't have to wonder if p must be projected or not. p has been already found
	// as a boundary grid point, that is a point in the nearby of the surface, so for sure
	// we have to project it. The only caveat is that for numerical reasons you could not be able to project
	// it. In this last case the safest thing to do is to set as the projection the (grid) point itself;
	// in this case the normal cannot be provided
	
	double proj[3],norm[3],dist;	
	projectToSphere(p,center,radius,proj,dist);
	getNormalToSphere(proj,center,radius,norm);

	(*proj1)=proj[0];
	(*proj2)=proj[1];
	(*proj3)=proj[2];

	(*normal1)=norm[0];
	(*normal2)=norm[1];
	(*normal3)=norm[2];

	return true;
}

void ExampleSurface::printSummary()
{
	cout << endl << INFO << "Example surface has a radius of " << radius;
}

void ExampleSurface::preProcessPanel()
{
	// if we used an acceleration data structure we could build it here
	// note that this function is called for each of the three coordinated panels
	// so it is called three times. 
	// This is done because, as for the current Connolly, Skin, and Mesh surfaces
	// the acceleration data structure dependes on the panel from which we are tracing
	// rays. 
}
void ExampleSurface::postRayCasting()
{
	// if we used an acceleration data structure we could clean it here
}

bool ExampleSurface::preBoundaryProjection()
{
	// if we used an acceleration data structure we could build it here
	return true;
}

bool ExampleSurface::build()
{
	// here we build the analytical description of the surface
	// in this case the only parameter to compute is the center.
	// once we find the center, we just check that the user radius is smaller
	// than the maximal distance from the baricenter to any atom just to assure
	// we don't go out of the embedding grid. If the radius is too big
	// we set it to the this maximal distance

	center[0]=0;
	center[1]=0;
	center[2]=0;

	int na = delphi->numAtoms;

	for (int i=0;i<na;i++)
	{
		double x = delphi->atoms[i]->pos[0];
		double y = delphi->atoms[i]->pos[1];
		double z = delphi->atoms[i]->pos[2];

		center[0]+=x;
		center[1]+=y;
		center[2]+=z;
	}

	center[0]/=na;
	center[1]/=na;
	center[2]/=na;

	double maxDist = 0, squaredDist;

	for (int i=0;i<na;i++)
	{
		DIST2(squaredDist,(delphi->atoms[i]->pos),center)
		if (squaredDist>maxDist)
			maxDist = squaredDist;
	}
	maxDist = sqrt(maxDist);

	if (maxDist<radius)
	{
		cout << endl << WARN << "Max distance is lower than radius. Setting radius = max distance";
		radius = maxDist;
	}
	return true;
}

bool ExampleSurface::save(char* fileName)
{
	// you should save the surface in your prefered format
	return true;
}


bool ExampleSurface::load(char* fileName) 
{
	// you should load the surface in your prefered format
	return true;
}