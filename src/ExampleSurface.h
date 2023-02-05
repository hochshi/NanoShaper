//---------------------------------------------------------
/**    @file	ExampleSurface.h
*     @brief	ExampleSurface.h is the header for CLASS
*               ExampleSurface.cpp								*/
//---------------------------------------------------------

#ifndef ExampleSurface_h
#define ExampleSurface_h

#include "Surface.h"

// 1 Angstrom default radius
#define DEFAULT_RADIUS 1.0

/** @brief This class is a simple example that shows the main concepts needed to
introduce a new surface in NanoShaper. This class builds a sphere centered on the baricenter
of input points whose radius is specified by the user. 
This surface is ray-casting based and it is a trivial Surface definition. On a real scenario
high performance is needed; to get performance the Developer can use the same strategies used in 
other surfaces such as in classes ConnollySurface and SkinSurface. 
The acceleration strategies are 2D/3D grids that allow fast retrieval of points or in general
of simplices and patches.
To add a new surface it is sufficient that you define a class withi this style (see registration class)
and then it will be automatically loaded; you don't have to change any make file or cmake file, just
write the source and put it the src folder.

@author Sergio Decherchi 
@date 29/10/2013
*/
class ExampleSurface: public Surface
{

private:

	double radius;
	double center[3];

public:
	/** Default constructor*/
	ExampleSurface();
	/** set DelPhi environment*/
	ExampleSurface(DelPhiShared* ds);			
	/** set configuration and DelPhi environment*/
	ExampleSurface(ConfigFile* cf,DelPhiShared* ds);			

	//////////////////////// INTERFACE MANDATORY METHODS /////////////////////////////////
	/** Compute the surface. Call it after load*/
	virtual bool build();
	/** Save it in a simple ASCII format (.example). It is not implemented*/
	virtual bool save(char* fileName);
	/** Load the surface from a file in .example format. It is not implemented*/
	virtual bool load(char* fileName);
	/** Print a summary of the surface*/
	virtual void printSummary();		
	/** Get a projection of a point on the surface. Return projection and normal*/
	virtual bool getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3);
	/** Get all the intersections of a ray that goes from P1 to P2 over the surface.
	The interesctions are returned with increasing distance order. 
	the first double in the vector is the t parameter for the intersection of the parametric 
	line and the surface, the double pointer is the normal vector.*/
	virtual void getRayIntersection(double p1[3],double p2[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals);
	/** function for the constructor without arguments*/
	virtual void init();
	/** functions for the constructor with config file argument*/
	virtual void init(ConfigFile* cf);
	/**function for the denstructor*/
	virtual void clear();
	/** here you should put your code to prepare a 2D/3D acceleration structure for ray tracing*/
	virtual void preProcessPanel();
	/** here you should put your code clean memory after ray casting*/
	virtual void postRayCasting();
	/** here you should put your code to prepare a 2D/3D acceleration structure for point to surface projection*/
	virtual bool preBoundaryProjection();
	/////////////////////////////////////////////////////////////

	void setRadius(double r) 
	{ 
		if (r<=0)
		{
			cout << endl << WARN << "Cannot set a radius <0. Setting the default value";
			radius = DEFAULT_RADIUS;
		}
		else
			radius = r;
	}
	
	double getRadius() 
	{
		return radius;
	}
	
	virtual ~ExampleSurface();

};


/** This is class used to register the surface, that is to automatically
notify its existence to NanoShaper. For every surface that you define
you have to follow the naming convention that follows in the class. The following variable
declaration, (please look at it) works perfectly and builds automatically this class but 
Swig seems to have problem with it so we write it explicitly*/

/*
static class ExampleSurfaceRegister{ 
	static Surface* createSurface(ConfigFile* conf,DelPhiShared* ds) 
	{ 
		return new ExampleSurface(conf,ds); 
	} 
	public: 
		ExampleSurfaceRegister() 
		{ 
			surfaceFactory().add("example",createSurface); 
		} 
} ExampleSurfaceRegisterObject;
*/

static SurfaceRecorder<ExampleSurface> exampleRecorder("example");

#endif