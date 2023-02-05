
//---------------------------------------------------------
/**    @file		ExternalSurface.h
*     @brief		ExternalSurface.h is the header for CLASS
*					ExternalSurface.cpp								*/
//---------------------------------------------------------

#ifndef ExternalSurface_h
#define ExternalSurface_h

#include "Surface.h"


/** @brief This class is wrapper and loads an external surface using the files epsmapx.txt, epsmapy.txt,
epsmapz.txt, status.txt and projections.txt that is the list of the boundary grid points indexes,values and normals

@author Sergio Decherchi 
@date 29/10/2011
*/
class ExternalSurface: public Surface
{
public:
	/** Default constructor*/
	ExternalSurface();
	/** set DelPhi environment*/
	ExternalSurface(DelPhiShared* ds);		
	/** set configuration and DelPhi environment*/
	ExternalSurface(ConfigFile* cf,DelPhiShared* ds);	

	//////////////////////// INTERFACE MANDATORY METHODS /////////////////////////////////
	/** Do nothing*/
	virtual bool build();
	/** Save the epsmaps and projections in a simple ascii format*/
	virtual bool save(char* fileName);
	/**Do nothing*/
	virtual bool load(char* fileName);
	/**Load the external surface and perform cavity detection if requested*/
	virtual bool getSurf(bool fill=false,double vol=0);		
	/** Print number of boundary grid points*/
	virtual void printSummary();		
	/** Do nothing*/
	virtual bool getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3);
	/** Do nothing*/
	virtual void getRayIntersection(double p1[3],double p2[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals);
	/** Do nothing*/
	virtual void getLocalArea(double gridPoint[3],double* area);
	/** Do nothing*/
	virtual double getSurfaceArea();
	/** Approximate the volume by summing the inner boxes of the grid plus the partial volumes
	of the boundary boxes.*/
	virtual double getVolume();
	/** function for the constructor without arguments*/
	virtual void init();
	/** functions for the constructor with config file argument*/
	virtual void init(ConfigFile* cf);
	/**function for the denstructor*/
	virtual void clear();
	/////////////////////////////////////////////////////////////


	virtual ~ExternalSurface();
};


//	REGISTER_SURFACE(ExternalSurface,"external")

// expand it explicitly because Swig is not able to expand it
static class ExternalSurfaceRegister{ 
	static Surface* createSurface(ConfigFile* conf,DelPhiShared* ds) 
	{ 
		return new ExternalSurface(conf,ds); 
	} 
	public: 
		ExternalSurfaceRegister() 
		{ 
			surfaceFactory().add("external",createSurface); 
		} 
} ExternalSurfaceRegisterObject;

#endif