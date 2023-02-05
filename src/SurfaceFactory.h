
#ifndef SurfaceFactory_h
#define SurfaceFactory_h

#include <map>
#include <string>
#include <iostream>
#include "ConfigFile.h"
#include "globals.h"

class Surface;
class DelPhiShared;
using namespace std;

typedef Surface* (*surface_instantiator)(ConfigFile* conf,DelPhiShared* ds);

class SurfaceFactory{
	
private:
	map<string,surface_instantiator> surfRegister;
	map<string,surface_instantiator>::iterator it;
public:
	
	void add(string surfaceName,surface_instantiator si)
	{
		surfRegister.insert(pair<string,surface_instantiator>(surfaceName,si));
	}
	
	Surface* create(ConfigFile* conf,DelPhiShared* ds)
	{		
		string surfName = conf->read<string>("Surface");
		if (!surfRegister.count(surfName))
		{
			cout << endl << surfName << " type is not registered!";
			return NULL;
		}		
		return surfRegister[surfName](conf,ds);
	}

	void print()
	{
		cout << endl << INFO << "Available surfaces:";
		for (it=surfRegister.begin();it!=surfRegister.end();it++)
		{
			cout << endl << INFO << "\t" << (*it).first;
		}
	}
};

SurfaceFactory& surfaceFactory();

template<class T> class SurfaceRecorder 
{

public:	
    SurfaceRecorder(string surface)
    {
		surfaceFactory().add(surface,createSurface);
    }

	static Surface* createSurface(ConfigFile* conf,DelPhiShared* ds) 
	{ 
		return new T(conf,ds); 
	} 
};

#endif