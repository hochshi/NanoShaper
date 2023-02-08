
#ifndef SurfaceFactory_h
#define SurfaceFactory_h

#include <map>
#include <memory>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string>
#include <iostream>
#include <ConfigFile.h>
#include <globals.h>

class Surface;
using SurfaceOP = std::shared_ptr<Surface>;
class DelPhiShared;
using DelPhiSharedOP = std::shared_ptr<DelPhiShared>;
using namespace std;

typedef SurfaceOP (*surface_instantiator)(ConfigFile* conf,DelPhiShared* ds);

class SurfaceFactory{
	
private:
	map<string,surface_instantiator> surfRegister;
	map<string,surface_instantiator>::iterator it;
public:
	
	void register_instantiator(string surfaceName,surface_instantiator si)
	{
		surfRegister.insert(pair<string,surface_instantiator>(surfaceName,si));
	}
	
	SurfaceOP create(ConfigFileOP conf,DelPhiSharedOP ds)
	{		
		string surfName = conf->read<string>("Surface");
		if (!surfRegister.count(surfName))
		{
      spdlog::error("{} type is not registered!", surfName);
			return NULL;
		}		
		return surfRegister[surfName](conf.get(),ds.get());
	}

	void print()
	{
    std::stringstream ss;
		ss << endl << INFO_STR << "Available surfaces:";
		for (it=surfRegister.begin();it!=surfRegister.end();it++)
		{
			ss << endl << INFO_STR << "\t" << (*it).first;
		}
    spdlog::info("{}", ss.str());
	}
};

SurfaceFactory& surfaceFactory();

template<class T> class SurfaceRecorder 
{

public:	
    SurfaceRecorder(string surface)
    {
		surfaceFactory().register_instantiator(surface,createSurface);
    }

	static SurfaceOP createSurface(ConfigFile* conf,DelPhiShared* ds) 
	{ 
    return std::make_shared<T>(conf, ds);
	} 
};

#endif
