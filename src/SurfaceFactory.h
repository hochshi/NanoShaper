
#ifndef SurfaceFactory_h
#define SurfaceFactory_h

#include <Configuration.h>
#include <globals.h>
#include <logging.h>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

class Surface;
using SurfaceOP = std::shared_ptr<Surface>;
class DelPhiShared;
using DelPhiSharedOP = std::shared_ptr<DelPhiShared>;
using namespace std;

typedef SurfaceOP (*surface_instantiator)(ConfigurationOP conf,
                                          DelPhiSharedOP ds);

class SurfaceFactory {

 private:
  SurfaceFactory() = default;
  ~SurfaceFactory() = default;

  map<string, surface_instantiator> surfRegister;
  map<string, surface_instantiator>::iterator it;

 public:
  SurfaceFactory(const SurfaceFactory&) = delete;
  SurfaceFactory& operator=(const SurfaceFactory&) = delete;

  static SurfaceFactory& getInstance() {
    static SurfaceFactory instance;
    return instance;
  }

  void register_instantiator(string surfaceName, surface_instantiator si) {
    surfRegister.insert(pair<string, surface_instantiator>(surfaceName, si));
  }

  SurfaceOP create(ConfigurationOP conf, DelPhiSharedOP ds) {
    string surfName = conf->surfName;
    if (!surfRegister.count(surfName)) {
      logging::log<logging::level::err>("{} type is not registered!", surfName);
      return nullptr;
    }
    return surfRegister[surfName](conf, ds);
  }

  void print() {
    std::stringstream ss;
    ss << endl << INFO_STR << "Available surfaces:";
    for (it = surfRegister.begin(); it != surfRegister.end(); it++) {
      ss << endl << INFO_STR << "\t" << (*it).first;
    }
    logging::log<logging::level::info>("{}", ss.str());
  }
};

// SurfaceFactory& surfaceFactory();

template <class T>
class SurfaceRecorder {

 public:
  SurfaceRecorder(string surface) {
    SurfaceFactory::getInstance().register_instantiator(surface, createSurface);
  }

  static SurfaceOP createSurface(ConfigurationOP conf, DelPhiSharedOP ds) {
    return std::make_shared<T>(conf, ds);
  }
};

#endif
