
//---------------------------------------------------------
/**    @file	BlobbySurface.h
*     @brief	BlobbySurface.h is the header for CLASS
                BlobbySurface.cpp
*/
//---------------------------------------------------------

#ifndef BlobbySurface_h
#define BlobbySurface_h

#include <Configuration.h>
#include <MeshSurface.h>
#include <Surface.h>
#include <SurfaceFactory.h>
#include <globals.h>
#include <memory>

#ifdef DBGMEM_CRT
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW
#endif

#define DEFAULT_BLOBBYNESS -2.5
#define DEFAULT_CUTOFF 6
/** @brief This class builds the blobby surface, triangulate it and the use mesh
surface routines for the rest of the work.

@author Sergio Decherchi
@date 10/01/2012
*/
class BlobbySurface : public MeshSurface {
private:
  /** blobbyness*/
  double B;
  /** cut-off distance in Angstrom to speed-up computations. It is computed
  based on B value. The bigger the gaussian the less the cut-off*/
  double cutoff;

public:
  /** Default constructor*/
  BlobbySurface();
  /** set DelPhi environment*/
  BlobbySurface(DelPhiShared *ds);
  BlobbySurface(ConfigurationOP cf, DelPhiShared *ds);

  //////////////////////// INTERFACE MANDATORY METHODS
  ////////////////////////////////////
  /** Compute blobby mesh and run MeshSurface::build()*/
  virtual bool build();
  /** Print number a summary of the blobby*/
  virtual void printSummary();
  /** function for the constructor without arguments*/
  virtual void init();
  /** functions for the constructor with config file argument*/
  virtual void init(ConfigurationOP cf);
  /**function for the denstructor*/
  virtual void clear();
  /////////////////////////////////////////////////////////////

  //////////////////////////// BLOBBY SURFACE METHODS
  //////////////////////////////////
  void setBlobbyness(double b);
  double getBlobbyness();

  virtual ~BlobbySurface();
};

// expand it explicitly because Swig is not able to expand it
static class BlobbySurfaceRegister {
  static SurfaceOP createSurface(ConfigurationOP conf, DelPhiShared *ds) {
    return std::make_shared<BlobbySurface>(conf, ds);
  }

public:
  BlobbySurfaceRegister() {
    SurfaceFactory::getInstance().register_instantiator("blobby",
                                                        createSurface);
  }
} BlobbySurfaceRegisterObject;

// static SurfaceRecorder<BlobbySurface> blobbyRecorder("blobby");

#endif
