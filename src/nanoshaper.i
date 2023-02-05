/* NanoShaper Python-C++ Swig Interface */
%module NanoShaper
%{
#include <string>
#include "globals.h"
#include "MeshSurface.h"
#include "ExternalSurface.h"
#include "SkinSurface.h"
#include "tnt/tnt_array1d.h"
#include "tnt/tnt_array2d.h"
#include "tnt/tnt_math_utils.h"
#include "jama/jama_eig.h"
#include "BlobbySurface.h"
#include "ConnollySurface.h"
#include "Surface.h"
#include "tools.h"
#include "ConfigFile.h"
#include "DelphiShared.h"
#include "octree.h"
#include "SurfaceFactory.h"
ConfigFile* init(string a);
void dispose(ConfigFile* cf);
void stopDebug();
void restartDebug();
void cite();
void normalMode(Surface* surf,DelPhiShared* dg);
void pocketMode(bool hasAtomInfo,ConfigFile* cf);
%}
int num_cores;
%include "std_string.i"
%include "DelphiShared.h"
%include "Surface.h"
%include "globals.h"
%include "MeshSurface.h"
%include "tnt/tnt_array1d.h"
%include "tnt/tnt_array2d.h"
%include "tnt/tnt_math_utils.h"
%include "jama/jama_eig.h"
%include "ExternalSurface.h"
%include "SkinSurface.h"
%include "BlobbySurface.h"
%include "ConnollySurface.h"
%include "tools.h"
%include "ConfigFile.h"
%include "SurfaceFactory.h"
ConfigFile* init(string a);
void dispose(ConfigFile* cf);
void stopDebug();
void restartDebug();
void cite();
void normalMode(Surface* surf,DelPhiShared* dg);
void pocketMode(bool hasAtomInfo,ConfigFile* cf);