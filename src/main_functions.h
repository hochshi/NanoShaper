//---------------------------------------------------------
/*    @file		main_functions.h
 *     @brief	main_functions.h Includes functions used in main.cpp
 *							    							*/
//---------------------------------------------------------

#ifndef main_functions_h
#define main_functions_h

#include <ConfigFile.h>
#include <DelphiShared.h>
#include <Surface.h>
#include <globals.h>
#include <memory>

ConfigFileOP load(std::string argv);
ConfigurationOP parse(ConfigFileOP cf);
// void dispose(ConfigFileOP cf);
// void stopDebug();
// void restartDebug();
void cite();
void normalMode(SurfaceOP surf, DelPhiSharedOP dg, ConfigurationOP conf);
void pocketMode(bool hasAtomInfo, ConfigFileOP cf, ConfigurationOP conf);

#endif
