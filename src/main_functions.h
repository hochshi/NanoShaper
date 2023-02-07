//---------------------------------------------------------
/*    @file		main_functions.h
*     @brief	main_functions.h Includes functions used in main.cpp
*							    							*/
//---------------------------------------------------------

#ifndef main_functions_h
#define main_functions_h

#include <globals.h>
#include <DelphiShared.h>
#include <Surface.h>

ConfigFile* init(std::string argv);
void dispose(ConfigFile* cf);
void stopDebug();
void restartDebug();
void cite();
void normalMode(Surface* surf,DelPhiShared* dg);
void pocketMode(bool hasAtomInfo,ConfigFile* cf);

#endif
