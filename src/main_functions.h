//---------------------------------------------------------
/*    @file		main_functions.h
*     @brief	main_functions.h Includes functions used in main.cpp
*							    							*/
//---------------------------------------------------------

#ifndef main_functions_h
#define main_functions_h

#include <DelphiShared.h>
#include <Surface.h>
#include <globals.h>

ConfigFile* init(string argv);
void dispose(ConfigFile* cf);
void stopDebug();
void restartDebug();
void cite();
void normalMode(Surface* surf,DelPhiShared* dg);
void pocketMode(bool hasAtomInfo,ConfigFile* cf);

#endif
