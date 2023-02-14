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

namespace nanoshaper {
ConfigFileOP load(std::string argv, string delimiter = "=",
                  string comment = "#", string sentry = "EndConfigFile",
                  std::string format = "plain");
ConfigurationOP parse(ConfigFileOP cf);
// void dispose(ConfigFileOP cf);
// void stopDebug();
// void restartDebug();
void cite();
void normalMode(SurfaceOP surf, DelPhiSharedOP dg, ConfigurationOP conf);
void pocketMode(bool hasAtomInfo, ConfigurationOP cf, ConfigurationOP conf);

}  // namespace nanoshaper
#endif
