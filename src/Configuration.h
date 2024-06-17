#ifndef Configuration_h
#define Configuration_h

#ifdef JSON_ENABLED
#include <nlohmann/json.hpp>
#endif  // JSON_ENABLED

#include <string>

namespace nanoshaper {
struct Configuration {
 public:
  void stopDebug() { debugStatus = false; }

  void restartDebug() { debugStatus = true; }

  double cavVol;
  int numMol;

  // grid (DelPhi) params
  double scale;
  double perfill;
  bool isAvailableAtomInfo = false;
  // mol file name
  std::string molFile;
  // sys name
  std::string sysName;

  bool multi_diel;

  // actions
  bool fillCavities;
  bool buildEpsmaps;
  bool buildStatus;
  bool tri;
  bool accTri;
  bool smoothing;
  bool tri2balls;
  bool projBGP;

  // save data
  bool saveEpsmaps;
  bool saveIdebmap;
  bool saveBgps;
  bool saveStatusMap;
  bool saveCavities;

  // global parameters
  std::string operativeMode;
  int numthd;
  bool printAvailSurf;
  int currentSeed;

  // pocket detection
  bool cavAndPockets;
  bool linkPockets;
  double pocketRadiusBig;
  double pocketRadiusSmall;
  double pocketRadiusLink;
  bool debug;
  bool debugStatus;

  // Rest of the values came from the code
  double blobby_B = -2.5;
  unsigned int maxSESDim2D = 50;
  unsigned int maxSESPatches2D = 400;
  unsigned int maxSESDim = 100;
  unsigned int maxSESPatches = 400;
  int mp = 100;
  double si_perfil = 1.5;
  double radius = 1.0;

  std::string sfname = "mesh.off";
  unsigned int maxMeshDim = 100;
  unsigned int maxMeshPatches = 250;
  unsigned int maxMeshDim2D = 100;
  unsigned int maxMeshPatches2D = 250;
  int NumMSMSfiles = 1;

  double skin_s = 0.45;
  unsigned int maxSkinDim = 100;
  unsigned int maxSkinPatches = 400;
  unsigned int maxSkinDim2D = 50;
  unsigned int maxSkinPatches2D = 400;
  bool useFastProjection = false;
  bool savePovRay = false;

  bool checkDuplicatedVertices = true;
  bool wellShaped = false;
  double probeRadius = 1.4;
  bool lb = true;
  bool vaFlag = false;
  bool computeNormals = false;
  bool saveMSMS = false;
  bool savePLY = false;
  double sternLayer = -1.;
  int Max_Atoms_Multi_Grid = 100;
  std::string surfName;
  // string surfName = "ses";
};

#ifdef JSON_ENABLED
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(
    Configuration, cavVol, numMol, scale, perfill, molFile, sysName, multi_diel,
    fillCavities, buildEpsmaps, buildStatus, tri, accTri, smoothing, tri2balls,
    projBGP, saveEpsmaps, saveIdebmap, saveBgps, saveStatusMap, saveCavities,
    operativeMode, numthd, printAvailSurf, currentSeed, cavAndPockets,
    linkPockets, pocketRadiusBig, pocketRadiusSmall, pocketRadiusLink, debug,
    debugStatus, blobby_B, maxSESDim2D, maxSESPatches2D, maxSESDim,
    maxSESPatches, mp, si_perfil, radius, sfname, maxMeshDim, maxMeshPatches,
    maxMeshDim2D, maxMeshPatches2D, NumMSMSfiles, skin_s, maxSkinDim,
    maxSkinPatches, maxSkinDim2D, maxSkinPatches2D, useFastProjection,
    savePovRay, checkDuplicatedVertices, wellShaped, probeRadius, lb, vaFlag,
    computeNormals, saveMSMS, savePLY, sternLayer, Max_Atoms_Multi_Grid, surfName)
#endif

using ConfigurationOP = std::shared_ptr<Configuration>;
}  // namespace nanoshaper
#endif
