
#ifndef DelphiShared_h
#define DelphiShared_h

#include <globals.h>
#include <memory>
#include <vector>

#ifdef DBGMEM_CRT
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW
#endif

#include <tools.h>

#define DEFAULT_PERFIL 80
#define DEFAULT_SCALE 2.0

// Status map codes
#define STATUS_CAVITY_POINT_SHAPE_UNDER_CHECK -3
#define STATUS_CAVITY_POINT_SHAPE_OK -2
#define STATUS_POINT_INSIDE -1
#define STATUS_POINT_TEMPORARY_OUT 2
// over this index are all cavities
#define STATUS_POINT_OUT 3
// normal cavity point
#define STATUS_FIRST_CAV 4
// support points are those which have in absolute number a cavity number (>=4)
// but an inverted sign. These points represent grid point in which the probe
// can stay. In a generic cavity point probe could not fit.
#define STATUS_FIRST_SUPPORT_CAV -4

// extern double debug_stern;

/** @brief This class is the DelPhi Shared variables environment which emulates
 * the DelPhi environment.
 */
class DelPhiShared {
public:
  ///////////////////////////////////// methods
  /////////////////////////////////////////////////////////
  /** default constructor*/
  DelPhiShared();

  /** default init*/
  void init();

  /** constructor. map indicates epsmap and idebmap, status is status map and
   * multi mean multidielectric*/
  DelPhiShared(double scale, double perfill, string fn, bool map, bool status,
               bool multi, bool atinfo = false);

  /** init function-constructor*/
  void init(double scale, double perfill, string fn, bool eps_flag,
            bool stat_flag, bool multi, bool atinfo);

  /** Fortran DelPhi binding function. Allows to link DelPhi old code to
      the new surface framework */
  void DelPhiBinding(double xmin, double ymin, double zmin, double xmax,
                     double ymax, double zmax, double c1, double c2, double c3,
                     double rmax, double perf, int *local_i_epsmap, int igrid,
                     double scale, double *local_i_scspos,
                     double *local_i_scsnor, bool *local_i_idebmap,
                     int *local_i_ibgp, int maxbgp, bool status, int *atsurf);

  /**If in DelPhi binding mode, call it before denstructor. It will be change
  the variables in a Fortran/DelPhi compatible mode and will avoid the
  denstruction of the variables passed to DelPhi. Always call this function ONLY
  just before the denstructor of the DelPhi shared environment*/
  void finalizeBinding(int *ibnum);
  /** deletion function*/
  void clear();
  /** Save the epsmaps in a format such that DelPhi can read them*/
  void saveEpsMaps(char *fname);
  /** Save the boundary grid points in a format such that DelPhi can read them*/
  void saveBGP(char *fname);
  /** Save internal status map for cavity detection*/
  void saveStatus(char *fname);
  /** Save only cavities info*/
  void saveCavities(bool onlySaveNonFilled = false);
  /** mimicks ProShape output*/
  void saveCavities2(bool onlySaveNonFilled, string sysName);
  /** Convert grid cavity points in atoms files for NanoShaper.
  This function is usefull to triangulate by NanoShaper the envelope of the
  cavities*/
  int cavitiesToAtoms(double rad);
  /** given the status map of a 1.4 surf, mark for all the active cavites if
   * they are pockets or cavities*/
  void markPockets(short *status, vector<bool> &isPocket);
  /** clear reference atoms for each cavity*/
  void clearCav2Atoms();
  /** reinit cav2atoms vector to current number of cavities. If not empty it
   * gets cleared*/
  void initCav2Atoms();
  /** parse atinfo vector from DelPhi*/
  static void parseAtomInfo(char *atinfo, string mol, int na, double *xn1,
                            double *rad);
  /** Save idebmap*/
  void saveIdebMap(char *fname);
  bool clearAndAllocEpsMaps();
  void clearEpsMaps();
  int loadAtoms(string fn);
  int loadAtoms(int na, double *pos, double *rad, double *charge, int *dielec, char *atinf);
  int loadAtoms(int na, const std::vector<double> &x,
                 const std::vector<double> &y, const std::vector<double> &z,
                 const std::vector<double> &r, const std::vector<double> &q,
                 const std::vector<int> &d, const std::vector<AtomInfo> &ai);
  int loadAtoms(int na, const double *x, const double *y, const double *z,
                 const double *r, const double *q, const int *d,
                 const AtomInfo *ai);
  /** Emulates DelPhi grid construction*/
  bool buildGrid(double scale, double perfill);
  bool getDelphiBinding() { return delphiBinding; }

  void buildEpsmap(bool a) { buildEpsMap = a; }

  void buildStatusmap(bool a) { buildStatus = a; }

  bool getMultiDiel() { return multi_diel; }

  int getNumAtoms() { return numAtoms; }
  virtual ~DelPhiShared();

  ///////////////////////////////////////////////////// variables
  /////////////////////////////////////////

  // multi dielectric mode
  bool multi_diel;
  double scale;
  double perfill;
  bool delphiBinding;
  bool buildEpsMap;
  bool buildStatus;
  // epsmap is a four dimensional matrix
  int *epsmap, *ibgp;
  bool *idebmap;
  int maxbgp;
  double *scsnor, *scsarea, *scspos;
  bool isAvailableAtomInfo;
  int *atsurf;
  // face area = side*side;
  double A;
  double *x;
  double *y;
  double *z;
  int nbgp;
  double side;
  int nx, ny, nz;
  // the status matrix is the status for each grid point:
  // -5 -> cavity point is inside a well shaped cavity (temporary)
  // -10 -> cavity point under shape check (temporary)
  // -1 -> point is inside
  // 2 -> point is temporary outside (may be really out or in cavity)
  // 3 -> point is really out
  // 4 -> point is in cavity number 1
  // 5 -> point is in cavity number 2
  // 6 -> point is in cavity number 3
  // and so on....
  // Thus 65536-4 cavities can be detected at best; it should be enough...
  short *status;
  // This temporary variable is introduce to detect if a cavity has been split
  // in two subcavities during Connolly filtering. Is used by the difference
  // function
  short *tempStatus;

  double rmaxdim;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double baricenter[3];
  // atoms vector. Pointer to atom objects
  Atom **atoms;
  // num atoms
  int numAtoms;
  char file[BUFLEN];
  double hside;
  /** dynamical vector of cavities*/
  vector<vector<int *> *> *cavitiesVec;
  /** dynamical vector of cavities size*/
  vector<double> cavitiesSize;
  /** dynamical vector of cavities filling. true if cav is filled*/
  vector<bool> cavitiesFlag;
  /** dynamical vector whose index is the cavity and contains the set of atoms
   * that produce that cavity*/
  vector<set<int> *> cav2atoms;
};

using DelPhiSharedOP = std::shared_ptr<DelPhiShared>;

#endif
