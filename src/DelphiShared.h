
#ifndef  DelphiShared_h
#define  DelphiShared_h

#include "globals.h"

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif

#include "tools.h"

// it now contains the Atom class that serves as interface
#include "nanoshaper.h"
using namespace NS;

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
// support points are those which have in absolute number a cavity number (>=4) but
// an inverted sign. These points represent grid point in which the probe can stay.
// In a generic cavity point probe could not fit.
#define STATUS_FIRST_SUPPORT_CAV -4

//extern double debug_stern;

/** @brief This class is the DelPhi Shared variables environment which emulates the DelPhi environment.
*/
class DelPhiShared {   
public:
	///////////////////////////////////// methods //////////////////////////////////////////////////////
	/** default constructor*/
	DelPhiShared();

	/** default init*/
	void init();

	/** constructor. map indicates epsmap and idebmap, status is status map and multi means multidielectric. */
	DelPhiShared(int max_num_atoms,double domain_shrinkage,bool optimize_grids,double scale,double perfill,string fn,bool map,bool status,bool multi,bool atinfo=false);

	/** constructor for API usage. Gives the atoms list directly instead that passing a file name*/
	DelPhiShared(double scale,double perfill,const std::vector<Atom> &atoms,int max_num_atoms,double domain_shrinkage,bool optimize_grids,bool map=false,bool status=true,bool multi=false);
	
	// SD: PB_NEW
	/** constructor for API Usage. Coupled with buildGrid(double scale,double *cmin,double *cmax). */
	DelPhiShared(double scale,double *coord_min,double *coord_max,const std::vector<Atom> &ats,int max_num_atoms,double domain_shrinkage,bool optimize_grids,bool eps_flag=false,bool stat_flag=false,bool multi=false,double tol_on_integer_grid=1e-2);

	// SD: PB_NEW
	/** constructor for API Usage. Coupled with buildGrid(double scale,double *cmin,double *cmax)*/
	void init(double scale,double *coord_min,double *coord_max,const std::vector<Atom> &ats,int max_num_atoms,double domain_shrinkage,bool optimize_grids,bool eps_flag,bool stat_flag,bool multi,double tol_on_integer_grid);

    /** constructor associated init*/
    void init(double scale,double perfill,const std::vector<Atom> &ats,int max_num_atoms,double domain_shrinkage,bool optimize_grids,bool eps_flag,bool stat_flag,bool multi,bool atinfo);

	/** init function-constructor*/
	void init(int max_num_atoms,double domain_shrinkage,bool optimize_grids,double scale,double perfill,string fn,bool eps_flag,bool stat_flag,bool multi,bool atinfo);
	
	/** Fortran DelPhi binding function. Allows to link DelPhi old code to the new surface framework. */
	void DelPhiBinding (double xmin,double ymin,double zmin,
						double xmax,double ymax,double zmax,
						double c1,double c2,double c3,double rmax,double perf,
						int *local_i_epsmap,int igrid,double scale,
						double *local_i_scspos,double *local_i_scsnor,
						bool *local_i_idebmap,int *local_i_ibgp,
						int maxbgp,bool status,int *atsurf);

	/**If in DelPhi binding mode, call it before denstructor. It will change the
	variables in a Fortran/DelPhi compatible mode and will avoid the denstruction of the variables
	passed to DelPhi. Always call this function ONLY just before the denstructor of the DelPhi shared
	environment. */
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
	void saveCavities2(bool onlySaveNonFilled,string sysName);
	/** Convert grid cavity points in atoms files for NanoShaper.
	This function is usefull to triangulate by NanoShaper the envelope of the cavities. */
	int cavitiesToAtoms(double rad);
	int cavitiesToAtomsWithBilevelStatusMap(double rad);
	/** given the status map of a 1.4 surf, mark for all the active cavites if they are pockets or cavitie. s*/
	void markPockets(int *external_status, vector<bool> &isPocket);
	void markPockets(int **external_bilevel_status, vector<bool> &isPocket);
	/** clear reference atoms for each cavity*/
	void clearCav2Atoms();
	/** reinit cav2atoms vector to current number of cavities. If not empty it gets cleared*/
	void initCav2Atoms();
	/** parse atinfo vector from DelPhi*/
	static void parseAtomInfo(char *atinfo,string mol,int na,double *xn1,double *rad);
	/** Save idebmap*/
	void saveIdebMap(char *fname);
	bool clearAndAllocEpsMaps();
	void clearEpsMaps();
	bool loadAtoms(string fn);
	bool loadAtoms(int na,double *pos,double *r,double *q,int *d,char *atinf);
	/** Emulates DelPhi grid construction*/
	bool buildGrid(double scale,double perfill);
	
	// SD: PB_NEW
	/** It builds a grid with side=1/scale usual but the grid is dual to the input grid.
	"Dual" means that NS grid points are indeed on edges of input grid. One assumes cmin and cmax stay on the edges. */
	bool buildGrid(double scale,double *cmin,double *cmax,double tol_on_integer_grid);

	bool getDelphiBinding()
	{
		return delphiBinding;
	}

	void buildEpsmap(bool a)
	{
		buildEpsMap = a;
	}

	void buildStatusMap(bool a)
	{
		buildStatus = a;
	}

	bool getMultiDiel()
	{
		return multi_diel;
	}

	int getNumAtoms()
	{
		return atoms.size();
	}

	// returns the lower and upper bounds of the grid
	void getBounds(double *cmin, double *cmax);

	virtual ~DelPhiShared();

	///////////////////////////////////////////////////// variables //////////////////////////////////////

	// multi dielectric mode
	bool multi_diel;
	double scale;
	double perfill;
	bool delphiBinding;
	bool buildEpsMap;
	bool buildStatus;
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
	int64_t nx,ny,nz;
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
	int *status;
	int **bilevel_status;
	// This temporary variable is introduce to detect if a cavity has been split
	// in two subcavities during Connolly filtering. Is used by the difference
	// function
	int *tempStatus;
	int **bilevel_tempStatus;

	double rmaxdim;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double baricenter[3];
	// atoms vector. Pointer to atom objects
	// Atom **atoms;
	std::vector<Atom> atoms;
	// num atoms
	// int numAtoms;
	int maxNumAtoms;
	double domainShrinkage;
	int optimizeGrids;
	char file[BUFLEN];
	double hside;
	/** dynamical vector of cavities */
	vector<vector<int*>*> *cavitiesVec;
	/** dynamical vector of cavities size*/
	vector<double> cavitiesSize;
	/** dynamical vector of cavities filling. true if cav is filled*/
	vector<bool> cavitiesFlag;
	/** dynamical vector whose index is the cavity and contains the set of atoms that produce that cavity*/
	vector<set<int>*> cav2atoms;
}; 

#endif
