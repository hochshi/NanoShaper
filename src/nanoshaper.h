#ifndef nanoshaper_h__
#define nanoshaper_h__

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include <math.h>
#include "raytracer_datatype.h"


// #define FLOAT_VERTICES

#if !defined(FLOAT_VERTICES)
#define VERTEX_TYPE double
#else
#define VERTEX_TYPE float
#endif


namespace NS {
    
enum surface_type
{
	ses = 0,
	skin = 1,
	blobby = 2
};

    
/**@brief A point in R^3*/
class Point
{
public:
	double pos[3];

	Point()
	{
		pos[0] = 0;
		pos[1] = 0;
		pos[2] = 0;
	}

	Point(const double x,const double y,const double z)
	{
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
	}

	Point(const Point &p)
	{
		pos[0] = p.pos[0];
		pos[1] = p.pos[1];
		pos[2] = p.pos[2];
	}
	
	double squaredDist(const Point *const x)
	{
		double acc = 0;

		for (int i=0; i<3; i++)
		{
			double t = x->pos[i] - pos[i];
			acc += t*t;
		}
		return acc;
	}

	double dot(const Point *const x)
	{
		double acc = 0;

		for (int i=0; i<3; i++)
			acc += x->pos[i]*pos[i];

		return acc;
	}

	double norm()
	{
		double acc = 0;

		for (int i=0; i<3; i++)
			acc += pos[i]*pos[i];

		return sqrt(acc);
	}
	
	double dist(const Point *const x)
	{
		double acc = 0;

		for (int i=0; i<3; i++)
		{
			double t = x->pos[i] - pos[i];
			acc += t*t;
		}
		return sqrt(acc);
	}

	void sumInPlace(const Point *const x)
	{
		for (int i=0; i<3; i++)
		{
			pos[0] += x->pos[0];
			pos[1] += x->pos[1];
			pos[2] += x->pos[2];
		}
	}

	Point *sum(const Point *const x)
	{
		Point *p = new Point();

		for (int i=0; i<3; i++)
		{
			p->pos[0] = pos[0] + x->pos[0];
			p->pos[1] = pos[1] + x->pos[1];
			p->pos[2] = pos[2] + x->pos[2];
		}
		return p;
	}
};


class AtomInfo
{
public:

	std::string name;
	int serial;
	int resNum;
	std::string resName;
	std::string chain;

	AtomInfo()
	{
		name = "X";
		serial = -1;
		resNum = -1;
		resName = "UNK";
		chain = "X";
	}

	AtomInfo(const std::string name,const int resNum,const std::string resName,const std::string chain)
	{
		this->name = name;
		this->serial = -1;
		this->resNum = resNum;
		this->resName = resName;
		this->chain = chain;
	}

	const std::string &getName()
	{
		return name;
	}

	void setName(const std::string name)
	{
		this->name = name;
	}

	const int &getSerial()
	{
		return serial;
	}

	void setSerial(const int s)
	{
		serial = s;
	}

	const std::string &getResName()
	{
		return resName;
	}

	void setResName(const std::string res)
	{
		resName = res;
	}

	const int &getResNum()
	{
		return resNum;
	}

	void setResNum(const int res)
	{
		resNum = res;
	}

	void setChain(const std::string ch)
	{
		chain = ch;
	}

	const std::string &getChain()
	{
		return chain;
	}
	
	/** copy constructor */
	AtomInfo(const AtomInfo &p)
	{
		name = p.name;
		serial = p.serial;
		resNum = p.resNum;
		resName = p.resName;
		chain = p.chain;
	}
    
    /** copy constructor */
	AtomInfo(AtomInfo &p)
	{
		name = p.name;
		serial = p.serial;
		resNum = p.resNum;
		resName = p.resName;
		chain = p.chain;
	}
};

/**@brief An atom: an atom is a point in R^3 with radius, charge and dielectric value */
class Atom: public Point
{
public:
	double radius;
	double radius2;
	double charge;
	int dielectric;
	AtomInfo ai;

	Atom():Point()
	{
		radius = 0;
		charge = 0;
		dielectric = 0;
	}
	
	Atom(const double x,const double y,const double z,const double r,const double c,const int d):Point(x,y,z)
	{
		radius = r;
		radius2 = r*r;
		charge = c;
		dielectric = d;
	}

	Atom(const double x,const double y,const double z,const double r,const double c,const int d,
    const std::string name,const std::string resName,const int resNum,const std::string chain):Point(x,y,z)
	{
		radius = r;
		radius2 = r*r;
		charge = c;
		dielectric = d;
		ai.setName(name);
		ai.setResName(resName);
		ai.setChain(chain);
		ai.setResNum(resNum);
	}
	
	// copy constructor
	Atom(const Atom &p):Point(p)
	{
		radius = p.radius;
		charge = p.charge;
		dielectric = p.dielectric;        
		ai.setName(p.ai.name);
		ai.setResName(p.ai.resName);
		ai.setChain(p.ai.chain);
		ai.setResNum(p.ai.resNum);
		ai.setSerial(p.ai.serial);
	}

	void print()
	{
		std::cout << std::endl << std::endl << "Atom Info:" << std::endl;
		std::cout << "\tname -> " << ai.getName() << std::endl;
		std::cout << "\tresidue name -> " << ai.getResName() << std::endl;
		std::cout << "\tresidue number -> " << ai.getResNum() << std::endl;
		std::cout << "\tchain -> " << ai.getChain() << std::endl;
		std::cout << "\tposx -> " << pos[0] << std::endl ;
		std::cout << "\tposy -> " << pos[1] << std::endl ;
		std::cout << "\tposz -> " << pos[2] << std::endl ;
		std::cout << "\tcharge -> " << charge << std::endl ;
		std::cout << "\tradius -> " << radius << std::endl ;
		std::cout << "\tdielectric -> " << dielectric;
	}

};

class NanoShaper
{
private:
    /**  global data structure of configuration data */
    void *cf;
    /** delphi shared object */
    void *ds;
    /** surface object */
    void *surf;
    /**set of current atoms*/
    std::vector<Atom> atoms;
    int currentPanel;    
    int panelResolver(unsigned direction);
    void pointResolver(double *pb,double coord);
    /** if enabled NS returns also from get surf the set of grid rays intersections and normals */
    bool collectGridRays;
    /** map which contains the set of intersections, divided per panel */
    rays_t raysArrayMap;
    bool optimizeGrids;
    /* tolerance on the NextGenPB grid check*/
    double grid_tol;
    
    
public :

    /** This routine accepts the atoms set, the type of surface, 
    the surface parameter, the stern layer thickness and the number of threads. 
    Optionally it accepts also a config file. If the file 
    is found, params from the file are used. The parameters passed here
    however always override the one from the configuration file.
    This function will call and delegate to initConstructor */
    NanoShaper (const std::vector<Atom>& atoms,
                surface_type flag,
                const double _parameter,
                const double _stern_layer_thickness,
                const unsigned nt=1,
                const std::string *configFile=nullptr);

    /**The real "constructor". Allows to build later to declaration of the object */
    void initConstructor (const std::vector<Atom> &atoms,
                surface_type flag,
                const double _parameter,
                const double _stern_layer_thickness,
                const unsigned num_cores=1,
                const std::string* configFile=nullptr);

    /**Cast a ray from a starting to an ending point and returns the set of encountered
    intersections and normals along the ray path
    By default it does not compute the normals. Rays can be cast only along axis oriented and positive directions; 
    the current direction has be set through the direction argument.
    direction = 0 -> along x axis, coord is interpreted as the x value of the final point
    direction = 1 -> along y axis, coord is interpreted as the y value of the final point
    direction = 2 -> along z axis, coord is interpreted as the z value of the final point
    It returns a status:
    -1 no surf available or coord<=pa[direction] or not valid direction or other issues
    +1 intersections are returned, possibly 0.
    Please consider that no parity check is done here; the only clean-up done is that pathologically near
    intersections are removed as done internally by nanoshaper. No assumptions is made on the parity, is up
    to the caller to perfoms consistency checks where needed. */
    short castAxisOrientedRay(double *pa,double coord,std::vector<std::pair<VERTEX_TYPE,VERTEX_TYPE*> > &inters,unsigned direction,bool computeNormals=false);

    /**Computes the analytical rapresentation of the surface without performing grid colouring, nor triangulation*/
    bool buildAnalyticalSurface(double randDisplacement=-1);  
    
    /**Colour the NS internal grids with the in/out information*/
    bool colourGrid(void);
	bool colourGrid(double *surf_volume);
    
    /**Runs the triangulation*/
    bool triangulate(double *surf_area);
    
    /**Clean up memory*/
    void clean();
    
    /** read config var*/        
    template<class T> T readConfig(const std::string &keyword);
    /** set config var*/
    template<class T> void setConfig(std::string keyword,const T val);              
    
    /** nearestGridPoint (with saturation), returns false in case of failure or saturation.
    Saturation means p is outside the grid but indices are still the nearest ones */
    bool getNearestGridPoint(double p[3],std::vector<unsigned> &indices);
            
    /**it returns +1 for outside, -1 for inside, 0 for unknown (grid not computed) */
    int getGridPointStatus(unsigned ix,unsigned iy,unsigned iz);
    
    /**it returns false if the point cannot be accessed */
    bool getGridPointCoordinates(unsigned ix,unsigned iy,unsigned iz,std::vector<double> &p);
    
    /**it returns the grid size */
    bool getGridSize(std::vector<unsigned> &p);
    
    /**It parses a config file. If nullptr is passed it populates the
    conf variable with default values */
    void initConfig(const std::string *filePath=nullptr);
    
    /**call it everytime the raycasting direction is changed (included the first time) */
    bool setDirection(unsigned direction);

    /** set if during grid colouring rays intersections and normals have to be retrieved */
    void setCollectGridRays(bool cgr)
    {
        collectGridRays = cgr;
    }

    /** can be used to set a tolerance on the NextGenPB grid check. To have effected it has to be called
    before the buildAnalytical surface function which also builds the grid. */
    void setGridTol(double tol)
    {
        grid_tol = tol;
    }

    /** set if many grids will be optimized or not */
    void setOptimizeGrids(bool optimize_grids)
    {
		optimizeGrids = optimize_grids;
    }

    /** returns the map where we store the grid and normals intersections */
    rays_t &getRaysMap()
    {
        return raysArrayMap;
    }

    NanoShaper()
    {
        currentPanel = -1;        
    }    
    ~NanoShaper () = default;

};
}

#endif  // nanoshaper_h__
