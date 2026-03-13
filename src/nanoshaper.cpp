
#include "nanoshaper.h"
#include "DelphiShared.h"
#include "Surface.h"

using namespace NS;

#define SURF (static_cast<Surface*>(surf))
#define CONFIG (static_cast<ConfigFile*>(cf))
#define DS (static_cast<DelPhiShared*>(ds))

#ifdef ENABLE_CGAL 
	#ifdef CGAL_LINKED_WITH_TBB
		# include <tbb/global_control.h>
	#endif
#endif


NanoShaper::NanoShaper( const std::vector<Atom> &ats,
                        surface_type flag,
                        const double _parameter,
                        const double _stern_layer_thickness,
                        const unsigned nt,
                        const std::string *configFile)
{
    initConstructor(ats,flag,_parameter,_stern_layer_thickness,nt,configFile);
}


void NanoShaper::initConstructor( const std::vector<Atom> &ats,
                                  surface_type flag,
                                  const double _parameter,
                                  const double _stern_layer_thickness,
                                  const unsigned num_cores,
                                  const std::string *configFile)
{
    cf = nullptr;
    ds = nullptr;
    surf = nullptr;
    collectGridRays = false;
    optimizeGrids = true;
    grid_tol = 1e-2;
    
    initConfig(configFile);
    // copy the atoms internally
    atoms = ats;
    // unknown panel
    currentPanel = -1;
    // override with local info
    CONFIG->add<int>("Number_thread", num_cores);

    int numThreads = 1;

    if (num_cores <= 0)
    {
        numThreads = (int)std::thread::hardware_concurrency();
        cout << endl << INFO << "Detected " << numThreads << " logical cores";
        cout << endl << INFO << "Setting " << numThreads << " threads";
    }
    else
    {
        cout << endl << INFO << "User selected num threads " << num_cores;
        numThreads = num_cores;
    }
    CONFIG->add<int>("Number_thread", numThreads);
    // We use this instead of previous global variable "num_cores"

    conf.numThreads = numThreads;

    #if defined(ENABLE_CGAL) && defined(CGAL_LINKED_WITH_TBB)
    initTBB (numThreads);
    #endif

    switch(flag)
    {
        case skin:
        {
            CONFIG->add<std::string>("Surface","skin");
            CONFIG->add<double>( "Skin_Surface_Parameter",_parameter);
            break;
        }
        case ses:
        {
            CONFIG->add<std::string>("Surface","ses");
            CONFIG->add<double>( "Probe_Radius",_parameter);
            break;
        }
        default:
            std::cout << ERR << "Unknown surface type. Existing\n"; std::exit(EXIT_FAILURE);
    }
}


/**it returns +1 for outside, -1 for inside, 0 for unknown (grid not computed)*/
int NanoShaper::getGridPointStatus(unsigned ix,unsigned iy,unsigned iz)
{      
    if (!optimizeGrids)
    {
        if (ix<0 || iy<0 || iz<0 || ix>=DS->nx || iy>=DS->ny || iz>=DS->nz || ds==nullptr || DS->status==nullptr)
            return 0;
        // -1 -> point is inside
        // 2 -> point is temporary outside (may be really out or in cavity)
        // 3 -> point is really out
        int val = read3DVector<int>(DS->status,ix,iy,iz,DS->nx,DS->ny,DS->nz);

        if (val < 0)
            return -1;
        else
            return +1;
    }
    else
    {
        if (ix<0 || iy<0 || iz<0 || ix>=DS->nx || iy>=DS->ny || iz>=DS->nz || ds==nullptr || DS->bilevel_status==nullptr)
            return 0;
        // -1 -> point is inside
        // 2 -> point is temporary outside (may be really out or in cavity)
        // 3 -> point is really out
        int val = readBilevelGrid<int>(DS->bilevel_status,STATUS_POINT_TEMPORARY_OUT,ix,iy,iz,DS->nx,DS->ny,DS->nz);

        if (val < 0)
            return -1;
        else
            return +1;
    }
}


/**it returns the grid size*/
bool NanoShaper::getGridSize(std::vector<unsigned> &p)
{
    if (!optimizeGrids)
    {
        if (ds == nullptr || DS->status == nullptr)
            return false;
    }
    else
    {
        if (ds == nullptr || DS->bilevel_status == nullptr)
            return false;
    }
    p.clear();
    p.push_back(DS->nx);
    p.push_back(DS->ny);
    p.push_back(DS->nz);

    return true;
}


/**it returns false if the point cannot be accessed*/
bool NanoShaper::getGridPointCoordinates(unsigned ix,unsigned iy,unsigned iz,std::vector<double> &p)
{
    if (!optimizeGrids)
    {
        if (ix<0 || iy<0 || iz<0 || ix>=DS->nx || iy>=DS->ny || iz>=DS->nz || ds==nullptr || DS->status==nullptr)
            return false;
    }
    else
    {
        if (ix<0 || iy<0 || iz<0 || ix>=DS->nx || iy>=DS->ny || iz>=DS->nz || ds==nullptr || DS->bilevel_status==nullptr)
            return false;
    }
    p.clear();
    p.push_back(DS->x[ix]);
    p.push_back(DS->y[iy]);
    p.push_back(DS->z[iz]);
    return true;
}


bool NanoShaper::setDirection(unsigned direction)
{
    if (direction > 2)
        return false;
        
    if (surf == nullptr)
        return false;
                       
    currentPanel = panelResolver(direction);
    // set the panel
    SURF->panel = currentPanel;
    //every time we change direction we need to prepare
    //the right data structures.
    SURF->preProcessPanel();

    return true;
}


int NanoShaper::panelResolver(unsigned direction)
{
    // X axis
    if (direction == 0)
        return 0;
    // Z axis
    else if (direction == 2)
        return 1;
    // Y axis
    else if (direction == 1)
        return 2;
    else
    {
        std::cout << ERR << "Please select a valid direction to cast the ray 0=x,1=y,2=z";
        return -1;
    }
}


void NanoShaper::pointResolver(double *pb,double coord)
{
    if (currentPanel == 0)
        pb[0] = coord;
    else if (currentPanel == 1)
        pb[2] = coord;
    else if (currentPanel=2)
        pb[1] = coord;
}


short NanoShaper::castAxisOrientedRay(double *pa, double coord, std::vector<std::pair<VERTEX_TYPE,VERTEX_TYPE*>> &intersections, unsigned direction, bool computeNormals)
{
    double pb[3];

    int desiredPanel = panelResolver(direction);

    if (desiredPanel < 0)
    {
        std::cout << ERR << "Wrong requested ray direction";
        return -1;
    }

    if (currentPanel != desiredPanel)
    {
        std::cout << ERR << "Cannot cast a ray in a direction that is not the current one. Please use setDirection to set current direction before casting rays";
        return -1;
    }

    // copy pa to pb
    pb[0] = pa[0];
    pb[1] = pa[1];
    pb[2] = pa[2];

    // update the correct coordinate based on the current panel
    pointResolver(pb, coord);

    double diff[3];
    diff[0] = pb[0] - pa[0];
    diff[1] = pb[1] - pa[1];
    diff[2] = pb[2] - pa[2];

    if (pa[direction] >= coord)
    {
        std::cout << ERR << "Please cast rays only in progressive directions";
        return -1;
    }

    if (surf == nullptr)
    {
        std::cout << ERR << "Please compute a surface before casting rays";
        return -1;
    }

    std::vector<std::pair<VERTEX_TYPE,VERTEX_TYPE*>> intersecs;

	SURF->getRayIntersection(pa,pb,intersecs,computeNormals,0);

    if (intersecs.size() == 0)
		return -1;

	double lastValid = INFINITY;

    // performs a cleanup of the intersections
    // too near ones are removed
    for (std::vector<std::pair<VERTEX_TYPE,VERTEX_TYPE*>>::iterator itt = intersecs.begin(); itt!=intersecs.end(); itt++)
    {
        double a = itt->first;
        a = pa[direction] + diff[direction]*a;

        if (fabs(a-lastValid) >= EPS_INT)
        {
            std::pair<VERTEX_TYPE,VERTEX_TYPE*> pp;
            // put the intersection coordinate
            pp.first = a;
            pp.second = itt->second;
            // keep the intersection
            intersections.push_back(pp);
            lastValid = a;
        }
        else
        {
            // discard the intersection and delete the normal if present
            if (computeNormals)
                delete itt->second;
        }
    }
    return 1;
}


bool NanoShaper::buildAnalyticalSurface(double randDisplacement)
{   
    if (ds != nullptr)
        delete DS;
 
    if (surf != nullptr)
        delete SURF;
                
    DelPhiShared *dsl = nullptr;

    if (!CONFIG->keyExists("PB_grid_mode") || !CONFIG->read<bool>("PB_grid_mode"))
    {
        dsl = new DelPhiShared(CONFIG->read<double>("Grid_scale"),
                               CONFIG->read<double>("Grid_perfil"),
                               atoms,
                               CONFIG->read<int>("Max_Num_Atoms"),
                               CONFIG->read<double>("Domain_Shrinkage"),
                               CONFIG->read<bool>("Optimize_Grids"),
                               CONFIG->read<bool>("Build_epsilon_maps"),
                               CONFIG->read<bool>("Build_status_map"),
                               CONFIG->read<bool>("Multi_Dielectric"));
    }
    else 
    {
        double coord_min[3],coord_max[3];

        coord_min[0] = CONFIG->read<double>("xmin");
        coord_min[1] = CONFIG->read<double>("ymin");
        coord_min[2] = CONFIG->read<double>("zmin");

        coord_max[0] = CONFIG->read<double>("xmax");
        coord_max[1] = CONFIG->read<double>("ymax");
        coord_max[2] = CONFIG->read<double>("zmax");

        dsl = new DelPhiShared(CONFIG->read<double>("Grid_scale"),
                               coord_min,coord_max,
                               atoms,
                               CONFIG->read<int>("Max_Num_Atoms"),
                               CONFIG->read<double>("Domain_Shrinkage"),
                               CONFIG->read<bool>("Optimize_Grids"),
                               CONFIG->read<bool>("Build_epsilon_maps"),
                               CONFIG->read<bool>("Build_status_map"),
                               CONFIG->read<bool>("Multi_Dielectric"),this->grid_tol);
    }

    ds = dsl;

    /*
    ns4.setConfig<bool>("Save_Status_map",true);

    vector<double> p;
    getGridPointCoordinates(0,0,0,p);

    std::cout << "\n" << p[0];
    std::cout << "\n" << p[1];
    std::cout << "\n" << p[2];

    exit(-1);
    */

    // Get surface from the factory
    Surface *ss = surfaceFactory().create(CONFIG,DS);
	surf = ss;

    if (randDisplacement >= 0)
        ss->setRandDisplacement(randDisplacement);

    return ss->build();
}


bool NanoShaper::triangulate(double *surf_area)
{
    if (surf != nullptr && ds != nullptr)
    {
        if (CONFIG->read<bool>("Smooth_Mesh"))
            *surf_area = SURF->triangulateSurface(false,false);
        else
            *surf_area = SURF->triangulateSurface(true,true);

        if (*surf_area == 0)
            return false;

        if (CONFIG->read<bool>("Smooth_Mesh"))
            SURF->smoothSurface(true,true);
        return true;
    }
    return false;
}


bool NanoShaper::getNearestGridPoint(double p[3],std::vector<unsigned> &indices)
{
    if (ds == nullptr)
        return false;
    
    if (!optimizeGrids)
    {
        if (DS->status == nullptr)
            return false;
    }
    else
    {
        if (DS->bilevel_status == nullptr)
            return false;
    }
        
    indices.clear();
    
    int ix = (int)rintp((p[0]-DS->xmin)/DS->side);
    int iy = (int)rintp((p[1]-DS->ymin)/DS->side);
    int iz = (int)rintp((p[2]-DS->zmin)/DS->side);
    
    bool outside = false;
    
    if (ix < 0)
    {
        ix = 0;
        outside = true;
    }
    if (ix >= DS->nx)
    {
        ix = DS->nx-1;
        outside = true;
    }

    if (iy < 0)
    {
        iy = 0;
        outside = true;
    }

    if (iy >= DS->ny)
    {
        iy = DS->ny-1;
        outside = true;
    }

    if (iz < 0)
    {
        iz = 0;
        outside = true;
    }
    if (iz >= DS->nz)
    {
        iz = DS->nz-1;
        outside = true;
    }
        
    indices.push_back((unsigned)ix);
    indices.push_back((unsigned)iy);
    indices.push_back((unsigned)iz);
                
    return outside;
}


/** read config var*/        
template<class T> T NanoShaper::readConfig(const std::string &keyword)
{
    return CONFIG->read<T>(keyword);
}


/** set config var*/
template<class T> void NanoShaper::setConfig(std::string keyword,const T val)
{
    CONFIG->add<T>(keyword,val);
}


template bool NanoShaper::readConfig<bool>(const std::string &keyword);
template int NanoShaper::readConfig<int>(const std::string &keyword);
template double NanoShaper::readConfig<double>(const std::string &keyword);
template float NanoShaper::readConfig<float>(const std::string &keyword);
template unsigned NanoShaper::readConfig<unsigned>(const std::string &keyword);
template short NanoShaper::readConfig<short>(const std::string &keyword);
template long NanoShaper::readConfig<long>(const std::string &keyword);
template std::string NanoShaper::readConfig<std::string>(const std::string &keyword);
    
template void NanoShaper::setConfig<bool>(std::string keyword,const bool val);
template void NanoShaper::setConfig<int>(std::string keyword,const int val);
template void NanoShaper::setConfig<double>(std::string keyword,const double val);
template void NanoShaper::setConfig<float>(std::string keyword,const float val);
template void NanoShaper::setConfig<unsigned>(std::string keyword,const unsigned val);
template void NanoShaper::setConfig<short>(std::string keyword,const short val);
template void NanoShaper::setConfig<long>(std::string keyword,const long val);
template void NanoShaper::setConfig<std::string>(std::string keyword,const std::string val);


bool NanoShaper::colourGrid(void)
{
    double surf_volume;

    return colourGrid(&surf_volume);
}


bool NanoShaper::colourGrid(double *surf_volume)
{
    // conventional execution mode
    if (!collectGridRays)
	    return SURF->getSurf(surf_volume,
                             CONFIG->read<bool>("Optimize_Grids"),
                             CONFIG->read<bool>("Cavity_Detection_Filling"),
                             CONFIG->read<double>("Conditional_Volume_Filling_Value"));
    // new execution mode, keeps track of intersections and normals at grid points
    else
    {
        vector<packet> intersectionsInfo;
        bool cvn = CONFIG->read<bool>("Compute_Vertex_Normals");
        // impose to compute normals in order to collect them
        SURF->setComputeNormals(true);

    	bool retVal = SURF->getSurf(surf_volume,
                                    CONFIG->read<bool>("Optimize_Grids"),
                                    CONFIG->read<bool>("Cavity_Detection_Filling"),
                                    CONFIG->read<double>("Conditional_Volume_Filling_Value"),
                                    &intersectionsInfo);

        if (!retVal)
            return false;

        // clean up in case it is needed
        for (unsigned k=0; k<3; k++)
            raysArrayMap[k].clear(); 

        // it initializes the map with empty crossing_t objects such that the guest code can verify 
        // that the rays have been cast and no intersections have been found.
        // later we fill where intersections have been found

        std::vector<unsigned> grid;

        getGridSize(grid);

        //panel 0, dir x
        for (int iy = 0; iy < grid[1]; ++iy) 
        {
            for (int iz = 0; iz < grid[2]; ++iz) 
            {
                std::vector<double> p;
                std::array<double,2> point;
                unsigned dir = 0;

                getGridPointCoordinates(0, iy, iz, p);

                // y,z panel, x direction      
                point[0] = p[1];
                point[1] = p[2];
                
                raysArrayMap[dir][point] = crossings_t();
                crossings_t &temp = raysArrayMap[dir][point];
                temp.point[0] = point[0];
                temp.point[1] = point[1];
                temp.dir = dir;  
                temp.init = 1;
            }
        }

        //panel 2, dir y
        for (int ix = 0; ix < grid[0]; ++ix) 
        {
            for (int iz = 0; iz < grid[2]; ++iz) 
            {
                std::vector<double> p;
                std::array<double,2> point;
                unsigned dir = 1;

                getGridPointCoordinates(ix, 0, iz, p);

                // y,z panel, x direction      
                point[0] = p[0];
                point[1] = p[2];
                
                raysArrayMap[dir][point] = crossings_t();
                crossings_t &temp = raysArrayMap[dir][point];
                temp.point[0] = point[0];
                temp.point[1] = point[1];
                temp.dir = dir;  
                temp.init = 1;
            }
        }

        //panel 1, dir z
        for (int ix = 0; ix < grid[0]; ++ix) 
        {
            for (int iy = 0; iy < grid[1]; ++iy) 
            {
                std::vector<double> p;
                std::array<double,2> point;
                unsigned dir = 2;

                getGridPointCoordinates(ix, iy, 0, p);

                // y,z panel, x direction      
                point[0] = p[0];
                point[1] = p[1];
                
                raysArrayMap[dir][point] = crossings_t();
                crossings_t &temp = raysArrayMap[dir][point];
                temp.point[0] = point[0];
                temp.point[1] = point[1];
                temp.dir = dir;  
                temp.init = 1;
            }
        }

        unsigned num_packs = intersectionsInfo.size();
        // printf("\npacks %d", num_packs);
        auto it = intersectionsInfo.begin();

        // reorganize the data according to the required data structure
		for (unsigned i_pack=0; i_pack<num_packs; i_pack++, it++)
		{      
            packet &pack = *it;
            auto it_intersections = pack.first->begin();
            #if !defined(COORD_NORM_PACKING)
            auto it_normals = pack.second->begin();
            #endif
            unsigned psize = pack.first->size();
		    cout << endl << INFO << "Unpacking rays packet " << i_pack << " of size " << psize;
            
            // unpacking and reshaping
            for (unsigned k=0; k<psize; k++)
            {
                int dir;
                #if !defined(COORD_NORM_PACKING)
                coordVec &cv_int = *(it_intersections++);
                coordVec &cv_norm = *(it_normals++);
                #else
                #if !defined(COMPRESS_INTERSECTION_COORDS)
                coordNormPacket &cv_int = *(it_intersections++);

                dir = cv_int.dir;
                #else
                compressedCoordNormPacket &cv_int = *(it_intersections++);

                int dummy[3];

                cv_int.getCompressedCoords(&dummy[0], &dummy[1], &dummy[2], &dir);
                #endif
                #endif

                std::array<double,2> point;

                // y,z panel, x direction
                if (dir == 0)
                {    
                   // point = mappa(array(cv_int.iy,cv_int.iz));
                   point[0] = cv_int.vec[1];
                   point[1] = cv_int.vec[2];
                }

                // x,z panel, y direction
                else if (dir == 1)
                {                                                
                    point[0] = cv_int.vec[0];
                    point[1] = cv_int.vec[2];
                }

                // x,y panel, z direction
                else if (dir == 2)
                {                    
                    point[0] = cv_int.vec[0];
                    point[1] = cv_int.vec[1];
                }

                // printf("\n%f %f %f , dir %d",cv_int.vec[0],cv_int.vec[1],cv_int.vec[2],dir);

                crossings_t &crossings = raysArrayMap[dir][point];
    
                // store the coordinate of the intersection only, not the entire vector
                crossings.inters.push_back(cv_int.vec[dir]);      
                // store the normal
                #if !defined(COORD_NORM_PACKING)
                crossings.normals.push_back(cv_norm.vec[0]);
                crossings.normals.push_back(cv_norm.vec[1]);
                crossings.normals.push_back(cv_norm.vec[2]);
                #else
                crossings.normals.push_back(cv_int.nor[0]);
                crossings.normals.push_back(cv_int.nor[1]);
                crossings.normals.push_back(cv_int.nor[2]);
                #endif

                // TODO cleanup the intersection vector and normal, as by now is no more useful
            } // for intersections
        } // for packet
        // restore the previous configuration on normals computation
        SURF->setComputeNormals(cvn);
        return true;
    }
}


void NanoShaper::clean()
{
    if (CONFIG != nullptr)
        delete CONFIG;

    if (DS != nullptr)
        delete DS;

    if (SURF != nullptr)
        delete SURF;
}


void NanoShaper::initConfig(const std::string *confFile)
{
	NanoShaper::clean();
    ConfigFile *cfl;
    
    // apply default configuration values. These values are particularly
    // well suite to colour a coarse grid
	if (confFile == nullptr)
	{
        cfl = new ConfigFile();
        cfl->add<std::string>("Surface", "skin");
        cfl->add<std::string>("Output_Folder", "");
        cfl->add<bool>("Build_epsilon_maps", false);
        cfl->add<double>("Probe_Radius", 1.4);
        cfl->add<double>("Skin_Surface_Parameter", 0.45);
        cfl->add<bool>("Save_eps_maps", false);
        cfl->add<bool>("Save_Status_map", false);
        cfl->add<bool>("Build_status_map", true);
        cfl->add<bool>("Save_ideb_map", false);
        cfl->add<bool>("Cavity_Detection_Filling", false);
        cfl->add<bool>("Project_boundary_grid_points", false);
        cfl->add<bool>("Accurate_Triangulation", false);
        cfl->add<bool>("Triangulation", false);
        cfl->add<std::string>("Operative_Mode", "normal");
        cfl->add<bool>("Debug_Internals", false);
        cfl->add<bool>("Patch_Based_Algorithm", true);
        cfl->add<bool>("Analytical_Ray_Vs_Torus_Intersection", true);
        cfl->add<bool>("Collect_face_intersections", false);
        cfl->add<double>("Parallel_Buildup_Halo_Thickness", 0.0);
        cfl->add<bool>("Force_Serial_Build", false);
        cfl->add<bool>("Optimize_Grids", true);
        cfl->add<int>("Max_Num_Atoms", -1);
        cfl->add<double>("Domain_Shrinkage", 0.);
        cfl->add<double>("Conditional_Volume_Filling_Value", 11.4);
        cfl->add<int>("Num_Wat_Pocket", 2);
        cfl->add<double>("Grid_scale", 1.0);
        cfl->add<double>("Grid_perfil", 80.0);
        cfl->add<std::string>("XYZR_FileName", "null");
        cfl->add<bool>("Multi_Dielectric", false);
        cfl->add<bool>("Smooth_Mesh", true);
        cfl->add<bool>("Tri2Balls", false);
        cfl->add<bool>("Save_eps_maps", false);
        cfl->add<bool>("Save_bgps", false);
        cfl->add<bool>("Save_Cavities", false);
        cfl->add<std::string>("Sys_Name", "mol_from_api");
        cfl->add<int>("Number_thread", 1);
        cfl->add<bool>("Print_Available_Surfaces", false);

        int seed = 1;
        cfl->add<int>("Seed", seed);
        srand(seed);

        cfl->add<bool>("Pockets_And_Cavities", true);
        cfl->add<bool>("Link_Pockets", false);
        cfl->add<double>("Pocket_Radius_Big", 3.0);
        cfl->add<double>("Pocket_Radius_Small", 1.4);
        cfl->add<double>("Pocket_Radius_Link", 1.0);
        cfl->add<double>("Memb_Height", 25.0);
        cfl->add<double>("Memb_Shift", 1.0);
        cfl->add<std::string>("Root_FileName", "");
        cfl->add<bool>("Compute_Vertex_Normals", false);
        cfl->add<bool>("Save_Mesh_MSMS_Format", false);
        cfl->add<bool>("Save_Mesh_PLY_Format", false);
        cfl->add<bool>("Load_Balancing", true);
        cfl->add<double>("Blobbyness", -2.5);
        cfl->add<std::string>("Surface_File_Name", "triangulatedSurf.off");
        cfl->add<bool>("Keep_Water_Shaped_Cavities", false);
        cfl->add<int>("Max_Probes_Self_Intersections", 100);
        cfl->add<int>("Self_Intersections_Grid_Coefficient", 1.5);
        cfl->add<bool>("Check_duplicated_vertices", true);
        cfl->add<bool>("Save_PovRay", false);
        cfl->add<int>("Max_mesh_auxiliary_grid_size", 100);
        cfl->add<int>("Max_mesh_patches_per_auxiliary_grid_cell", 250);
        cfl->add<int>("Max_mesh_auxiliary_grid_2d_size", 100);
        cfl->add<int>("Max_mesh_patches_per_auxiliary_grid_2d_cell",  250);
        cfl->add<int>("Max_ses_patches_auxiliary_grid_size", 100);
        cfl->add<int>("Max_ses_patches_per_auxiliary_grid_cell", 400);
        cfl->add<int>("Max_ses_patches_auxiliary_grid_2d_size", 50);
        cfl->add<int>("Max_ses_patches_per_auxiliary_grid_2d_cell", 400);
        cfl->add<int>("Max_skin_patches_auxiliary_grid_size", 100);
        cfl->add<int>("Max_skin_patches_per_auxiliary_grid_cell", 400);
        cfl->add<int>("Max_skin_patches_auxiliary_grid_2d_size", 50);
        cfl->add<int>("Max_skin_patches_per_auxiliary_grid_2d_cell", 400);

        // new grid mode keywords
        cfl->add<bool>("PB_grid_mode",false);
        cfl->add<double>("xmin",-1.0);
        cfl->add<double>("ymin",-1.0);
        cfl->add<double>("zmin",-1.0);
        cfl->add<double>("xmax",-1.0);
        cfl->add<double>("ymax",-1.0);
        cfl->add<double>("zmax",-1.0);
	}
	else
	{
        // get data from configuration file
        try
        {
            cfl = new ConfigFile(confFile->c_str());
        }
        catch(...)
        {
            cout << endl << ERR << "Cannot read " << confFile << endl;
            exit(-1);
        }    
    }   
    cf = (void *)cfl;
}
