#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include<vector>

#define CLIENT_INFO " <<CLIENT INFO>>"

#include "../src/nanoshaper.h"

/**utility function to print intersection map*/
void print_map(const std::array<std::map<std::array<double, 2>, crossings_t, map_compare>, 3>& r);

// in this main code we present a set of example usage of the NS api

int main(int argc, char* argv[])
{
    NS::surface_type surf_type = NS::skin;
    NS::surface_type surf_type2 = NS::ses;
    double skin_param = 0.45;
    double probe = 0.01;
    double stern_layer = 2.;
    double radius = 2.0;
    double charge = 1;
    double dielectric=0;
    unsigned numberOfThreads = 1;
    // analytical test with a single sphere of radius 1 and scale 1
    // disable random displacement
    double rand_displ = 0.1;
    unsigned X_DIR=0,Y_DIR=1,Z_DIR=2;

    double surf_volume = 0.0;
    double surf_area = 0.0;

    // EXAMPLE 1, generic usage example
    /*
    std::vector<NS::Atom> atoms;

    //4 spheres, actually one
    atoms.push_back(NS::Atom(5.1, 10.1,1.0, radius,charge,dielectric));
    atoms.push_back(NS::Atom(10.2,10.2,1.7, radius,charge,dielectric));
    atoms.push_back(NS::Atom(5.3, 5.3,-0.1, radius,charge,dielectric));
    atoms.push_back(NS::Atom(10.4,5.4,-0.5, radius,charge,dielectric));


    NS::NanoShaper ns(atoms,surf_type,skin_param,stern_layer,numberOfThreads);

    ns.setConfig<double>("Grid_scale", 3.0 );
    // 50 is needed for small number of atoms
    ns.setConfig<double>("Grid_perfil", 50.0 );

    ns.buildAnalyticalSurface(rand_displ);

    unsigned ix = 2;
    unsigned iy = 1;
    unsigned iz = 3;
    std::vector<unsigned> gridSize;
    std::vector<double> coords;

    bool gridReady = ns.getGridSize(gridSize);
    std::cout << std::endl << CLIENT_INFO << " Grid is ready " << gridReady;

    short val =  ns.getGridPointStatus(ix,iy,iz);
    std::cout << std::endl << CLIENT_INFO << " Grid val " << val;

    // colour coarse grained grid (status map)
    ns.colourGrid(&surf_volume);

    gridReady = ns.getGridSize(gridSize);

    std::cout << std::endl << CLIENT_INFO << " Grid is ready2 " << gridReady;

    if (gridReady)
    {
        std::cout << std::endl << CLIENT_INFO << " Grid size " << gridSize[0] << " " << gridSize[1] << " " << gridSize[2];

        val =  ns.getGridPointStatus(ix,iy,iz);
        std::cout << std::endl << CLIENT_INFO << " Grid val " << val;

        ns.getGridPointCoordinates(ix,iy,iz,coords);
        std::cout << std::endl << CLIENT_INFO << " Coords " << coords[0] << " " << coords[1] << " " << coords[2];

        ns.triangulate(&surf_area);
        // few rays
        system("cp triangulatedSurf.off lowRes.off");
    }

    // same grid, more rays, more accurate triangulation
    NS::NanoShaper ns2(atoms,surf_type,skin_param,stern_layer,numberOfThreads);
    ns2.setConfig<double>("Grid_scale",3.0);
    ns2.setConfig<double>("Grid_perfil",50.0);
    ns2.setConfig<bool>("Accurate_Triangulation",true);
    ns2.buildAnalyticalSurface(rand_displ);
    ns2.colourGrid(&surf_volume);
    ns2.triangulate(&surf_area);
    system("cp triangulatedSurf.off highRes.off");

    // play with grid and rays
    double startp[3]={7.5,2,0};
    double offset = 11;
    double endp[3]={startp[0],startp[1]+offset,startp[2]};
    unsigned y_direction = 1;
    bool computeNormals = false;
    std::vector<unsigned> indices;
    ns2.getNearestGridPoint(startp,indices);
    std::cout << std::endl << CLIENT_INFO << " Start point on the grid " << indices[0] << " " << indices[1] << " " << indices[2];
    short stat = ns2.getGridPointStatus(indices[0],indices[1],indices[2]);
    int startx = indices[0];
    int starty = indices[1];
    int startz = indices[2];

    std::cout << std::endl << CLIENT_INFO << " Status " << stat;
    ns2.getNearestGridPoint(endp,indices);
    int endy = indices[1];

    std::cout << std::endl << CLIENT_INFO << " End point on the grid " << indices[0] << " " << indices[1] << " " << indices[2];
    stat = ns2.getGridPointStatus(indices[0],indices[1],indices[2]);
    std::cout << std::endl << CLIENT_INFO << " Status " << stat;

    for (int i=starty;i<=endy;i++)
    {
        stat = ns2.getGridPointStatus(startx,i,startz);
        std::cout << std::endl << CLIENT_INFO << " Status " << stat << " index " << i;
    }

    std::vector<std::pair<double,double*> > inters;
    // this call prepare the datastuctures. Call it if you change direction of ray casting
    // this is done outside as it is an expensive function hence if multiple rays are cast
    // in the same direction it is convenient to call it separately
    ns2.setDirection(y_direction);
    ns2.castAxisOrientedRay(startp,endp[y_direction],inters,y_direction,computeNormals);
    std::cout << std::endl << CLIENT_INFO <<"Ray path:" ;
    for (unsigned i=0;i<inters.size();i++)
        std::cout << std::endl << CLIENT_INFO << "Hit at " << inters[i].first;

    inters.clear();
    startp[0]=startp[0]-2.5;
    endp[0]=endp[0]-2.5;
    // same direction, we don't need to call setDirection again
    ns2.castAxisOrientedRay(startp,endp[y_direction],inters,y_direction,computeNormals);

    std::cout << std::endl << CLIENT_INFO << "Ray path:" ;
    for (unsigned i=0;i<inters.size();i++)
        std::cout << std::endl << CLIENT_INFO << "Hit at " << inters[i].first;

    std::cout << std::endl;

    //////////////////////////////////////////////
    // end EXAMPLE 1


    // EXAMPLE 2 : here we show how to use the NS api not to colour the grid only,
    // but also to collect all the intersections and the normals.
    // This can be used to interface NS to a PB solver which takes advantage of the analytical
    // intersections and normals or in general to any guest code which can be interested in getting
    // the intersections.

    std::vector<NS::Atom> atoms;
    atoms.push_back(NS::Atom(0,0,0, +1,charge,dielectric));
    atoms.push_back(NS::Atom(0,0,0, 0,charge,dielectric));
    atoms.push_back(NS::Atom(0,0,0, 0,charge,dielectric));
    atoms.push_back(NS::Atom(0,0,0, 0,charge,dielectric));

    NS::NanoShaper ns_grid_rays(atoms,surf_type,skin_param,stern_layer,numberOfThreads);
    ns_grid_rays.setConfig<double>("Grid_scale",2.0);
    ns_grid_rays.setConfig<double>("Grid_perfil",50.0);
    ns_grid_rays.setConfig<bool>("Accurate_Triangulation",true);
    ns_grid_rays.setConfig<bool>("Build_epsilon_maps",true);
    // we use a random diplacement equal to 0 to get do fully analytical test
    ns_grid_rays.buildAnalyticalSurface(rand_displ);
    // remember to set this to true in order to collect rays intersections data
    ns_grid_rays.setCollectGridRays(true);
    ns_grid_rays.colourGrid(&surf_volume);
    auto raysMap = ns_grid_rays.getRaysMap();

    std::cout << "\n" << CLIENT_INFO << " Along x found " << raysMap[X_DIR].size() << " active panels";
    std::cout << "\n" << CLIENT_INFO << " Along y found " << raysMap[Y_DIR].size() << " active panels";
    std::cout << "\n" << CLIENT_INFO << " Along z found " << raysMap[Z_DIR].size() << " active panels";


    // iterate through the map and print the intersections and normals
    // shows all the intersections in the x direction
    //auto x_crossings = raysMap[X_DIR];
    //auto it = x_crossings.begin();
    //while (it != x_crossings.end())
    //{
    //    auto itt = it->second;
    //    std::cout << "\n" << CLIENT_INFO << " Cast ray along x from pos y " << itt->point[0] << " pos z " << itt->point[1];
    //    for (int i=0;i<itt->inters.size();i++)
    //        std::cout << "\n\t" << CLIENT_INFO << " Intersection x val " << itt->inters[i] << " normal " << itt->normals[i][0] << "," << itt->normals[i][1] << "," << itt->normals[i][2];
    //    ++it;
    //}
    //auto y_crossings = raysMap[Y_DIR];
    //auto it = y_crossings.begin();
    //while (it != y_crossings.end())
    //{
    //    auto itt = it->second;
    //    std::cout << "\n" << CLIENT_INFO << " Cast ray along y from pos x " << itt->point[0] << " pos z " << itt->point[1];
    //    for (int i=0;i<itt->inters.size();i++)
    //        std::cout << "\n\t" << CLIENT_INFO << " Intersection y val " << itt->inters[i] << " normal " << itt->normals[i][0] << "," << itt->normals[i][1] << "," << itt->normals[i][2];
    //    ++it;
    //}

    //auto z_crossings = raysMap[Z_DIR];
    //it = z_crossings.begin();
    //while (it != z_crossings.end())
    //{
    //    auto itt = it->second;
    //    std::cout << "\n" << CLIENT_INFO << " Cast ray along z from pos x " << itt->point[0] << " pos y " << itt->point[1];
    //    for (int i=0;i<itt->inters.size();i++)
    //        std::cout << "\n\t" << CLIENT_INFO << " Intersection z val " << itt->inters[i] << " normal " << itt->normals[i][0] << "," << itt->normals[i][1] << "," << itt->normals[i][2];
    //    ++it;
    //}
    */
    // end EXAMPLE 2

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    // EXAMPLE 3 : here we test the collection of intersections + the new grid build up mode.
    // in this mode we embed an arbitrary parallelepipe grid into the NS cubic grid.
    // Then we collect and return intersections. This allows to interface to solver which
    // builds an arbitray grid. We require to know the lower and upper corner of the grid.
    /*
    std::cout << "\n\n\nNew grid mode + intersections collection";

    std::vector<NS::Atom> atoms2;

    atoms2.push_back(NS::Atom(1.000,-1.000,6.500, 2.0,charge,dielectric));
    atoms2.push_back(NS::Atom(2.000,0.000,1.000,1.50,charge,dielectric));
    atoms2.push_back(NS::Atom(0.000,1.000,-2.800,2.00,charge,dielectric));
    atoms2.push_back(NS::Atom(0.750,-1.500,-7.000,1.50,charge,dielectric));
    atoms2.push_back(NS::Atom(-2.000,-1.000,3.000,1.50,charge,dielectric));

    NS::NanoShaper ns4(atoms2,surf_type2,probe,stern_layer,numberOfThreads);

    // set here a consistent grid scale
    ns4.setConfig<double>("Grid_scale", 2.0 );
    ns4.setConfig<bool>("Accurate_Triangulation",true);
    ns4.setConfig<bool>("Build_epsilon_maps",true);

    // build the grid in the new mode
    ns4.setConfig<bool>("PB_grid_mode",true);

    // impose here the min and max of the cube
    ns4.setConfig<double>("xmin",-4.5);
    ns4.setConfig<double>("ymin",-4.0);
    ns4.setConfig<double>("zmin",-10.5);

    ns4.setConfig<double>("xmax",4.5);
    ns4.setConfig<double>("ymax",4.0);
    ns4.setConfig<double>("zmax",10.5);

    ns4.buildAnalyticalSurface(rand_displ);
    // remember to set this to true in order to collect rays intersections data
    ns4.setCollectGridRays(true);
    ns4.colourGrid(&surf_volume);
    // retrieve intersections data from the map
    auto raysMap2 = ns4.getRaysMap();

    std::cout << "\n" << CLIENT_INFO << " Along x found " << raysMap2[X_DIR].size() << " active panels";
    std::cout << "\n" << CLIENT_INFO << " Along y found " << raysMap2[Y_DIR].size() << " active panels";
    std::cout << "\n" << CLIENT_INFO << " Along z found " << raysMap2[Z_DIR].size() << " active panels";

    // iterate through the map and print the intersections and normals
    // shows all the intersections in the x direction
    auto x_crossings2 = raysMap2[X_DIR];
    auto it2 = x_crossings2.begin();
    while (it2 != x_crossings2.end())
    {
        auto itt = it2->second;
        std::cout << "\n" << CLIENT_INFO << " Cast ray along x from pos y " << itt.point[0] << " pos z " << itt.point[1];
        for (int i=0;i<itt.inters.size();i++)
            std::cout << "\n\t" << CLIENT_INFO << " Intersection x val " << itt.inters[i] << " normal " << itt.normals[3*i+0] << "," << itt.normals[3*i+1] << "," << itt.normals[3*i+2];
        ++it2;
    }

    ns4.triangulate(&surf_area);
    */
    // end EXAMPLE 3

    // EXAMPLE 3 bis
    std::cout << "\n\n\nNew grid mode + intersections collection";
    std::vector<NS::Atom> atoms2;

    atoms2.push_back(NS::Atom( 0.0, 0.0, 0.0, 2.0,charge,dielectric));
    atoms2.push_back(NS::Atom( 0.0, 0.0, 0.1, 0.0,charge,dielectric));
    atoms2.push_back(NS::Atom( 0.0, 0.1, 0.0, 0.0,charge,dielectric));
    atoms2.push_back(NS::Atom( 0.1, 0.0, 0.0, 0.0,charge,dielectric));
    atoms2.push_back(NS::Atom( 0.0, 0.0,-0.1, 0.0,charge,dielectric));
    atoms2.push_back(NS::Atom( 0.0,-0.1, 0.0, 0.0,charge,dielectric));
    atoms2.push_back(NS::Atom(-0.1, 0.0, 0.0, 0.0,charge,dielectric));

    NS::NanoShaper ns4(atoms2,surf_type2,probe,stern_layer,numberOfThreads);

    // set here a consistent grid scale
    ns4.setConfig<double>("Grid_scale", 2.0 );
    ns4.setConfig<bool>("Accurate_Triangulation",true);
    ns4.setConfig<bool>("Build_epsilon_maps",false);

    // build the grid in the new mode
    ns4.setConfig<bool>("PB_grid_mode",true);

    ns4.setConfig<bool>("Patch_Based_Algorithm",true);
    ns4.setConfig<bool>("Analytical_Ray_Vs_Torus_Intersection",true);

    // impose here the min and max of the cube
    ns4.setConfig<double>("xmin",-4.0);
    ns4.setConfig<double>("ymin",-4.0);
    ns4.setConfig<double>("zmin",-4.0);

    ns4.setConfig<double>("xmax",4.0);
    ns4.setConfig<double>("ymax",4.0);
    ns4.setConfig<double>("zmax",4.0);
    // you can put here the tollerance for the PB grid check.
    // always use before buildAnalyticalSurface
    // ns4.setGridToll(1e-1);
    ns4.buildAnalyticalSurface(rand_displ);
    // remember to set this to true in order to collect rays intersections data
    ns4.setCollectGridRays(true);

    ns4.colourGrid(&surf_volume);

    std::cout << "\n" << CLIENT_INFO << " Estimated volume " << surf_volume << " [A^3]";

    std::cout << "\n" << CLIENT_INFO << " Estimated area " << surf_area << " [A^2]";

    // retrieve intersections data from the map
    auto raysMap2 = ns4.getRaysMap();
    print_map (raysMap2);
    std::cout << "\n" << CLIENT_INFO << " Along x found " << raysMap2[X_DIR].size() << " active panels";
    std::cout << "\n" << CLIENT_INFO << " Along y found " << raysMap2[Y_DIR].size() << " active panels";
    std::cout << "\n" << CLIENT_INFO << " Along z found " << raysMap2[Z_DIR].size() << " active panels";

    // iterate through the map and print the intersections and normals
    // shows all the intersections in the x direction
    auto x_crossings2 = raysMap2[Y_DIR];
    auto it2 = x_crossings2.begin();
    while (it2 != x_crossings2.end())
    {
        auto itt = it2->second;
        std::cout << "\n" << CLIENT_INFO << " Cast ray along x from pos y " << itt.point[0] << " pos z " << itt.point[1];
        for (int i=0;i<itt.inters.size();i++)
            std::cout << "\n\t" << CLIENT_INFO << " Intersection x val " << itt.inters[i] << " normal " << itt.normals[3*i+0] << "," << itt.normals[3*i+1] << "," << itt.normals[3*i+2];
        ++it2;
    }

    ns4.triangulate(&surf_area);

    std::cout << "\n";
    return 0;
}


void print_map(const std::array<std::map<std::array<double, 2>, crossings_t, map_compare>, 3> &r)
{
  std::ofstream ray_cached_file;
  // ray_cached_file.open ("ray_cache_ns.txt");
  for (int i = 0; i < 3; ++i)
  {
    std::string filename = "ray_cache_ns_";
    std::string extension = ".txt";

    filename += std::to_string(i);
    filename += extension;
    ray_cached_file.open (filename.c_str ());

    if (ray_cached_file.is_open ())
    {
      int count = 0;
      for (auto it : r[i])
      {
        ray_cached_file << "[[" << it.first.at(0) << ", " << it.first.at(1) << "]";
        if (it.second.inters.size() > 0)
        {
          count++;
        }
        for (int i = 0; i < it.second.inters.size (); i++)
          ray_cached_file << ", " << it.second.inters[i];
        ray_cached_file << "]" << std::endl;
      }
      ray_cached_file << std::endl;
      ray_cached_file << "num intersected rays: " << count << std::endl;
    }
    ray_cached_file.close ();
  }
}
