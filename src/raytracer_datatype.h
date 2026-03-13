#ifndef RAYTRACER_DATATYPE_H
#define RAYTRACER_DATATYPE_H

#include <map>
#include <vector>
#include <array>


struct crossings_t
{
  static double start[3], end[3];

  double point[2]; // point x and z coords: the ones that prescribe the ray
  unsigned dir;
  bool init = 0; // 1 when the ray is initilized

  std::vector<double> inters; // intersections
  std::vector<double> normals; // normals
};


struct map_compare
{
  bool operator ()(const std::array<double,2> &arr1, const std::array<double,2> &arr2) const
  {
    double tol = 1.e-3;

    if ((arr1.at(0) - arr2.at(0)) < -tol)
      return 1;
    else if (std::fabs(arr1.at(0) - arr2.at(0)) <= tol)
    {
      if ((arr1.at(1) - arr2.at(1)) < -tol)
        return 1;
      else
        return 0;
    }
    else
      return 0;
  }
};


using rays_t = std::array<std::map<std::array<double, 2>, crossings_t, map_compare>, 3>;

#endif
