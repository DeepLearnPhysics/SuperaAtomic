#ifndef __SUPERA_POINT_CXX__
#define __SUPERA_POINT_CXX__

#include "Point.h"
#include "SuperaType.h"
#include <sstream>
#include <math.h>
namespace supera{

  Point3D::Point3D()
    : x(kINVALID_DOUBLE), y(kINVALID_DOUBLE), z(kINVALID_DOUBLE)
  {}

  Point3D::Point3D(double xv, double yv, double zv)
    : x(xv), y(yv), z(zv)
  {}

  Vertex::Vertex() : pos(), time(kINVALID_DOUBLE)
  {}

  Vertex::Vertex(double xv, double yv, double zv, double tv)
   : pos(xv,yv,zv), time(tv)
  {}

  void Vertex::reset()
  { pos.x = pos.y = pos.z = time = kINVALID_DOUBLE;}

 std::string Vertex::dump() const
  {
    std::stringstream ss;
    ss << "(" << pos.x << ", " << pos.y << ", " << pos.z << ", " << time << ")" << std::endl;
    return ss.str();
  }

  void Vertex::approx(unsigned int power)
  {
    double factor = std::pow(10, power);
    pos.x = (double)( ((double)((signed long long)(pos.x * factor)) / factor ));
    pos.y = (double)( ((double)((signed long long)(pos.y * factor)) / factor ));
    pos.z = (double)( ((double)((signed long long)(pos.z * factor)) / factor ));
    time  = (double)( ((double)((signed long long)(time  * factor)) / factor ));
  }


}

#endif