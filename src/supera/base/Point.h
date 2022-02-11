/**
 * \file Point.h
 *
 * \ingroup base
 * 
 * \brief Class def header for a class Point
 *
 * @author kazuhiro
 */

/** \addtogroup base
    @{*/
#ifndef __SUPERA_POINT_H__
#define __SUPERA_POINT_H__

#include <iostream>
#include <cmath>

namespace supera{

  /**
     \class Point3D
     Simple 3D point struct (unit of "x", "y" and "z" are not defined here and app specific)
  */
  class Point3D {
  public:
    Point3D();
    Point3D(double xv, double yv, double zv);
    ~Point3D() {}

    Point3D(const Point3D& pt) : x(pt.x), y(pt.y), z(pt.z) {}

    double x, y, z;
    
    inline bool operator== (const Point3D& rhs) const
    { return (x == rhs.x && y == rhs.y && z == rhs.z); }
    inline bool operator!= (const Point3D& rhs) const
    { return !(rhs == (*this)); }

    inline Point3D& operator/= (const double rhs)
    { x /= rhs; y /= rhs; z /= rhs; return (*this); }
    inline Point3D& operator*= (const double rhs)
    { x *= rhs; y *= rhs; z *= rhs; return (*this); }
    inline Point3D& operator+= (const Point3D& rhs)
    { x += rhs.x; y += rhs.y; z += rhs.z; return (*this); }
    inline Point3D& operator-= (const Point3D& rhs)
    { x -= rhs.x; y -= rhs.y; z -= rhs.z; return (*this); }

    inline Point3D operator/ (const double rhs) const
    { return Point3D(x/rhs,y/rhs,z/rhs); }
    inline Point3D operator* (const double rhs) const
    { return Point3D(x*rhs,y*rhs,z*rhs); }
    inline Point3D operator+ (const Point3D& rhs) const
    { return Point3D(x+rhs.x,y+rhs.y,z+rhs.z); }
    inline Point3D operator- (const Point3D& rhs) const
    { return Point3D(x-rhs.x,y-rhs.y,z-rhs.z); }

    inline double squared_distance(const Point3D& pt) const
    { return pow(x-pt.x,2)+pow(y-pt.y,2)+pow(z-pt.z,2); }
    inline double distance(const Point3D& pt) const
    { return sqrt(squared_distance(pt)); }
    inline Point3D direction(const Point3D& pt) const
    { Point3D res(pt.x - x, pt.y - y, pt.z - z); res /= distance(pt); return res; }

  };

  class Vertex {
  public:
    /// Particle ID default constructor
    Vertex();
    Vertex(double x, double y, double z, double t);

    /// Reset function
    void reset();

    /// Default destructor
    virtual ~Vertex(){};

    inline bool operator== (const Vertex& rhs) const
    {return ( pos == rhs.pos && time == rhs.time );}

    inline bool operator!= (const Vertex& rhs) const
    {return !((*this) == rhs);}

    inline bool operator< (const Vertex& rhs) const
    {
      if( pos.x     < rhs.pos.x ) return true;
      if( rhs.pos.x < pos.x     ) return false;
      if( pos.y     < rhs.pos.y ) return true;
      if( rhs.pos.y < pos.y     ) return false;
      if( pos.z     < rhs.pos.z ) return true;
      if( rhs.pos.z < pos.z     ) return false;
      if( time      < rhs.time  ) return true;
      if( rhs.time  < time      ) return false;

      return false;
    }
    void approx(unsigned int power=6);
    std::string dump() const;

    Point3D pos;
    double time;
  };
}
#endif