#ifndef __IMAGEMETA3D_CXX__
#define __IMAGEMETA3D_CXX__

#include "ImageMeta3D.h"
#include "supera/base/meatloaf.h"
#include <algorithm>
#include <cstdint>
#include <sstream>

namespace supera {

  ImageMeta3D::ImageMeta3D()
  { clear(); }

  void ImageMeta3D::clear()
  {
    (BBox3D)(*this) = BBox3D();
    _valid = false;
    _xnum = _ynum = _znum = 0;
    _xlen = _ylen = _zlen = -1.;
  }

  void ImageMeta3D::update(size_t xnum, size_t ynum, size_t znum) {

    if(empty())
      throw meatloaf("Empty voxel volume definition cannot be updated with new voxel counts!");

    if(xnum == kINVALID_SIZE || xnum == 0)
      throw meatloaf("ImageMeta3D::Set x voxel count not set!");

    if(ynum == kINVALID_SIZE || ynum == 0)
      throw meatloaf("ImageMeta3D::Set y voxel count not set!");

    if(znum == kINVALID_SIZE || znum == 0)
      throw meatloaf("ImageMeta3D::Set z voxel count not set!");

    _xlen = (max_x() - min_x()) / ((double)xnum);
    _xnum = xnum;

    _ylen = (max_y() - min_y()) / ((double)ynum);
    _ynum = ynum;

    _zlen = (max_z() - min_z()) / ((double)znum);
    _znum = znum;

    _num_element = _xnum * _ynum * _znum;
    _valid = true;
  }

  VoxelID_t ImageMeta3D::id(double x, double y, double z) const
  {
    if(!_valid) throw meatloaf("ImageMeta3D::ID cannot be called on invalid meta!");
    if(x > max_x() || x < min_x()) return kINVALID_VOXELID;
    if(y > max_y() || y < min_y()) return kINVALID_VOXELID;
    if(z > max_z() || z < min_z()) return kINVALID_VOXELID;

    VoxelID_t xindex = (x - min_x()) / _xlen;
    VoxelID_t yindex = (y - min_y()) / _ylen;
    VoxelID_t zindex = (z - min_z()) / _zlen;

    if(xindex == _xnum) xindex -= 1;
    if(yindex == _ynum) yindex -= 1;
    if(zindex == _znum) zindex -= 1;

    return (zindex * (_xnum * _ynum) + yindex * _xnum + xindex);
  }
  VoxelID_t ImageMeta3D::index(const size_t i_x, const size_t i_y, const size_t i_z) const
  {
    if (!_valid) throw meatloaf("ImageMeta3D::ID cannot be called on invalid meta!");
    if (i_x >= _xnum) return kINVALID_VOXELID;
    if (i_y >= _ynum) return kINVALID_VOXELID;
    if (i_z >= _znum) return kINVALID_VOXELID;

    return (i_z * (_xnum * _ynum) + i_y * _xnum + i_x);

  }

  VoxelID_t ImageMeta3D::shift(const VoxelID_t origin_id,
                             const int shift_x,
                             const int shift_y,
                             const int shift_z) const
  {
    // Validate the grid and origin before decoding unsigned voxel coordinates.
    if (!_valid)
      throw meatloaf("ImageMeta3D::shift cannot be called on invalid meta!");
    if (origin_id >= _num_element)
      return kINVALID_VOXELID;

    size_t x, y, z;
    id_to_xyz_index(origin_id, x, y, z);

    // Check each signed displacement before unsigned arithmetic. Widen first so
    // negating INT_MIN is defined, and never narrow the flattened voxel ID.
    auto shift_axis = [](size_t &index, size_t count, int displacement) {
      const auto delta = static_cast<std::int64_t>(displacement);
      if (delta < 0) {
        const auto distance = static_cast<std::uint64_t>(-delta);
        if (distance > index)
          return false;
        index -= distance;
      } else {
        const auto distance = static_cast<std::uint64_t>(delta);
        if (distance >= count - index)
          return false;
        index += distance;
      }
      return true;
    };

    // A failed axis shift leaves the caller's grid unchanged and yields no voxel.
    if (!shift_axis(x, _xnum, shift_x) || !shift_axis(y, _ynum, shift_y) ||
        !shift_axis(z, _znum, shift_z))
      return kINVALID_VOXELID;

    return index(x, y, z);
  }

  Point3D ImageMeta3D::position(VoxelID_t id) const
  {
    if(!_valid) throw meatloaf("ImageMeta3D::pos cannot be called on invalid meta!");
    if(id >= _num_element) throw meatloaf("ImageMeta3D::pos invalid VoxelID_t!");

    VoxelID_t zid = id / (_xnum * _ynum);
    id -= zid * (_xnum * _ynum);
    VoxelID_t yid = id / _xnum;
    VoxelID_t xid = (id - yid * _xnum);

    return Point3D(min_x() + ((double)xid + 0.5) * _xlen,
		   min_y() + ((double)yid + 0.5) * _ylen,
		   min_z() + ((double)zid + 0.5) * _zlen);
  }

  double ImageMeta3D::pos_x(VoxelID_t id) const
  {
    if(!_valid) throw meatloaf("ImageMeta3D::pos_x cannot be called on invalid meta!");
    if(id >= _num_element) throw meatloaf("ImageMeta3D::pos_x invalid VoxelID_t!");

    VoxelID_t zid = id / (_xnum * _ynum);
    id -= zid * (_xnum * _ynum);
    VoxelID_t yid = id / _xnum;
    VoxelID_t xid = (id - yid * _xnum);

    return min_x() + ((double)xid + 0.5) * _xlen;
  }

  double ImageMeta3D::pos_y(VoxelID_t id) const
  {
    if(!_valid) throw meatloaf("ImageMeta3D::pos_y cannot be called on invalid meta!");
    if(id >= _num_element) throw meatloaf("ImageMeta3D::pos_y invalid VoxelID_t!");

    VoxelID_t zid = id / (_xnum * _ynum);
    id -= zid * (_xnum * _ynum);
    VoxelID_t yid = id / _xnum;
    return min_y() + ((double)yid + 0.5) * _ylen;
  }

  double ImageMeta3D::pos_z(VoxelID_t id) const
  {
    if(!_valid) throw meatloaf("ImageMeta3D::pos_z cannot be called on invalid meta!");
    if(id >= _num_element) throw meatloaf("ImageMeta3D::pos_z invalid VoxelID_t!");

    VoxelID_t zid = id / (_xnum * _ynum);
    return min_z() + ((double)zid + 0.5) * _zlen;
  }

  size_t ImageMeta3D::id_to_x_index(VoxelID_t id) const
  {
    if(id >= _num_element) throw meatloaf("ImageMeta3D::pos invalid VoxelID_t!");

    VoxelID_t zid = id / (_xnum * _ynum);
    id -= zid * (_xnum * _ynum);
    VoxelID_t yid = id / _xnum;
    VoxelID_t xid = (id - yid * _xnum);

    return xid;
  }

  size_t ImageMeta3D::id_to_y_index(VoxelID_t id) const
  {
    if(id >= _num_element) throw meatloaf("ImageMeta3D::pos invalid VoxelID_t!");

    VoxelID_t zid = id / (_xnum * _ynum);
    id -= zid * (_xnum * _ynum);
    VoxelID_t yid = id / _xnum;

    return yid;
  }


  size_t ImageMeta3D::id_to_z_index(VoxelID_t id) const
  {
    if(id >= _num_element) throw meatloaf("ImageMeta3D::pos invalid VoxelID_t!");

    VoxelID_t zid = id / (_xnum * _ynum);

    return zid;
  }

  void ImageMeta3D::id_to_xyz_index(VoxelID_t id, size_t& x, size_t& y, size_t& z) const
  {
    if(id >= _num_element) throw meatloaf("ImageMeta3D::pos invalid VoxelID_t!");
    z = id / (_xnum * _ynum);
    id -= z * (_xnum * _ynum);
    y = id / _xnum;
    x = (id - y * _xnum);
  }

  std::string  ImageMeta3D::dump() const
  {
    std::stringstream ss;
    ss << "X range: " << min_x() << " => " << max_x() << " ... " << _xnum << " bins" << std::endl
       << "Y range: " << min_y() << " => " << max_y() << " ... " << _ynum << " bins" << std::endl
       << "Z range: " << min_z() << " => " << max_z() << " ... " << _znum << " bins" << std::endl;
    return std::string(ss.str());
  }

  std::string ImageMeta3D::dump2cpp(const std::string &instanceName) const
  {
     std::stringstream ss;

     ss << "supera::ImageMeta3D " << instanceName << ";\n";
     ss << instanceName << ".set("
        << min_x() << ", " << min_y() << ", " << min_z() << ", "
        << max_x() << ", " << max_y() << ", " << max_z() << ", "
        << num_voxel_x() << ", " << num_voxel_y() << ", " << num_voxel_z()
        << ");\n";

    return ss.str();

  }


  VoxelSet ImageMeta3D::edep2voxelset(const std::vector<supera::EDep> edeps) const
  {
    VoxelSet result;
    result.reserve(edeps.size());
    for(auto const& edep : edeps){
      auto vox_id = this->id(edep.x,edep.y,edep.z);
      if(vox_id == supera::kINVALID_VOXELID)
        continue;
      result.emplace(vox_id, edep.e, true);
    }
    return result;
  }

};

#endif
