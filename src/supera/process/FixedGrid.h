#ifndef SUPERA_PROCESS_FIXEDGRID_H
#define SUPERA_PROCESS_FIXEDGRID_H

#include "supera/data/Grid3D.h"
#include "supera/data/ImageMeta3D.h"

namespace supera {

/// Checked v2 boundary around Atomic's existing voxel indexer.
/// Construction validates metadata once; queries never mutate the grid.
class FixedGrid {
public:
  explicit FixedGrid(Grid3D metadata);

  /// Truncate requested length/width to counts, then derive effective widths.
  /// This preserves the legacy fixed-box convention without depending on FHiCL.
  static FixedGrid from_lengths(const std::array<double, 3> &minimum_cm,
                                const std::array<double, 3> &length_cm,
                                const std::array<double, 3> &requested_width_cm);

  const Grid3D &metadata() const {
    return metadata_;
  }

  RecordIndex size() const {
    return size_;
  }

  std::array<double, 3> widths_cm() const;

  /// Finite points outside the volume return missing_index. Nonfinite input fails.
  RecordIndex id(const std::array<double, 3> &position_cm) const;

  /// Discrete coordinates and origins must be valid, not missing or out of range.
  RecordIndex index(const std::array<RecordIndex, 3> &coordinates) const;
  std::array<double, 3> position(RecordIndex id) const;

  /// A valid origin shifted outside the volume returns missing_index.
  RecordIndex shift(RecordIndex origin, const std::array<int, 3> &displacement) const;

private:
  Grid3D metadata_;
  ImageMeta3D indexer_;
  RecordIndex size_{};
};

} // namespace supera
#endif
