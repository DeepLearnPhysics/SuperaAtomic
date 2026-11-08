# Fixed-grid contract

Status: validated grid foundation for the reconstructed space-point checkpoint.
The coordinated adapter now includes extraction, shared collision selection, and
LArCV2 projection; this document specifies the grid contract.

## Data and computation

`supera/data/Grid3D.h` is a standard-library-only record of minimum/maximum bounds
in centimeters and integer voxel counts. It is available through `supera::data`.
`supera/process/FixedGrid.h` provides validated construction and queries through
`supera::labeling`. It wraps the existing `ImageMeta3D` implementation; it does not
introduce a second coordinate-to-voxel formula or a LArCV dependency. The old
unsigned invalid-voxel sentinel is translated only at this explicit bridge.

The grid is Cartesian and x-fastest. Both outer faces are included, with an exact
upper-bound coordinate assigned to the final voxel. This deliberately preserves
the existing fixed-box behavior. It is not a detector geometry or an active-TPC
acceptance mask. Callers must supply points in the same coordinate frame as the
bounds; this step performs no frame transformation.

`Grid3D::counts` is authoritative. `FixedGrid::from_lengths` accepts minimum bounds,
lengths and requested widths, truncates length/width to integer counts, then derives
effective widths from the actual bounds/counts. Requested and effective widths
can differ. The ICARUS configuration gives 6144 bins per axis and 231,928,233,984
voxels; no voxel ID is narrowed to 32 bits.

## Validity and failure behavior

Construction rejects nonfinite bounds, nonpositive/overflowing lengths, zero
counts, overflowing count products, and nonfinite or underflowed effective widths.
Per-axis counts are limited to 2^53 - 1 for the existing double-based indexer. The
product must fit both the canonical signed 64-bit range and the underlying native
index types. Requested-width quotients are checked before integer conversion.
The record defaults to unspecified bounds and zero counts, not a usable grid.

| Operation | Input/result contract |
| --- | --- |
| `id(position_cm)` | Requires three finite coordinates; a finite out-of-volume point returns `missing_index` |
| `index(coordinates)` | Requires three nonnegative, in-range discrete indices; malformed indices throw |
| `position(id)` | Requires an existing voxel ID; returns its center in centimeters |
| `shift(origin, displacement)` | Requires a valid origin; an out-of-volume neighbor returns `missing_index` |
| `widths_cm()` | Returns effective widths derived by the established indexer |

A missing measurement cannot silently become an ordinary geometric rejection.
Invalid metadata throws `invalid_argument`; invalid required discrete indices
throw `out_of_range`. Shift offsets currently use the legacy signed `int` domain,
including safe handling of INT_MIN. A constructed grid cannot be reconfigured;
queries are const and validate before calling legacy operations.

These guards protect this new entry point. Direct legacy `ImageMeta3D` use has not
been globally rewritten or declared safe for arbitrary malformed metadata.

## Validation and reproducibility

- `FixedGrid_test.cxx` exercises non-cubic flattening, independently computed voxel
  centers, all outer boundaries and their neighboring doubles, internal boundaries,
  truncation, large IDs, invalid positions/indices, count overflow and underflow.
- All six shared/labeling tests pass with address, undefined-behavior and
  float-to-integer overflow sanitizers.
- `FixedGridReference_test.cxx` uses the actual LArCV2 implementation. Its optional
  filename arguments inspect the twenty ICARUS SpacePoint product trees, checking
  every event's bounds, counts, widths and size exactly. It samples stored voxel
  centers across each product's ID range and includes independent boundary probes.
- Against both immutable references, the local native run checked 400 product-events
  and 13,181 coordinate probes with no differences. This establishes metadata and
  sampled grid-query agreement, not a new twenty-product rasterization result.
- An installed-package consumer checks `find_package(Supera COMPONENTS processing)`
  and the exported target without ROOT or LArCV present. The CI workflow includes
  this consumer; data-only builds also compile the new metadata record.

When both processing and LArCV2 are enabled, CMake builds
`fixed_grid_reference_test`. Without filenames, its CTest invocation performs
synthetic installed-LArCV boundary checks. Supply the staged and standard reference
paths to reproduce the full fixture metadata probe. In the pinned ICARUS lab, the
same actual source files were also compiled directly with ROOT/LArCV2 for this
probe; commands and logs live under `research/fixed-grid/` and
`results/fixed-grid-reference.log`.

The subsequent SpacePoint/Hit extraction, shared representative selection, and
LArCV2 projection are now implemented. The coordinated adapter output review
records exact agreement for every field of all twenty configured products.
