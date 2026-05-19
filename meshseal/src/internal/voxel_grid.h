#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace meshseal::internal {

// 3D occupancy grid used as a coarse inside/outside oracle for the
// volumetric-fallback / soup-reconstruction stages. Three jobs:
//
//   1. Rasterize input triangles → mark intersected cells as Surface.
//   2. Flood-fill 6-connected from boundary cells → mark reached cells as
//      Outside. Unreached non-Surface cells remain Inside (interior pockets
//      enclosed by the input geometry).
//   3. Sample(p) returns the label at world point p so callers can decide
//      whether a triangle / vertex / fragment lies on the boundary of the
//      inferred volume.
//
// The grid does NOT generate output geometry — it's purely an oracle.
// Output triangles come from the original input, classified by sampling
// this grid (typically at small offsets above and below their centroids).
//
// Design choices:
//   - First-pass rasterization uses the triangle's AABB, not exact triangle-
//     AABB intersection (Akenine-Möller SAT). This is conservative — it can
//     over-mark cells along thin slivers, but never misses a triangle. A SAT
//     refinement is a future improvement; AABB is sufficient at the
//     resolutions we use.
//   - A 2-cell padding ring around the input bbox guarantees that flood-fill
//     seeds (the entire grid boundary) are exterior even when input triangles
//     touch the bbox of the input.
class VoxelGrid {
public:
    enum class Label : std::uint8_t {
        Inside  = 0,  // not reached by flood fill, no surface — interior pocket
        Outside = 1,  // reached by 6-connected flood fill from grid boundary
        Surface = 2,  // contains at least one input triangle's AABB
    };

    // Build empty grid spanning [lo, hi] padded by 2 cells on each side.
    // cell size = max-extent / target_cells_per_long_axis.
    VoxelGrid(const std::array<double, 3>& lo,
              const std::array<double, 3>& hi,
              std::uint32_t target_cells_per_long_axis);

    // Mark all cells overlapping the triangle's AABB as Surface.
    void rasterize(const std::array<double, 3>& a,
                   const std::array<double, 3>& b,
                   const std::array<double, 3>& c);

    // 6-connected flood fill seeded from every cell on the grid's outer
    // boundary. Cells that are still Inside afterwards are interior pockets.
    // Call AFTER all rasterize() calls.
    void flood_fill_exterior();

    // Sample the cell at world point p. Out-of-bounds samples return Outside
    // (the grid is padded so this only happens far from the input).
    Label sample(const std::array<double, 3>& p) const;

    // Cell side length. Useful for choosing sample-offset distances.
    double cell_size() const { return cell_; }

    // Grid resolution. (mostly for diagnostics)
    std::array<std::uint32_t, 3> dims() const { return dim_; }

private:
    std::array<double, 3>        lo_;     // padded grid origin (world coords)
    std::array<std::uint32_t, 3> dim_;
    double                       cell_;
    std::vector<Label>           labels_;

    inline std::size_t idx(std::uint32_t i, std::uint32_t j, std::uint32_t k) const {
        return (static_cast<std::size_t>(k) * dim_[1] + j) * dim_[0] + i;
    }
};

} // namespace meshseal::internal
