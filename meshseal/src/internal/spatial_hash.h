#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace meshseal::internal {

// Maps quantized (ix, iy, iz) cell coordinates to a list of vertex indices.
// Used for O(1) average-case nearest-vertex lookup during welding.
class SpatialHash {
public:
    // cell_size: the quantization step (= weld_tolerance)
    explicit SpatialHash(double cell_size) : cell_size_(cell_size) {}

    // Insert vertex at position p with original index idx.
    void insert(const std::array<double, 3>& p, uint32_t idx) {
        const int64_t ix = static_cast<int64_t>(std::floor(p[0] / cell_size_));
        const int64_t iy = static_cast<int64_t>(std::floor(p[1] / cell_size_));
        const int64_t iz = static_cast<int64_t>(std::floor(p[2] / cell_size_));
        cells_[encode(ix, iy, iz)].push_back(idx);
    }

    // Return all vertex indices in cells within radius r of p (3×3×3 = 27-cell neighborhood search).
    std::vector<uint32_t> query_neighbors(const std::array<double, 3>& p) const {
        const int64_t cx = static_cast<int64_t>(std::floor(p[0] / cell_size_));
        const int64_t cy = static_cast<int64_t>(std::floor(p[1] / cell_size_));
        const int64_t cz = static_cast<int64_t>(std::floor(p[2] / cell_size_));

        std::vector<uint32_t> result;
        for (int64_t dx = -1; dx <= 1; ++dx) {
            for (int64_t dy = -1; dy <= 1; ++dy) {
                for (int64_t dz = -1; dz <= 1; ++dz) {
                    const auto it = cells_.find(encode(cx + dx, cy + dy, cz + dz));
                    if (it != cells_.end()) {
                        result.insert(result.end(), it->second.begin(), it->second.end());
                    }
                }
            }
        }
        return result;
    }

private:
    double cell_size_;
    std::unordered_map<uint64_t, std::vector<uint32_t>> cells_;

    static uint64_t encode(int64_t ix, int64_t iy, int64_t iz) {
        return static_cast<uint64_t>(ix & 0x1FFFFF)
             | (static_cast<uint64_t>(iy & 0x1FFFFF) << 21)
             | (static_cast<uint64_t>(iz & 0x1FFFFF) << 42);
    }
};

} // namespace meshseal::internal
