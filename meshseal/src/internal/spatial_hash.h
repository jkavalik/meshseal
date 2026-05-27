#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace meshseal::internal {

// Maps quantized (ix, iy, iz) cell coordinates to a list of vertex indices.
// Used for O(1) average-case nearest-vertex lookup during welding.
//
// The cell coords are stored in a struct key (CellKey) so cells at indices
// further apart than 2^21 do NOT alias — the previous 21-bit-per-axis
// truncation wrapped at ~2M cells/axis, causing spurious bucket collisions
// on extreme-aspect-ratio inputs (lego/untitled2.stl: 34km long, 5M cells
// along X). Callers metric-filter so the aliasing was correctness-safe
// but added perf overhead on large meshes.
class SpatialHash {
public:
    // cell_size: the quantization step (= weld_tolerance)
    explicit SpatialHash(double cell_size) : cell_size_(cell_size) {}

    // Insert vertex at position p with original index idx.
    void insert(const std::array<double, 3>& p, uint32_t idx) {
        cells_[cell_of(p)].push_back(idx);
    }

    // Return all vertex indices in cells within radius r of p (3×3×3 = 27-cell neighborhood search).
    std::vector<uint32_t> query_neighbors(const std::array<double, 3>& p) const {
        const CellKey c = cell_of(p);
        std::vector<uint32_t> result;
        for (int64_t dx = -1; dx <= 1; ++dx) {
            for (int64_t dy = -1; dy <= 1; ++dy) {
                for (int64_t dz = -1; dz <= 1; ++dz) {
                    const CellKey q{c.x + dx, c.y + dy, c.z + dz};
                    const auto it = cells_.find(q);
                    if (it != cells_.end()) {
                        result.insert(result.end(), it->second.begin(), it->second.end());
                    }
                }
            }
        }
        return result;
    }

private:
    struct CellKey {
        int64_t x, y, z;
        bool operator==(const CellKey& o) const noexcept {
            return x == o.x && y == o.y && z == o.z;
        }
    };
    struct CellKeyHash {
        std::size_t operator()(const CellKey& k) const noexcept {
            // boost::hash_combine-style mix of three int64s. The unordered_map
            // resolves collisions via chaining; the hash quality just needs
            // to be reasonable, not perfect.
            std::size_t h = std::hash<int64_t>{}(k.x);
            h ^= std::hash<int64_t>{}(k.y) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            h ^= std::hash<int64_t>{}(k.z) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            return h;
        }
    };

    CellKey cell_of(const std::array<double, 3>& p) const {
        return CellKey{
            static_cast<int64_t>(std::floor(p[0] / cell_size_)),
            static_cast<int64_t>(std::floor(p[1] / cell_size_)),
            static_cast<int64_t>(std::floor(p[2] / cell_size_))
        };
    }

    double cell_size_;
    std::unordered_map<CellKey, std::vector<uint32_t>, CellKeyHash> cells_;
};

} // namespace meshseal::internal
