#include "voxel_grid.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>

namespace meshseal::internal {

VoxelGrid::VoxelGrid(const std::array<double, 3>& lo,
                     const std::array<double, 3>& hi,
                     std::uint32_t target) {
    if (target < 4) target = 4;

    const double ext_x = hi[0] - lo[0];
    const double ext_y = hi[1] - lo[1];
    const double ext_z = hi[2] - lo[2];
    double maxext = std::max({ext_x, ext_y, ext_z});
    if (maxext <= 0.0) maxext = 1.0; // degenerate input — pretend unit
    cell_ = maxext / static_cast<double>(target);
    if (cell_ <= 0.0) cell_ = 1e-9;

    // 2-cell padding so the entire grid boundary is exterior.
    // Guard each axis: if hi < lo (degenerate / NaN-derived), substitute
    // a unit extent. Without this, sz can be negative and `ceil(neg/cell)`
    // casts to a huge uint32_t → enormous allocation.
    for (int k = 0; k < 3; ++k) {
        const double ext_k = (hi[k] > lo[k]) ? (hi[k] - lo[k]) : 0.0;
        lo_[k] = lo[k] - 2.0 * cell_;
        const double sz = ext_k + 4.0 * cell_;
        std::uint32_t d = std::max<std::uint32_t>(
            4, static_cast<std::uint32_t>(std::ceil(sz / cell_)));
        // Per-axis safety cap. 4096 cells/axis is more than any practical
        // mesh needs; combined with the product cap below this prevents
        // runaway allocation on pathological aspect ratios.
        if (d > 4096u) d = 4096u;
        dim_[k] = d;
    }
    // Product cap: a degenerate aspect ratio could still combine three
    // 4096s into 68 GB of Label-bytes. Cap the product at 1 GiB cells.
    // Falls back to a uniform coarser grid by halving the largest dim
    // until under cap.
    constexpr std::size_t kMaxCells = 1ull << 30; // 1 GiB cells
    std::size_t total = static_cast<std::size_t>(dim_[0]) * dim_[1] * dim_[2];
    while (total > kMaxCells) {
        int kmax = 0;
        if (dim_[1] > dim_[kmax]) kmax = 1;
        if (dim_[2] > dim_[kmax]) kmax = 2;
        if (dim_[kmax] <= 8u) break; // can't go below floor
        dim_[kmax] /= 2u;
        cell_ *= 2.0;
        total = static_cast<std::size_t>(dim_[0]) * dim_[1] * dim_[2];
    }
    labels_.assign(total, Label::Inside);
}

void VoxelGrid::rasterize(const std::array<double, 3>& a,
                          const std::array<double, 3>& b,
                          const std::array<double, 3>& c) {
    double tlo[3], thi[3];
    for (int k = 0; k < 3; ++k) {
        tlo[k] = std::min({a[k], b[k], c[k]});
        thi[k] = std::max({a[k], b[k], c[k]});
    }
    int ilo[3], ihi[3];
    for (int k = 0; k < 3; ++k) {
        const int lo_idx = static_cast<int>(std::floor((tlo[k] - lo_[k]) / cell_));
        const int hi_idx = static_cast<int>(std::floor((thi[k] - lo_[k]) / cell_));
        ilo[k] = std::max(0, lo_idx);
        ihi[k] = std::min(static_cast<int>(dim_[k]) - 1, hi_idx);
    }
    for (int i = ilo[0]; i <= ihi[0]; ++i)
    for (int j = ilo[1]; j <= ihi[1]; ++j)
    for (int k = ilo[2]; k <= ihi[2]; ++k) {
        labels_[idx(static_cast<std::uint32_t>(i),
                    static_cast<std::uint32_t>(j),
                    static_cast<std::uint32_t>(k))] = Label::Surface;
    }
}

void VoxelGrid::flood_fill_exterior() {
    std::queue<std::array<int, 3>> q;
    auto try_seed = [&](std::uint32_t i, std::uint32_t j, std::uint32_t k) {
        if (labels_[idx(i, j, k)] == Label::Inside) {
            labels_[idx(i, j, k)] = Label::Outside;
            q.push({static_cast<int>(i), static_cast<int>(j), static_cast<int>(k)});
        }
    };

    // Seed from all six faces of the grid boundary.
    for (std::uint32_t j = 0; j < dim_[1]; ++j)
    for (std::uint32_t k = 0; k < dim_[2]; ++k) {
        try_seed(0, j, k);
        try_seed(dim_[0] - 1, j, k);
    }
    for (std::uint32_t i = 0; i < dim_[0]; ++i)
    for (std::uint32_t k = 0; k < dim_[2]; ++k) {
        try_seed(i, 0, k);
        try_seed(i, dim_[1] - 1, k);
    }
    for (std::uint32_t i = 0; i < dim_[0]; ++i)
    for (std::uint32_t j = 0; j < dim_[1]; ++j) {
        try_seed(i, j, 0);
        try_seed(i, j, dim_[2] - 1);
    }

    static const int neigh[6][3] = {
        { 1, 0, 0}, {-1, 0, 0},
        { 0, 1, 0}, { 0,-1, 0},
        { 0, 0, 1}, { 0, 0,-1}
    };
    while (!q.empty()) {
        const auto p = q.front(); q.pop();
        for (auto& n : neigh) {
            const int ni = p[0] + n[0];
            const int nj = p[1] + n[1];
            const int nk = p[2] + n[2];
            if (ni < 0 || nj < 0 || nk < 0 ||
                ni >= static_cast<int>(dim_[0]) ||
                nj >= static_cast<int>(dim_[1]) ||
                nk >= static_cast<int>(dim_[2])) continue;
            auto& lab = labels_[idx(static_cast<std::uint32_t>(ni),
                                    static_cast<std::uint32_t>(nj),
                                    static_cast<std::uint32_t>(nk))];
            if (lab != Label::Inside) continue;
            lab = Label::Outside;
            q.push({ni, nj, nk});
        }
    }
}

VoxelGrid::Label VoxelGrid::sample(const std::array<double, 3>& p) const {
    const int i = static_cast<int>(std::floor((p[0] - lo_[0]) / cell_));
    const int j = static_cast<int>(std::floor((p[1] - lo_[1]) / cell_));
    const int k = static_cast<int>(std::floor((p[2] - lo_[2]) / cell_));
    if (i < 0 || j < 0 || k < 0 ||
        i >= static_cast<int>(dim_[0]) ||
        j >= static_cast<int>(dim_[1]) ||
        k >= static_cast<int>(dim_[2])) {
        return Label::Outside;
    }
    return labels_[idx(static_cast<std::uint32_t>(i),
                       static_cast<std::uint32_t>(j),
                       static_cast<std::uint32_t>(k))];
}

} // namespace meshseal::internal
