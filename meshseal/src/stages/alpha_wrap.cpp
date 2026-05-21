#include "alpha_wrap.h"
#include "../internal/manifold_weld.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>
#include <manifold/common.h>
#include <manifold/manifold.h>

namespace meshseal::stages {

namespace {

using Vec3 = std::array<double, 3>;
constexpr float kBig = 1e18f;

inline double vd2(const Vec3& a, const Vec3& b) {
    const double x=a[0]-b[0], y=a[1]-b[1], z=a[2]-b[2];
    return x*x + y*y + z*z;
}

// 1D squared-distance transform (Felzenszwalb & Huttenlocher 2012).
// On entry f[q] = 0 at seeds, large elsewhere. On exit f[q] = min_p (f[p] +
// (q-p)^2) — the squared distance to the nearest seed along the line.
void edt1d(std::vector<double>& f) {
    const int n = static_cast<int>(f.size());
    if (n == 0) return;
    std::vector<int>    v(n);
    std::vector<double> z(n + 1);
    std::vector<double> d(n);
    int k = 0;
    v[0] = 0;
    z[0] = -1e30;
    z[1] =  1e30;
    for (int q = 1; q < n; ++q) {
        double s;
        for (;;) {
            s = ((f[q] + double(q) * q) - (f[v[k]] + double(v[k]) * v[k]))
                / (2.0 * q - 2.0 * v[k]);
            if (k > 0 && s <= z[k]) --k;
            else break;
        }
        ++k;
        v[k]     = q;
        z[k]     = s;
        z[k + 1] = 1e30;
    }
    k = 0;
    for (int q = 0; q < n; ++q) {
        while (z[k + 1] < q) ++k;
        const double dq = q - v[k];
        d[q] = dq * dq + f[v[k]];
    }
    f.swap(d);
}

// 3D Euclidean distance transform by separable 1D passes. `g` holds 0 at
// seed voxels and kBig elsewhere; on exit it holds the squared distance to
// the nearest seed, in voxel units.
void edt3d(std::vector<float>& g, int nx, int ny, int nz) {
    auto at = [&](int i, int j, int k) -> std::size_t {
        return (std::size_t(k) * ny + j) * nx + i;
    };
    std::vector<double> line;
    line.resize(nx);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) line[i] = g[at(i, j, k)];
            edt1d(line);
            for (int i = 0; i < nx; ++i) g[at(i, j, k)] = float(line[i]);
        }
    line.resize(ny);
    for (int k = 0; k < nz; ++k)
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) line[j] = g[at(i, j, k)];
            edt1d(line);
            for (int j = 0; j < ny; ++j) g[at(i, j, k)] = float(line[j]);
        }
    line.resize(nz);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            for (int k = 0; k < nz; ++k) line[k] = g[at(i, j, k)];
            edt1d(line);
            for (int k = 0; k < nz; ++k) g[at(i, j, k)] = float(line[k]);
        }
}

} // anonymous namespace

AlphaWrapResult alpha_wrap(const Mesh& mesh, double alpha, int voxel_res) {
    AlphaWrapResult result;
    if (mesh.faces.empty() || mesh.vertices.empty()) {
        result.reason = "empty input";
        return result;
    }

    // ---- bounding box ----
    Vec3 lo = mesh.vertices[0], hi = mesh.vertices[0];
    for (const auto& v : mesh.vertices)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], v[k]);
            hi[k] = std::max(hi[k], v[k]);
        }
    const double ext[3] = { hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2] };
    const double max_ext = std::max({ ext[0], ext[1], ext[2] });
    const double diag = std::sqrt(ext[0]*ext[0] + ext[1]*ext[1] + ext[2]*ext[2]);
    if (max_ext <= 0.0 || diag <= 0.0) {
        result.reason = "degenerate bbox";
        return result;
    }

    if (voxel_res < 32) voxel_res = 32;
    double cell = max_ext / double(voxel_res);

    // alpha: auto = 2% of bbox diagonal, floored so the ball spans >= 4
    // voxels (an under-resolved ball produces a faceted, unreliable wrap).
    if (alpha < 0.0) alpha = 0.02 * diag;
    const double min_alpha = 4.0 * cell;
    if (alpha < min_alpha) alpha = min_alpha;
    result.alpha_used = alpha;

    // ---- grid: pad enough that the boundary is clear exterior ----
    const int pad = int(std::ceil(alpha / cell)) + 3;
    int nx = int(std::ceil(ext[0] / cell)) + 2 * pad + 1;
    int ny = int(std::ceil(ext[1] / cell)) + 2 * pad + 1;
    int nz = int(std::ceil(ext[2] / cell)) + 2 * pad + 1;

    // cap total voxels — coarsen uniformly if the grid is too large.
    while (double(nx) * ny * nz > 18.0e6 && voxel_res > 48) {
        voxel_res = voxel_res * 3 / 4;
        cell = max_ext / double(voxel_res);
        if (alpha < 4.0 * cell) alpha = 4.0 * cell;
        result.alpha_used = alpha;
        const int p2 = int(std::ceil(alpha / cell)) + 3;
        nx = int(std::ceil(ext[0] / cell)) + 2 * p2 + 1;
        ny = int(std::ceil(ext[1] / cell)) + 2 * p2 + 1;
        nz = int(std::ceil(ext[2] / cell)) + 2 * p2 + 1;
    }
    const Vec3 glo = { lo[0] - pad * cell, lo[1] - pad * cell, lo[2] - pad * cell };
    const std::size_t ncell = std::size_t(nx) * ny * nz;
    auto idx = [&](int i, int j, int k) -> std::size_t {
        return (std::size_t(k) * ny + j) * nx + i;
    };

    // ---- seed Surface voxels: subdivide each triangle by longest-edge
    //      bisection until every sub-triangle fits inside a voxel, then mark
    //      its corners. Work is O(surface area / cell^2 + perimeter / cell)
    //      regardless of triangle shape — a model-spanning or sliver triangle
    //      cannot blow it up (unlike an AABB scan). ----
    std::vector<float> field(ncell, kBig);
    auto mark = [&](const Vec3& p) {
        const int i = int((p[0] - glo[0]) / cell);
        const int j = int((p[1] - glo[1]) / cell);
        const int k = int((p[2] - glo[2]) / cell);
        if (i >= 0 && j >= 0 && k >= 0 && i < nx && j < ny && k < nz)
            field[idx(i, j, k)] = 0.0f;
    };
    {
        const double cell2 = cell * cell;
        std::vector<std::array<Vec3, 3>> work;
        for (const auto& f : mesh.faces) {
            work.push_back({ mesh.vertices[f[0]], mesh.vertices[f[1]],
                             mesh.vertices[f[2]] });
            while (!work.empty()) {
                const std::array<Vec3, 3> t = work.back();
                work.pop_back();
                const Vec3& a = t[0]; const Vec3& b = t[1]; const Vec3& c = t[2];
                const double eab = vd2(a,b), ebc = vd2(b,c), eca = vd2(c,a);
                const double m = std::max({ eab, ebc, eca });
                if (m <= cell2) {
                    mark(a); mark(b); mark(c);
                    mark({ (a[0]+b[0]+c[0])/3, (a[1]+b[1]+c[1])/3,
                           (a[2]+b[2]+c[2])/3 });
                    continue;
                }
                if (eab >= ebc && eab >= eca) {
                    const Vec3 mid{ (a[0]+b[0])/2, (a[1]+b[1])/2, (a[2]+b[2])/2 };
                    work.push_back({ a, mid, c });
                    work.push_back({ mid, b, c });
                } else if (ebc >= eca) {
                    const Vec3 mid{ (b[0]+c[0])/2, (b[1]+c[1])/2, (b[2]+c[2])/2 };
                    work.push_back({ a, b, mid });
                    work.push_back({ a, mid, c });
                } else {
                    const Vec3 mid{ (c[0]+a[0])/2, (c[1]+a[1])/2, (c[2]+a[2])/2 };
                    work.push_back({ a, b, mid });
                    work.push_back({ mid, b, c });
                }
            }
        }
    }

    // ---- distance transform #1: dist_to_surface (squared, voxel units) ----
    edt3d(field, nx, ny, nz);
    const double alpha_vox = alpha / cell;
    const float  alpha_vox_sq = float(alpha_vox * alpha_vox);

    // ---- reachable probe centres C: flood-fill from the grid boundary
    //      through voxels the radius-alpha ball fits in (dist_to_surface
    //      squared > alpha_vox^2). ----
    std::vector<std::uint8_t> reach(ncell, 0);
    std::vector<std::size_t>  stack;
    auto seed = [&](int i, int j, int k) {
        const std::size_t s = idx(i, j, k);
        if (!reach[s] && field[s] > alpha_vox_sq) { reach[s] = 1; stack.push_back(s); }
    };
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) { seed(i,j,0); seed(i,j,nz-1); }
    for (int k = 0; k < nz; ++k) for (int i = 0; i < nx; ++i) { seed(i,0,k); seed(i,ny-1,k); }
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) { seed(0,j,k); seed(nx-1,j,k); }
    while (!stack.empty()) {
        const std::size_t s = stack.back(); stack.pop_back();
        const int i = int(s % nx);
        const int j = int((s / nx) % ny);
        const int k = int(s / (std::size_t(nx) * ny));
        const int di[6]={1,-1,0,0,0,0}, dj[6]={0,0,1,-1,0,0}, dk[6]={0,0,0,0,1,-1};
        for (int n = 0; n < 6; ++n) {
            const int ii=i+di[n], jj=j+dj[n], kk=k+dk[n];
            if (ii<0||jj<0||kk<0||ii>=nx||jj>=ny||kk>=nz) continue;
            const std::size_t t = idx(ii,jj,kk);
            if (!reach[t] && field[t] > alpha_vox_sq) { reach[t]=1; stack.push_back(t); }
        }
    }

    // ---- distance transform #2: dist_to_C, converted to WORLD units ----
    std::vector<float> distC(ncell, kBig);
    for (std::size_t s = 0; s < ncell; ++s) if (reach[s]) distC[s] = 0.0f;
    edt3d(distC, nx, ny, nz);
    for (std::size_t s = 0; s < ncell; ++s)
        distC[s] = float(std::sqrt(double(distC[s])) * cell);

    // ---- signed field  s = alpha - dist_to_C  (s < 0 inside the solid) ----
    auto sdf = [&](manifold::vec3 mp) -> double {
        double gx = (mp.x - glo[0]) / cell - 0.5;
        double gy = (mp.y - glo[1]) / cell - 0.5;
        double gz = (mp.z - glo[2]) / cell - 0.5;
        int i0 = int(std::floor(gx)), j0 = int(std::floor(gy)), k0 = int(std::floor(gz));
        double tx = gx - i0, ty = gy - j0, tz = gz - k0;
        auto clampi = [](int a, int n){ return a<0?0:(a>=n?n-1:a); };
        const int i1=clampi(i0+1,nx), j1=clampi(j0+1,ny), k1=clampi(k0+1,nz);
        i0=clampi(i0,nx); j0=clampi(j0,ny); k0=clampi(k0,nz);
        auto D=[&](int i,int j,int k){ return double(distC[idx(i,j,k)]); };
        const double c00=D(i0,j0,k0)*(1-tx)+D(i1,j0,k0)*tx;
        const double c10=D(i0,j1,k0)*(1-tx)+D(i1,j1,k0)*tx;
        const double c01=D(i0,j0,k1)*(1-tx)+D(i1,j0,k1)*tx;
        const double c11=D(i0,j1,k1)*(1-tx)+D(i1,j1,k1)*tx;
        const double c0=c00*(1-ty)+c10*ty, c1=c01*(1-ty)+c11*ty;
        const double dist_to_C = c0*(1-tz)+c1*tz;
        return dist_to_C - alpha;   // >0 inside the solid (manifold convention)
    };

    // ---- marching cubes via manifold::LevelSet ----
    const manifold::Box bounds(
        manifold::vec3(glo[0], glo[1], glo[2]),
        manifold::vec3(glo[0]+nx*cell, glo[1]+ny*cell, glo[2]+nz*cell));
    const double edge_length = cell * 2.0;
    manifold::Manifold m;
    try {
        m = manifold::Manifold::LevelSet(sdf, bounds, edge_length, /*level=*/0.0);
    } catch (const std::exception& e) {
        result.reason = std::string("LevelSet threw: ") + e.what();
        return result;
    }

    std::string fail;
    Mesh out = internal::manifold_to_welded_mesh(m, fail);
    if (out.faces.empty()) {
        result.reason = fail.empty() ? "alpha_wrap produced empty mesh" : fail;
        return result;
    }
    result.mesh    = std::move(out);
    result.success = true;
    return result;
}

} // namespace meshseal::stages
