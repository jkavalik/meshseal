#include "winding_number.h"
#include <cmath>

namespace meshseal::internal {

double tri_solid_angle(const std::array<double, 3>& a,
                       const std::array<double, 3>& b,
                       const std::array<double, 3>& c,
                       const std::array<double, 3>& p)
{
    const double ax = a[0]-p[0], ay = a[1]-p[1], az = a[2]-p[2];
    const double bx = b[0]-p[0], by = b[1]-p[1], bz = b[2]-p[2];
    const double cx = c[0]-p[0], cy = c[1]-p[1], cz = c[2]-p[2];

    const double la = std::sqrt(ax*ax + ay*ay + az*az);
    const double lb = std::sqrt(bx*bx + by*by + bz*bz);
    const double lc = std::sqrt(cx*cx + cy*cy + cz*cz);

    // Degenerate: query point coincides with a vertex.
    if (la < 1e-30 || lb < 1e-30 || lc < 1e-30) return 0.0;

    // Numerator: det([a-p, b-p, c-p]) = a·(b×c).
    const double det =
        ax * (by*cz - bz*cy) -
        ay * (bx*cz - bz*cx) +
        az * (bx*cy - by*cx);

    // Denominator: product term as per the robust atan2 form.
    const double dot_ab = ax*bx + ay*by + az*bz;
    const double dot_bc = bx*cx + by*cy + bz*cz;
    const double dot_ca = cx*ax + cy*ay + cz*az;
    const double denom = la*lb*lc + dot_ab*lc + dot_bc*la + dot_ca*lb;

    return 2.0 * std::atan2(det, denom);
}

double generalized_winding_number(const Mesh& mesh,
                                  const std::array<double, 3>& p)
{
    constexpr double kInv4Pi = 1.0 / (4.0 * 3.14159265358979323846);
    const size_t nv = mesh.vertices.size();
    double sum = 0.0;
    for (const auto& f : mesh.faces) {
        // Defensive: skip any face whose indices are out of range. Public
        // repair() Phase 0 validates this for normal callers, but this
        // function is also reachable via `orient_wn` which calls it O(F²)
        // times — one tainted face would silently UB across every query.
        if (f[0] >= nv || f[1] >= nv || f[2] >= nv) continue;
        sum += tri_solid_angle(mesh.vertices[f[0]],
                               mesh.vertices[f[1]],
                               mesh.vertices[f[2]], p);
    }
    return sum * kInv4Pi;
}

} // namespace meshseal::internal
