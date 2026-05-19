#include "arrangement.h"
#include <algorithm>
#include <cmath>

namespace meshseal::internal {

namespace {

constexpr double kPlaneEps = 1e-12;
constexpr double kIntervalEps = 1e-12;

// Linear interpolation Pa→Pb at the plane crossing implied by signed
// distances dPa, dPb (assumed opposite-signed).
inline Vec3d edge_plane_intersect(const Vec3d& Pa, const Vec3d& Pb,
                                  double dPa, double dPb) {
    const double t = dPa / (dPa - dPb);
    return add(scale(Pa, 1.0 - t), scale(Pb, t));
}

// Given a triangle (V0,V1,V2) with signed distances (d0,d1,d2) to some
// reference plane, return the two points where its edges cross that plane.
// Caller must have verified that the triangle does cross the plane (not all
// same sign, not all ~zero).
bool tri_plane_cross_segment(const Vec3d& V0, const Vec3d& V1, const Vec3d& V2,
                             double d0, double d1, double d2,
                             Vec3d& p_out, Vec3d& q_out) {
    const int npos = (d0 > 0) + (d1 > 0) + (d2 > 0);
    const int nneg = (d0 < 0) + (d1 < 0) + (d2 < 0);
    if (npos == 0 || nneg == 0) return false;

    // The "lone" vertex is the single vertex on one side of the plane.
    const bool lone_pos = (npos == 1);
    auto matches_lone = [lone_pos](double d) {
        return lone_pos ? (d > 0) : (d < 0);
    };

    int a;
    if (matches_lone(d0))      a = 0;
    else if (matches_lone(d1)) a = 1;
    else                       a = 2;
    const int b = (a + 1) % 3;
    const int c = (a + 2) % 3;

    const Vec3d V[3] = {V0, V1, V2};
    const double  d[3] = {d0, d1, d2};

    p_out = edge_plane_intersect(V[a], V[b], d[a], d[b]);
    q_out = edge_plane_intersect(V[a], V[c], d[a], d[c]);
    return true;
}

// Clip the infinite line A + t·dir (dir assumed unit, in T's plane modulo
// fp noise) against triangle T = (V0,V1,V2). Returns the t-interval
// [t_lo, t_hi] over which the line is inside T (within tolerance). If the
// line does not cross T, returns an inverted interval (t_lo > t_hi).
//
// Used by split_triangle_by_cuts to extend cuts with endpoints in T's
// interior — the line is parametrically extended along its direction
// until it exits T through an edge, converting interior endpoints into
// edge endpoints the splitter can consume.
std::pair<double, double> clip_line_in_triangle_plane(
    const Vec3d& V0, const Vec3d& V1, const Vec3d& V2,
    const Vec3d& A, const Vec3d& unit_dir)
{
    // Project to 2D along the axis nearest T's normal (best conditioning).
    const Vec3d N = cross(sub(V1, V0), sub(V2, V0));
    int na = 0;
    if (std::abs(N[1]) > std::abs(N[na])) na = 1;
    if (std::abs(N[2]) > std::abs(N[na])) na = 2;
    const int u = (na + 1) % 3;
    const int v = (na + 2) % 3;

    const double V0u = V0[u], V0v = V0[v];
    const double V1u = V1[u], V1v = V1[v];
    const double V2u = V2[u], V2v = V2[v];
    const double Au  = A[u],  Av  = A[v];
    const double Du  = unit_dir[u], Dv = unit_dir[v];

    std::vector<double> ts;
    ts.reserve(3);

    auto try_edge = [&](double E0u_, double E0v_, double E1u_, double E1v_) {
        // Solve E0 + s·(E1−E0) = A + t·D for s ∈ [0,1], any t.
        const double a = E1u_ - E0u_, b = -Du;
        const double c = E1v_ - E0v_, d = -Dv;
        const double det = a * d - b * c;
        if (std::abs(det) < 1e-18) return; // parallel
        const double r0 = Au - E0u_, r1 = Av - E0v_;
        const double s  = (d * r0 - b * r1) / det;
        const double t  = (a * r1 - c * r0) / det;
        if (s < -1e-9 || s > 1.0 + 1e-9) return;
        ts.push_back(t);
    };
    try_edge(V0u, V0v, V1u, V1v);
    try_edge(V1u, V1v, V2u, V2v);
    try_edge(V2u, V2v, V0u, V0v);

    if (ts.size() < 2) return {1.0, -1.0}; // line missed the triangle
    std::sort(ts.begin(), ts.end());
    return {ts.front(), ts.back()};
}

// Is point P within `tol` of segment (A,B)? If yes, also report its
// barycentric parameter t (0=A, 1=B), clamped to [0,1].
bool point_on_segment(const Vec3d& P, const Vec3d& A, const Vec3d& B,
                      double tol, double& t_out) {
    const Vec3d AB = sub(B, A);
    const double ab_sq = dot(AB, AB);
    if (ab_sq < tol * tol) return false; // degenerate edge
    const Vec3d AP = sub(P, A);
    double t = dot(AP, AB) / ab_sq;
    if (t < -tol || t > 1.0 + tol) return false;
    if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
    const Vec3d proj = add(A, scale(AB, t));
    const Vec3d diff = sub(P, proj);
    if (std::sqrt(dot(diff, diff)) > tol) return false;
    t_out = t;
    return true;
}

// Identify which edge of triangle (V0,V1,V2) a point lies on.
//   0 → edge V0→V1
//   1 → edge V1→V2
//   2 → edge V2→V0
// Returns -1 if none. Tolerance applied via point_on_segment.
int find_edge(const Vec3d& P, const Vec3d& V0, const Vec3d& V1, const Vec3d& V2,
              double tol, double& t_out) {
    double t;
    if (point_on_segment(P, V0, V1, tol, t)) { t_out = t; return 0; }
    if (point_on_segment(P, V1, V2, tol, t)) { t_out = t; return 1; }
    if (point_on_segment(P, V2, V0, tol, t)) { t_out = t; return 2; }
    return -1;
}

// Triangulate a quad (q0,q1,q2,q3) in CCW order by picking the shorter
// diagonal (more numerically stable, avoids thin slivers).
void triangulate_quad(const Vec3d& q0, const Vec3d& q1,
                      const Vec3d& q2, const Vec3d& q3,
                      std::vector<std::array<Vec3d, 3>>& out) {
    const Vec3d d02 = sub(q2, q0);
    const Vec3d d13 = sub(q3, q1);
    const double l02 = dot(d02, d02);
    const double l13 = dot(d13, d13);
    if (l02 <= l13) {
        out.push_back({q0, q1, q2});
        out.push_back({q0, q2, q3});
    } else {
        out.push_back({q0, q1, q3});
        out.push_back({q1, q2, q3});
    }
}

// Split triangle (V0,V1,V2) by a single cut (P→Q) whose endpoints sit on
// edges ep (for P) and eq (for Q). Pre: ep != eq, both in {0,1,2}.
//
// Encoding: edge 0 is V0→V1, edge 1 is V1→V2, edge 2 is V2→V0. The
// "common vertex" is the triangle vertex shared between the two cut edges:
//   ep,eq in {0,1} → common V1, far V0 and V2
//   ep,eq in {1,2} → common V2, far V1 and V0
//   ep,eq in {0,2} → common V0, far V2 and V1
//
// Yields one small triangle (containing the common vertex) plus a quad
// (containing the two far vertices), which is then triangulated.
void split_one_cut(int ep, int eq,
                   const Vec3d& V0, const Vec3d& V1, const Vec3d& V2,
                   const Vec3d& P, const Vec3d& Q,
                   std::vector<std::array<Vec3d, 3>>& out) {
    auto have = [](int e1, int e2, int a, int b) {
        return (e1 == a && e2 == b) || (e1 == b && e2 == a);
    };

    if (have(ep, eq, 0, 1)) {
        // Common vertex V1; P on V0V1 (or Q), Q on V1V2 (or P).
        const Vec3d on01 = (ep == 0) ? P : Q;
        const Vec3d on12 = (eq == 1) ? Q : P;
        out.push_back({on01, V1, on12});
        triangulate_quad(V0, on01, on12, V2, out);
    } else if (have(ep, eq, 1, 2)) {
        // Common vertex V2.
        const Vec3d on12 = (ep == 1) ? P : Q;
        const Vec3d on20 = (eq == 2) ? Q : P;
        out.push_back({on12, V2, on20});
        triangulate_quad(V0, V1, on12, on20, out);
    } else { // have(ep, eq, 0, 2)
        // Common vertex V0.
        const Vec3d on01 = (ep == 0) ? P : Q;
        const Vec3d on20 = (eq == 2) ? Q : P;
        out.push_back({on20, V0, on01});
        triangulate_quad(on01, V1, V2, on20, out);
    }
}

} // anonymous

// ────────────────────────────────────────────────────────────────────────
// tri_tri_intersect — Möller's algorithm, non-coplanar case only.
// ────────────────────────────────────────────────────────────────────────
std::vector<Vec3d> tri_tri_intersect(
    const Vec3d& a1, const Vec3d& b1, const Vec3d& c1,
    const Vec3d& a2, const Vec3d& b2, const Vec3d& c2)
{
    // Plane of T1.
    const Vec3d e1 = sub(b1, a1);
    const Vec3d e2 = sub(c1, a1);
    const Vec3d N1 = cross(e1, e2);
    const double n1_sq = dot(N1, N1);
    if (n1_sq < kPlaneEps) return {}; // degenerate T1
    const double d1 = -dot(N1, a1);

    // Signed distances of T2's vertices to T1's plane.
    const double dU0 = dot(N1, a2) + d1;
    const double dU1 = dot(N1, b2) + d1;
    const double dU2 = dot(N1, c2) + d1;

    // Reject if T2 is entirely on one side of T1's plane.
    if ((dU0 >  kPlaneEps && dU1 >  kPlaneEps && dU2 >  kPlaneEps) ||
        (dU0 < -kPlaneEps && dU1 < -kPlaneEps && dU2 < -kPlaneEps)) return {};
    // Coplanar: skip in v1.
    if (std::abs(dU0) < kPlaneEps && std::abs(dU1) < kPlaneEps &&
        std::abs(dU2) < kPlaneEps) return {};

    // Plane of T2.
    const Vec3d f1 = sub(b2, a2);
    const Vec3d f2 = sub(c2, a2);
    const Vec3d N2 = cross(f1, f2);
    const double n2_sq = dot(N2, N2);
    if (n2_sq < kPlaneEps) return {};
    const double d2 = -dot(N2, a2);

    const double dV0 = dot(N2, a1) + d2;
    const double dV1 = dot(N2, b1) + d2;
    const double dV2 = dot(N2, c1) + d2;

    if ((dV0 >  kPlaneEps && dV1 >  kPlaneEps && dV2 >  kPlaneEps) ||
        (dV0 < -kPlaneEps && dV1 < -kPlaneEps && dV2 < -kPlaneEps)) return {};

    // Segment of T1 crossing T2's plane.
    Vec3d p1, q1;
    if (!tri_plane_cross_segment(a1, b1, c1, dV0, dV1, dV2, p1, q1)) return {};

    // Segment of T2 crossing T1's plane.
    Vec3d p2, q2;
    if (!tri_plane_cross_segment(a2, b2, c2, dU0, dU1, dU2, p2, q2)) return {};

    // Both segments lie on the intersection line of the two planes.
    // Parametrize along D = N1 × N2 and intersect the two intervals.
    const Vec3d D = cross(N1, N2);

    const double tp1 = dot(p1, D);
    const double tq1 = dot(q1, D);
    const double tp2 = dot(p2, D);
    const double tq2 = dot(q2, D);

    const double t1_min = std::min(tp1, tq1);
    const double t1_max = std::max(tp1, tq1);
    const double t2_min = std::min(tp2, tq2);
    const double t2_max = std::max(tp2, tq2);

    const double t_min = std::max(t1_min, t2_min);
    const double t_max = std::min(t1_max, t2_max);
    if (t_min > t_max + kIntervalEps) return {};

    // Map t_min, t_max back to 3D by linear interpolation along p1→q1.
    const double denom = tq1 - tp1;
    auto point_at = [&](double t) {
        const double s = (std::abs(denom) < kIntervalEps)
            ? 0.0 : (t - tp1) / denom;
        return add(p1, scale(sub(q1, p1), s));
    };

    Vec3d pt_lo = point_at(t_min);
    Vec3d pt_hi = point_at(t_max);

    // Reject degenerate (zero-length) segments — they're vertex-on-edge or
    // edge-on-edge touching, not actual interpenetration.
    const Vec3d d = sub(pt_hi, pt_lo);
    if (dot(d, d) < kIntervalEps * kIntervalEps) return {};

    return {pt_lo, pt_hi};
}

// ────────────────────────────────────────────────────────────────────────
// split_triangle_by_cuts — recursive edge-to-edge splitting.
// ────────────────────────────────────────────────────────────────────────
std::vector<std::array<Vec3d, 3>> split_triangle_by_cuts(
    const Vec3d& V0, const Vec3d& V1, const Vec3d& V2,
    const std::vector<std::pair<Vec3d, Vec3d>>& cuts,
    double edge_tol)
{
    std::vector<std::array<Vec3d, 3>> pieces;
    pieces.push_back({V0, V1, V2});

    for (const auto& cut : cuts) {
        // Direction along the cut line.
        const Vec3d dvec = sub(cut.second, cut.first);
        const double dlen = std::sqrt(dot(dvec, dvec));
        if (dlen < 1e-12) continue; // degenerate cut
        const Vec3d unit_dir = scale(dvec, 1.0 / dlen);

        std::vector<std::array<Vec3d, 3>> next;
        next.reserve(pieces.size() * 2);
        for (const auto& tri : pieces) {
            // Fast path: original endpoints already on this piece's edges.
            double t_dummy;
            int ep = find_edge(cut.first,  tri[0], tri[1], tri[2], edge_tol, t_dummy);
            int eq = find_edge(cut.second, tri[0], tri[1], tri[2], edge_tol, t_dummy);
            Vec3d use_P = cut.first;
            Vec3d use_Q = cut.second;

            if (ep < 0 || eq < 0) {
                // At least one endpoint is in the piece's interior.
                // Extend the cut line to the piece's edges.
                auto interval = clip_line_in_triangle_plane(
                    tri[0], tri[1], tri[2], cut.first, unit_dir);
                const double t_lo = interval.first;
                const double t_hi = interval.second;
                if (t_hi <= t_lo + 1e-12) {
                    // Line doesn't enter this piece.
                    next.push_back(tri);
                    continue;
                }
                // Honour the cut's original extent where it's the tighter
                // constraint (i.e. don't extend BEYOND the actual physical
                // intersection segment when that endpoint is already on a
                // piece edge). Only extend endpoints that are interior.
                const double tP_orig = 0.0;
                const double tQ_orig = dlen;
                double tP = (ep < 0) ? t_lo  : tP_orig;
                double tQ = (eq < 0) ? t_hi  : tQ_orig;
                // Clamp to where the line is inside the piece.
                tP = std::max(tP, t_lo);
                tQ = std::min(tQ, t_hi);
                if (tQ - tP < 1e-12) { next.push_back(tri); continue; }
                use_P = add(cut.first, scale(unit_dir, tP));
                use_Q = add(cut.first, scale(unit_dir, tQ));

                // Re-test against piece edges.
                ep = find_edge(use_P, tri[0], tri[1], tri[2], edge_tol, t_dummy);
                eq = find_edge(use_Q, tri[0], tri[1], tri[2], edge_tol, t_dummy);
            }

            if (ep < 0 || eq < 0 || ep == eq) {
                // Extension didn't yield distinct edges (likely the cut runs
                // along an edge of this piece). Skip; keep piece intact.
                next.push_back(tri);
                continue;
            }
            split_one_cut(ep, eq, tri[0], tri[1], tri[2], use_P, use_Q, next);
        }
        pieces.swap(next);
    }
    return pieces;
}

} // namespace meshseal::internal
