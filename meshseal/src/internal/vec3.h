#pragma once

#include <array>
#include <cmath>
#include <cstdint>

namespace meshseal::internal {

using Vec3d = std::array<double, 3>;

inline Vec3d add(Vec3d a, Vec3d b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline Vec3d sub(Vec3d a, Vec3d b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline Vec3d scale(Vec3d a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

inline double dot(Vec3d a, Vec3d b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline Vec3d cross(Vec3d a, Vec3d b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

inline double norm(Vec3d a) {
    return std::sqrt(dot(a, a));
}

inline Vec3d normalize(Vec3d a) {
    const double n = norm(a);
    if (n < 1e-15) {
        return {0.0, 0.0, 0.0};
    }
    return scale(a, 1.0 / n);
}

inline double triangle_area(Vec3d a, Vec3d b, Vec3d c) {
    return 0.5 * norm(cross(sub(b, a), sub(c, a)));
}

inline double signed_tet_volume(Vec3d a, Vec3d b, Vec3d c) {
    return dot(a, cross(b, c)) / 6.0;
}

} // namespace meshseal::internal
