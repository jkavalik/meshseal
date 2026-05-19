#pragma once
#include <meshseal/meshseal.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <utility>

namespace meshseal::test {

// A valid closed cube [0,1]^3: 8 vertices, 12 triangles (outward normals)
inline Mesh unit_cube() {
    Mesh m;
    m.vertices = {
        {0.0, 0.0, 0.0}, // 0
        {1.0, 0.0, 0.0}, // 1
        {1.0, 1.0, 0.0}, // 2
        {0.0, 1.0, 0.0}, // 3
        {0.0, 0.0, 1.0}, // 4
        {1.0, 0.0, 1.0}, // 5
        {1.0, 1.0, 1.0}, // 6
        {0.0, 1.0, 1.0}, // 7
    };
    m.faces = {
        // Bottom (-Z): normal = (0,0,-1)
        {0, 2, 1}, {0, 3, 2},
        // Top (+Z): normal = (0,0,+1)
        {4, 5, 6}, {4, 6, 7},
        // Front (-Y): normal = (0,-1,0)
        {0, 1, 5}, {0, 5, 4},
        // Back (+Y): normal = (0,+1,0)
        {2, 3, 7}, {2, 7, 6},
        // Left (-X): normal = (-1,0,0)
        {0, 4, 7}, {0, 7, 3},
        // Right (+X): normal = (+1,0,0)
        {1, 2, 6}, {1, 6, 5},
    };
    return m;
}

// Open cube — unit cube with top +Z face removed: 8 vertices, 10 triangles
inline Mesh open_cube() {
    Mesh m = unit_cube();
    // Remove the two top (+Z) faces (indices 2 and 3)
    m.faces.erase(m.faces.begin() + 2, m.faces.begin() + 4);
    return m;
}

// Cube with face 0 having flipped winding (inward normal):
// 8 vertices, 12 triangles
inline Mesh cube_one_flipped() {
    Mesh m = unit_cube();
    // Face 0 is {0,2,1}; swap v1 and v2 → {0,1,2} (inward normal on bottom)
    m.faces[0] = {0, 1, 2};
    return m;
}

// Valid tetrahedron: 4 vertices, 4 triangles, outward normals.
// Base triangle sits on z=0 plane; apex at positive z.
inline Mesh unit_tetrahedron() {
    Mesh m;
    const double s3 = std::sqrt(3.0);
    const double s6 = std::sqrt(6.0);
    m.vertices = {
        {0.0,       0.0,     0.0},      // 0
        {1.0,       0.0,     0.0},      // 1
        {0.5,  s3 / 2.0,     0.0},      // 2
        {0.5,  s3 / 6.0, s6 / 3.0},    // 3  apex
    };
    // Outward-normal winding (verified: each face normal points away from
    // the opposite vertex and away from the centroid):
    //   {0,2,1} — base, normal -Z (away from apex at z > 0)
    //   {0,1,3} — side, normal in -Y direction (away from v2 at large y)
    //   {1,2,3} — side, normal in +X direction (away from v0 at x = 0)
    //   {0,3,2} — side, normal in -X direction (away from v1 at x = 1)
    m.faces = {
        {0, 2, 1},
        {0, 1, 3},
        {1, 2, 3},
        {0, 3, 2},
    };
    return m;
}

// Valid tetrahedron + one zero-area degenerate face (all 3 indices = 0)
inline Mesh tetrahedron_with_degenerate() {
    Mesh m = unit_tetrahedron();
    m.faces.push_back({0, 0, 0}); // degenerate: all three vertices identical
    return m;
}

// Valid tetrahedron + an exact duplicate of face 0
inline Mesh tetrahedron_with_duplicate() {
    Mesh m = unit_tetrahedron();
    m.faces.push_back(m.faces[0]);
    return m;
}

// Tetrahedron with all faces' winding inverted (normals pointing inward)
inline Mesh unit_tetrahedron_inverted() {
    auto m = unit_tetrahedron();
    // Flip all faces by swapping v1 and v2
    for (auto& f : m.faces)
        std::swap(f[1], f[2]);
    return m;
}

// Globally inverted cube (all normals inward)
inline Mesh unit_cube_inverted() {
    auto m = unit_cube();
    for (auto& f : m.faces)
        std::swap(f[1], f[2]);
    return m;
}

// Bowtie: two tetrahedra sharing only vertex 0
// Vertices 0-3: first tetrahedron; vertices 4-6: second tetrahedron; vertex 0 shared
inline Mesh bowtie_tetrahedra() {
    // First tet: v0, v1, v2, v3
    // Second tet: v0, v4, v5, v6
    Mesh m;
    m.vertices = {
        {0.0, 0.0, 0.0},  // shared bowtie vertex
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {-1.0, 0.0, 0.0},
        {0.0, -1.0, 0.0},
        {0.0, 0.0, -1.0},
    };
    // First tet faces (outward normals): v0,v1,v2 / v0,v3,v1 / v0,v2,v3 / v1,v3,v2
    m.faces = {
        {0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3},
        // Second tet: v0,v4,v5 / v0,v6,v4 / v0,v5,v6 / v4,v6,v5
        {0, 4, 5}, {0, 6, 4}, {0, 5, 6}, {4, 6, 5}
    };
    return m;
}

// Cube with a fin: unit_cube plus one extra triangle sharing edge (0,1),
// making that edge non-manifold (shared by 3 faces instead of 2).
// unit_cube face {0,2,1} has directed edge 1->0 on (0,1).
// unit_cube face {0,1,5} has directed edge 0->1 on (0,1).
// The fin face {0,1,new_v} also has directed edge 0->1 on (0,1) -> non-manifold.
inline Mesh cube_with_fin() {
    auto m = unit_cube();
    // New vertex outside the cube bottom face
    m.vertices.push_back({0.5, -0.5, 0.0});
    uint32_t new_v = static_cast<uint32_t>(m.vertices.size()) - 1; // index 8
    m.faces.push_back({0, 1, new_v}); // fin face: directed edge 0->1 duplicates face {0,1,5}
    return m;
}

// Two separate unit cubes: one at origin, one at (10, 10, 10) — far apart, 2 shells
inline Mesh two_separate_cubes() {
    auto m1 = unit_cube();
    auto m2 = unit_cube();
    // Offset m2 vertices by (10, 10, 10)
    const uint32_t offset = static_cast<uint32_t>(m1.vertices.size());
    for (auto& v : m2.vertices)
        v = {v[0] + 10.0, v[1] + 10.0, v[2] + 10.0};
    for (auto& f : m2.faces)
        f = {f[0] + offset, f[1] + offset, f[2] + offset};
    m1.vertices.insert(m1.vertices.end(), m2.vertices.begin(), m2.vertices.end());
    m1.faces.insert(m1.faces.end(), m2.faces.begin(), m2.faces.end());
    return m1;
}

// One big cube + one tiny cube (tiny should be dropped by volume threshold)
inline Mesh big_and_tiny_cube() {
    auto big = unit_cube(); // volume ~1
    Mesh tiny;
    // Tiny cube at (5,5,5) with side 0.001
    double s = 0.001;
    tiny.vertices = {
        {5.0,   5.0,   5.0},   {5.0+s, 5.0,   5.0},
        {5.0+s, 5.0+s, 5.0},   {5.0,   5.0+s, 5.0},
        {5.0,   5.0,   5.0+s}, {5.0+s, 5.0,   5.0+s},
        {5.0+s, 5.0+s, 5.0+s}, {5.0,   5.0+s, 5.0+s}
    };
    // Same winding structure as unit_cube (normals all outward) — face indices differ
    tiny.faces = {
        {0,3,1},{1,3,2}, {4,5,7},{5,6,7},
        {0,1,4},{1,5,4}, {2,3,6},{3,7,6},
        {0,4,3},{3,4,7}, {1,2,5},{2,6,5}
    };
    const uint32_t off = static_cast<uint32_t>(big.vertices.size());
    for (auto& f : tiny.faces)
        f = {f[0]+off, f[1]+off, f[2]+off};
    big.vertices.insert(big.vertices.end(), tiny.vertices.begin(), tiny.vertices.end());
    big.faces.insert(big.faces.end(), tiny.faces.begin(), tiny.faces.end());
    return big;
}

// Tetrahedron missing one face (open, 1 triangular hole)
inline Mesh open_tetrahedron() {
    Mesh m;
    m.vertices = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
        {0.5, 0.866025, 0.0}, {0.5, 0.288675, 0.816497}
    };
    m.faces = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}}; // base missing
    return m;
}

// Two overlapping tetrahedra — self-intersecting mesh
inline Mesh two_overlapping_tetrahedra() {
    Mesh m;
    // Tet 1: (0,0,0), (2,0,0), (1,2,0), (1,1,2)
    m.vertices = {
        {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {1.0, 2.0, 0.0}, {1.0, 1.0, 2.0},
        // Tet 2: offset so it overlaps tet 1
        {0.5, 0.5, 0.5}, {2.5, 0.5, 0.5}, {1.5, 2.5, 0.5}, {1.5, 1.5, 2.5}
    };
    m.faces = {
        // Tet 1 faces (outward normals):
        {0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2},
        // Tet 2 faces (outward normals):
        {4, 6, 5}, {4, 5, 7}, {5, 6, 7}, {4, 7, 6}
    };
    return m;
}

// Two back-to-back triangles: a fin (ultra-thin feature)
// Triangle at z=0 facing up and triangle at z=0.0001 facing down, sharing same area
inline Mesh back_to_back_triangles() {
    Mesh m;
    m.vertices = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, 1.0, 0.0},           // face 0
        {0.0, 0.0, 0.0005}, {1.0, 0.0, 0.0005}, {0.5, 1.0, 0.0005} // face 1: z=5e-4 > weld tol, ≤ thin threshold
    };
    m.faces = {
        {0, 1, 2},   // normal pointing up
        {3, 5, 4}    // normal pointing down (reversed winding)
    };
    return m;
}

// Triangle soup: a cube's 12 faces with no shared vertices (each face is independent).
// Every triangle has correct outward winding; the non-manifold property comes solely
// from open boundary edges (no two adjacent faces share a vertex index).
// Manifold's soup constructor welds the duplicate positions and produces a valid solid.
inline Mesh cube_triangle_soup() {
    auto m = unit_cube();
    // Expand: each face gets its own 3 vertices (no vertex sharing)
    Mesh soup;
    for (const auto& f : m.faces) {
        uint32_t base = static_cast<uint32_t>(soup.vertices.size());
        soup.vertices.push_back(m.vertices[f[0]]);
        soup.vertices.push_back(m.vertices[f[1]]);
        soup.vertices.push_back(m.vertices[f[2]]);
        soup.faces.push_back({base, base+1, base+2});
    }
    return soup;
}

// ─── Multi-component and overlap mesh fixtures ────────────────────────────────

// Generic axis-aligned box from (x0,y0,z0) to (x1,y1,z1), outward normals.
// Uses the same vertex-ordering and face winding as unit_cube().
inline Mesh make_box(double x0, double y0, double z0,
                     double x1, double y1, double z1) {
    Mesh m;
    m.vertices = {
        {x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},  // 0-3 bottom ring
        {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1},  // 4-7 top ring
    };
    m.faces = {
        {0, 2, 1}, {0, 3, 2},  // -Z (bottom)
        {4, 5, 6}, {4, 6, 7},  // +Z (top)
        {0, 1, 5}, {0, 5, 4},  // -Y (front)
        {2, 3, 7}, {2, 7, 6},  // +Y (back)
        {0, 4, 7}, {0, 7, 3},  // -X (left)
        {1, 2, 6}, {1, 6, 5},  // +X (right)
    };
    return m;
}

// Two overlapping closed cubes: A=[0,1]^3, B=[0.5,1.5]×[0,1]×[0,1].
// The overlap region [0.5,1]×[0,1]×[0,1] has volume 0.5; the boolean union
// is a 1.5×1×1 box with volume 1.5.
// Both shells are valid closed manifolds; the combined mesh is a 2-component
// non-planar arrangement with geometrically intersecting faces.
// Represents kitchen_sink, multi_defect_box, overlapping_shells_with_hole
// after their holes have been filled.
inline Mesh two_overlapping_cubes() {
    auto a = make_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    auto b = make_box(0.5, 0.0, 0.0, 1.5, 1.0, 1.0);
    const uint32_t off = static_cast<uint32_t>(a.vertices.size());
    for (auto& f : b.faces)
        f = {f[0] + off, f[1] + off, f[2] + off};
    a.vertices.insert(a.vertices.end(), b.vertices.begin(), b.vertices.end());
    a.faces.insert(a.faces.end(),       b.faces.begin(),    b.faces.end());
    return a;
}

// Small cube [1,3]^3 sitting completely inside a large cube [0,4]^3.
// Both shells are valid closed manifolds; their boolean union equals the
// outer cube (volume 64) because every point of the inner is already in the outer.
// Represents nested real-world shell patterns (matryoshka_nested, nested_shells).
inline Mesh nested_cubes() {
    auto outer = make_box(0.0, 0.0, 0.0, 4.0, 4.0, 4.0);
    auto inner = make_box(1.0, 1.0, 1.0, 3.0, 3.0, 3.0);
    const uint32_t off = static_cast<uint32_t>(outer.vertices.size());
    for (auto& f : inner.faces)
        f = {f[0] + off, f[1] + off, f[2] + off};
    outer.vertices.insert(outer.vertices.end(), inner.vertices.begin(), inner.vertices.end());
    outer.faces.insert(outer.faces.end(),       inner.faces.begin(),    inner.faces.end());
    return outer;
}

// Three 2×2×2 boxes with pairwise chain overlaps:
//   A=[0,2]×[0,2]×[0,2],  B=[1,3]×[0,2]×[0,2],  C=[2,4]×[0,2]×[0,2]
// A∩B=[1,2]×... (vol 4), B∩C=[2,3]×... (vol 4), A∩C={} (zero volume).
// Boolean union = [0,4]×[0,2]×[0,2] = volume 16.
// Represents multi-shell real-world models (clown, doorman, human, o3d_monkey)
// where several components overlap pairwise and must all be merged.
inline Mesh three_chain_overlapping_boxes() {
    auto a = make_box(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
    auto b = make_box(1.0, 0.0, 0.0, 3.0, 2.0, 2.0);
    auto c = make_box(2.0, 0.0, 0.0, 4.0, 2.0, 2.0);
    auto append = [](Mesh& dst, Mesh src) {
        const uint32_t off = static_cast<uint32_t>(dst.vertices.size());
        for (auto& f : src.faces)
            f = {f[0] + off, f[1] + off, f[2] + off};
        dst.vertices.insert(dst.vertices.end(), src.vertices.begin(), src.vertices.end());
        dst.faces.insert(dst.faces.end(),       src.faces.begin(),    src.faces.end());
    };
    append(a, b);
    append(a, c);
    return a;
}

// Two overlapping cubes where one has a hole: an open_cube [0,1]^3 (top face
// removed) plus a complete cube [0.5,1.5]×[0,1]×[0,1].
// The shells overlap in [0.5,1]×[0,1]×[0,1] (vol 0.5).
// After hole-filling the open cube, the ideal repair result is a single merged
// solid with volume 1.5.
// Represents overlapping_shells_with_hole, kitchen_sink, clown, doorman —
// combined-defect cases where orig_has_open_boundaries causes the intersection
// stage guard to skip merging the now-closed overlapping shells.
inline Mesh overlapping_cubes_one_holed() {
    Mesh a = open_cube();  // [0,1]^3 with top (+Z) face removed: 10 faces, 8 vertices
    Mesh b = make_box(0.5, 0.0, 0.0, 1.5, 1.0, 1.0);
    const uint32_t off = static_cast<uint32_t>(a.vertices.size());
    for (auto& f : b.faces)
        f = {f[0] + off, f[1] + off, f[2] + off};
    a.vertices.insert(a.vertices.end(), b.vertices.begin(), b.vertices.end());
    a.faces.insert(a.faces.end(),       b.faces.begin(),    b.faces.end());
    return a;
}

} // namespace meshseal::test

// ─── Shared topology / geometry helpers ──────────────────────────────────────
// These are available to every test file via #include "helpers.h".
// NOTE: All helpers below assume no degenerate faces (no repeated vertex indices
// within a single face). Results are undefined on degenerate input.

namespace meshseal::test {

// Compute signed volume via the divergence theorem.
// Positive  → faces have outward normals.
// Negative  → faces have inward normals.
// Meaningful only for a closed, consistently-wound mesh.
inline double signed_volume(const Mesh& m) {
    double vol = 0.0;
    for (const auto& f : m.faces) {
        const auto& a = m.vertices[f[0]];
        const auto& b = m.vertices[f[1]];
        const auto& c = m.vertices[f[2]];
        vol += (a[0]*(b[1]*c[2]-b[2]*c[1])
              + a[1]*(b[2]*c[0]-b[0]*c[2])
              + a[2]*(b[0]*c[1]-b[1]*c[0]));
    }
    return vol / 6.0;
}

// Returns true iff every undirected edge is used by exactly two faces with
// *opposite* directed orientations (a→b in one, b→a in the other).
//
// This is a strictly stronger check than signed_volume > 0:
//   • signed_volume is a global aggregate — misoriented faces can be hidden
//     by the surrounding majority (e.g. 1 flipped face on a cube still gives vol ≈ +1).
//   • is_consistently_wound catches every individual winding fault because a single
//     flipped face breaks the directed-edge pairing on each of its 3 edges.
//
// Requires: closed, 2-manifold mesh. (Open/non-manifold meshes will return false. Returns true vacuously for an empty mesh.)
inline bool is_consistently_wound(const Mesh& m) {
    std::map<std::pair<uint32_t,uint32_t>, int> directed;
    for (const auto& f : m.faces)
        for (int i = 0; i < 3; ++i)
            directed[{f[i], f[(i+1)%3]}]++;
    for (const auto& [e, cnt] : directed) {
        if (cnt != 1) return false;
        auto rev = std::make_pair(e.second, e.first);
        auto it  = directed.find(rev);
        if (it == directed.end() || it->second != 1) return false;
    }
    return true;
}

// Returns the maximum number of faces incident to any single undirected edge.
// • 2 = 2-manifold (every edge shared by exactly 2 faces).
// • ≥3 = non-manifold edge present.
// • 1 = boundary edge present (open mesh).
inline uint32_t max_edge_valence(const Mesh& m) {
    std::map<std::pair<uint32_t,uint32_t>, uint32_t> cnt;
    for (const auto& f : m.faces)
        for (int i = 0; i < 3; ++i) {
            auto a = f[i], b = f[(i+1)%3];
            if (a > b) std::swap(a, b);
            cnt[{a,b}]++;
        }
    uint32_t mx = 0;
    for (const auto& [e, c] : cnt) mx = std::max(mx, c);
    return mx;
}

// Returns the number of boundary edges (edges incident to exactly 1 face).
// Zero means the mesh is closed.
inline uint32_t boundary_edge_count(const Mesh& m) {
    std::map<std::pair<uint32_t,uint32_t>, uint32_t> cnt;
    for (const auto& f : m.faces)
        for (int i = 0; i < 3; ++i) {
            auto a = f[i], b = f[(i+1)%3];
            if (a > b) std::swap(a, b);
            cnt[{a,b}]++;
        }
    uint32_t bdy = 0;
    for (const auto& [e, c] : cnt) if (c == 1) ++bdy;
    return bdy;
}

} // namespace meshseal::test (topology helpers)
