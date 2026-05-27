#include "stl_io.h"
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace meshseal {

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// The binary STL format is defined as little-endian (4-byte floats + uint32
// counts). The implementations below are explicit byte-order — works on big-
// endian hosts as well as little-endian. Modern compilers (GCC, Clang, MSVC)
// recognize the byte-OR / byte-store patterns and emit a single load/store
// plus a conditional bswap, so there is no perf cost on x86/ARM (LE) and the
// code stays correct on the rare BE host (POWER BE, MIPS BE).

static uint32_t read_u32_le(const uint8_t* p) {
    return  static_cast<uint32_t>(p[0])
         | (static_cast<uint32_t>(p[1]) <<  8)
         | (static_cast<uint32_t>(p[2]) << 16)
         | (static_cast<uint32_t>(p[3]) << 24);
}

static void write_u32_le(uint8_t* p, uint32_t v) {
    p[0] = static_cast<uint8_t>( v        & 0xffu);
    p[1] = static_cast<uint8_t>((v >>  8) & 0xffu);
    p[2] = static_cast<uint8_t>((v >> 16) & 0xffu);
    p[3] = static_cast<uint8_t>((v >> 24) & 0xffu);
}

static float read_f32_le(const uint8_t* p) {
    const uint32_t u = read_u32_le(p);
    float v;
    std::memcpy(&v, &u, sizeof(v));  // type-pun without UB
    return v;
}

static void write_f32_le(uint8_t* p, float v) {
    uint32_t u;
    std::memcpy(&u, &v, sizeof(u));
    write_u32_le(p, u);
}

// ---------------------------------------------------------------------------
// ASCII STL parsing
// ---------------------------------------------------------------------------

// Parse an ASCII STL: scan for `vertex x y z` records; every 3 form a face.
// Tolerant of arbitrary whitespace/indentation and ignores the normal,
// `facet`/`loop`/`endloop`/`endfacet`/`endsolid` keywords.
static Mesh read_stl_ascii(const uint8_t* data, size_t size) {
    // Match the binary STL cap so an attacker can't bypass the binary cap
    // by including a "solid" prefix and routing into the ASCII path.
    // 100M triangles = 300M vertex records.
    constexpr size_t kMaxAsciiVerts = 300'000'000u;
    std::string text(reinterpret_cast<const char*>(data), size);
    std::istringstream ss(text);
    std::string tok;
    std::vector<std::array<double, 3>> verts;
    while (ss >> tok) {
        if (tok == "vertex") {
            double x, y, z;
            if (!(ss >> x >> y >> z)) {
                throw StlError("malformed ASCII STL: bad vertex record");
            }
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
                throw StlError("NaN/Inf vertex in ASCII STL");
            }
            if (verts.size() >= kMaxAsciiVerts) {
                throw StlError("ASCII STL vertex count exceeds 300M cap");
            }
            verts.push_back({x, y, z});
        }
    }
    if (verts.empty()) {
        throw StlError("ASCII STL contains no vertices");
    }
    if (verts.size() % 3 != 0) {
        throw StlError("malformed ASCII STL: vertex count not a multiple of 3");
    }
    Mesh mesh;
    mesh.vertices = std::move(verts);
    const uint32_t nt = static_cast<uint32_t>(mesh.vertices.size() / 3);
    mesh.faces.reserve(nt);
    for (uint32_t i = 0; i < nt; ++i) {
        mesh.faces.push_back({i*3, i*3 + 1, i*3 + 2});
    }
    return mesh;
}

// Decide whether the buffer is an ASCII STL. A binary STL is exactly
// 84 + 50·count bytes (count = u32 at offset 80); an ASCII file's offset-80
// bytes are text, so the size will not match. The `solid` prefix alone is
// not decisive (binary headers may also start with "solid"), but combined
// with a size mismatch — or the literal `facet` keyword appearing early —
// it is reliable.
static bool looks_like_ascii_stl(const uint8_t* data, size_t size) {
    size_t s = 0;
    while (s < size && std::isspace(static_cast<unsigned char>(data[s]))) ++s;
    if (size - s < 5 || std::memcmp(data + s, "solid", 5) != 0) return false;
    if (size < 84) return true;  // too small to be a sane binary STL
    const uint32_t tc = read_u32_le(data + 80);
    const size_t bin_expected = 84 + static_cast<size_t>(tc) * 50;
    if (size != bin_expected) return true;
    // Sizes coincide (extremely rare) — fall back to scanning for `facet`.
    const size_t scan = size < 512 ? size : 512;
    for (size_t i = s; i + 5 <= scan; ++i)
        if (std::memcmp(data + i, "facet", 5) == 0) return true;
    return false;
}

// ---------------------------------------------------------------------------
// read_stl_bytes — core parsing with full validation
// ---------------------------------------------------------------------------

Mesh read_stl_bytes(const uint8_t* data, size_t size) {
    if (looks_like_ascii_stl(data, size)) {
        return read_stl_ascii(data, size);
    }

    if (size < 84) {
        throw StlError("file too small");
    }

    const uint32_t tri_count = read_u32_le(data + 80);

    if (tri_count > 100'000'000u) {
        std::ostringstream oss;
        oss << "implausible triangle count: " << tri_count;
        throw StlError(oss.str());
    }

    const size_t expected = 84 + static_cast<size_t>(tri_count) * 50;
    if (size != expected) {
        std::ostringstream oss;
        oss << "size mismatch: expected " << expected << " bytes, got " << size;
        throw StlError(oss.str());
    }

    Mesh mesh;
    mesh.vertices.reserve(static_cast<size_t>(tri_count) * 3);
    mesh.faces.reserve(tri_count);

    const uint8_t* p = data + 84;
    for (uint32_t i = 0; i < tri_count; ++i, p += 50) {
        // Skip normal (bytes 0-11), read 3 vertices (bytes 12-47)
        for (int v = 0; v < 3; ++v) {
            const uint8_t* vp = p + 12 + v * 12;
            float x = read_f32_le(vp);
            float y = read_f32_le(vp + 4);
            float z = read_f32_le(vp + 8);
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
                std::ostringstream oss;
                oss << "NaN/Inf vertex in triangle " << i;
                throw StlError(oss.str());
            }
            mesh.vertices.push_back({
                static_cast<double>(x),
                static_cast<double>(y),
                static_cast<double>(z)
            });
        }
        const uint32_t base = i * 3;
        mesh.faces.push_back({base, base + 1, base + 2});
    }

    return mesh;
}

// ---------------------------------------------------------------------------
// read_stl — file-based wrapper
// ---------------------------------------------------------------------------

Mesh read_stl(const std::filesystem::path& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) {
        throw StlError("cannot open file: " + path.string());
    }
    const std::streamoff sz = f.tellg();
    if (sz < 0) {
        throw StlError("cannot determine file size: " + path.string());
    }
    const auto file_size = static_cast<size_t>(sz);
    f.seekg(0, std::ios::beg);

    std::vector<uint8_t> buf(file_size);
    if (!f.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(file_size))) {
        throw StlError("read error: " + path.string());
    }

    return read_stl_bytes(buf.data(), buf.size());
}

// ---------------------------------------------------------------------------
// write_stl_bytes — serialize to buffer
// ---------------------------------------------------------------------------

std::vector<uint8_t> write_stl_bytes(const Mesh& mesh) {
    if (mesh.faces.size() > static_cast<size_t>(UINT32_MAX)) {
        throw StlError("mesh has too many faces for STL format");
    }
    const uint32_t tri_count = static_cast<uint32_t>(mesh.faces.size());
    const size_t total = 84 + static_cast<size_t>(tri_count) * 50;

    std::vector<uint8_t> buf(total, 0);

    write_u32_le(buf.data() + 80, tri_count);

    uint8_t* p = buf.data() + 84;
    for (uint32_t i = 0; i < tri_count; ++i, p += 50) {
        const auto& face = mesh.faces[i];

        // Compute face normal from vertex order (cross product)
        const auto& a = mesh.vertices[face[0]];
        const auto& b = mesh.vertices[face[1]];
        const auto& c = mesh.vertices[face[2]];

        const double ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
        const double vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2];
        const double nx = uy * vz - uz * vy;
        const double ny = uz * vx - ux * vz;
        const double nz = ux * vy - uy * vx;
        const double len = std::sqrt(nx * nx + ny * ny + nz * nz);

        if (len > 0.0) {
            write_f32_le(p,     static_cast<float>(nx / len));
            write_f32_le(p + 4, static_cast<float>(ny / len));
            write_f32_le(p + 8, static_cast<float>(nz / len));
        } // else degenerate — buf already zero

        for (int v = 0; v < 3; ++v) {
            const auto& vert = mesh.vertices[face[static_cast<size_t>(v)]];
            uint8_t* vp = p + 12 + v * 12;
            write_f32_le(vp,     static_cast<float>(vert[0]));
            write_f32_le(vp + 4, static_cast<float>(vert[1]));
            write_f32_le(vp + 8, static_cast<float>(vert[2]));
        }
        // attribute byte count already zero (buf initialised to 0)
    }

    return buf;
}

// ---------------------------------------------------------------------------
// write_stl — file-based wrapper
// ---------------------------------------------------------------------------

void write_stl(const Mesh& mesh, const std::filesystem::path& path) {
    const auto buf = write_stl_bytes(mesh);
    std::ofstream f(path, std::ios::binary);
    if (!f) {
        throw StlError("cannot open file for writing: " + path.string());
    }
    if (!f.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size()))) {
        throw StlError("write error: " + path.string());
    }
}

} // namespace meshseal
