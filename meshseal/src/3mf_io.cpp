#include "3mf_io.h"
#include <miniz.h>
#include <fast_float/fast_float.h>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace meshseal {

// miniz declares MZ_DEFAULT_COMPRESSION etc. in an anonymous enum (signed),
// but its API expects mz_uint (unsigned). Cast once into a typed constant
// so the call sites don't trip MSVC's C4245 signed/unsigned-mismatch warning.
static constexpr mz_uint kMzDefaultCompression =
    static_cast<mz_uint>(MZ_DEFAULT_COMPRESSION);

// Zip-bomb guard: reject 3MF archives whose declared uncompressed size for
// the model XML or slicer config exceeds these caps. A 1KB malicious 3MF
// can otherwise advertise a 10GB uncompressed entry and force allocation
// of the full size during extract_to_heap.
static constexpr mz_uint64 kMaxModelXmlUncompressed     = 256ull * 1024 * 1024; // 256 MB
static constexpr mz_uint64 kMaxSlicerConfigUncompressed =        1 * 1024 * 1024; //   1 MB

// ---------------------------------------------------------------------------
// Locale-independent numeric parsers.
//
// `std::stod` / `std::stoi` / `std::stoul` honour the C locale: on Czech /
// German Windows, decimal separator is `,` and `std::stod("1.5")` parses
// as `1`. The 3MF spec mandates `.` as the decimal separator, so the
// reader MUST be locale-insensitive. `std::from_chars` is the right tool.
// ---------------------------------------------------------------------------

static bool parse_double(const std::string& s, double& out) {
    if (s.empty()) return false;
    const char* first = s.data();
    const char* last  = s.data() + s.size();
    // Use fast_float instead of std::from_chars: the C++17 fp overload is
    // missing on Apple Clang's libc++ (explicitly = deleted) and on
    // libstdc++ < 11 (GCC 9/10). fast_float ships the polyfill with the
    // exact same API.
    auto r = fast_float::from_chars(first, last, out);
    return r.ec == std::errc() && r.ptr == last;
}

static bool parse_uint32(const std::string& s, uint32_t& out) {
    if (s.empty()) return false;
    const char* first = s.data();
    const char* last  = s.data() + s.size();
    unsigned long long tmp = 0;
    auto r = std::from_chars(first, last, tmp);
    if (r.ec != std::errc() || r.ptr != last) return false;
    if (tmp > std::numeric_limits<uint32_t>::max()) return false;
    out = static_cast<uint32_t>(tmp);
    return true;
}

static bool parse_int(const std::string& s, int& out) {
    if (s.empty()) return false;
    const char* first = s.data();
    const char* last  = s.data() + s.size();
    auto r = std::from_chars(first, last, out);
    return r.ec == std::errc() && r.ptr == last;
}

// ---------------------------------------------------------------------------
// Internal XML helpers
// ---------------------------------------------------------------------------

// Extract the value of a named attribute from a single XML tag string.
// Looks for: ' name="value"' (with leading space to avoid partial matches like fv1).
static std::string get_attr(const std::string& tag, const std::string& name) {
    const std::string key = " " + name + "=\"";
    size_t pos = tag.find(key);
    if (pos != std::string::npos) {
        pos += key.size();
        const size_t end = tag.find('"', pos);
        if (end != std::string::npos) {
            return tag.substr(pos, end - pos);
        }
    }
    // try single-quote variant
    const std::string key2 = " " + name + "='";
    pos = tag.find(key2);
    if (pos != std::string::npos) {
        pos += key2.size();
        const size_t end = tag.find('\'', pos);
        if (end != std::string::npos) {
            return tag.substr(pos, end - pos);
        }
    }
    return {};
}

// ---------------------------------------------------------------------------
// Content strings used when writing
// ---------------------------------------------------------------------------

static const char* k_content_types =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\n"
    "  <Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\n"
    "  <Default Extension=\"model\" ContentType=\"application/vnd.ms-package.3dmanufacturing-3dmodel+xml\"/>\n"
    "</Types>\n";

static const char* k_rels =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n"
    "  <Relationship Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel\""
    " Target=\"/3D/3dmodel.model\" Id=\"rel0\"/>\n"
    "</Relationships>\n";

// ---------------------------------------------------------------------------
// ZIP packing — shared by file and in-memory writers
// ---------------------------------------------------------------------------

// Add the standard 3MF entries ([Content_Types].xml, _rels/.rels,
// 3D/3dmodel.model, and optionally Metadata/Slic3r_PE_model.config) to an
// already-initialised writer. Returns false on any miniz failure.
static bool add_3mf_entries(mz_zip_archive& zip,
                            const std::string& model_xml,
                            const std::string* slic3r_cfg) {
    bool ok = true;
    ok = ok && mz_zip_writer_add_mem(&zip, "[Content_Types].xml",
        k_content_types, std::strlen(k_content_types), kMzDefaultCompression);
    ok = ok && mz_zip_writer_add_mem(&zip, "_rels/.rels",
        k_rels, std::strlen(k_rels), kMzDefaultCompression);
    ok = ok && mz_zip_writer_add_mem(&zip, "3D/3dmodel.model",
        model_xml.data(), model_xml.size(), kMzDefaultCompression);
    if (slic3r_cfg && !slic3r_cfg->empty()) {
        ok = ok && mz_zip_writer_add_mem(&zip, "Metadata/Slic3r_PE_model.config",
            slic3r_cfg->data(), slic3r_cfg->size(), kMzDefaultCompression);
    }
    return ok;
}

static void pack_3mf_to_file(const std::string& model_xml,
                             const std::string* slic3r_cfg,
                             const std::filesystem::path& path) {
    mz_zip_archive zip{};
    if (!mz_zip_writer_init_file(&zip, path.string().c_str(), 0)) {
        throw ThreeMfError("cannot create 3MF file: " + path.string());
    }
    bool ok = add_3mf_entries(zip, model_xml, slic3r_cfg);
    ok = ok && mz_zip_writer_finalize_archive(&zip);
    mz_zip_writer_end(&zip);
    if (!ok) {
        throw ThreeMfError("error writing 3MF archive: " + path.string());
    }
}

static std::vector<uint8_t> pack_3mf_to_buffer(const std::string& model_xml,
                                               const std::string* slic3r_cfg) {
    mz_zip_archive zip{};
    if (!mz_zip_writer_init_heap(&zip, 0, 0)) {
        throw ThreeMfError("cannot initialise in-memory 3MF writer");
    }
    bool ok = add_3mf_entries(zip, model_xml, slic3r_cfg);
    void* buf = nullptr;
    size_t buf_size = 0;
    ok = ok && mz_zip_writer_finalize_heap_archive(&zip, &buf, &buf_size);
    std::vector<uint8_t> out;
    if (ok && buf) {
        const uint8_t* p = static_cast<const uint8_t*>(buf);
        out.assign(p, p + buf_size);   // copy out before end() frees the heap buffer
    }
    mz_zip_writer_end(&zip);           // frees the internal heap buffer
    if (!ok) {
        throw ThreeMfError("error writing 3MF archive to buffer");
    }
    return out;
}

// ---------------------------------------------------------------------------
// model.xml / config builders — shared by file and in-memory writers
// ---------------------------------------------------------------------------

static std::string build_single_model_xml(const Mesh& mesh) {
    std::string model;
    model.reserve(64 + mesh.vertices.size() * 60 + mesh.faces.size() * 40);

    model +=
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<model unit=\"millimeter\" xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\">\n"
        "  <resources>\n"
        "    <object id=\"1\" type=\"model\">\n"
        "      <mesh>\n"
        "        <vertices>\n";

    char buf[128];
    for (const auto& v : mesh.vertices) {
        std::snprintf(buf, sizeof(buf),
            "          <vertex x=\"%.17g\" y=\"%.17g\" z=\"%.17g\"/>\n",
            v[0], v[1], v[2]);
        model += buf;
    }

    model += "        </vertices>\n        <triangles>\n";

    for (const auto& f : mesh.faces) {
        std::snprintf(buf, sizeof(buf),
            "          <triangle v1=\"%u\" v2=\"%u\" v3=\"%u\"/>\n",
            f[0], f[1], f[2]);
        model += buf;
    }

    model +=
        "        </triangles>\n"
        "      </mesh>\n"
        "    </object>\n"
        "  </resources>\n"
        "  <build>\n"
        "    <item objectid=\"1\"/>\n"
        "  </build>\n"
        "</model>\n";
    return model;
}

static std::string build_volumes_model_xml(const std::vector<ThreeMfVolume>& volumes) {
    std::string model;
    // Rough pre-allocation: one <vertex .../> per vertex is ~60 chars,
    // one <triangle .../> per face is ~40 chars; plus per-object overhead.
    size_t total_v = 0, total_f = 0;
    for (const auto& vol : volumes) {
        total_v += vol.mesh.vertices.size();
        total_f += vol.mesh.faces.size();
    }
    model.reserve(256 + total_v * 60 + total_f * 40 + volumes.size() * 200);

    model +=
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<model unit=\"millimeter\" xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\">\n"
        "  <resources>\n";

    char buf[128];
    for (size_t i = 0; i < volumes.size(); ++i) {
        const auto& vol = volumes[i];
        // Object IDs start at 1 (3MF convention).
        std::snprintf(buf, sizeof(buf),
            "    <object id=\"%zu\" type=\"model\">\n"
            "      <mesh>\n"
            "        <vertices>\n", i + 1);
        model += buf;
        for (const auto& v : vol.mesh.vertices) {
            std::snprintf(buf, sizeof(buf),
                "          <vertex x=\"%.17g\" y=\"%.17g\" z=\"%.17g\"/>\n",
                v[0], v[1], v[2]);
            model += buf;
        }
        model += "        </vertices>\n        <triangles>\n";
        for (const auto& f : vol.mesh.faces) {
            std::snprintf(buf, sizeof(buf),
                "          <triangle v1=\"%u\" v2=\"%u\" v3=\"%u\"/>\n",
                f[0], f[1], f[2]);
            model += buf;
        }
        model +=
            "        </triangles>\n"
            "      </mesh>\n"
            "    </object>\n";
    }

    model +=
        "  </resources>\n"
        "  <build>\n";
    for (size_t i = 0; i < volumes.size(); ++i) {
        std::snprintf(buf, sizeof(buf),
            "    <item objectid=\"%zu\"/>\n", i + 1);
        model += buf;
    }
    model +=
        "  </build>\n"
        "</model>\n";
    return model;
}

// Slic3r_PE_model.config preserving the volume partition. The PrusaSlicer
// slicer reader expects ALL volumes packed into one <object> with
// firstid/lastid pointing at concatenated triangle indices. This is the
// inverse of the read_3mf_volumes parsing.
static std::string build_slic3r_cfg(const std::vector<ThreeMfVolume>& volumes) {
    std::string slic3r_cfg;
    slic3r_cfg.reserve(256 + volumes.size() * 200);
    slic3r_cfg += "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<config>\n";
    // Single <object id="1"> entry — concatenated face indices.
    // Note: this is provided as a hint for downstream PrusaSlicer round-
    // trip; our PRIMARY output is the multi-<object> 3MF above (each
    // slicer object is a separate physical body). PrusaSlicer can read
    // either form; multi-<object> matches what e.g. write_3mf gives.
    slic3r_cfg += "  <object id=\"1\" instances_count=\"1\">\n";
    char buf[128];
    uint32_t cur_first = 0;
    for (const auto& vol : volumes) {
        const uint32_t nf = static_cast<uint32_t>(vol.mesh.faces.size());
        if (nf == 0) continue;
        std::snprintf(buf, sizeof(buf),
            "    <volume firstid=\"%u\" lastid=\"%u\">\n",
            cur_first, cur_first + nf - 1);
        slic3r_cfg += buf;
        // Escape the name minimally (no quotes / angle brackets in
        // typical Slic3r names; PrusaSlicer's names look like
        // "Bee_v3.stl_1_4"). Pass-through for now.
        slic3r_cfg += "      <metadata type=\"volume\" key=\"name\" value=\"";
        for (char c : vol.name) {
            if (c == '"' || c == '<' || c == '>' || c == '&') continue;
            slic3r_cfg += c;
        }
        slic3r_cfg += "\"/>\n";
        slic3r_cfg += "      <metadata type=\"volume\" key=\"volume_type\" value=\"ModelPart\"/>\n";
        std::snprintf(buf, sizeof(buf),
            "      <metadata type=\"volume\" key=\"extruder\" value=\"%d\"/>\n",
            vol.extruder);
        slic3r_cfg += buf;
        slic3r_cfg += "    </volume>\n";
        cur_first += nf;
    }
    slic3r_cfg += "  </object>\n</config>\n";
    return slic3r_cfg;
}

// ---------------------------------------------------------------------------
// write_3mf / write_3mf_bytes
// ---------------------------------------------------------------------------

void write_3mf(const Mesh& mesh, const std::filesystem::path& path) {
    pack_3mf_to_file(build_single_model_xml(mesh), nullptr, path);
}

std::vector<uint8_t> write_3mf_bytes(const Mesh& mesh) {
    return pack_3mf_to_buffer(build_single_model_xml(mesh), nullptr);
}

// ---------------------------------------------------------------------------
// write_3mf_volumes / write_3mf_volumes_bytes — multi-object 3MF
// ---------------------------------------------------------------------------

void write_3mf_volumes(const std::vector<ThreeMfVolume>& volumes,
                       const std::filesystem::path& path,
                       bool write_slic3r_config) {
    if (volumes.empty()) {
        throw ThreeMfError("write_3mf_volumes: no volumes provided");
    }
    const std::string model = build_volumes_model_xml(volumes);
    std::string cfg;
    if (write_slic3r_config) cfg = build_slic3r_cfg(volumes);
    pack_3mf_to_file(model, (write_slic3r_config && !cfg.empty()) ? &cfg : nullptr, path);
}

std::vector<uint8_t> write_3mf_volumes_bytes(const std::vector<ThreeMfVolume>& volumes,
                                             bool write_slic3r_config) {
    if (volumes.empty()) {
        throw ThreeMfError("write_3mf_volumes_bytes: no volumes provided");
    }
    const std::string model = build_volumes_model_xml(volumes);
    std::string cfg;
    if (write_slic3r_config) cfg = build_slic3r_cfg(volumes);
    return pack_3mf_to_buffer(model, (write_slic3r_config && !cfg.empty()) ? &cfg : nullptr);
}

// ---------------------------------------------------------------------------
// read_3mf
// ---------------------------------------------------------------------------

// Parse one <mesh>...</mesh> block (assumed to live inside the substring
// [begin..end)) and return the resulting Mesh. The mesh vertex indices are
// LOCAL to this mesh block; the caller is responsible for any concatenation.
static Mesh parse_3mf_mesh_block(const std::string& xml,
                                 size_t begin, size_t end) {
    // Defensive: only search within [begin, end). The previous implementation
    // searched the full xml string with a post-hoc `> end` reject; on
    // misordered or back-to-back multi-object inputs that pattern could
    // (in principle) match a sibling object's tags. Constraining the
    // upper bound up front eliminates the cross-object aliasing risk.
    if (end > xml.size()) end = xml.size();
    if (begin > end) begin = end;
    const size_t span = end - begin;
    auto find_bounded = [&](const char* needle) -> size_t {
        const size_t p = xml.find(needle, begin);
        if (p == std::string::npos || p >= end) return std::string::npos;
        return p;
    };
    const size_t verts_open  = find_bounded("<vertices>");
    const size_t verts_close = find_bounded("</vertices>");
    const size_t tris_open   = find_bounded("<triangles>");
    const size_t tris_close  = find_bounded("</triangles>");
    if (verts_open  == std::string::npos ||
        verts_close == std::string::npos ||
        tris_open   == std::string::npos ||
        tris_close  == std::string::npos) {
        throw ThreeMfError("missing <vertices>/<triangles> in 3MF <mesh> block");
    }
    (void)span;

    // Caps on count to bound parser work and downstream allocations even
    // when the model XML is within the 256 MB zip-bomb extraction cap.
    // 100M each is generous (every meshseal-known fixture is <2M).
    constexpr size_t kMaxParseVerts = 100'000'000u;
    constexpr size_t kMaxParseTris  = 100'000'000u;

    Mesh mesh;
    // <vertex>
    {
        size_t pos = verts_open + std::strlen("<vertices>");
        while (pos < verts_close) {
            const size_t tag_start = xml.find("<vertex", pos);
            if (tag_start == std::string::npos || tag_start >= verts_close) break;
            const size_t tag_end = xml.find('>', tag_start);
            if (tag_end == std::string::npos) break;
            const std::string tag = xml.substr(tag_start, tag_end - tag_start + 1);
            const std::string sx = get_attr(tag, "x");
            const std::string sy = get_attr(tag, "y");
            const std::string sz = get_attr(tag, "z");
            if (sx.empty() || sy.empty() || sz.empty())
                throw ThreeMfError("vertex tag missing x/y/z attribute");
            double dx = 0, dy = 0, dz = 0;
            if (!parse_double(sx, dx) || !parse_double(sy, dy) || !parse_double(sz, dz))
                throw ThreeMfError("invalid vertex coordinate (non-numeric)");
            if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz))
                throw ThreeMfError("NaN/Inf vertex coordinate");
            if (mesh.vertices.size() >= kMaxParseVerts)
                throw ThreeMfError("3MF vertex count exceeds 100M cap");
            mesh.vertices.push_back({ dx, dy, dz });
            pos = tag_end + 1;
        }
    }
    // <triangle>
    {
        size_t pos = tris_open + std::strlen("<triangles>");
        while (pos < tris_close) {
            const size_t tag_start = xml.find("<triangle", pos);
            if (tag_start == std::string::npos || tag_start >= tris_close) break;
            const size_t tag_end = xml.find('>', tag_start);
            if (tag_end == std::string::npos) break;
            const std::string tag = xml.substr(tag_start, tag_end - tag_start + 1);
            const std::string sv1 = get_attr(tag, "v1");
            const std::string sv2 = get_attr(tag, "v2");
            const std::string sv3 = get_attr(tag, "v3");
            if (sv1.empty() || sv2.empty() || sv3.empty())
                throw ThreeMfError("triangle tag missing v1/v2/v3 attribute");
            uint32_t i1 = 0, i2 = 0, i3 = 0;
            if (!parse_uint32(sv1, i1) || !parse_uint32(sv2, i2) || !parse_uint32(sv3, i3))
                throw ThreeMfError("invalid triangle index (non-numeric or out-of-range)");
            const uint32_t vcount = static_cast<uint32_t>(mesh.vertices.size());
            if (i1 >= vcount || i2 >= vcount || i3 >= vcount)
                throw ThreeMfError("triangle index out of range");
            if (mesh.faces.size() >= kMaxParseTris)
                throw ThreeMfError("3MF triangle count exceeds 100M cap");
            mesh.faces.push_back({ i1, i2, i3 });
            pos = tag_end + 1;
        }
    }
    return mesh;
}

// Locate each <object id="N" type="model"> ... </object> block in the 3MF
// XML, parse its <mesh>, and return one (id, name, mesh) tuple per object.
// IDs and names are LOCAL to the XML (caller is responsible for any
// remapping). For a 3MF with a single <object>, returns one entry. For
// a multi-object 3MF (as produced by write_3mf_volumes), returns one
// entry per object.
struct ParsedObject {
    int         id = 0;
    std::string name;
    Mesh        mesh;
};
static std::vector<ParsedObject> parse_3mf_objects(const std::string& xml) {
    std::vector<ParsedObject> out;
    size_t pos = 0;
    while (pos < xml.size()) {
        const size_t obj_open = xml.find("<object", pos);
        if (obj_open == std::string::npos) break;
        const size_t obj_tag_end = xml.find('>', obj_open);
        if (obj_tag_end == std::string::npos) break;
        const size_t obj_close = xml.find("</object>", obj_tag_end);
        if (obj_close == std::string::npos) break;
        const std::string tag = xml.substr(obj_open, obj_tag_end - obj_open + 1);
        // Skip object refs (type="other" or components-only); we only want
        // objects with a <mesh> inside.
        const size_t mesh_open = xml.find("<mesh>", obj_tag_end);
        if (mesh_open != std::string::npos && mesh_open < obj_close) {
            ParsedObject po;
            const std::string sid = get_attr(tag, "id");
            if (!sid.empty()) {
                int idv = 0;
                if (parse_int(sid, idv)) po.id = idv;
            }
            po.name = get_attr(tag, "name");
            po.mesh = parse_3mf_mesh_block(xml, obj_tag_end, obj_close);
            out.push_back(std::move(po));
        }
        pos = obj_close + 1;
    }
    if (out.empty()) {
        throw ThreeMfError("no <object> blocks with <mesh> in 3MF model");
    }
    return out;
}

// Concatenate all parsed objects into a single mesh, offsetting vertex
// indices. Shared by read_3mf / read_3mf_bytes.
static Mesh concat_objects(const std::vector<ParsedObject>& objects) {
    Mesh out;
    size_t total_v = 0, total_f = 0;
    for (const auto& o : objects) {
        total_v += o.mesh.vertices.size();
        total_f += o.mesh.faces.size();
    }
    out.vertices.reserve(total_v);
    out.faces.reserve(total_f);
    for (const auto& o : objects) {
        const uint32_t base = static_cast<uint32_t>(out.vertices.size());
        out.vertices.insert(out.vertices.end(),
            o.mesh.vertices.begin(), o.mesh.vertices.end());
        for (const auto& f : o.mesh.faces) {
            out.faces.push_back({ f[0] + base, f[1] + base, f[2] + base });
        }
    }
    return out;
}

// RAII closer for a reader archive so an exception during extraction/parse
// can't leak the miniz state.
namespace {
struct ZipReaderGuard {
    mz_zip_archive* z;
    ~ZipReaderGuard() { if (z) mz_zip_reader_end(z); }
};
} // namespace

// Extract one entry from an initialised reader archive. Returns nullopt if
// the entry is absent. Enforces the uncompressed-size cap (zip-bomb guard);
// throws ThreeMfError on cap violation or extraction failure.
static std::optional<std::string> extract_entry(mz_zip_archive& zip,
                                                 const char* name,
                                                 mz_uint64 cap) {
    int idx = mz_zip_reader_locate_file(&zip, name, nullptr, 0);
    if (idx < 0) return std::nullopt;
    mz_zip_archive_file_stat stat{};
    if (!mz_zip_reader_file_stat(&zip, static_cast<mz_uint>(idx), &stat)) {
        throw ThreeMfError(std::string("cannot stat 3MF entry: ") + name);
    }
    if (stat.m_uncomp_size > cap) {
        throw ThreeMfError(std::string("3MF entry exceeds size cap: ") + name);
    }
    size_t out_size = 0;
    void* raw = mz_zip_reader_extract_to_heap(&zip,
        static_cast<mz_uint>(idx), &out_size, 0);
    if (!raw) {
        throw ThreeMfError(std::string("failed to extract 3MF entry: ") + name);
    }
    std::string s(static_cast<const char*>(raw), out_size);
    mz_free(raw);
    return s;
}

// Pull the model XML (required) and the Slic3r config (optional) out of an
// initialised reader archive.
static void load_3mf_parts(mz_zip_archive& zip,
                           std::string& model_xml,
                           std::string& slic3r_cfg) {
    auto m = extract_entry(zip, "3D/3dmodel.model", kMaxModelXmlUncompressed);
    if (!m) throw ThreeMfError("3D/3dmodel.model not found in archive");
    model_xml = std::move(*m);
    auto c = extract_entry(zip, "Metadata/Slic3r_PE_model.config",
                           kMaxSlicerConfigUncompressed);
    slic3r_cfg = c ? std::move(*c) : std::string();
}

Mesh read_3mf(const std::filesystem::path& path) {
    mz_zip_archive zip{};
    if (!mz_zip_reader_init_file(&zip, path.string().c_str(), 0)) {
        throw ThreeMfError("cannot open 3MF file: " + path.string());
    }
    ZipReaderGuard guard{&zip};
    std::string xml, cfg;
    load_3mf_parts(zip, xml, cfg);
    return concat_objects(parse_3mf_objects(xml));
}

Mesh read_3mf_bytes(const uint8_t* data, size_t size) {
    mz_zip_archive zip{};
    if (!mz_zip_reader_init_mem(&zip, data, size, 0)) {
        throw ThreeMfError("cannot open 3MF buffer");
    }
    ZipReaderGuard guard{&zip};
    std::string xml, cfg;
    load_3mf_parts(zip, xml, cfg);
    return concat_objects(parse_3mf_objects(xml));
}

// ---------------------------------------------------------------------------
// read_3mf_volumes — read 3MF and split into shells
// ---------------------------------------------------------------------------
//
// Two paths handled:
//
//   1) Multi-<object> 3MF (e.g. what write_3mf_volumes produces, or any
//      conforming multi-object 3MF): each <object> becomes its own
//      ThreeMfVolume. The name is the object's `name` attribute (if any).
//      If Metadata/Slic3r_PE_model.config exists with a matching number
//      of <volume> blocks, the names and extruders are taken from there.
//
//   2) Single-<object> 3MF with Slic3r_PE_model.config (PrusaSlicer's
//      packing): everything is in one <object> with the per-shell
//      triangle-index partition in `Metadata/Slic3r_PE_model.config` as
//      `<volume firstid="N" lastid="M">` entries.
//
// For 3MFs without either form, falls back to a single volume covering
// the whole mesh.
//
// volumes_from_xml is the pure-string core, shared by the file and
// in-memory entry points: the model XML and the (possibly empty) slicer
// config string are already extracted from the archive.
static std::vector<ThreeMfVolume> volumes_from_xml(const std::string& xml,
                                                   const std::string& cfg) {
    auto objects = parse_3mf_objects(xml);

    // For the multi-object case, build volumes directly from objects.
    if (objects.size() > 1) {
        std::vector<ThreeMfVolume> out;
        out.reserve(objects.size());
        for (auto& o : objects) {
            ThreeMfVolume tv;
            tv.mesh = std::move(o.mesh);
            tv.name = o.name;
            tv.extruder = 0;
            out.push_back(std::move(tv));
        }
        // Best-effort: if the file also has Slic3r_PE_model.config with the
        // same number of <volume> entries, take names + extruders from
        // there (so write_3mf_volumes round-trip preserves metadata).
        if (!cfg.empty()) {
            std::vector<std::pair<std::string, int>> meta;
            size_t pos = 0;
            while (pos < cfg.size()) {
                const size_t vo = cfg.find("<volume ", pos);
                if (vo == std::string::npos) break;
                const size_t ve = cfg.find("</volume>", vo);
                if (ve == std::string::npos) break;
                const std::string body = cfg.substr(vo, ve - vo);
                std::string nm; int ex = 0;
                size_t mp = 0;
                while (mp < body.size()) {
                    const size_t mo = body.find("<metadata", mp);
                    if (mo == std::string::npos) break;
                    const size_t me = body.find("/>", mo);
                    if (me == std::string::npos) break;
                    const std::string mtag = body.substr(mo, me - mo + 2);
                    const std::string k = get_attr(mtag, "key");
                    const std::string v = get_attr(mtag, "value");
                    if (k == "name") nm = v;
                    else if (k == "extruder") { int ev = 0; if (parse_int(v, ev)) ex = ev; }
                    mp = me + 2;
                }
                meta.emplace_back(nm, ex);
                pos = ve + 1;
            }
            if (meta.size() == out.size()) {
                for (size_t i = 0; i < out.size(); ++i) {
                    if (out[i].name.empty()) out[i].name = meta[i].first;
                    out[i].extruder = meta[i].second;
                }
            }
        }
        return out;
    }

    // Single-object path: concat the (only) object to one Mesh and look
    // for Slic3r_PE_model.config volume ranges within it.
    Mesh full = std::move(objects.front().mesh);
    const uint32_t total_faces = static_cast<uint32_t>(full.faces.size());

    struct VolRange { uint32_t firstid, lastid; std::string name; int extruder; };
    std::vector<VolRange> ranges;

    if (!cfg.empty()) {
        // Parse each <volume firstid="N" lastid="M"> ... </volume>
        size_t pos = 0;
        while (pos < cfg.size()) {
            const size_t vol_open = cfg.find("<volume ", pos);
            if (vol_open == std::string::npos) break;
            const size_t open_end = cfg.find('>', vol_open);
            if (open_end == std::string::npos) break;
            const size_t vol_close = cfg.find("</volume>", open_end);
            if (vol_close == std::string::npos) break;
            const std::string open_tag =
                cfg.substr(vol_open, open_end - vol_open + 1);
            const std::string body =
                cfg.substr(open_end + 1, vol_close - open_end - 1);
            const std::string sf = get_attr(open_tag, "firstid");
            const std::string sl = get_attr(open_tag, "lastid");
            if (!sf.empty() && !sl.empty()) {
                VolRange vr;
                if (!parse_uint32(sf, vr.firstid) ||
                    !parse_uint32(sl, vr.lastid)) {
                    pos = vol_close + 1; continue;
                }
                // name + extruder inside <metadata key="name|extruder" value="...">
                size_t mp = 0;
                while (mp < body.size()) {
                    const size_t mo = body.find("<metadata", mp);
                    if (mo == std::string::npos) break;
                    const size_t me = body.find("/>", mo);
                    if (me == std::string::npos) break;
                    const std::string mtag = body.substr(mo, me - mo + 2);
                    const std::string k = get_attr(mtag, "key");
                    const std::string v = get_attr(mtag, "value");
                    if (k == "name") vr.name = v;
                    else if (k == "extruder") {
                        int ev = 0;
                        if (parse_int(v, ev)) vr.extruder = ev;
                    }
                    mp = me + 2;
                }
                // Sanity: bound to actual face count; skip degenerate
                if (vr.lastid < total_faces && vr.firstid <= vr.lastid)
                    ranges.push_back(std::move(vr));
            }
            pos = vol_close + 1;
        }
    }

    // Fallback: no metadata or no parseable volumes → single volume.
    if (ranges.empty()) {
        ThreeMfVolume tv;
        tv.mesh = std::move(full);
        tv.name = "";
        tv.extruder = 0;
        return { std::move(tv) };
    }

    // Build one Mesh per volume by extracting the face range and
    // re-mapping vertex indices (each volume gets its own vertex pool).
    std::vector<ThreeMfVolume> out;
    out.reserve(ranges.size());
    for (const auto& vr : ranges) {
        ThreeMfVolume tv;
        tv.name = vr.name;
        tv.extruder = vr.extruder;
        std::vector<uint32_t> remap(full.vertices.size(),
                                    std::numeric_limits<uint32_t>::max());
        tv.mesh.faces.reserve(vr.lastid - vr.firstid + 1);
        for (uint32_t fi = vr.firstid; fi <= vr.lastid; ++fi) {
            const auto& src = full.faces[fi];
            std::array<uint32_t, 3> nf;
            for (int k = 0; k < 3; ++k) {
                uint32_t sv = src[k];
                if (remap[sv] == std::numeric_limits<uint32_t>::max()) {
                    remap[sv] = static_cast<uint32_t>(tv.mesh.vertices.size());
                    tv.mesh.vertices.push_back(full.vertices[sv]);
                }
                nf[k] = remap[sv];
            }
            tv.mesh.faces.push_back(nf);
        }
        out.push_back(std::move(tv));
    }
    return out;
}

std::vector<ThreeMfVolume> read_3mf_volumes(const std::filesystem::path& path) {
    mz_zip_archive zip{};
    if (!mz_zip_reader_init_file(&zip, path.string().c_str(), 0)) {
        throw ThreeMfError("cannot open 3MF file: " + path.string());
    }
    ZipReaderGuard guard{&zip};
    std::string xml, cfg;
    load_3mf_parts(zip, xml, cfg);
    return volumes_from_xml(xml, cfg);
}

std::vector<ThreeMfVolume> read_3mf_volumes_bytes(const uint8_t* data, size_t size) {
    mz_zip_archive zip{};
    if (!mz_zip_reader_init_mem(&zip, data, size, 0)) {
        throw ThreeMfError("cannot open 3MF buffer");
    }
    ZipReaderGuard guard{&zip};
    std::string xml, cfg;
    load_3mf_parts(zip, xml, cfg);
    return volumes_from_xml(xml, cfg);
}

} // namespace meshseal
