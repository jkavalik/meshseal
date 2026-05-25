#include "3mf_io.h"
#include <miniz.h>
#include <array>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace meshseal {

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

// ---------------------------------------------------------------------------
// write_3mf
// ---------------------------------------------------------------------------

void write_3mf(const Mesh& mesh, const std::filesystem::path& path) {
    // Build the 3dmodel.model XML
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

    // Pack into ZIP
    mz_zip_archive zip{};
    if (!mz_zip_writer_init_file(&zip, path.string().c_str(), 0)) {
        throw ThreeMfError("cannot create 3MF file: " + path.string());
    }

    bool ok = true;

    ok = ok && mz_zip_writer_add_mem(&zip, "[Content_Types].xml",
        k_content_types, std::strlen(k_content_types),
        MZ_DEFAULT_COMPRESSION);

    static const char* k_rels =
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n"
        "  <Relationship Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel\""
        " Target=\"/3D/3dmodel.model\" Id=\"rel0\"/>\n"
        "</Relationships>\n";

    ok = ok && mz_zip_writer_add_mem(&zip, "_rels/.rels",
        k_rels, std::strlen(k_rels),
        MZ_DEFAULT_COMPRESSION);

    ok = ok && mz_zip_writer_add_mem(&zip, "3D/3dmodel.model",
        model.data(), model.size(),
        MZ_DEFAULT_COMPRESSION);

    ok = ok && mz_zip_writer_finalize_archive(&zip);
    mz_zip_writer_end(&zip);

    if (!ok) {
        throw ThreeMfError("error writing 3MF archive: " + path.string());
    }
}

// ---------------------------------------------------------------------------
// read_3mf
// ---------------------------------------------------------------------------

Mesh read_3mf(const std::filesystem::path& path) {
    mz_zip_archive zip{};
    if (!mz_zip_reader_init_file(&zip, path.string().c_str(), 0)) {
        throw ThreeMfError("cannot open 3MF file: " + path.string());
    }

    // Locate 3D/3dmodel.model (case-sensitive per spec)
    int idx = mz_zip_reader_locate_file(&zip, "3D/3dmodel.model", nullptr, 0);
    if (idx < 0) {
        mz_zip_reader_end(&zip);
        throw ThreeMfError("3D/3dmodel.model not found in archive");
    }

    size_t out_size = 0;
    void* raw = mz_zip_reader_extract_to_heap(&zip, static_cast<mz_uint>(idx), &out_size, 0);
    mz_zip_reader_end(&zip);

    if (!raw) {
        throw ThreeMfError("failed to extract 3D/3dmodel.model");
    }

    std::string xml(static_cast<const char*>(raw), out_size);
    mz_free(raw);

    // --- Parse vertices ---
    const size_t verts_open = xml.find("<vertices>");
    const size_t verts_close = xml.find("</vertices>");
    if (verts_open == std::string::npos || verts_close == std::string::npos) {
        throw ThreeMfError("missing <vertices> section in 3MF model");
    }

    const size_t tris_open = xml.find("<triangles>");
    const size_t tris_close = xml.find("</triangles>");
    if (tris_open == std::string::npos || tris_close == std::string::npos) {
        throw ThreeMfError("missing <triangles> section in 3MF model");
    }

    Mesh mesh;

    // Parse <vertex ...> tags
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
            if (sx.empty() || sy.empty() || sz.empty()) {
                throw ThreeMfError("vertex tag missing x/y/z attribute");
            }
            try {
                mesh.vertices.push_back({
                    std::stod(sx),
                    std::stod(sy),
                    std::stod(sz)
                });
            } catch (const std::exception& e) {
                throw ThreeMfError(std::string("invalid vertex coordinate: ") + e.what());
            }
            pos = tag_end + 1;
        }
    }

    // Parse <triangle ...> tags
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
            if (sv1.empty() || sv2.empty() || sv3.empty()) {
                throw ThreeMfError("triangle tag missing v1/v2/v3 attribute");
            }
            uint32_t i1, i2, i3;
            try {
                i1 = static_cast<uint32_t>(std::stoul(sv1));
                i2 = static_cast<uint32_t>(std::stoul(sv2));
                i3 = static_cast<uint32_t>(std::stoul(sv3));
            } catch (const std::exception& e) {
                throw ThreeMfError(std::string("invalid triangle index: ") + e.what());
            }
            const uint32_t vcount = static_cast<uint32_t>(mesh.vertices.size());
            if (i1 >= vcount || i2 >= vcount || i3 >= vcount) {
                throw ThreeMfError("triangle index out of range");
            }
            mesh.faces.push_back({i1, i2, i3});
            pos = tag_end + 1;
        }
    }

    return mesh;
}

// ---------------------------------------------------------------------------
// read_3mf_volumes — read 3MF and split into shells using Slic3r metadata
// ---------------------------------------------------------------------------
//
// Multi-color / multi-part 3MFs from PrusaSlicer pack several physical
// shells into one <object> in the main model and put the per-shell
// triangle-index partition into Metadata/Slic3r_PE_model.config (a
// Slic3r-specific extension). Reading this metadata gives a CLEAN shell
// partition without any geometric guessing — the slicer already
// identified the shells when the user authored the multi-part print.
//
// For 3MFs without Slic3r_PE_model.config or without <volume> blocks,
// this falls back to a single volume covering the whole mesh.
std::vector<ThreeMfVolume> read_3mf_volumes(const std::filesystem::path& path) {
    Mesh full = read_3mf(path);
    const uint32_t total_faces = static_cast<uint32_t>(full.faces.size());

    // Re-open the archive to look for Slic3r_PE_model.config.
    struct VolRange { uint32_t firstid, lastid; std::string name; int extruder; };
    std::vector<VolRange> ranges;

    mz_zip_archive zip{};
    if (mz_zip_reader_init_file(&zip, path.string().c_str(), 0)) {
        int idx = mz_zip_reader_locate_file(&zip,
            "Metadata/Slic3r_PE_model.config", nullptr, 0);
        if (idx >= 0) {
            size_t out_size = 0;
            void* raw = mz_zip_reader_extract_to_heap(&zip,
                static_cast<mz_uint>(idx), &out_size, 0);
            if (raw) {
                std::string cfg(static_cast<const char*>(raw), out_size);
                mz_free(raw);

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
                        try {
                            vr.firstid = static_cast<uint32_t>(std::stoul(sf));
                            vr.lastid  = static_cast<uint32_t>(std::stoul(sl));
                        } catch (...) {
                            pos = vol_close + 1; continue;
                        }
                        // name + extruder inside <metadata key="name|extruder" value="...">
                        size_t mp = 0;
                        while (mp < body.size()) {
                            const size_t mo = body.find("<metadata", mp);
                            if (mo == std::string::npos) break;
                            const size_t me = body.find("/>", mo);
                            if (me == std::string::npos) break;
                            const std::string mtag =
                                body.substr(mo, me - mo + 2);
                            const std::string k = get_attr(mtag, "key");
                            const std::string v = get_attr(mtag, "value");
                            if (k == "name") vr.name = v;
                            else if (k == "extruder") {
                                try { vr.extruder = std::stoi(v); }
                                catch (...) {}
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
        }
        mz_zip_reader_end(&zip);
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

} // namespace meshseal
