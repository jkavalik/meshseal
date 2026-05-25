#pragma once
#include "../include/meshseal/meshseal.h"
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace meshseal {

struct ThreeMfError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Read a minimal 3MF file (single mesh object, no materials).
// Preserves full double-precision coordinates.
// Throws ThreeMfError on malformed input.
Mesh read_3mf(const std::filesystem::path& path);

// One named volume parsed from a 3MF file. Multi-color / multi-part 3MFs
// from PrusaSlicer (and others) pack several shells into a single
// <object> in the main model and put the per-shell triangle-index
// partition into `Metadata/Slic3r_PE_model.config` as `<volume firstid="N"
// lastid="M">` entries. read_3mf_volumes returns one entry per volume,
// each carrying the sub-mesh extracted from the indicated face range.
// For 3MFs without volume metadata, returns a single entry covering the
// whole mesh.
struct ThreeMfVolume {
    Mesh        mesh;       // sub-mesh for this volume (own vertex pool, faces re-indexed)
    std::string name;       // volume name from metadata (e.g. "black.stl", "Bee_v3.stl_1_4")
    int         extruder = 0;
};
std::vector<ThreeMfVolume> read_3mf_volumes(const std::filesystem::path& path);

// Write mesh as minimal 3MF. Coordinates written as full-precision decimal strings.
void write_3mf(const Mesh& mesh, const std::filesystem::path& path);

// Write a multi-volume 3MF: separate <object> per ThreeMfVolume, each
// with its own <vertices>/<triangles>. The <build> section lists all
// objects so a slicer sees them as physically distinct bodies.
//
// This is the lossless multi-volume output format — preserves the per-
// shell separation the way PrusaSlicer's "Fix by Windows repair" round-
// trip does. STL output of multi-volume meshes is lossy (positional
// welding on re-read merges shared boundaries); this is not.
//
// When write_slic3r_config is true, also writes
// Metadata/Slic3r_PE_model.config so re-reading the file with
// read_3mf_volumes recovers the same partition. The triangle ranges
// (firstid/lastid) are computed by concatenating the volumes in order.
void write_3mf_volumes(const std::vector<ThreeMfVolume>& volumes,
                       const std::filesystem::path& path,
                       bool write_slic3r_config = true);

} // namespace meshseal
