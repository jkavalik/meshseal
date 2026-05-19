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

// Write mesh as minimal 3MF. Coordinates written as full-precision decimal strings.
void write_3mf(const Mesh& mesh, const std::filesystem::path& path);

} // namespace meshseal
