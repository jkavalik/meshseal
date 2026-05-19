#pragma once
#include "../include/meshseal/meshseal.h"
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace meshseal {

struct StlError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Read binary STL from file. Throws StlError on malformed input.
Mesh read_stl(const std::filesystem::path& path);

// Read binary STL from memory buffer. Throws StlError on malformed input.
Mesh read_stl_bytes(const uint8_t* data, size_t size);

// Write mesh as binary STL. Face normals computed from vertex order.
// The Mesh stores double-precision vertices; STL stores float32.
// Precision loss is expected and documented — use 3MF for lossless round-trips.
void write_stl(const Mesh& mesh, const std::filesystem::path& path);

// Write mesh as binary STL to a byte buffer (returned).
std::vector<uint8_t> write_stl_bytes(const Mesh& mesh);

} // namespace meshseal
