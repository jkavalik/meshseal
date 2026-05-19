#include <meshseal/meshseal.h>
#include "../meshseal/src/stl_io.h"
#include "../meshseal/src/3mf_io.h"
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: meshseal_cli <input> <output> [--quiet]\n";
        return 1;
    }

    const std::filesystem::path in_path(argv[1]);
    const std::filesystem::path out_path(argv[2]);
    bool quiet = false;
    for (int i = 3; i < argc; ++i) {
        if (std::string(argv[i]) == "--quiet") quiet = true;
    }

    meshseal::Mesh mesh;
    try {
        const auto ext = in_path.extension().string();
        if (ext == ".stl" || ext == ".STL") {
            mesh = meshseal::read_stl(in_path);
        } else if (ext == ".3mf" || ext == ".3MF") {
            mesh = meshseal::read_3mf(in_path);
        } else {
            std::cerr << "Unsupported input format: " << ext << "\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error reading input: " << e.what() << "\n";
        return 1;
    }

    const auto result = meshseal::repair(mesh);

    if (result.partial_failure) {
        std::cerr << "warning: repair partially failed\n";
    }

    if (!quiet) {
        std::cout << "watertight=" << (result.watertight ? "true" : "false")
                  << " is_volume=" << (result.is_volume ? "true" : "false")
                  << " components=" << result.component_count
                  << " confidence=" << result.confidence << "\n";
        for (const auto& note : result.notes) {
            std::cout << "  " << note << "\n";
        }
    }

    try {
        const auto ext = out_path.extension().string();
        if (ext == ".stl" || ext == ".STL") {
            meshseal::write_stl(result.mesh, out_path);
        } else if (ext == ".3mf" || ext == ".3MF") {
            meshseal::write_3mf(result.mesh, out_path);
        } else {
            std::cerr << "Unsupported output format: " << ext << "\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error writing output: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
