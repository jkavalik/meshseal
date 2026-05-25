#include <meshseal/meshseal.h>
#include "../meshseal/src/stl_io.h"
#include "../meshseal/src/3mf_io.h"
#include <cctype>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

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

    auto str_lower = [](std::string s) {
        for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    };
    const std::string in_ext  = str_lower(in_path.extension().string());
    const std::string out_ext = str_lower(out_path.extension().string());

    // Multi-volume 3MF → 3MF: parse all volumes, repair each independently,
    // write a multi-<object> 3MF preserving the per-volume separation.
    // This matches PrusaSlicer's "Fix by Windows repair algorithm"
    // round-trip semantics (each ModelVolume repaired in isolation, written
    // back into its own 3MF object). Other input/output combos fall through
    // to the single-mesh path.
    const bool in_3mf  = (in_ext  == ".3mf");
    const bool out_3mf = (out_ext == ".3mf");

    if (in_3mf && out_3mf) {
        std::vector<meshseal::ThreeMfVolume> volumes;
        try {
            volumes = meshseal::read_3mf_volumes(in_path);
        } catch (const std::exception& e) {
            std::cerr << "Error reading input: " << e.what() << "\n";
            return 1;
        }
        if (volumes.size() > 1) {
            if (!quiet) {
                std::cout << "3MF: " << volumes.size()
                          << " volumes detected — per-volume repair\n";
            }
            std::vector<meshseal::ThreeMfVolume> out_vols;
            out_vols.reserve(volumes.size());
            bool any_partial = false;
            for (size_t i = 0; i < volumes.size(); ++i) {
                if (!quiet) {
                    std::cout << "  [volume " << i << "] name=" << volumes[i].name
                              << " extruder=" << volumes[i].extruder
                              << " F=" << volumes[i].mesh.faces.size() << "\n";
                }
                auto sub = meshseal::repair(volumes[i].mesh);
                if (sub.partial_failure) any_partial = true;
                if (!quiet) {
                    for (const auto& n : sub.notes)
                        std::cout << "    " << n << "\n";
                }
                meshseal::ThreeMfVolume tv;
                tv.mesh = std::move(sub.mesh);
                tv.name = volumes[i].name;
                tv.extruder = volumes[i].extruder;
                if (!tv.mesh.faces.empty())
                    out_vols.push_back(std::move(tv));
            }
            if (any_partial) std::cerr << "warning: repair partially failed\n";
            try {
                meshseal::write_3mf_volumes(out_vols, out_path);
            } catch (const std::exception& e) {
                std::cerr << "Error writing output: " << e.what() << "\n";
                return 1;
            }
            return 0;
        }
        // 3MF with no Slic3r volume metadata — fall through to single mesh.
    }

    meshseal::Mesh mesh;
    try {
        if (in_ext == ".stl") {
            mesh = meshseal::read_stl(in_path);
        } else if (in_ext == ".3mf") {
            mesh = meshseal::read_3mf(in_path);
        } else {
            std::cerr << "Unsupported input format: " << in_ext << "\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error reading input: " << e.what() << "\n";
        return 1;
    }

    meshseal::RepairResult result = meshseal::repair(mesh);

    if (result.partial_failure) {
        std::cerr << "warning: repair partially failed\n";
    }

    if (!quiet) {
        std::cout << "watertight=" << (result.watertight ? "true" : "false")
                  << " is_volume=" << (result.is_volume ? "true" : "false")
                  << " components=" << result.component_count
                  << " self_intersections=" << result.self_intersections
                  << " confidence=" << result.confidence << "\n";
        for (const auto& note : result.notes) {
            std::cout << "  " << note << "\n";
        }
    }

    try {
        if (out_ext == ".stl") {
            meshseal::write_stl(result.mesh, out_path);
        } else if (out_ext == ".3mf") {
            meshseal::write_3mf(result.mesh, out_path);
        } else {
            std::cerr << "Unsupported output format: " << out_ext << "\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error writing output: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
