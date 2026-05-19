#pragma once

#include "meshseal.h"
#include <string>
#include <vector>

namespace meshseal {

class Pipeline {
public:
    Pipeline& weld();
    Pipeline& degenerate();
    Pipeline& orient();
    Pipeline& non_manifold();
    Pipeline& holes();
    Pipeline& shells();
    Pipeline& intersections();
    Pipeline& thin_features();
    Pipeline& soup_reconstruct();

    RepairResult run(const Mesh& mesh, const RepairOptions& opts = RepairOptions{}) const;

private:
    std::vector<std::string> stages_;
};

} // namespace meshseal
