#include "../include/meshseal/pipeline.h"
#include "../include/meshseal/meshseal.h"
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace meshseal {

// Stage prerequisite table
static const std::vector<std::pair<std::string, std::vector<std::string>>> k_prereqs = {
    {"degenerate",      {"weld"}},
    {"orient",          {"degenerate"}},
    {"non_manifold",    {"orient"}},
    {"holes",           {"non_manifold"}},
    {"shells",          {"holes"}},
    {"intersections",   {"shells"}},
    {"thin_features",   {"shells"}},
    {"soup_reconstruct",{"intersections"}},
};

Pipeline& Pipeline::weld()            { stages_.push_back("weld");            return *this; }
Pipeline& Pipeline::degenerate()      { stages_.push_back("degenerate");      return *this; }
Pipeline& Pipeline::orient()          { stages_.push_back("orient");          return *this; }
Pipeline& Pipeline::non_manifold()    { stages_.push_back("non_manifold");    return *this; }
Pipeline& Pipeline::holes()           { stages_.push_back("holes");           return *this; }
Pipeline& Pipeline::shells()          { stages_.push_back("shells");          return *this; }
Pipeline& Pipeline::intersections()   { stages_.push_back("intersections");   return *this; }
Pipeline& Pipeline::thin_features()   { stages_.push_back("thin_features");   return *this; }
Pipeline& Pipeline::soup_reconstruct(){ stages_.push_back("soup_reconstruct");return *this; }

RepairResult Pipeline::run(const Mesh& mesh, const RepairOptions& base_opts) const {
    // Build a custom RepairOptions with only the selected stages enabled.
    RepairOptions opts;
    opts.weld_tolerance   = base_opts.weld_tolerance;
    opts.weld             = false;
    opts.degenerate       = false;
    opts.orient           = false;
    opts.non_manifold     = false;
    opts.holes            = false;
    opts.shells           = false;
    opts.intersections    = false;
    opts.thin_features    = false;
    opts.soup_reconstruct = false;

    RepairResult result;

    for (const auto& stage : stages_) {
        // Warn about missing prerequisites via notes
        for (const auto& [s, prereqs] : k_prereqs) {
            if (s != stage) continue;
            for (const auto& req : prereqs) {
                if (std::find(stages_.begin(), stages_.end(), req) == stages_.end()) {
                    result.notes.push_back(
                        "warning: stage '" + stage + "' declared prerequisite '" +
                        req + "' is not in this pipeline");
                }
            }
        }

        if (stage == "weld")             opts.weld             = true;
        else if (stage == "degenerate")  opts.degenerate       = true;
        else if (stage == "orient")      opts.orient           = true;
        else if (stage == "non_manifold")opts.non_manifold     = true;
        else if (stage == "holes")       opts.holes            = true;
        else if (stage == "shells")      opts.shells           = true;
        else if (stage == "intersections")opts.intersections   = true;
        else if (stage == "thin_features")opts.thin_features   = true;
        else if (stage == "soup_reconstruct")opts.soup_reconstruct = true;
    }

    // Delegate to the standard repair() with the customised opts
    auto r = repair(mesh, opts);
    // Preserve any prerequisite warnings we added
    r.notes.insert(r.notes.begin(), result.notes.begin(), result.notes.end());
    return r;
}

} // namespace meshseal
