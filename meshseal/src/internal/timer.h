#pragma once
#include <chrono>
#include <string>
#include <map>

namespace meshseal::internal {

// RAII scoped timer. On destruction, adds elapsed ms to the given map.
class ScopedTimer {
public:
    ScopedTimer(const std::string& name, std::map<std::string, double>& times)
        : name_(name), times_(times),
          start_(std::chrono::steady_clock::now()) {}

    ~ScopedTimer() {
        auto end = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(end - start_).count();
        times_[name_] += ms;  // += so repeated runs accumulate
    }

private:
    std::string name_;
    std::map<std::string, double>& times_;
    std::chrono::time_point<std::chrono::steady_clock> start_;
};

} // namespace meshseal::internal
