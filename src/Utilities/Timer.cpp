#include "Timer.h"
#include <json-c/json.h>
#include <cstdio>
#include <numeric>
#include <stdexcept>

namespace {
// Function-local static (Meyer's singleton) instead of a plain global --
// avoids static-initialization-order issues across translation units
// while still giving every caller in this process the same registry.
std::map<std::string, std::vector<double>> &registry() {
    static std::map<std::string, std::vector<double>> instance;
    return instance;
}
}  // namespace

void Timer::record(const std::string &category, double seconds) {
    registry()[category].push_back(seconds);
}

const std::map<std::string, std::vector<double>> &Timer::all() {
    return registry();
}

void Timer::printSummary() {
    printf("=========================================================================== \n");
    printf("                      TIMING SUMMARY (seconds)                       \n");
    printf("=========================================================================== \n");
    printf("%-28s %10s %12s %12s\n", "category", "calls", "total(s)", "mean(s)");
    for (const auto &entry : registry()) {
        const std::string &name = entry.first;
        const std::vector<double> &samples = entry.second;
        double total = std::accumulate(samples.begin(), samples.end(), 0.0);
        double mean = samples.empty() ? 0.0 : total / static_cast<double>(samples.size());
        printf("%-28s %10zu %12.4f %12.6f\n", name.c_str(), samples.size(), total, mean);
    }
}

void Timer::saveJson(const std::string &path) {
    json_object *root = json_object_new_object();
    for (const auto &entry : registry()) {
        json_object *arr = json_object_new_array();
        for (double s : entry.second) json_object_array_add(arr, json_object_new_double(s));
        json_object_object_add(root, entry.first.c_str(), arr);
    }

    if (json_object_to_file_ext(path.c_str(), root, JSON_C_TO_STRING_PRETTY) != 0) {
        json_object_put(root);
        throw std::runtime_error("Timer::saveJson: failed to write " + path);
    }
    json_object_put(root);
}

ScopedTimer::~ScopedTimer() {
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
    Timer::record(category_, elapsed);
}
