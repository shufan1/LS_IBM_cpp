#include "Timer.h"
#include <json-c/json.h>
#include <cstdio>
#include <limits>
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

    // total_runtime is excluded from the table and reported below instead,
    // matching RunADRE_no_LS.m:117 -- it wraps all the others, so summing
    // it in with them would double-count.
    double instrumentedTotal = 0.0;
    for (const auto &entry : registry()) {
        const std::string &name = entry.first;
        if (name == "total_runtime") continue;
        const std::vector<double> &samples = entry.second;
        double total = std::accumulate(samples.begin(), samples.end(), 0.0);
        // MATLAB's mean([]) is NaN, and it prints that for never-called
        // categories (e.g. solve_p_piso with PISO off). Matched here so
        // the two tables diff cleanly.
        double mean = samples.empty() ? std::numeric_limits<double>::quiet_NaN()
                                      : total / static_cast<double>(samples.size());
        printf("%-28s %10zu %12.4f %12.6f\n", name.c_str(), samples.size(), total, mean);
        instrumentedTotal += total;
    }

    auto it = registry().find("total_runtime");
    if (it != registry().end()) {
        const std::vector<double> &samples = it->second;
        double total = std::accumulate(samples.begin(), samples.end(), 0.0);
        printf("--------------------------------------------------------------------------- \n");
        printf("%-28s %10zu %12.4f\n", "total_runtime", samples.size(), total);
        if (total > 0.0)
            printf("instrumented categories cover %.1f%% of total_runtime\n",
                   100.0 * instrumentedTotal / total);
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
