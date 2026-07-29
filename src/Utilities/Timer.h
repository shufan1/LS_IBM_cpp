#pragma once
#include <chrono>
#include <map>
#include <string>
#include <vector>

// Mirrors the MATLAB `global Timer` struct-of-arrays pattern used
// throughout beta_AD_correctu (SolveUVP.m, SolveTransportADRE.m,
// modelSimulation.m, etc.): one named category per measured code
// section, each holding one elapsed-seconds sample per call. Category
// names should match their MATLAB counterparts where there's a direct
// equivalent (e.g. "diff_flux_u"/"diff_flux_v"), so C++ vs. MATLAB
// timings can be compared directly.
//
// Not rank-aware: each MPI rank has its own independent registry (this
// is a plain process-local static, same as MATLAB's `global`). Callers
// are responsible for only printing/saving from rank 0 if that's what's
// wanted (see main.cpp).
class Timer {
public:
    // Appends one elapsed-seconds sample to `category`'s array. The
    // category is created on first use -- no separate registration step
    // needed, matching MATLAB's `if ~isfield(...)` guard.
    static void record(const std::string &category, double seconds);

    // All recorded categories and their raw per-call sample arrays --
    // for saving/plotting/comparing against MATLAB's timer_results.mat.
    static const std::map<std::string, std::vector<double>> &all();

    // Prints a category/calls/total(s)/mean(s) table, matching
    // RunADRE.m's TIMING SUMMARY block.
    static void printSummary();

    // Dumps all() to a JSON file (one array of per-call samples per
    // category), for direct comparison against MATLAB's
    // timer_results.mat via the same plot_timer_results.m-style
    // tooling. Throws std::runtime_error if the file can't be written.
    static void saveJson(const std::string &path);
};

// RAII scoped timer: construct at the top of a block you want timed; it
// records the elapsed time into `category` when it goes out of scope
// (including on early return/exception), mirroring MATLAB's
// `t=tic; ...; Timer.category(end+1)=toc(t);` pattern without needing a
// matching manual "stop" call to remember.
//
// Usage: { ScopedTimer t("diff_flux_u"); ... code to time ...; }
class ScopedTimer {
public:
    explicit ScopedTimer(std::string category)
        : category_(std::move(category)), start_(std::chrono::steady_clock::now()) {}
    ~ScopedTimer();

    ScopedTimer(const ScopedTimer &) = delete;
    ScopedTimer &operator=(const ScopedTimer &) = delete;

private:
    std::string category_;
    std::chrono::steady_clock::time_point start_;
};
