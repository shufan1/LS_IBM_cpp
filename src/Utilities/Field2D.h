#pragma once
#include <vector>

// A 2D view over a flat, contiguous std::vector<double> -- gives U(i,j)
// syntax without giving up contiguous storage (needed for cache
// locality, and for eventually copying into/out of PETSc Vecs/Mats).
// Row-major: linear index = i*ny + j. No bounds checking, by design --
// this is meant for hot assembly loops (COEFFU/COEFFV/COEFFP-style code)
// where per-access checks would add up.
class Field2D {
public:
    Field2D() = default;
    Field2D(int nx, int ny, double fill = 0.0)
        : nx_(nx), ny_(ny), data_(static_cast<size_t>(nx) * ny, fill) {}

    double &operator()(int i, int j) { return data_[static_cast<size_t>(i) * ny_ + j]; }
    double operator()(int i, int j) const { return data_[static_cast<size_t>(i) * ny_ + j]; }

    int nx() const { return nx_; }
    int ny() const { return ny_; }

    // Escape hatch to the raw contiguous buffer -- for PETSc
    // MatSetValues/VecSetValues, json-c array serialization, etc.
    std::vector<double> &data() { return data_; }
    const std::vector<double> &data() const { return data_; }

private:
    int nx_ = 0, ny_ = 0;
    std::vector<double> data_;
};
