#include "Interp.h"
#include <algorithm>

namespace {
// Largest index i in [0, n-2] with a[i] <= q (so i and i+1 both index a
// valid cell). Binary search (O(log n)) -- the same O(log n)-not-O(n)
// reasoning as the update_A1g.m ghost-point fix, since this is called
// once per U-/V-grid node. Clamps to the end cells for out-of-range q.
int bracketIndex(const std::vector<double> &a, double q) {
    int n = static_cast<int>(a.size());
    if (q <= a[0]) return 0;
    if (q >= a[n - 1]) return n - 2;
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (a[mid] <= q) lo = mid;
        else hi = mid;
    }
    return lo;
}
}  // namespace

double bilinearInterp(const std::vector<double> &x, const std::vector<double> &y,
                      const Field2D &V, double xq, double yq) {
    int i = bracketIndex(x, xq);
    int j = bracketIndex(y, yq);

    double tx = (xq - x[i]) / (x[i + 1] - x[i]);
    double ty = (yq - y[j]) / (y[j + 1] - y[j]);
    // Constant extrapolation outside the grid (see Interp.h) -- clamp the
    // blend weights rather than letting a query past the edge produce a
    // negative/over-unity weight.
    tx = std::min(1.0, std::max(0.0, tx));
    ty = std::min(1.0, std::max(0.0, ty));

    double f00 = V(i, j);
    double f10 = V(i + 1, j);
    double f01 = V(i, j + 1);
    double f11 = V(i + 1, j + 1);

    return f00 * (1.0 - tx) * (1.0 - ty) + f10 * tx * (1.0 - ty) +
           f01 * (1.0 - tx) * ty + f11 * tx * ty;
}
