#include "mwpm/algorithms.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace mwpm {

namespace {

// Standard interval DP for non-crossing matchings on points lying in convex
// position. Let cost[i][j] be the minimum total weight of a perfect matching
// of points[i..j] (inclusive). The recurrence is
//
//   cost[i][j] = min over k in {i+1, i+3, ..., j} of
//                  d(points[i], points[k]) + cost[i+1..k-1] + cost[k+1..j]
//
// (with k - i odd, so that both sub-ranges have even length). Empty ranges
// contribute cost 0.
//
// Optimisations relative to a textbook implementation:
//
//   * `cost` and `partner` are stored as single contiguous N x N arrays (one
//     allocation each) instead of nested `vector<vector<>>`. This trims a
//     factor of ~N allocations and gives the inner loop a predictable stride.
//   * The Euclidean distance d(points[i], points[k]) appears in the inner
//     loop indexed only by (i, k), not j. The naive nested-loop version
//     therefore recomputes the same sqrt O(N) times per (i, k) pair, i.e.
//     O(N^3) sqrt calls total. We precompute a triangular distance table
//     once -- O(N^2) sqrt calls -- and then the triple loop becomes pure
//     adds + compares. For N = 1000 that drops the sqrt count from ~5 * 10^8
//     to ~5 * 10^5; DP becomes roughly an order of magnitude faster.
//   * `partner` only needs to discriminate among indices 0..N-1; we keep it
//     as `int` (cache-friendlier than the previous `vector<vector<int>>`).

}  // namespace

MatchingResult dp_matching(const std::vector<Point>& points) {
    const int n = static_cast<int>(points.size());
    if (n % 2 != 0) {
        throw std::invalid_argument("dp_matching requires an even number of points");
    }
    MatchingResult result;
    if (n == 0) return result;

    const std::size_t nn = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);

    // Single flat allocation each; cost(i, j) lives at i * n + j.
    std::vector<double> cost(nn, 0.0);
    std::vector<int> partner(nn, -1);

    // Distance table, same flat layout. Only the (i, k) entries with k > i
    // and k - i odd are ever read from this table; we still fill the full
    // upper triangle to keep indexing branch-free (writes to k - i even
    // entries are essentially free in the streaming write pattern below).
    std::vector<double> dist(nn, 0.0);
    for (int i = 0; i < n; ++i) {
        const double xi = points[i].x;
        const double yi = points[i].y;
        const std::size_t row_base = static_cast<std::size_t>(i) * static_cast<std::size_t>(n);
        for (int j = i + 1; j < n; ++j) {
            const double dx = points[j].x - xi;
            const double dy = points[j].y - yi;
            dist[row_base + static_cast<std::size_t>(j)] = std::sqrt(dx * dx + dy * dy);
        }
    }

    // Fill by interval length. Lengths 2, 4, 6, ... up to n.
    double* const cost_data = cost.data();
    int* const partner_data = partner.data();
    const double* const dist_data = dist.data();

    for (int len = 2; len <= n; len += 2) {
        for (int i = 0; i + len - 1 < n; ++i) {
            const int j = i + len - 1;
            const std::size_t row_i = static_cast<std::size_t>(i) * static_cast<std::size_t>(n);
            const std::size_t row_i1 = row_i + static_cast<std::size_t>(n);  // (i+1) * n

            double best = std::numeric_limits<double>::infinity();
            int best_k = -1;
            // k = i + 1: left sub-range is empty (k - 1 < i + 1).
            {
                const int k = i + 1;
                const double right = (j >= k + 1)
                    ? cost_data[(static_cast<std::size_t>(k) + 1) *
                                static_cast<std::size_t>(n) + static_cast<std::size_t>(j)]
                    : 0.0;
                const double total = dist_data[row_i + static_cast<std::size_t>(k)] + right;
                if (total < best) {
                    best = total;
                    best_k = k;
                }
            }
            // k = i + 3, i + 5, ..., j. Both sub-ranges non-empty.
            for (int k = i + 3; k <= j; k += 2) {
                const double left = cost_data[row_i1 + static_cast<std::size_t>(k - 1)];
                const double right = (j >= k + 1)
                    ? cost_data[(static_cast<std::size_t>(k) + 1) *
                                static_cast<std::size_t>(n) + static_cast<std::size_t>(j)]
                    : 0.0;
                const double total =
                    dist_data[row_i + static_cast<std::size_t>(k)] + left + right;
                if (total < best) {
                    best = total;
                    best_k = k;
                }
            }
            cost_data[row_i + static_cast<std::size_t>(j)] = best;
            partner_data[row_i + static_cast<std::size_t>(j)] = best_k;
        }
    }

    result.cost = cost_data[static_cast<std::size_t>(n - 1)];  // cost(0, n-1)

    // Reconstruct pairs by walking the partner table.
    std::vector<Pair> pairs;
    pairs.reserve(static_cast<std::size_t>(n / 2));

    std::vector<std::pair<int, int>> stack;
    stack.reserve(static_cast<std::size_t>(n / 2));
    stack.emplace_back(0, n - 1);
    while (!stack.empty()) {
        auto [i, j] = stack.back();
        stack.pop_back();
        if (i > j) continue;
        const int k = partner_data[static_cast<std::size_t>(i) *
                                   static_cast<std::size_t>(n) +
                                   static_cast<std::size_t>(j)];
        pairs.emplace_back(static_cast<std::size_t>(i), static_cast<std::size_t>(k));
        if (i + 1 <= k - 1) stack.emplace_back(i + 1, k - 1);
        if (k + 1 <= j) stack.emplace_back(k + 1, j);
    }
    result.pairs = std::move(pairs);
    return result;
}

}  // namespace mwpm
