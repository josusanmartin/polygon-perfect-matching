#include "mwpm/algorithms.hpp"

#include <limits>
#include <stdexcept>

namespace mwpm {

namespace {

// Standard interval DP for non-crossing matchings on points lying in convex
// position. Let cost[i][j] be the minimum total weight of a perfect matching
// of points[i..j] (inclusive). The recurrence is
//
//   cost[i][j] = min over k in {i+1, i+3, ..., j} of
//                  d(points[i], points[k]) + cost[i+1..k-1] + cost[k+1..j]
//
// (with k - i odd, so that both sub-ranges have even length). For convenience
// we treat empty ranges as having cost 0 and an empty matching.

}  // namespace

MatchingResult dp_matching(const std::vector<Point>& points) {
    const int n = static_cast<int>(points.size());
    if (n % 2 != 0) {
        throw std::invalid_argument("dp_matching requires an even number of points");
    }
    MatchingResult result;
    if (n == 0) return result;

    // cost[i][j] for 0 <= i <= j < n with (j - i + 1) even.
    std::vector<std::vector<double>> cost(n, std::vector<double>(n, 0.0));
    // partner[i][j] = the index k that point i is paired with for the optimal
    // matching of [i..j]. Used to reconstruct the matching pairs.
    std::vector<std::vector<int>> partner(n, std::vector<int>(n, -1));

    // Fill by interval length. Lengths 2, 4, 6, ... up to n.
    for (int len = 2; len <= n; len += 2) {
        for (int i = 0; i + len - 1 < n; ++i) {
            const int j = i + len - 1;
            double best = std::numeric_limits<double>::infinity();
            int best_k = -1;
            for (int k = i + 1; k <= j; k += 2) {
                // Pair i with k, leaving sub-intervals [i+1..k-1] and [k+1..j].
                const double left = (k - 1 >= i + 1) ? cost[i + 1][k - 1] : 0.0;
                const double right = (j >= k + 1) ? cost[k + 1][j] : 0.0;
                const double total = distance(points[i], points[k]) + left + right;
                if (total < best) {
                    best = total;
                    best_k = k;
                }
            }
            cost[i][j] = best;
            partner[i][j] = best_k;
        }
    }

    result.cost = cost[0][n - 1];

    // Reconstruct pairs by walking the partner table.
    std::vector<Pair> pairs;
    pairs.reserve(static_cast<std::size_t>(n / 2));

    // Stack-based traversal over intervals.
    std::vector<std::pair<int, int>> stack;
    stack.emplace_back(0, n - 1);
    while (!stack.empty()) {
        auto [i, j] = stack.back();
        stack.pop_back();
        if (i > j) continue;
        const int k = partner[i][j];
        pairs.emplace_back(static_cast<std::size_t>(i), static_cast<std::size_t>(k));
        if (i + 1 <= k - 1) stack.emplace_back(i + 1, k - 1);
        if (k + 1 <= j) stack.emplace_back(k + 1, j);
    }
    result.pairs = std::move(pairs);
    return result;
}

}  // namespace mwpm
