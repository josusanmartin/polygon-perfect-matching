#include "mwpm/algorithms.hpp"

#include <limits>
#include <stdexcept>

namespace mwpm {

namespace {

void brute_recurse(const std::vector<Point>& points,
                   std::vector<bool>& used,
                   std::vector<Pair>& current,
                   double current_cost,
                   double& best_cost,
                   std::vector<Pair>& best_pairs) {
    if (current_cost >= best_cost) {
        // Cheap pruning: any partial matching whose accumulated cost has
        // already exceeded the current best cannot improve it.
        return;
    }

    // Find the lowest-indexed unmatched point; pair it with every other
    // unmatched candidate. This canonicalizes the recursion (no double-counting)
    // and yields all (2k - 1)!! pairings on 2k points.
    std::size_t i = 0;
    while (i < points.size() && used[i]) ++i;
    if (i == points.size()) {
        best_cost = current_cost;
        best_pairs = current;
        return;
    }

    used[i] = true;
    for (std::size_t j = i + 1; j < points.size(); ++j) {
        if (used[j]) continue;
        used[j] = true;
        current.emplace_back(i, j);
        brute_recurse(points, used, current, current_cost + distance(points[i], points[j]),
                      best_cost, best_pairs);
        current.pop_back();
        used[j] = false;
    }
    used[i] = false;
}

}  // namespace

MatchingResult brute_force_matching(const std::vector<Point>& points) {
    if (points.size() % 2 != 0) {
        throw std::invalid_argument("brute_force_matching requires an even number of points");
    }
    MatchingResult result;
    if (points.empty()) {
        return result;
    }

    std::vector<bool> used(points.size(), false);
    std::vector<Pair> current;
    current.reserve(points.size() / 2);
    double best_cost = std::numeric_limits<double>::infinity();
    std::vector<Pair> best_pairs;

    brute_recurse(points, used, current, 0.0, best_cost, best_pairs);

    result.cost = best_cost;
    result.pairs = std::move(best_pairs);
    return result;
}

}  // namespace mwpm
