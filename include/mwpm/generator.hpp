#pragma once

#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include "mwpm/point.hpp"

namespace mwpm {

struct Triangle {
    std::array<Point, 3> vertices;
    // Number of points lying on edge i, where edge i goes from vertices[i] to
    // vertices[(i + 1) % 3]. The three counts always sum to an even number.
    std::array<int, 3> per_edge_counts{};
};

// Generates a random triangle whose vertex coordinates are bounded so that the
// resulting Euclidean distances stay numerically well behaved. The number of
// points placed on each edge is randomized independently, then adjusted so the
// total is even and does not exceed `max_total_points`.
Triangle generate_random_triangle(std::mt19937_64& rng, int max_total_points);

// Like `generate_random_triangle` but forces the total number of points to be
// exactly `exact_total` (which must be a positive even integer). The total is
// partitioned across the three edges by a random composition. This is the
// generator used for benchmark runs that want a clean N value rather than a
// distribution.
Triangle generate_random_triangle_exact(std::mt19937_64& rng, int exact_total);

// Generates `total` points distributed on the edges of `tri` (using
// `tri.per_edge_counts`) and returns them in cyclic counterclockwise order
// around the triangle. The placement along each edge is randomized.
std::vector<Point> generate_points_on_triangle(std::mt19937_64& rng, const Triangle& tri);

}  // namespace mwpm
