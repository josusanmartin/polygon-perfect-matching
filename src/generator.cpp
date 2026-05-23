#include "mwpm/generator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mwpm {

namespace {

// Returns the signed area of triangle (a, b, c). Positive if CCW.
double signed_area(const Point& a, const Point& b, const Point& c) {
    return 0.5 * ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

// Sample a non-degenerate, CCW-oriented triangle whose vertices are bounded
// in [-100, 100]^2. Shared by both generator entry points.
void sample_triangle_vertices(std::mt19937_64& rng, std::array<Point, 3>& v) {
    std::uniform_real_distribution<double> coord(-100.0, 100.0);
    while (true) {
        for (auto& p : v) p = Point{coord(rng), coord(rng)};
        const double area = signed_area(v[0], v[1], v[2]);
        if (std::abs(area) < 1.0) continue;          // reject (near-)degenerate
        if (area < 0.0) std::swap(v[1], v[2]);       // force CCW
        break;
    }
}

}  // namespace

Triangle generate_random_triangle(std::mt19937_64& rng, int max_total_points) {
    Triangle tri;
    sample_triangle_vertices(rng, tri.vertices);

    // Distribute points on the three edges. Each edge gets between 0 and
    // ceil(max_total_points / 3) points; we then bump one count by 1 if the
    // total is odd to satisfy the parity constraint. We force at least 1 here
    // so the resampling loop can terminate even when max_total_points is very
    // small (e.g. max_total_points == 2).
    const int max_per_edge = std::max(1, (max_total_points + 2) / 3);
    std::uniform_int_distribution<int> count_dist(0, max_per_edge);

    while (true) {
        for (int i = 0; i < 3; ++i) {
            tri.per_edge_counts[i] = count_dist(rng);
        }
        int total = tri.per_edge_counts[0] + tri.per_edge_counts[1] + tri.per_edge_counts[2];
        if (total % 2 != 0) {
            // Bump a uniformly chosen edge to fix parity.
            std::uniform_int_distribution<int> idx(0, 2);
            ++tri.per_edge_counts[idx(rng)];
            ++total;
        }
        if (total > max_total_points) {
            continue;  // resample to respect the cap
        }
        if (total >= 2) {
            break;  // need at least one pair to match
        }
    }

    return tri;
}

Triangle generate_random_triangle_exact(std::mt19937_64& rng, int exact_total) {
    if (exact_total < 2 || exact_total % 2 != 0) {
        throw std::invalid_argument(
            "generate_random_triangle_exact: exact_total must be even and >= 2");
    }
    Triangle tri;
    sample_triangle_vertices(rng, tri.vertices);

    // Random composition of exact_total into three non-negative integers via
    // two i.i.d. uniform draws. The resulting marginal distribution is
    // slightly biased (composition (a, b, c) is reached with probability
    // proportional to 1 / (1 + (exact_total - a))) but for benchmarking we
    // only need *some* random spread across the edges.
    const int e1 = std::uniform_int_distribution<int>(0, exact_total)(rng);
    const int e2 = std::uniform_int_distribution<int>(0, exact_total - e1)(rng);
    const int e3 = exact_total - e1 - e2;
    tri.per_edge_counts = {e1, e2, e3};
    return tri;
}

std::vector<Point> generate_points_on_triangle(std::mt19937_64& rng, const Triangle& tri) {
    std::vector<Point> result;
    const int total = tri.per_edge_counts[0] + tri.per_edge_counts[1] + tri.per_edge_counts[2];
    result.reserve(static_cast<std::size_t>(total));

    std::uniform_real_distribution<double> param(0.0, 1.0);

    // For each CCW edge i: a = vertices[i], b = vertices[(i+1) % 3]. We sample
    // n_i parameters in (0, 1), sort them, then emit a + t * (b - a) for each
    // parameter. Sorting ensures the points appear in CCW order along the edge
    // (and therefore in CCW order around the triangle perimeter).
    for (int i = 0; i < 3; ++i) {
        const int n_i = tri.per_edge_counts[i];
        if (n_i == 0) continue;

        std::vector<double> ts(static_cast<std::size_t>(n_i));
        for (auto& t : ts) {
            // Strictly interior placements (avoid 0 and 1 to prevent points
            // landing exactly on triangle vertices, which would make distance
            // bookkeeping ambiguous).
            do {
                t = param(rng);
            } while (t <= 1e-9 || t >= 1.0 - 1e-9);
        }
        std::sort(ts.begin(), ts.end());

        const Point& a = tri.vertices[i];
        const Point& b = tri.vertices[(i + 1) % 3];
        for (double t : ts) {
            result.push_back(Point{a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)});
        }
    }

    return result;
}

}  // namespace mwpm
