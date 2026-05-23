#pragma once

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace mwpm {

struct Point {
    double x{0.0};
    double y{0.0};

    Point() = default;
    constexpr Point(double x_, double y_) : x(x_), y(y_) {}
};

inline double distance(const Point& a, const Point& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

// A matched pair represented by indices into the original cyclic point vector.
using Pair = std::pair<std::size_t, std::size_t>;

struct MatchingResult {
    double cost{0.0};
    std::vector<Pair> pairs;
};

}  // namespace mwpm
