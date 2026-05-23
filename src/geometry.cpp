#include "mwpm/geometry.hpp"

#include <cmath>
#include <sstream>

namespace mwpm {

namespace {

double signed_area(const std::vector<Point>& vertices) {
    const std::size_t n = vertices.size();
    double area = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = vertices[i];
        const Point& b = vertices[(i + 1) % n];
        area += a.x * b.y - b.x * a.y;
    }
    return 0.5 * area;
}

Point interpolate(const Point& a, const Point& b, double s) {
    return Point{a.x + s * (b.x - a.x), a.y + s * (b.y - a.y)};
}

double cross(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x);
}

// Returns +1 for CCW turn, -1 for CW, 0 for collinear.
int turn_sign(const Point& a, const Point& b, const Point& c, double eps = 1e-9) {
    const double cr = cross(a, b, c);
    if (cr > eps) return 1;
    if (cr < -eps) return -1;
    return 0;
}

void validate_convex_ccw_polygon(const std::vector<Point>& vertices) {
    const std::size_t n = vertices.size();
    int sign = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = vertices[i];
        const Point& b = vertices[(i + 1) % n];
        const Point& c = vertices[(i + 2) % n];
        const int s = turn_sign(a, b, c);
        if (s == 0) {
            throw GeometryError(
                "polygon vertices are collinear or degenerate at a corner");
        }
        if (sign == 0) {
            sign = s;
        } else if (sign != s) {
            throw GeometryError("polygon vertices must form a convex polygon");
        }
    }
    if (sign != 1) {
        throw GeometryError("polygon vertices must be listed in counterclockwise order");
    }
}

}  // namespace

BuiltGeometry build_geometry(const PolygonInput& input) {
    const std::size_t n_edges = input.vertices.size();
    if (n_edges < 3) {
        throw GeometryError("vertices must contain at least 3 points");
    }
    if (input.edge_points.size() != n_edges) {
        throw GeometryError("edge_points length must match vertices length");
    }
    if (std::abs(signed_area(input.vertices)) < 1e-9) {
        throw GeometryError("polygon vertices are degenerate (zero area)");
    }
    validate_convex_ccw_polygon(input.vertices);

    BuiltGeometry geom;
    std::size_t total = 0;
    for (const auto& edge : input.edge_points) {
        total += edge.size();
    }
    if (total < 2) {
        throw GeometryError("need at least 2 boundary points in total");
    }
    if (total % 2 != 0) {
        throw GeometryError("total number of boundary points must be even");
    }

    geom.points.reserve(total);
    geom.refs.reserve(total);

    for (std::size_t e = 0; e < n_edges; ++e) {
        const auto& params = input.edge_points[e];
        for (std::size_t i = 1; i < params.size(); ++i) {
            if (params[i] < params[i - 1]) {
                std::ostringstream oss;
                oss << "edge_points[" << e << "] must be sorted ascending";
                throw GeometryError(oss.str());
            }
        }
        const Point& a = input.vertices[e];
        const Point& b = input.vertices[(e + 1) % n_edges];
        for (std::size_t j = 0; j < params.size(); ++j) {
            const double s = params[j];
            if (s < 0.0 || s > 1.0) {
                std::ostringstream oss;
                oss << "edge_points[" << e << "][" << j << "] = " << s
                    << " is outside [0, 1]";
                throw GeometryError(oss.str());
            }
            geom.points.push_back(interpolate(a, b, s));
            geom.refs.push_back(PointRef{static_cast<int>(e), static_cast<int>(j)});
        }
    }

    return geom;
}

}  // namespace mwpm
