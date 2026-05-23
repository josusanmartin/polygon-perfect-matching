#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "mwpm/point.hpp"

namespace mwpm {

// A reference to one input point: edge index and position within that edge's
// sorted barycentric-coordinate list.
struct PointRef {
    int edge{0};
    int index{0};
};

// Polygon geometry read from JSON: CCW vertices and per-edge point parameters.
struct PolygonInput {
    std::vector<Point> vertices;
    // edge_points[e] holds sorted barycentric coordinates s in [0, 1] along edge
    // e (from vertices[e] toward vertices[(e + 1) % n]).
    std::vector<std::vector<double>> edge_points;
};

// Flat CCW point list plus parallel (edge, index) references for I/O round-trip.
struct BuiltGeometry {
    std::vector<Point> points;
    std::vector<PointRef> refs;
};

class GeometryError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

BuiltGeometry build_geometry(const PolygonInput& input);

}  // namespace mwpm
