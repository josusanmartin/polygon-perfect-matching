#include "mwpm/io.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

namespace mwpm {

namespace {

using json = nlohmann::json;

class IoError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

json point_ref_json(const PointRef& ref) {
    return json::array({ref.edge, ref.index});
}

struct Viewport {
    double scale{1.0};
    double offset_x{0.0};
    double offset_y{0.0};

    double map_x(double x) const { return offset_x + scale * x; }
    double map_y(double y) const { return offset_y + scale * y; }
};

Viewport compute_viewport(const std::vector<Point>& points, double canvas = 100.0,
                          double margin = 6.0) {
    double min_x = points.front().x;
    double max_x = points.front().x;
    double min_y = points.front().y;
    double max_y = points.front().y;
    for (const Point& p : points) {
        min_x = std::min(min_x, p.x);
        max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y);
        max_y = std::max(max_y, p.y);
    }

    const double width = std::max(max_x - min_x, 1e-9);
    const double height = std::max(max_y - min_y, 1e-9);
    const double inner = canvas - 2.0 * margin;
    const double scale = inner / std::max(width, height);

    Viewport vp;
    vp.scale = scale;
    vp.offset_x = margin + 0.5 * (inner - scale * width) - scale * min_x;
    vp.offset_y = margin + 0.5 * (inner - scale * height) - scale * min_y;
    return vp;
}

std::string svg_header() {
    return R"(<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100" width="100" height="100">
  <rect x="0" y="0" width="100" height="100" fill="#fafafa"/>
)";
}

}  // namespace

PolygonInput read_input_json(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw IoError("failed to open input file: " + path);
    }

    json doc;
    try {
        in >> doc;
    } catch (const json::parse_error& e) {
        throw IoError(std::string("JSON parse error: ") + e.what());
    }

    if (!doc.contains("vertices") || !doc["vertices"].is_array()) {
        throw IoError("input must contain a \"vertices\" array");
    }
    if (!doc.contains("edge_points") || !doc["edge_points"].is_array()) {
        throw IoError("input must contain an \"edge_points\" array");
    }

    PolygonInput input;
    for (const auto& v : doc["vertices"]) {
        if (!v.is_array() || v.size() != 2 || !v[0].is_number() || !v[1].is_number()) {
            throw IoError("each vertex must be a pair of numbers [x, y]");
        }
        input.vertices.emplace_back(v[0].get<double>(), v[1].get<double>());
    }

    for (const auto& edge : doc["edge_points"]) {
        if (!edge.is_array()) {
            throw IoError("each edge_points entry must be an array of numbers");
        }
        std::vector<double> params;
        params.reserve(edge.size());
        for (const auto& s : edge) {
            if (!s.is_number()) {
                throw IoError("edge point coordinates must be numbers");
            }
            params.push_back(s.get<double>());
        }
        input.edge_points.push_back(std::move(params));
    }

    return input;
}

void write_output_json(const std::string& path,
                       const MatchingResult& result,
                       const BuiltGeometry& geometry) {
    json segments = json::array();
    for (const auto& [i, j] : result.pairs) {
        segments.push_back(json::array({point_ref_json(geometry.refs[i]),
                                        point_ref_json(geometry.refs[j])}));
    }

    json out{{"cost", result.cost}, {"segments", std::move(segments)}};

    std::ofstream os(path);
    if (!os) {
        throw IoError("failed to open output file: " + path);
    }
    os << std::setprecision(17) << out.dump(2) << '\n';
}

void write_solution_svg(const std::string& path,
                        const PolygonInput& input,
                        const BuiltGeometry& geometry,
                        const MatchingResult& result) {
    std::vector<Point> all;
    all.reserve(input.vertices.size() + geometry.points.size());
    all.insert(all.end(), input.vertices.begin(), input.vertices.end());
    all.insert(all.end(), geometry.points.begin(), geometry.points.end());
    const Viewport vp = compute_viewport(all);

    constexpr double kPolygonStroke = 0.45;
    constexpr double kSegmentStroke = 2.0 * kPolygonStroke;

    std::ostringstream svg;
    svg << svg_header();

    svg << "  <polygon fill=\"none\" stroke=\"#111111\" stroke-width=\""
        << kPolygonStroke << "\" stroke-linejoin=\"round\" points=\"";
    const std::size_t n = input.vertices.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Point& p = input.vertices[i];
        if (i > 0) svg << ' ';
        svg << vp.map_x(p.x) << ',' << vp.map_y(p.y);
    }
    svg << "\"/>\n";

    svg << "  <g stroke=\"#2563eb\" stroke-width=\"" << kSegmentStroke
        << "\" stroke-linecap=\"round\">\n";
    for (const auto& [i, j] : result.pairs) {
        const Point& a = geometry.points[i];
        const Point& b = geometry.points[j];
        svg << "    <line x1=\"" << vp.map_x(a.x) << "\" y1=\"" << vp.map_y(a.y)
            << "\" x2=\"" << vp.map_x(b.x) << "\" y2=\"" << vp.map_y(b.y) << "\"/>\n";
    }
    svg << "  </g>\n";

    svg << "  <g fill=\"#111111\">\n";
    for (const Point& p : geometry.points) {
        svg << "    <circle cx=\"" << vp.map_x(p.x) << "\" cy=\"" << vp.map_y(p.y)
            << "\" r=\"0.95\"/>\n";
    }
    svg << "  </g>\n";

    svg << "</svg>\n";

    std::ofstream out(path);
    if (!out) {
        throw IoError("failed to open SVG file: " + path);
    }
    out << svg.str();
}

}  // namespace mwpm
