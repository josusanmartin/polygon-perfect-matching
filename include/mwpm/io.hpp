#pragma once

#include <string>

#include "mwpm/geometry.hpp"
#include "mwpm/point.hpp"

namespace mwpm {

PolygonInput read_input_json(const std::string& path);

void write_output_json(const std::string& path,
                       const MatchingResult& result,
                       const BuiltGeometry& geometry);

void write_solution_svg(const std::string& path,
                        const PolygonInput& input,
                        const BuiltGeometry& geometry,
                        const MatchingResult& result);

}  // namespace mwpm
