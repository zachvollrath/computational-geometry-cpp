#pragma once

#include "cg/points.hpp"

#include <string>

namespace cg {

struct PolygonArea {
  double area = 0.0;
  double signed_area = 0.0;
};

// Signed area is positive for a counter-clockwise ring. A repeated closing vertex is ignored.
PolygonArea shoelace(const Points& polygon, const std::string& method = "auto");

}  // namespace cg
