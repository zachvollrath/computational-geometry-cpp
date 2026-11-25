#pragma once

#include "cg/points.hpp"

#include <array>
#include <string>
#include <vector>

namespace cg {

struct Triangulation {
  std::vector<std::array<int, 3>> triangles;
};

// Delaunay triangulation of a planar point set. Triangles are wound counter-clockwise.
// method: auto, bowyer_watson.
Triangulation delaunay(const Points& points, const std::string& method = "auto");

}  // namespace cg
