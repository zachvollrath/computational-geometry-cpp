#pragma once

#include "cg/points.hpp"

#include <array>
#include <string>
#include <vector>

namespace cg {

struct ConvexHull {
  std::vector<int> indices;
  std::vector<double> vertices;
  std::vector<std::array<int, 3>> faces;
};

// 2D hulls are counter-clockwise, starting at the lexicographically smallest vertex.
// 3D results are sorted vertex indices plus the hull triangles.
// method: auto, monotone, incremental.
ConvexHull convex_hull(const Points& points, const std::string& method = "auto");

}  // namespace cg
