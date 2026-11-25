#pragma once

#include "cg/points.hpp"

#include <string>
#include <vector>

namespace cg {

struct OrientedBox {
  std::size_t dim = 0;
  std::vector<double> center;
  std::vector<double> axes;
  std::vector<double> extents;
  std::vector<double> corners;
  double volume = 0.0;

  double axis(std::size_t col, std::size_t row) const { return axes[col * dim + row]; }
  bool contains(const Points& points, double tol = 1e-6) const;
};

// Minimum-area rectangle in 2D (one side lies on a hull edge) and the minimum box
// over frames built from 3D hull triangles. method: auto, calipers.
OrientedBox minimum_bounding_box(const Points& points, const std::string& method = "auto");

}  // namespace cg
