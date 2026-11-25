#pragma once

#include "cg/points.hpp"

#include <string>
#include <vector>

namespace cg {

struct ClosestPair {
  std::size_t index_a = 0;
  std::size_t index_b = 0;
  std::vector<double> point_a;
  std::vector<double> point_b;
  double distance = 0.0;
};

// method: auto, brute, divide_conquer, kdtree
// auto uses divide-and-conquer in 2D and a KD-tree otherwise.
ClosestPair closest_pair(const Points& points, const std::string& method = "auto");

}  // namespace cg
