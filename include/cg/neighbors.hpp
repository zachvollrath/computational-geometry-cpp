#pragma once

#include "cg/points.hpp"

#include <string>
#include <vector>

namespace cg {

struct NeighborResult {
  std::size_t queries = 0;
  std::size_t k = 0;
  std::vector<std::size_t> indices;
  std::vector<double> distances;

  std::size_t index(std::size_t row, std::size_t slot) const { return indices[row * k + slot]; }
  double distance(std::size_t row, std::size_t slot) const { return distances[row * k + slot]; }
};

// k nearest neighbors. queries == nullptr searches the point set against itself
// and drops the zero-distance self match. method: auto, brute, kdtree.
NeighborResult nearest_neighbors(const Points& points, std::size_t k = 1,
                                 const Points* queries = nullptr,
                                 const std::string& method = "auto");

}  // namespace cg
