#include "cg/laplacian.hpp"

#include "cg/delaunay.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cg {
namespace {

std::vector<std::array<int, 2>> edges_of(const std::vector<std::array<int, 3>>& faces) {
  std::vector<std::array<int, 2>> edges;
  for (const auto& face : faces) {
    const int loop[3][2] = {{face[0], face[1]}, {face[1], face[2]}, {face[2], face[0]}};
    for (const auto& edge : loop) {
      if (edge[0] == edge[1]) {
        continue;
      }
      const int a = std::min(edge[0], edge[1]);
      const int b = std::max(edge[0], edge[1]);
      const std::array<int, 2> key{a, b};
      if (std::find(edges.begin(), edges.end(), key) == edges.end()) {
        edges.push_back(key);
      }
    }
  }
  return edges;
}

double weight_of(const Points& pts, int i, int j, const std::string& weight) {
  if (weight == "uniform") {
    return 1.0;
  }
  double sum = 0.0;
  for (std::size_t axis = 0; axis < pts.d; ++axis) {
    const double delta = pts(static_cast<std::size_t>(i), axis) - pts(static_cast<std::size_t>(j), axis);
    sum += delta * delta;
  }
  return 1.0 / (std::sqrt(sum) + 1e-12);
}

Points step(const Points& pts, const std::vector<std::array<int, 2>>& edges, double lambda,
            const std::string& weight) {
  std::vector<double> acc(pts.data.size(), 0.0);
  std::vector<double> wsum(pts.n, 0.0);
  for (const auto& edge : edges) {
    const double w = weight_of(pts, edge[0], edge[1], weight);
    for (std::size_t axis = 0; axis < pts.d; ++axis) {
      acc[static_cast<std::size_t>(edge[0]) * pts.d + axis] +=
          w * pts(static_cast<std::size_t>(edge[1]), axis);
      acc[static_cast<std::size_t>(edge[1]) * pts.d + axis] +=
          w * pts(static_cast<std::size_t>(edge[0]), axis);
    }
    wsum[static_cast<std::size_t>(edge[0])] += w;
    wsum[static_cast<std::size_t>(edge[1])] += w;
  }
  Points updated = pts;
  for (std::size_t i = 0; i < pts.n; ++i) {
    if (wsum[i] <= 0.0) {
      continue;
    }
    for (std::size_t axis = 0; axis < pts.d; ++axis) {
      const double average = acc[i * pts.d + axis] / wsum[i];
      updated(i, axis) = pts(i, axis) + lambda * (average - pts(i, axis));
    }
  }
  return updated;
}

}  // namespace

Points laplacian_smooth(const Points& points, const std::vector<std::array<int, 3>>* faces, int iterations,
                        double lambda, const std::string& weight, const std::string& method) {
  require_method(method, {"auto", "direct"});
  require_method(weight, {"uniform", "inverse_distance"});
  if (iterations < 0) {
    throw std::invalid_argument("iterations must be non-negative");
  }
  if (lambda < 0.0 || lambda > 1.0) {
    throw std::invalid_argument("lamb must lie in [0, 1]");
  }
  check_points(points, 1, -1, 2, 3);
  std::vector<std::array<int, 3>> owned;
  const std::vector<std::array<int, 3>>* tris = faces;
  if (faces == nullptr) {
    if (points.d != 2) {
      throw std::invalid_argument("faces are required for a 3D mesh");
    }
    owned = delaunay(points).triangles;
    tris = &owned;
  } else {
    for (const auto& face : *faces) {
      for (int index : face) {
        if (index < 0 || static_cast<std::size_t>(index) >= points.n) {
          throw std::invalid_argument("face indices are out of range");
        }
      }
    }
  }
  const auto edges = edges_of(*tris);
  Points current = points;
  for (int i = 0; i < iterations; ++i) {
    current = step(current, edges, lambda, weight);
  }
  return current;
}

}  // namespace cg
