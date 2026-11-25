#include "cg/closest_pair.hpp"

#include "cg/neighbors.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace cg {
namespace {

ClosestPair pack(const Points& points, std::size_t i, std::size_t j) {
  if (i == j) {
    throw std::runtime_error("closest pair search returned a single index");
  }
  if (i > j) {
    std::swap(i, j);
  }
  ClosestPair pair;
  pair.index_a = i;
  pair.index_b = j;
  pair.point_a.assign(points.data.begin() + static_cast<std::ptrdiff_t>(i * points.d),
                      points.data.begin() + static_cast<std::ptrdiff_t>((i + 1) * points.d));
  pair.point_b.assign(points.data.begin() + static_cast<std::ptrdiff_t>(j * points.d),
                      points.data.begin() + static_cast<std::ptrdiff_t>((j + 1) * points.d));
  pair.distance = distance(points, i, j);
  return pair;
}

ClosestPair brute(const Points& points) {
  double best = std::numeric_limits<double>::infinity();
  std::size_t bi = 0;
  std::size_t bj = 1;
  for (std::size_t i = 0; i < points.n; ++i) {
    for (std::size_t j = i + 1; j < points.n; ++j) {
      const double dist = distance(points, i, j);
      if (dist < best) {
        best = dist;
        bi = i;
        bj = j;
      }
    }
  }
  return pack(points, bi, bj);
}

ClosestPair divide_conquer(const Points& points) {
  check_points(points, 2, 2);
  std::vector<int> order(points.n);
  for (std::size_t i = 0; i < points.n; ++i) {
    order[i] = static_cast<int>(i);
  }
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    if (points(static_cast<std::size_t>(a), 0) != points(static_cast<std::size_t>(b), 0)) {
      return points(static_cast<std::size_t>(a), 0) < points(static_cast<std::size_t>(b), 0);
    }
    return a < b;
  });

  double best = std::numeric_limits<double>::infinity();
  int bi = 0;
  int bj = 1;
  auto consider = [&](int i, int j) {
    if (i == j) {
      return;
    }
    const double dx = points(static_cast<std::size_t>(i), 0) - points(static_cast<std::size_t>(j), 0);
    const double dy = points(static_cast<std::size_t>(i), 1) - points(static_cast<std::size_t>(j), 1);
    const double dist = std::hypot(dx, dy);
    if (dist < best) {
      best = dist;
      if (i < j) {
        bi = i;
        bj = j;
      } else {
        bi = j;
        bj = i;
      }
    }
  };

  auto solve = [&](auto&& self, int begin, int end) -> void {
    const int count = end - begin;
    if (count <= 3) {
      for (int a = begin; a < end; ++a) {
        for (int b = a + 1; b < end; ++b) {
          consider(order[static_cast<std::size_t>(a)], order[static_cast<std::size_t>(b)]);
        }
      }
      return;
    }
    const int mid = begin + count / 2;
    const double split_x = points(static_cast<std::size_t>(order[static_cast<std::size_t>(mid)]), 0);
    self(self, begin, mid);
    self(self, mid, end);
    std::vector<int> strip;
    strip.reserve(static_cast<std::size_t>(count));
    for (int i = begin; i < end; ++i) {
      const int index = order[static_cast<std::size_t>(i)];
      if (std::abs(points(static_cast<std::size_t>(index), 0) - split_x) <= best) {
        strip.push_back(index);
      }
    }
    if (strip.size() < 2) {
      return;
    }
    std::stable_sort(strip.begin(), strip.end(), [&](int a, int b) {
      if (points(static_cast<std::size_t>(a), 1) != points(static_cast<std::size_t>(b), 1)) {
        return points(static_cast<std::size_t>(a), 1) < points(static_cast<std::size_t>(b), 1);
      }
      return a < b;
    });
    for (std::size_t a = 0; a < strip.size(); ++a) {
      const std::size_t last = std::min(a + 8, strip.size());
      for (std::size_t b = a + 1; b < last; ++b) {
        consider(strip[a], strip[b]);
      }
    }
  };
  solve(solve, 0, static_cast<int>(points.n));
  return pack(points, static_cast<std::size_t>(bi), static_cast<std::size_t>(bj));
}

ClosestPair by_tree(const Points& points) {
  const NeighborResult neighbors = nearest_neighbors(points, 1, nullptr, "kdtree");
  double best = std::numeric_limits<double>::infinity();
  std::size_t bi = 0;
  std::size_t bj = 1;
  for (std::size_t i = 0; i < points.n; ++i) {
    const std::size_t j = neighbors.index(i, 0);
    const double dist = neighbors.distance(i, 0);
    if (dist < best) {
      best = dist;
      bi = std::min(i, j);
      bj = std::max(i, j);
    }
  }
  return pack(points, bi, bj);
}

}  // namespace

ClosestPair closest_pair(const Points& points, const std::string& method) {
  require_method(method, {"auto", "brute", "divide_conquer", "kdtree"});
  check_points(points, 2);
  std::string chosen = method;
  if (chosen == "auto") {
    chosen = points.d == 2 ? "divide_conquer" : "kdtree";
  }
  if (chosen == "brute") {
    return brute(points);
  }
  if (chosen == "divide_conquer") {
    return divide_conquer(points);
  }
  return by_tree(points);
}

}  // namespace cg
