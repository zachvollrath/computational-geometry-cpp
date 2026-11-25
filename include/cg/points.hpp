#pragma once

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cg {

struct Points {
  std::size_t n = 0;
  std::size_t d = 0;
  std::vector<double> data;

  double operator()(std::size_t i, std::size_t axis) const { return data[i * d + axis]; }
  double& operator()(std::size_t i, std::size_t axis) { return data[i * d + axis]; }
};

inline Points make_points(std::size_t d, std::initializer_list<std::initializer_list<double>> rows) {
  Points p;
  p.d = d;
  p.n = rows.size();
  p.data.reserve(p.n * d);
  for (const auto& row : rows) {
    if (row.size() != d) {
      throw std::invalid_argument("row width does not match dimension");
    }
    p.data.insert(p.data.end(), row.begin(), row.end());
  }
  return p;
}

inline void check_points(const Points& p, std::size_t min_count, int dim = -1, int min_dim = 1,
                         int max_dim = -1) {
  if (p.d == 0 || p.data.size() != p.n * p.d) {
    throw std::invalid_argument("points must have shape (n, d)");
  }
  if (p.n < min_count) {
    throw std::invalid_argument("expected at least " + std::to_string(min_count) + " points, got " +
                                std::to_string(p.n));
  }
  if (dim >= 0 && p.d != static_cast<std::size_t>(dim)) {
    throw std::invalid_argument("expected points with dimension " + std::to_string(dim) +
                                ", got " + std::to_string(p.d));
  }
  if (static_cast<int>(p.d) < min_dim || (max_dim >= 0 && static_cast<int>(p.d) > max_dim)) {
    throw std::invalid_argument("unsupported point dimension " + std::to_string(p.d));
  }
  for (double v : p.data) {
    if (!std::isfinite(v)) {
      throw std::invalid_argument("points must be finite");
    }
  }
}

inline void require_method(const std::string& method, std::initializer_list<const char*> known) {
  for (const char* name : known) {
    if (method == name) {
      return;
    }
  }
  std::string list;
  for (const char* name : known) {
    if (!list.empty()) {
      list += ", ";
    }
    list += name;
  }
  throw std::invalid_argument("unknown method '" + method + "'; expected one of " + list);
}

inline double distance(const Points& p, std::size_t i, std::size_t j) {
  double sum = 0.0;
  for (std::size_t axis = 0; axis < p.d; ++axis) {
    const double delta = p(i, axis) - p(j, axis);
    sum += delta * delta;
  }
  return std::sqrt(sum);
}

inline double dot(const double* a, const double* b, std::size_t d) {
  double sum = 0.0;
  for (std::size_t i = 0; i < d; ++i) {
    sum += a[i] * b[i];
  }
  return sum;
}

}  // namespace cg
