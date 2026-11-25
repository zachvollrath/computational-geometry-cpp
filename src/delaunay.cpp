#include "cg/delaunay.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cg {
namespace {

struct I256 {
  std::uint64_t d[4] = {};
  bool neg = false;

  int sign() const {
    if ((d[0] | d[1] | d[2] | d[3]) == 0) {
      return 0;
    }
    return neg ? -1 : 1;
  }
};

int cmp_mag(const I256& a, const I256& b) {
  for (int i = 3; i >= 0; --i) {
    if (a.d[i] != b.d[i]) {
      return a.d[i] < b.d[i] ? -1 : 1;
    }
  }
  return 0;
}

void add_u128_at(I256& value, unsigned __int128 part, int limb) {
  unsigned __int128 carry = part;
  for (int i = limb; i < 4 && carry != 0; ++i) {
    carry += value.d[i];
    value.d[i] = static_cast<std::uint64_t>(carry);
    carry >>= 64;
  }
  if (carry != 0) {
    throw std::runtime_error("exact predicate overflow");
  }
}

I256 mul_i128(__int128 a, __int128 b) {
  const bool neg = (a < 0) != (b < 0);
  const unsigned __int128 ua = a < 0 ? static_cast<unsigned __int128>(-a) : static_cast<unsigned __int128>(a);
  const unsigned __int128 ub = b < 0 ? static_cast<unsigned __int128>(-b) : static_cast<unsigned __int128>(b);
  const std::uint64_t a0 = static_cast<std::uint64_t>(ua);
  const std::uint64_t a1 = static_cast<std::uint64_t>(ua >> 64);
  const std::uint64_t b0 = static_cast<std::uint64_t>(ub);
  const std::uint64_t b1 = static_cast<std::uint64_t>(ub >> 64);
  I256 result;
  add_u128_at(result, static_cast<unsigned __int128>(a0) * b0, 0);
  add_u128_at(result, static_cast<unsigned __int128>(a0) * b1, 1);
  add_u128_at(result, static_cast<unsigned __int128>(a1) * b0, 1);
  add_u128_at(result, static_cast<unsigned __int128>(a1) * b1, 2);
  result.neg = neg && result.sign() != 0;
  return result;
}

I256 add_mag(const I256& a, const I256& b) {
  I256 result;
  unsigned __int128 carry = 0;
  for (int i = 0; i < 4; ++i) {
    carry += static_cast<unsigned __int128>(a.d[i]) + b.d[i];
    result.d[i] = static_cast<std::uint64_t>(carry);
    carry >>= 64;
  }
  if (carry != 0) {
    throw std::runtime_error("exact predicate overflow");
  }
  return result;
}

I256 sub_mag(const I256& a, const I256& b) {
  I256 result;
  unsigned __int128 borrow = 0;
  for (int i = 0; i < 4; ++i) {
    const unsigned __int128 current = a.d[i];
    const unsigned __int128 sub = static_cast<unsigned __int128>(b.d[i]) + borrow;
    if (current < sub) {
      result.d[i] = static_cast<std::uint64_t>((current + (static_cast<unsigned __int128>(1) << 64)) - sub);
      borrow = 1;
    } else {
      result.d[i] = static_cast<std::uint64_t>(current - sub);
      borrow = 0;
    }
  }
  return result;
}

I256 add_signed(const I256& a, const I256& b) {
  if (a.neg == b.neg) {
    I256 result = add_mag(a, b);
    result.neg = a.neg && result.sign() != 0;
    return result;
  }
  const int order = cmp_mag(a, b);
  if (order == 0) {
    return {};
  }
  if (order > 0) {
    I256 result = sub_mag(a, b);
    result.neg = a.neg && result.sign() != 0;
    return result;
  }
  I256 result = sub_mag(b, a);
  result.neg = b.neg && result.sign() != 0;
  return result;
}

struct GridPoint {
  long long x = 0;
  long long y = 0;
};

bool in_circumcircle(const std::vector<GridPoint>& points, const std::array<int, 3>& triangle,
                     const GridPoint& query) {
  const GridPoint& a = points[static_cast<std::size_t>(triangle[0])];
  const GridPoint& b = points[static_cast<std::size_t>(triangle[1])];
  const GridPoint& c = points[static_cast<std::size_t>(triangle[2])];
  const __int128 orient = (static_cast<__int128>(b.x) - a.x) * (static_cast<__int128>(c.y) - a.y) -
                          (static_cast<__int128>(b.y) - a.y) * (static_cast<__int128>(c.x) - a.x);
  if (orient == 0) {
    return true;
  }
  const __int128 ax = static_cast<__int128>(a.x) - query.x;
  const __int128 ay = static_cast<__int128>(a.y) - query.y;
  const __int128 bx = static_cast<__int128>(b.x) - query.x;
  const __int128 by = static_cast<__int128>(b.y) - query.y;
  const __int128 cx = static_cast<__int128>(c.x) - query.x;
  const __int128 cy = static_cast<__int128>(c.y) - query.y;
  I256 det = mul_i128(ax * ax + ay * ay, bx * cy - cx * by);
  det = add_signed(det, mul_i128(bx * bx + by * by, cx * ay - ax * cy));
  det = add_signed(det, mul_i128(cx * cx + cy * cy, ax * by - bx * ay));
  if (orient < 0) {
    if (det.sign() != 0) {
      det.neg = !det.neg;
    }
  }
  return det.sign() > 0;
}

void require_distinct(const Points& pts) {
  std::set<std::pair<long long, long long>> seen;
  for (std::size_t i = 0; i < pts.n; ++i) {
    const double sx = std::round(pts(i, 0) * 1e10);
    const double sy = std::round(pts(i, 1) * 1e10);
    if (!std::isfinite(sx) || !std::isfinite(sy)) {
      throw std::invalid_argument("Delaunay input points must be distinct");
    }
    if (!seen.emplace(static_cast<long long>(std::llround(sx)), static_cast<long long>(std::llround(sy)))
             .second) {
      throw std::invalid_argument("Delaunay input points must be distinct");
    }
  }
}

std::vector<std::array<int, 3>> bowyer_watson(const Points& pts) {
  const int count = static_cast<int>(pts.n);
  double min_x = pts(0, 0);
  double min_y = pts(0, 1);
  double max_x = min_x;
  double max_y = min_y;
  for (std::size_t i = 1; i < pts.n; ++i) {
    min_x = std::min(min_x, pts(i, 0));
    min_y = std::min(min_y, pts(i, 1));
    max_x = std::max(max_x, pts(i, 0));
    max_y = std::max(max_y, pts(i, 1));
  }
  const double span = std::max(max_x - min_x, max_y - min_y);
  const double scale = span == 0.0 ? 1.0 : span;
  std::vector<GridPoint> cloud(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    const double lx = (pts(static_cast<std::size_t>(i), 0) - min_x) / scale;
    const double ly = (pts(static_cast<std::size_t>(i), 1) - min_y) / scale;
    cloud[static_cast<std::size_t>(i)] = {std::llround(lx * 1000000000.0), std::llround(ly * 1000000000.0)};
  }
  const long long bound = 1000000000000000LL;
  cloud.push_back({-bound, -bound});
  cloud.push_back({bound, -bound});
  cloud.push_back({0, bound});

  std::vector<std::array<int, 3>> triangles{{count, count + 1, count + 2}};
  for (int index = 0; index < count; ++index) {
    std::vector<std::array<int, 3>> bad;
    for (const auto& triangle : triangles) {
      if (in_circumcircle(cloud, triangle, cloud[static_cast<std::size_t>(index)])) {
        bad.push_back(triangle);
      }
    }
    std::map<std::pair<int, int>, int> boundary;
    for (const auto& triangle : bad) {
      const std::pair<int, int> edges[3] = {
          std::minmax(triangle[0], triangle[1]),
          std::minmax(triangle[1], triangle[2]),
          std::minmax(triangle[2], triangle[0]),
      };
      for (const auto& edge : edges) {
        ++boundary[edge];
      }
    }
    std::vector<std::array<int, 3>> kept;
    for (const auto& triangle : triangles) {
      if (std::find(bad.begin(), bad.end(), triangle) == bad.end()) {
        kept.push_back(triangle);
      }
    }
    for (const auto& [edge, times] : boundary) {
      if (times == 1) {
        kept.push_back({edge.first, edge.second, index});
      }
    }
    triangles.swap(kept);
  }

  std::vector<std::array<int, 3>> result;
  for (auto triangle : triangles) {
    if (triangle[0] >= count || triangle[1] >= count || triangle[2] >= count) {
      continue;
    }
    const double ax = pts(static_cast<std::size_t>(triangle[0]), 0);
    const double ay = pts(static_cast<std::size_t>(triangle[0]), 1);
    const double bx = pts(static_cast<std::size_t>(triangle[1]), 0);
    const double by = pts(static_cast<std::size_t>(triangle[1]), 1);
    const double cx = pts(static_cast<std::size_t>(triangle[2]), 0);
    const double cy = pts(static_cast<std::size_t>(triangle[2]), 1);
    const double area = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
    if (area < 0.0) {
      std::swap(triangle[1], triangle[2]);
    }
    result.push_back(triangle);
  }
  if (result.empty()) {
    throw std::invalid_argument("points do not have a Delaunay triangulation");
  }
  return result;
}

}  // namespace

Triangulation delaunay(const Points& points, const std::string& method) {
  require_method(method, {"auto", "bowyer_watson"});
  check_points(points, 3, 2);
  require_distinct(points);
  Triangulation mesh;
  mesh.triangles = bowyer_watson(points);
  return mesh;
}

}  // namespace cg
