#include "cg/convex_hull.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

namespace cg {
namespace {

double cross2(const Points& pts, int origin, int a, int b) {
  const double ax = pts(static_cast<std::size_t>(a), 0) - pts(static_cast<std::size_t>(origin), 0);
  const double ay = pts(static_cast<std::size_t>(a), 1) - pts(static_cast<std::size_t>(origin), 1);
  const double bx = pts(static_cast<std::size_t>(b), 0) - pts(static_cast<std::size_t>(origin), 0);
  const double by = pts(static_cast<std::size_t>(b), 1) - pts(static_cast<std::size_t>(origin), 1);
  return ax * by - ay * bx;
}

std::vector<int> chain(const Points& pts, const std::vector<int>& order) {
  std::vector<int> built;
  for (int index : order) {
    while (built.size() >= 2 &&
           cross2(pts, built[built.size() - 2], built.back(), index) <= 0.0) {
      built.pop_back();
    }
    built.push_back(index);
  }
  return built;
}

std::vector<int> dedup_ring(const std::vector<int>& indices, const Points& pts) {
  std::vector<int> cleaned;
  for (int index : indices) {
    if (cleaned.empty() ||
        std::hypot(pts(static_cast<std::size_t>(index), 0) - pts(static_cast<std::size_t>(cleaned.back()), 0),
                   pts(static_cast<std::size_t>(index), 1) - pts(static_cast<std::size_t>(cleaned.back()), 1)) >
            1e-12) {
      cleaned.push_back(index);
    }
  }
  if (cleaned.size() > 1 &&
      std::hypot(pts(static_cast<std::size_t>(cleaned.front()), 0) - pts(static_cast<std::size_t>(cleaned.back()), 0),
                 pts(static_cast<std::size_t>(cleaned.front()), 1) - pts(static_cast<std::size_t>(cleaned.back()), 1)) <=
          1e-12) {
    cleaned.pop_back();
  }
  if (cleaned.empty()) {
    cleaned.push_back(indices.empty() ? 0 : indices.front());
  }
  return cleaned;
}

std::vector<int> canonicalize_2d(std::vector<int> indices, const Points& pts) {
  std::vector<int> unique;
  for (int index : indices) {
    if (std::find(unique.begin(), unique.end(), index) == unique.end()) {
      unique.push_back(index);
    }
  }
  indices.swap(unique);
  if (indices.size() >= 3) {
    double area = 0.0;
    for (std::size_t pos = 0; pos < indices.size(); ++pos) {
      const int a = indices[pos];
      const int b = indices[(pos + 1) % indices.size()];
      area += pts(static_cast<std::size_t>(a), 0) * pts(static_cast<std::size_t>(b), 1) -
              pts(static_cast<std::size_t>(a), 1) * pts(static_cast<std::size_t>(b), 0);
    }
    if (area < 0.0) {
      std::reverse(indices.begin(), indices.end());
    }
  }
  if (indices.empty()) {
    return indices;
  }
  std::size_t start = 0;
  for (std::size_t pos = 1; pos < indices.size(); ++pos) {
    const int cur = indices[pos];
    const int best = indices[start];
    const double cx = pts(static_cast<std::size_t>(cur), 0);
    const double cy = pts(static_cast<std::size_t>(cur), 1);
    const double bx = pts(static_cast<std::size_t>(best), 0);
    const double by = pts(static_cast<std::size_t>(best), 1);
    if (std::tie(cx, cy, cur) < std::tie(bx, by, best)) {
      start = pos;
    }
  }
  std::rotate(indices.begin(), indices.begin() + static_cast<std::ptrdiff_t>(start), indices.end());
  return indices;
}

std::vector<int> order_2d(std::vector<int> indices, const Points& pts) {
  std::vector<int> ordered;
  for (int index : indices) {
    if (std::find(ordered.begin(), ordered.end(), index) == ordered.end()) {
      ordered.push_back(index);
    }
  }
  if (ordered.size() >= 3) {
    double cx = 0.0;
    double cy = 0.0;
    for (int index : ordered) {
      cx += pts(static_cast<std::size_t>(index), 0);
      cy += pts(static_cast<std::size_t>(index), 1);
    }
    cx /= static_cast<double>(ordered.size());
    cy /= static_cast<double>(ordered.size());
    std::vector<double> angles(ordered.size());
    for (std::size_t i = 0; i < ordered.size(); ++i) {
      angles[i] = std::atan2(pts(static_cast<std::size_t>(ordered[i]), 1) - cy,
                             pts(static_cast<std::size_t>(ordered[i]), 0) - cx);
    }
    std::vector<std::size_t> perm(ordered.size());
    for (std::size_t i = 0; i < perm.size(); ++i) {
      perm[i] = i;
    }
    std::stable_sort(perm.begin(), perm.end(), [&](std::size_t a, std::size_t b) {
      if (angles[a] != angles[b]) {
        return angles[a] < angles[b];
      }
      return ordered[a] < ordered[b];
    });
    std::vector<int> sorted;
    sorted.reserve(ordered.size());
    for (std::size_t i : perm) {
      sorted.push_back(ordered[i]);
    }
    ordered.swap(sorted);
  }
  return canonicalize_2d(std::move(ordered), pts);
}

std::vector<int> monotone_indices(const Points& pts) {
  std::vector<int> order(pts.n);
  for (std::size_t i = 0; i < pts.n; ++i) {
    order[i] = static_cast<int>(i);
  }
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    const double ax = pts(static_cast<std::size_t>(a), 0);
    const double bx = pts(static_cast<std::size_t>(b), 0);
    if (ax != bx) {
      return ax < bx;
    }
    const double ay = pts(static_cast<std::size_t>(a), 1);
    const double by = pts(static_cast<std::size_t>(b), 1);
    if (ay != by) {
      return ay < by;
    }
    return a < b;
  });
  const std::vector<int> lower = chain(pts, order);
  std::vector<int> reversed(order.rbegin(), order.rend());
  const std::vector<int> upper = chain(pts, reversed);
  std::vector<int> ring;
  if (!lower.empty()) {
    ring.insert(ring.end(), lower.begin(), lower.end() - 1);
  }
  if (!upper.empty()) {
    ring.insert(ring.end(), upper.begin(), upper.end() - 1);
  }
  return order_2d(dedup_ring(ring, pts), pts);
}

ConvexHull pack(const Points& pts, std::vector<int> indices, std::vector<std::array<int, 3>> faces) {
  ConvexHull hull;
  hull.indices = std::move(indices);
  hull.faces = std::move(faces);
  hull.vertices.resize(hull.indices.size() * pts.d);
  for (std::size_t i = 0; i < hull.indices.size(); ++i) {
    for (std::size_t axis = 0; axis < pts.d; ++axis) {
      hull.vertices[i * pts.d + axis] = pts(static_cast<std::size_t>(hull.indices[i]), axis);
    }
  }
  return hull;
}

void sub3(const Points& pts, int index, const double origin[3], double out[3]) {
  for (int axis = 0; axis < 3; ++axis) {
    out[axis] = pts(static_cast<std::size_t>(index), static_cast<std::size_t>(axis)) - origin[axis];
  }
}

void cross3(const double a[3], const double b[3], double out[3]) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

double dot3(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm3(const double a[3]) { return std::sqrt(dot3(a, a)); }

double orient_point(const Points& pts, const std::array<int, 3>& face, const double q[3]) {
  double ab[3], ac[3], ad[3];
  const double origin[3] = {pts(static_cast<std::size_t>(face[0]), 0),
                            pts(static_cast<std::size_t>(face[0]), 1),
                            pts(static_cast<std::size_t>(face[0]), 2)};
  sub3(pts, face[1], origin, ab);
  sub3(pts, face[2], origin, ac);
  for (int axis = 0; axis < 3; ++axis) {
    ad[axis] = q[axis] - origin[axis];
  }
  double normal[3];
  cross3(ab, ac, normal);
  return dot3(normal, ad);
}

double coord_scale(const Points& pts) {
  double scale = 1.0;
  for (double value : pts.data) {
    scale = std::max(scale, std::abs(value));
  }
  return scale;
}

bool point_in_triangle(const Points& pts, const std::array<int, 3>& face, const double q[3],
                       double eps) {
  double ab[3], ac[3];
  const double origin[3] = {pts(static_cast<std::size_t>(face[0]), 0),
                            pts(static_cast<std::size_t>(face[0]), 1),
                            pts(static_cast<std::size_t>(face[0]), 2)};
  sub3(pts, face[1], origin, ab);
  sub3(pts, face[2], origin, ac);
  double normal[3];
  cross3(ab, ac, normal);
  const int drop = std::abs(normal[0]) >= std::abs(normal[1]) && std::abs(normal[0]) >= std::abs(normal[2])
                       ? 0
                       : (std::abs(normal[1]) >= std::abs(normal[2]) ? 1 : 2);
  const int u = (drop + 1) % 3;
  const int v = (drop + 2) % 3;
  const double ax = origin[u];
  const double ay = origin[v];
  const double bx = pts(static_cast<std::size_t>(face[1]), static_cast<std::size_t>(u));
  const double by = pts(static_cast<std::size_t>(face[1]), static_cast<std::size_t>(v));
  const double cx = pts(static_cast<std::size_t>(face[2]), static_cast<std::size_t>(u));
  const double cy = pts(static_cast<std::size_t>(face[2]), static_cast<std::size_t>(v));
  const double px = q[u];
  const double py = q[v];
  auto side = [&](double x0, double y0, double x1, double y1, double x2, double y2) {
    return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
  };
  const double c1 = side(ax, ay, bx, by, px, py);
  const double c2 = side(bx, by, cx, cy, px, py);
  const double c3 = side(cx, cy, ax, ay, px, py);
  const bool neg = c1 < -eps || c2 < -eps || c3 < -eps;
  const bool pos = c1 > eps || c2 > eps || c3 > eps;
  return !(neg && pos);
}

struct Mesh {
  std::vector<std::array<int, 3>> faces;
  double inside[3] = {};
  double eps = 1e-9;
  double area_eps = 1e-8;
};

void flip_against_inside(const Points& pts, std::array<int, 3>& face, const double inside[3],
                         double eps) {
  if (orient_point(pts, face, inside) > eps) {
    std::swap(face[1], face[2]);
  }
}

void add_strict(Mesh& mesh, const Points& pts, int index) {
  double q[3] = {pts(static_cast<std::size_t>(index), 0), pts(static_cast<std::size_t>(index), 1),
                 pts(static_cast<std::size_t>(index), 2)};
  std::vector<int> visible;
  for (int face = 0; face < static_cast<int>(mesh.faces.size()); ++face) {
    if (orient_point(pts, mesh.faces[static_cast<std::size_t>(face)], q) > mesh.eps) {
      visible.push_back(face);
    }
  }
  if (visible.empty()) {
    return;
  }
  std::map<std::pair<int, int>, int> count;
  std::map<std::pair<int, int>, std::pair<int, int>> directed;
  for (int face_index : visible) {
    const auto& face = mesh.faces[static_cast<std::size_t>(face_index)];
    const int edges[3][2] = {{face[0], face[1]}, {face[1], face[2]}, {face[2], face[0]}};
    for (const auto& edge : edges) {
      const std::pair<int, int> key = std::minmax(edge[0], edge[1]);
      if (count[key] == 0) {
        directed[key] = {edge[0], edge[1]};
      }
      ++count[key];
    }
  }
  std::vector<char> dead(mesh.faces.size(), 0);
  for (int face_index : visible) {
    dead[static_cast<std::size_t>(face_index)] = 1;
  }
  std::vector<std::array<int, 3>> next;
  for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
    if (!dead[i]) {
      next.push_back(mesh.faces[i]);
    }
  }
  for (const auto& [key, times] : count) {
    if (times != 1) {
      continue;
    }
    const auto [u, v] = directed[key];
    std::array<int, 3> face{u, v, index};
    flip_against_inside(pts, face, mesh.inside, mesh.eps);
    next.push_back(face);
  }
  mesh.faces.swap(next);
}

bool absorb_coplanar(Mesh& mesh, const Points& pts, std::vector<char>& skip) {
  for (int index = 0; index < static_cast<int>(pts.n); ++index) {
    if (skip[static_cast<std::size_t>(index)]) {
      continue;
    }
    double q[3] = {pts(static_cast<std::size_t>(index), 0), pts(static_cast<std::size_t>(index), 1),
                   pts(static_cast<std::size_t>(index), 2)};
    bool strict_out = false;
    bool on_triangle = false;
    std::vector<int> coplanar;
    for (int face = 0; face < static_cast<int>(mesh.faces.size()); ++face) {
      const double volume = orient_point(pts, mesh.faces[static_cast<std::size_t>(face)], q);
      if (volume > mesh.eps) {
        strict_out = true;
        break;
      }
      if (std::abs(volume) <= mesh.eps) {
        coplanar.push_back(face);
        if (point_in_triangle(pts, mesh.faces[static_cast<std::size_t>(face)], q, mesh.area_eps)) {
          on_triangle = true;
        }
      }
    }
    if (strict_out) {
      add_strict(mesh, pts, index);
      return true;
    }
    if (coplanar.empty() || on_triangle) {
      continue;
    }
    std::vector<int> used;
    for (int face_index : coplanar) {
      for (int vertex : mesh.faces[static_cast<std::size_t>(face_index)]) {
        if (std::find(used.begin(), used.end(), vertex) == used.end()) {
          used.push_back(vertex);
        }
      }
    }
    used.push_back(index);
    double ab[3], ac[3], normal[3];
    const auto& seed = mesh.faces[static_cast<std::size_t>(coplanar.front())];
    const double origin[3] = {pts(static_cast<std::size_t>(seed[0]), 0),
                              pts(static_cast<std::size_t>(seed[0]), 1),
                              pts(static_cast<std::size_t>(seed[0]), 2)};
    sub3(pts, seed[1], origin, ab);
    sub3(pts, seed[2], origin, ac);
    cross3(ab, ac, normal);
    const int drop =
        std::abs(normal[0]) >= std::abs(normal[1]) && std::abs(normal[0]) >= std::abs(normal[2])
            ? 0
            : (std::abs(normal[1]) >= std::abs(normal[2]) ? 1 : 2);
    const int u = (drop + 1) % 3;
    const int v = (drop + 2) % 3;
    Points flat;
    flat.n = used.size();
    flat.d = 2;
    flat.data.resize(flat.n * 2);
    for (std::size_t i = 0; i < used.size(); ++i) {
      if (used[i] == index) {
        flat(i, 0) = q[u];
        flat(i, 1) = q[v];
      } else {
        flat(i, 0) = pts(static_cast<std::size_t>(used[i]), static_cast<std::size_t>(u));
        flat(i, 1) = pts(static_cast<std::size_t>(used[i]), static_cast<std::size_t>(v));
      }
    }
    const std::vector<int> ring = monotone_indices(flat);
    bool includes = false;
    for (int local : ring) {
      if (used[static_cast<std::size_t>(local)] == index) {
        includes = true;
      }
    }
    if (!includes || ring.size() < 3) {
      skip[static_cast<std::size_t>(index)] = 1;
      continue;
    }
    std::vector<char> dead(mesh.faces.size(), 0);
    for (int face_index : coplanar) {
      dead[static_cast<std::size_t>(face_index)] = 1;
    }
    std::vector<std::array<int, 3>> next;
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
      if (!dead[i]) {
        next.push_back(mesh.faces[i]);
      }
    }
    const int fan = used[static_cast<std::size_t>(ring.front())];
    for (std::size_t k = 1; k + 1 < ring.size(); ++k) {
      std::array<int, 3> face{fan, used[static_cast<std::size_t>(ring[k])],
                              used[static_cast<std::size_t>(ring[k + 1])]};
      flip_against_inside(pts, face, mesh.inside, mesh.eps);
      next.push_back(face);
    }
    mesh.faces.swap(next);
    return true;
  }
  return false;
}

std::vector<int> extremes_along_farthest(const Points& pts) {
  if (pts.n == 1) {
    return {0};
  }
  int a = 0;
  int b = 1;
  double best = -1.0;
  for (int i = 0; i < static_cast<int>(pts.n); ++i) {
    for (int j = i + 1; j < static_cast<int>(pts.n); ++j) {
      double sum = 0.0;
      for (std::size_t axis = 0; axis < 3; ++axis) {
        const double delta = pts(static_cast<std::size_t>(i), axis) - pts(static_cast<std::size_t>(j), axis);
        sum += delta * delta;
      }
      if (sum > best) {
        best = sum;
        a = i;
        b = j;
      }
    }
  }
  std::vector<int> indices{a, b};
  std::sort(indices.begin(), indices.end());
  return indices;
}

ConvexHull hull3d(const Points& pts) {
  const double scale = coord_scale(pts);
  const double vol_eps = 1e-9 * scale * scale * scale;
  const double area_eps = 1e-9 * scale * scale;
  int i0 = 0;
  for (int i = 1; i < static_cast<int>(pts.n); ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      const double av = pts(static_cast<std::size_t>(i), static_cast<std::size_t>(axis));
      const double bv = pts(static_cast<std::size_t>(i0), static_cast<std::size_t>(axis));
      if (av != bv) {
        if (av < bv) {
          i0 = i;
        }
        break;
      }
    }
  }
  int i1 = -1;
  double best_d = -1.0;
  for (int i = 0; i < static_cast<int>(pts.n); ++i) {
    if (i == i0) {
      continue;
    }
    double sum = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
      const double delta = pts(static_cast<std::size_t>(i), static_cast<std::size_t>(axis)) -
                           pts(static_cast<std::size_t>(i0), static_cast<std::size_t>(axis));
      sum += delta * delta;
    }
    if (sum > best_d) {
      best_d = sum;
      i1 = i;
    }
  }
  int i2 = -1;
  double best_area = area_eps;
  double axis_v[3] = {};
  if (i1 >= 0) {
    for (int axis = 0; axis < 3; ++axis) {
      axis_v[axis] = pts(static_cast<std::size_t>(i1), static_cast<std::size_t>(axis)) -
                     pts(static_cast<std::size_t>(i0), static_cast<std::size_t>(axis));
    }
  }
  for (int i = 0; i < static_cast<int>(pts.n); ++i) {
    if (i == i0 || i == i1) {
      continue;
    }
    double edge[3];
    for (int axis = 0; axis < 3; ++axis) {
      edge[axis] = pts(static_cast<std::size_t>(i), static_cast<std::size_t>(axis)) -
                   pts(static_cast<std::size_t>(i0), static_cast<std::size_t>(axis));
    }
    double normal[3];
    cross3(axis_v, edge, normal);
    const double area = norm3(normal);
    if (area > best_area) {
      best_area = area;
      i2 = i;
    }
  }
  int i3 = -1;
  double best_vol = vol_eps;
  std::array<int, 3> seed{i0, i1, i2};
  if (i2 >= 0) {
    for (int i = 0; i < static_cast<int>(pts.n); ++i) {
      if (i == i0 || i == i1 || i == i2) {
        continue;
      }
      double q[3] = {pts(static_cast<std::size_t>(i), 0), pts(static_cast<std::size_t>(i), 1),
                     pts(static_cast<std::size_t>(i), 2)};
      const double volume = std::abs(orient_point(pts, seed, q));
      if (volume > best_vol) {
        best_vol = volume;
        i3 = i;
      }
    }
  }
  if (i1 < 0) {
    return pack(pts, {0}, {});
  }
  if (i2 < 0 || i3 < 0) {
    if (i2 < 0) {
      return pack(pts, extremes_along_farthest(pts), {});
    }
    Points flat;
    flat.n = pts.n;
    flat.d = 2;
    flat.data.resize(flat.n * 2);
    double ab[3], ac[3], normal[3];
    const double origin[3] = {pts(static_cast<std::size_t>(i0), 0), pts(static_cast<std::size_t>(i0), 1),
                              pts(static_cast<std::size_t>(i0), 2)};
    for (int axis = 0; axis < 3; ++axis) {
      ab[axis] = pts(static_cast<std::size_t>(i1), static_cast<std::size_t>(axis)) - origin[axis];
      ac[axis] = pts(static_cast<std::size_t>(i2), static_cast<std::size_t>(axis)) - origin[axis];
    }
    cross3(ab, ac, normal);
    const double nlen = std::max(norm3(normal), 1e-15);
    const double alen = std::max(norm3(ab), 1e-15);
    for (int axis = 0; axis < 3; ++axis) {
      normal[axis] /= nlen;
      ab[axis] /= alen;
    }
    double yaxis[3];
    cross3(normal, ab, yaxis);
    for (std::size_t i = 0; i < pts.n; ++i) {
      double delta[3];
      for (int axis = 0; axis < 3; ++axis) {
        delta[axis] = pts(i, static_cast<std::size_t>(axis)) - origin[axis];
      }
      flat(i, 0) = dot3(delta, ab);
      flat(i, 1) = dot3(delta, yaxis);
    }
    auto indices = monotone_indices(flat);
    std::sort(indices.begin(), indices.end());
    return pack(pts, std::move(indices), {});
  }

  Mesh mesh;
  mesh.eps = vol_eps;
  mesh.area_eps = area_eps;
  const int tet[4] = {i0, i1, i2, i3};
  for (int axis = 0; axis < 3; ++axis) {
    mesh.inside[axis] = 0.0;
    for (int vertex : tet) {
      mesh.inside[axis] += pts(static_cast<std::size_t>(vertex), static_cast<std::size_t>(axis));
    }
    mesh.inside[axis] /= 4.0;
  }
  const int combos[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
  for (const auto& combo : combos) {
    std::array<int, 3> face{tet[combo[0]], tet[combo[1]], tet[combo[2]]};
    flip_against_inside(pts, face, mesh.inside, mesh.eps);
    mesh.faces.push_back(face);
  }
  for (int i = 0; i < static_cast<int>(pts.n); ++i) {
    if (i == i0 || i == i1 || i == i2 || i == i3) {
      continue;
    }
    add_strict(mesh, pts, i);
  }
  std::vector<char> skip(pts.n, 0);
  skip[static_cast<std::size_t>(i0)] = 1;
  skip[static_cast<std::size_t>(i1)] = 1;
  skip[static_cast<std::size_t>(i2)] = 1;
  skip[static_cast<std::size_t>(i3)] = 1;
  for (int guard = 0; guard < static_cast<int>(pts.n); ++guard) {
    if (!absorb_coplanar(mesh, pts, skip)) {
      break;
    }
  }
  std::vector<int> indices;
  for (const auto& face : mesh.faces) {
    for (int vertex : face) {
      if (std::find(indices.begin(), indices.end(), vertex) == indices.end()) {
        indices.push_back(vertex);
      }
    }
  }
  std::sort(indices.begin(), indices.end());
  return pack(pts, std::move(indices), mesh.faces);
}

}  // namespace

ConvexHull convex_hull(const Points& points, const std::string& method) {
  require_method(method, {"auto", "monotone", "incremental"});
  check_points(points, 1, -1, 2, 3);
  std::string chosen = method;
  if (chosen == "auto") {
    chosen = points.d == 2 ? "monotone" : "incremental";
  }
  if (chosen == "monotone") {
    if (points.d != 2) {
      throw std::invalid_argument("expected points with dimension 2, got " + std::to_string(points.d));
    }
    return pack(points, monotone_indices(points), {});
  }
  if (points.d != 3) {
    throw std::invalid_argument("expected points with dimension 3, got " + std::to_string(points.d));
  }
  return hull3d(points);
}

}  // namespace cg
