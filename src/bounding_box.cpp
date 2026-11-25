#include "cg/bounding_box.hpp"

#include "cg/convex_hull.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace cg {
namespace {

OrientedBox freeze_box(std::vector<double> center, std::vector<double> axes, std::vector<double> extents,
                       std::vector<double> corners) {
  OrientedBox box;
  box.dim = center.size();
  for (double& extent : extents) {
    extent = std::max(extent, 0.0);
  }
  box.volume = 1.0;
  for (double extent : extents) {
    box.volume *= extent;
  }
  box.center = std::move(center);
  box.axes = std::move(axes);
  box.extents = std::move(extents);
  box.corners = std::move(corners);
  return box;
}

std::vector<double> corners_of(const std::vector<double>& center, const std::vector<double>& axes,
                               const std::vector<double>& extents) {
  const std::size_t dim = center.size();
  const std::size_t count = std::size_t{1} << dim;
  std::vector<double> corners(count * dim);
  for (std::size_t mask = 0; mask < count; ++mask) {
    for (std::size_t row = 0; row < dim; ++row) {
      double value = center[row];
      for (std::size_t col = 0; col < dim; ++col) {
        const double sign = ((mask >> col) & 1U) ? 0.5 : -0.5;
        value += axes[col * dim + row] * sign * extents[col];
      }
      corners[mask * dim + row] = value;
    }
  }
  return corners;
}

OrientedBox aabb_in_frame(const Points& pts, const std::vector<double>& origin,
                          const std::vector<double>& axes) {
  const std::size_t dim = pts.d;
  std::vector<double> lo(dim, std::numeric_limits<double>::infinity());
  std::vector<double> hi(dim, -std::numeric_limits<double>::infinity());
  for (std::size_t i = 0; i < pts.n; ++i) {
    for (std::size_t col = 0; col < dim; ++col) {
      double coord = 0.0;
      for (std::size_t row = 0; row < dim; ++row) {
        coord += (pts(i, row) - origin[row]) * axes[col * dim + row];
      }
      lo[col] = std::min(lo[col], coord);
      hi[col] = std::max(hi[col], coord);
    }
  }
  std::vector<double> extents(dim);
  std::vector<double> mid(dim);
  for (std::size_t col = 0; col < dim; ++col) {
    extents[col] = hi[col] - lo[col];
    mid[col] = (lo[col] + hi[col]) / 2.0;
  }
  std::vector<double> center(dim);
  for (std::size_t row = 0; row < dim; ++row) {
    center[row] = origin[row];
    for (std::size_t col = 0; col < dim; ++col) {
      center[row] += axes[col * dim + row] * mid[col];
    }
  }
  auto corners = corners_of(center, axes, extents);
  return freeze_box(std::move(center), axes, std::move(extents), std::move(corners));
}

OrientedBox degenerate(const Points& pts) {
  std::vector<double> center(pts.d, 0.0);
  for (std::size_t i = 0; i < pts.n; ++i) {
    for (std::size_t axis = 0; axis < pts.d; ++axis) {
      center[axis] += pts(i, axis);
    }
  }
  for (double& value : center) {
    value /= static_cast<double>(pts.n);
  }
  std::vector<double> axes(pts.d * pts.d, 0.0);
  for (std::size_t axis = 0; axis < pts.d; ++axis) {
    axes[axis * pts.d + axis] = 1.0;
  }
  std::vector<double> extents(pts.d, 0.0);
  auto corners = corners_of(center, axes, extents);
  return freeze_box(std::move(center), std::move(axes), std::move(extents), std::move(corners));
}

void jacobi(double a[3][3], double v[3][3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      v[i][j] = i == j ? 1.0 : 0.0;
    }
  }
  for (int iter = 0; iter < 32; ++iter) {
    int p = 0;
    int q = 1;
    double largest = std::abs(a[0][1]);
    if (std::abs(a[0][2]) > largest) {
      largest = std::abs(a[0][2]);
      p = 0;
      q = 2;
    }
    if (std::abs(a[1][2]) > largest) {
      largest = std::abs(a[1][2]);
      p = 1;
      q = 2;
    }
    if (largest < 1e-15) {
      break;
    }
    const double app = a[p][p];
    const double aqq = a[q][q];
    const double apq = a[p][q];
    const double tau = (aqq - app) / (2.0 * apq);
    const double t = std::copysign(1.0, tau) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
    const double c = 1.0 / std::sqrt(1.0 + t * t);
    const double s = t * c;
    a[p][p] = app - t * apq;
    a[q][q] = aqq + t * apq;
    a[p][q] = a[q][p] = 0.0;
    for (int r = 0; r < 3; ++r) {
      if (r == p || r == q) {
        continue;
      }
      const double arp = a[r][p];
      const double arq = a[r][q];
      a[r][p] = a[p][r] = c * arp - s * arq;
      a[r][q] = a[q][r] = s * arp + c * arq;
    }
    for (int r = 0; r < 3; ++r) {
      const double vrp = v[r][p];
      const double vrq = v[r][q];
      v[r][p] = c * vrp - s * vrq;
      v[r][q] = s * vrp + c * vrq;
    }
  }
}

OrientedBox pca_box(const Points& pts) {
  std::vector<double> center(pts.d, 0.0);
  for (std::size_t i = 0; i < pts.n; ++i) {
    for (std::size_t axis = 0; axis < pts.d; ++axis) {
      center[axis] += pts(i, axis);
    }
  }
  for (double& value : center) {
    value /= static_cast<double>(pts.n);
  }
  bool flat = pts.n == 1;
  if (!flat) {
    flat = true;
    for (std::size_t i = 0; i < pts.n && flat; ++i) {
      for (std::size_t axis = 0; axis < pts.d; ++axis) {
        if (std::abs(pts(i, axis) - center[axis]) > 1e-12) {
          flat = false;
        }
      }
    }
  }
  if (flat || pts.d != 3) {
    if (pts.d == 2 && !flat) {
      double cov[2][2] = {};
      for (std::size_t i = 0; i < pts.n; ++i) {
        const double x = pts(i, 0) - center[0];
        const double y = pts(i, 1) - center[1];
        cov[0][0] += x * x;
        cov[0][1] += x * y;
        cov[1][1] += y * y;
      }
      cov[1][0] = cov[0][1];
      const double trace = cov[0][0] + cov[1][1];
      const double det = cov[0][0] * cov[1][1] - cov[0][1] * cov[0][1];
      const double gap = std::sqrt(std::max(0.0, trace * trace / 4.0 - det));
      const double lambda = trace / 2.0 + gap;
      double x = cov[0][1];
      double y = lambda - cov[0][0];
      if (std::abs(x) + std::abs(y) < 1e-15) {
        x = 1.0;
        y = 0.0;
      }
      const double length = std::hypot(x, y);
      x /= length;
      y /= length;
      std::vector<double> axes{x, y, -y, x};
      return aabb_in_frame(pts, center, axes);
    }
    return degenerate(pts);
  }
  double cov[3][3] = {};
  for (std::size_t i = 0; i < pts.n; ++i) {
    double delta[3];
    for (int axis = 0; axis < 3; ++axis) {
      delta[axis] = pts(i, static_cast<std::size_t>(axis)) - center[static_cast<std::size_t>(axis)];
    }
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        cov[row][col] += delta[row] * delta[col];
      }
    }
  }
  double basis[3][3];
  jacobi(cov, basis);
  int order[3] = {0, 1, 2};
  std::sort(order, order + 3, [&](int a, int b) { return cov[a][a] > cov[b][b]; });
  std::vector<double> axes(9);
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      axes[static_cast<std::size_t>(col) * 3 + static_cast<std::size_t>(row)] = basis[row][order[col]];
    }
  }
  const double det = axes[0] * (axes[4] * axes[8] - axes[5] * axes[7]) -
                     axes[3] * (axes[1] * axes[8] - axes[2] * axes[7]) +
                     axes[6] * (axes[1] * axes[5] - axes[2] * axes[4]);
  if (det < 0.0) {
    axes[6] *= -1.0;
    axes[7] *= -1.0;
    axes[8] *= -1.0;
  }
  return aabb_in_frame(pts, center, axes);
}

OrientedBox min_area_rect(const Points& pts, const std::vector<double>& hull) {
  const std::size_t count = hull.size() / 2;
  if (count <= 1) {
    return degenerate(pts);
  }
  OrientedBox chosen;
  bool have = false;
  for (std::size_t index = 0; index < count; ++index) {
    const std::size_t next = (index + 1) % count;
    const double ex = hull[next * 2] - hull[index * 2];
    const double ey = hull[next * 2 + 1] - hull[index * 2 + 1];
    const double length = std::hypot(ex, ey);
    if (length <= 1e-15) {
      continue;
    }
    const double xx = ex / length;
    const double xy = ey / length;
    std::vector<double> axes{xx, xy, -xy, xx};
    std::vector<double> origin{hull[index * 2], hull[index * 2 + 1]};
    OrientedBox candidate = aabb_in_frame(pts, origin, axes);
    if (!have || candidate.volume < chosen.volume) {
      chosen = std::move(candidate);
      have = true;
    }
  }
  if (!have) {
    return degenerate(pts);
  }
  return chosen;
}

void cross3(const double a[3], const double b[3], double out[3]) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

OrientedBox min_volume_box(const Points& pts, const std::vector<std::array<int, 3>>& faces) {
  if (faces.empty()) {
    return pca_box(pts);
  }
  bool have = false;
  OrientedBox chosen;
  for (const auto& tri : faces) {
    double a[3], b[3], c[3];
    for (int axis = 0; axis < 3; ++axis) {
      a[axis] = pts(static_cast<std::size_t>(tri[0]), static_cast<std::size_t>(axis));
      b[axis] = pts(static_cast<std::size_t>(tri[1]), static_cast<std::size_t>(axis)) - a[axis];
      c[axis] = pts(static_cast<std::size_t>(tri[2]), static_cast<std::size_t>(axis)) - a[axis];
    }
    double normal[3];
    cross3(b, c, normal);
    const double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    if (length <= 1e-15) {
      continue;
    }
    for (double& value : normal) {
      value /= length;
    }
    const double helper[3] = {std::abs(normal[0]) < 0.9 ? 1.0 : 0.0, std::abs(normal[0]) < 0.9 ? 0.0 : 1.0, 0.0};
    double x_axis[3];
    cross3(helper, normal, x_axis);
    const double xlen = std::sqrt(x_axis[0] * x_axis[0] + x_axis[1] * x_axis[1] + x_axis[2] * x_axis[2]);
    for (double& value : x_axis) {
      value /= xlen;
    }
    double y_axis[3];
    cross3(normal, x_axis, y_axis);
    std::vector<double> frame(9);
    for (int row = 0; row < 3; ++row) {
      frame[row] = x_axis[row];
      frame[3 + row] = y_axis[row];
      frame[6 + row] = normal[row];
    }
    Points flat;
    flat.n = pts.n;
    flat.d = 2;
    flat.data.resize(flat.n * 2);
    std::vector<double> local_z(pts.n);
    for (std::size_t i = 0; i < pts.n; ++i) {
      double delta[3];
      for (int axis = 0; axis < 3; ++axis) {
        delta[axis] = pts(i, static_cast<std::size_t>(axis)) - a[axis];
      }
      for (int col = 0; col < 2; ++col) {
        double coord = 0.0;
        for (int row = 0; row < 3; ++row) {
          coord += delta[row] * frame[static_cast<std::size_t>(col) * 3 + static_cast<std::size_t>(row)];
        }
        flat(i, static_cast<std::size_t>(col)) = coord;
      }
      local_z[i] = delta[0] * frame[6] + delta[1] * frame[7] + delta[2] * frame[8];
    }
    const ConvexHull hull = convex_hull(flat);
    const OrientedBox rect = min_area_rect(flat, hull.vertices);
    const double z_lo = *std::min_element(local_z.begin(), local_z.end());
    const double z_hi = *std::max_element(local_z.begin(), local_z.end());
    const double center_local[3] = {rect.center[0], rect.center[1], (z_lo + z_hi) / 2.0};
    std::vector<double> center(3);
    for (int row = 0; row < 3; ++row) {
      center[static_cast<std::size_t>(row)] =
          a[row] + frame[row] * center_local[0] + frame[3 + row] * center_local[1] + frame[6 + row] * center_local[2];
    }
    std::vector<double> world(9);
    const double local_axes[2][3] = {{rect.axis(0, 0), rect.axis(0, 1), 0.0},
                                     {rect.axis(1, 0), rect.axis(1, 1), 0.0}};
    for (int col = 0; col < 2; ++col) {
      for (int row = 0; row < 3; ++row) {
        world[static_cast<std::size_t>(col) * 3 + static_cast<std::size_t>(row)] =
            frame[row] * local_axes[col][0] + frame[3 + row] * local_axes[col][1] + frame[6 + row] * local_axes[col][2];
      }
    }
    world[6] = normal[0];
    world[7] = normal[1];
    world[8] = normal[2];
    std::vector<double> extents{rect.extents[0], rect.extents[1], z_hi - z_lo};
    auto corners = corners_of(center, world, extents);
    OrientedBox candidate = freeze_box(std::move(center), std::move(world), std::move(extents), std::move(corners));
    if (!have || candidate.volume < chosen.volume) {
      chosen = std::move(candidate);
      have = true;
    }
  }
  if (!have) {
    return pca_box(pts);
  }
  return chosen;
}

}  // namespace

bool OrientedBox::contains(const Points& points, double tol) const {
  if (points.d != dim) {
    return false;
  }
  for (std::size_t i = 0; i < points.n; ++i) {
    for (std::size_t col = 0; col < dim; ++col) {
      double coord = 0.0;
      for (std::size_t row = 0; row < dim; ++row) {
        coord += (points(i, row) - center[row]) * axis(col, row);
      }
      if (std::abs(coord) > extents[col] / 2.0 + tol) {
        return false;
      }
    }
  }
  return true;
}

OrientedBox minimum_bounding_box(const Points& points, const std::string& method) {
  require_method(method, {"auto", "calipers"});
  check_points(points, 1, -1, 2, 3);
  if (points.d == 2) {
    const ConvexHull hull = convex_hull(points);
    return min_area_rect(points, hull.vertices);
  }
  const ConvexHull hull = convex_hull(points);
  return min_volume_box(points, hull.faces);
}

}  // namespace cg
