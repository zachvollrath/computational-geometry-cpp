#include "cg/svg.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cg {
namespace {

struct Pen {
  std::ostringstream out;
  int width = 0;
  int height = 0;

  explicit Pen(int w, int h) : width(w), height(h) {
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w << "\" height=\"" << h
        << "\" viewBox=\"0 0 " << w << " " << h << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
  }

  std::string finish() {
    out << "</svg>\n";
    return out.str();
  }
};

std::vector<std::array<double, 2>> project(const std::vector<std::array<double, 2>>& points, double min_x,
                                           double min_y, double max_x, double max_y, int size, int margin) {
  if (size <= 2 * margin) {
    throw std::invalid_argument("size must be larger than twice the margin");
  }
  const double span_x = std::max(max_x - min_x, 1e-12);
  const double span_y = std::max(max_y - min_y, 1e-12);
  const double scale = static_cast<double>(size - 2 * margin) / std::max(span_x, span_y);
  const double used_x = span_x * scale;
  const double used_y = span_y * scale;
  const double origin_x = margin + (static_cast<double>(size - 2 * margin) - used_x) / 2.0;
  const double origin_y = margin + (static_cast<double>(size - 2 * margin) - used_y) / 2.0;
  std::vector<std::array<double, 2>> pixels;
  pixels.reserve(points.size());
  for (const auto& point : points) {
    const double x = origin_x + (point[0] - min_x) * scale;
    const double y = origin_y + (used_y - (point[1] - min_y) * scale);
    pixels.push_back({x, y});
  }
  return pixels;
}

std::vector<std::array<double, 2>> project_cloud(const std::vector<std::array<double, 2>>& points, int size,
                                                 int margin) {
  double min_x = points[0][0];
  double min_y = points[0][1];
  double max_x = min_x;
  double max_y = min_y;
  for (const auto& point : points) {
    min_x = std::min(min_x, point[0]);
    min_y = std::min(min_y, point[1]);
    max_x = std::max(max_x, point[0]);
    max_y = std::max(max_y, point[1]);
  }
  return project(points, min_x, min_y, max_x, max_y, size, margin);
}

std::vector<std::array<double, 2>> rows2(const Points& points) {
  if (points.d != 2) {
    throw std::invalid_argument("visualization expects points of shape (n, 2)");
  }
  std::vector<std::array<double, 2>> rows(points.n);
  for (std::size_t i = 0; i < points.n; ++i) {
    rows[i] = {points(i, 0), points(i, 1)};
  }
  return rows;
}

void label(Pen& pen, int, const std::string& text, int x) {
  pen.out << "<text x=\"" << x << "\" y=\"" << (pen.height - 16)
          << "\" font-family=\"sans-serif\" font-size=\"16\" fill=\"#141414\">" << text << "</text>\n";
}

void dots(Pen& pen, const std::vector<std::array<double, 2>>& pixels, const char* color, double radius) {
  for (const auto& pixel : pixels) {
    pen.out << "<circle cx=\"" << pixel[0] << "\" cy=\"" << pixel[1] << "\" r=\"" << radius
            << "\" fill=\"" << color << "\"/>\n";
  }
}

}  // namespace

std::string render_closest_pair(const Points& points, const ClosestPair& pair, int size, int margin) {
  const auto pixels = project_cloud(rows2(points), size, margin);
  Pen pen(size, size);
  dots(pen, pixels, "#282828", 4);
  pen.out << "<line x1=\"" << pixels[pair.index_a][0] << "\" y1=\"" << pixels[pair.index_a][1]
          << "\" x2=\"" << pixels[pair.index_b][0] << "\" y2=\"" << pixels[pair.index_b][1]
          << "\" stroke=\"#d22828\" stroke-width=\"2\"/>\n";
  for (std::size_t index : {pair.index_a, pair.index_b}) {
    pen.out << "<circle cx=\"" << pixels[index][0] << "\" cy=\"" << pixels[index][1]
            << "\" r=\"7\" fill=\"#d22828\"/>\n";
  }
  std::ostringstream text;
  text << "d = " << pair.distance << "  (" << pair.index_a << ", " << pair.index_b << ")";
  label(pen, margin, text.str(), margin);
  return pen.finish();
}

std::string render_convex_hull(const Points& points, const ConvexHull& hull, int size, int margin) {
  const auto pixels = project_cloud(rows2(points), size, margin);
  Pen pen(size, size);
  dots(pen, pixels, "#282828", 4);
  pen.out << "<polygon points=\"";
  for (int index : hull.indices) {
    pen.out << pixels[static_cast<std::size_t>(index)][0] << ","
            << pixels[static_cast<std::size_t>(index)][1] << " ";
  }
  pen.out << "\" fill=\"none\" stroke=\"#145ab4\" stroke-width=\"2\"/>\n";
  for (int index : hull.indices) {
    pen.out << "<circle cx=\"" << pixels[static_cast<std::size_t>(index)][0] << "\" cy=\""
            << pixels[static_cast<std::size_t>(index)][1] << "\" r=\"6\" fill=\"#145ab4\"/>\n";
  }
  label(pen, margin, "hull vertices = " + std::to_string(hull.indices.size()), margin);
  return pen.finish();
}

std::string render_segments(const std::vector<Segment>& segments, const std::vector<SegmentIntersection>& hits,
                            int size, int margin) {
  std::vector<std::array<double, 2>> cloud;
  for (const Segment& segment : segments) {
    cloud.push_back(segment[0]);
    cloud.push_back(segment[1]);
  }
  for (const auto& hit : hits) {
    cloud.push_back({hit.x, hit.y});
  }
  const auto pixels = project_cloud(cloud, size, margin);
  Pen pen(size, size);
  for (std::size_t i = 0; i < segments.size(); ++i) {
    pen.out << "<line x1=\"" << pixels[2 * i][0] << "\" y1=\"" << pixels[2 * i][1] << "\" x2=\""
            << pixels[2 * i + 1][0] << "\" y2=\"" << pixels[2 * i + 1][1]
            << "\" stroke=\"#3c3c3c\" stroke-width=\"2\"/>\n";
  }
  for (std::size_t i = 2 * segments.size(); i < pixels.size(); ++i) {
    pen.out << "<circle cx=\"" << pixels[i][0] << "\" cy=\"" << pixels[i][1]
            << "\" r=\"6\" fill=\"#d22828\"/>\n";
  }
  label(pen, margin, "intersections = " + std::to_string(hits.size()), margin);
  return pen.finish();
}

std::string render_bounding_box(const Points& points, const OrientedBox& box, int size, int margin) {
  auto cloud = rows2(points);
  const std::size_t corner_count = box.corners.size() / box.dim;
  for (std::size_t i = 0; i < corner_count; ++i) {
    cloud.push_back({box.corners[i * box.dim], box.corners[i * box.dim + 1]});
  }
  const auto pixels = project_cloud(cloud, size, margin);
  Pen pen(size, size);
  std::vector<std::array<double, 2>> point_pixels(pixels.begin(), pixels.begin() + static_cast<std::ptrdiff_t>(points.n));
  dots(pen, point_pixels, "#282828", 4);
  std::vector<std::size_t> order(corner_count);
  for (std::size_t i = 0; i < corner_count; ++i) {
    order[i] = i;
  }
  double cx = 0.0;
  double cy = 0.0;
  for (std::size_t i = 0; i < corner_count; ++i) {
    cx += box.corners[i * box.dim];
    cy += box.corners[i * box.dim + 1];
  }
  cx /= static_cast<double>(corner_count);
  cy /= static_cast<double>(corner_count);
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    const double aa = std::atan2(box.corners[a * box.dim + 1] - cy, box.corners[a * box.dim] - cx);
    const double bb = std::atan2(box.corners[b * box.dim + 1] - cy, box.corners[b * box.dim] - cx);
    return aa < bb;
  });
  pen.out << "<polygon points=\"";
  for (std::size_t index : order) {
    pen.out << pixels[points.n + index][0] << "," << pixels[points.n + index][1] << " ";
  }
  pen.out << "\" fill=\"none\" stroke=\"#148c28\" stroke-width=\"2\"/>\n";
  std::ostringstream text;
  text << "area = " << box.volume;
  label(pen, margin, text.str(), margin);
  return pen.finish();
}

std::string render_neighbors(const Points& points, const NeighborResult& neighbors, int size, int margin) {
  const auto pixels = project_cloud(rows2(points), size, margin);
  if (neighbors.queries != points.n) {
    throw std::invalid_argument("neighbor rows must match the point count");
  }
  Pen pen(size, size);
  std::vector<std::pair<std::size_t, std::size_t>> seen;
  for (std::size_t i = 0; i < points.n; ++i) {
    for (std::size_t slot = 0; slot < neighbors.k; ++slot) {
      std::size_t j = neighbors.index(i, slot);
      if (i == j) {
        continue;
      }
      const auto edge = std::minmax(i, j);
      if (std::find(seen.begin(), seen.end(), edge) != seen.end()) {
        continue;
      }
      seen.push_back(edge);
      pen.out << "<line x1=\"" << pixels[i][0] << "\" y1=\"" << pixels[i][1] << "\" x2=\"" << pixels[j][0]
              << "\" y2=\"" << pixels[j][1] << "\" stroke=\"#3c96be\" stroke-width=\"1\"/>\n";
    }
  }
  dots(pen, pixels, "#282828", 4);
  label(pen, margin, "k = " + std::to_string(neighbors.k), margin);
  return pen.finish();
}

std::string render_polygon(const Points& polygon, const PolygonArea& area, int size, int margin) {
  Points ring = polygon;
  if (ring.n >= 2 && std::hypot(ring(0, 0) - ring(ring.n - 1, 0), ring(0, 1) - ring(ring.n - 1, 1)) <= 1e-12) {
    --ring.n;
    ring.data.resize(ring.n * ring.d);
  }
  const auto pixels = project_cloud(rows2(ring), size, margin);
  Pen pen(size, size);
  pen.out << "<polygon points=\"";
  for (const auto& pixel : pixels) {
    pen.out << pixel[0] << "," << pixel[1] << " ";
  }
  pen.out << "\" fill=\"#c4e8d6\" stroke=\"#286e28\" stroke-width=\"2\"/>\n";
  dots(pen, pixels, "#1e461e", 5);
  std::ostringstream text;
  text << "area = " << area.area;
  label(pen, margin, text.str(), margin);
  return pen.finish();
}

std::string render_delaunay(const Points& points, const Triangulation& mesh, int size, int margin) {
  const auto pixels = project_cloud(rows2(points), size, margin);
  Pen pen(size, size);
  for (const auto& triangle : mesh.triangles) {
    pen.out << "<polygon points=\"";
    for (int index : triangle) {
      pen.out << pixels[static_cast<std::size_t>(index)][0] << ","
              << pixels[static_cast<std::size_t>(index)][1] << " ";
    }
    pen.out << "\" fill=\"#cee1e6\" stroke=\"#284650\" stroke-width=\"1\"/>\n";
  }
  dots(pen, pixels, "#282828", 4);
  label(pen, margin, "triangles = " + std::to_string(mesh.triangles.size()), margin);
  return pen.finish();
}

std::string render_smoothing(const Points& original, const Points& smoothed,
                             const std::vector<std::array<int, 3>>& faces, int size, int margin) {
  if (original.d != 3 || smoothed.n != original.n || smoothed.d != 3) {
    throw std::invalid_argument("smoothing visualization expects matching (n, 3) meshes");
  }
  double min_x = original(0, 0);
  double min_y = original(0, 1);
  double max_x = min_x;
  double max_y = min_y;
  double min_z = std::min(original(0, 2), smoothed(0, 2));
  double max_z = std::max(original(0, 2), smoothed(0, 2));
  for (const Points* cloud : {&original, &smoothed}) {
    for (std::size_t i = 0; i < cloud->n; ++i) {
      min_x = std::min(min_x, (*cloud)(i, 0));
      min_y = std::min(min_y, (*cloud)(i, 1));
      max_x = std::max(max_x, (*cloud)(i, 0));
      max_y = std::max(max_y, (*cloud)(i, 1));
      min_z = std::min(min_z, (*cloud)(i, 2));
      max_z = std::max(max_z, (*cloud)(i, 2));
    }
  }
  const double span = std::max(max_z - min_z, 1e-12);
  Pen pen(size * 2, size);
  const char* titles[2] = {"original", "smoothed"};
  const Points* clouds[2] = {&original, &smoothed};
  for (int panel = 0; panel < 2; ++panel) {
    std::vector<std::array<double, 2>> xy(clouds[panel]->n);
    for (std::size_t i = 0; i < clouds[panel]->n; ++i) {
      xy[i] = {(*clouds[panel])(i, 0), (*clouds[panel])(i, 1)};
    }
    auto pixels = project(xy, min_x, min_y, max_x, max_y, size, margin);
    for (auto& pixel : pixels) {
      pixel[0] += panel * size;
    }
    for (const auto& face : faces) {
      double height = 0.0;
      for (int index : face) {
        height += (*clouds[panel])(static_cast<std::size_t>(index), 2);
      }
      height /= 3.0;
      const double t = std::min(1.0, std::max(0.0, (height - min_z) / span));
      const int red = static_cast<int>(255 * t);
      const int blue = static_cast<int>(255 * (1.0 - t));
      pen.out << "<polygon points=\"";
      for (int index : face) {
        pen.out << pixels[static_cast<std::size_t>(index)][0] << ","
                << pixels[static_cast<std::size_t>(index)][1] << " ";
      }
      pen.out << "\" fill=\"rgb(" << red << ",70," << blue << ")\" stroke=\"#282828\" stroke-width=\"1\"/>\n";
    }
    label(pen, margin, titles[panel], margin + panel * size);
  }
  return pen.finish();
}

void save_svg(const std::string& path, const std::string& svg) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("failed to write " + path);
  }
  file << svg;
  if (!file) {
    throw std::runtime_error("failed to write " + path);
  }
}

}  // namespace cg
