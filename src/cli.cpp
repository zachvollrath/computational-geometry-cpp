#include "cg/bounding_box.hpp"
#include "cg/closest_pair.hpp"
#include "cg/convex_hull.hpp"
#include "cg/delaunay.hpp"
#include "cg/laplacian.hpp"
#include "cg/neighbors.hpp"
#include "cg/points.hpp"
#include "cg/segments.hpp"
#include "cg/shoelace.hpp"
#include "cg/svg.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <numbers>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Options {
  std::string example;
  int n = 200;
  int dim = 2;
  int seed = 0;
  int k = 3;
  int iterations = 8;
  std::string out;
  bool has_angle = false;
  double angle = 0.0;
};

void usage() {
  std::cerr
      << "usage: geometry <example> [options]\n"
      << "examples: closest-pair, convex-hull, smooth, intersections, obb, neighbors, shoelace, delaunay\n"
      << "options: --n N  --dim 2|3  --seed N  --k K  --iterations N  --angle DEG  --out FILE.svg\n";
}

bool take_int(int& i, int argc, char** argv, int& slot) {
  if (i + 1 >= argc) {
    return false;
  }
  slot = std::stoi(argv[++i]);
  return true;
}

bool take_double(int& i, int argc, char** argv, double& slot) {
  if (i + 1 >= argc) {
    return false;
  }
  slot = std::stod(argv[++i]);
  return true;
}

bool take_string(int& i, int argc, char** argv, std::string& slot) {
  if (i + 1 >= argc) {
    return false;
  }
  slot = argv[++i];
  return true;
}

Options parse(int argc, char** argv) {
  Options options;
  int index = 1;
  if (argc < 2 || argv[1][0] == '-') {
    options.example = "closest-pair";
  } else {
    options.example = argv[1];
    index = 2;
  }
  if (options.example == "closest-pair") {
    options.n = 200;
  } else if (options.example == "smooth") {
    options.n = 6;
  } else if (options.example == "intersections") {
    options.n = 20;
  } else if (options.example == "shoelace") {
    options.n = 8;
  } else if (options.example == "delaunay") {
    options.n = 40;
  }
  for (int i = index; i < argc; ++i) {
    const std::string flag = argv[i];
    if (flag == "--n" && take_int(i, argc, argv, options.n)) {
      continue;
    }
    if (flag == "--dim" && take_int(i, argc, argv, options.dim)) {
      continue;
    }
    if (flag == "--seed" && take_int(i, argc, argv, options.seed)) {
      continue;
    }
    if (flag == "--k" && take_int(i, argc, argv, options.k)) {
      continue;
    }
    if (flag == "--iterations" && take_int(i, argc, argv, options.iterations)) {
      continue;
    }
    if (flag == "--angle" && take_double(i, argc, argv, options.angle)) {
      options.has_angle = true;
      continue;
    }
    if (flag == "--out" && take_string(i, argc, argv, options.out)) {
      continue;
    }
    if (flag == "--help" || flag == "-h") {
      usage();
      std::exit(0);
    }
    std::cerr << "unknown argument " << flag << "\n";
    usage();
    std::exit(1);
  }
  return options;
}

cg::Points cloud(const Options& options) {
  if (options.n < 2) {
    throw std::invalid_argument("--n must be at least 2");
  }
  std::mt19937 rng(static_cast<std::uint32_t>(options.seed));
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  cg::Points points;
  points.n = static_cast<std::size_t>(options.n);
  points.d = static_cast<std::size_t>(options.dim);
  points.data.resize(points.n * points.d);
  for (double& value : points.data) {
    value = uniform(rng);
  }
  return points;
}

std::pair<cg::Points, double> rotated_cloud(const Options& options) {
  if (options.n < 2) {
    throw std::invalid_argument("--n must be at least 2");
  }
  std::mt19937 rng(static_cast<std::uint32_t>(options.seed));
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);
  cg::Points points;
  points.n = static_cast<std::size_t>(options.n);
  points.d = static_cast<std::size_t>(options.dim);
  points.data.resize(points.n * points.d);
  const double scale[3] = {3.0, 1.0, 0.4};
  for (std::size_t i = 0; i < points.n; ++i) {
    for (std::size_t axis = 0; axis < points.d; ++axis) {
      points(i, axis) = (uniform(rng) - 0.5) * scale[axis];
    }
  }
  double angle = options.has_angle ? options.angle * std::numbers::pi / 180.0
                                  : std::uniform_real_distribution<double>(0.2, std::numbers::pi - 0.2)(rng);
  if (points.d == 2) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    for (std::size_t i = 0; i < points.n; ++i) {
      const double x = points(i, 0);
      const double y = points(i, 1);
      points(i, 0) = c * x - s * y;
      points(i, 1) = s * x + c * y;
    }
  } else {
    double rotation[9];
    if (options.has_angle) {
      const double c = std::cos(angle);
      const double s = std::sin(angle);
      rotation[0] = c;
      rotation[1] = -s;
      rotation[2] = 0;
      rotation[3] = s;
      rotation[4] = c;
      rotation[5] = 0;
      rotation[6] = 0;
      rotation[7] = 0;
      rotation[8] = 1;
    } else {
      double axis[3] = {normal(rng), normal(rng), normal(rng)};
      const double length = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
      for (double& value : axis) {
        value /= length;
      }
      const double skew[9] = {0, -axis[2], axis[1], axis[2], 0, -axis[0], -axis[1], axis[0], 0};
      double skew2[9] = {};
      for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
          for (int k = 0; k < 3; ++k) {
            skew2[row * 3 + col] += skew[row * 3 + k] * skew[k * 3 + col];
          }
        }
      }
      const double sine = std::sin(angle);
      const double cosine = 1.0 - std::cos(angle);
      for (int i = 0; i < 9; ++i) {
        rotation[i] = (i % 4 == 0 ? 1.0 : 0.0) + sine * skew[i] + cosine * skew2[i];
      }
    }
    for (std::size_t i = 0; i < points.n; ++i) {
      double in[3] = {points(i, 0), points(i, 1), points(i, 2)};
      for (int row = 0; row < 3; ++row) {
        points(i, static_cast<std::size_t>(row)) =
            rotation[row * 3] * in[0] + rotation[row * 3 + 1] * in[1] + rotation[row * 3 + 2] * in[2];
      }
    }
  }
  return {points, angle};
}

template <typename Fn>
double milliseconds(Fn&& fn) {
  const auto start = std::chrono::steady_clock::now();
  fn();
  const auto stop = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::milli>(stop - start).count();
}

void demo_closest(const Options& options) {
  const cg::Points points = cloud(options);
  std::vector<std::string> methods = {"brute", "kdtree"};
  if (options.dim == 2) {
    methods.insert(methods.begin() + 1, "divide_conquer");
  }
  for (const std::string& method : methods) {
    cg::ClosestPair pair;
    const double elapsed = milliseconds([&] { pair = cg::closest_pair(points, method); });
    std::printf("%-16s  d=%.8f  pair=(%zu, %zu)  %8.2f ms\n", method.c_str(), pair.distance, pair.index_a,
                pair.index_b, elapsed);
    if (!options.out.empty() && options.dim == 2 && method == "divide_conquer") {
      cg::save_svg(options.out, cg::render_closest_pair(points, pair));
      std::cout << "wrote " << options.out << "\n";
    }
  }
}

void demo_hull(const Options& options) {
  const cg::Points points = cloud(options);
  const std::string method = options.dim == 2 ? "monotone" : "incremental";
  cg::ConvexHull hull;
  const double elapsed = milliseconds([&] { hull = cg::convex_hull(points, method); });
  std::printf("%-16s  vertices=%4zu  %8.2f ms\n", method.c_str(), hull.indices.size(), elapsed);
  if (!options.out.empty() && options.dim == 2) {
    cg::save_svg(options.out, cg::render_convex_hull(points, hull));
    std::cout << "wrote " << options.out << "\n";
  }
}

void demo_smooth(const Options& options) {
  if (options.n < 2) {
    throw std::invalid_argument("--n must be at least 2");
  }
  const int n = options.n;
  cg::Points points;
  points.n = static_cast<std::size_t>(n * n);
  points.d = 3;
  points.data.resize(points.n * 3);
  std::mt19937 rng(static_cast<std::uint32_t>(options.seed));
  std::normal_distribution<double> noise(0.0, 0.05);
  for (int y = 0; y < n; ++y) {
    for (int x = 0; x < n; ++x) {
      const std::size_t i = static_cast<std::size_t>(y * n + x);
      points(i, 0) = n == 1 ? 0.0 : static_cast<double>(x) / (n - 1);
      points(i, 1) = n == 1 ? 0.0 : static_cast<double>(y) / (n - 1);
      points(i, 2) = noise(rng);
    }
  }
  std::vector<std::array<int, 3>> faces;
  for (int y = 0; y < n - 1; ++y) {
    for (int x = 0; x < n - 1; ++x) {
      const int i = y * n + x;
      faces.push_back({i, i + 1, i + n});
      faces.push_back({i + 1, i + n + 1, i + n});
    }
  }
  const cg::Points smoothed =
      cg::laplacian_smooth(points, &faces, options.iterations, 0.5, "inverse_distance");
  double shift = 0.0;
  for (std::size_t i = 0; i < points.data.size(); ++i) {
    const double delta = smoothed.data[i] - points.data[i];
    shift += delta * delta;
  }
  shift = std::sqrt(shift) / static_cast<double>(points.n);
  std::printf("grid %dx%d  iterations=%d  mean shift=%.6f\n", n, n, options.iterations, shift);
  if (!options.out.empty()) {
    cg::save_svg(options.out, cg::render_smoothing(points, smoothed, faces));
    std::cout << "wrote " << options.out << "\n";
  }
}

void demo_segments(const Options& options) {
  if (options.n < 2) {
    throw std::invalid_argument("--n must be at least 2");
  }
  std::mt19937 rng(static_cast<std::uint32_t>(options.seed));
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::vector<cg::Segment> segments(static_cast<std::size_t>(options.n));
  for (auto& segment : segments) {
    for (auto& point : segment) {
      point = {uniform(rng), uniform(rng)};
    }
  }
  std::vector<cg::SegmentIntersection> hits;
  const double elapsed = milliseconds([&] { hits = cg::segment_intersections(segments); });
  std::printf("%-16s  intersections=%4zu  %8.2f ms\n", "sweep", hits.size(), elapsed);
  if (!options.out.empty()) {
    cg::save_svg(options.out, cg::render_segments(segments, hits));
    std::cout << "wrote " << options.out << "\n";
  }
}

void demo_box(const Options& options) {
  const auto [points, angle] = rotated_cloud(options);
  std::printf("rotated the point set by %.1f deg\n", angle * 180.0 / std::numbers::pi);
  cg::OrientedBox box;
  const double elapsed = milliseconds([&] { box = cg::minimum_bounding_box(points); });
  std::printf("%-16s  volume=%.6f  %8.2f ms\n", "calipers", box.volume, elapsed);
  if (!options.out.empty() && options.dim == 2) {
    cg::save_svg(options.out, cg::render_bounding_box(points, box));
    std::cout << "wrote " << options.out << "\n";
  }
}

void demo_neighbors(const Options& options) {
  if (options.k < 1 || options.k >= options.n) {
    throw std::invalid_argument("--k must be between 1 and n-1");
  }
  const cg::Points points = cloud(options);
  for (const char* method : {"brute", "kdtree"}) {
    cg::NeighborResult result;
    const double elapsed =
        milliseconds([&] { result = cg::nearest_neighbors(points, static_cast<std::size_t>(options.k), nullptr, method); });
    double mean = 0.0;
    for (std::size_t i = 0; i < points.n; ++i) {
      mean += result.distance(i, 0);
    }
    mean /= static_cast<double>(points.n);
    std::printf("%-16s  mean d=%.6f  %8.2f ms\n", method, mean, elapsed);
    if (!options.out.empty() && options.dim == 2 && std::strcmp(method, "kdtree") == 0) {
      cg::save_svg(options.out, cg::render_neighbors(points, result));
      std::cout << "wrote " << options.out << "\n";
    }
  }
}

void demo_shoelace(const Options& options) {
  if (options.n < 3) {
    throw std::invalid_argument("--n must be at least 3");
  }
  std::mt19937 rng(static_cast<std::uint32_t>(options.seed));
  std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * std::numbers::pi);
  std::uniform_real_distribution<double> radius_dist(0.35, 1.0);
  std::vector<double> angles(static_cast<std::size_t>(options.n));
  for (double& angle : angles) {
    angle = angle_dist(rng);
  }
  std::sort(angles.begin(), angles.end());
  cg::Points polygon;
  polygon.n = angles.size();
  polygon.d = 2;
  polygon.data.resize(polygon.n * 2);
  for (std::size_t i = 0; i < polygon.n; ++i) {
    const double radius = radius_dist(rng);
    polygon(i, 0) = std::cos(angles[i]) * radius;
    polygon(i, 1) = std::sin(angles[i]) * radius;
  }
  const cg::PolygonArea area = cg::shoelace(polygon);
  std::printf("%-16s  area=%.6f  signed=%.6f\n", "shoelace", area.area, area.signed_area);
  if (!options.out.empty()) {
    cg::save_svg(options.out, cg::render_polygon(polygon, area));
    std::cout << "wrote " << options.out << "\n";
  }
}

void demo_delaunay(const Options& options) {
  if (options.n < 3) {
    throw std::invalid_argument("--n must be at least 3");
  }
  Options flat = options;
  flat.dim = 2;
  const cg::Points points = cloud(flat);
  cg::Triangulation mesh;
  const double elapsed = milliseconds([&] { mesh = cg::delaunay(points); });
  std::printf("%-16s  triangles=%4zu  %8.2f ms\n", "bowyer_watson", mesh.triangles.size(), elapsed);
  if (!options.out.empty()) {
    cg::save_svg(options.out, cg::render_delaunay(points, mesh));
    std::cout << "wrote " << options.out << "\n";
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse(argc, argv);
    if (options.example == "closest-pair") {
      demo_closest(options);
    } else if (options.example == "convex-hull") {
      demo_hull(options);
    } else if (options.example == "smooth") {
      demo_smooth(options);
    } else if (options.example == "intersections") {
      demo_segments(options);
    } else if (options.example == "obb") {
      demo_box(options);
    } else if (options.example == "neighbors") {
      demo_neighbors(options);
    } else if (options.example == "shoelace") {
      demo_shoelace(options);
    } else if (options.example == "delaunay") {
      demo_delaunay(options);
    } else {
      usage();
      return 1;
    }
  } catch (const std::exception& exc) {
    std::cerr << exc.what() << "\n";
    return 1;
  }
  return 0;
}
