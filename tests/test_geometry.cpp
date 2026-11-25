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
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace {

int g_failed = 0;

void check(bool cond, const char* expr, int line) {
  if (!cond) {
    std::cerr << "FAIL line " << line << ": " << expr << "\n";
    ++g_failed;
  }
}

#define CHECK(cond) check(static_cast<bool>(cond), #cond, __LINE__)

bool near(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol; }

void test_closest_pair() {
  const cg::Points square = cg::make_points(2, {{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  const auto divide = cg::closest_pair(square, "divide_conquer");
  const auto brute = cg::closest_pair(square, "brute");
  const auto tree = cg::closest_pair(square, "kdtree");
  CHECK(near(divide.distance, 1.0));
  CHECK(near(brute.distance, divide.distance));
  CHECK(near(tree.distance, divide.distance));
  CHECK(divide.index_a < divide.index_b);

  const cg::Points two = cg::make_points(2, {{0, 0}, {3, 4}});
  const auto pair = cg::closest_pair(two);
  CHECK(pair.index_a == 0 && pair.index_b == 1);
  CHECK(near(pair.distance, 5.0));

  const cg::Points dup = cg::make_points(2, {{1, 2}, {4, 5}, {1, 2}, {9, 0}});
  for (const char* method : {"brute", "divide_conquer", "kdtree"}) {
    const auto found = cg::closest_pair(dup, method);
    CHECK(near(found.distance, 0.0));
    CHECK(found.index_a == 0 && found.index_b == 2);
  }

  std::mt19937 rng(7);
  std::normal_distribution<double> normal;
  cg::Points cloud;
  cloud.n = 80;
  cloud.d = 2;
  cloud.data.resize(cloud.n * cloud.d);
  for (double& value : cloud.data) {
    value = normal(rng);
  }
  const auto a = cg::closest_pair(cloud, "brute");
  const auto b = cg::closest_pair(cloud, "divide_conquer");
  const auto c = cg::closest_pair(cloud, "kdtree");
  CHECK(near(a.distance, b.distance));
  CHECK(near(a.distance, c.distance));
  CHECK(near(a.distance, cg::distance(cloud, a.index_a, a.index_b)));

  cg::Points volume = cloud;
  volume.d = 3;
  volume.n = 40;
  volume.data.resize(volume.n * 3);
  for (double& value : volume.data) {
    value = normal(rng);
  }
  CHECK(near(cg::closest_pair(volume, "brute").distance, cg::closest_pair(volume, "kdtree").distance));
  bool threw = false;
  try {
    cg::closest_pair(volume, "divide_conquer");
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);

  const cg::Points line = cg::make_points(2, {{0, 0}, {2, 0}, {5, 0}, {9, 0}});
  CHECK(near(cg::closest_pair(line, "divide_conquer").distance, 2.0));

  const cg::Points cluster =
      cg::make_points(2, {{0, 0}, {10, 10}, {10, 10.05}, {-4, 3}, {8, -6}});
  const auto close = cg::closest_pair(cluster, "divide_conquer");
  CHECK(close.index_a == 1 && close.index_b == 2);
  CHECK(near(close.distance, 0.05));

  threw = false;
  try {
    cg::closest_pair(cg::make_points(2, {{0, 0}}));
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
  const std::string picture = cg::render_closest_pair(square, divide, 200, 20);
  CHECK(picture.find("width=\"200\"") != std::string::npos);
  CHECK(picture.find("svg") != std::string::npos);
}

void test_hull() {
  const cg::Points square = cg::make_points(2, {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.25, 0.25}});
  const auto hull = cg::convex_hull(square);
  CHECK(hull.indices == std::vector<int>({0, 1, 2, 3}));

  const cg::Points line = cg::make_points(2, {{0, 0}, {2, 0}, {5, 0}, {9, 0}});
  const auto ends = cg::convex_hull(line);
  CHECK(ends.indices == std::vector<int>({0, 3}));

  std::mt19937 rng(3);
  std::normal_distribution<double> normal;
  cg::Points cloud;
  cloud.n = 60;
  cloud.d = 2;
  cloud.data.resize(cloud.n * 2);
  for (double& value : cloud.data) {
    value = normal(rng);
  }
  const auto random_hull = cg::convex_hull(cloud);
  CHECK(random_hull.indices.size() >= 3);
  for (std::size_t i = 0; i < cloud.n; ++i) {
    bool inside = true;
    for (std::size_t e = 0; e < random_hull.indices.size(); ++e) {
      const int a = random_hull.indices[e];
      const int b = random_hull.indices[(e + 1) % random_hull.indices.size()];
      const double cross = (cloud(static_cast<std::size_t>(b), 0) - cloud(static_cast<std::size_t>(a), 0)) *
                               (cloud(i, 1) - cloud(static_cast<std::size_t>(a), 1)) -
                           (cloud(static_cast<std::size_t>(b), 1) - cloud(static_cast<std::size_t>(a), 1)) *
                               (cloud(i, 0) - cloud(static_cast<std::size_t>(a), 0));
      if (cross < -1e-8) {
        inside = false;
      }
    }
    CHECK(inside);
  }

  const cg::Points tet = cg::make_points(3, {{0, 0, 0}, {1, 0, 0}, {0.2, 1, 0}, {0.2, 0.3, 1}, {0.3, 0.3, 0.2}});
  const auto solid = cg::convex_hull(tet);
  CHECK(std::set<int>(solid.indices.begin(), solid.indices.end()) == std::set<int>({0, 1, 2, 3}));
  bool threw = false;
  try {
    cg::convex_hull(tet, "monotone");
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
  CHECK(cg::render_convex_hull(square, hull, 200, 20).find("svg") != std::string::npos);
}

std::vector<int> brute_hull_vertices(const cg::Points& pts) {
  double scale = 1.0;
  for (double value : pts.data) {
    scale = std::max(scale, std::abs(value));
  }
  const double eps = 1e-8 * scale * scale * scale;
  const double area_eps = 1e-8 * scale * scale;
  std::vector<char> on(pts.n, 0);
  auto volume = [&](int a, int b, int c, int d) {
    const double ab[3] = {pts(static_cast<std::size_t>(b), 0) - pts(static_cast<std::size_t>(a), 0),
                          pts(static_cast<std::size_t>(b), 1) - pts(static_cast<std::size_t>(a), 1),
                          pts(static_cast<std::size_t>(b), 2) - pts(static_cast<std::size_t>(a), 2)};
    const double ac[3] = {pts(static_cast<std::size_t>(c), 0) - pts(static_cast<std::size_t>(a), 0),
                          pts(static_cast<std::size_t>(c), 1) - pts(static_cast<std::size_t>(a), 1),
                          pts(static_cast<std::size_t>(c), 2) - pts(static_cast<std::size_t>(a), 2)};
    const double ad[3] = {pts(static_cast<std::size_t>(d), 0) - pts(static_cast<std::size_t>(a), 0),
                          pts(static_cast<std::size_t>(d), 1) - pts(static_cast<std::size_t>(a), 1),
                          pts(static_cast<std::size_t>(d), 2) - pts(static_cast<std::size_t>(a), 2)};
    const double n0 = ab[1] * ac[2] - ab[2] * ac[1];
    const double n1 = ab[2] * ac[0] - ab[0] * ac[2];
    const double n2 = ab[0] * ac[1] - ab[1] * ac[0];
    return n0 * ad[0] + n1 * ad[1] + n2 * ad[2];
  };
  for (int i = 0; i < static_cast<int>(pts.n); ++i) {
    for (int j = i + 1; j < static_cast<int>(pts.n); ++j) {
      for (int k = j + 1; k < static_cast<int>(pts.n); ++k) {
        const double n0 = (pts(static_cast<std::size_t>(j), 1) - pts(static_cast<std::size_t>(i), 1)) *
                              (pts(static_cast<std::size_t>(k), 2) - pts(static_cast<std::size_t>(i), 2)) -
                          (pts(static_cast<std::size_t>(j), 2) - pts(static_cast<std::size_t>(i), 2)) *
                              (pts(static_cast<std::size_t>(k), 1) - pts(static_cast<std::size_t>(i), 1));
        const double n1 = (pts(static_cast<std::size_t>(j), 2) - pts(static_cast<std::size_t>(i), 2)) *
                              (pts(static_cast<std::size_t>(k), 0) - pts(static_cast<std::size_t>(i), 0)) -
                          (pts(static_cast<std::size_t>(j), 0) - pts(static_cast<std::size_t>(i), 0)) *
                              (pts(static_cast<std::size_t>(k), 2) - pts(static_cast<std::size_t>(i), 2));
        const double n2 = (pts(static_cast<std::size_t>(j), 0) - pts(static_cast<std::size_t>(i), 0)) *
                              (pts(static_cast<std::size_t>(k), 1) - pts(static_cast<std::size_t>(i), 1)) -
                          (pts(static_cast<std::size_t>(j), 1) - pts(static_cast<std::size_t>(i), 1)) *
                              (pts(static_cast<std::size_t>(k), 0) - pts(static_cast<std::size_t>(i), 0));
        if (n0 * n0 + n1 * n1 + n2 * n2 <= area_eps * area_eps) {
          continue;
        }
        bool pos = false;
        bool neg = false;
        for (int t = 0; t < static_cast<int>(pts.n); ++t) {
          if (t == i || t == j || t == k) {
            continue;
          }
          const double vol = volume(i, j, k, t);
          pos = pos || vol > eps;
          neg = neg || vol < -eps;
        }
        if (!(pos && neg)) {
          on[static_cast<std::size_t>(i)] = 1;
          on[static_cast<std::size_t>(j)] = 1;
          on[static_cast<std::size_t>(k)] = 1;
        }
      }
    }
  }
  std::vector<int> indices;
  for (int i = 0; i < static_cast<int>(pts.n); ++i) {
    if (on[static_cast<std::size_t>(i)]) {
      indices.push_back(i);
    }
  }
  return indices;
}

void test_hull_3d_matches_supporting_planes() {
  std::mt19937 rng(12);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  cg::Points cloud;
  cloud.n = 18;
  cloud.d = 3;
  cloud.data.resize(cloud.n * 3);
  for (double& value : cloud.data) {
    value = uniform(rng);
  }
  const auto hull = cg::convex_hull(cloud);
  const auto brute = brute_hull_vertices(cloud);
  CHECK(std::set<int>(hull.indices.begin(), hull.indices.end()) == std::set<int>(brute.begin(), brute.end()));
}

void test_box() {
  const double theta = 27.0 * std::acos(-1.0) / 180.0;
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  cg::Points points = cg::make_points(2, {{0, 0}, {3, 0}, {3, 1}, {0, 1}, {1, 0.4}});
  for (std::size_t i = 0; i < points.n; ++i) {
    const double x = points(i, 0);
    const double y = points(i, 1);
    points(i, 0) = c * x - s * y;
    points(i, 1) = s * x + c * y;
  }
  const auto box = cg::minimum_bounding_box(points);
  CHECK(near(box.volume, 3.0, 1e-6));
  CHECK(box.contains(points));
  CHECK(cg::render_bounding_box(points, box, 200, 20).find("svg") != std::string::npos);

  cg::Points cube;
  cube.d = 3;
  for (int x = 0; x < 2; ++x) {
    for (int y = 0; y < 2; ++y) {
      for (int z = 0; z < 2; ++z) {
        cube.data.insert(cube.data.end(), {x - 0.5, y - 0.5, z - 0.5});
      }
    }
  }
  cube.data.insert(cube.data.end(), {0, 0, 0});
  cube.n = cube.data.size() / 3;
  const double angle = 0.4;
  const double cc = std::cos(angle);
  const double ss = std::sin(angle);
  for (std::size_t i = 0; i < cube.n; ++i) {
    const double x = cube(i, 0);
    const double y = cube(i, 1);
    cube(i, 0) = cc * x - ss * y;
    cube(i, 1) = ss * x + cc * y;
  }
  const auto solid = cg::minimum_bounding_box(cube);
  CHECK(near(solid.volume, 1.0, 1e-3));
  CHECK(solid.contains(cube, 1e-4));

  const auto point = cg::minimum_bounding_box(cg::make_points(2, {{2, -3}}));
  CHECK(near(point.volume, 0.0));
  CHECK(near(point.center[0], 2.0) && near(point.center[1], -3.0));
}

void test_neighbors() {
  std::mt19937 rng(9);
  std::uniform_real_distribution<double> uniform;
  cg::Points points;
  points.n = 35;
  points.d = 3;
  points.data.resize(points.n * 3);
  for (double& value : points.data) {
    value = uniform(rng);
  }
  const auto brute = cg::nearest_neighbors(points, 4, nullptr, "brute");
  const auto tree = cg::nearest_neighbors(points, 4, nullptr, "kdtree");
  for (std::size_t i = 0; i < brute.distances.size(); ++i) {
    CHECK(near(brute.distances[i], tree.distances[i], 1e-9));
    CHECK(brute.indices[i] == tree.indices[i]);
  }

  const cg::Points plane = cg::make_points(2, {{0, 0}, {1, 0}, {0, 2}});
  const auto found = cg::nearest_neighbors(plane, 1);
  CHECK(found.index(0, 0) == 1);
  CHECK(found.index(1, 0) == 0);
  CHECK(found.index(2, 0) == 0);
  CHECK(near(found.distance(0, 0), 1.0));
  CHECK(near(found.distance(2, 0), 2.0));

  const cg::Points queries = cg::make_points(3, {{0.1, 0, 0}});
  const cg::Points data = cg::make_points(3, {{0, 0, 0}, {1, 0, 0}});
  const auto explicit_query = cg::nearest_neighbors(data, 1, &queries, "brute");
  CHECK(explicit_query.index(0, 0) == 0);
  CHECK(near(explicit_query.distance(0, 0), 0.1));

  bool threw = false;
  try {
    cg::nearest_neighbors(cg::make_points(2, {{0, 0}, {1, 0}, {2, 0}}), 3);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
  CHECK(cg::render_neighbors(plane, found, 200, 20).find("svg") != std::string::npos);
}

void test_shoelace() {
  const auto square = cg::shoelace(cg::make_points(2, {{0, 0}, {1, 0}, {1, 1}, {0, 1}}));
  CHECK(near(square.signed_area, 1.0));
  CHECK(near(square.area, 1.0));
  const auto clockwise = cg::shoelace(cg::make_points(2, {{0, 0}, {0, 1}, {1, 1}, {1, 0}}));
  CHECK(near(clockwise.signed_area, -1.0));
  CHECK(near(clockwise.area, 1.0));
  const auto closed = cg::shoelace(cg::make_points(2, {{0, 0}, {2, 0}, {2, 3}, {0, 3}, {0, 0}}));
  CHECK(near(closed.area, 6.0));
  const cg::Points triangle = cg::make_points(2, {{0, 0}, {1, 0}, {0, 1}});
  CHECK(cg::render_polygon(triangle, cg::shoelace(triangle), 200, 20).find("svg") != std::string::npos);
}

std::vector<std::array<double, 2>> pairwise(const std::vector<cg::Segment>& segments) {
  std::vector<std::array<double, 2>> found;
  for (std::size_t i = 0; i < segments.size(); ++i) {
    for (std::size_t j = i + 1; j < segments.size(); ++j) {
      const auto& a = segments[i][0];
      const auto& b = segments[i][1];
      const auto& c = segments[j][0];
      const auto& d = segments[j][1];
      const double denom = (b[0] - a[0]) * (d[1] - c[1]) - (b[1] - a[1]) * (d[0] - c[0]);
      if (std::abs(denom) <= 1e-12) {
        continue;
      }
      const double t = ((c[0] - a[0]) * (d[1] - c[1]) - (c[1] - a[1]) * (d[0] - c[0])) / denom;
      const double u = ((c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0])) / denom;
      if (0.0 <= t && t <= 1.0 && 0.0 <= u && u <= 1.0) {
        found.push_back({a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])});
      }
    }
  }
  return found;
}

void test_segments() {
  const std::vector<cg::Segment> cross{{{{0, 0}, {1, 1}}}, {{{0, 1}, {1, 0}}}};
  const auto hits = cg::segment_intersections(cross);
  CHECK(hits.size() == 1);
  CHECK(near(hits[0].x, 0.5) && near(hits[0].y, 0.5));
  CHECK(hits[0].segment_a == 0 && hits[0].segment_b == 1);

  const auto shared = cg::segment_intersections({{{{0, 0}, {1, 0}}}, {{{1, 0}, {1, 1}}}});
  CHECK(shared.size() == 1);
  CHECK(near(shared[0].x, 1.0) && near(shared[0].y, 0.0));

  const auto overlap = cg::segment_intersections({{{{0, 0}, {4, 0}}}, {{{2, 0}, {6, 0}}}});
  CHECK(overlap.size() == 2);
  CHECK(near(overlap[0].x, 2.0) && near(overlap[0].y, 0.0));
  CHECK(near(overlap[1].x, 4.0) && near(overlap[1].y, 0.0));

  const auto vertical = cg::segment_intersections({{{{0, -1}, {0, 1}}}, {{{-1, 0}, {1, 0}}}});
  CHECK(vertical.size() == 1);
  CHECK(near(vertical[0].x, 0.0) && near(vertical[0].y, 0.0));

  std::mt19937 rng(5);
  std::normal_distribution<double> normal;
  std::vector<cg::Segment> random(25);
  for (auto& segment : random) {
    for (auto& point : segment) {
      point = {normal(rng), normal(rng)};
    }
  }
  const auto swept = cg::segment_intersections(random);
  auto expected = pairwise(random);
  CHECK(swept.size() == expected.size());
  for (const auto& hit : swept) {
    const auto match = std::find_if(expected.begin(), expected.end(), [&](const auto& point) {
      return near(point[0], hit.x, 1e-6) && near(point[1], hit.y, 1e-6);
    });
    CHECK(match != expected.end());
    if (match != expected.end()) {
      expected.erase(match);
    }
  }

  const auto none = cg::segment_intersections({{{{0, 0}, {1, 0}}}, {{{0, 1}, {1, 1}}}, {{{2, 0}, {3, 1}}}});
  CHECK(none.empty());
  CHECK(cg::render_segments(cross, hits, 200, 20).find("svg") != std::string::npos);
}

void test_laplacian_and_delaunay() {
  const cg::Points bump = cg::make_points(3, {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0.5, 0.5, 0.4}});
  const std::vector<std::array<int, 3>> faces{{0, 1, 4}, {1, 3, 4}, {3, 2, 4}, {2, 0, 4}};
  const auto smoothed = cg::laplacian_smooth(bump, &faces, 3, 0.4);
  const std::vector<double> expected{0.392, 0.392, 0.08983703703703703, 0.608, 0.392, 0.08983703703703703,
                                     0.392, 0.608, 0.08983703703703703, 0.608, 0.608, 0.08983703703703703,
                                     0.5,   0.5,   0.13048888888888888};
  CHECK(smoothed.data.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    CHECK(near(smoothed.data[i], expected[i], 1e-12));
  }
  const auto inverse = cg::laplacian_smooth(bump, &faces, 2, 0.5, "inverse_distance");
  CHECK(near(inverse(4, 2), 0.13809802932482335, 1e-12));
  const auto copy = cg::laplacian_smooth(bump, &faces, 0);
  CHECK(copy.data == bump.data);

  std::mt19937 rng(1);
  std::uniform_real_distribution<double> uniform;
  cg::Points flat;
  flat.n = 25;
  flat.d = 2;
  flat.data.resize(flat.n * 2);
  for (double& value : flat.data) {
    value = uniform(rng);
  }
  const auto moved = cg::laplacian_smooth(flat, nullptr, 5, 0.5);
  double shift = 0.0;
  for (std::size_t i = 0; i < flat.data.size(); ++i) {
    shift += std::abs(moved.data[i] - flat.data[i]);
  }
  CHECK(shift > 0.0);
  CHECK(std::isfinite(shift));

  bool threw = false;
  try {
    cg::laplacian_smooth(bump, nullptr);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
  const auto picture = cg::render_smoothing(bump, inverse, faces, 160, 16);
  CHECK(picture.find("width=\"320\"") != std::string::npos);

  const cg::Points center = cg::make_points(2, {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}});
  std::set<std::array<int, 3>> got;
  for (auto tri : cg::delaunay(center).triangles) {
    std::sort(tri.begin(), tri.end());
    got.insert(tri);
  }
  const std::set<std::array<int, 3>> expected_faces{{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {0, 3, 4}};
  CHECK(got == expected_faces);

  cg::Points samples;
  samples.n = 30;
  samples.d = 2;
  samples.data.resize(samples.n * 2);
  std::mt19937 sample_rng(4);
  std::uniform_real_distribution<double> unit;
  for (double& value : samples.data) {
    value = unit(sample_rng);
  }
  const auto mesh = cg::delaunay(samples);
  const auto hull = cg::convex_hull(samples);
  CHECK(mesh.triangles.size() == 2 * samples.n - 2 - hull.indices.size());
  for (const auto& tri : mesh.triangles) {
    const double ax = samples(static_cast<std::size_t>(tri[0]), 0);
    const double ay = samples(static_cast<std::size_t>(tri[0]), 1);
    const double bx = samples(static_cast<std::size_t>(tri[1]), 0);
    const double by = samples(static_cast<std::size_t>(tri[1]), 1);
    const double cx = samples(static_cast<std::size_t>(tri[2]), 0);
    const double cy = samples(static_cast<std::size_t>(tri[2]), 1);
    const double d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    CHECK(std::abs(d) > 1e-12);
    const double aa = ax * ax + ay * ay;
    const double bb = bx * bx + by * by;
    const double cc = cx * cx + cy * cy;
    const double ux = (aa * (by - cy) + bb * (cy - ay) + cc * (ay - by)) / d;
    const double uy = (aa * (cx - bx) + bb * (ax - cx) + cc * (bx - ax)) / d;
    const double radius2 = (ux - ax) * (ux - ax) + (uy - ay) * (uy - ay);
    for (std::size_t i = 0; i < samples.n; ++i) {
      if (static_cast<int>(i) == tri[0] || static_cast<int>(i) == tri[1] || static_cast<int>(i) == tri[2]) {
        continue;
      }
      const double dx = samples(i, 0) - ux;
      const double dy = samples(i, 1) - uy;
      CHECK(dx * dx + dy * dy >= radius2 - 1e-8);
    }
  }
  CHECK(cg::render_delaunay(center, cg::delaunay(center), 200, 20).find("svg") != std::string::npos);
}

}  // namespace

int main() {
  test_closest_pair();
  test_hull();
  test_hull_3d_matches_supporting_planes();
  test_box();
  test_neighbors();
  test_shoelace();
  test_segments();
  test_laplacian_and_delaunay();
  if (g_failed != 0) {
    std::cerr << g_failed << " check(s) failed\n";
    return 1;
  }
  std::cout << "all checks passed\n";
  return 0;
}
