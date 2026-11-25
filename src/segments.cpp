#include "cg/segments.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace cg {
namespace {

constexpr double kEps = 1e-9;

using Vec2 = std::array<double, 2>;

double orient(const Vec2& a, const Vec2& b, const Vec2& c) {
  return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

bool on_segment(const Vec2& a, const Vec2& b, const Vec2& p) {
  return std::min(a[0], b[0]) - 1e-8 <= p[0] && p[0] <= std::max(a[0], b[0]) + 1e-8 &&
         std::min(a[1], b[1]) - 1e-8 <= p[1] && p[1] <= std::max(a[1], b[1]) + 1e-8;
}

int dominant_axis(const Segment& segment) {
  const double dx = std::abs(segment[1][0] - segment[0][0]);
  const double dy = std::abs(segment[1][1] - segment[0][1]);
  return dx >= dy ? 0 : 1;
}

bool collinear(const Segment& a, const Segment& b) {
  double scale = 1.0;
  for (const Segment& segment : {a, b}) {
    for (const Vec2& point : segment) {
      scale = std::max(scale, std::max(std::abs(point[0]), std::abs(point[1])));
    }
  }
  const double tol = kEps * scale;
  return std::abs(orient(a[0], a[1], b[0])) <= tol && std::abs(orient(a[0], a[1], b[1])) <= tol &&
         std::abs(orient(b[0], b[1], a[0])) <= tol && std::abs(orient(b[0], b[1], a[1])) <= tol;
}

Vec2 point_at(const Segment& segment, int axis, double value) {
  const double delta = segment[1][static_cast<std::size_t>(axis)] - segment[0][static_cast<std::size_t>(axis)];
  if (std::abs(delta) <= kEps) {
    return segment[0];
  }
  const double t = (value - segment[0][static_cast<std::size_t>(axis)]) / delta;
  return {segment[0][0] + t * (segment[1][0] - segment[0][0]),
          segment[0][1] + t * (segment[1][1] - segment[0][1])};
}

std::vector<Vec2> unique_points(const std::vector<Vec2>& points) {
  std::vector<Vec2> unique;
  for (const Vec2& point : points) {
    bool seen = false;
    for (const Vec2& kept : unique) {
      if (std::hypot(point[0] - kept[0], point[1] - kept[1]) <= 1e-8) {
        seen = true;
        break;
      }
    }
    if (!seen) {
      unique.push_back(point);
    }
  }
  return unique;
}

std::vector<Vec2> overlap_points(const Segment& a, const Segment& b) {
  const int axis = dominant_axis(a);
  const auto bounds = [&](const Segment& segment) {
    return std::minmax(segment[0][static_cast<std::size_t>(axis)], segment[1][static_cast<std::size_t>(axis)]);
  };
  const auto [a_lo, a_hi] = bounds(a);
  const auto [b_lo, b_hi] = bounds(b);
  const double lo = std::max(a_lo, b_lo);
  const double hi = std::min(a_hi, b_hi);
  if (lo > hi + kEps) {
    return {};
  }
  if (std::abs(hi - lo) <= kEps) {
    return {point_at(a, axis, lo)};
  }
  return {point_at(a, axis, lo), point_at(a, axis, hi)};
}

std::optional<Vec2> proper_point(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d) {
  const double dx = b[0] - a[0];
  const double dy = b[1] - a[1];
  const double ox = d[0] - c[0];
  const double oy = d[1] - c[1];
  const double denom = dx * oy - dy * ox;
  if (std::abs(denom) <= kEps) {
    return std::nullopt;
  }
  const double t = ((c[0] - a[0]) * oy - (c[1] - a[1]) * ox) / denom;
  if (t < -1e-8 || t > 1.0 + 1e-8) {
    return std::nullopt;
  }
  return Vec2{a[0] + t * dx, a[1] + t * dy};
}

std::vector<Vec2> intersect_pair(const Segment& seg_a, const Segment& seg_b) {
  const Vec2& a = seg_a[0];
  const Vec2& b = seg_a[1];
  const Vec2& c = seg_b[0];
  const Vec2& d = seg_b[1];
  const double o1 = orient(a, b, c);
  const double o2 = orient(a, b, d);
  const double o3 = orient(c, d, a);
  const double o4 = orient(c, d, b);
  double scale = 1.0;
  for (const Vec2& point : {a, b, c, d}) {
    scale = std::max(scale, std::max(std::abs(point[0]), std::abs(point[1])));
  }
  const double tol = kEps * scale;
  if (std::abs(o1) <= tol && std::abs(o2) <= tol && std::abs(o3) <= tol && std::abs(o4) <= tol) {
    return overlap_points(seg_a, seg_b);
  }
  std::vector<Vec2> points;
  if (o1 * o2 < 0.0 && o3 * o4 < 0.0) {
    if (const auto point = proper_point(a, b, c, d)) {
      points.push_back(*point);
    }
    return points;
  }
  for (const auto& [corner, side] : {std::pair{c, o1}, std::pair{d, o2}, std::pair{a, o3}, std::pair{b, o4}}) {
    if (std::abs(side) <= tol && on_segment(a, b, corner) && on_segment(c, d, corner)) {
      points.push_back(corner);
    }
  }
  return unique_points(points);
}

std::vector<SegmentIntersection> hits_of(const std::vector<Segment>& segments, int i, int j) {
  if (i > j) {
    std::swap(i, j);
  }
  std::vector<SegmentIntersection> hits;
  for (const Vec2& point : intersect_pair(segments[static_cast<std::size_t>(i)],
                                          segments[static_cast<std::size_t>(j)])) {
    hits.push_back(SegmentIntersection{point[0], point[1], i, j});
  }
  return hits;
}

std::vector<std::vector<int>> collinear_families(const std::vector<Segment>& segments) {
  std::vector<std::vector<int>> families;
  for (int index = 0; index < static_cast<int>(segments.size()); ++index) {
    bool placed = false;
    for (auto& family : families) {
      if (collinear(segments[static_cast<std::size_t>(family.front())],
                    segments[static_cast<std::size_t>(index)])) {
        family.push_back(index);
        placed = true;
        break;
      }
    }
    if (!placed) {
      families.push_back({index});
    }
  }
  return families;
}

std::vector<SegmentIntersection> interval_sweep(const std::vector<Segment>& segments,
                                                const std::vector<int>& family) {
  const int axis = dominant_axis(segments[static_cast<std::size_t>(family.front())]);
  const int other = 1 - axis;
  struct Event {
    double x;
    double y;
    int kind;
    int index;
    bool operator<(const Event& rhs) const {
      return std::tie(x, y, kind, index) < std::tie(rhs.x, rhs.y, rhs.kind, rhs.index);
    }
  };
  std::vector<Event> events;
  for (int index : family) {
    const Segment& segment = segments[static_cast<std::size_t>(index)];
    const bool forward = segment[0][static_cast<std::size_t>(axis)] < segment[1][static_cast<std::size_t>(axis)] ||
                         (segment[0][static_cast<std::size_t>(axis)] == segment[1][static_cast<std::size_t>(axis)] &&
                          segment[0][static_cast<std::size_t>(other)] <= segment[1][static_cast<std::size_t>(other)]);
    const Vec2& lo = forward ? segment[0] : segment[1];
    const Vec2& hi = forward ? segment[1] : segment[0];
    events.push_back(Event{lo[static_cast<std::size_t>(axis)], lo[static_cast<std::size_t>(other)], 0, index});
    events.push_back(Event{hi[static_cast<std::size_t>(axis)], hi[static_cast<std::size_t>(other)], 1, index});
  }
  std::sort(events.begin(), events.end());
  std::vector<int> active;
  std::vector<SegmentIntersection> found;
  for (const Event& event : events) {
    if (event.kind == 0) {
      for (int other_index : active) {
        auto hits = hits_of(segments, other_index, event.index);
        found.insert(found.end(), hits.begin(), hits.end());
      }
      active.push_back(event.index);
    } else {
      active.erase(std::find(active.begin(), active.end(), event.index));
    }
  }
  return found;
}

std::vector<SegmentIntersection> bentley_ottmann(const std::vector<Segment>& segments,
                                                 const std::set<std::pair<int, int>>& skip) {
  std::vector<Segment> ordered = segments;
  for (Segment& segment : ordered) {
    if (std::tie(segment[1][0], segment[1][1]) < std::tie(segment[0][0], segment[0][1])) {
      std::swap(segment[0], segment[1]);
    }
  }

  struct Payload {
    bool is_pair = false;
    int a = 0;
    int b = 0;
  };
  std::vector<Payload> payloads;
  using Event = std::tuple<double, double, int, int>;
  std::vector<Event> heap;
  auto push = [&](double x, double y, int kind, Payload data) {
    const int token = static_cast<int>(payloads.size());
    payloads.push_back(data);
    heap.emplace_back(x, y, kind, token);
    std::push_heap(heap.begin(), heap.end(), std::greater<Event>{});
  };
  for (int index = 0; index < static_cast<int>(ordered.size()); ++index) {
    push(ordered[static_cast<std::size_t>(index)][0][0], ordered[static_cast<std::size_t>(index)][0][1], 0,
         Payload{false, index, 0});
    push(ordered[static_cast<std::size_t>(index)][1][0], ordered[static_cast<std::size_t>(index)][1][1], 2,
         Payload{false, index, 0});
  }

  std::vector<int> status;
  std::vector<SegmentIntersection> found;
  std::set<std::pair<int, int>> seen;
  double sweep_x = -std::numeric_limits<double>::infinity();

  auto y_at = [&](int index, double x) {
    const Segment& segment = ordered[static_cast<std::size_t>(index)];
    const double dx = segment[1][0] - segment[0][0];
    if (std::abs(dx) <= kEps) {
      return segment[0][1];
    }
    const double t = std::min(1.0, std::max(0.0, (x - segment[0][0]) / dx));
    return segment[0][1] + t * (segment[1][1] - segment[0][1]);
  };
  auto slope = [&](int index) {
    const Segment& segment = ordered[static_cast<std::size_t>(index)];
    const double dx = segment[1][0] - segment[0][0];
    const double dy = segment[1][1] - segment[0][1];
    if (std::abs(dx) <= kEps) {
      return dy >= 0.0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
    }
    return dy / dx;
  };
  auto reorder = [&]() {
    std::sort(status.begin(), status.end(), [&](int a, int b) {
      const double ya = y_at(a, sweep_x);
      const double yb = y_at(b, sweep_x);
      if (ya != yb) {
        return ya < yb;
      }
      const double sa = slope(a);
      const double sb = slope(b);
      if (sa != sb) {
        return sa < sb;
      }
      return a < b;
    });
  };
  auto locate = [&](int index) {
    const auto it = std::find(status.begin(), status.end(), index);
    return it == status.end() ? -1 : static_cast<int>(it - status.begin());
  };

  auto consider = [&](int i, int j) {
    if (i > j) {
      std::swap(i, j);
    }
    const std::pair<int, int> pair{i, j};
    if (skip.count(pair) || seen.count(pair)) {
      return;
    }
    auto hits = hits_of(segments, i, j);
    if (hits.empty()) {
      return;
    }
    std::vector<SegmentIntersection> ahead;
    for (const auto& hit : hits) {
      if (hit.x >= sweep_x - 1e-7) {
        ahead.push_back(hit);
      }
    }
    if (ahead.empty()) {
      return;
    }
    seen.insert(pair);
    found.insert(found.end(), hits.begin(), hits.end());
    const auto point = std::min_element(ahead.begin(), ahead.end(), [](const auto& a, const auto& b) {
      return std::tie(a.x, a.y) < std::tie(b.x, b.y);
    });
    push(point->x, point->y, 1, Payload{true, i, j});
  };

  while (!heap.empty()) {
    std::pop_heap(heap.begin(), heap.end(), std::greater<Event>{});
    const auto [x, y, kind, token] = heap.back();
    heap.pop_back();
    (void)y;
    sweep_x = x;
    const Payload data = payloads[static_cast<std::size_t>(token)];
    if (kind == 0) {
      status.push_back(data.a);
      reorder();
      const int pos = locate(data.a);
      if (pos > 0) {
        consider(status[static_cast<std::size_t>(pos - 1)], status[static_cast<std::size_t>(pos)]);
      }
      if (pos + 1 < static_cast<int>(status.size())) {
        consider(status[static_cast<std::size_t>(pos)], status[static_cast<std::size_t>(pos + 1)]);
      }
    } else if (kind == 2) {
      const int pos = locate(data.a);
      if (pos < 0) {
        continue;
      }
      const int below = pos > 0 ? status[static_cast<std::size_t>(pos - 1)] : -1;
      const int above = pos + 1 < static_cast<int>(status.size()) ? status[static_cast<std::size_t>(pos + 1)] : -1;
      status.erase(status.begin() + pos);
      if (below >= 0 && above >= 0) {
        consider(below, above);
      }
    } else {
      int pi = locate(data.a);
      int pj = locate(data.b);
      if (pi < 0 || pj < 0 || std::abs(pi - pj) != 1) {
        continue;
      }
      if (pi > pj) {
        std::swap(pi, pj);
      }
      std::swap(status[static_cast<std::size_t>(pi)], status[static_cast<std::size_t>(pj)]);
      if (pi > 0) {
        consider(status[static_cast<std::size_t>(pi - 1)], status[static_cast<std::size_t>(pi)]);
      }
      if (pj + 1 < static_cast<int>(status.size())) {
        consider(status[static_cast<std::size_t>(pj)], status[static_cast<std::size_t>(pj + 1)]);
      }
    }
  }
  return found;
}

std::vector<SegmentIntersection> unique_hits(const std::vector<SegmentIntersection>& found) {
  std::vector<SegmentIntersection> unique;
  for (const auto& hit : found) {
    bool seen = false;
    for (const auto& kept : unique) {
      if (kept.segment_a == hit.segment_a && kept.segment_b == hit.segment_b &&
          std::hypot(kept.x - hit.x, kept.y - hit.y) <= 1e-7) {
        seen = true;
        break;
      }
    }
    if (!seen) {
      unique.push_back(hit);
    }
  }
  return unique;
}

}  // namespace

std::vector<SegmentIntersection> segment_intersections(const std::vector<Segment>& segments) {
  if (segments.size() < 2) {
    throw std::invalid_argument("intersection requires at least 2 segments");
  }
  for (const Segment& segment : segments) {
    for (const auto& point : segment) {
      if (!std::isfinite(point[0]) || !std::isfinite(point[1])) {
        throw std::invalid_argument("segments must be finite");
      }
    }
  }
  const auto families = collinear_families(segments);
  std::set<std::pair<int, int>> collinear_pairs;
  std::vector<SegmentIntersection> found;
  for (const auto& family : families) {
    if (family.size() < 2) {
      continue;
    }
    for (std::size_t i = 0; i < family.size(); ++i) {
      for (std::size_t j = i + 1; j < family.size(); ++j) {
        const int a = family[i];
        const int b = family[j];
        collinear_pairs.emplace(std::min(a, b), std::max(a, b));
      }
    }
    auto hits = interval_sweep(segments, family);
    found.insert(found.end(), hits.begin(), hits.end());
  }
  auto swept = bentley_ottmann(segments, collinear_pairs);
  found.insert(found.end(), swept.begin(), swept.end());
  found = unique_hits(found);
  std::sort(found.begin(), found.end(), [](const SegmentIntersection& a, const SegmentIntersection& b) {
    return std::tie(a.x, a.y, a.segment_a, a.segment_b) < std::tie(b.x, b.y, b.segment_a, b.segment_b);
  });
  return found;
}

}  // namespace cg
