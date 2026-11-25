#pragma once

#include <array>
#include <vector>

namespace cg {

using Segment = std::array<std::array<double, 2>, 2>;

struct SegmentIntersection {
  double x = 0.0;
  double y = 0.0;
  int segment_a = 0;
  int segment_b = 0;
};

// Bentley-Ottmann sweep, plus a 1D sweep for overlapping collinear segments.
// Overlaps contribute the endpoints of the overlap. A shared endpoint is reported once.
std::vector<SegmentIntersection> segment_intersections(const std::vector<Segment>& segments);

}  // namespace cg
