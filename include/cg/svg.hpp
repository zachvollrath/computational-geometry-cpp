#pragma once

#include "cg/bounding_box.hpp"
#include "cg/closest_pair.hpp"
#include "cg/convex_hull.hpp"
#include "cg/delaunay.hpp"
#include "cg/neighbors.hpp"
#include "cg/points.hpp"
#include "cg/segments.hpp"
#include "cg/shoelace.hpp"

#include <array>
#include <string>
#include <vector>

namespace cg {

std::string render_closest_pair(const Points& points, const ClosestPair& pair, int size = 800,
                                int margin = 40);
std::string render_convex_hull(const Points& points, const ConvexHull& hull, int size = 800,
                               int margin = 40);
std::string render_segments(const std::vector<Segment>& segments,
                            const std::vector<SegmentIntersection>& hits, int size = 800,
                            int margin = 40);
std::string render_bounding_box(const Points& points, const OrientedBox& box, int size = 800,
                                int margin = 40);
std::string render_neighbors(const Points& points, const NeighborResult& neighbors, int size = 800,
                             int margin = 40);
std::string render_smoothing(const Points& original, const Points& smoothed,
                             const std::vector<std::array<int, 3>>& faces, int size = 800,
                             int margin = 40);
std::string render_polygon(const Points& polygon, const PolygonArea& area, int size = 800,
                           int margin = 40);
std::string render_delaunay(const Points& points, const Triangulation& mesh, int size = 800,
                            int margin = 40);

void save_svg(const std::string& path, const std::string& svg);

}  // namespace cg
