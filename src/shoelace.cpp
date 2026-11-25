#include "cg/shoelace.hpp"

#include <cmath>
#include <stdexcept>

namespace cg {

PolygonArea shoelace(const Points& polygon, const std::string& method) {
  require_method(method, {"auto", "shoelace"});
  check_points(polygon, 3, 2);
  std::size_t count = polygon.n;
  if (std::hypot(polygon(0, 0) - polygon(count - 1, 0), polygon(0, 1) - polygon(count - 1, 1)) <=
      1e-12) {
    --count;
  }
  if (count < 3) {
    throw std::invalid_argument("a polygon ring needs at least 3 vertices");
  }
  double sum = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    const std::size_t j = (i + 1) % count;
    sum += polygon(i, 0) * polygon(j, 1) - polygon(i, 1) * polygon(j, 0);
  }
  PolygonArea area;
  area.signed_area = 0.5 * sum;
  area.area = std::abs(area.signed_area);
  return area;
}

}  // namespace cg
