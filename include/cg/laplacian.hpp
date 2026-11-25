#pragma once

#include "cg/points.hpp"

#include <array>
#include <string>
#include <vector>

namespace cg {

// One iteration: v <- v + lambda * (sum(w v_n) / sum(w) - v).
// weight is "uniform" or "inverse_distance".
// When faces is null, a 2D point set is triangulated with Bowyer-Watson first.
Points laplacian_smooth(const Points& points, const std::vector<std::array<int, 3>>* faces = nullptr,
                        int iterations = 10, double lambda = 0.5,
                        const std::string& weight = "uniform",
                        const std::string& method = "auto");

}  // namespace cg
