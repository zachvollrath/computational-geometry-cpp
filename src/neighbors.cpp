#include "cg/neighbors.hpp"

#include <algorithm>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cg {
namespace {

struct Item {
  double d2 = 0.0;
  int index = 0;
  bool operator<(const Item& other) const {
    if (d2 != other.d2) {
      return d2 < other.d2;
    }
    return index < other.index;
  }
};

struct KDTree {
  const Points* points = nullptr;
  struct Node {
    int index = 0;
    int axis = 0;
    int left = -1;
    int right = -1;
  };
  std::vector<Node> nodes;

  int build(std::vector<int>& order, int begin, int end, int depth) {
    if (begin >= end) {
      return -1;
    }
    const int axis = depth % static_cast<int>(points->d);
    const int mid = begin + (end - begin) / 2;
    std::nth_element(order.begin() + begin, order.begin() + mid, order.begin() + end,
                     [&](int a, int b) {
                       const double av = (*points)(static_cast<std::size_t>(a), axis);
                       const double bv = (*points)(static_cast<std::size_t>(b), axis);
                       if (av != bv) {
                         return av < bv;
                       }
                       return a < b;
                     });
    const int id = static_cast<int>(nodes.size());
    nodes.push_back(Node{order[static_cast<std::size_t>(mid)], axis, -1, -1});
    const int left = build(order, begin, mid, depth + 1);
    const int right = build(order, mid + 1, end, depth + 1);
    nodes[static_cast<std::size_t>(id)].left = left;
    nodes[static_cast<std::size_t>(id)].right = right;
    return id;
  }

  void search(int node, const double* query, int query_index, std::size_t k,
              std::priority_queue<Item>& best) const {
    if (node < 0) {
      return;
    }
    const Node& current = nodes[static_cast<std::size_t>(node)];
    if (current.index != query_index) {
      double d2 = 0.0;
      for (std::size_t axis = 0; axis < points->d; ++axis) {
        const double delta = query[axis] - (*points)(static_cast<std::size_t>(current.index), axis);
        d2 += delta * delta;
      }
      if (best.size() < k || d2 < best.top().d2 ||
          (d2 == best.top().d2 && current.index < best.top().index)) {
        if (best.size() == k) {
          best.pop();
        }
        best.push(Item{d2, current.index});
      }
    }
    const double delta = query[current.axis] - (*points)(static_cast<std::size_t>(current.index),
                                                         static_cast<std::size_t>(current.axis));
    const int first = delta < 0.0 ? current.left : current.right;
    const int second = delta < 0.0 ? current.right : current.left;
    search(first, query, query_index, k, best);
    const double limit = best.size() < k ? std::numeric_limits<double>::infinity() : best.top().d2;
    if (delta * delta <= limit) {
      search(second, query, query_index, k, best);
    }
  }
};

std::vector<Item> knn_row(const Points& data, const double* query, int query_index, std::size_t k,
                          bool use_tree, const KDTree* tree) {
  std::vector<Item> found;
  found.reserve(k);
  if (use_tree) {
    std::priority_queue<Item> best;
    tree->search(tree->nodes.empty() ? -1 : 0, query, query_index, k, best);
    while (!best.empty()) {
      found.push_back(best.top());
      best.pop();
    }
  } else {
    std::priority_queue<Item> best;
    for (int index = 0; index < static_cast<int>(data.n); ++index) {
      if (index == query_index) {
        continue;
      }
      double d2 = 0.0;
      for (std::size_t axis = 0; axis < data.d; ++axis) {
        const double delta = query[axis] - data(static_cast<std::size_t>(index), axis);
        d2 += delta * delta;
      }
      if (best.size() < k || d2 < best.top().d2 ||
          (d2 == best.top().d2 && index < best.top().index)) {
        if (best.size() == k) {
          best.pop();
        }
        best.push(Item{d2, index});
      }
    }
    while (!best.empty()) {
      found.push_back(best.top());
      best.pop();
    }
  }
  std::sort(found.begin(), found.end(), [](const Item& a, const Item& b) {
    if (a.d2 != b.d2) {
      return a.d2 < b.d2;
    }
    return a.index < b.index;
  });
  return found;
}

}  // namespace

NeighborResult nearest_neighbors(const Points& points, std::size_t k, const Points* queries,
                                 const std::string& method) {
  require_method(method, {"auto", "brute", "kdtree"});
  if (k < 1) {
    throw std::invalid_argument("k must be at least 1");
  }
  check_points(points, 1);
  const bool exclude_self = queries == nullptr;
  Points query_storage;
  const Points* query = queries;
  if (exclude_self) {
    query = &points;
  } else {
    check_points(*queries, 1, static_cast<int>(points.d));
  }
  const std::size_t available = points.n - (exclude_self ? 1 : 0);
  if (k > available) {
    throw std::invalid_argument("k=" + std::to_string(k) + " exceeds the " +
                                std::to_string(available) + " available neighbors");
  }

  const std::string chosen = method == "auto" ? "kdtree" : method;
  KDTree tree;
  if (chosen == "kdtree") {
    tree.points = &points;
    std::vector<int> order(points.n);
    for (std::size_t i = 0; i < points.n; ++i) {
      order[i] = static_cast<int>(i);
    }
    tree.build(order, 0, static_cast<int>(points.n), 0);
  }

  NeighborResult result;
  result.queries = query->n;
  result.k = k;
  result.indices.resize(query->n * k);
  result.distances.resize(query->n * k);
  std::vector<double> row(points.d);
  for (std::size_t r = 0; r < query->n; ++r) {
    for (std::size_t axis = 0; axis < points.d; ++axis) {
      row[axis] = (*query)(r, axis);
    }
    const int skip = exclude_self ? static_cast<int>(r) : -1;
    const std::vector<Item> found =
        knn_row(points, row.data(), skip, k, chosen == "kdtree", &tree);
    if (found.size() < k) {
      throw std::runtime_error("nearest-neighbor query did not return enough neighbors");
    }
    for (std::size_t slot = 0; slot < k; ++slot) {
      result.indices[r * k + slot] = static_cast<std::size_t>(found[slot].index);
      result.distances[r * k + slot] = std::sqrt(std::max(0.0, found[slot].d2));
    }
  }
  return result;
}

}  // namespace cg
