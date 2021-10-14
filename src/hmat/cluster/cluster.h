#ifndef HMATRIX_CLUSTER_H
#define HMATRIX_CLUSTER_H

#include <limits>

#include <il/Array.h>
#include <il/Array2C.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/Tree.h>
#include <il/algorithmArray2D.h>

#include <hmat/hmatrix/HMatrixType.h>

namespace bie {

struct Cluster {
  il::Tree<il::Range, 2> partition;
  il::Array<il::int_t> permutation;
};

Cluster cluster(il::int_t leaf_size, il::io_t, il::Array2D<double>& node);

void cluster_rec(il::spot_t s, il::int_t leaf_size, il::io_t,
                 il::Tree<il::Range, 2>& tree, il::Array2D<double>& node,
                 il::Array<il::int_t>& permutation);

// block-cluster Tree I*I
il::Tree<il::SubHMatrix, 4> hmatrixTreeIxI(const il::Array2D<double>& node,
                                        const il::Tree<il::Range, 2>& tree,
                                        double eta);

void hmatrixTreeIxI_rec(const il::Array2D<double>& node,
                     const il::Tree<il::Range, 2>& range_tree, double eta,
                     il::spot_t s, il::spot_t s0, il::spot_t s1, il::io_t,
                     il::Tree<il::SubHMatrix, 4>& hmatrix_tree);

// block-cluster Tree I*J
il::Tree<il::SubHMatrix, 4> hmatrixTreeIxJ(const il::Array2D<double>& node0,
                                        const il::Tree<il::Range, 2>& tree0,
                                        const il::Array2D<double>& node1,
                                        const il::Tree<il::Range, 2>& tree1,
                                        double eta);

void hmatrixTreeIxJ_rec(const il::Array2D<double>& node0,
                     const il::Tree<il::Range, 2>& range_tree0,
                     const il::Array2D<double>& node1,
                     const il::Tree<il::Range, 2>& range_tree1, double eta,
                     il::spot_t s, il::spot_t s0, il::spot_t s1, il::io_t,
                     il::Tree<il::SubHMatrix, 4>& hmatrix_tree);



////////////////////////////////////////////////////////////////////////////////
//// utilities functions below.
inline double distance(const il::Array2D<double>& node, il::Range range_0,
                       il::Range range_1) {
  IL_EXPECT_FAST(range_0.begin < range_0.end);
  IL_EXPECT_FAST(range_1.begin < range_1.end);
  IL_EXPECT_FAST(
      (range_0.begin == range_1.begin && range_0.end == range_1.end) ||
      range_0.end <= range_1.begin || range_1.end <= range_0.begin);

  const il::int_t dim = node.size(1);

  double distance = 0.0;
  for (il::int_t d = 0; d < dim; ++d) {
    const il::MinMax<double> bound_0 = il::minMax(node, range_0, d);
    const il::MinMax<double> bound_1 = il::minMax(node, range_1, d);
    distance += il::ipow<2>(il::max(0.0, bound_1.min - bound_0.max)) +
                il::ipow<2>(il::max(0.0, bound_0.min - bound_1.max));
  }
  distance = std::sqrt(distance);
  return distance;
}

// distance - for 2 set of points
inline double distance(const il::Array2D<double>& node0,const il::Array2D<double>& node1, il::Range range_0,
                       il::Range range_1) {
  IL_EXPECT_FAST(range_0.begin < range_0.end);
  IL_EXPECT_FAST(range_1.begin < range_1.end);
 // IL_EXPECT_FAST(
   //   (range_0.begin == range_1.begin && range_0.end == range_1.end) ||
    //  range_0.end <= range_1.begin || range_1.end <= range_0.begin);
  IL_EXPECT_FAST(node0.size(1)==node1.size(1));

  const il::int_t dim = node0.size(1);

  double distance = 0.0;
  for (il::int_t d = 0; d < dim; ++d) {
    const il::MinMax<double> bound_0 = il::minMax(node0, range_0, d);
    const il::MinMax<double> bound_1 = il::minMax(node1, range_1, d);
    distance += il::ipow<2>(il::max(0.0, bound_1.min - bound_0.max)) +
                il::ipow<2>(il::max(0.0, bound_0.min - bound_1.max));
  }
  distance = std::sqrt(distance);
  return distance;
}

inline double diameter(const il::Array2D<double>& node, il::Range range) {
  IL_EXPECT_FAST(range.begin < range.end);

  const il::int_t dim = node.size(1);

  double diameter = 0.0;
  for (il::int_t d = 0; d < dim; ++d) {
    const il::MinMax<double> bound = il::minMax(node, range, d);
    diameter += il::ipow<2>(bound.max - bound.min);
  }
  diameter = std::sqrt(diameter);

  return diameter;
}

inline bool isAdmissible(const il::Array2D<double>& node, double eta,
                         il::Range range_0, il::Range range_1) {
  if (range_0.begin == range_1.begin && range_0.end == range_1.end) {
    return false;
  } else {
    const double diam_0 = diameter(node, range_0);
    const double diam_1 = diameter(node, range_1);
    const double dist = distance(node, range_0, range_1);

    return il::max(diam_0, diam_1) <= eta * dist;
  }
}

// case of 2 set of nodes (source / receivers) possibly different
inline bool isAdmissible(const il::Array2D<double>& node0, const il::Array2D<double>& node1,double eta,
                         il::Range range_0, il::Range range_1) {
//  if (range_0.begin == range_1.begin && range_0.end == range_1.end) {
//    return false;
//  } else {
//
    const double diam_0 = diameter(node0, range_0);
    const double diam_1 = diameter(node1, range_1);
    const double dist = distance(node0,node1, range_0, range_1);

    return il::max(diam_0, diam_1) <= eta * dist;
}


// k: is the node number. The tree is numbered as
//
//                 0            -- Level 0
//               /   \
//             1       2        -- Level 1
//           /  |    /   |
//         3     4  5     6     -- Level 2
//                   \
//                    12        -- Level 3 -- Depth 3
//
// At construction we reserve memory as if the tree would be balanced. Then, at
// runtime, we grow the different nodes of the tree.
//
class BinaryTree {
 private:
  il::int_t initial_depth_;
  il::int_t depth_;
  il::Array<il::int_t> data_;

 public:
  BinaryTree(il::int_t nb_nodes, il::int_t leaf_size);
  il::Range operator[](il::int_t) const;
  void addLeft(il::int_t k, il::int_t i_begin, il::int_t i_end);
  void addRight(il::int_t k, il::int_t i_begin, il::int_t i_end);
  bool hasLeftNode(il::int_t k) const;
  bool hasRightNode(il::int_t k) const;
  il::int_t leftNode(il::int_t k) const;
  il::int_t rightNode(il::int_t k) const;
  il::int_t begin(il::int_t k) const;
  il::int_t end(il::int_t k) const;
  il::int_t root() const;
  il::int_t initialDepth() const;
  il::int_t depth() const;
};

inline BinaryTree::BinaryTree(il::int_t nb_nodes, il::int_t leaf_size) {
  IL_EXPECT_FAST(nb_nodes >= 0);
  IL_EXPECT_FAST(leaf_size > 0);

  il::int_t depth = 0;
  il::int_t power = 1;
  while (power * leaf_size < nb_nodes) {
    depth += 1;
    power *= 2;
  }
  const il::int_t n = 2 * (2 * power - 1);

  initial_depth_ = depth;
  depth_ = depth;
  data_.Resize(n);
  for (il::int_t i = 0; i < n; ++i) {
    data_[i] = -1;
  }

  data_[0] = 0;
  data_[1] = nb_nodes;
}

inline il::Range BinaryTree::operator[](il::int_t k) const {
  return il::Range{begin(k), end(k)};
}

inline void BinaryTree::addLeft(il::int_t k, il::int_t i_begin,
                                il::int_t i_end) {
  IL_EXPECT_FAST(2 * k + 1 < data_.size());

  const il::int_t k_left = 2 * k + 1;
  if (2 * (2 * k_left + 1) >= data_.size()) {
    data_.Resize(2 * data_.size() + 1, -1);
    ++depth_;
  }
  data_[2 * k_left] = i_begin;
  data_[2 * k_left + 1] = i_end;
}

inline void BinaryTree::addRight(il::int_t k, il::int_t i_begin,
                                 il::int_t i_end) {
  IL_EXPECT_FAST(2 * k + 1 < data_.size());

  const il::int_t k_right = 2 * k + 2;
  if (2 * (2 * k_right + 1) >= data_.size()) {
    data_.Resize(2 * data_.size() + 1, -1);
    ++depth_;
  }
  data_[2 * k_right] = i_begin;
  data_[2 * k_right + 1] = i_end;
}

inline bool BinaryTree::hasLeftNode(il::int_t k) const {
  return (2 * (2 * k + 1) < data_.size()) && data_[2 * (2 * k + 1)] >= 0;
}

inline bool BinaryTree::hasRightNode(il::int_t k) const {
  return (2 * (2 * k + 2) < data_.size()) && data_[2 * (2 * k + 2)] >= 0;
}

inline il::int_t BinaryTree::leftNode(il::int_t k) const { return 2 * k + 1; }

inline il::int_t BinaryTree::rightNode(il::int_t k) const { return 2 * k + 2; }

inline il::int_t BinaryTree::begin(il::int_t k) const { return data_[2 * k]; }

inline il::int_t BinaryTree::end(il::int_t k) const { return data_[2 * k + 1]; }

inline il::int_t BinaryTree::root() const { return 0; }

inline il::int_t BinaryTree::initialDepth() const { return initial_depth_; }

inline il::int_t BinaryTree::depth() const { return depth_; }

}

#endif  // HMATRIX_CLUSTER_H
