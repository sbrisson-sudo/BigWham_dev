
#pragma once

#include <il/StaticArray.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/Map.h>

#include <hmat/hmatrix/HMatrixNode.h>
#include <hmat/hmatrix/HMatrixType.h>

namespace bie {

template <typename T>
class HMatrix {
 private:
  struct HNode {
   public:
    bie::HMatrixNode<T> value;
    il::int_t parent;
    il::StaticArray<il::int_t, 4> children;
   public:
    HNode() : value{}, parent{-1}, children{-1} {};
  };
  il::Array<HNode> tree_;

 public:
  HMatrix();
  HMatrix(il::int_t d);

  il::int_t size(il::int_t d) const;
  il::int_t size(il::int_t d, il::spot_t s) const;

  il::spot_t root() const;
  il::spot_t parent(il::spot_t s) const;
  il::spot_t child(il::spot_t s, il::int_t i0, il::int_t i1) const;

  bool isFullRank(il::spot_t s) const;
  bool isLowRank(il::spot_t s) const;
  bool isHierarchical(il::spot_t s) const;
  bool isFullLu(il::spot_t s) const;
  bie::HMatrixType type(il::spot_t s) const;

  il::int_t rankOfLowRank(il::spot_t s) const;
  void UpdateRank(il::spot_t s, il::int_t r);

  void SetHierarchical(il::spot_t s);
  void SetFullRank(il::spot_t s, il::int_t n0, il::int_t n1);
  void SetLowRank(il::spot_t s, il::int_t n0, il::int_t n1, il::int_t r);
  void ConvertToFullLu(il::spot_t s);

  il::Array2DView<T> asFullRank(il::spot_t s) const;
  il::Array2DEdit<T> AsFullRank(il::spot_t s);

  il::Array2DView<T> asLowRankA(il::spot_t s) const;
  il::Array2DEdit<T> AsLowRankA(il::spot_t s);
  il::Array2DView<T> asLowRankB(il::spot_t s) const;
  il::Array2DEdit<T> AsLowRankB(il::spot_t s);

  il::ArrayView<int> asFullLuPivot(il::spot_t s) const;
  il::ArrayEdit<int> AsFullLuPivot(il::spot_t s);
  il::Array2DView<T> asFullLu(il::spot_t s) const;
  il::Array2DEdit<T> AsFullLu(il::spot_t s);

  bool isBuilt() const;

  void printCompression(il::spot_t s) const;
  il::int_t memorySize(il::spot_t s) const;
  il::int_t memorySize() const;

 private:
  bool isBuilt(il::spot_t s) const;
};

template <typename T>
HMatrix<T>::HMatrix() : tree_{1} {}

template <typename T>
il::int_t HMatrix<T>::size(il::int_t d) const {
  return size(d, il::spot_t{0});
}

template <typename T>
il::spot_t HMatrix<T>::root() const {
  return il::spot_t{0};
}

template <typename T>
il::spot_t HMatrix<T>::parent(il::spot_t s) const {
  return il::spot_t{tree_[s.index].parent};
}

template <typename T>
il::spot_t HMatrix<T>::child(il::spot_t s, il::int_t i0, il::int_t i1) const {
  return il::spot_t{tree_[s.index].children[i1 * 2 + i0]};
}

template <typename T>
bool HMatrix<T>::isFullRank(il::spot_t s) const {
  return tree_[s.index].value.isFullRank();
}

template <typename T>
bool HMatrix<T>::isLowRank(il::spot_t s) const {
  return tree_[s.index].value.isLowRank();
}

template <typename T>
bool HMatrix<T>::isHierarchical(il::spot_t s) const {
  return tree_[s.index].value.isHierarchical();
}

template <typename T>
bool HMatrix<T>::isFullLu(il::spot_t s) const {
  return tree_[s.index].value.isFullLu();
}

template <typename T>
bie::HMatrixType HMatrix<T>::type(il::spot_t s) const {
  return tree_[s.index].value.type();
}

template <typename T>
il::int_t HMatrix<T>::rankOfLowRank(il::spot_t s) const {
  return tree_[s.index].value.rankOfLowRank();
}

template <typename T>
void HMatrix<T>::UpdateRank(il::spot_t s, il::int_t r) {
  return tree_[s.index].value.UpdateRank(r);
}

template <typename T>
void HMatrix<T>::SetHierarchical(il::spot_t s) {
  tree_[s.index].value.SetHierarchical();
  il::int_t n = tree_.size();
  tree_.Resize(n + 4);
  tree_[s.index].children[0] = n;
  tree_[s.index].children[1] = n + 1;
  tree_[s.index].children[2] = n + 2;
  tree_[s.index].children[3] = n + 3;
  tree_[n].parent = s.index;
  tree_[n + 1].parent = s.index;
  tree_[n + 2].parent = s.index;
  tree_[n + 3].parent = s.index;
}

template <typename T>
void HMatrix<T>::SetFullRank(il::spot_t s, il::int_t n0, il::int_t n1) {
  tree_[s.index].value.SetFullRank(n0, n1);
}

template <typename T>
void HMatrix<T>::SetLowRank(il::spot_t s, il::int_t n0, il::int_t n1,
                            il::int_t r) {
  tree_[s.index].value.SetLowRank(n0, n1, r);
}

template <typename T>
void HMatrix<T>::ConvertToFullLu(il::spot_t s) {
  return tree_[s.index].value.ConvertToFullLu();
}

template <typename T>
il::Array2DView<T> HMatrix<T>::asFullRank(il::spot_t s) const {
  return tree_[s.index].value.asFullRank().view();
}

template <typename T>
il::Array2DEdit<T> HMatrix<T>::AsFullRank(il::spot_t s) {
    return tree_[s.index].value.AsFullRank().Edit();
}

template <typename T>
il::Array2DView<T> HMatrix<T>::asLowRankA(il::spot_t s) const {
  return tree_[s.index].value.asLowRankA().view();
}

template <typename T>
il::Array2DEdit<T> HMatrix<T>::AsLowRankA(il::spot_t s) {
  return tree_[s.index].value.AsLowRankA().Edit();
}

template <typename T>
il::Array2DView<T> HMatrix<T>::asLowRankB(il::spot_t s) const {
  return tree_[s.index].value.asLowRankB().view();
}

template <typename T>
il::Array2DEdit<T> HMatrix<T>::AsLowRankB(il::spot_t s) {
  return tree_[s.index].value.AsLowRankB().Edit();
}

template <typename T>
il::ArrayView<int> HMatrix<T>::asFullLuPivot(il::spot_t s) const {
  return tree_[s.index].value.asFullLuPivot().view();
}

template <typename T>
il::ArrayEdit<int> HMatrix<T>::AsFullLuPivot(il::spot_t s) {
  return tree_[s.index].value.AsFullLuPivot().Edit();
}

template <typename T>
il::Array2DView<T> HMatrix<T>::asFullLu(il::spot_t s) const {
  return tree_[s.index].value.asFullLu().view();
}

template <typename T>
il::Array2DEdit<T> HMatrix<T>::AsFullLu(il::spot_t s) {
  return tree_[s.index].value.AsFullLu().Edit();
}

template <typename T>
bool HMatrix<T>::isBuilt() const {
  return isBuilt(root());
}

template <typename T>
void HMatrix<T>::printCompression(il::spot_t s) const {
  tree_[s.index].value.printCompression();
}

template <typename T>
il::int_t HMatrix<T>::memorySize(il::spot_t s) const {
  return tree_[s.index].value.memorySize();
}

template <typename T>
il::int_t HMatrix<T>::memorySize() const {
  il::int_t ans = 0;
  for (il::int_t k = 0; k < tree_.size(); ++k) {
    ans += tree_[k].value.memorySize();
  }
  return ans;
}

template <typename T>
il::int_t HMatrix<T>::size(il::int_t d, il::spot_t s) const {
  if (tree_[s.index].value.isFullRank() || tree_[s.index].value.isLowRank() ||
      tree_[s.index].value.isFullLu()) {
    return tree_[s.index].value.size(d);
  } else if (tree_[s.index].value.isEmpty()) {
    IL_UNREACHABLE;
  } else {
    const il::int_t n00 = size(d, child(s, 0, 0));
    const il::int_t n10 = size(d, child(s, 1, 0));
    const il::int_t n01 = size(d, child(s, 0, 1));
    const il::int_t n11 = size(d, child(s, 1, 1));
    if (d == 0) {
      IL_EXPECT_MEDIUM(n00 == n01);
      IL_EXPECT_MEDIUM(n10 == n11);
      return n00 + n10;
    } else {
      IL_EXPECT_MEDIUM(n00 == n10);
      IL_EXPECT_MEDIUM(n01 == n11);
      return n00 + n01;
    }
  }
  IL_UNREACHABLE;
  return -1;
}

template <typename T>
bool HMatrix<T>::isBuilt(il::spot_t s) const {
  if (tree_[s.index].value.isFullRank() || tree_[s.index].value.isLowRank()) {
    return true;
  } else if (tree_[s.index].value.isEmpty()) {
    return false;
  } else {
    const il::spot_t s00 = child(s, 0, 0);
    const il::spot_t s10 = child(s, 1, 0);
    const il::spot_t s01 = child(s, 0, 1);
    const il::spot_t s11 = child(s, 1, 1);
    return isBuilt(s00) && isBuilt(s10) && isBuilt(s01) && isBuilt(s11);
  }
}

}  // namespace bie