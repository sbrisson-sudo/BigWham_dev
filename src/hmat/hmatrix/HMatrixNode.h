//
// This file is part of BigWham.
//
// Created by Francois Fayard - 2018
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#pragma once

#include <hmat/hmatrix/HMatrixType.h>
#include <il/Array2D.h>
#include <il/linearAlgebra/Matrix.h>

namespace il {

template <typename T>
class HMatrixNode {
 private:
  bool empty_;
  il::HMatrixType matrix_type_;
  il::Array2D<T> A_;
  il::Array2D<T> B_;
  il::Array<int> pivot_;

 public:
  HMatrixNode();
  HMatrixNode(il::Array2D<T> A);
  HMatrixNode(il::Array2D<T> A, il::Array2D<T> B);
  il::int_t size(il::int_t d) const;
  bool isEmpty() const;
  bool isFullRank() const;
  bool isLowRank() const;
  bool isHierarchical() const;
  bool isFullLu() const;
  il::HMatrixType type() const;
  il::int_t rankOfLowRank() const;
  void SetEmpty();
  void SetHierarchical();
  void SetFullRank(il::int_t n0, il::int_t n1);
  void SetFullRank(il::Array2D<T> A);
  void SetLowRank(il::int_t n0, il::int_t n1, il::int_t r);
  void SetLowRank(il::Array2D<T> A, il::Array2D<T> B);
  void ConvertToFullLu();
  void UpdateRank(il::int_t r);
  const il::Array2D<T>& asFullRank() const;
  il::Array2D<T>& AsFullRank();
  const il::Array2D<T>& asLowRankA() const;
  il::Array2D<T>& AsLowRankA();
  const il::Array2D<T>& asLowRankB() const;
  il::Array2D<T>& AsLowRankB();
  const il::Array2D<T>& asFullLu() const;
  il::Array2D<T>& AsFullLu();
  const il::Array<int>& asFullLuPivot() const;
  il::Array<int>& AsFullLuPivot();
};

template <typename T>
HMatrixNode<T>::HMatrixNode() : A_{}, B_{}, pivot_{} {
  empty_ = true;
};

template <typename T>
HMatrixNode<T>::HMatrixNode(il::Array2D<T> A)
    : A_{std::move(A)}, B_{}, pivot_{} {
  empty_ = false;
  matrix_type_ = il::HMatrixType::FullRank;
};

template <typename T>
HMatrixNode<T>::HMatrixNode(il::Array2D<T> A, il::Array2D<T> B)
    : A_{std::move(A)}, B_{std::move(B)}, pivot_{} {
  empty_ = false;
  matrix_type_ = il::HMatrixType::LowRank;
}

template <typename T>
il::int_t HMatrixNode<T>::size(il::int_t d) const {
  IL_EXPECT_MEDIUM(!empty_);
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::LowRank ||
                   matrix_type_ == il::HMatrixType::FullRank ||
                   matrix_type_ == il::HMatrixType::FullLu);
  IL_EXPECT_MEDIUM(d == 0 || d == 1);

  if (matrix_type_ == il::HMatrixType::LowRank) {
    if (d == 0) {
      return A_.size(0);
    } else {
      return B_.size(0);
    }
  } else if (matrix_type_ == il::HMatrixType::FullRank ||
             matrix_type_ == il::HMatrixType::FullLu) {
    return A_.size(d);
  }
  IL_UNREACHABLE;
  return -1;
}

template <typename T>
bool HMatrixNode<T>::isEmpty() const {
  return empty_;
}

template <typename T>
bool HMatrixNode<T>::isFullRank() const {
  return !empty_ && matrix_type_ == il::HMatrixType::FullRank;
}

template <typename T>
bool HMatrixNode<T>::isLowRank() const {
  return !empty_ && matrix_type_ == il::HMatrixType::LowRank;
}

template <typename T>
bool HMatrixNode<T>::isHierarchical() const {
  return !empty_ && matrix_type_ == il::HMatrixType::Hierarchical;
}

template <typename T>
bool HMatrixNode<T>::isFullLu() const {
  return !empty_ && matrix_type_ == il::HMatrixType::FullLu;
}

template <typename T>
il::HMatrixType HMatrixNode<T>::type() const {
  IL_EXPECT_MEDIUM(!empty_);

  return matrix_type_;
}

template <typename T>
il::int_t HMatrixNode<T>::rankOfLowRank() const {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::LowRank);
  IL_EXPECT_MEDIUM(A_.size(1) == B_.size(1));

  return A_.size(1);
}

template <typename T>
void HMatrixNode<T>::SetEmpty() {
  empty_ = true;
  A_ = il::Array2D<T>{};
  B_ = il::Array2D<T>{};
  pivot_ = il::Array<int>{};
}

template <typename T>
void HMatrixNode<T>::SetHierarchical() {
  empty_ = false;
  matrix_type_ = il::HMatrixType::Hierarchical;
  A_ = il::Array2D<T>{};
  B_ = il::Array2D<T>{};
  pivot_ = il::Array<int>{};
}

template <typename T>
void HMatrixNode<T>::SetFullRank(il::int_t n0, il::int_t n1) {
  empty_ = false;
  matrix_type_ = il::HMatrixType::FullRank;
  A_.Resize(n0, n1);
}

template <typename T>
void HMatrixNode<T>::SetFullRank(il::Array2D<T> A) {
  empty_ = false;
  matrix_type_ = il::HMatrixType::FullRank;
  A_ = std::move(A);
  B_ = il::Array2D<T>{};
  pivot_ = il::Array<int>{};
}

template <typename T>
void HMatrixNode<T>::SetLowRank(il::int_t n0, il::int_t n1, il::int_t r) {
  empty_ = false;
  matrix_type_ = il::HMatrixType::LowRank;
  A_.Resize(n0, r);
  B_.Resize(n1, r);
}

template <typename T>
void HMatrixNode<T>::SetLowRank(il::Array2D<T> A, il::Array2D<T> B) {
  empty_ = false;
  matrix_type_ = il::HMatrixType::LowRank;
  A_ = std::move(A);
  B_ = std::move(B);
  pivot_ = il::Array<int>{};
}

template <typename T>
void HMatrixNode<T>::ConvertToFullLu() {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::FullRank);
  IL_EXPECT_MEDIUM(A_.size(0) == A_.size(1));

  matrix_type_ = il::HMatrixType::FullLu;
  pivot_.Resize(A_.size(0));
}

template <typename T>
void HMatrixNode<T>::UpdateRank(il::int_t r) {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::LowRank);
  IL_EXPECT_MEDIUM(A_.size(1) == B_.size(1));

  A_.Resize(A_.size(0), r);
  B_.Resize(B_.size(0), r);
}

template <typename T>
const il::Array2D<T>& HMatrixNode<T>::asFullRank() const {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::FullRank);

  return A_;
}

template <typename T>
il::Array2D<T>& HMatrixNode<T>::AsFullRank() {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::FullRank);

  return A_;
}

template <typename T>
const il::Array2D<T>& HMatrixNode<T>::asLowRankA() const {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::LowRank);

  return A_;
}

template <typename T>
il::Array2D<T>& HMatrixNode<T>::AsLowRankA() {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::LowRank);

  return A_;
}

template <typename T>
const il::Array2D<T>& HMatrixNode<T>::asLowRankB() const {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::LowRank);

  return B_;
}

template <typename T>
il::Array2D<T>& HMatrixNode<T>::AsLowRankB() {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::LowRank);

  return B_;
}

template <typename T>
const il::Array2D<T>& HMatrixNode<T>::asFullLu() const {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::FullLu);

  return A_;
}

template <typename T>
il::Array2D<T>& HMatrixNode<T>::AsFullLu() {
  IL_EXPECT_MEDIUM(matrix_type_ == il::HMatrixType::FullLu);

  return A_;
}

template <typename T>
const il::Array<int>& HMatrixNode<T>::asFullLuPivot() const {
  return pivot_;
}

template <typename T>
il::Array<int>& HMatrixNode<T>::AsFullLuPivot() {
  return pivot_;
}

}  // namespace il
