//
// This file is part of HFP.
//
// Created by Brice Lecampion on  12.01.2023
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//
#pragma once

#include <string>

#include <hmat/arrayFunctor/MatrixGenerator.h>
#include <src/core/BEMesh.h>
#include <src/core/ElasticProperties.h>

namespace bie {

/*

 Class Square Matrix generator for BIE - note that the element type of the source and receiver elements should be the same!

*/
template <typename T,class El,class BIE_Kernel>
class SquareMatrixGenerator : public bie::MatrixGenerator<T> {
 private:

  il::Array<il::int_t> permutation_;
  bie::BEMesh<El> mesh_;

  BIE_Kernel bie_kernel_;

  il::int_t dof_dim_; // unknowns per nodes
  il::int_t size_; // total square matrix of size_*size_
  il::int_t number_points_; // size_ / block_size_

  public:
    SquareMatrixGenerator(bie::BEMesh<El> & mesh, BIE_Kernel & bie_kernel,il::Array<il::int_t>& permutation);
    il::int_t size(il::int_t d) const override;
    il::int_t blockSize() const override;
    il::int_t sizeAsBlocks(il::int_t d) const override;
    void set(il::int_t b0, il::int_t b1, il::io_t,il::Array2DEdit<T> M) const override;

};

template <typename T,class El,class Bie_def>
SquareMatrixGenerator<T,El,Bie_def>::SquareMatrixGenerator(bie::BEMesh<El> & mesh,Bie_def & bie_kernel,il::Array<il::int_t>& permutation) :
      mesh_{mesh},bie_kernel_{bie_kernel},permutation_{permutation}
      {
          number_points_= mesh_.numberCollocationPoints();
          dof_dim_ =bie_kernel.getDofDim() ;
          size_ = number_points_ * dof_dim_ ;
      };



template <typename T,class El,class Bie_def>
il::int_t SquareMatrixGenerator<T,El,Bie_def>::size(il::int_t d) const {
    IL_EXPECT_MEDIUM(d == 0 || d == 1);
  return size_;
};

template <typename T,class El,class Bie_def>
il::int_t SquareMatrixGenerator<T,El,Bie_def>::blockSize() const {
  return dof_dim_;
}

template  <typename T,class El,class Bie_def>
il::int_t SquareMatrixGenerator<T,El,Bie_def>::sizeAsBlocks(il::int_t d) const {
    IL_EXPECT_MEDIUM(d == 0 || d == 1);
  return number_points_;
}

template <typename T,class El,class Bie_def>
void SquareMatrixGenerator<T,El,Bie_def>::set(il::int_t b0, il::int_t b1, il::io_t,
                                il::Array2DEdit<T> M) const {
  IL_EXPECT_MEDIUM(M.size(0) % blockSize() == 0);
  IL_EXPECT_MEDIUM(M.size(1) % blockSize() == 0);
  IL_EXPECT_MEDIUM(b0 + M.size(0) / blockSize() <= number_points_);
  IL_EXPECT_MEDIUM(b1 + M.size(1) / blockSize() <= number_points_);

  il::int_t jj=M.size(1) / blockSize();
#pragma omp parallel if(M.size(1) / blockSize()>200)
  {
#pragma omp for
    for (il::int_t j1 = 0; j1 < M.size(1) / blockSize(); ++j1) {

      El source_elt{},receiver_elt{};
      il::int_t k1 = b1 + j1;
      // j1 source node
      // from k1 - permute back to original mesh ordering using permutation of the clusters.
      il::int_t old_k1 = permutation_[k1];
      il::int_t e_k1 = mesh_.inElement(old_k1);
      il::int_t is_l = mesh_.localCollocationPointId(old_k1);

      il::Array2D<double> xv=mesh_.getVertices(e_k1);
      source_elt.setElement(xv);

      for (il::int_t j0 = 0; j0 < M.size(0) / blockSize(); ++j0) {
        il::int_t k0 = b0 + j0;
        il::int_t old_k0 = permutation_[k0];
        il::int_t e_k0 = mesh_.inElement(old_k0);  //  receiver element
        il::int_t ir_l = mesh_.localCollocationPointId(old_k0);
        receiver_elt.setElement(mesh_.getVertices(e_k0));
        std::vector<double> st=bie_kernel_.influence(source_elt,is_l,receiver_elt,ir_l);
        il::int_t k=0;
        for (il::int_t j = 0; j < dof_dim_; j++) {
          for (il::int_t i = 0; i <  dof_dim_; i++) {
            M(j0 * dof_dim_ + i, j1 * dof_dim_ + j) = st[k];
            k++;
          }
        }
      }
    }
  }

}

}
