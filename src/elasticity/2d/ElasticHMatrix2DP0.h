//
// This file is part of HFP.
//
// Created by Brice Lecampion on 23.10.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//
#pragma once

#include <string>

#include <hmat/arrayFunctor/MatrixGenerator.h>

#include <src/core/Mesh2D.h>
#include <src/core/ElasticProperties.h>
#include <elasticity/2d/ElasticS3DP0_element.h>
#include <elasticity/2d/FullMatrixAssembly2D.h>


namespace bie {

template <typename T>
class ElasticHMatrix2DP0 : public il::MatrixGenerator<T> {
 private:
  il::Array2D<double> point_;
  double ker_opts_;

  il::Array<il::int_t> permutation_;
  bie::Mesh mesh_;
  bie::ElasticProperties elas_;

 public:
  ElasticHMatrix2DP0(const il::Array2D<double> &point,
                     const il::Array<il::int_t> &permutation, bie::Mesh &mesh,
                     bie::ElasticProperties &elas, double ker_opts);
  il::int_t size(il::int_t d) const override;
  il::int_t blockSize() const override;
  il::int_t sizeAsBlocks(il::int_t d) const override;
  void set(il::int_t b0, il::int_t b1, il::io_t,
           il::Array2DEdit<T> M) const override;
};

template <typename T>
ElasticHMatrix2DP0<T>::ElasticHMatrix2DP0(
    const il::Array2D<double> &point, const il::Array<il::int_t> &permutation,
    bie::Mesh &mesh, bie::ElasticProperties &elas, double ker_opts)
    : point_{point},
      permutation_{permutation},
      mesh_{mesh},
      elas_{elas},
      ker_opts_{ker_opts} {
  IL_EXPECT_FAST(point_.size(1) == 2);
};

template <typename T>
il::int_t ElasticHMatrix2DP0<T>::size(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);

  return mesh_.numberDDDofs();
};

template <typename T>
il::int_t ElasticHMatrix2DP0<T>::blockSize() const {
  return 2;
}

template <typename T>
il::int_t ElasticHMatrix2DP0<T>::sizeAsBlocks(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);

  return (mesh_.numberDDDofs() / 2);
}

template <typename T>
void ElasticHMatrix2DP0<T>::set(il::int_t b0, il::int_t b1, il::io_t,
                                il::Array2DEdit<T> M) const {
  IL_EXPECT_MEDIUM(M.size(0) % blockSize() == 0);
  IL_EXPECT_MEDIUM(M.size(1) % blockSize() == 0);
  IL_EXPECT_MEDIUM(b0 + M.size(0) / blockSize() <= point_.size(0));
  IL_EXPECT_MEDIUM(b1 + M.size(1) / blockSize() <= point_.size(0));


#ifndef NUMBEROFTHREADS
#define NUMBEROFTHREADS 4
#endif
#pragma omp parallel for num_threads(NUMBEROFTHREADS)
  for (il::int_t j1 = 0; j1 < M.size(1) / blockSize(); ++j1) {
    il::int_t old_k1;
    il::int_t old_k0;
    il::int_t e_k1, e_k0, is_l, ir_l;
    il::StaticArray2D<double, 2, 2> stnl;

    il::int_t k1 = b1 + j1;
    // j1 source node
    // from k1 - permute back to original mesh ordering using permutation of the
    // clusters.
    old_k1 = permutation_[k1];
    e_k1 = il::floor(old_k1 / (mesh_.interpolationOrder() + 1));  // element
    is_l = 1;
    if (old_k1 % (mesh_.interpolationOrder() + 1) == 0) {
      // will not work for quadratic element
      is_l = 0;
    }

    bie::SegmentData seg_s = mesh_.getElementData(e_k1);

    for (il::int_t j0 = 0; j0 < M.size(0) / blockSize(); ++j0) {
      il::int_t k0 = b0 + j0;
      old_k0 = permutation_[k0];
      e_k0 = il::floor(old_k0 / (mesh_.interpolationOrder() + 1));  // element
      ir_l = 1;
      if (old_k0 % (mesh_.interpolationOrder() + 1) ==
          0) {  // will not work for quadratic element
        ir_l = 0;
      }
      bie::SegmentData seg_r = mesh_.getElementData(e_k0);
      il::int_t const p1 = 1;
      stnl = normal_shear_stress_kernel_s3d_dp0_dd_nodal(
          seg_s, seg_r, is_l, ir_l, elas_, ker_opts_);

      for (il::int_t j = 0; j < 2; j++) {
        for (il::int_t i = 0; i < 2; i++) {
          M(j0 * 2 + i, j1 * 2 + j) = stnl(i, j);
        }
      }
    }
  }

}

}
