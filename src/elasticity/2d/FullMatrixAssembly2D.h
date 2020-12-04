//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 20.05.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#ifndef BIGWHAM_FULLMATRIXASSEMBLY2D_H
#define BIGWHAM_FULLMATRIXASSEMBLY2D_H

#pragma once

#include <tbb/tbb.h>

#include <src/core/ElasticProperties.h>
#include <src/core/Mesh2D.h>
#include <src/core/SegmentData.h>

#include <src/elasticity/jacobi_prec_assembly.h>

#include <utility>

namespace bie {



typedef il::StaticArray2D<double, 2, 2> (*vKernelCallNodal)(
    SegmentData &source_elt, SegmentData &receiver_elt, il::int_t s_col,  il::int_t i_col,
    const ElasticProperties &Elas, double ker_options);


il::Array2D<double> serialFullMatrix2d(
    Mesh &mesh, const ElasticProperties &elas, vKernelCallNodal KernelCall,
    double ker_options);

il::Array2D<double> parallelFullMatrix2d(Mesh &mesh,
                                         const ElasticProperties &elas,
                                         il::Array<il::int_t> &permutation,
                                         double ker_options,
                                         vKernelCallNodal KernelCall);


//////////////////////////////////////////////////////////////////////////////
// a structure to allow parallel assembly via TBB
// can be used also for the H-mat creation
struct BlockMat {
  il::Array2DEdit<double>  &Mat_;

  const  bie::Mesh mesh_;
  const bie::ElasticProperties elas_;
  const il::Array<il::int_t> permutation_;
  double opts_;
  il::int_t b0_,b1_;

  vKernelCallNodal KernelCall_;

  BlockMat(vKernelCallNodal KernelCall, double ker_opts,
           const bie::ElasticProperties &elas, const bie::Mesh &mesh,
           const il::Array<il::int_t> &permutation, il::int_t b0, il::int_t b1,
           il::Array2DEdit<double> &M)
      :
      Mat_{M},mesh_{mesh},elas_{elas},KernelCall_{KernelCall},opts_{ker_opts},
      permutation_{permutation},b0_{b0},b1_{b1} {};

  void operator()(const tbb::blocked_range<il::int_t>& r) const {

    for( size_t j1=r.begin(); j1<r.end(); ++j1 ){
      il::int_t old_k1;
      il::int_t old_k0;
      il::int_t e_k1, e_k0, is_l, ir_l;
      il::StaticArray2D<double, 2, 2> stnl;

      il::int_t k1 =b1_+ j1;
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

      for (il::int_t j0 = 0; j0 < Mat_.size(0)/2 ; ++j0) {
        il::int_t k0 =b0_+ j0;
        old_k0 = permutation_[k0];
        e_k0 = il::floor(old_k0 / (mesh_.interpolationOrder() + 1));  // element
        ir_l = 1;
        if (old_k0 % (mesh_.interpolationOrder() + 1) ==
            0) {
          // will not work for quadratic element
          ir_l = 0;
        }
        bie::SegmentData seg_r = mesh_.getElementData(e_k0);
        il::int_t const p1 = 1;

        stnl=KernelCall_(seg_s, seg_r, is_l, ir_l, elas_, opts_);

        for (il::int_t j = 0; j < 2; j++) {
          for (il::int_t i = 0; i < 2; i++) {
            Mat_(j0 * 2 + i, j1 * 2 + j) = stnl(i, j);
          }
        }
      }
    }
  }
};

}

#endif  // BIGWHAM_FULLMATRIXASSEMBLY2D_H
