//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 20.05.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <tbb/tbb.h>
#include <elasticity/2d/FullMatrixAssembly2D.h>

namespace  bie{

il::Array2D<double> serialFullMatrix2d(
        Mesh2D &mesh, const ElasticProperties &elas, vKernelCallNodal KernelCall,
        double ker_options) {

  il::Array2D<double> M{mesh.numberCollocationPoints()*2,mesh.numberCollocationPoints()*2,0.};


  for (il::int_t j1 = 0; j1 < M.size(1)/2; j1++) {

    il::int_t old_k1;
    il::int_t old_k0;
    il::int_t e_k1, e_k0, is_l, ir_l;
    il::StaticArray2D<double, 2, 2> stnl;

    il::int_t k1 =  j1;
    // j1 source node
    // from k1 - permute back to original mesh ordering using permutation of the
    // clusters.
    old_k1 = k1;
    e_k1 = il::floor(old_k1 / (mesh.interpolationOrder() + 1));  // element
    is_l = 1;
    if (old_k1 % (mesh.interpolationOrder() + 1) == 0) {
      // will not work for quadratic element
      is_l = 0;
    }

    bie::SegmentData seg_s = mesh.getElementData(e_k1);

    for (il::int_t j0 = 0; j0 < M.size(0)/2 ; ++j0) {
      il::int_t k0 = j0;
      old_k0 = k0;
      e_k0 = il::floor(old_k0 / (mesh.interpolationOrder() + 1));  // element
      ir_l = 1;
      if (old_k0 % (mesh.interpolationOrder() + 1) ==
          0) {
        // will not work for quadratic element
        ir_l = 0;
      }
      bie::SegmentData seg_r = mesh.getElementData(e_k0);
      il::int_t const p1 = 1;

      stnl=KernelCall(seg_s, seg_r, is_l, ir_l, elas, ker_options);

      for (il::int_t j = 0; j < 2; j++) {
        for (il::int_t i = 0; i < 2; i++) {
          M(j0 * 2 + i, j1 * 2 + j) = stnl(i, j);
        }
      }
    }
  }
  return M;
}

/// parallelize assembly

il::Array2D<double> parallelFullMatrix2d(Mesh2D& mesh,
                                         const ElasticProperties& elas,
                                         il::Array<il::int_t>& permutation,
                                         double ker_options,
                                         vKernelCallNodal KernelCall) {

  il::Array2D<double> M{mesh.numberCollocationPoints() * 2,
                        mesh.numberCollocationPoints() * 2, 0.};

  il::Array2DEdit<double> Me=M.Edit();

  tbb::parallel_for(tbb::blocked_range<il::int_t>(0,M.size(1)/2),
                    bie::BlockMat(KernelCall, ker_options, elas, mesh,
                                  permutation, 0, 0, Me));

  return  M;
}


}