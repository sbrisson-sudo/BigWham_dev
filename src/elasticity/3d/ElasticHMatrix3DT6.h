//
// This file is_n_a part of HFP.
//
// Created by Carlo Peruzzo on 08.08.19.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2019.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#pragma once

#include <Hmat-lib/arrayFunctor/MatrixGenerator.h>
#include <elasticity/3d/Elastic3DT6_element.h>
#include <src/core/ElasticProperties.h>
#include <src/core/TriangularElementData.h>
#include <string>

namespace bie {

template <typename T>
class ElasticHMatrix3DT6 : public il::MatrixGenerator<T> {
 private:
  il::Array2D<double> point_;

  il::Array<il::int_t> permutation_;
  const bie::Mesh3D mesh_;
  //  il::Array2D<il::int_t> binary_ind_pts_at_front_;

 public:
  il::int_t I_want_global_DD;
  il::int_t I_want_global_traction;
  bie::ElasticProperties const elas_;
  ElasticHMatrix3DT6(il::Array2D<double> &point, const il::Array<il::int_t> &permutation,
      bie::Mesh3D &i_meshtools, bie::ElasticProperties &elas,
      il::int_t I_want_global_DD,
      il::int_t I_want_global_traction);

  il::int_t size(il::int_t d) const override;
  il::int_t blockSize() const override;
  il::int_t sizeAsBlocks(il::int_t d) const override;
  void set(il::int_t b0, il::int_t b1, il::io_t,
           il::Array2DEdit<T> M) const override;
};

template <typename T>
ElasticHMatrix3DT6<T>::ElasticHMatrix3DT6(
    il::Array2D<double> &point, const il::Array<il::int_t> &permutation,
    bie::Mesh3D &i_meshtools, bie::ElasticProperties &elas,
    il::int_t I_want_global_DD,
    il::int_t I_want_global_traction)  // il::Array2D<il::int_t>
                                       // &binary_ind_pts_at_front
    : point_{point}, //std::move(point) never fucking do that !
      permutation_{permutation},
      mesh_{i_meshtools},
      elas_{elas},
      I_want_global_DD{I_want_global_DD},
      I_want_global_traction{I_want_global_traction}
//      binary_ind_pts_at_front_{binary_ind_pts_at_front}

{
  IL_EXPECT_FAST(point_.size(1) == 3);  // size(1) point==3 in 3D
};

template <typename T>
il::int_t ElasticHMatrix3DT6<T>::size(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);

  return mesh_
      .numberCollPts() * 3;  // num of nodes * (# of degree of freedom per node)
};

template <typename T>
il::int_t ElasticHMatrix3DT6<T>::blockSize() const {
  return 3;  // # of degree of freedom per node
}

template <typename T>
il::int_t ElasticHMatrix3DT6<T>::sizeAsBlocks(il::int_t d) const {
  IL_EXPECT_MEDIUM(d == 0 || d == 1);

  return (mesh_.numberCollPts());
}

template <typename T>
void ElasticHMatrix3DT6<T>::set(il::int_t b0, il::int_t b1, il::io_t,
                                il::Array2DEdit<T> M) const {
  //    This function, should fill a submatrix of our matrix. The top left
  //    corner of this submatrix should be at (block)-row b0 and (block)-column
  //    b1 and our object should fill the given matrix.
  IL_EXPECT_MEDIUM(M.size(0) % blockSize() == 0);
  IL_EXPECT_MEDIUM(M.size(1) % blockSize() == 0);
  IL_EXPECT_MEDIUM(b0 + M.size(0) / blockSize() <= point_.size(0));
  IL_EXPECT_MEDIUM(b1 + M.size(1) / blockSize() <= point_.size(0));

  //    il::int_t old_k1;
  //    il::int_t old_k0;
  //    il::int_t e_k1, e_k0, is_l, ir_l;
  //    il::Array2D<double> stnl{3,3,0.0};

#ifndef NUMBEROFTHREADS
#define NUMBEROFTHREADS 4
#endif
#pragma omp parallel for num_threads(NUMBEROFTHREADS)
  for (il::int_t j1 = 0; j1 < M.size(1) / blockSize();
       ++j1)  // Loop over a subset of source nodes
  {
    il::int_t old_k1;
    il::int_t old_k0;
    il::int_t e_k1, e_k0, is_l, ir_l;
    il::Array2D<double> stnl{3, 3, 0.0};

    // j1 source node number (it isn't in the global numeration but is local:
    // from 0 to M.size(1) / blockSize()) b1 is the (block)-column number k1 is
    // the index in the full Hmat of the the (block)-column number that we will
    // write
    il::int_t k1 = b1 + j1;

    // Now:
    // from k1 - permute back to original mesh ordering using permutation of the
    // clusters.
    // old_k1 is the (block)-column number in the full (never assembled)
    // elasticity matrix in the global numeration
    old_k1 = permutation_[k1];

    // we import element #0,1,... with connectivity to 3 vertexes but then
    // within each element we have 6 nodes for each node 3 DD we assume that
    // element # 0 has nodes from 0 to 5
    //                        # 1 has nodes from 6 to 11
    //                        ... and so on ...
    // the range of nodes for one element is: el_number*6+0 to el_number*6+5
    // todo: make it general now is only for -> 6 nodes per triangular element
    // (with 3 vertexes)
    e_k1 = il::floor(
        old_k1 / 6.);  // element number where we have the node old_k1 inside
                       // floor gives the integer part of the division

    // is_l : index from 0 to 5 looping over the nodes of the source element (we
    // give it to the kernel)
    // only for triagular P2 elements
    // % gives the remainder of the division
    is_l = old_k1 % 6;


    il::StaticArray2D<double, 3, 3> xv = mesh_.getVerticesElt(e_k1); // get vertices' coordinates of source element
    bie::TriangularElementData elem_data_s(xv, 2); // 2 = interpolation order
    
    // Loop over a subset of collocation points

    for (il::int_t j0 = 0; j0 < M.size(0) / blockSize(); ++j0) {
      il::int_t k0 = b0 + j0;     // analogously as the outer loop
      old_k0 = permutation_[k0];  // analogously as the outer loop
      e_k0 = il::floor(
          old_k0 /
          6.);  // element number where we have the coll point old_k0 inside
                // I got an error here: it was e_k0 = il::floor(old_k1 / 6.);

      // ir_l : index from 0 to 5 looping over the collocation points of the
      // receiver element (we give it to the kernel)
      ir_l = old_k0 % 6;  // only for triagular P2 elements

      xv = mesh_.getVerticesElt(e_k0); // get vertices' coordinates of receiver element
      bie::TriangularElementData elem_data_r(xv, 2); // 2 = interpolation order

      // call to the kernel
      stnl = NodeDDtriplet_to_CPtraction_influence(
          elem_data_s, elem_data_r, is_l, ir_l, elas_, I_want_global_DD,
          I_want_global_traction);

      /* stnl is a matrix 3x3 like that if the source node is NOT at the
         boundary: t_dir_x_node(ir_l)_dd1_on_node_(is_l)
         t_dir_x_node(ir_l)_dd2_on_node_(is_l)
         t_dir_x_node(ir_l)_dd3_on_node_(is_l)
         t_dir_y_node(ir_l)_dd1_on_node_(is_l)
         t_dir_y_node(ir_l)_dd2_on_node_(is_l)
         t_dir_y_node(ir_l)_dd3_on_node_(is_l)
         t_dir_z_node(ir_l)_dd1_on_node_(is_l)
         t_dir_z_node(ir_l)_dd2_on_node_(is_l)
         t_dir_z_node(ir_l)_dd3_on_node_(is_l) */

      /* stnl is a matrix 3x3 like that if the source node is at the boundary:
         t_dir_x_node(ir_l)_dd1_on_node_(is_l)                    0 0 0
         t_dir_y_node(ir_l)_dd2_on_node_(is_l)                    0 0 0
         t_dir_z_node(ir_l)_dd3_on_node_(is_l) */

      for (il::int_t j = 0; j < 3; j++) {
        for (il::int_t i = 0; i < 3; i++) {
          M(j0 * 3 + i, j1 * 3 + j) = stnl(i, j);

          // I'm writing on
          // M( direction , number of DD )
        }
      }
    }
  }
}
}  // namespace bie
