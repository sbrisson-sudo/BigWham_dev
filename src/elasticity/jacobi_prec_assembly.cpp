//
// This file is part of HFP.
//
// Created by Brice Lecampion on 12.05.20.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

// Inclusion from Inside Loop library
#include <il/linearAlgebra.h>

// Inclusion from the project
#include <src/elasticity/jacobi_prec_assembly.h>

#include <elasticity/3d/Elastic3DT6_element.h>

namespace bie {

///////////////////////////////////////////////////////////////////////////////
//// 2D kernels with virtual call
il::Array<double> self_influence_elastic(
    Mesh &mesh, const bie::ElasticProperties &elas, vKernelCall KernelCall,
    double ker_options) {
  // Mesh :: 2D mesh object
  // elas :: ElasticProperties object
  // KernelCall :: virtual kernel function

  il::int_t p = mesh.interpolationOrder();

  il::StaticArray2D<double, 2, 2> R;
  il::Array<il::int_t> dofe{2 * (p + 1), 0};

  il::StaticArray2D<double, 2, 4> stnl;

  il::Array<double> diag_Kmat{mesh.numberDDDofs(), 0.};

  for (il::int_t e = 0; e < mesh.numberOfElts(); ++e) {  // loop on all elements

    //   get characteristic of element # e
    bie::SegmentData mysege = mesh.getElementData(e);
    // Rotation matrix of the element w.r. to x-axis.
    R = bie::rotationMatrix2D(mysege.theta());

    // vector of dof id of  element e
    for (il::int_t i = 0; i < 2 * (p + 1); ++i) {
      dofe[i] =e * 2 * (p + 1) + i;// mesh.dofDD(e, i); // todo needs to be changed
    };

    // loop on collocation points of element e
    for (il::int_t ic = 0; ic < p + 1; ++ic) {
      stnl = KernelCall(mysege, mysege, ic, elas, ker_options);

      for (il::int_t j0 = 0; j0 < 2; ++j0) {
        diag_Kmat[dofe[2 * ic] + j0] = stnl(j0, (2*ic) + j0);
      }
    }
  }
  return diag_Kmat;
};



/////////////////////////////////////////////////////////////////////////////////
///// 3D functions -
// note: not generalized - just for T2 , not commented

 il::Array<double> T2_self_influence_elastic(const il::Array<il::int_t> &permutation, bie::Mesh3D &i_meshtools,
    bie::ElasticProperties &elas,
                                              il::int_t I_want_global_DD,
                                              il::int_t I_want_global_traction)
{
il::int_t old_i,e_k1, is_l;
il::Array2D<double> stnl{3, 3, 0.0};
il::Array<double> jacobi_prec{3 * i_meshtools.numberCollPts(),0.0};

for (il::int_t i = 0;i<i_meshtools.numberCollPts();++i) // source nodes
{
// Now:
// from i - permute back to original mesh ordering using permutation of the clusters.
// old_k1 is the (block)-column number in the full (never assembled) elasticity matrix in the global numeration
old_i = permutation[i];

// we import element #0,1,... with connectivity to 3 vertexes but then within each element we have 6 nodes
// for each node 3 DD
// we assume that element # 0 has nodes from 0 to 5
//                        # 1 has nodes from 6 to 11
//                        ... and so on ...
// the range of nodes for one element is: el_number*6+0 to el_number*6+5
// todo: make it general now is only for -> 6 nodes per triangular element (with 3 vertexes)
e_k1 = il::floor(old_i / 6.);  // element number where we have the node old_k1 inside
// floor gives the integer part of the division


// is_l : index from 0 to 5 looping over the nodes of the source element  (we give it to the kernel)
//only for triagular P2 elements
// % gives the remainder of the division
is_l = old_i % 6;

il::Array2D<double> xv = i_meshtools.getVerticesElt(e_k1);
bie::FaceData elem_data_s(xv, 2); // 2 = interpolation order

// call to the kernel
stnl = NodeDDtriplet_to_CPtraction_influence(elem_data_s, elem_data_s, is_l, is_l, elas, I_want_global_DD,I_want_global_traction);

// stnl is a matrix 3x3 like that:
// t_dir_x_node(ir_l)_dd1_on_node_(is_l)  t_dir_x_node(ir_l)_dd2_on_node_(is_l)  t_dir_x_node(ir_l)_dd3_on_node_(is_l)
// t_dir_y_node(ir_l)_dd1_on_node_(is_l)  t_dir_y_node(ir_l)_dd2_on_node_(is_l)  t_dir_y_node(ir_l)_dd3_on_node_(is_l)
// t_dir_z_node(ir_l)_dd1_on_node_(is_l)  t_dir_z_node(ir_l)_dd2_on_node_(is_l)  t_dir_z_node(ir_l)_dd3_on_node_(is_l)

for (il::int_t j = 0;j < 3; j++) {jacobi_prec[i*3 + j] =stnl(j, j); }
}
return jacobi_prec;
}
}
