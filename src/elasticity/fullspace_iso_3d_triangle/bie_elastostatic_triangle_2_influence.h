//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 11.07.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2024.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_BIE_ELASTOSTATIC_TRIANGLE_2_INFLUENCE_H
#define BIGWHAM_BIE_ELASTOSTATIC_TRIANGLE_2_INFLUENCE_H


#pragma once

#include <il/blas.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/StaticArray3D.h>
#include <il/StaticArray4D.h>
#include <il/linearAlgebra.h>


#include "elasticity/bie_elastostatic.h"
#include "elements/triangle.h"
//#include "elements/point.h"
#include "core/oldies/FaceData.h"
#include "elasticity/fullspace_iso_3d_triangle/constants.h"
#include "elasticity/fullspace_iso_3d_triangle/elastic_3dT6_element.h"
#include "elasticity/fullspace_iso_3d_triangle/tensor_utilities_3DT6.h"
#include "elasticity/fullspace_iso_3d_triangle/h_potential_3DT6.h"
#include "elasticity/fullspace_iso_3d_triangle/element_utilities_3DT6.h"

namespace bigwham {

    template
    class BieElastostatic<Triangle<2>, Triangle<2>, ElasticKernelType::H>;
    //   template class BieElastostatic<Triangle<2>, Point<3>, ElasticKernelType::W>;
//    template class BieElastostatic<Triangle<0>, Point<3>, ElasticKernelType::T>;

// Hyper-singular Kernel !
    template<>
    std::vector<double>
    BieElastostatic<Triangle<2>, Triangle<2>, ElasticKernelType::H>::influence(
            const BoundaryElement &source_elt, il::int_t i_s,
            const BoundaryElement &receiver_elt, il::int_t i_r) const {

  //      NodeDDtriplet_to_CPtraction_influence_matrix{3, 3, 0.0};

        auto el_vert_s_aux = source_elt.vertices();
        auto el_vert_r_aux = receiver_elt.vertices();

        bigwham::FaceData elem_data_r(el_vert_r_aux, 2);  // 2 = interpolation order
        bigwham::FaceData elem_data_s(el_vert_s_aux, 2);  // 2 = interpolation order

        il::Array2D<double> NodeDDtriplet_to_CPtraction_influence_matrix=traction_influence_3DT6(elem_data_s, elem_data_r,i_s,i_r, elas_,0,0);
        //         NodeDDtriplet_to_CPtraction_influence_matrix;
        // t_dir_x_node(i_r)_dd1_on_node_(i_s)  t_dir_x_node(i_r)_dd2_on_node_(i_s)  t_dir_x_node(i_r)_dd3_on_node_(i_s)
        // t_dir_y_node(i_r)_dd1_on_node_(i_s)  t_dir_y_node(i_r)_dd2_on_node_(i_s)  t_dir_y_node(i_r)_dd3_on_node_(i_s)
        // t_dir_z_node(i_r)_dd1_on_node_(i_s)  t_dir_z_node(i_r)_dd2_on_node_(i_s)  t_dir_z_node(i_r)_dd3_on_node_(i_s)

        // column major format
        std::vector<double> stnl(9, 0.);
        int k = 0;
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
                // t1/D1,t2/D1,t3/D1, t1/D2,t2/D2,t3/D2 , t1/D3,t2/D3,t3/D3
                stnl[k] = NodeDDtriplet_to_CPtraction_influence_matrix(i, j);
                k++;
            }
        }
        return stnl;

    }
}
#endif //BIGWHAM_BIE_ELASTOSTATIC_TRIANGLE_2_INFLUENCE_H
