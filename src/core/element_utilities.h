//
// This file is part of HFP.
//
// Created by D. Nikolski on 1/10/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

// todo - move this file to elasticity/3d ?

#ifndef BIE_ELEM_UTILITIES_H
#define BIE_ELEM_UTILITIES_H

#include <complex>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

namespace bie {

// position of a point with respect to an element
// (complex local coordinate representation)
    struct HZ {
        double h;
        std::complex<double> z;
    };

// element properties in one structure
    struct Element_Struct_T {
        // vertices' coordinates
        il::StaticArray2D<double, 3, 3> vert;
        // vertices' "weights" (defining the positions of edge nodes)
        //il::StaticArray<double, 3> vert_wts;
        // rotation tensor (reference coordinates to el-t local coordinates)
        il::StaticArray2D<double, 3, 3> r_tensor;
        // collocation points' coordinates
        il::StaticArray<il::StaticArray<double, 3>, 6> cp_crd;
        // coefficients of basis (shape) functions of the el-t
        il::StaticArray2D<std::complex<double>, 6, 6> sf_m;
        // values of nodal SF at collocation points
        il::StaticArray<il::StaticArray<double, 6>, 6> sf_cp;
    };

/////// the utilities ///////

// Element's local coordinate system manipulations

    il::StaticArray2D<double, 3, 3> make_el_r_tensor
            (const il::StaticArray2D<double, 3, 3> &el_vert);

    il::StaticArray<std::complex<double>, 3> make_el_tau_crd
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray2D<double, 3, 3> &r_tensor);

    HZ make_el_pt_hz
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray<double, 3> &X0,
             const il::StaticArray2D<double, 3, 3> &r_tensor);

    il::StaticArray2D<std::complex<double>, 2, 2> make_el_tau_2_mc
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray2D<double, 3, 3> &r_tensor);

// Element's basis (shape) functions

    il::StaticArray2D<std::complex<double>, 6, 6> make_el_sfm_uniform
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             il::io_t, il::StaticArray2D<double, 3, 3> &r_tensor);

    il::StaticArray2D<std::complex<double>, 6, 6> make_el_sfm_nonuniform
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray<double, 3> &vertex_wts,
             il::io_t, il::StaticArray2D<double, 3, 3> &r_tensor);

    //il::StaticArray2D<std::complex<double>, 6, 6> make_el_sfm_beta
    // (const il::StaticArray2D<double, 3, 3> &el_vert,
    // const il::StaticArray<double, 3> &vertex_wts,
    // double beta,
    // il::io_t, il::StaticArray2D<double, 3, 3> &r_tensor);

    il::StaticArray2D<std::complex<double>, 6, 6> shift_el_sfm
            (std::complex<double> z);

// Collocation points

    il::StaticArray<il::StaticArray<double, 3>, 6> el_cp_uniform
            (const il::StaticArray2D<double, 3, 3> &el_vert, double beta, double betaMiddlenode);

    il::StaticArray<il::StaticArray<double, 3>, 6> el_cp_nonuniform
            (const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray<double, 3> &vertex_wts,
             double beta, double betaMiddlenode);

    // This function defines the whole set of element properties:
    // vertex coordinates, rotational tensor, collocation points,
    // coefficients of nodal shape functions, and their values for each CP
    Element_Struct_T set_ele_struct(il::StaticArray2D<double, 3, 3> &el_vert,
                        //il::StaticArray<double, 3> %vert_wts,
                        double beta, double betaMiddlepoint);

// Integration over one element

    il::StaticArray<std::complex<double>, 6> el_p2_cbp_integral
            (std::complex<double> a, std::complex<double> b);

/*
    il::StaticArray<std::complex<double>, 15> el_p4_cbp_integral
            (std::complex<double> a, std::complex<double> b);
*/
    il::StaticArray<double, 6> el_p2_sf_integral
            (il::StaticArray2D<std::complex<double>, 6, 6> el_sfm,
             il::StaticArray<std::complex<double>, 3> el_tau);

// auxiliary functions (norm, cross product)

    double l2norm(const il::StaticArray<double, 3> &a);
    il::StaticArray<double, 3> normalize(const il::StaticArray<double, 3> &a);
    il::StaticArray<double, 3> cross(const il::StaticArray<double, 3> &a,
                                     const il::StaticArray<double, 3> &b);

}

#endif //BIE_ELEM_UTILITIES_H
