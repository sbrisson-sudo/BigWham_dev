//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 1/24/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <cmath>
#include <complex>
#include <il/math.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/linearAlgebra.h>
#include <il/linearAlgebra/dense/norm.h>


#include "elasticity/fullspace_iso_3d_triangle/constants.h"
#include "element_utilities.h"

namespace bie {

// auxiliary functions (norm, cross product)

    double l2norm(const il::StaticArray<double, 3> &a) {
// L2 norm of a vector
        double n_a = 0.0;
        for (int k = 0; k < a.size(); ++k) {
            n_a += a[k] * a[k];
        }
        n_a = std::sqrt(n_a);
        return n_a;
    }

    il::StaticArray<double, 3> normalize(const il::StaticArray<double, 3> &a) {
// normalized 3D vector
        il::StaticArray<double, 3> e;
        // double n_a = l2norm(a);
        double n_a = il::norm(a, il::Norm::L2);
        for (int k = 0; k < a.size(); ++k) {
            e[k] = a[k] / n_a;
        }
        return e;
    }

    il::StaticArray<double, 3> cross(const il::StaticArray<double, 3> &a,const il::StaticArray<double, 3> &b) {
// cross product of two 3D vectors
        IL_EXPECT_FAST(a.size() == 3);
        IL_EXPECT_FAST(b.size() == 3);
        il::StaticArray<double, 3> c;
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
        return c;
    }

// Element's local coordinate system manipulations

    il::StaticArray2D<double, 3, 3> make_el_r_tensor(const il::StaticArray2D<double, 3, 3> &el_vert) {
// This function calculates the rotation tensor -
// coordinate transform from "global" (reference) 
// Cartesian coordinate system to the element's local 
// Cartesian coordinate system with origin at the 
// first vertex of the element (el_vert(j, 0))
        il::StaticArray2D<double, 3, 3> r_tensor;
        il::StaticArray<double, 3> a1{}, a2{}, a3{},
                e1{}, e2{}, e3{};
        for (int j = 0; j < 3; ++j) {
            a1[j] = el_vert(j, 1) - el_vert(j, 0);
            a2[j] = el_vert(j, 2) - el_vert(j, 0);
        }
        e1 = normalize(a1); // 1_st tangent <=> the first edge on the triangle
        // a3 = cross(e1, a2);
        a3 = il::cross(e1, a2);
        e3 = normalize(a3); //normal
        e2 = normalize(il::cross(e3, e1)); // 2_nd tangent
        for (int j = 0; j < 3; ++j) {
            r_tensor(0, j) = e1[j];
            r_tensor(1, j) = e2[j];
            r_tensor(2, j) = e3[j];
        }
        return r_tensor;
    }

    il::StaticArray<std::complex<double>, 3> make_el_tau_crd(const il::StaticArray2D<double, 3, 3> &el_vert,const il::StaticArray2D<double, 3, 3> &r_tensor) {
// This function calculates the tau-coordinates
// of the element's vertices
// r-tensor is the rotation tensor for the
// coordinate transform from "global" (reference) 
// Cartesian coordinate system to the element's local one
// Note: use il::dot(r_tensor, a) for R.a
// and il::dot(r_tensor, il::Blas::kTranspose, a) for R^T.a

//  el_vert coordinates in the global system
//  x0 x1 x2
//  y0 y1 y2
//  z0 z1 z2
        il::StaticArray<std::complex<double>, 3> tau{0.0};
        il::StaticArray<double, 3> loc_origin;
        il::StaticArray2D<double, 3, 3> b_vects{0.0};
        for (int k = 0; k < 3; ++k) {
            // Here the 1st vertex is chosen as the origin
            loc_origin[k] = el_vert(k, 0);
        }
        for (int k = 0; k < 3; ++k) {
            for (int j = 0; j < 3; ++j) {
                // Basis vectors (rows)
                b_vects(k, j) = el_vert(k, j) - loc_origin[k];
            }
        }
        // Rotated b_vects
        il::StaticArray2D<double, 3, 3> b_v_r = il::dot(r_tensor, b_vects);
        for (int k = 0; k < 3; ++k) {
            tau[k] = std::complex<double>(b_v_r(0, k), b_v_r(1, k));
        }
        return tau;
    }

    il::StaticArray2D<std::complex<double>, 2, 2> make_el_tau_2_mc(const il::StaticArray2D<double, 3, 3> &el_vert,
             const il::StaticArray2D<double, 3, 3> &r_tensor) {
// This function calculates the coordinate transform
// from local Cartesian coordinates to "master element"
// ([tau, conj(tau)] to [x,y] <[eta1,eta2] according
// to the figure B.1 of Dmitry's thesys>) with origin at the
// first vertex of the element (el_vert(j, 0))
        il::StaticArray2D<std::complex<double>, 2, 2> transform_mtr{0.0};
        il::StaticArray<std::complex<double>, 2> z23{0.0};
        il::StaticArray<double, 3> xsi{0.0}, vv{0.0};
        std::complex<double> m_det;
        //  el_vert coordinates in the global system
        //  x0 x1 x2
        //  y0 y1 y2
        //  z0 z1 z2
        for (int k = 0; k < 2; ++k) {
            for (int n = 0; n < 3; ++n) {
                vv[n] = el_vert(n, k + 1) - el_vert(n, 0);
            }
            xsi = il::dot(r_tensor, vv);
            z23[k] = std::complex<double>(xsi[0], xsi[1]);
        }
        // common denominator (determinant)
        m_det = z23[0] * std::conj(z23[1]) - z23[1] * std::conj(z23[0]);
        // inverse transform
        transform_mtr(0, 0) = std::conj(z23[1]) / m_det;
        transform_mtr(0, 1) = -z23[1] / m_det;
        transform_mtr(1, 0) = -std::conj(z23[0]) / m_det;
        transform_mtr(1, 1) = z23[0] / m_det;
        return transform_mtr;
    }

    HZ make_el_pt_hz(const il::StaticArray2D<double, 3, 3> &el_vert,const il::StaticArray<double, 3> &m_pt_crd,const il::StaticArray2D<double, 3, 3> &r_tensor) {
// This function calculates the h- and tau-coordinates 
// of the point m_pt_cr with respect to the element's 
// local Cartesian coordinate system with origin at 
// the first vertex of the element (el_vert(j, 0))
        HZ h_z;
        il::StaticArray<double, 3> loc_origin, m_pt_c_r;
        for (int k = 0; k < 3; ++k) {
            loc_origin[k] = el_vert(k, 0);
            m_pt_c_r[k] = m_pt_crd[k] - loc_origin[k];
        }
        m_pt_c_r = il::dot(r_tensor, m_pt_c_r);
        h_z.h = -m_pt_c_r[2];
        h_z.z = std::complex<double>(m_pt_c_r[0], m_pt_c_r[1]);
        return h_z;
    }

// Element's basis (shape) functions

    il::StaticArray2D<std::complex<double>, 6, 6> make_el_sfm_uniform(const il::StaticArray2D<double, 3, 3> &el_vert,il::io_t, il::StaticArray2D<double, 3, 3> &r_tensor) {
// This function calculates the basis (shape) functions' 
// coefficients (rows of sfm) for a triangular boundary 
// element with 2nd order polynomial approximation of unknowns
// in terms of complex (tau, conj(tau)) representation of 
// local element's coordinates with trivial (middle) 
// edge partitioning;
// returns the same as 
// make_el_sfm_nonuniform(el_vert, {1.0, 1.0, 1.0})

        // tau_2_mc defines inverse coordinate transform 
        // [tau, conj(tau)] to [x,y] (see make_el_tau_2_mc)
        r_tensor = make_el_r_tensor(el_vert); //creating the rotation tensor from global to local
                                              // in cartesian coord. system (if you multiply on the left of a point)
                                              // I  e_global->e_local      I*v_global=v_local
                                              // I^T*v_local=v_global or v_local^T*I=v_global^T

        il::StaticArray2D<std::complex<double>, 2, 2> tau_2_mc =make_el_tau_2_mc(el_vert, r_tensor);
        il::StaticArray2D<std::complex<double>, 6, 6> sfm{0.0}, sfm_mc{0.0};

        // coefficients of shape functions (rows) 
        // for master element (0<=x,y<=1); ~[1, x, y, x^2, y^2, x*y]
        // Matrix M(s) in D.Nikolskyi notation
        sfm_mc(0, 0) = 1.0;
        sfm_mc(0, 1) = -3.0;
        sfm_mc(0, 2) = -3.0;
        sfm_mc(0, 3) = 2.0;
        sfm_mc(0, 4) = 2.0;
        sfm_mc(0, 5) = 4.0;
        sfm_mc(1, 1) = -1.0;
        sfm_mc(1, 3) = 2.0;
        sfm_mc(2, 2) = -1.0;
        sfm_mc(2, 4) = 2.0;
        sfm_mc(3, 5) = 4.0;
        sfm_mc(4, 2) = 4.0;
        sfm_mc(4, 4) = -4.0;
        sfm_mc(4, 5) = -4.0;
        sfm_mc(5, 1) = 4.0;
        sfm_mc(5, 3) = -4.0;
        sfm_mc(5, 5) = -4.0;

        // inverse coordinate transform 
        // [1, tau, tau_c, tau^2, tau_c^2, tau*tau_c] 
        // to [1, x, y, x^2, y^2, x*y]
        il::StaticArray2D<std::complex<double>, 3, 3> cq{0.0};
        for (int j = 0; j <= 1; ++j) {
            for (int k = 0; k <= 1; ++k) {
                cq(j, k) = tau_2_mc(j, k) * tau_2_mc(j, k);
            }
            cq(j, 2) = 2.0 * tau_2_mc(j, 0) * tau_2_mc(j, 1);
            cq(2, j) = tau_2_mc(0, j) * tau_2_mc(1, j);
        }
        cq(2, 2) = tau_2_mc(0, 0) * tau_2_mc(1, 1) + tau_2_mc(1, 0) * tau_2_mc(0, 1);

        // matrix S in D.Nikolskiy notation
        il::StaticArray2D<std::complex<double>, 6, 6> tau_sq_2_mc{0.0};
        tau_sq_2_mc(0, 0) = 1.0;
        for (int j = 0; j <= 1; ++j) {
            for (int k = 0; k <= 1; ++k) {
                tau_sq_2_mc(j + 1, k + 1) = tau_2_mc(j, k);
            }
        }
        for (int j = 0; j <= 2; ++j) {
            for (int k = 0; k <= 2; ++k) {
                tau_sq_2_mc(j + 3, k + 3) = cq(j, k);
            }
        }
        // assembly of sfm M.S
        sfm = il::dot(sfm_mc, tau_sq_2_mc);
        return sfm;
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    il::StaticArray2D<std::complex<double>, 6, 6> make_el_sfm_nonuniform(const il::StaticArray2D<double, 3, 3> &el_vert,const il::StaticArray<double, 3> &vertex_wts,il::io_t, il::StaticArray2D<double, 3, 3> &r_tensor) {
// This function calculates the basis (shape) functions' 
// coefficients (rows of sfm) for a triangular boundary 
// element with 2nd order polynomial approximation of unknowns
// in terms of complex (tau, conj(tau)) representation of 
// local element's coordinates with non-trivial edge partitioning 
// defined by "weights" vertex_wts

        double p12 = vertex_wts[0] / vertex_wts[1],
                p13 = vertex_wts[0] / vertex_wts[2],
                p23 = vertex_wts[1] / vertex_wts[2],
                c122 = p12 + 1.0,     // (vertex_wts[0]+w(2))/vertex_wts[1];
                c121 = 1.0 / p12 + 1.0,    // (vertex_wts[0]+w(2))/vertex_wts[0];
                c12q = c121 + c122,
                c233 = p23 + 1.0,     // (vertex_wts[1]+w(3))/vertex_wts[2];
                c232 = 1.0 / p23 + 1.0,    // (vertex_wts[1]+w(3))/vertex_wts[1];
                c23q = c232 + c233,
                c133 = p13 + 1.0,     // (vertex_wts[0]+w(3))/vertex_wts[2];
                c131 = 1.0 / p13 + 1.0,    // (vertex_wts[0]+w(3))/vertex_wts[0];
                c13q = c131 + c133;

        r_tensor = make_el_r_tensor(el_vert);
        // tau_2_mc defines inverse coordinate transform 
        // [tau, conj(tau)] to [x,y] (see make_el_tau_2_mc)
        il::StaticArray2D<std::complex<double>, 2, 2> tau_2_mc =make_el_tau_2_mc(el_vert, r_tensor);
        il::StaticArray2D<std::complex<double>, 6, 6> sfm{0.0}, sfm_mc{0.0};

        // coefficients of shape functions (rows) 
        // for master element (0<=x,y<=1); ~[1, x, y, x^2, y^2, x*y]
        sfm_mc(0, 0) = 1.0;
        sfm_mc(0, 1) = -p12 - 2.0;
        sfm_mc(0, 2) = -p13 - 2.0;
        sfm_mc(0, 3) = c122;
        sfm_mc(0, 4) = c133;
        sfm_mc(0, 5) = p13 + p12 + 2.0;
        sfm_mc(1, 1) = -1.0 / p12;
        sfm_mc(1, 3) = c121;
        sfm_mc(1, 5) = 1.0 / p12 - p23;
        sfm_mc(2, 2) = -1.0 / p13;
        sfm_mc(2, 4) = c131;
        sfm_mc(2, 5) = 1.0 / p13 - 1.0 / p23;
        sfm_mc(3, 5) = c23q;
        sfm_mc(4, 2) = c13q;
        sfm_mc(4, 4) = -c13q;
        sfm_mc(4, 5) = -c13q;
        sfm_mc(5, 1) = c12q;
        sfm_mc(5, 3) = -c12q;
        sfm_mc(5, 5) = -c12q;

        // inverse coordinate transform 
        // [1, tau, tau_c, tau^2, tau_c^2, tau*tau_c] 
        // to [1, x, y, x^2, y^2, x*y]
        il::StaticArray2D<std::complex<double>, 3, 3> cq{0.0};
        for (int j = 0; j <= 1; ++j) {
            for (int k = 0; k <= 1; ++k) {
                cq(j, k) = tau_2_mc(j, k) * tau_2_mc(j, k);
            }
            cq(j, 2) = 2.0 * tau_2_mc(j, 0) * tau_2_mc(j, 1);
            cq(2, j) = tau_2_mc(0, j) * tau_2_mc(1, j);
        }
        cq(2, 2) = tau_2_mc(0, 0) * tau_2_mc(1, 1) + 
                tau_2_mc(1, 0) * tau_2_mc(0, 1);

        il::StaticArray2D<std::complex<double>, 6, 6> tau_sq_2_mc{0.0};
        tau_sq_2_mc(0, 0) = 1.0;
        for (int j = 0; j <= 1; ++j) {
            for (int k = 0; k <= 1; ++k) {
                tau_sq_2_mc(j + 1, k + 1) = tau_2_mc(j, k);
            }
        }
        for (int j = 0; j <= 2; ++j) {
            for (int k = 0; k <= 2; ++k) {
                tau_sq_2_mc(j + 3, k + 3) = cq(j, k);
            }
        }
        // assembly of sfm
        sfm = il::dot(sfm_mc, tau_sq_2_mc);
        return sfm;
    }

    //il::StaticArray2D<std::complex<double>, 6, 6> make_el_sfm_beta
    // (const il::StaticArray2D<double, 3, 3> &el_vert,
    // const il::StaticArray<double,3> &vertex_wts
    // double beta,
    // il::io_t, il::StaticArray2D<double, 3, 3> &r_tensor) {
// This function calculates the basis (shape) functions'
// coefficients (rows of sfm) for a triangular boundary
// element with 2nd order polynomial approximation of unknowns
// in terms of complex (tau, conj(tau)) representation of
// local element's coordinates with nodes' offset to the
// centroid (e.g. at collocation points) defined by beta
// and non-trivial edge partitioning defined by
// "weights" vertex_wts
    //};

    il::StaticArray2D<std::complex<double>, 6, 6> shift_el_sfm(std::complex<double> z) {
        // "shifted" sfm from z, tau[m], and local sfm
        // in D.Nikolskiy notation - the Z(s) matrix
        std::complex<double> zc = std::conj(z);
        il::StaticArray2D<std::complex<double>, 6, 6> shift_2_z{0.0};
        shift_2_z(0, 0) = 1.0;
        shift_2_z(1, 0) = z;
        shift_2_z(1, 1) = 1.0;
        shift_2_z(2, 0) = zc;
        shift_2_z(2, 2) = 1.0;
        shift_2_z(3, 0) = z * z;
        shift_2_z(3, 1) = 2.0 * z;
        shift_2_z(3, 3) = 1.0;
        shift_2_z(4, 0) = zc * zc;
        shift_2_z(4, 2) = 2.0 * zc;
        shift_2_z(4, 4) = 1.0;
        shift_2_z(5, 0) = z * zc;
        shift_2_z(5, 1) = zc;
        shift_2_z(5, 2) = z;
        shift_2_z(5, 5) = 1.0;
        return shift_2_z;
    }

// Collocation points

    il::StaticArray<il::StaticArray<double, 3>, 6> el_cp_uniform(const il::StaticArray2D<double, 3, 3> &el_vert, double beta,double betaMiddlenode) {
// This function calculates the coordinates
// of the collocation points on a triangular boundary element
// with 2nd order polynomial approximation of unknowns
// and trivial (middle) edge partitioning;
// offset of the points to the centroid is defined by beta;
// returns the same as el_cp_nonuniform(el_vert, {1.0, 1.0, 1.0}, beta)
        il::StaticArray<il::StaticArray<double, 3>, 6> coll_pt_crd;
        il::StaticArray<double, 3> el_centroid{0.0};

        for (int j = 0; j < 3; ++j) {
            for (int v = 0; v < 3; ++v) {
                el_centroid[j] += el_vert(j, v) / 3.0;
            }
        }
        for (int v = 0; v < 3; ++v) {
            // the edge across the v-th vertex
                                 // v will be 0, 1, 2
            int m = (v + 1) % 3; // m will be 1, 2, 0
            int l = (m + 1) % 3; // l will be 2, 0, 1
            for (int j = 0; j < 3; ++j)
            {  // j is the component of the coord. x y z
                // from the vertex
                (coll_pt_crd[v])[j] =(1.0 - beta) * el_vert(j, v) +beta * el_centroid[j];
                // from the mid side node
                (coll_pt_crd[v + 3])[j] =0.5 * (1.0 - betaMiddlenode) *(el_vert(j, m) + el_vert(j, l)) +betaMiddlenode * el_centroid[j];
            }
        }
        return coll_pt_crd;
        // the ordering of the collocation points in this vector is as follows:
        // [0,1,2,3,4,5]
        //
        //                 0
        //               / + \
        //              /     \
        //             5 +   + 4
        //            /   (+)   \
        //           / +   +   + \
        //          1------3------2
    }

    il::StaticArray<il::StaticArray<double, 3>, 6> el_cp_nonuniform(const il::StaticArray2D<double, 3, 3> &el_vert,const il::StaticArray<double, 3> &vertex_wts,
             double beta,double betaMiddlenode) {
// This function calculates the coordinates
// of the collocation points on a triangular boundary element
// with 2nd order polynomial approximation of unknowns
// and non-trivial (middle) edge partitioning;
// offset of the points to the centroid is defined by beta;
// returns the same as el_cp_nonuniform(el_vert, {1.0, 1.0, 1.0}, beta)
        il::StaticArray<il::StaticArray<double, 3>, 6> coll_pt_crd;
        il::StaticArray<double, 3> el_centroid{0.0};

        for (int j = 0; j < 3; ++j) {
            for (int v = 0; v < 3; ++v) {
                el_centroid[j] += el_vert(j, v) / 3.0;
            }
        }
        for (int v = 0; v < 3; ++v) {
            // the edge across the v-th vertex
            int m = (v + 1) % 3;
            int l = (m + 1) % 3;
            for (int j = 0; j < 3; ++j) {
                // from the vertex
                (coll_pt_crd[v])[j] =
                        (1.0 - beta) * el_vert(j, v) + beta * el_centroid[j];
                // from the mid side node
                (coll_pt_crd[v + 3])[j] =
                        (1.0 - betaMiddlenode) *  (vertex_wts[m] * el_vert(j, m) + vertex_wts[l] * el_vert(j, l)) /
                        (vertex_wts[m] + vertex_wts[l]) + betaMiddlenode * el_centroid[j];
            }
        }
        return coll_pt_crd;
    }

//////////////////////////////////////////////////////////////
    Element_Struct_T set_ele_struct(il::StaticArray2D<double, 3, 3> &el_vert,
             //il::StaticArray<double, 3> &vert_wts,
             double beta,double betaMiddlenode) {
// This function defines the whole set of element properties:
// vertex coordinates, rotational tensor, collocation points,
// coefficients of nodal shape functions, and their values for each CP
        Element_Struct_T ele_s;

        // set vertices' coordinates
        for (il::int_t j = 0; j < 3; ++j) {
            for (il::int_t k = 0; k < 3; ++k) {
                ele_s.vert(k, j) = el_vert(k, j);
            }
            // set vert_wts[j]
        }

        // Basis (shape) functions and rotation tensor of the el-t
        ele_s.sf_m = make_el_sfm_uniform(ele_s.vert, il::io, ele_s.r_tensor);
        //ele_s.sf_m = make_el_sfm_nonuniform
        // (ele_s.r_tensor, ele_s.el_vert, ele_s.vert_wts);

        // Collocation points' coordinates
        ele_s.cp_crd = el_cp_uniform(ele_s.vert, beta, betaMiddlenode);
        //ele_s.cp_crd = el_cp_nonuniform(ele_s.vert, ele_s.vert_wts, beta);

        // values of nodal SF at CP
        for (int n = 0; n < 6; ++n) {
            // position of CP
            HZ hz = make_el_pt_hz(ele_s.vert, ele_s.cp_crd[n], ele_s.r_tensor);
            std::complex<double> tau = hz.z;
            il::StaticArray<std::complex<double>, 6> tau_v {0.0};
            tau_v[0] = 1.0;
            tau_v[1] = tau;
            tau_v[2] = std::conj(tau);
            tau_v[3] = tau_v[1] * tau_v[1];
            tau_v[4] = tau_v[2] * tau_v[2];
            tau_v[5] = tau_v[1] * tau_v[2];
            il::StaticArray<std::complex<double>, 6> sf_cp_c =
                    il::dot(ele_s.sf_m, tau_v);
            for (int k = 0; k < 6; ++k) {
                (ele_s.sf_cp[n])[k] = std::real(sf_cp_c[k]);
            }
        }

        return ele_s;
    }

// Integration of shape functions over the element

    il::StaticArray<std::complex<double>, 6> el_p2_cbp_integral(std::complex<double> a, std::complex<double> b) {
// This function calculates
// the surface-to-contour integral conversion
// (Cauchy-Borel-Pompeiu) for monomials of 2nd order
        std::complex<double> ca = std::conj(a), cb = std::conj(b),
                a2 = a * a, b2 = b * b, ab = a * b,
                ca2 = std::conj(a2), cb2 = std::conj(b2), cab = std::conj(ab);
        std::complex<double> c = 0.25 * ii * (a - b);
        il::StaticArray<std::complex<double>, 6> l_int;

        // const part
        l_int[0] = (ca + cb) * c;

        // linear part; ~\tau , ~\tau\conj
        l_int[1] = (ca * ( 2.0 * a + b ) + cb * ( 2.0 * b + a)) * c / 3.0;
        l_int[2] = (ca2 + cab + cb2) * c / 3.0;

        // quadratic part; ~\tau^2 , ~\tau\conj^2, ~\tau*\tau\conj
        l_int[3] = (ca * (3.0 * a2 + 2.0 * ab + b2) +
                cb * (3.0 * b2 + 2.0 * ab + a2)) * c / 6.0;
        l_int[4] = (ca + cb) * (ca2 + cb2) * c / 6.0;
        l_int[5] = (ca2 * (3.0 * a + b) +
                2.0 * cab * (a + b) +
                cb2 * (3.0 * b + a)) * c / 12.0;
        return l_int;
    }

/*
    il::StaticArray<std::complex<double>, 15> el_p4_cbp_integral
            (std::complex<double> a, std::complex<double> b) {
// This function calculates
// the surface-to-contour integral conversion
// (Cauchy-Borel-Pompeiu) for monomials of 4th order
        std::complex<double> ca = std::conj(a), cb = std::conj(b),
        a2 = a * a, b2 = b * b, ab = a * b, a3 = a * a2, b3 = b * b2,
        a4 = a2 * a2, b4 = b2 * b2, a2b2 = a2 * b2,
        //a5 = a3 * a2, b5 = b3 * b2, a6 = a3 * a3, b6 = b3 * b3,
        amb2 = (a - b) * (a - b),
        ca2 = std::conj(a2), cb2 = std::conj(b2), cab = std::conj(ab),
        ca3 = std::conj(a3), cb3 = std::conj(b3),
        ca4 = std::conj(a4), cb4 = std::conj(b4), ca2b2 = std::conj(a2b2),
        c = 0.25 * ii * (a - b);
        il::StaticArray<std::complex<double>, 15> l_int;

        // const part
        l_int[0] = (ca + cb) * c;

        // linear part; ~\tau , ~\tau\conj
        l_int[1] = (ca * ( 2.0 * a + b ) + cb * ( 2.0 * b + a)) * c / 3.0;
        l_int[2] = (ca2 + cab + cb2) * c / 3.0;

        // quadratic part; ~\tau^2 , ~\tau\conj^2, ~\tau*\tau\conj
        l_int[3] = (ca * (3.0 * a2 + 2.0 * ab + b2) +
                cb * (3.0 * b2 + 2.0 * ab + a2)) * c / 6.0;
        l_int[4] = (ca + cb) * (ca2 + cb2) * c / 6.0;
        l_int[5] = (ca2 * (3.0 * a + b) +
                2.0 * cab * (a + b) +
                cb2 * (3.0 * b + a)) * c / 12.0;

        // cubic part;
        l_int[6] = ((4.0 * a5 - 5.0 * a4 * b + b5) * ca +
                    (a5 - 5.0 * a * b4 + 4.0 * b5) * cb) / amb2 * c / 10.0;
        l_int[7] = ((3.0 * a2 + 4.0 * ab + 3.0 * b2) * cab +
                6.0 * a2 * ca2 + (3.0 * ab + b2) * ca2 +
                (a2 + 3.0 * ab) * cb2 + 6.0 * b2 * cb2) * c / 30.0;
        l_int[8] = ((4.0 * a + b) * ca3 + (3.0 * a + 2.0 * b) * ca2 * cb +
                (2.0 * a + 3.0 * b) * ca * cb2 +
                (a + 4.0 * b) * cb3) * c / 30.0;
        l_int[9] = (ca5 - cb5) / (ca - cb) * c / 10.0;

        // 4-ic part;
        l_int[10] = ((5.0 * a6 - 6.0 * a5 * b + b6) * ca +
                    (a6 - 6.0 * a * b5 + 5.0 * b6) * cb) / amb2 * c / 15.0;
        l_int[11] = ((10.0 * a3 + 6.0 * a2 * b + 3.0 * a * b2 + b3) * ca2 +
                2.0 * (a + b) * (2.0 * a2 + ab + 2.0 b2) * cab +
                (a3 + 3.0 * a2 * b + 6.0 * a * b2 + 10.0 * b3) * cb2) * c / 60.0;
        l_int[12] = ((10.0 * a2 + 4.0 * ab + b2) * ca3 +
                3.0 * (2.0 * a2 + 2.0 * ab + b2) * ca2 * cb +
                3.0 * (a2 + 2.0 * ab + 2.0 * b2) * ca * cb2 +
                (a2 + 4.0 * ab + 10.0 * b2) * c3) * c / 90.0;
        l_int[13] = ((5.0 * a + b) * ca4 + 2.0 * (2.0 * a + b) * ca3 * cb +
                3.0 * (a + b) * ca2b2 + 2.0 (a + 2.0 * b) * ca * cb3 +
                (a + 5.0 * b) * cb4) * c / 60.0;
        l_int[14] = (ca6 - cb6) / (ca - cb) * c / 15.0;
        return l_int;
    }
*/

    il::StaticArray<double, 6> el_p2_sf_integral(il::StaticArray2D<std::complex<double>, 6, 6> el_sfm,il::StaticArray<std::complex<double>, 3> el_tau) {
// This function calculates the integral
// of 2nd order polynomial shape functions
// over an element
        il::StaticArray<std::complex<double>, 6> l_int, m_int{0.0};
        for (int m = 0; m < 3; ++m) {
            int n = (m + 1) % 3;
            std::complex<double> a = el_tau[m], b = el_tau[n];
            l_int = el_p2_cbp_integral(a, b);
            // il::blas(1.0, l_int, 1.0, il::io, m_int);
            for (int k = 0; k < 6; ++k) {
                m_int[k] += l_int[k];
            }
        }
        il::blas(1.0, el_sfm, m_int, 0.0, il::io, l_int);
        il::StaticArray<double, 6> el_sf_int;
        for (int k = 0; k < 6; ++k) {
            el_sf_int[k] = std::real(l_int[k]);
        }
        return el_sf_int;
    }

}
