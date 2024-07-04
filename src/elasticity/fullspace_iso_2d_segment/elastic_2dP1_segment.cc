//
// This file part of BigWham
//
// Created by Brice Lecampion on 29.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// contains fundamental plane-strain elasticity kernels.
// for fracture segment with piece-wise linear variation of
// displacement discontinuities (DD)
// convention: Stress positive in tension
//      Displacement discontinuity positive in overlap dd=u^- - u^+
//

#include <il/Array2D.h>
#include <il/linearAlgebra.h>
#include <il/math.h>

#include "core/elastic_properties.h"
#include "elastic_2dP1_segment.h"


namespace bigwham {

il::StaticArray2D<double, 2, 4> stresses_kernel_dp1_dd(double h, double Ep,
                                                       double x, double y) {
  // function computing the stresses at (x,y) induced by a linear DD segment (of
  // total length h) centered on the origin [-h/2,h/2]
  // it returns stresses due to a linear variation from an unit value at the
  // left node (node 1)  to zero at the right node (coordinates 2) for both
  // shear and
  // opening displacement discontinuity
  // and stresses due to a linear variation from an unit value at the right
  // coordinates (node 2) to zero at the right node (coordinates 1)
  // for both shear and opening displacement discontinuity
  //
  // Ep is the plane strain Young'smodulus
  // notation of stress component:
  // 1 : left node influence, 2 : right coordinates influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress
  // e.g.:  sxxs1 - sxx component due to a shear DD with a linear variation with
  // unit value at coordinates 1 .
  // note that we have the following relations : sxxn = sxys,  sxyn = syys such
  // that we don t output all values.

  //  should have a check for h>0
  double overhalfh = 2. / h;
  // minus sign here for convention of positive DD in overlap
  double Elascoef = -1. * Ep / (4. * il::pi);

  double xp = x * overhalfh;  // change of variable to perform the computation
  // for a unit segment [-1,1]
  double yp = y * overhalfh;  // change of variable to perform the computation
  // for a unit segment [-1,1]

  // auxiliary variables.
  double xpm1 = xp - 1.;
  double xpp1 = xp + 1.;
  double dxpp3 = 2. * xp + 3.;
  double xpm2 = xp - 2.;
  double r1 = (xpm1 * xpm1) + (yp * yp);
  double r2 = (xpp1 * xpp1) + (yp * yp);
  double yp2m1 = yp * yp - 1;
  double yp2p1 = yp * yp + 1;

  double AtanAux1 = atanh((2. * xp) / ((xp * xp) + yp2p1)) * 0.5;
  double AtanAux2 = atanh((2. * xp) / ((xp * xp) + yp2p1)) * 0.5;
  //
  double AtanAuxxpm1 = atan(xpm1 / yp);
  double AtanAuxxpp1 = atan(xpp1 / yp);

  double sxxs1, sxxs2;

  // 1 : left node influence, 2 : right coordinates influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress

  double sxys1 = AtanAux1 - (1. / r1) * (1. / (r2 * r2)) *
                                ((xpm1 * xpm1) * (xpp1 * xpp1 * xpp1) +
                                 2. * xp * (3. + xp) * xpp1 * (yp * yp) +
                                 xpm1 * (yp * yp * yp * yp));

  double syys1 =
      (1. / r1) * (1. / (r2 * r2)) *
      (2 * xpm1 * yp * (xpp1 * xpp1) - 2 * (1 + 3 * xp) * (yp * yp * yp));

  double syyn1 =
      AtanAux1 - xpp1 * (1. / (r2 * r2)) * ((xpp1 * xpp1) + 3 * (yp * yp)) +
      2 * xp * (yp * yp) *
          (1. /
           (2 * yp2m1 * (xp * xp) + (xp * xp * xp * xp) + (yp2p1 * yp2p1)));

  double sxys2 = -AtanAux2 + (1. / (r1 * r1)) * (1. / r2) *
                                 ((xpm1 * xpm1 * xpm1) * (xpp1 * xpp1) +
                                  2 * (-3 + xp) * xp * xpm1 * (yp * yp) +
                                  xpp1 * (yp * yp * yp * yp));

  double syys2 =
      (1. / (r1 * r1)) * (1. / r2) *
      (2 * xpp1 * yp * (xpm1 * xpm1) + (2 - 6 * xp) * (yp * yp * yp));

  double syyn2 = -AtanAux2 + (1. / (r1 * r1)) * (1. / r2) *
                                 ((xpm1 * xpm1 * xpm1) * (xpp1 * xpp1) +
                                  2 * (2 + xp) * xpm1 * xpp1 * (yp * yp) +
                                  (-3 + xp) * (yp * yp * yp * yp));

  if (yp != 0.)  // case of observation point on the unit segment line
  {
    sxxs1 = AtanAuxxpm1 - AtanAuxxpp1 +
            2 * yp * (1. / r1) * (1. / (r2 * r2)) *
                (xpm1 * xpm2 * (xpp1 * xpp1) + (3 + dxpp3 * xp) * (yp * yp) +
                 (yp * yp * yp * yp));  //

    sxxs2 = -AtanAuxxpm1 + (1. / (r1 * r1)) * (1. / r2) *
                               (AtanAuxxpp1 * r2 * (r1 * r1) -
                                2 * (2 + xp) * xpp1 * yp * (xpm1 * xpm1) -
                                2 * (3 + xp * (-3 + 2 * xp)) * (yp * yp * yp) -
                                2 * (yp * yp * yp * yp * yp));
  } else {
    sxxs1 =
        ((il::pi) * (xpm1 / sqrt((xpm1 * xpm1)) - xpp1 / sqrt((xpp1 * xpp1)))) /
        2.;  // could be simplified further as 2 or 0 depending on abs(xp)
    sxxs2 = ((il::pi) *
             (-(xpm1 / sqrt((xpm1 * xpm1))) + xpp1 / sqrt((xpp1 * xpp1)))) /
            2.;
  }

  // note that we have the following relations : sxxn = sxys,  sxyn = syys
  // we return a matrix with 4 columns and 2 rows
  // row 1: effect of node 1, row 2 effect of coordinates 2
  // columns sxxs, sxys, syys, syyn    (knowing that we sxxn and sxyn are
  // respectively equal to sxys and syys )

  il::StaticArray2D<double, 2, 4> Stress;

  // switch back to the segment [-h/2,h/2]
  //  we put directly here the constant -E'/(4Pi)
  Stress(0, 0) = Elascoef * overhalfh * sxxs1;
  Stress(0, 1) = Elascoef * overhalfh * sxys1;
  Stress(0, 2) = Elascoef * overhalfh * syys1;
  Stress(0, 3) = Elascoef * overhalfh * syyn1;

  Stress(1, 0) = Elascoef * overhalfh * sxxs2;
  Stress(1, 1) = Elascoef * overhalfh * sxys2;
  Stress(1, 2) = Elascoef * overhalfh * syys2;
  Stress(1, 3) = Elascoef * overhalfh * syyn2;

  return Stress;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//
il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_dp1_dd(
    SegmentData &source_elt, SegmentData &receiver_elt, il::int_t i_col,
    const ElasticProperties &Elas, double ker_options) {
  //   Function to get the normal and shear stress at a point on a surface
  // (with given normal and shear vector) induced by
  // a linear DD segment of size h
  // (the linear DD segment is assumed to be the unit element
  // along the x axis [-h/2,h/2]
  // Material of Plane-strain Young's modulus Ep
  // INPUTS:
  // St : returned stress
  // source_elt : element data structure of the source element
  // receiver_elt  : element data structure of the receiver element
  // i_col : integer for the collocation point number where to compute the
  // normal andshear stress in the  receiver element (0, 1)
  // Elas :: elastic properties object
  // ker_options : dummy argument here (double) - needed for agnostic call ..

  // switch to the frame of the source element....
  il::StaticArray2D<double, 2, 2> R = bigwham::rotationMatrix2D(source_elt.theta());

  il::StaticArray2D<double, 2, 2> Rt = R;
  Rt(0, 1) = R(1, 0);
  Rt(1, 0) = R(0, 1);

  il::StaticArray<double, 2> xe;
  for (int i = 0; i < 2; ++i) {
    xe[i] = receiver_elt.CollocationPoints(i_col, i) - source_elt.Xmid(i);
  }

  xe = il::dot(Rt, xe);

  il::StaticArray<double, 2> n = il::dot(Rt, receiver_elt.n());
  il::StaticArray<double, 2> s = il::dot(Rt, receiver_elt.s());

  double h = source_elt.size();

  double n1n1 = n[0] * n[0];
  double n2n2 = n[1] * n[1];
  double n1n2 = n[0] * n[1];
  double n1s1 = n[0] * s[0];
  double n2s2 = n[1] * s[1];

  double n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

  // columns sxxs, sxys, syys, syyn
  // (knowing that sxxn and sxyn are respectively equal to sxys and syys )
  il::StaticArray2D<double, 2, 4> stress_l =
      stresses_kernel_dp1_dd(h, Elas.young_modulus_plane_strain(), xe[0], xe[1]);

  // shear stress
  // coordinates 1
  // shear dd
  double sh1s = n1s1 * stress_l(0, 0) + n1s2pn2s1 * stress_l(0, 1) +
                n2s2 * stress_l(0, 2);
  // normal dd
  double sh1n = n1s1 * stress_l(0, 1) + n1s2pn2s1 * stress_l(0, 2) +
                n2s2 * stress_l(0, 3);
  // coordinates 2
  // shear dd
  double sh2s = n1s1 * stress_l(1, 0) + n1s2pn2s1 * stress_l(1, 1) +
                n2s2 * stress_l(1, 2);
  // normal dd
  double sh2n = n1s1 * stress_l(1, 1) + n1s2pn2s1 * stress_l(1, 2) +
                n2s2 * stress_l(1, 3);

  // normal stress
  // coordinates 1
  // shear dd
  double sn1s =
      n1n1 * stress_l(0, 0) + 2 * n1n2 * stress_l(0, 1) + n2n2 * stress_l(0, 2);
  // normal dd
  double sn1n =
      n1n1 * stress_l(0, 1) + 2 * n1n2 * stress_l(0, 2) + n2n2 * stress_l(0, 3);
  // coordinates 2
  // shear dd
  double sn2s =
      n1n1 * stress_l(1, 0) + 2 * n1n2 * stress_l(1, 1) + n2n2 * stress_l(1, 2);
  // normal dd
  double sn2n =
      n1n1 * stress_l(1, 1) + 2 * n1n2 * stress_l(1, 2) + n2n2 * stress_l(1, 3);

  // output the desired stress components induced either by normal or shear dd
  // of the 2 nodes of the linear DD segment.
  il::StaticArray2D<double, 2, 4> St;
  // shear stress (node 1 shear dd, normal dd ; coordinates 2 shear dd , normal
  // dd)
  St(0, 0) = sh1s;
  St(0, 1) = sh1n;
  St(0, 2) = sh2s;
  St(0, 3) = sh2n;

  // normal stress (coordinates 1 shear dd, normal dd ; node 2 shear dd , normal
  // dd)
  St(1, 0) = sn1s;
  St(1, 1) = sn1n;
  St(1, 2) = sn2s;
  St(1, 3) = sn2n;
  // St dimensions 2 * 4

  return St;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//      KERNEL FUNCTION WORKING BY NODES - for hmat

il::StaticArray<double, 4> We_segment_1(il::int_t local_node_i,
                                        double h, double Ep,
                                        double x, double y) {
  // function computing the stresses at (x,y) induced by a linear DD segment (of
  // total length h) centered on the origin [-h/2,h/2]
  // it returns stresses due to a linear variation from an unit value at the
  // left node (node 1)  to zero at the right node (coordinates 2) for both
  // shear and opening displacement discontinuity
  // OR stresses due to a linear variation from an unit value at the right
  // coordinates (node 2) to zero at the right node (coordinates 1) for both
  // shear and opening displacement discontinuity
  //
  // Inputs
  // local_node_i :: integer describing if we need to compute the effect of the
  // left (1) or right (2) nodes
  // h :: element size
  // Ep ::  is the plane strain Young'smodulus
  // x :: x position of the point where to compute the stress (element frame)
  // y :: y position of the point where to compute the stress (element frame)

  // Notation of stress component:
  // 1 : left node influence, 2 : right coordinates influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress
  // e.g.:  sxxs1 - sxx component due to a shear DD with a linear variation with
  // unit value at coordinates 1 .
  // note that we have the following relations : sxxn = sxys,  sxyn = syys such
  // that we don t output all values.

  IL_EXPECT_FAST(local_node_i == 0 || local_node_i == 1);

  //  should have a check for h>0
  double overhalfh = 2. / h;
  // minus sign here for convention of positive DD in overlap
  double Elascoef = -1. * Ep / (4. * il::pi);

  double xp = x * overhalfh;  // change of variable to perform the computation
  // for a unit segment [-1,1]
  double yp = y * overhalfh;  // change of variable to perform the computation
  // for a unit segment [-1,1]

  // auxiliary variables.
  double xpm1 = xp - 1.;
  double xpp1 = xp + 1.;
  double dxpp3 = 2. * xp + 3.;
  double xpm2 = xp - 2.;
  double r1 = (xpm1 * xpm1) + (yp * yp);
  double r2 = (xpp1 * xpp1) + (yp * yp);
  double yp2m1 = (yp * yp) - 1;
  double yp2p1 = (yp * yp) + 1;

  double AtanAux1 = atanh((2. * xp) / ((xp * xp) + yp2p1)) * 0.5;
  double AtanAux2 = atanh((2. * xp) / ((xp * xp) + yp2p1)) * 0.5;
  //
  double AtanAuxxpm1 = atan(xpm1 / yp);
  double AtanAuxxpp1 = atan(xpp1 / yp);

  double sxxs1, sxxs2;

  // 1 : left node influence, 2 : right coordinates influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress

  il::StaticArray<double, 4> Stress;

  if (local_node_i == 0) {  // node 1 (in C convention;))

    double sxys1 = AtanAux1 - (1. / r1) * (1. / (r2 * r2)) *
                                  ((xpm1 * xpm1) * (xpp1 * xpp1 * xpp1) +
                                   2. * xp * (3. + xp) * xpp1 * (yp * yp) +
                                   xpm1 * (yp * yp * yp * yp));

    double syys1 =
        (1. / r1) * (1. / (r2 * r2)) *
        (2 * xpm1 * yp * (xpp1 * xpp1) - 2 * (1 + 3 * xp) * (yp * yp * yp));

    double syyn1 =
        AtanAux1 - xpp1 * (1. / (r2 * r2)) * ((xpp1 * xpp1) + 3 * (yp * yp)) +
        2 * xp * yp * yp *
            (1. /
             (2 * yp2m1 * (xp * xp) + (xp * xp * xp * xp) + (yp2p1 * yp2p1)));

    if (yp != 0.)  // case of observation point on the unit sgement line
    {
      sxxs1 = AtanAuxxpm1 - AtanAuxxpp1 +
              2 * yp * (1. / r1) * (1. / (r2 * r2)) *
                  (xpm1 * xpm2 * (xpp1 * xpp1) + (3 + dxpp3 * xp) * (yp * yp) +
                   (yp * yp * yp * yp));  //

    } else {
      sxxs1 = ((il::pi) *
               (xpm1 / sqrt((xpm1 * xpm1)) - xpp1 / sqrt((xpp1 * xpp1)))) /
              2.;  // could be simplified further as 2 or 0 depending on abs(xp)
    }

    // switch back to the segment [-h/2,h/2]
    //  we put directly here the constant -E'/(4Pi)
    Stress[0] = Elascoef * overhalfh * sxxs1;  // Sxx due to shear  dd
    Stress[1] = Elascoef * overhalfh * sxys1;  // Sxx due to normal  dd
    Stress[2] = Elascoef * overhalfh * syys1;  // Syy due to shear  dd
    Stress[3] = Elascoef * overhalfh * syyn1;  // Syy due to normal  dd
  };

  if (local_node_i == 1) {  // node 2 in C convention
    double sxys2 = -AtanAux2 + (1. / (r1 * r1)) * (1. / r2) *
                                   ((xpm1 * xpm1 * xpm1) * (xpp1 * xpp1) +
                                    2 * (-3 + xp) * xp * xpm1 * (yp * yp) +
                                    xpp1 * (yp * yp * yp * yp));

    double syys2 =
        (1. / (r1 * r1)) * (1. / r2) *
        (2 * xpp1 * yp * (xpm1 * xpm1) + (2 - 6 * xp) * (yp * yp * yp));
    double syyn2 = -AtanAux2 + (1. / (r1 * r1)) * (1. / r2) *
                                   ((xpm1 * xpm1 * xpm1) * (xpp1 * xpp1) +
                                    2 * (2 + xp) * xpm1 * xpp1 * (yp * yp) +
                                    (-3 + xp) * (yp * yp * yp * yp));

    if (yp != 0.)  // case of observation point on the unit sgement line
    {
      sxxs2 =
          -AtanAuxxpm1 + (1. / (r1 * r1)) * (1. / r2) *
                             (AtanAuxxpp1 * r2 * (r1 * r1) -
                              2 * (2 + xp) * xpp1 * yp * (xpm1 * xpm1) -
                              2 * (3 + xp * (-3 + 2 * xp)) * (yp * yp * yp) -
                              2 * (yp * yp * yp * yp * yp));
    } else {
      sxxs2 = ((il::pi) *
               (-(xpm1 / sqrt((xpm1 * xpm1))) + xpp1 / sqrt((xpp1 * xpp1)))) /
              2.;
    }

    // switch back to the segment [-h/2,h/2]
    //  we put directly here the constant -E'/(4Pi)

    Stress[0] = Elascoef * overhalfh * sxxs2;
    Stress[1] = Elascoef * overhalfh * sxys2;
    Stress[2] = Elascoef * overhalfh * syys2;
    Stress[3] = Elascoef * overhalfh * syyn2;
  }

  // note that we have the following relations : sxxn = sxys,  sxyn = syys
  // we thus return a vector with 4 entries (instead of 6 or an array of 2*3)
  //   sxxs, sxys, syys, syyn    (knowing that  sxxn and sxyn are respectively equal to sxys and syys )

  return Stress;
}

////////////////////////////////////////////////////////////////////////////////
il::StaticArray2D<double, 2, 2> normal_shear_stress_kernel_dp1_dd_nodal(
    SegmentData &source_elt, SegmentData &receiver_elt, il::int_t s_col,
    il::int_t i_col, const ElasticProperties &Elas, double ker_options) {
  //   Function to get the normal and shear stress at a point on a surface
  // (with given normal and shear vector) induced by
  // a the local node s_col of a linear DD segment of size h
  // (the linear DD segment is assumed to be the unit element
  // along the x axis [-h/2,h/2]
  // Material of Plane-strain Young's modulus Ep
  // INPUTS:
  // St : returned stress
  // source_elt : element data structure of the source element
  // receiver_elt : element data structure of the receiver element
  // s_col : integer describing if we compute here the effect of node 0
  //        or 1 of the source element  (local node number here, i.e. 0 or 1)
  // i_col : integer for the collocation point number of the receiver element
  //        where to compute the normal and shear stress in the receiver element
  //        value either (0, 1)
  // Elas :: elastic properties object  (local node number here, i.e. 0 or 1)
  // ker_options : dummy argument here (double) - needed for agnostic call  due
  // to other kernel function
  // OUTPUT
  // 2*2 matrix containing the effect of node s_col of the source element, on
  // i_col
  // first row shear traction  (effect of shear and opening DD)
  // second row normal traction (effect of shear and opening DD)

  // switch to the frame of the source element....
  il::StaticArray2D<double, 2, 2> R = bigwham::rotationMatrix2D(source_elt.theta());

  il::StaticArray2D<double, 2, 2> Rt = R;
  Rt(0, 1) = R(1, 0);
  Rt(1, 0) = R(0, 1);

  // receiver collocation point in the reference frame of the source element.
  il::StaticArray<double, 2> xe;
  for (int i = 0; i < 2; ++i) {
    xe[i] = receiver_elt.CollocationPoints(i_col, i) - source_elt.Xmid(i);
  }
  xe = il::dot(Rt, xe);

  il::StaticArray<double, 2> n = il::dot(Rt, receiver_elt.n());
  il::StaticArray<double, 2> s = il::dot(Rt, receiver_elt.s());

  double h = source_elt.size();

  double n1n1 = n[0] * n[0];
  double n2n2 = n[1] * n[1];
  double n1n2 = n[0] * n[1];
  double n1s1 = n[0] * s[0];
  double n2s2 = n[1] * s[1];

  double n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

  // columns sxxs, sxys, syys, syyn
  // (knowing that sxxn and sxyn are respectively equal to sxys and syys )
  il::StaticArray<double, 4> stress_l =
          We_segment_1(s_col, h, Elas.young_modulus_plane_strain(), xe[0], xe[1]);

  // shear stress
  // shear dd
  double shs =
      n1s1 * stress_l[0] + n1s2pn2s1 * stress_l[1] + n2s2 * stress_l[2];
  // normal dd
  double shn =
      n1s1 * stress_l[1] + n1s2pn2s1 * stress_l[2] + n2s2 * stress_l[3];

  // normal stress
  // shear dd
  double sns = n1n1 * stress_l[0] + 2 * n1n2 * stress_l[1] + n2n2 * stress_l[2];
  // normal dd
  double snn = n1n1 * stress_l[1] + 2 * n1n2 * stress_l[2] + n2n2 * stress_l[3];

  // output the desired stress components induced either by normal or shear dd
  // of the 2 nodes of the linear DD segment.
  il::StaticArray2D<double, 2, 2> St;

  // shear stress (node 1 shear dd, normal dd ; coordinates 2 shear dd , normal
  // dd)
  St(0, 0) = shs;
  St(0, 1) = shn;

  // normal stress (coordinates 1 shear dd, normal dd ; node 2 shear dd , normal
  // dd)
  St(1, 0) = sns;
  St(1, 1) = snn;

  return St;
}

il::StaticArray<double, 3> point_stress_2d_dp1_dd(
    il::StaticArray<double, 2> &observ_pt, SegmentData &source_elt,
    il::StaticArray<double, 4> &nodal_dd, ElasticProperties &Elas,
    double ker_options) {
  // Function to get the stress at a point
  // induced by a piece-wise linear DD segment of size h
  //
  // Used to calculate stresses at given points
  //
  // INPUTS:
  // source_elt : element data structure of the source element
  // elt_dd : (shear, normal)
  // observ_pt : observation point coordinates (x, y)
  // Elas :: elastic properties object
  //
  // OUTPUTS:
  // stress in global (reference) coordinates (xx, xy, yy)
  // Make sure to add in-situ stress later
  //
  // CONVENTIONS:
  // positive compression / positive opening DD;
  // positive tension / positive overlap DD;

  // the tensor to switch to the frame of the source element....
  il::StaticArray2D<double, 2, 2> R = bigwham::rotationMatrix2D(source_elt.theta());

  // the inverse rotation...
  il::StaticArray2D<double, 2, 2> Rt = R;
  Rt(0, 1) = R(1, 0);
  Rt(1, 0) = R(0, 1);

  // directional vector from the element to the observation point
  il::StaticArray<double, 2> xe;
  for (int i = 0; i < 2; ++i) {
    xe[i] = observ_pt[i] - source_elt.Xmid(i);
  }
  // -"- in source element's coordinates
  xe = il::dot(Rt, xe);

  // stress kernel at the observation point in source element's coordinates
  il::StaticArray2D<double, 2, 4> stress_l =
      stresses_kernel_dp1_dd(source_elt.size(), Elas.young_modulus_plane_strain(), xe[0], xe[1]);

  // note that we have the following relations : sxxn = sxys,  sxyn = syys
  // we return a vector with 4 entries
  //   sxxs, sxys, syys, syyn    (knowing that we sxxn and sxyn are
  // respectively equal to sxys and syys )

  il::StaticArray2D<double, 2, 3> stress_tot_shear, stress_tot_normal;

  for (il::int_t i = 0; i < 2; i++) {
    stress_tot_shear(i, 0) = stress_l(i, 0);
    stress_tot_shear(i, 1) = stress_l(i, 1);
    stress_tot_shear(i, 2) = stress_l(i, 2);

    stress_tot_normal(i, 0) = stress_l(i, 1);
    stress_tot_normal(i, 1) = stress_l(i, 2);
    stress_tot_normal(i, 2) = stress_l(i, 3);
  }

  // stresses at the observation point (in source element's coordinates)
  il::StaticArray<double, 2> elt_shear_dd;
  elt_shear_dd[0] = nodal_dd[0];
  elt_shear_dd[1] = nodal_dd[2];
  il::StaticArray<double, 2> elt_normal_dd;
  elt_normal_dd[0] = nodal_dd[1];
  elt_normal_dd[1] = nodal_dd[3];

  il::StaticArray<double, 3> stress_due_to_shear_g =
      il::dot(elt_shear_dd,stress_tot_shear);
  il::StaticArray<double, 3> stress_due_to_normal_g =
      il::dot(elt_normal_dd,stress_tot_normal);

  // stress tensor in matrix form (in source element's coordinates)
  il::StaticArray2D<double, 2, 2> stress_m, stress_i;
  stress_m(0, 0) = stress_due_to_shear_g[0] + stress_due_to_normal_g[0];
  stress_m(0, 1) = stress_due_to_shear_g[1] + stress_due_to_normal_g[1];
  stress_m(1, 0) = stress_due_to_shear_g[1] + stress_due_to_normal_g[1];
  stress_m(1, 1) = stress_due_to_shear_g[2] + stress_due_to_normal_g[2];

  // stress rotation to global (reference) coordinates
  stress_i = il::dot(stress_m, Rt);
  stress_m = il::dot(R, stress_i);

  il::StaticArray<double, 3> stress_g;
  stress_g[0] = stress_m(0, 0);
  stress_g[1] = 0.5 * (stress_m(0, 1) + stress_m(1, 0));
  stress_g[2] = stress_m(1, 1);

  return stress_g;
};

//------------------------------------------------------------------------------
}  // namespace bigwham


