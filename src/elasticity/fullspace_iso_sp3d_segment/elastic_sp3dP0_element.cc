//
// This file part of BigWham
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include <limits>

#include <il/linearAlgebra.h>
#include <il/math.h>

#include "core/elastic_properties.h"
#include "elastic_sp3dP0_element.h"

namespace bigwham {
//------------------------------------------------------------------------------
// implementation of the Wu & Olson (2015) Simplified 3D kernels for constant
// height
// fractures

// as per their notation, the height is along direction 2
// e_1 (x) corresponds to e_x in 2D plane-strain
// e_3 (z) corresponds to e_y in 2D plane-strain

// RONGVED SOLUTION FOR A P0 Rectangular dislocation in a full space
// dislocation is centered on the origin in the plane e_3=0 (z) , (-a,a) in e_1
// (x) (-b,b)
// in e_2 (z

// AUXILIARY FUNCTIONS NEEDED
// second order derivatives of I(x,y,z,xi,eta)
double ip11(double x, double y, double z, double xi, double eta) {
  //  (x - \[Xi])/((y - \[Eta] + Sqrt(Power(z,2) + Power(y - \[Eta],2) + Power(x
  //  - \[Xi],2)))*
  //      Sqrt(Power(z,2) + Power(y - \[Eta],2) + Power(x - \[Xi],2)))
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (x - xi) / ((R + y - eta) * R);
}

double ip12(double x, double y, double z, double xi, double eta) {
  // double R ;

  return 1. / sqrt(x * x + y * y + z * z - 2 * y * eta + eta * eta -
                   2 * x * xi + xi * xi);
}

double ip13(double x, double y, double z, double xi, double eta) {
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return z / ((R + y - eta) * R);
}

double ip22(double x, double y, double z, double xi, double eta) {
  //  (y - \[Eta])/(Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] (x + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] - \[Xi]))
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (y - eta) / ((R + x - xi) * R);
}

double ip23(double x, double y, double z, double xi, double eta) {
  //  z/(Sqrt[z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] (x + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] - \[Xi]))
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return z / ((R + x - xi) * R);
}

double ip33(double x, double y, double z, double xi, double eta) {
  //  ((y - \[Eta]) (2 z^2 + (y - \[Eta])^2 + (x - \[Xi])^2) (x - \
//\[Xi]))/((z^2 + (y - \[Eta])^2) (z^2 + (x - \[Xi])^2) Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2])
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (x - xi) * (y - eta) *
         (2 * z * z + (y - eta) * (y - eta) + (xi - x) * (xi - x)) /
         (R * (z * z + (x - xi) * (x - xi)) * (z * z + (y - eta) * (y - eta)));
}

//// third order derivatives of I(x,y,z,xi,eta)

double ip111(double x, double y, double z, double xi, double eta) {
  //  (R2 (Sqrt[R2] + y - \[Eta]) -
  //      Sqrt[R2] (x - \[Xi])^2 - (Sqrt[R2] + y - \[Eta]) (x - \[Xi])^2)/(R2^(
  //      3/2) (Sqrt[R2] + y - \[Eta])^2)

  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return (R * (R2 - 2. * ((x - xi) * (x - xi))) +
          (y - eta) * (R2 - ((x - xi) * (x - xi)))) /
         (pow(R2, 3. / 2.) * ((R + y - eta) * (R + y - eta)));
}

double ip112(double x, double y, double z, double xi, double eta) {
  //  (-x + \[Xi])/(z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^(3/2)
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return (xi - x) / pow(R2, 3. / 2.);
}

double ip113(double x, double y, double z, double xi, double eta) {
  //-((z (y - \[Eta] +
  //  2 Sqrt[z^2 + (y - \[Eta])^2 + (x - \[Xi])^2]) (x - \[Xi]))/((y - \
//\[Eta] + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2])^2 (z^2 + (y - \[Eta])^2 + \
//(x - \[Xi])^2)^(3/2)))

  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return z * (xi - x) * (2. * R + y - eta) /
         (pow(R2, 3. / 2.) * ((R + y - eta) * (R + y - eta)));
}

double ip122(double x, double y, double z, double xi, double eta) {
  //(-y + \[Eta])/(z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^(3/2)
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return (eta - y) / pow(R2, 3. / 2.);
}

double ip123(double x, double y, double z, double xi, double eta) {
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return -z / pow(R2, 3. / 2.);
}

double ip133(double x, double y, double z, double xi, double eta) {
  //  (R (R2 - 2 z^2) + (R2 - z^2) (y - \[Eta]))/(R2^(
  //      3/2) (R + y - \[Eta])^2)

  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return (R * (R2 - 2. * z * z) + (R2 - z * z) * (y - eta)) /
         (pow(R2, 3. / 2.) * ((R + y - eta) * (R + y - eta)));
}

double ip222(double x, double y, double z, double xi, double eta) {
  //  (R (R2 - 2 (y - \[Eta])^2) + (R2 - (y - \[Eta])^2) (x - \[Xi]))/(R2^(
  //      3/2) (R + x - \[Xi])^2)
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return (R * (R2 - 2. * ((y - eta) * (y - eta))) +
          (x - xi) * (R2 - (y - eta) * (y - eta))) /
         (pow(R2, 3. / 2.) * ((R + x - xi) * (R + x - xi)));
}

double ip223(double x, double y, double z, double xi, double eta) {
  //  -((z (y - \[Eta]) (2 R + x - \[Xi]))/(R2^(3/2) (R + x - \[Xi])^2))

  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return z * (eta - y) * (2 * R + x - xi) /
         (pow(R2, 3. / 2.) * ((R + x - xi) * (R + x - xi)));
}

double ip233(double x, double y, double z, double xi, double eta) {
  //  (R (R2 - 2 z^2) + (R2 - z^2) (x - \[Xi]))/(R2^(3/2) (R + x - \[Xi])^2)
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return (R * (R2 - 2. * z * z) + (x - xi) * (R2 - z * z)) /
         (pow(R2, 3. / 2.) * ((R + x - xi) * (R + x - xi)));
}

double ip333(double x, double y, double z, double xi, double eta) {
  //    (-(2 R2^2 + R2 z^2 + 3 z^4) (z^2 + (y - \[Eta])^2) (y - \[Eta]) (z^2 + (x - \
    //    \[Xi])^2) (x - \[Xi]) + 2 (y - \[Eta])^3 (2 z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^2 (x - \
    //    \[Xi])^3)/(R2^(3/2) z (z^2 + (y - \[Eta])^2)^2 (z^2 + (x - \[Xi])^2)^2)

  double R2, R, z2, xmxi, ymeta;
  xmxi = (x - xi);
  ymeta = (y - eta);
  R2 = xmxi * xmxi + ymeta * ymeta + z * z;
  z2 = z * z;

  return (-(2 * R2 * R2 + R2 * z2 + 3 * z2 * z2) * (z2 + pow(ymeta, 2)) *
              ymeta * (z2 + pow(xmxi, 2)) * xmxi +
          2 * pow(ymeta, 3) * pow(2 * z2 + pow(ymeta, 2) + pow(xmxi, 2), 2) *
              pow(xmxi, 3)) /
         (pow(R2, 1.5) * z * pow(z2 + pow(ymeta, 2), 2) *
          pow(z2 + pow(xmxi, 2), 2));
}

// CHEMMERY Integration function

typedef double (*vFunctionCall)(double x, double y, double z, double xi,
                                double eta);

double rectangular_integration(double x, double y, double z, double a, double b,
                               vFunctionCall Func) {
  double ma, mb;
  ma = -a;
  mb = -b;
  return (Func(x, y, z, a, b) - Func(x, y, z, a, mb) - Func(x, y, z, ma, b) +
          Func(x, y, z, ma, mb));
}

//
//// Fundamental stress kernel - only  effect of DD_x(e_1) and DD_z (e3)
//// this function output the 3D stress in the 3D frame for completeness ?
// il::StaticArray2D<double, 2, 6> all_3d_stresses_kernel_s3d_p0_dd(
//    double x, double y, double z, double a, double b, double G, double nu) {
//  // x , y , z location where to compute stress
//  //  a,b  1/2 size of the rectangular DD
//  //  G Shear modulus, nu Poisson's ratio'
//  // Ep instead of G ?
//  //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
//  //  DDx (shear), DDy (shear), DDz (normal)
//
//  //  CONVENTION: Positive DDs in Overlap, Stresses positive in tension
//
//------------------------------------------------------------------------------

//  this function ouptuts the stress in the plane (x2=0) s11,s33,s13
//  simplified 3D kernel
il::StaticArray2D<double, 2, 3> stresses_kernel_s3d_p0_dd(double a, double b,
                                                          double G, double nu,
                                                          double xx,
                                                          double yy) {
  //  Simplified 3D Kernel - 3D kernel of Rongved,
  // with corrections from Wu & Olson 2015
  //
  //  Rectangular DD on plan z=0   x \in [-a,a], y \in [-b,b]
  //  DDx (shear),  DDz (normal)
  //  CONVENTION: Positive DD in Overlap, stress positive in tension
  //
  // Inputs:
  // a 1/2 size of the rectangular DD,
  // b corresponds to the 1/2 FRACTURE HEIGHT
  // G Shear modulus, nu Poisson's ratio'
  // xx, yy location where to compute stress in the 2D elastic plane
  //
  // Outputs:
  // 2 by 3 array  - stress components in  the 2D xy frame
  // 1st row contains the stresses due to the shear DD : sxx, sxy, syy
  // 2nd row contains the stresses due to the normal DD : sxx, sxy, syy
  //
  // CONVENTIONS:
  // positive compression / positive opening DD;
  // positive tension / positive overlap DD;

  // Note
  // switch to the 3D cartesian frame x,y,z (of the DD element)
  // with the plane of the DD being z=0
  // the simplified 3D kernel then looks at stresses
  // in the mid plane of the element y=0.

  double x = xx;
  double z = yy;
  double y = 0.;

  double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;
  double Ip111, Ip122, Ip133, Ip112, Ip113, Ip123, Ip222, Ip223, Ip233, Ip333;
  double Ce = 1. * G / (4 * il::pi * (1. - nu));
  double C11, C12, D, H;

  D = sqrt(xx * xx + yy * yy); // distance to DD center
  H = 2 * b; // height -> may want to modify that & input H directly ?

  // correction factors
  C11 =
      1. - pow(D, 1.5) * pow(H / 1., 4.) / pow(D * D + H * H, (1.5 + 4.) / 2.);
  C12 =
      1. - pow(D, 1.8) * pow(H / 1., 4.) / pow(D * D + H * H, (1.8 + 4.) / 2.);

  //
  il::StaticArray2D<double, 2, 3> StressInPlane;
  // compute the Is function derivatives....

  //  Ip11 = rectangular_integration(x, y, z, a, b, ip11);
  Ip22 = rectangular_integration(x, y, z, a, b, ip22);
  Ip33 = rectangular_integration(x, y, z, a, b, ip33);
  //  Ip23 = rectangular_integration(x, y, z, a, b, ip23);
  //  Ip12 = rectangular_integration(x, y, z, a, b, ip12);
  Ip13 = rectangular_integration(x, y, z, a, b, ip13);

  Ip111 = rectangular_integration(x, y, z, a, b, ip111);
  //  Ip122 = rectangular_integration(x, y, z, a, b, ip122);
  Ip133 = rectangular_integration(x, y, z, a, b, ip133);
  //  Ip112 = rectangular_integration(x, y, z, a, b, ip112);
  Ip113 = rectangular_integration(x, y, z, a, b, ip113);
  //  Ip123 = rectangular_integration(x, y, z, a, b, ip123);
  //  Ip222 = rectangular_integration(x, y, z, a, b, ip222);
  //  Ip233 = rectangular_integration(x, y, z, a, b, ip233);
  //  Ip223 = rectangular_integration(x, y, z, a, b, ip223);
  Ip333 = rectangular_integration(x, y, z, a, b, ip333);

  // StressInPlane row is dof (DDshear,DDnormal),
  // columns are sxx,sxy,syy (in 2D plane frame)

  // stress due to displacement discontinuity DDx (shear)
  // sxx
  StressInPlane(0, 0) = Ce * C11 * (2. * Ip13 - z * Ip111);
  //  Stress(0, 1) = Ce * (2. * nu * Ip13 - z * Ip122);   // syy
  //  Stress(0, 3) = Ce * ((1. - nu) * Ip23 - z * Ip112);  // sxy
  // sxz (3) -> Correspond to sxy in 2d Plan elasticity frame
  StressInPlane(0, 1) = Ce * C12 * (Ip33 + nu * Ip22 - z * Ip113);
  //  Stress(0, 5) = Ce * (-nu * Ip12 - z * Ip123);        // syz
  // szz (3D)  -> correspond to syy in 2d Plan elasticity frame
  StressInPlane(0, 2) = Ce * C11 * (-z * Ip133);

  // stress du to displacement discontinuity DDz (normal)
  StressInPlane(1, 0) =
      Ce * C11 * (Ip33 + (1. - 2. * nu) * Ip22 - z * Ip113); // sxx
  //  Stress(1, 3) = Ce * (-(1. - 2. * nu) * Ip12 - z * Ip123);        // sxy
  StressInPlane(1, 1) = Ce * C12 * (-z * Ip133); // sxz[which is sxy in the 2D]
  //  Stress(1, 5) = Ce * (z * Ip233);                                 // syz
  //  Stress(1, 1) = Ce * (Ip33 + (1. - 2. * nu) * Ip11 - z * Ip223);  // syy
  double EPSILON;
  EPSILON = std::numeric_limits<double>::epsilon();
  if (abs(z) <= EPSILON) {
    StressInPlane(1, 2) = Ce * C11 * (Ip33);
  } else {
    StressInPlane(1, 2) = Ce * C11 * (Ip33 - z * Ip333); // modified 2021
  }
  // szz[which is syy in the 2D]
  return StressInPlane;
}

////------------------------------------------------------------------------------
//il::StaticArray<double, 3> point_stress_s3d_dp0_dd(
//    il::StaticArray<double, 2> &observ_pt, SegmentData &source_elt,
//    il::StaticArray<double, 4> &nodal_dd, // we use only [0] and [1]
//    ElasticProperties &Elas, double ker_options) {
//  // Function to get the stress at a point
//  // induced by a piece-wise constant DD segment of size h and height
//  // ker_options
//  //
//  // Used to calculate stresses at given points
//  //
//  // INPUTS:
//  // source_elt : element data structure of the source element
//  // elt_dd : (shear, normal)
//  // observ_pt : observation point coordinates (x, y)
//  // Elas :: elastic properties object
//  // ker_options : dummy argument here for the height of the simplified 3D elt
//  //
//  // OUTPUTS:
//  // stress in global (reference) coordinates (xx, xy, yy)
//  // Make sure to add in-situ stress later
//  //
//  // CONVENTIONS:
//  // positive compression / positive opening DD;
//  // positive tension / positive overlap DD;
//
//  // the tensor to switch to the frame of the source element....
//  il::StaticArray2D<double, 2, 2> R = bigwham::rotationMatrix2D(source_elt.theta());
//
//  // the inverse rotation...
//  il::StaticArray2D<double, 2, 2> Rt = R;
//  Rt(0, 1) = R(1, 0);
//  Rt(1, 0) = R(0, 1); // il::transpose() ?
//
//  // directional vector from the element to the observation point
//  il::StaticArray<double, 2> xe;
//  for (int i = 0; i < 2; ++i) {
//    xe[i] = observ_pt[i] - source_elt.Xmid(i);
//  }
//  // -"- in source element's coordinates
//  xe = il::dot(Rt, xe);
//
//  // stress kernel at the observation point in source element's coordinates
//  il::StaticArray2D<double, 2, 3> stress_l = stresses_kernel_s3d_p0_dd(
//      source_elt.size() / 2., ker_options / 2., Elas.shear_modulus(),
//      Elas.poisson_ratio(), xe[0], xe[1]);
//
//  // stresses at the observation point (in source element's coordinates)
//  il::StaticArray<double, 2> elt_dd;
//  elt_dd[0] = nodal_dd[0];
//  elt_dd[1] = nodal_dd[1]; // here we ignore [2] & [3] components (for P1)
//  il::StaticArray<double, 3> stress_g = il::dot(elt_dd, stress_l);
//
//  // stress tensor in matrix form (in source element's coordinates)
//  il::StaticArray2D<double, 2, 2> stress_m, stress_i;
//  stress_m(0, 0) = stress_g[0];
//  stress_m(0, 1) = stress_g[1];
//  stress_m(1, 0) = stress_g[1];
//  stress_m(1, 1) = stress_g[2];
//
//  // stress rotation to global (reference) coordinates
//  stress_i = il::dot(stress_m, Rt);
//  stress_m = il::dot(R, stress_i);
//
//  stress_g[0] = stress_m(0, 0);
//  stress_g[1] = 0.5 * (stress_m(0, 1) + stress_m(1, 0));
//  stress_g[2] = stress_m(1, 1);
//
//  return stress_g;
//};
} // namespace bigwham
