//
// This file is part of HFPx3D.
//
// Created by Brice Lecampion on 04.02.19.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2019.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//


#include <il/math.h>

#include "R0_element.h"


// RONGVED SOLUTION FOR A P0 Rectangular dislocation in a full space
// dislocation is centered on the origin in the plane z=0 , (-a,a) in x (-b,b)
// in y

// first order derivatives of I(x,y,z,xi,eta)
double ip1(double& x, double& y, double& z, double& xi, double& eta) {
    double R;
    R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

    return log(R + y - eta);
}

double ip2(double& x, double& y, double& z, double& xi, double& eta) {
    double R;
    R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

    return log(R + x - xi);
}

double ip3(double& x, double& y, double& z, double& xi, double& eta) {
    double R;
    R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

    return -atan((x - xi) * (y - eta) / (z * R));
}

// second order derivatives of I(x,y,z,xi,eta)
double ip11(double& x, double& y, double& z, double& xi, double& eta) {
//    (x - \[Xi])/((y - \[Eta] + Sqrt(Power(z,2) + Power(y - \[Eta],2) + Power(x
//    - \[Xi],2)))*
//        Sqrt(Power(z,2) + Power(y - \[Eta],2) + Power(x - \[Xi],2)))
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (x - xi) / ((R + y - eta) * R);
}

double ip12(double& x, double& y, double& z, double& xi, double& eta) {
  // double R ;

  return 1. / sqrt(x * x + y * y + z * z - 2 * y * eta + eta * eta -
      2 * x * xi + xi * xi);
}

double ip13(double& x, double& y, double& z, double& xi, double& eta) {
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return z / ((R + y - eta) * R);
}

double ip22(double& x, double& y, double& z, double& xi, double& eta) {
  //  (y - \[Eta])/(Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] (x + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] - \[Xi]))
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (y - eta) / ((R + x - xi) * R);
}

double ip23(double& x, double& y, double& z, double& xi, double& eta) {
  //  z/(Sqrt[z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] (x + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] - \[Xi]))
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return z / ((R + x - xi) * R);
}

double ip33(double& x, double& y, double& z, double& xi, double& eta) {
  //  ((y - \[Eta]) (2 z^2 + (y - \[Eta])^2 + (x - \[Xi])^2) (x - \
	//\[Xi]))/((z^2 + (y - \[Eta])^2) (z^2 + (x - \[Xi])^2) Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2])
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (x - xi) * (y - eta) *
      (2 * z * z + (y - eta) * (y - eta) + (xi - x) * (xi - x)) /
      (R * (z * z + (x - xi) * (x - xi)) * (z * z + (y - eta) * (y - eta)));
}

//// third order derivatives of I(x,y,z,xi,eta)

double ip111(double& x, double& y, double& z, double& xi, double& eta) {
  //  (R2 (Sqrt[R2] + y - \[Eta]) -
  //      Sqrt[R2] (x - \[Xi])^2 - (Sqrt[R2] + y - \[Eta]) (x - \[Xi])^2)/(R2^(
  //      3/2) (Sqrt[R2] + y - \[Eta])^2)

  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return (R * (R2 - 2. * pow(x - xi, 2)) + (y - eta) * (R2 - pow(x - xi, 2))) /
      (pow(R2, 3. / 2.) * pow(R + y - eta, 2.));
}

double ip112(double& x, double& y, double& z, double& xi, double& eta) {
  //  (-x + \[Xi])/(z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^(3/2)
  double R2;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return (xi - x) / pow(R2, 3. / 2.);
}

double ip113(double& x, double& y, double& z, double& xi, double& eta) {
  //-((z (y - \[Eta] +
  //  2 Sqrt[z^2 + (y - \[Eta])^2 + (x - \[Xi])^2]) (x - \[Xi]))/((y - \
	//\[Eta] + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2])^2 (z^2 + (y - \[Eta])^2 + \
	//(x - \[Xi])^2)^(3/2)))

  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return z * (xi - x) * (2. * R + y - eta) /
      (pow(R2, 3. / 2.) * pow(R + y - eta, 2));
}

double ip122(double& x, double& y, double& z, double& xi, double& eta) {
  //(-y + \[Eta])/(z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^(3/2)
  double R2;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return (eta - y) / pow(R2, 3. / 2.);
}

double ip123(double& x, double& y, double& z, double& xi, double& eta) {
  double R2;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return -z / pow(R2, 3. / 2.);
}

double ip133(double& x, double& y, double& z, double& xi, double& eta) {
  //  (R (R2 - 2 z^2) + (R2 - z^2) (y - \[Eta]))/(R2^(
  //      3/2) (R + y - \[Eta])^2)

  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return (R * (R2 - 2. * z * z) + (R2 - z * z) * (y - eta)) /
      (pow(R2, 3. / 2.) * pow(R + y - eta, 2.));
}

double ip222(double& x, double& y, double& z, double& xi, double& eta) {
  //  (R (R2 - 2 (y - \[Eta])^2) + (R2 - (y - \[Eta])^2) (x - \[Xi]))/(R2^(
  //      3/2) (R + x - \[Xi])^2)
  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return (R * (R2 - 2. * pow(y - eta, 2.)) +
      (x - xi) * (R2 - (y - eta) * (y - eta))) /
      (pow(R2, 3. / 2.) * pow(R + x - xi, 2.));
}

double ip223(double& x, double& y, double& z, double& xi, double& eta) {
  //  -((z (y - \[Eta]) (2 R + x - \[Xi]))/(R2^(3/2) (R + x - \[Xi])^2))

  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return z * (eta - y) * (2 * R + x - xi) /
      (pow(R2, 3. / 2.) * pow(R + x - xi, 2.));
}

double ip233(double& x, double& y, double& z, double& xi, double& eta) {
  //  (R (R2 - 2 z^2) + (R2 - z^2) (x - \[Xi]))/(R2^(3/2) (R + x - \[Xi])^2)
  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return (R * (R2 - 2. * z * z) + (x - xi) * (R2 - z * z)) /
      (pow(R2, 3. / 2.) * pow(R + x - xi, 2.));
}

// CHEMMERY Integration function

typedef double (*vFunctionCall)(double& x, double& y, double& z, double& xi,
                                double& eta);

double rectangular_integration(double& x, double& y, double& z, double& a,
                               double& b, vFunctionCall Func) {
  double ma, mb;
  ma = -a;
  mb = -b;
  return (Func(x, y, z, a, b) - Func(x, y, z, a, mb) - Func(x, y, z, ma, b) +
      Func(x, y, z, ma, mb));
}

// Fundamental stress kernel
il::StaticArray2D<double, 3, 6> StressesKernelRectangularP0DD(
    double& x, double& y, double& z, double& a, double& b, double& G,
    double& nu) {
  // x , y , z location where to compute stress
  //  a,b  1/2 size of the rectangular DD
  //  G Shear modulus, nu Poisson's ratio'
  //  Ep instead of G ?
  //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
  //  DDx (shear), DDy (shear), DDz (normal)

  double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;

  double Ip111, Ip122, Ip133, Ip112, Ip113, Ip123, Ip222, Ip223, Ip233;

  double Ce = G / (4 * il::pi * (1. - nu));
  //  double sxx, sxy, sxz, syy, syz, szz;
  //
  il::StaticArray2D<double, 3, 6> Stress;
  // compute the Is function derivatives....

  Ip11 = rectangular_integration(x, y, z, a, b, ip11);
  Ip22 = rectangular_integration(x, y, z, a, b, ip22);
  Ip33 = rectangular_integration(x, y, z, a, b, ip33);
  Ip23 = rectangular_integration(x, y, z, a, b, ip23);
  Ip12 = rectangular_integration(x, y, z, a, b, ip12);
  Ip13 = rectangular_integration(x, y, z, a, b, ip13);

  Ip111 = rectangular_integration(x, y, z, a, b, ip111);
  Ip122 = rectangular_integration(x, y, z, a, b, ip122);
  Ip133 = rectangular_integration(x, y, z, a, b, ip133);
  Ip112 = rectangular_integration(x, y, z, a, b, ip112);
  Ip113 = rectangular_integration(x, y, z, a, b, ip113);
  Ip123 = rectangular_integration(x, y, z, a, b, ip123);
  Ip222 = rectangular_integration(x, y, z, a, b, ip222);
  Ip233 = rectangular_integration(x, y, z, a, b, ip233);
  Ip223 = rectangular_integration(x, y, z, a, b, ip223);

  // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

  // stress due to displacement discontinuity DDx (shear)
  Stress(0, 0) = Ce * (2. * Ip13 - z * Ip111);         // sxx
  Stress(0, 1) = Ce * (2. * nu * Ip13 - z * Ip122);    // syy
  Stress(0, 2) = Ce * (-z * Ip133);                    // szz
  Stress(0, 3) = Ce * ((1. - nu) * Ip23 - z * Ip112);  // sxy
  Stress(0, 4) = Ce * (Ip33 + nu * Ip22 - z * Ip113);  // sxz
  Stress(0, 5) = Ce * (-nu * Ip12 - z * Ip123);        // syz

  // stress due to displacement discontinuity  DDy (shear)
  Stress(1, 0) = Ce * (2. * nu * Ip23 - z * Ip112);    // sxx
  Stress(1, 1) = Ce * (2. * Ip23 - z * Ip222);         // syy
  Stress(1, 2) = Ce * (-z * Ip233);                    // szz
  Stress(1, 3) = Ce * ((1. - nu) * Ip13 - z * Ip122);  // sxy
  Stress(1, 4) = Ce * (-nu * Ip12 - z * Ip123);        // sxz
  Stress(1, 5) = Ce * (Ip33 + nu * Ip11 - z * Ip223);  // syz

  // stress due to displacement discontinuity DDz (normal)
  Stress(2, 0) = Ce * (Ip33 + (1. - 2. * nu) * Ip22 - z * Ip113);  // sxx
  Stress(2, 1) = Ce * (Ip33 + (1. - 2. * nu) * Ip11 - z * Ip223);  // syy
  Stress(2, 2) = Ce * (Ip33 - z * Ip113);                          // szz
  Stress(2, 3) = Ce * (-(1. - 2. * nu) * Ip12 - z * Ip123);        // sxy
  Stress(2, 4) = Ce * (-z * Ip133);                                // sxz
  Stress(2, 5) = Ce * (z * Ip233);                                 // syz

  return Stress;
}


// todo : implement mormalshearTraction routine to be used in BIE assembly

