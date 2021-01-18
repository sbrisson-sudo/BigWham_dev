//
// This file is part of HFPx3D.
//
// Created by Brice Lecampion on 04.02.19.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2019.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//


#include <cmath>
#include <il/linearAlgebra.h>
#include <il/blas.h>
#include "Elastic3DR0_element.h"
#include <iostream>

namespace bie{
    // RONGVED SOLUTION FOR A P0 Rectangular dislocation in a full space
    // dislocation is centered on the origin in the plane z=0 , (-a,a) in x (-b,b)
    // in y

    //--------------------------------------------------------------------------------------------------//
    //                                                                                                  //
    //  Special case to treat any case when z = 0.                                                      //
    //  i.e. when the point where we ask for the coordinates lies on the plane of the element z == 0    //
    //                                                                                                  //
    //  PLEASE: do note that any stress is theoretically singular when evaluated at the edges of an     //
    //  element or its prolongation.                                                                    //
    //                                                                                                  //
    //--------------------------------------------------------------------------------------------------//

    double ip11_lim_z_to_0(double& x, double& y, double& z, double& xi, double& eta) {
        //(x - \[Xi])/((y - \[Eta] + Sqrt[(y - \[Eta])^2 + (x - \[Xi])^2])
        //   Sqrt[(y - \[Eta])^2 + (x - \[Xi])^2])
        double y_m_eta, x_m_xi, my_sqrt_;
        y_m_eta = y - eta ;
        x_m_xi = x - xi ;
        my_sqrt_ = sqrt(y_m_eta * y_m_eta + x_m_xi * x_m_xi);

        return x_m_xi / (my_sqrt_ * (y + my_sqrt_ - eta));
    }

    double ip22_lim_z_to_0(double& x, double& y, double& z, double& xi, double& eta) {
        // (y - \[Eta])/(Sqrt[(y - \[Eta])^2 + (x - \[Xi])^2] (x +
        //   Sqrt[(y - \[Eta])^2 + (x - \[Xi])^2] - \[Xi]))
        double y_m_eta, x_m_xi, my_sqrt_;
        y_m_eta = y - eta ;
        x_m_xi = x - xi ;
        my_sqrt_ = sqrt(y_m_eta * y_m_eta + x_m_xi * x_m_xi);

        return y_m_eta / (my_sqrt_ * (x + my_sqrt_ - xi));
    }

    double ip33_lim_z_to_0(double& x, double& y, double& z, double& xi, double& eta) {
        // Sqrt[(y - \[Eta])^2 + (x - \[Xi])^2]/((y - \[Eta]) (x - \[Xi]))
        double y_m_eta, x_m_xi;
        y_m_eta = y - eta ;
        x_m_xi = x - xi ;

        return sqrt(y_m_eta * y_m_eta + x_m_xi * x_m_xi) / (y_m_eta * x_m_xi);
    }

    double Ip33_lim_z_to_0_and_x_to_a(double& x, double& y, double& a, double& b) {
        // -(Sqrt[(a + x)^2 + (-b + y)^2]/((a + x) (-b + y))) + Sqrt[(a +
        //    x)^2 + (b + y)^2]/((a + x) (b + y))
        double a_plus_x, b_plus_y, y_minus_b, sqrt_1st, sqrt_2nd;
        a_plus_x = a + x ;
        b_plus_y = b + y ;
        y_minus_b = y - b ;
        sqrt_1st = sqrt(a_plus_x * a_plus_x + y_minus_b * y_minus_b);
        sqrt_2nd = sqrt(a_plus_x * a_plus_x + b_plus_y * b_plus_y);

        return - sqrt_1st / (a_plus_x * y_minus_b) + sqrt_2nd / (a_plus_x * b_plus_y);
    }

    double Ip33_lim_z_to_0_and_y_to_b(double& x, double& y, double& a, double& b) {
        // -(Sqrt[(-a + x)^2 + (b + y)^2]/((-a + x) (b + y))) + Sqrt[(a +
        //    x)^2 + (b + y)^2]/((a + x) (b + y))
        double a_plus_x, b_plus_y, x_minus_a, sqrt_1st, sqrt_2nd;
        a_plus_x = a + x ;
        b_plus_y = b + y ;
        x_minus_a = x - a ;
        sqrt_1st = sqrt(b_plus_y * b_plus_y + x_minus_a * x_minus_a);
        sqrt_2nd = sqrt(a_plus_x * a_plus_x + b_plus_y * b_plus_y);

        return - sqrt_1st / (b_plus_y * x_minus_a) + sqrt_2nd / (a_plus_x * b_plus_y);
    }

    double ip12_lim_z_to_0(double& x, double& y, double& z, double& xi, double& eta) {
        return 1. / sqrt( (y - eta) * (y - eta) + (x - xi) * (x - xi));
    }

    //--------------------------------------------------------------------------------------------------//

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

    double ip333(double& x, double& y, double& z, double& xi, double& eta) {
    //    (-(2 R2^2 + R2 z^2 + 3 z^4) (z^2 + (y - \[Eta])^2) (y - \[Eta]) (z^2 + (x - \
    //    \[Xi])^2) (x - \[Xi]) + 2 (y - \[Eta])^3 (2 z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^2 (x - \
    //    \[Xi])^3)/(R2^(3/2) z (z^2 + (y - \[Eta])^2)^2 (z^2 + (x - \[Xi])^2)^2)

        double R2, z2, xmxi, ymeta;
        xmxi = (x - xi);
        ymeta = (y - eta);
        R2 = xmxi * xmxi + ymeta * ymeta + z * z;
        z2 = z * z;

        return (-(2 * R2 * R2 + R2 * z2 + 3 * z2 * z2) *
                 (z2 + pow(ymeta,2)) * ymeta *
                 (z2 + pow(xmxi,2)) * xmxi + 2 * pow(ymeta,3) *
                 pow(2 * z2 + pow(ymeta,2) + pow(xmxi,2),2) *
                 pow(xmxi,3)) / (pow(R2,1.5) * z * pow(z2 + pow(ymeta,2),2) *
                 pow(z2 + pow(xmxi,2),2));
    }

    // Chinnery Integration function - M. A. Chinnery, The deformations of the ground around surface faults, Bulletin of the Seismological Society of America (1961), 51 355-372

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
    il::StaticArray2D<double, 3, 6> StressesKernelR0_NON_NULL_z(
        double& x, double& y, double& z, double& a, double& b, double& G,
        double& nu) {
      // x , y , z location where to compute stress
      //  a,b  1/2 size of the rectangular DD
      //  G Shear modulus, nu Poisson's ratio'
      //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
      //  DDx (shear), DDy (shear), DDz (normal)

      double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;

      double Ip111, Ip122, Ip133, Ip112, Ip113, Ip123, Ip222, Ip223, Ip233, Ip333;

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
      Ip333 = rectangular_integration(x, y, z, a, b, ip333);

      // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

      // stress due to displacement discontinuity DDx (shear)
      Stress(0, 0) = Ce * (2. * Ip13 - z * Ip111);         // sxx
      Stress(0, 1) = Ce * (2.   * nu * Ip13 - z * Ip122);  // syy
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
      Stress(2, 2) = Ce * (Ip33 - z * Ip333);                          // szz  (fixed by CP 2020)
      Stress(2, 3) = Ce * (-(1. - 2. * nu) * Ip12 - z * Ip123);        // sxy
      Stress(2, 4) = Ce * (-z * Ip133);                                // sxz
      Stress(2, 5) = Ce * (-z * Ip233);                                // syz (a minus by CP 2020)

      return Stress;
            // DDx (shear)  -> | sxx, syy, szz, sxy, sxz, syz  |
            // DDy (shear)  -> | sxx, syy, szz, sxy, sxz, syz  |
            // DDz (normal) -> | sxx, syy, szz, sxy, sxz, syz  |
    }

    bool is_stress_singular_at_given_location(double& x, double& y, double& a, double& b, bool verbose)
    {   double EPSILON;
        EPSILON = std::numeric_limits<double>::epsilon();
        if (
                ( abs(abs(x)/a - 1.) <= EPSILON && (abs(y) <= b) )  ||
                ( abs(abs(y)/b - 1.) <= EPSILON && (abs(x) <= a ))
                )
        {   if (verbose){
                std::cout << "WARNING: \n " \
                      << " you are computing the stress along the edge of a 3DR0 element. \n" \
                      << " At that location some components of the stress tensor are theoretically infinite. \n" \
                      << " Suggestions: \n" \
                      << "    1) check that your collocation points do not lie on the edge of another element \n" \
                      << "    2) check that you are not computing the stress on the edge of an element \n";
                }
            return true ;
        }
        else return false;
    }


    il::StaticArray2D<double, 3, 6> StressesKernelR0_NULL_z(
            double& x, double& y, double& a, double& b, double& G,
            double& nu) {
        // x , y , (z = 0) location where to compute stress
        //  a,b  1/2 size of the rectangular DD
        //  G Shear modulus, nu Poisson's ratio'
        //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
        //  DDx (shear), DDy (shear), DDz (normal)

        bool answer;
        il::StaticArray2D<double, 3, 6> Stress;
        double EPSILON;
        EPSILON = std::numeric_limits<double>::epsilon();

        if (!is_stress_singular_at_given_location(x, y, a, b))
        {
            double z = 0.;

            double Ce = G / (4 * il::pi * (1. - nu));
            //  double sxx, sxy, sxz, syy, syz, szz;
            //

            // compute the Is function derivatives....

            double Ip11_lim_z_to_0, Ip22_lim_z_to_0, Ip33_lim_z_to_0, Ip12_lim_z_to_0;
            Ip11_lim_z_to_0 = rectangular_integration(x, y, z, a, b, ip11_lim_z_to_0);
            Ip22_lim_z_to_0 = rectangular_integration(x, y, z, a, b, ip22_lim_z_to_0);

            if (abs(abs(x) / a - 1.) <= EPSILON && (abs(y) > b))
                Ip33_lim_z_to_0 = Ip33_lim_z_to_0_and_x_to_a(x, y, a, b);
            else if (abs(abs(y) / b - 1.) <= EPSILON && (abs(x) > a))
                Ip33_lim_z_to_0 = Ip33_lim_z_to_0_and_y_to_b(x, y, a, b);
            else
                Ip33_lim_z_to_0 = rectangular_integration(x, y, z, a, b, ip33_lim_z_to_0);

            Ip12_lim_z_to_0 = rectangular_integration(x, y, z, a, b, ip12_lim_z_to_0);

            // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

            // stress due to displacement discontinuity DDx (shear)
            Stress(0, 0) = 0.;                                              // sxx
            Stress(0, 1) = 0.;                                              // syy
            Stress(0, 2) = 0.;                                              // szz
            Stress(0, 3) = 0.;                                              // sxy
            Stress(0, 4) = Ce * (Ip33_lim_z_to_0 + nu * Ip22_lim_z_to_0);   // sxz
            Stress(0, 5) = Ce * (-nu) * Ip12_lim_z_to_0;                    // syz

            // stress due to displacement discontinuity  DDy (shear)
            Stress(1, 0) = 0.;                                              // sxx
            Stress(1, 1) = 0.;                                              // syy
            Stress(1, 2) = 0.;                                              // szz
            Stress(1, 3) = 0.;                                              // sxy
            Stress(1, 4) = Ce * (-nu) * Ip12_lim_z_to_0;                    // sxz
            Stress(1, 5) = Ce * (Ip33_lim_z_to_0 + nu * Ip11_lim_z_to_0);   // syz

            // stress due to displacement discontinuity DDz (normal)
            Stress(2, 0) = Ce * (Ip33_lim_z_to_0 + (1 - 2 * nu) * Ip22_lim_z_to_0);  // sxx
            Stress(2, 1) = Ce * (Ip33_lim_z_to_0 + (1 - 2 * nu) * Ip11_lim_z_to_0);  // syy
            Stress(2, 2) = Ce * Ip33_lim_z_to_0;                                     // szz
            Stress(2, 3) = Ce * (-1 + 2 * nu) * Ip12_lim_z_to_0;                     // sxy
            Stress(2, 4) = 0.;                                                       // sxz
            Stress(2, 5) = 0.;                                                       // syz
        }
        else {  for (il::int_t i = 0; i < 3; i++) {
                    for (il::int_t j = 0; j < 6; j++) {
                        Stress(i,j) = NAN;
                     }
                 }
        }
        return Stress;
        // DDx (shear)  -> | sxx, syy, szz, sxy, sxz, syz  |
        // DDy (shear)  -> | sxx, syy, szz, sxy, sxz, syz  |
        // DDz (normal) -> | sxx, syy, szz, sxy, sxz, syz  |
    }

    il::StaticArray2D<double, 3, 6> StressesKernelR0(
            double& x, double& y, double& z, double& a, double& b, double& G,
            double& nu)
    {
        double EPSILON;
        EPSILON = std::numeric_limits<double>::epsilon();
        if ( abs(z) <= EPSILON ) {return StressesKernelR0_NULL_z(x, y, a, b, G, nu);}
        else                     {return StressesKernelR0_NON_NULL_z(x, y, z, a, b, G, nu);}
    }

    // Fundamental displacement kernel
    il::Array2D<double> DisplacementKernelR0(
            double& x, double& y, double& z, double& a, double& b,
            double& nu) {
        //  x , y , z location where to compute displacement
        //  a,b  1/2 size of the rectangular DD
        //  nu Poisson's ratio'
        //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
        //  DDx (shear), DDy (shear), DDz (normal)

        double Ip1, Ip2, Ip3;

        double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;

        double Ce = -1 / (8 * il::pi * (1. - nu));

        il::Array2D<double> Displacement{3,3,0.};
        // compute the Is function derivatives....
        Ip1 = rectangular_integration(x, y, z, a, b, ip1);
        Ip2 = rectangular_integration(x, y, z, a, b, ip2);
        Ip3 = rectangular_integration(x, y, z, a, b, ip3);

        Ip11 = rectangular_integration(x, y, z, a, b, ip11);
        Ip22 = rectangular_integration(x, y, z, a, b, ip22);
        Ip33 = rectangular_integration(x, y, z, a, b, ip33);
        Ip23 = rectangular_integration(x, y, z, a, b, ip23);
        Ip12 = rectangular_integration(x, y, z, a, b, ip12);
        Ip13 = rectangular_integration(x, y, z, a, b, ip13);

        // Displacement row is dof (DDx,DDy,DDx), columns are Ux,Uy,Uz in the local reference system

        // displacement due to displacement discontinuity DDx (shear)
        Displacement(0, 0) = Ce * (z * Ip11 - 2 * (1 - nu) * Ip3);  // Ux
        Displacement(1, 0) = Ce * (z * Ip12);                       // Uy
        Displacement(2, 0) = Ce * (z * Ip13 - (1 - 2 * nu) * Ip1);  // Uz

        // displacement due to displacement discontinuity  DDy (shear)
        Displacement(0, 1) = Ce * (z * Ip12);                       // Ux
        Displacement(1, 1) = Ce * (z * Ip22 - 2 * (1 - nu) * Ip3);  // Uy
        Displacement(2, 1) = Ce * (z * Ip23 - (1 - 2 * nu) * Ip2);  // Uz

        // displacement due to displacement discontinuity DDz (normal)
        Displacement(0, 2) = Ce * (z * Ip13 + (1 - 2 * nu) * Ip1);  // Ux
        Displacement(1, 2) = Ce * (z * Ip23 + (1 - 2 * nu) * Ip2);  // Uy
        Displacement(2, 2) = Ce * (z * Ip33 - 2 * (1 - nu) * Ip3);  // Uz

        return Displacement; // expressed in the reference system of the DD element
        // index        ->    DDx (shear)    DDy (shear)     DDz (normal)
        //   0      -> |       Ux,            Ux,             Ux            |
        //   1      -> |       Uy,            Uy,             Uy            |
        //   2      -> |       Uz,            Uz,             Uz            |
    }

    il::Array2D<double> transpose(il::Array2D<double> M){
        il::Array2D<double> MT{3,3,0.};
        for (il::int_t i = 0; i < 3; i++) {
            for (il::int_t j = 0; j < 3; j++) {
                MT(j,i) = M(i,j);
            }
        }
        return MT;
    }

    il::Array2D<double> change_ref_system (const il::Array2D<double>& linearApplication,il::int_t change_domain, il::int_t change_codomain, const il::Array2D<double>& RglobalTOlocal_domain, const il::Array2D<double>& RglobalTOlocal_codomain){
        // Description:
        // A linear application takes values from a domain and outputs values in a codomain.
        // This function changes the base in the domain, in the codomain or in bonth.
        //
        // Input:
        // linearApplication is a matrix that thought to be in a local reference system both in the domain and in the codomain
        // change_domain can be 0 (false) or 1 (true)
        //      - if false the domain will be expressed with respect to the local reference system to the source element
        //      - if true the domain will be expressed with respect to the global reference system
        // change_codomain can be 0 (false) or 1 (true)
        //      - if false the domain will be expressed with respect to the local reference system to the receiver element
        //      - if true the codomain will be expressed with respect to the global reference system
        //
        // RglobalTOlocal_domain it is a matrix that rotates a vector from the global reference system (r.s.) to the local r.s. of the source element
        // RglobalTOlocal_codomain it is a matrix that rotates a vector from the global reference system (r.s.) to the local r.s. of the receiver element
        //
        // Output:
        // A rotated matrix i.e. linear application
        il::Array2D<double> rotatedLinearApplication{3,3,0.};

        if (change_domain == 0 && change_codomain == 0)
        {   // False - False
            // Explanation with regard to the problem at hand:
            //    R(global to local codomain)*R(from local to global domain)*(M_localDD_localTraction_source) = (M_globalDD_localTraction_receiver)
            //    R(global to local codomain)*R(from local to global domain)*(M_localDD_localTraction_source) = (M_globalDD_localDisplacements_receiver)
            il::Array2D<double> RT = transpose(RglobalTOlocal_domain);
            rotatedLinearApplication = il::dot(RT,linearApplication);
            rotatedLinearApplication = il::dot(RglobalTOlocal_codomain,rotatedLinearApplication);
        }
        else if (change_domain == 0 && change_codomain == 1)
        {   // False - True
            // Explanation with regard to the problem at hand:
            //    R(from local domain to global)*(M_localDD_localTraction) = (M_globalDD_localTraction)
            //    R(from local domain to global)*(M_localDD_localDisplacement) = (M_globalDD_localDisplacement)
            il::Array2D<double> RT = transpose(RglobalTOlocal_domain);
            rotatedLinearApplication = il::dot(RT,linearApplication);

        }
        else if (change_domain == 1 && change_codomain == 1)
        {   // True - True
            // Explanation with regard to the problem at hand:
            //    R(from local domain to global)*(M_localDD_localTraction)*R(from global to local domain) = (M_globalDD_globalTraction)
            //    R(from local domain to global)*(M_localDD_localDisplacement)*R(from global to local domain) = (M_globalDD_globalDisplacement)
            rotatedLinearApplication = il::dot(linearApplication,RglobalTOlocal_domain);
            il::Array2D<double> RT = transpose(RglobalTOlocal_domain);
            rotatedLinearApplication = il::dot(RT,rotatedLinearApplication);
        }
        else if (change_domain == 1 && change_codomain == 0)
        {   // True - False
            // Explanation with regard to the problem at hand:
            //    R(global to local codomain)*R(from local domain to global)*(M_localDD_localTraction)*R(from global to local domain) = (M_globalDD_localTraction)
            //    R(global to local codomain)*R(from local domain to global)*(M_localDD_localDisplacement)*R(from global to local domain) = (M_globalDD_localDisplacement)

            rotatedLinearApplication = il::dot(linearApplication,RglobalTOlocal_domain);
            il::Array2D<double> RT = transpose(RglobalTOlocal_domain);
            rotatedLinearApplication = il::dot(RT,rotatedLinearApplication);
            rotatedLinearApplication = il::dot(RglobalTOlocal_codomain,rotatedLinearApplication);
        }
        else { std::cout << "ERROR: bad options given for switch in routine: change_ref_system = " << change_domain << "\n";}

        return rotatedLinearApplication;
    }

    il::StaticArray<double,2> get_a_and_b(il::Array2D <double>xv, double NoV) {
        /*
         * This function returns two values:
         *
         * a := half length of the 1st edge of an element
         * b := half length of the last edge of an element
         */
        il::StaticArray<double,2> a_and_b;

        // vec01 goes from vertex 0 to vertex 1
        // vec02 goes from vertex 0 to vertex NoV_ (vertex 2 in case of triangular element)
        il::StaticArray<double, 3> vec01, vec02;
        vec01[0] = xv(1, 0) - xv(0, 0);
        vec01[1] = xv(1, 1) - xv(0, 1);
        vec01[2] = xv(1, 2) - xv(0, 2);
        vec02[0] = xv(NoV - 1, 0) - xv(0, 0);
        vec02[1] = xv(NoV - 1, 1) - xv(0, 1);
        vec02[2] = xv(NoV - 1, 2) - xv(0, 2);

        double vec01norm=sqrt(il::dot(vec01, vec01)), vec02norm=sqrt(il::dot(vec02, vec02));

        a_and_b[0] =vec01norm/2.; // a
        a_and_b[1] =vec02norm/2.; // b

        return a_and_b ;}

    il::Array2D<double> NodeDDtriplet_to_CPtraction_influence(
        FaceData &elem_data_s, // source element
        FaceData &elem_data_r, // receiver element
        ElasticProperties const &elas_, // elastic properties
        il::int_t I_want_global_DD = 0,
        il::int_t I_want_global_traction = 0)
    {

        double G = elas_.getG(), nu = elas_.getNu();
        il::StaticArray<double,2> a_and_b = get_a_and_b(elem_data_s.getVertices(),elem_data_s.getNoV());
        double a = a_and_b[0], b = a_and_b[1];

        il::Array2D<double> el_cp_s;
        el_cp_s = elem_data_s.getCollocationPoints();

        il::Array2D<double> el_cp_r;
        el_cp_r = elem_data_r.getCollocationPoints();

        il::Array2D<double> R = elem_data_s.rotationMatrix();

        il::Array<double> dsr{3};
        for (int i = 0; i < 3; ++i) { dsr[i] = el_cp_r(0,i) - el_cp_s(0,i); }

        // dsr contains the component of the distance between the source and the receiver
        dsr = il::dot(R, dsr);

        il::StaticArray2D<double, 3, 6> Stress;
        double EPSILON;
        EPSILON = std::numeric_limits<double>::epsilon();
        if (abs(dsr[2]) > EPSILON){
            Stress = StressesKernelR0(dsr[0],
                                      dsr[1],
                                      dsr[2],
                                      a, b,
                                      G,nu);}
        else{
            Stress = StressesKernelR0_NULL_z(dsr[0],
                                              dsr[1],
                                              a, b,
                                              G,nu);
        }
        // in the reference system of the source element both in the domain and in the codomain
        // index        ->    0    1    2    3    4    5
        // DDx (shear)  -> | sxx, syy, szz, sxy, sxz, syz  |
        // DDy (shear)  -> | sxx, syy, szz, sxy, sxz, syz  |
        // DDz (normal) -> | sxx, syy, szz, sxy, sxz, syz  |

        // normal vector at the receiver location in the reference system of the source element
        il::Array<double> nr = elem_data_r.getNormal();
        nr = il::dot(R,nr);

        il::Array<double> traction_temp;
        il::Array2D<double> DDs_to_traction_local_local{3,3,0.0}, sigma_temp{3,3,0.0};

        for (int i = 0; i < 3; ++i)
        { //loop over the rows of Stress
                sigma_temp(0,0) = Stress(i,0); // sig_xx
                sigma_temp(0,1) = Stress(i,3); // sig_xy
                sigma_temp(0,2) = Stress(i,4); // sig_xz
                sigma_temp(1,0) = Stress(i,3); // sig_yx
                sigma_temp(1,1) = Stress(i,1); // sig_yy
                sigma_temp(1,2) = Stress(i,5); // sig_yz
                sigma_temp(2,0) = Stress(i,4); // sig_zx
                sigma_temp(2,1) = Stress(i,5); // sig_zy
                sigma_temp(2,2) = Stress(i,2); // sig_zz

                traction_temp = il::dot(sigma_temp, nr);
                for (int j = 0; j < 3; ++j) {
                    DDs_to_traction_local_local(j,i) = traction_temp[j];
                    // | t1/Dshear1   t1/Dshear2  t1/Dnormal |
                    // | t2/Dshear1   t2/Dshear2  t2/Dnormal |
                    // | t3/Dshear1   t3/Dshear2  t3/Dnormal |
                    // localDD & local traction
                    // in the reference system of the source element both in the domain and in the codomain
                }
        }

        return change_ref_system(DDs_to_traction_local_local, I_want_global_DD, I_want_global_traction, R, elem_data_r.rotationMatrix());

        // | t1/Dshear1   t1/Dshear2  t1/Dnormal |
        // | t2/Dshear1   t2/Dshear2  t2/Dnormal |
        // | t3/Dshear1   t3/Dshear2  t3/Dnormal |
        //
        // directions 1, 2 or 3
    }

    il::Array2D<double> NodeDDtriplet_to_CPdisplacement_influence(
            FaceData &elem_data_s, // source element
            FaceData &elem_data_r, // receiver element
            ElasticProperties const &elas_, // elastic properties
            il::int_t I_want_global_DD,
            il::int_t I_want_global_displacement)
    {

        il::StaticArray<double,2> a_and_b = get_a_and_b(elem_data_s.getVertices(),elem_data_s.getNoV());
        double a = a_and_b[0], b = a_and_b[1];
        double nu = elas_.getNu();

        il::Array2D<double> el_cp_s;
        el_cp_s = elem_data_s.getCollocationPoints();

        il::Array2D<double> el_cp_r;
        el_cp_r = elem_data_r.getCollocationPoints();

        il::Array2D<double> R;
        R = elem_data_s.rotationMatrix();

        il::Array<double> dsr{3};
        for (int i = 0; i < 3; ++i) { dsr[i] = el_cp_r(0,i) - el_cp_s(0,i); }

        // dsr contains the component of the distance between the source and the receiver
        dsr = il::dot(R, dsr);

        il::Array2D<double> DDs_to_Displacement_local_local =DisplacementKernelR0(dsr[0],
                                                                               dsr[1],
                                                                               dsr[2],
                                                                               a, b,
                                                                               nu);
        // index        ->    DDx (shear)    DDy (shear)     DDz (normal)
        //   0      -> |       Ux,            Ux,             Ux            |
        //   1      -> |       Uy,            Uy,             Uy            |
        //   2      -> |       Uz,            Uz,             Uz            |

        return change_ref_system(DDs_to_Displacement_local_local, I_want_global_DD, I_want_global_displacement, R, elem_data_r.rotationMatrix());

        // | U1/Dshear1   U1/Dshear2  U1/Dnormal |
        // | U2/Dshear1   U2/Dshear2  U2/Dnormal |
        // | U3/Dshear1   U3/Dshear2  U3/Dnormal |
        //
        // directions 1, 2 or 3
    }
}