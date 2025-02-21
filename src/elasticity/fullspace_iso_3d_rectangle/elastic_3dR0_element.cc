//
// This file part of BigWham
//
// Created by Brice Lecampion on 04.02.19.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE file for more details.
//

#include <cmath>
#include <limits>

#include "elasticity/fullspace_iso_3d_rectangle//elastic_3dR0_common.h"
#include "elasticity/fullspace_iso_3d_rectangle/elastic_3dR0_element.h"

namespace bigwham{

    // RONGVED SOLUTION FOR A P0 Rectangular dislocation in a full space
    // dislocation is centered on the origin in the plane z=0 , (-a,a) in x (-b,b)
    // in y

    il::StaticArray2D<double, 3, 6> StressesKernelR0(
            double& x, double& y, double& z, double& a, double& b, double& G,
            double& nu) {
        //  assume a reference system x' y' z' with origin at the center of a rectangular DD and with axis x', y' parallel to the element's edges
        //  x , y , z location (with respect to x', y', z') where to compute stress.
        //  a,b  1/2 size of the rectangular DD
        //  G Shear modulus, nu Poisson's ratio'
        //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
        //  DDx (shear), DDy (shear), DDz (normal)

        il::StaticArray2D<double, 3, 6> Stress;
        double EPSILON;
        EPSILON = 100000 * std::numeric_limits<double>::epsilon();

        if (!is_stress_singular_at_given_location(x, y, z ,a, b, true))
        {

            double Ce = G / (4. * il::pi * (1. - nu));
            //  double sxx, sxy, sxz, syy, syz, szz;
            //

            // compute the Is function derivatives....

            double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;

            double Ip111, Ip122, Ip133, Ip112, Ip113, Ip123, Ip222, Ip223, Ip233, Ip333 = 0.;
            if (il::abs(z)<0.01)
            {
                if ((il::abs(x)/a)<10.){
                    double abs_y = il::abs(y);
                    Ip11 = rectangular_integration(x, abs_y, z, a, b, ip11);
                    Ip33 = rectangular_integration(x, abs_y, z, a, b, ip33);

                    if (y < 0.){Ip12 = - rectangular_integration(x, abs_y, z, a, b, ip12);} else { Ip12 = rectangular_integration(x, abs_y, z, a, b, ip12);}

                    Ip13 = rectangular_integration(x, abs_y, z, a, b, ip13);

                    Ip111 = rectangular_integration(x, abs_y, z, a, b, ip111);
                    Ip122 = rectangular_integration(x, abs_y, z, a, b, ip122);
                    Ip133 = rectangular_integration(x, abs_y, z, a, b, ip133);

                    if (y < 0.){Ip112 = -rectangular_integration(x, abs_y, z, a, b, ip112);} else { Ip112 = rectangular_integration(x, abs_y, z, a, b, ip112);}

                    Ip113 = rectangular_integration(x, abs_y, z, a, b, ip113);

                    if (y < 0.){Ip123 = -rectangular_integration(x, abs_y, z, a, b, ip123);} else { Ip123 = rectangular_integration(x, abs_y, z, a, b, ip123);}

                    // Ip333 <--- this is non trivial that it goes to 0 when z->0 i.e. limit to be taken
                } else{

                    Ip11 = rectangular_integration(x, y, z, a, b, ip11);
                    Ip33 = rectangular_integration(x, y, z, a, b, ip33);
                    Ip12 = rectangular_integration(x, y, z, a, b, ip12);
                    Ip13 = rectangular_integration(x, y, z, a, b, ip13);

                    Ip111 = rectangular_integration(x, y, z, a, b, ip111);
                    Ip122 = rectangular_integration(x, y, z, a, b, ip122);
                    Ip133 = rectangular_integration(x, y, z, a, b, ip133);

                    Ip112 = rectangular_integration(x, y, z, a, b, ip112);
                    Ip113 = rectangular_integration(x, y, z, a, b, ip113);
                    Ip123 = rectangular_integration(x, y, z, a, b, ip123);

                    // Ip333 <--- this is non trivial that it goes to 0 when z->0 i.e. limit to be taken
                }


                if ((il::abs(y)/b)<10.){
                    double abs_x = il::abs(x);
                    Ip22 = rectangular_integration(abs_x , y, z, a, b, ip22);
                    Ip23 = rectangular_integration(abs_x , y, z, a, b, ip23);
                    Ip222 = rectangular_integration(abs_x, y, z, a, b, ip222);
                    Ip233 = rectangular_integration(abs_x, y, z, a, b, ip233);
                    Ip223 = rectangular_integration(abs_x, y, z, a, b, ip223);
                }
                else{
                    Ip22 = rectangular_integration(x, y, z, a, b, ip22);
                    Ip23 = rectangular_integration(x, y, z, a, b, ip23);
                    Ip222 = rectangular_integration(x, y, z, a, b, ip222);
                    Ip233 = rectangular_integration(x, y, z, a, b, ip233);
                    Ip223 = rectangular_integration(x, y, z, a, b, ip223);
                }
            }
            else{
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
                // Ip333 <--- this is non trivial that it goes to 0 when z->0 i.e. limit to be taken
            }
//            Ip11 = rectangular_integration(x, y, z, a, b, ip11);
//            Ip22 = rectangular_integration(x, y, z, a, b, ip22);
//            Ip33 = rectangular_integration(x, y, z, a, b, ip33);
//            Ip23 = rectangular_integration(x, y, z, a, b, ip23);
//            Ip12 = rectangular_integration(x, y, z, a, b, ip12);
//            Ip13 = rectangular_integration(x, y, z, a, b, ip13);
//
//            Ip111 = rectangular_integration(x, y, z, a, b, ip111);
//            Ip122 = rectangular_integration(x, y, z, a, b, ip122);
//            Ip133 = rectangular_integration(x, y, z, a, b, ip133);
//            Ip112 = rectangular_integration(x, y, z, a, b, ip112);
//            Ip113 = rectangular_integration(x, y, z, a, b, ip113);
//            Ip123 = rectangular_integration(x, y, z, a, b, ip123);
//            Ip222 = rectangular_integration(x, y, z, a, b, ip222);
//            Ip233 = rectangular_integration(x, y, z, a, b, ip233);
//            Ip223 = rectangular_integration(x, y, z, a, b, ip223);
            // Ip333 <--- this is non trivial that it goes to 0 when z->0 i.e. limit to be taken

            double Ip33_lim_z_to_0, z_times_Ip333;

            if (il::abs(z) <= EPSILON )
                z_times_Ip333 = 0.;
            else{ Ip333 = rectangular_integration(x, y, z, a, b, ip333);
                  z_times_Ip333 = z * Ip333;}

            if (il::abs(il::abs(x) / a - 1.) <= EPSILON && (il::abs(y) > b) && il::abs(z) <= EPSILON)
                Ip33_lim_z_to_0 = Ip33_lim_z_to_0_and_x_to_a(x, y, a, b);
            else if (il::abs(il::abs(y) / b - 1.) <= EPSILON && (il::abs(x) > a) && il::abs(z) <= EPSILON)
                Ip33_lim_z_to_0 = Ip33_lim_z_to_0_and_y_to_b(x, y, a, b);
            else
                Ip33_lim_z_to_0 = Ip33;

            // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

            // stress due to displacement discontinuity DDx (shear)
            Stress(0, 0) = Ce * (2. * Ip13 - z * Ip111);                    // sxx -> if z=0. it will be 0. (unchanged expresssion)
            Stress(0, 1) = Ce * (2.   * nu * Ip13 - z * Ip122);             // syy -> if z=0. it will be 0. (unchanged expresssion)
            Stress(0, 2) = Ce * (-z * Ip133);                               // szz -> if z=0. it will be 0. (unchanged expresssion)
            Stress(0, 3) = Ce * ((1. - nu) * Ip23 - z * Ip112);             // sxy -> if z=0. it will be 0. (unchanged expresssion)
            Stress(0, 4) = Ce * (Ip33_lim_z_to_0 + nu * Ip22 - z * Ip113);  // sxz
            Stress(0, 5) = Ce * (-nu * Ip12 - z * Ip123);                   // syz if z=0. (unchanged expresssion)

            // stress due to displacement discontinuity  DDy (shear)

            Stress(1, 0) = Ce * (2. * nu * Ip23 - z * Ip112);               // sxx -> if z=0. it will be 0. (unchanged expresssion)
            Stress(1, 1) = Ce * (2. * Ip23 - z * Ip222);                    // syy -> if z=0. it will be 0. (unchanged expresssion)
            Stress(1, 2) = Ce * (-z * Ip233);                               // szz -> if z=0. it will be 0. (unchanged expresssion)
            Stress(1, 3) = Ce * ((1. - nu) * Ip13 - z * Ip122);             // sxy -> if z=0. it will be 0. (unchanged expresssion)
            Stress(1, 4) = Stress(0, 5);                             // sxz -> if z=0. (unchanged expresssion)
            Stress(1, 5) = Ce * (Ip33_lim_z_to_0 + nu * Ip11 - z * Ip223);  // syz

            // stress due to displacement discontinuity DDz (normal)
            Stress(2, 0) = Ce * (Ip33_lim_z_to_0 + (1 - 2 * nu) * Ip22 - z * Ip113); // sxx
            Stress(2, 1) = Ce * (Ip33_lim_z_to_0 + (1 - 2 * nu) * Ip11 - z * Ip223); // syy
            Stress(2, 2) = Ce * (Ip33_lim_z_to_0 - z_times_Ip333);                   // szz
            Stress(2, 3) = Ce * ((-1 + 2 * nu) * Ip12 - z * Ip123);                  // sxy if z=0. (unchanged expresssion)
            Stress(2, 4) = Ce * (-z * Ip133);                                        // sxz -> if z=0. it will be 0. (unchanged expresssion)
            Stress(2, 5) = Ce * (-z * Ip233);                                        // syz -> if z=0. it will be 0. (unchanged expresssion)

//            for (il::int_t i = 0; i < 3; i++) {
//                for (il::int_t j = 0; j < 6; j++) {
//                    if (std::isnan(Stress(i,j))){
//                        printf("found NAN");
//                    }
//                }
//            }

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


    // Fundamental displacement kernel
    il::Array2D<double> DisplacementKernelR0(
            double& x, double& y, double& z, double& a, double& b,
            double& nu) {
        //  x , y , z location where to compute displacement
        //  a,b  1/2 size of the rectangular DD
        //  nu Poisson's ratio'
        //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
        //  DDx (shear), DDy (shear), DDz (normal)

        double Ip1, Ip2, Ip3_out_plane_z_EQ_0;

        double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;

        double Ce = (-1. / (8. * il::pi * (1. - nu)));

        il::Array2D<double> Displacement{3,3,0.};

        // compute the Is function derivatives....
        Ip3_out_plane_z_EQ_0 = rectangular_integration(x, y, z, a, b, ip3);


        if (il::abs(z)<0.01)
        {
            if ((il::abs(x)/a)<10.){
                double abs_y = il::abs(y);
                Ip1 = rectangular_integration(x, abs_y, z, a, b, ip1);
                Ip11 = rectangular_integration(x, abs_y, z, a, b, ip11);
                Ip33 = rectangular_integration(x, abs_y, z, a, b, ip33);
                if (y < 0.){Ip12 = - rectangular_integration(x, abs_y, z, a, b, ip12);} else { Ip12 = rectangular_integration(x, abs_y, z, a, b, ip12);}
                Ip13 = rectangular_integration(x, abs_y, z, a, b, ip13);
            }else{
                Ip1 = rectangular_integration(x, y, z, a, b, ip1);
                Ip11 = rectangular_integration(x, y, z, a, b, ip11);
                Ip33 = rectangular_integration(x, y, z, a, b, ip33);
                Ip12 = rectangular_integration(x, y, z, a, b, ip12);
                Ip13 = rectangular_integration(x, y, z, a, b, ip13);
            }

            if ((il::abs(y)/b)<10.){
                double abs_x = il::abs(x);
                Ip2 = rectangular_integration(abs_x , y, z, a, b, ip2);
                Ip22 = rectangular_integration(abs_x , y, z, a, b, ip22);
                Ip23 = rectangular_integration(abs_x , y, z, a, b, ip23);
            }else{
                Ip2 = rectangular_integration(x, y, z, a, b, ip2);
                Ip22 = rectangular_integration(x, y, z, a, b, ip22);
                Ip23 = rectangular_integration(x, y, z, a, b, ip23);
            }

        }
        else{
            Ip1 = rectangular_integration(x, y, z, a, b, ip1);
            Ip2 = rectangular_integration(x, y, z, a, b, ip2);
            Ip11 = rectangular_integration(x, y, z, a, b, ip11);
            Ip22 = rectangular_integration(x, y, z, a, b, ip22);
            Ip33 = rectangular_integration(x, y, z, a, b, ip33);
            Ip23 = rectangular_integration(x, y, z, a, b, ip23);
            Ip12 = rectangular_integration(x, y, z, a, b, ip12);
            Ip13 = rectangular_integration(x, y, z, a, b, ip13);
        }

//        Ip1 = rectangular_integration(x, y, z, a, b, ip1);
//        Ip2 = rectangular_integration(x, y, z, a, b, ip2);
//        Ip11 = rectangular_integration(x, y, z, a, b, ip11);
//        Ip22 = rectangular_integration(x, y, z, a, b, ip22);
//        Ip33 = rectangular_integration(x, y, z, a, b, ip33);
//        Ip23 = rectangular_integration(x, y, z, a, b, ip23);
//        Ip12 = rectangular_integration(x, y, z, a, b, ip12);
//        Ip13 = rectangular_integration(x, y, z, a, b, ip13);

        double Ip3 = get_Ip3 ( x, y, z, a, b, Ip3_out_plane_z_EQ_0);

        // Displacement row is dof (DDx,DDy,DDx), columns are Ux,Uy,Uz in the local reference system

        // displacement due to displacement discontinuity DDx (shear)
        Displacement(0, 0) = Ce * (z * Ip11 - 2 * (1 - nu) * Ip3);  // Ux
        Displacement(1, 0) = Ce * (z * Ip12);                       // Uy -> if z=0. it will be 0. (unchanged expresssion)
        Displacement(2, 0) = Ce * (z * Ip13 - (1 - 2 * nu) * Ip1);  // Uz

        // displacement due to displacement discontinuity  DDy (shear)
        Displacement(0, 1) = Displacement(1, 0);             // Ux  -> if z=0. it will be 0. (unchanged expresssion)
        Displacement(1, 1) = Ce * (z * Ip22 - 2 * (1 - nu) * Ip3);  // Uy
        Displacement(2, 1) = Ce * (z * Ip23 - (1 - 2 * nu) * Ip2);  // Uz

        // displacement due to displacement discontinuity DDz (normal)
        Displacement(0, 2) = Ce * (z * Ip13 + (1 - 2 * nu) * Ip1);  // Ux
        Displacement(1, 2) = Ce * (z * Ip23 + (1 - 2 * nu) * Ip2);  // Uy
        Displacement(2, 2) = Ce * (z * Ip33 - 2 * (1 - nu) * Ip3);  // Uz

//        Uncomment to check nan
//        for (il::int_t i = 0; i < 3; i++) {
//            for (il::int_t j = 0; j < 3; j++) {
//                if (std::isnan(Displacement(i,j))){
//                    printf("found NAN");
//                }
//            }
//        }
        return Displacement; // expressed in the reference system of the DD element
        // index        ->    DDx (shear)    DDy (shear)     DDz (normal)
        //   0      -> |       Ux,            Ux,             Ux            |
        //   1      -> |       Uy,            Uy,             Uy            |
        //   2      -> |       Uz,            Uz,             Uz            |
    }
//
}
