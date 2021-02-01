//
// This file is part of HFPx3D.
//
// Created by Carlo Peruzzo on 01.01.21.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <iostream>
#include <cmath.h>


namespace bie{
    double pi_() {return 3.14159265358979323846}

    double r_(double &x, double &y){return sqrt(x * x + y * y)}

    double theta_(double &x, double &y){return atan(y / x)}

    double f1_(double &x, double &y, double &z, double &a){
        double zz = z * z;
        double f1 = sqrt((a + r) * (a + r) + zz);
        return f1
    }

    double f2_(double &x, double &y, double &z, double &a){
        double zz = z * z;
        double f2 = sqrt((a - r) * (a - r) + zz);
        return f2
    }

    double l1_(double &x, double &y, double &z, double &a){
        double f1 = f1_(double &x, double &y, double &z, double &a);
        double f2 = f2_(double &x, double &y, double &z, double &a);
        return 0.5 * (f1 - f2)
    }

    double l2_(double &x, double &y, double &z, double &a){
        double f1 = f1_(double &x, double &y, double &z, double &a);
        double f2 = f2_(double &x, double &y, double &z, double &a);
        return 0.5 * (f1 + f2)
    }

    double l1_l1_(double &x, double &y, double &z, double &a){
        double l1 = l1_(double &x, double &y, double &z, double &a)
        return l1 * l1
    }

    double l2_l2_(double &x, double &y, double &z, double &a){
        double l2 = l2_(double &x, double &y, double &z, double &a)
        return l2 * l2
    }

    double one_m_2nu( double& nu){return (1 - 2. * nu)}

    double f3_(double &x, double &y, double &z, double &a){
        double l2l2 = l2_l2_(double &x, double &y, double &z, double &a);
        return sqrt(l2l2 - a*a)
    }

    double f4_(double &x, double &y, double &z, double &a){
        double l1l1 = l1_l1_(double &x, double &y, double &z, double &a);
        return sqrt(a*a - l1l1)
    }

    double f5_(double &x, double &y, double &z, double &a){
        double l2l2 = l2_l2_(double &x, double &y, double &z, double &a);
        double f3 = f3_(double &x, double &y, double &z, double &a);
        return a * f3 / l2l2
    }

    double ux_ux_nl_base_function(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double l1 = l1_(double &x, double &y, double &z, double &a);
        double l2 = l2_(double &x, double &y, double &z, double &a);
        double l1l1 = l1_l1_(double &x, double &y, double &z, double &a);
        double l2l2 = l2_l2_(double &x, double &y, double &z, double &a);
        const double pi = pi_();
        double c1 = pz / (2. * pi * G);
        double c2 = one_m_2nu(nu);
        double c3 = f5_(x, y, z, a);
        double c4 = 2. * a * a * abs(z) * f4_(x, y, z, a);
        double c5 = l2l2 * (l2l2 - l1l1);
        return c1 * ( c2 * (c3 - asin(a / l2)) + c4 / c5)
    }

    double sign_(double &x){
        if (x>=0){ return 1}
        else { return -1}
    }

    double f6_(double &x, double &y, double &z, double &a, double &nu){
        double l2 = l2_(double &x, double &y, double &z, double &a);
        double c1 = 4 * nu - 5;
        double c2 = 4 * (1 - nu);
        return c1 * z * asin(a / l2) + c2 * sign_(z) * f4_(x, y, z, a)
    }

    double f7_(double &x, double &y, double &z, double &a){
        double l2l2 = l2_l2_(double &x, double &y, double &z, double &a);
        double l1l1 = l1_l1_(double &x, double &y, double &z, double &a);

        return a * z * f3_(x, y, z, a) / (l2l2 -l1l1)
    }

    /*
     *
     * DISPLACEMENTS DUE TO AN UNIFORM NORMAL LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for z >=0
     *
     */

    double ux_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return x * ux_ux_nl_base_function(x, y, z, a, G, nu, pz)
    }

    double uy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return y * ux_ux_nl_base_function(x, y, z, a, G, nu, pz)
    }

    double uz_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double l1 = l1_(double &x, double &y, double &z, double &a);
        double l2 = l2_(double &x, double &y, double &z, double &a);
        double l1l1 = l1_l1_(double &x, double &y, double &z, double &a);
        double l2l2 = l2_l2_(double &x, double &y, double &z, double &a);
        const double pi = pi_();
        double c1 = pz / (pi * G);
        double c2 = 2. * (1 - nu);
        double f4 = f4_(x, y, z, a)
        double c3 = sign_(z) * f4;
        double c4 = asin(a / l2);
        double c5 = (a * f4) / (l2l2 - l1l1)

        return c1 * ( c2 * (c3 - z * c4) + z * (c4 - c5))
    }

    /*
     *
     * DISPLACEMENTS DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for z >=0
     *
     */
    double ux_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        const double pi = pi_();
        double c1 = 1. / (pi * G * (2. - nu));
        double f6 = f6_(x, y, z, a, nu);
        double f7 = f7_(x, y, z, a);
        double theta = theta_(x, y);
        double l1l1 = l1_l1_(double &x, double &y, double &z, double &a);
        double l2l2 = l2_l2_(double &x, double &y, double &z, double &a);
        double c2 = (px + (px * sin(2. * theta) + py * cos(2. * theta)) * l1l1 / l2l2);

        return c1 * (f6 * px + f7 * c2)
    }

    double uy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        const double pi = pi_();
        double c1 = 1. / (pi * G * (2. - nu));
        double f6 = f6_(x, y, z, a, nu);
        double f7 = f7_(x, y, z, a);
        double theta = theta_(x, y);
        double l1l1 = l1_l1_(double &x, double &y, double &z, double &a);
        double l2l2 = l2_l2_(double &x, double &y, double &z, double &a);
        double c2 = (py + (px * sin(2. * theta) - py * cos(2. * theta)) * l1l1 / l2l2);

        return c1 * (f6 * py + f7 * c2)
    }

    double uz_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        const double pi = pi_();
        double c1 = (2. * x * px + 2. * y * py) / (pi * G * (2. - nu));
        double c2 = (1 - 2. * nu) / 2.;
        double l2 = l2_(double &x, double &y, double &z, double &a);
        double l2l2 = l2 * l2;
        double f3 = f3_(x, y, z, a);
        double c3 = (asin(a / l2) - a * f3 / l2l2);
        double f4 = f4_(x, y, z, a);
        double l1l1 = l1_l1_(double &x, double &y, double &z, double &a);
        double c4 = z * a * a * f4 / (l2l2 * (l2l2 - l1l1))

        return c1 * (c2 * c3 + c4 )
    }

    std::Array2D<double> symmetry_(double& x, double& y, double& z){
        std::Array2D<double> coord(3,0.);
        if (z < 0.){
            coord[0] = -x; coord[1] = -y; coord[2] = +z;
        }
        else {
            coord[0] = +x; coord[1] = +y; coord[2] = +z;
        }
        return coord
    }

    typedef double (*vFunctionCall_switch_z)double &x_, double &y_, double &z_, double &a_, double&G_, double&nu_, double&pz_);
    typedef double (*vFunctionCall_switch_xy)double &x_, double &y_, double &z_, double &a_, double&G_, double&nu_, double&px_, double&py_);

    double Func_call(double& x, double& y, double& z, double& a, double& G, double& nu, double& px, double& py, vFunctionCall_switch_xy Func) {
        std::Array2D<double> coord = symmetry_(x, y, z)
        return Func(coord[0], coord[1], coord[2], a, G, nu, px, py);
    }

    double Func_call(double& x, double& y, double& z, double& a, double& G, double& nu, double& pz, vFunctionCall_switch_z Func) {
        std::Array2D<double> coord = symmetry_(x, y, z)
        return Func(coord[0], coord[1], coord[2], a, G, nu, pz);
    }
    // ------------------------------------------------------------------------------------


    /*
     *
     * STRESSES DUE TO AN UNIFORM NORMAL LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for z >=0
     *
     */
    double sig_1_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_2R_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_xx_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return 0.5 * (sig_1_nl_(x, y, z, a, G, nu, pz) + sig_2R_nl_(x, y, z, a, G, nu, pz))
    }
    double sig_yy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return 0.5 * (sig_1_nl_(x, y, z, a, G, nu, pz) - sig_2R_nl_(x, y, z, a, G, nu, pz))
    }
    double sig_xy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_zy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_zx_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_zz_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}

    /*
     *
     * STRESSES DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for z >=0
     *
     */
    double sig_1_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_2R_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_xx_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_yy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_xy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_zy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_zx_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}
    double sig_zz_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){}

    // ------------------------------------------------------------------------------------
    /*
     *
     * DISPLACEMENTS
     * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
     * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for any z
     *
     */

    // normal
    double Ux_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return symmetry_switch_z(x, y, z, a, G, nu, pz, ux_nl_)
    }

    double Uy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return symmetry_switch_z(x, y, z, a, G, nu, pz, uy_nl_)
    }

    double Uz_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return symmetry_switch_z(x, y, z, a, G, nu, pz, uz_nl_)
    }

    // shear
    double Ux_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        return symmetry_switch_xy(x, y, z, a, G, nu, px, py, ux_sl_)
    }

    double Uy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        return symmetry_switch_xy(x, y, z, a, G, nu, px, py, uy_sl_)
    }

    double Uz_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        return symmetry_switch_xy(x, y, z, a, G, nu, px, py, uz_sl_)
    }

    /*
     *
     * STRESSES
     * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
     * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for any z
     *
     */
}