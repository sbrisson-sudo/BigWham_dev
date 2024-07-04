//
// This file is part of HFPx3D.
//
// Created by Carlo Peruzzo on 01.01.21.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2021.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <vector>

#include "penny_shaped_crack_analytical_sol.h"

namespace bigwham{
    double pi_() {return 3.14159265358979323846;}

    double r_(double &x, double &y){return sqrt(x * x + y * y);}

    double theta_(double &x, double &y){return atan(y / x);}

    double f1_(double &x, double &y, double &z, double &a){
        double r = r_(x, y);
        double zz = z * z;
        double f1 = sqrt((a + r) * (a + r) + zz);
        return f1;
    }

    double f2_(double &x, double &y, double &z, double &a){
        double r = r_(x, y);
        double zz = z * z;
        double f2 = sqrt((a - r) * (a - r) + zz);
        return f2;
    }

    double l1_(double &x, double &y, double &z, double &a){
        double f1 = f1_(x, y, z, a);
        double f2 = f2_(x, y, z, a);
        return 0.5 * (f1 - f2);
    }

    double l2_(double &x, double &y, double &z, double &a){
        double f1 = f1_(x, y, z, a);
        double f2 = f2_(x, y, z, a);
        return 0.5 * (f1 + f2);
    }

    double l1_l1_(double &x, double &y, double &z, double &a){
        double l1 = l1_(x, y, z, a);
        return l1 * l1;
    }

    double l2_l2_(double &x, double &y, double &z, double &a){
        double l2 = l2_(x, y, z, a);
        return l2 * l2;
    }

    double one_m_2nu( double& nu){return (1 - 2. * nu);}

    double f3_(double &x, double &y, double &z, double &a){
        double l2l2 = l2_l2_(x, y, z, a);
        return sqrt(l2l2 - a*a);
    }

    double f4_(double &x, double &y, double &z, double &a){
        double l1l1 = l1_l1_(x, y, z, a);
        return sqrt(a*a - l1l1);
    }

    double f5_(double &x, double &y, double &z, double &a){
        double l2l2 = l2_l2_(x, y, z, a);
        double f3 = f3_(x, y, z, a);
        return a * f3 / l2l2;
    }

    double ux_ux_nl_base_function(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double l2 = l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double l2l2 = l2_l2_(x, y, z, a);
        const double pi = pi_();
        double c1 = pz / (2. * pi * G);
        double c2 = one_m_2nu(nu);
        double c3 = f5_(x, y, z, a);
        double c4 = 2. * a * a * abs(z) * f4_(x, y, z, a);
        double c5 = l2l2 * (l2l2 - l1l1);
        return c1 * ( c2 * (c3 - asin(a / l2)) + c4 / c5);
    }

    double sign_(double &x){
        if (x>=0){ return 1;}
        else { return -1;}
    }

    double f6_(double &x, double &y, double &z, double &a, double &nu){
        double l2 = l2_(x, y, z, a);
        double c1 = 4 * nu - 5;
        double c2 = 4 * (1 - nu);
        return c1 * z * asin(a / l2) + c2 * sign_(z) * f4_(x, y, z, a);
    }

    double f7_(double &x, double &y, double &z, double &a){
        double l2l2 = l2_l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);

        return a * z * f3_(x, y, z, a) / (l2l2 -l1l1);
    }

    /*
     *
     * DISPLACEMENTS DUE TO AN UNIFORM NORMAL LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for z >=0
     *
     */

    double ux_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return x * ux_ux_nl_base_function(x, y, z, a, G, nu, pz);
    }

    double uy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return y * ux_ux_nl_base_function(x, y, z, a, G, nu, pz);
    }

    double uz_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double l2 = l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double l2l2 = l2_l2_(x, y, z, a);
        const double pi = pi_();
        double c1 = pz / (pi * G);
        double c2 = 2. * (1 - nu);
        double f3 = f3_(x, y, z, a); //sqrt(l2l2 - a*a)
        double f4 = f4_(x, y, z, a); //sqrt(a*a - l1l1)
        double c3 = sign_(z) * f4;
        double c4 = asin(a / l2);
        double c5 = (a * f3) / (l2l2 - l1l1);

        return c1 * ( c2 * (c3 - z * c4) + z * (c4 - c5));
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
        double l1l1 = l1_l1_(x, y, z, a);
        double l2l2 = l2_l2_(x, y, z, a);
        double c2 = (px + (px * cos(2. * theta) + py * sin(2. * theta)) * l1l1 / l2l2);

        return c1 * (f6 * px + f7 * c2);
    }

    double uy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        const double pi = pi_();
        double c1 = 1. / (pi * G * (2. - nu));
        double f6 = f6_(x, y, z, a, nu);
        double f7 = f7_(x, y, z, a);
        double theta = theta_(x, y);
        double l1l1 = l1_l1_(x, y, z, a);
        double l2l2 = l2_l2_(x, y, z, a);
        double c2 = (py + (px * sin(2. * theta) - py * cos(2. * theta)) * l1l1 / l2l2);

        return c1 * (f6 * py + f7 * c2);
    }

    double uz_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        const double pi = pi_();
        double c1 = (2. * x * px + 2. * y * py) / (pi * G * (2. - nu));
        double c2 = (1 - 2. * nu) / 2.;
        double l2 = l2_(x, y, z, a);
        double l2l2 = l2 * l2;
        double f3 = f3_(x, y, z, a);
        double c3 = (asin(a / l2) - a * f3 / l2l2);
        double f4 = f4_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double c4 = z * a * a * f4 / (l2l2 * (l2l2 - l1l1));

        return c1 * (c2 * c3 + c4 );
    }


    /*
     *
     * STRESSES DUE TO AN UNIFORM NORMAL LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for z >=0
     *
     */

    double f8_(double &x, double &y, double &z, double &a){
        double f3 = f3_(x, y, z, a); // sqrt(l2l2 - a*a)
        double l2l2 = l2_l2_(x, y, z, a);
        double l2 = l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        return a * f3 / (l2l2 - l1l1) - asin(a / l2);
    }

    double common_nl_sigxy_sig2R(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        const double pi = pi_();
        double c1 =  pz / pi;

        double l1l1 = l1_l1_(x, y, z, a);
        double f3 = f3_(x, y, z, a); // sqrt(l2l2 - a*a)
        double c2 = a * l1l1 * f3;

        double l2l2 = l2_l2_(x, y, z, a);
        double c3 = l2l2 * (l2l2 - l1l1);
        double r = r_(x,y);
        double c4 = z * z * (a * a * (6. * l2l2 - 2. * l1l1 + r * r ) - 5. * l2l2 * l2l2);
        double c5 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - a * a);

        return c1 * c2 / c3 * (1. - 2. * nu + c4 / c5);
    }

    double common_nl_sigzx_sigzy(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        const double pi = pi_();
        double c1 = - 2. * pz / pi;
        double l1l1 = l1_l1_(x, y, z, a);
        double l1 = l1_(x, y, z, a);
        double l2 = l2_(x, y, z, a);
        double f3 = f3_(x, y, z, a); // sqrt(l2l2 - a*a)
        double l2l2 = l2_l2_(x, y, z, a);
        double r = r_(x, y);
        double c2 = z * l1 * f3 * (a * a * (4. * l2l2 - 5. * r * r) + l1l1 * l1l1);

        double c3 = l2 * (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1);
        return c1 * c2 / c3;
    }

    double sig_1_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        const double pi = pi_();
        double c1 = 2. * pz / pi;
        double c2 = 1 + 2 * nu;
        double c3 = f8_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double l2l2 = l2_l2_(x, y, z, a);
        double f3 = f3_(x, y, z, a); // sqrt(l2l2 - a*a)
        double r = r_(x, y);
        double aa = a * a, zz = z * z, rr = r * r;
        double c4 = a * zz * (l1l1 * l1l1 + aa * (2. * aa + 2. * zz - 3. * rr));
        double c5 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1) * f3;

        return c1 * (c2 * c3 + c4 / c5);
    }
    double sig_2R_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double theta = theta_(x,y);
        return 2 * cos(2 * theta) * common_nl_sigxy_sig2R(x, y, z, a, G, nu, pz);
    }
    double sig_xx_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return 0.5 * (sig_1_nl_(x, y, z, a, G, nu, pz) + sig_2R_nl_(x, y, z, a, G, nu, pz));
    }
    double sig_yy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return 0.5 * (sig_1_nl_(x, y, z, a, G, nu, pz) - sig_2R_nl_(x, y, z, a, G, nu, pz));
    }
    double sig_xy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double theta = theta_(x,y);
        return sin(2 * theta) * common_nl_sigxy_sig2R(x, y, z, a, G, nu, pz);
    }
    double sig_zy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double theta = theta_(x,y);
        return sin(theta) * common_nl_sigzx_sigzy(x, y, z, a, G, nu, pz);
    }
    double sig_zx_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double theta = theta_(x,y);
        return cos(theta) * common_nl_sigzx_sigzy(x, y, z, a, G, nu, pz);
    }
    double sig_zz_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        const double pi = pi_();
        double c1 = 2. * pz / pi;
        double c3 = f8_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double l2l2 = l2_l2_(x, y, z, a);
        double f3 = f3_(x, y, z, a); // sqrt(l2l2 - a*a)
        double r = r_(x, y);
        double aa = a * a, zz = z * z, rr = r * r;
        double c4 = a * zz * (l1l1 * l1l1 + aa * (2. * aa + 2. * zz - 3. * rr));
        double c5 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1) * f3;

        return c1 * ( c3 - c4 / c5);
    }

    /*
     *
     * STRESSES DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for z >=0
     *
     */

    double f9_(double &x, double &y, double &z, double &a){
        double l1 = l1_(x, y, z, a);
        double r = r_(x,y);
        double f3 = f3_(x, y, z, a); // sqrt(l2l2 - a*a)
        double l2l2 = l2_l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double l2 = l2_(x, y, z, a);

        double c1 = z * l1 * f3 * (a * a * (4. * l2l2 - 5. * r * r) + l1l1 * l1l1);
        double c2 = l2 * (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1);
        return c1 / c2;
    }

    double f10_(double &x, double &y, double &z, double &a){
        double r = r_(x,y);
        double f4 = f4_(x, y, z, a); // sqrt(a*a - l1l1)
        double l2l2 = l2_l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double aa = a * a;
        double c1 = z * f4 * (l1l1 * l1l1 + aa * (2. * aa + 2. * z * z - 3. * r * r));
        double c2 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1);
        return c1 / c2;
    }

    double f11_(double &x, double &y, double &z, double &a, double &nu){
        double r = r_(x,y);
        double aa = a * a, rr = r * r;
        double l2l2 = l2_l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double f3 = f3_(x, y, z, a); // sqrt(l2l2 - a*a)
        double f4 = f4_(x, y, z, a); // sqrt(a*a - l1l1)
        double c1 = nu * a * f3;
        double c2 = z * f4 * (aa * (6. * l2l2 - 2 * l1l1 + rr) - 5. * l2l2 * l2l2);
        double c3 = (l2l2 - l1l1) * (l2l2 - l1l1);
        return c1 + c2 / c3;
    }

    double f12_(double &x, double &y, double &z, double &a){
        double l2 = l2_(x, y, z, a);
        double l1 = l1_(x, y, z, a);
        double l2l2 = l2_l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double f4 = f4_(x, y, z, a); // sqrt(a*a - l1l1)
        return  a * l1 * f4 / (l2 * (l2l2 - l1l1));
    }

    double f13_(double &x, double &y, double &z, double &a){
        double r = r_(x,y);
        double f3 = f3_(x, y, z, a); // sqrt(l2l2 - a*a)
        double l2l2 = l2_l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double l2 = l2_(x, y, z, a);
        double l1 = l1_(x, y, z, a);

        double c1 = a * a * (4. * l2l2 - 5. * r * r) + l1l1 * l1l1;
        double c2 = (l2l2 - l1l1) * (l2l2 - l1l1) * (l2l2 - l1l1);
        return z * l1 * f3 * c1 / (l2 * c2);
    }

    double f14_(double &x, double &y, double &z, double &a){
        double l2l2 = l2_l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        double l1 = l1_(x, y, z, a);
        double l2 = l2_(x, y, z, a);
        double f3 = f3_(x, y, z, a); // sqrt(l2l2 - a*a)
        return 4. * z * a * a * l1 * f3 / (l2l2 * l2 * (l2l2 - l1l1));
    }

    double sig_1_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        const double pi = pi_();
        double theta = theta_(x, y);
        double c1 = 4. * (px * cos(theta) + py * sin(theta)) / (pi * (2. - nu));
        double c2 = f12_(x, y, z, a);
        double c3 = f13_(x, y, z, a);
        return c1 * (-2. * (1. + nu) * c2 + c3);
    }
    double sig_2R_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        const double pi = pi_();
        double theta = theta_(x, y);
        double c1 = -2./(pi * (2. - nu));
        double c2 = 4. * (1. - nu) * f12_(x, y, z, a);
        double c3 = f13_(x, y, z, a);
        double c4 = px * cos(theta) - py * sin(theta);
        double c5 = f14_(x, y, z, a);
        double c6 = px * cos(3. * theta) + py * sin(3. * theta);

        return c1 * ((c2 - c3) * c4 + (c5 - c3) * c6);
    }
    double sig_xx_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        return 0.5 * (sig_1_sl_(x, y, z, a, G, nu, px, py) + sig_2R_sl_(x, y, z, a, G, nu, px, py));
    }
    double sig_yy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        return 0.5 * (sig_1_sl_(x, y, z, a, G, nu, px, py) - sig_2R_sl_(x, y, z, a, G, nu, px, py));
    }
    double sig_xy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        const double pi = pi_();
        double theta = theta_(x, y);
        double c1 = -1./(pi * (2. - nu));
        double c2 = 4. * (1. - nu) * f12_(x, y, z, a);
        double c3 = f13_(x, y, z, a);
        double c4 = px * sin(theta) + py * cos(theta);
        double c5 = f14_(x, y, z, a);
        double c6 = px * sin(3. * theta) - py * cos(3. * theta);

        return c1 * ((c2 - c3) * c4 + (c5 - c3) * c6);
    }

    double sig_zy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double c1 = 2. / (pi_() * (2. - nu));
        double theta = theta_(x, y);
        double c2 = f8_(x, y, z, a);
        double c3 = f10_(x, y, z, a);
        double c4 = f11_(x, y, z, a, nu);
        double l2l2 = l2_l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        return c1 * ( ((2. - nu) * c2 + c3) * py + c4 * l1l1 * (px * sin(2. * theta) - py * cos(2. * theta)) / (l2l2 * (l2l2 - l1l1)) );
    }

    double sig_zx_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double c1 = 2. / (pi_() * (2. - nu));
        double theta = theta_(x, y);
        double c2 = f8_(x, y, z, a);
        double c3 = f10_(x, y, z, a);
        double c4 = f11_(x, y, z, a, nu);
        double l2l2 = l2_l2_(x, y, z, a);
        double l1l1 = l1_l1_(x, y, z, a);
        return c1 * ( ((2. - nu) * c2 + c3) * px + c4 * l1l1 * (px * cos(2. * theta) + py * sin(2. * theta)) / (l2l2 * (l2l2 - l1l1)) );
    }

    double sig_zz_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double theta = theta_(x, y);
        double c1 = -4. * (px * cos(theta) + py * sin(theta)) / (pi_() * (2. - nu));
        return c1 * f9_(x, y, z, a);
    }

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
        return ux_nl_(x, y, z, a, G, nu, pz);
    }

    double Uy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return uy_nl_(x, y, z, a, G, nu, pz);
    }

    double Uz_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return uz_nl_(x, y, z, a, G, nu, pz);
    }

    // shear
    double Ux_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double value, abs_z;
        abs_z = abs(z);
        value = ux_sl_(x, y, abs_z, a, G, nu, px, py);
        if (z<0){
            value = - value;
        }
        return value;
    }

    double Uy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double value, abs_z;
        abs_z = abs(z);
        value = uy_sl_(x, y, abs_z, a, G, nu, px, py);
        if (z<0){
            value = - value;
        }
        return value;
    }

    double Uz_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double value, abs_z;
        abs_z = abs(z);
        value = uz_sl_(x, y, abs_z, a, G, nu, px, py);
        return value;
    }

    /*
     *
     * STRESSES
     * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
     * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
     * expressions valid for any z
     *
     */
// Normal
    double Sig_xx_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return sig_xx_nl_(x, y, z, a, G, nu, pz);
    }

    double Sig_yy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return sig_yy_nl_(x, y, z, a, G, nu, pz);
    }

    double Sig_zz_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return sig_zz_nl_(x, y, z, a, G, nu, pz);
    }

    double Sig_xy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        return sig_xy_nl_(x, y, z, a, G, nu, pz);
    }

    double Sig_zy_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double value;
        value = sig_zy_nl_(x, y, z, a, G, nu, pz);
        if (z<0 && !((x>0 && y>0) || (x>0 && y<0))){
            value = - value;
        }
        return value;
    }

    double Sig_zx_nl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& pz){
        double value;
        value = sig_zx_nl_(x, y, z, a, G, nu, pz);
        if (z<0 && !((x>0 && y>0) || (x>0 && y<0))){
            value = - value;
        }
        return value;
    }

// Shear
    double Sig_xx_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double value, abs_z;
        abs_z = abs(z);
        value = sig_xx_sl_(x, y, abs_z, a, G, nu, px, py);
        if (z<0 && ((x>0 &&y>0)||(x>0 && y<0))){
            value = - value;
        }
        return value;
    }

    double Sig_yy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double value, abs_z;
        abs_z = abs(z);
        value = sig_yy_sl_(x, y, abs_z, a, G, nu, px, py);
        if (z<0 && ((x>0 &&y>0)||(x>0 && y<0))){
            value = - value;
        }
        return value;
    }

    double Sig_zz_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double value, abs_z;
        abs_z = abs(z);
        value = sig_zz_sl_(x, y, abs_z, a, G, nu, px, py);
        if (z<0 && ((x>0 &&y>0)||(x>0 && y<0))){
            value = - value;
        }
        return value;
    }

    double Sig_xy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double value, abs_z;
        abs_z = abs(z);
        value = sig_xy_sl_(x, y, abs_z, a, G, nu, px, py);
        if (z<0 && ((x>0 &&y>0)||(x>0 && y<0))){
            value = - value;
        }
        return value;
    }

    double Sig_zy_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double value, abs_z;
        abs_z = abs(z);
        value = sig_zy_sl_(x, y, abs_z, a, G, nu, px, py);
        return value;
    }

    double Sig_zx_sl_(double &x, double &y, double &z, double &a, double& G, double& nu, double& px, double& py){
        double value, abs_z;
        abs_z = abs(z);
        value = sig_zx_sl_(x, y, abs_z, a, G, nu, px, py);
        return value;
    }
}
