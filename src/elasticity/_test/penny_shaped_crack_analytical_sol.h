//
// This file is part of BigWham.
//
// Created by Carlo Peruzzo on 01.01.21.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE file for more details.
//

#include <iostream>
#include <math.h>

namespace bigwham {

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
double Ux_nl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &pz);

double Uy_nl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &pz);

double Uz_nl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &pz);

// shear
double Ux_sl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &px, double &py);

double Uy_sl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &px, double &py);

double Uz_sl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &px, double &py);

/*
 *
 * STRESSES
 * - DUE TO AN UNIFORM SHEAR LOADING px and py OVER A CIRCULAR CRACK OF RADIUS a
 * - DUE TO AN UNIFORM SHEAR LOADING pz OVER A CIRCULAR CRACK OF RADIUS a
 * expressions valid for any z
 *
 */

double Sig_xx_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double Sig_yy_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double Sig_zz_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double Sig_xy_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double Sig_zy_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double Sig_zx_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double Sig_xx_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double Sig_yy_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double Sig_zz_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double Sig_xy_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double Sig_zy_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double Sig_zx_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

// ------------------------------------------------------------------------------------

// normal
double ux_nl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &pz);

double uy_nl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &pz);

double uz_nl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &pz);

// shear
double ux_sl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &px, double &py);

double uy_sl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &px, double &py);

double uz_sl_(double &x, double &y, double &z, double &a, double &G, double &nu,
              double &px, double &py);

double sig_xx_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double sig_yy_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double sig_zz_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double sig_xy_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double sig_zy_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double sig_zx_nl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &pz);

double sig_xx_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double sig_yy_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double sig_zz_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double sig_xy_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double sig_zy_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

double sig_zx_sl_(double &x, double &y, double &z, double &a, double &G,
                  double &nu, double &px, double &py);

} // namespace bigwham
