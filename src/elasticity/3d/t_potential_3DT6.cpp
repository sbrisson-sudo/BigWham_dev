//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 17.12.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2021.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
#include <elasticity/3d/constants.h>
#include <elasticity/3d/h_potential_3DT6.h>
#include <il/StaticArray.h>
#include <il/StaticArray3D.h>
#include <il/StaticArray4D.h>
#include <il/math.h>
#include <complex>

#include "t_potential_3DT6.h"

namespace bie {







// Limit case (h==0, plane) - all stress components
    il::StaticArray3D<std::complex<double>, 6, 4, 3> s_ij_lim_t(double nu, std::complex<double> eix,std::complex<double> d){
//  nu :: poisson's ratio
// iex :: position of source point in local coor
// d ::
//SnH
// hi
        double c_1_2nu = 1.0 + 2.0 * nu;
        double c_1_m2nu = 1.0 - 2.0 * nu;
        double c_2_mnu = 2.0 - nu;

        double cos_x = std::real(eix);
        double sin_x = std::imag(eix);
        // double tan_x = sin_x/cos_x;
        std::complex<double> e2x = eix * eix;
        double h0_lim = std::atanh(sin_x); // std::complex<double> ?
        // double h0_lim = 0.5*(std::log(1.0+sin_x)-std::log(1.0-sin_x));

        double d_1 = std::abs(d); // DM in matlab code
        // double d_2 = d_1*d_1; double d_4 = d_2*d_2;
        std::complex<double> e = std::polar(1.0, std::arg(d)), // = d/d_1
        e_2 = e * e,         e_3 = e * e_2,         e_4 = e_2 * e_2, p1, p2;
// du ==e
        std::complex<double> e_c = std::conj(e);

        il::StaticArray<std::complex<double>, 6> v1{}, v2{},v1c{},v2c{},v3{},v4{},v3c{},v5{},v6{};

        //T
        v1[0] = 2.*e * h0_lim;
        v1[1] =2.*e_2*d_1 * h0_lim;
        v1[2]= -2.*d_1 * h0_lim;
        v1[3]=  3.*e_3 * d_1 * d_1 * h0_lim;
        v1[4]= -e_c * d_1 * d_1  * h0_lim;
        v1[5]=    -e*(d_1 * d_1) * h0_lim;

        v1c[0]=std::conj(v1[0]);//2.*std::conj(e) * h0_lim;
        v1c[1]=-2.* d_1 * h0_lim;//v1[2]
        v1c[2]=2.*(e_c*e_c)* d_1 * h0_lim;
        v1c[3]=-e*(d_1 * d_1) * h0_lim;
        v1c[4]=3.*(e_c*e_c*e_c)*(d_1*d_1) * h0_lim;
        v1c[5]=std::conj(v1[5]);//-duc.*(d_1*d_1) * h0_lim;

        double ron = d_1 /cos_x;
        std::complex<double> PmnHi=d_1*(4.*eix + 1/cos_x);
        v2[0]=0;
        v2[1]=4.*e_2*d_1*eix;
        v2[2]=0.;
        v2[3]=(e_c*e_c)*PmnHi;
        v2[4]=e_c*d_1*ron;
        v2[5]=e*d_1*ron;

        std::complex<double> eix_c = std::conj(eix);
        v2c[0]=0.;
        v2c[1]=0;
        v2c[2]=4.*(e_c*e_c)*d_1*eix_c;
        v2c[3]=e*d_1*ron;
        v2c[4]=(e_c*e_c)*std::conj(PmnHi);
        v2c[5]=e_c*d_1*ron;

        v3[0]=5.*d_1*eix;
        v3[1]=24.*e_2*d_1*eix;
        v3[2]=0;
        v3[3]=6.*e_3*PmnHi;
        v3[4]=6.*e_c*d_1*ron;
        v3[5]=6.*e*d_1*ron;


        v3c[0]=5.*d_1*eix_c;
        v3c[1]=0;
        v3c[2]=24.*(e_c*e_c)*d_1*eix_c;
        v3c[3]=6.*e*d_1*ron;
        v3c[4]=6.*(e_c*e_c*e_c)*conj(PmnHi);
        v3c[5]=6.*e_c*d_1*ron;


        v4[0]=-2.*e_3;
        v4[1]=-6.*e_4*d_1;
        v4[2]=6.*d_1*e_2;
        v4[3]=-15.*e_3*(d_1 * d_1)*e_2;
        v4[4]=-3.*e_c*(d_1 *d_1)*e_2;
        v4[5]=9.*e*(d_1*d_1) *e_2;

        //E2Hi=EiHi.^2;
        v5[0]=e_3*(4.+e2x/2.);
        v5[0]=-4.*e_4*d_1*(3.-e2x);
        v5[0]=12.*e_2*d_1*eix;
        v5[0]=-e_4*e*ron*(15.+10.*e2x-2.*e2x*e2x);
        v5[0]=-3.*e*ron;
        v5[0]=3.*e_3*d_1*PmnHi;

        v6[0]=2.5*ii*e*eix;
        v6[1]=0;
        v6[2]=0;
        v6[3]=0;
        v6[4]=0;
        v6[5]=0;

        il::StaticArray3D<std::complex<double>, 6, 4, 3> c_array{0.0};

        int SnH =-1; // approaching from the opposite of the normal vector... to pass as arg
        int LS =-1;
        if (SnH==1) {
            LS=1;
        }

        for (int j = 0; j < c_array.size(0); ++j) {
            c_array(j, 0, 0)=(1./16.)*(6.0*v1[j]+ii*v3[j]);
            c_array(j, 0, 1)=(1./16.)*(6.0*v1c[j]-ii*v3c[j]);

            c_array(j, 1, 0) =(1./8.)*(7.0-8*nu)*(v1[j]+ii*v2[j]+v6[j]);
            c_array(j, 1, 1) =(1./8.)*(v4[j]+ii*v5[j]);

            c_array(j, 2, 2) = 0.25*(c_1_m2nu)*(v1[j]+ii*v2[j]);

            c_array(j, 3, 0) = -0.25*c_1_m2nu*(v1[j]+ii*v2[j]);
            c_array(j, 3, 1) = -0.25*c_1_m2nu*(v1c[j]-ii*v2c[j]);
        }

        c_array(0, 0, 2)= c_array(0, 0, 2)-LS*nu*Hi; //Hi ... h
        c_array(2, 0, 0)= c_array(1, 0, 0)+LS*(1-nu)*Hi/4.;// Hi...h to pass
        c_array(3, 0, 0)= c_array(3, 0, 2)+LS*(1-nu)*Hi/2.;

        return c_array;
// % divide everything by 4.*pi.*(1-Nu) if necessary
    }


}