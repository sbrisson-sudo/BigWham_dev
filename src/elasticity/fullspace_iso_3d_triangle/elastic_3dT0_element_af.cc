//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 03.02.2024.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne),
// Switzerland, Geo-Energy Laboratory, 2016-2025.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// implementation of  Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, artefact-free solution.
//  Geophysical Journal International
// with modifications (notably moving to the dislocation convention of positive DD in overlap)

#include <cmath>
#include <iostream>

#include <il/core.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linearAlgebra.h>
#include <il/linearAlgebra/dense/norm.h>

#include "elements/triangle.h"
#include "elastic_3dT0_element_af.h"

namespace bigwham{

    il::StaticArray2D<double,3,3> AngDisDisp(double x, double y, double z, double alpha, double nu) {
        // Displacement from an angular dislocation in an elastic isotropic full space
        //calculates the "incomplete" displacements (without the Burgers' function contribution) associated with an angular dislocation in
        // an elastic isotropic full-space.
        // the ADCS: Angular Dislocation Coordinate System is defined such that x is the normal ! and
        // the dislocation is in plane x=0 with one edge along z and the other at angle alpha from e_z
        // see Nikkko & Walker 2015, fig 1.
        //
        // note  - convention of positive DDs in opening.
        // Calculate intermediate values
        double cosA = std::cos(alpha);
        double sinA = std::sin(alpha);
        double eta = y * cosA - z * sinA;
        double zeta = y * sinA + z * cosA;
        double r = std::sqrt(x * x + y * y + z * z);
        double pre_factor = 1.0 / (8 * il::pi * (1 - nu));

        // Avoid complex results for logarithmic terms
        if (zeta >= r) {
            zeta = r-std::numeric_limits<double>::epsilon();
        }
        if (z >= r) {
            z = r-std::numeric_limits<double>::epsilon();
         }

        // Calculate displacement components
        double ux = pre_factor * (x * y / r / (r - z) - x * eta / r / (r - zeta));
        double vx = pre_factor * (eta * sinA / (r - zeta) - y * eta / r / (r - zeta) +
                                                  y * y / r / (r - z) + (1 - 2 * nu) * (cosA * std::log(r - zeta) - std::log(r - z)));
        double wx = pre_factor * (eta * cosA / (r - zeta) - y / r - eta * z / r / (r - zeta) -
                                                  (1 - 2 * nu) * sinA * std::log(r - zeta));

        double uy = pre_factor * (x * x * cosA / r / (r - zeta) - x * x / r / (r - z) - (1 - 2 * nu) * (cosA * std::log(r - zeta) - std::log(r - z)));
        double vy = pre_factor * x * (y * cosA / r / (r - zeta) - sinA * cosA / (r - zeta) - y / r / (r - z));
        double wy = pre_factor * x * (z * cosA / r / (r - zeta) - cosA * cosA / (r - zeta) + 1 / r);

        double uz = pre_factor * sinA  * ((1 - 2 * nu) * std::log(r - zeta) - x * x / r / (r - zeta));
        double vz =pre_factor * x * sinA   * (sinA / (r - zeta) - y / r / (r - zeta));
        double wz = pre_factor * x * sinA   * (cosA / (r - zeta) - z / r / (r - zeta));

        // We return the displacement influence tensor - in the Angular Dislocation Coordinate system.
        il::StaticArray2D<double,3,3> U_angular_dislocation{0.};
// index        ->    bx, by, bz
        //   0      -> |       Ux,            Ux,             Ux            |
        //   1      -> |       Uy,            Uy,             Uy            |
        //   2      -> |       Uz,            Uz,             Uz            |

        U_angular_dislocation(0,0)=ux;
        U_angular_dislocation(1,0)=vx;
        U_angular_dislocation(2,0)=wx;

        U_angular_dislocation(0,1)=uy;
        U_angular_dislocation(1,1)=vy;
        U_angular_dislocation(2,1)=wy;

        U_angular_dislocation(0,2)=uz;
        U_angular_dislocation(1,2)=vz;
        U_angular_dislocation(2,2)=wz;
//        for (il::int_t i=0;i<3;i++){
//            std::cout << " - " << U_angular_dislocation(0,i) << " - " << U_angular_dislocation(1,i)<< " - " << U_angular_dislocation(2,i) <<"\n";
//        }

        return U_angular_dislocation;
    }


    // Function to find trimode
    int trimodefinder(const il::StaticArray<double,3> x_obs, const il::Array<double> p1,const il::Array<double> p2,
                      const il::Array<double> p3) {
//       x_obs   x, y and z share the same size and
//       correspond to the y, z and x coordinates in the TDCS, respectively.
//       p1,p2 and p3 are two-component matrices representing the y and z coordinates
//       of the TD vertices in the TDCS, respectively.
//        y,z,x
        double x=x_obs[1];//y
        double y=x_obs[2];//z
        double z=x_obs[0];//x

        // Calculate barycentric coordinates
        double a, b, c;

        double denominator = (p2[1] - p3[1]) * (p1[0] - p3[0]) + (p3[0] - p2[0]) * (p1[1] - p3[1]);

        a= ((p2[1] - p3[1]) * (x - p3[0]) + (p3[0] - p2[0]) * (y - p3[1])) / denominator;
        b = ((p3[1] - p1[1]) * (x - p3[0]) + (p1[0] - p3[0]) * (y - p3[1])) / denominator;
        c = 1 - a - b;

        // Determine trimode values
        int trimode=1;
        if (a <= 0.0 && b > c && c  > a ) {
                trimode= -1;
        } else if (b <= 0.0 && c > a && a > b) {
                trimode = -1;
        } else if (c <= 0.0 && a > b && b > c) {
                trimode = -1;
        } else if (a == 0.0 && b >= 0.0 && c >= 0.0) {
                trimode = 0;
        } else if (a >= 0.0 && b == 0.0 && c >= 0.0) {
                trimode = 0;
        } else if (a >= 0.0 && b >= 0.0 && c == 0.0) {
                trimode = 0;
        }
        if (trimode == 0 && z != 0.0) {
                trimode = 1;
        }

        return trimode;
    }

// Function to perform TD setup and displacement calculations in TDCS
    il::StaticArray2D<double,3,3>  TDSetupD(const il::StaticArray<double, 3> &x_obs,
                 const double alpha,const double nu,const il::Array<double>& TriVertex, const il::Array<double>& SideVec) {
        // note  - convention of positive DDs in opening.

        // Transformation matrix
        il::StaticArray2D<double,2,2> A{0.},At{0.};
        A(0,0)=SideVec[2];
        A(0,1)=-SideVec[1];
        A(1,0)=SideVec[1];
        A(1,1)=SideVec[2];

        for (il::int_t i=0;i<2;i++){
            At(i,0)=A(0,i);
            At(i,1)=A(1,i);
        }

        il::StaticArray<double,2> x_aux{0.};
        x_aux[0]=x_obs[1]-TriVertex[1];
        x_aux[1]=x_obs[2]-TriVertex[2];
        auto r1 = il::dot(A,x_aux);
        double y1 = r1[0];
        double z1 = r1[1];

        // Transform in-plane slip vector components from TDCS to ADCS
//        il::Array<double> b_aux{2,1}; //#by , bz
//        vector<double> r2 = il::dot(A,b_aux)
//        by1 = r2[0];
//        bz1 = r2[1];

        il::StaticArray<double,2> b_aux_1{0.},b_aux_2{0.};
        b_aux_1[0]=1; // by in tdcs
        b_aux_2[1]=1; // bz in tdcs
        il::StaticArray<double,2> r2_1=il::dot(A,b_aux_1); // by in adcs
        il::StaticArray<double,2> r2_2=il::dot(A,b_aux_2); // bz in adcs

        // Calculate displacements in ADCS
        il::StaticArray2D<double,3,3> disp_adcs_aux= AngDisDisp(x_obs[0],y1, z1, -il::pi+alpha,nu); // in adcs
        il::StaticArray2D<double,3,3> disp_adcs{0.};
        // add in-plane effects
        for(il::int_t i=0;i<3;i++){
            for (il::int_t j=0;j<3;j++){
                if (j==0){
                    disp_adcs(i,j)=disp_adcs_aux(i,j);
                }
                if (j==1){  // U_i1 = U_i1(adcs)*by_acds[1]  + U_i2(adcs)*by_acds[2]
                    disp_adcs(i,j)= disp_adcs_aux(i,1)*r2_1[0]+disp_adcs_aux(i,2)*r2_1[1];
                }
                if (j==2){
                    disp_adcs(i,j)= disp_adcs_aux(i,1)*r2_2[0]+disp_adcs_aux(i,2)*r2_2[1];
                }
            }
        }

        il::StaticArray2D<double,3,3> disp_tdcs{0.};

        // Transform from ADCS to  TDCS
        for (il::int_t j=0;j<3;j++){
            il::StaticArray<double,2> vw;
            vw[0]=disp_adcs(1,j);
            vw[1]=disp_adcs(2,j);
            il::StaticArray<double,2> r3 = il::dot(At,vw);
            disp_tdcs(0,j)= disp_adcs(0,j);
            disp_tdcs(1,j)=r3[0];
            disp_tdcs(2,j)=r3[1];
        }
        return disp_tdcs;
    }
    //
    il::StaticArray2D<double,3,3> Displacements_TDCS(const il::StaticArray<double, 3> &x_obs,
                                                     const il::StaticArray2D<double, 3, 3> &P,
                                                     const double nu) {
        // compute displacements influence due to a triangular dislocation in an elastic full space
        // in the TDCS TRIANGULAR DISLOCATION COORDINATE SYSTEM... of Nikkho & Walker
        //
        // The x, y and z axes of the triangular dislocation coordinate system (TDCS) are normal to the plane,
        //  parallel to the strike direction and parallel to the dip vector, respectively.
        // The vertex P2 of the TD is the origin of the TDCS.
        // x_obs :: x observation point in the TDCS !!
        //
        // This routine for fundamental displacemnet WORKS only for a triangle with a normal = e_3 !!
        //
        // note  - convention of positive DDs in opening.
        IL_EXPECT_MEDIUM(P(0,2)==P(1,2) && P(0,2)==P(2,2) && P(1,2)==P(2,2) );
        il::Array<double> p1{3,0.0},p2{3,0.0},p3{3,0.0};
        for (il::int_t i=0;i<3;i++){
            p1[i]=P(0,i);
            p2[i]=P(1,i);
            p3[i]=P(2,i);
           // std::cout << " P1 " << p1[i] << "P2 " << p2[i] << "P3 "  << p3[i] << "\n";
        }
        // Calculate unit vectors TDCS
        il::Array<double> p2mp1{3,0.},p3mp1{3,0.},p3mp2{3,0.};
        for (il::int_t i=0;i<3;i++){
            p2mp1[i]=p2[i]-p1[i];
            p3mp1[i]=p3[i]-p1[i];
            p3mp2[i]=p3[i]-p2[i];
        }
        double n_p2mp1 = il::norm(p2mp1, il::Norm::L2);
        double n_p3mp1 = il::norm(p3mp1, il::Norm::L2);
        double n_p3mp2 = il::norm(p3mp2, il::Norm::L2);
        il::Array<double> e12{3},e13{3},e23{3};
        for (il::int_t i=0;i<3;i++){
            e12[i]=p2mp1[i]/n_p2mp1;
            e13[i]=p3mp1[i]/n_p3mp1;
            e23[i]=p3mp2[i]/n_p3mp2;
        }
        il::Array<double> m_e12{3},m_e13{3},m_e23{3};
        for (il::int_t i=0;i<3;i++){
            m_e12[i]=-e12[i];
            m_e13[i]=-e13[i];
            m_e23[i]=-e23[i];
        }
        // Calculate TD angles using dot product
        double A = acos(il::dot(e12, e13));
        double B = acos(-il::dot(e12, e23));
        double C = acos(il::dot(e23, e13));
       // std::cout << " A " << A << "B " << B << " C " << C <<"\n";
        il::Array<double> p1_a{2},p2_a{2},p3_a{2};
        for(il::int_t i=0;i<2;i++){ // y,z
            p1_a[i]=p1[i+1];
            p2_a[i]=p2[i+1];
            p3_a[i]=p3[i+1];
        }
        int Trimode = trimodefinder(x_obs,p1_a,p2_a,p3_a);
        //std::cout << " trimode !" << Trimode <<"\n";
        // if Trimode ==0 -> NAN
       // IL_EXPECT_FAST(Trimode != 0);
        if (Trimode ==0 ){
            std::cout << "Error obs point on a vertex !\n";
            il::abort();
        }

        il::StaticArray2D<double,3,3> disp_1Tp,disp_2Tp,disp_3Tp;
        if (Trimode == 1){
            // CONFIGURATION I
           // % Calculate first angular dislocation contribution
            disp_1Tp = TDSetupD(x_obs,A,  nu,p1, m_e13); //TDSetupD(xp,yp,zp,A,bx,by,bz,nu,p1,-e13);

          //  % Calculate second angular dislocation contribution
            disp_2Tp = TDSetupD(x_obs,B,nu,p2,e12); // TDSetupD(xp,yp,zp,B,bx,by,bz,nu,p2,e12);

            // % Calculate third angular dislocation contribution
            disp_3Tp = TDSetupD(x_obs,C,nu,p3,e23);// TDSetupD(xp,yp,zp,C,bx,by,bz,nu,p3,e23);
        }
        else {
            if (Trimode == -1){
                // CONFIGURATION II
               // % Calculate first angular dislocation contribution
              disp_1Tp =  TDSetupD(x_obs,A,  nu,p1, e13);
              //  % Calculate second angular dislocation contribution
              disp_2Tp = TDSetupD(x_obs,B,nu,p2,m_e12);
             //   % Calculate third angular dislocation contribution
              disp_3Tp = TDSetupD(x_obs,C,nu,p3,m_e23);
            }
        }
        // % Calculate the Burgers' function contribution corresponding to the TD
        il::Array<double> a{il::value,{-x_obs[0],p1[1]-x_obs[1],p1[2]-x_obs[2]}};
        il::Array<double> b{il::value,{-x_obs[0],-x_obs[1],-x_obs[2]}};
        il::Array<double> c{il::value,{-x_obs[0],p3[1]-x_obs[1],p3[2]-x_obs[2]}};
        double na = il::norm(a,il::Norm::L2);
        double nb = il::norm(b,il::Norm::L2);
        double nc = il::norm(c,il::Norm::L2);
      //  std::cout << "na " << na << " nb " << nb << "nc " << nc <<"\n";
        double aux = (na * nb  * nc  +  il::dot(a,b)*nc + il::dot(a,c)*nb + il::dot(b,c)*na) ;
        //std::cout << "aux " << aux << "\n";
        double aux_1 = (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]));
       // std::cout << "aux 1 " << aux_1 << "\n";
        double Fi = -2.0*atan2(aux_1,aux)/4./il::pi;
        //std::cout << " Fi :" << Fi << "\n" ;
        il::StaticArray2D<double,3,3> disp_all{0.};

        for (il::int_t j=0;j<3;j++){
            for(il::int_t i=0;i<3;i++){
                disp_all(i,j)=disp_1Tp(i,j)+disp_2Tp(i,j)+disp_3Tp(i,j);
            }
        }
        for (int j=0;j<3;j++){
            disp_all(j,j)= disp_all(j,j)+Fi;
        }

        return disp_all;
    }


    il::StaticArray2D<double,3,3> Displacements_EFCS(const il::StaticArray<double, 3> &x_obs,
                                                     const il::StaticArray2D<double, 3, 3> &P,
                                                     const double nu){
        // compute displacements influence due to a triangular dislocation in an elastic full space
        // in the EFCS TRIANGULAR DISLOCATION COORDINATE SYSTEM... of Nikkho & Walker
        //  X- East, Y-North, Z-UP
        // note that for some strange reason the displacement is uy,-ux, uz !! (see just below)
        //
        //% Calculate unit strike, dip and normal to TD vectors: For a horizontal TD
        //% as an exception, if the normal vector points upward, the strike and dip
        //% vectors point Northward and Westward, whereas if the normal vector points
        //% downward, the strike and dip vectors point Southward and Westward,
        //% respectively.
        //
        // note  - convention of positive DDs in opening.
        il::StaticArray<double,3> P2mP1{0.},P3mP1{0.};
        il::StaticArray<double,3> eY{il::value,{0.,1.,0.}},eZ{il::value,{0.,0.,1.}};
        for (int i=0;i<3;i++){
            P2mP1[i]=P(1,i)-P(0,i);
            P3mP1[i]=P(2,i)-P(0,i);
        }
        il::StaticArray<double,3> Vnorm = il::cross(P2mP1,P3mP1);
        for (int i=0;i<3;i++){
            Vnorm[i]=Vnorm[i]/il::norm(Vnorm,il::Norm::L2);
        }
        il::StaticArray<double,3> Vstrike = il::cross(eZ,Vnorm);
        double vstrike_norm = il::norm(Vstrike,il::Norm::L2);
        if (vstrike_norm==0.){
            for (int i=0;i<3;i++){
                Vstrike[i]=eY[i]*Vnorm[2];
            }
        }
        vstrike_norm = il::norm(Vstrike,il::Norm::L2);
        for (int i=0;i<3;i++){
            Vstrike[i]=Vstrike[i]/vstrike_norm;
        }
        il::StaticArray<double,3> Vdip = il::cross(Vnorm,Vstrike);
        il::StaticArray2D<double,3,3> At{0.};
        for (int j=0;j<3;j++){
            At(0,j)=Vnorm[j];
            At(1,j)=Vstrike[j];
            At(2,j)=Vdip[j];
        }
        // Transform coordinates to the TDCS
        il::StaticArray<double,3> p1,p2,p3,x_obs_tdsc;
        for (int i=0;i<3;i++){
            p1[i]=P(0,i)-P(1,i);
            p2[i]=P(1,i)-P(1,i);
            p3[i]=P(2,i)-P(1,i);
            x_obs_tdsc[i]=x_obs[i]-P(1,i);
        }
        p1 = il::dot(At,p1);
        p2 = il::dot(At,p2);
        p3 = il::dot(At,p3);
        x_obs_tdsc=il::dot(At,x_obs_tdsc);
        il::StaticArray2D<double,3,3> P_tdcs{0.};
        for (int j=0;j<3;j++){
            P_tdcs(0,j)=p1[j];
            P_tdcs(1,j)=p2[j];
            P_tdcs(2,j)=p3[j];
           // std::cout << " x obs tdcs " << x_obs_tdsc[j] << "\n";
            //std::cout << " P1 " << p1[j] << "P2 " << p2[j] << "P3 "  << p3[j] << "\n";
        }
        il::StaticArray2D<double, 3,3> disp_tdcs = Displacements_TDCS(x_obs_tdsc,P_tdcs,nu);
        // return to EFCS // be careful with switch of columns for source identity !
        //   disp_efcs_i^X = disp_tdcs_i^y, disp_efcs_i^Y = disp_tdcs_i^z, disp_efcs_i^z = disp_tdcs_i^x

        // then we need to change the displacement coordinates system
        il::StaticArray2D<double,3,3> A{0.};
        for (int j=0;j<3;j++){
            A(j,0)=Vnorm[j];
            A(j,1)=Vstrike[j];
            A(j,2)=Vdip[j];
        }
        // A.Uij.AT
        il::StaticArray2D<double, 3,3>  disp_efcs = il::dot(il::dot(A,disp_tdcs),At);

        return disp_efcs ;
    }



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Fundamental displacement kernel = displacement influence coefficients
    il::StaticArray2D<double, 3, 3> DisplacementKernelT0_af(const il::StaticArray<double, 3> &x,
                                                         const il::StaticArray2D<double, 3, 3> &xv,
                                                         const double nu) {
        // Compute displacement influence in the local T0 coordinates system of the element as per Fata convention
        ///  Usual T0 coordinate system is defined as:
        // origin on first vertex, and e_1 along first-second vertex,  and e_3 as normal to the element plane
        // inputs
        //   -x: source/receiver point's coordinates in the global reference system in
        //   the form (x,y,z) -xv: vertices' coordinates in the global reference
        //   system, in the following 2D array form:
        //     x0 y0 z0
        //     x1 y1 z1
        //     x2 y2 z2
        //   -nu: Poisson's ratio
        // (note the shear modulus does not enter here)
        // output
        //   Displacement = 3 x 3 matrix with the displacement influence coefficients in the local T0 coordinates system
        //   arranged as: U1, U2, U3 -> for rows
        //   DD1 (shear along e1 local), DD2 (shear along e2 local), DD3 (normal) -> for columns

        // note : here we switch to positive DD in overlap convention for consistency with other DD kernels

        bigwham::Triangle<0> elt;
        il::Array2D<double> xyz{3,3,0.};
        il::Array<double> x_obs{3,0.};
       // il::StaticArray2D<double, 3, 3> el_vertices_s_static{0};
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
                xyz(i, j) = xv(i, j);
            }
            x_obs[j]=x[j]-xv(0,j);
        }
        elt.SetElement(xyz);
        auto x_obs_local = elt.ConvertToLocal(x_obs);
        il::Array<double> y1{3},y2{3},y3{3};
        for (int i = 0; i < 3; i++) {
            // with origin at y1
            y1[i]=xyz(0,i)-xv(0,i);
            y2[i]=xyz(1,i)-xv(0,i);
            y3[i]=xyz(2,i)-xv(0,i);
        }

        //
        il::Array<double> P1=elt.ConvertToLocal(y1);
        il::Array<double> P2=elt.ConvertToLocal(y2);
        il::Array<double> P3=elt.ConvertToLocal(y3);

        il::StaticArray2D<double, 3, 3> Ps{0.};
        il::StaticArray<double, 3> x_obs_tdcs{0.};
        for (il::int_t j=0;j<3;j++){
            Ps(0,j)=P1[j];
            Ps(1,j)=P2[j];
            Ps(2,j)=P3[j];
            x_obs_tdcs[j]=x_obs_local[j];
          //  std::cout << " P1 " << Ps(0,j)  << " P2 " << Ps(1,j) << " P3 " << Ps(2,j) << " obs x"  << x_obs_tdcs[j]<< "\n";
        }

        // displacement in the Local coordinates system of the element !! needs to be converted back to the global frame ....
        il::StaticArray2D<double,3,3> disp=Displacements_EFCS(x_obs_tdcs,Ps, nu);
        // switch to positive DD in overlap
        for (int  j=0;j<3;j++){
            for (int  i=0;i<3;i++){
                disp(i,j)=-1.*  disp(i,j);
            }
        }

        return disp;
    }


}