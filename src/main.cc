//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 22.12.19.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <iostream>
#include <string>

#include "io/bigwham_io.h"

/////////////////////////////////////////////////////////////////////////////////
int main() {
    // create a simple mesh for a griffith crack -
    // use the bigwhamio interface.
    ///  simple mesh
    // check diag
    int n_elts = 4;
    std::vector<double> coor(2 * (n_elts + 1), 0.);
    double L = 1.;
    double h = 2. * L / n_elts;
    int k = 0;
    for (int i = 0; i < n_elts + 1; i++) {
        coor[k] = i * h - L;
        k = k + 2;
    }
    std::vector<int> conn(n_elts * 2, 0.);
    k = 0;
    for (int i = 0; i < n_elts; i++) {
        conn[k] = i;
        conn[k + 1] = i + 1;
        k = k + 2;
    }

    std::vector<double> properties{1., 0.};

    bool verbose = true;
    bool homegeneous_block_size = false;
    bool useCuda = false;
    BigWhamIO my_io{coor, conn, "2DS0-H", properties, verbose, homegeneous_block_size, useCuda};
    my_io.BuildHierarchicalMatrix(32, 2, 1.e-3);
    auto nt = my_io.GetOmpThreads() ;
    std::cout <<"numb of threads " << nt <<"\n";
    il::Array<double> x(my_io.MatrixSize(1),0.);
    for(il::int_t i=0;i<n_elts ;i++){
         x[2*i+1]=4.0*sqrt(L*L-0.25*(coor[2*i]+coor[2*(i+1)])*(coor[2*i]+coor[2*(i+1)]));
     }
    // obs mesh ....
    int n_obs=3;
    std::vector<double> obs_coor(2*n_obs,0.);
    k=0;
    for (int i=0;i<n_obs;i++) {
        obs_coor[k+1]=i*h+h; //  x=0, y-axis
        k=k+2;
    }
    for (int i=0;i<n_obs;i++){
        std::cout << "x " << obs_coor[2*i] << " y " << obs_coor[2*i+1] << "\n";
    }
    il::Array<double> test= my_io.ComputeDisplacements(obs_coor, x.view());
    il::Array<double> test_stress= my_io.ComputeStresses(obs_coor, x.view());

    for (int i=0;i<n_obs;i++){
        std::cout << "u_x " << test[i*2+0] <<"\n";
        std::cout << "u_y " << test[i*2+1] <<"\n";
    }

    for (int i=0;i<n_obs;i++){
        std::cout << "s_xx " << test_stress[i*3] <<"\n";
        std::cout << "s_yy " << test_stress[i*3+1] <<"\n";
        std::cout << "s_xy " << test_stress[i*3+2] <<"\n";
    }

    std::cout << "\n End of BigWham - exe "   << "\n";
}

/////////////////////////////////////////////////////////////////////////////////
