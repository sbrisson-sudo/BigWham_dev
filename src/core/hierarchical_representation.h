//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 26.01.23.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland, Geo-Energy Laboratory, 2016-2023.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef BIGWHAM_HIERARCHICAL_REPRESENTATION_H
#define BIGWHAM_HIERARCHICAL_REPRESENTATION_H

#pragma once
#include <il/Timer.h>
#include <il/Array.h>

#include <src/core/BEMesh.h>
#include <src/hmat/hmatrix/toHPattern.h>
#include <src/hmat/cluster/cluster.h>

namespace bie{

struct HRepresentation{
    bie::HPattern pattern_;
    il::Array<il::int_t> permutation_0_; // for rows
    il::Array<il::int_t> permutation_1_; // for columns
    bool is_square_=true; // if true - square matrix, and only permutation_0_ is stored
};

template<class El>
bie::HRepresentation h_representation_square_matrix(bie::BEMesh<El>& mesh,il::int_t max_leaf_size,double eta) {
    bie::HRepresentation hr ;
    hr.is_square_=true;
    // creation of the cluster
    // first get all collocation points in the mesh
    il::Timer tt;
    il::Array2D<double> Xcol=mesh.getCollocationPoints();
    tt.Start();
    bie::Cluster cluster = bie::cluster(max_leaf_size, il::io, Xcol);
    std::cout << "Cluster tree creation time :  " << tt.time() << "\n";
    tt.Stop();
    tt.Reset();
    hr.permutation_0_=cluster.permutation;

    tt.Start();
    const il::Tree<bie::SubHMatrix, 4> block_tree = bie::hmatrixTreeIxI(Xcol, cluster.partition,eta);
    tt.Stop();
    std::cout << "Time for binary cluster tree construction  " << tt.time() <<"\n";
    std::cout << " binary cluster tree depth =" << block_tree.depth() << "\n";
    hr.pattern_= bie::createPattern(block_tree);
    hr.pattern_.nr=mesh.numberCollocationPoints();
    hr.pattern_.nc=mesh.numberCollocationPoints();

    std::cout << " Number of blocks =" << hr.pattern_.n_B << "\n";
    std::cout << " Number of full blocks =" << hr.pattern_.n_FRB << "\n";
    std::cout << " Number of low rank blocks =" << hr.pattern_.n_LRB << "\n";

    return hr;
}


    template<class El_s,class El_r>
    bie::HRepresentation h_representation_rectangular_matrix(bie::BEMesh<El_s>& source_mesh,bie::BEMesh<El_r>& receiver_mesh,il::int_t max_leaf_size,double eta) {
        bie::HRepresentation hr ;
        hr.is_square_=false;
        // creation of the cluster
        // first get all collocation points in the mesh
        il::Timer tt;
        il::Array2D<double> Xcol_source=source_mesh.getCollocationPoints();
        il::Array2D<double> Xcol_receiver=receiver_mesh.getCollocationPoints();

        tt.Start();
        bie::Cluster cluster_s = bie::cluster(max_leaf_size, il::io, Xcol_source);
        std::cout << "Cluster tree creation time for the source mesh :  " << tt.time() << "\n";
        tt.Stop();
        tt.Reset();
        hr.permutation_1_=cluster_s.permutation;  // sources permutation

        tt.Start();
        bie::Cluster cluster_r = bie::cluster(max_leaf_size, il::io, Xcol_receiver);
        std::cout << "Cluster tree creation time for the source mesh :  " << tt.time() << "\n";
        tt.Stop();
        tt.Reset();
        hr.permutation_0_=cluster_r.permutation;  // receivers permutation

        tt.Start();
        const il::Tree<bie::SubHMatrix, 4> block_tree = bie::hmatrixTreeIxJ(Xcol_receiver, cluster_r.partition,Xcol_source,cluster_s.partition,eta);
        tt.Stop();
        std::cout << "Time for binary cluster tree construction  " << tt.time() <<"\n";
        std::cout << " binary cluster tree depth =" << block_tree.depth() << "\n";
        hr.pattern_= bie::createPattern(block_tree);
        hr.pattern_.nr=receiver_mesh.numberCollocationPoints();
        hr.pattern_.nc=source_mesh.numberCollocationPoints();

        std::cout << " Number of blocks =" << hr.pattern_.n_B << "\n";
        std::cout << " Number of full blocks =" << hr.pattern_.n_FRB << "\n";
        std::cout << " Number of low rank blocks =" << hr.pattern_.n_LRB << "\n";

        return hr;
    }





}

#endif //BIGWHAM_HIERARCHICAL_REPRESENTATION_H
