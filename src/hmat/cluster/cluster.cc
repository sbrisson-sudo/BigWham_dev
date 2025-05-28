//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 08.09.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
#include <limits>
#include <cmath>
#include <ctime>
#include <fstream>

#include "omp.h"

#include <il/Tree.h>
#include <json.hpp>

#include "cluster.h"

#define TIMING
// #define DEBUG

namespace bigwham {

// Block Cluster Tree Creation - IxI :: square H-matrix
il::Tree<bigwham::SubHMatrix, 4> hmatrixTreeIxI(
    const il::Array2D<double> &node, const il::Tree<il::Range, 2> &range_tree,
    double eta, const bool homogeneous_size, const int fixed_rank) {
  il::Tree<bigwham::SubHMatrix, 4> hmatrix_tree{};

  if (homogeneous_size){
    bigwham::hmatrixTreeIxI_rec_size_conservative(node, range_tree, eta, hmatrix_tree.root(),
    range_tree.root(), range_tree.root(), il::io,
    hmatrix_tree, fixed_rank);

  } else {
    bigwham::hmatrixTreeIxI_rec(node, range_tree, eta, hmatrix_tree.root(),
    range_tree.root(), range_tree.root(), il::io,
    hmatrix_tree);
  }

  hmatrix_tree.setDepth();
  return hmatrix_tree;
};

void hmatrixTreeIxI_rec(const il::Array2D<double> &node,
                     const il::Tree<il::Range, 2> &range_tree, double eta,
                     il::spot_t s, il::spot_t s0, il::spot_t s1, il::io_t,
                     il::Tree<bigwham::SubHMatrix, 4> &hmatrix_tree) {
  const bool is_admissible =
      isAdmissible(node, eta, range_tree.value(s0), range_tree.value(s1));
  if (is_admissible) {
    hmatrix_tree.Set(s,
                     bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
                                         bigwham::HMatrixType::LowRank});
  } else {
    if (range_tree.hasChild(s0, 0) && range_tree.hasChild(s0, 1) &&
        range_tree.hasChild(s1, 0) && range_tree.hasChild(s1, 1)) {
      hmatrix_tree.Set(
              s, bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
                                     bigwham::HMatrixType::Hierarchical});
      hmatrix_tree.AddChild(s, 0);
      hmatrixTreeIxI_rec(node, range_tree, eta, hmatrix_tree.child(s, 0),
                         range_tree.child(s0, 0), range_tree.child(s1, 0),
                         il::io, hmatrix_tree);
      hmatrix_tree.AddChild(s, 1);
      hmatrixTreeIxI_rec(node, range_tree, eta, hmatrix_tree.child(s, 1),
                         range_tree.child(s0, 1), range_tree.child(s1, 0),
                         il::io, hmatrix_tree);
      hmatrix_tree.AddChild(s, 2);
      hmatrixTreeIxI_rec(node, range_tree, eta, hmatrix_tree.child(s, 2),
                         range_tree.child(s0, 0), range_tree.child(s1, 1),
                         il::io, hmatrix_tree);
      hmatrix_tree.AddChild(s, 3);
      hmatrixTreeIxI_rec(node, range_tree, eta, hmatrix_tree.child(s, 3),
                         range_tree.child(s0, 1), range_tree.child(s1, 1),
                         il::io, hmatrix_tree);
    } else if ((range_tree.hasChild(s0, 0) && range_tree.hasChild(s0, 1)) ||
               (range_tree.hasChild(s1, 0) && range_tree.hasChild(s1, 1))) {
      hmatrix_tree.Set(
              s, bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
                                     bigwham::HMatrixType::FullRank});
    } else {
      hmatrix_tree.Set(
              s, bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
                                     bigwham::HMatrixType::FullRank});
    }
  }
}

void hmatrixTreeIxI_rec_size_conservative(const il::Array2D<double> &node,
  const il::Tree<il::Range, 2> &range_tree, double eta,
  il::spot_t s, il::spot_t s0, il::spot_t s1, il::io_t,
  il::Tree<bigwham::SubHMatrix, 4> &hmatrix_tree, const int fixed_rank) {

  const bool is_admissible =
  isAdmissible(node, eta, range_tree.value(s0), range_tree.value(s1));

  // If using fixed rank : we have to ensure that both dimensions of the block are 
  // at least twice the fixed rank
  int size_range_0 = range_tree.value(s0).end - range_tree.value(s0).begin;
  int size_range_1 = range_tree.value(s1).end - range_tree.value(s1).begin;
  bool is_bigger_than_rank = (fixed_rank < 0) || ((size_range_0 > fixed_rank) && (size_range_1 > fixed_rank));

  // 3rd condition to be admissible : the two clusters have abput the same size (if the biggest one has children)
  bool block_aspect_ratio_ok = true;
  int max_aspect_ratio = 1;
  if (size_range_0 > max_aspect_ratio*size_range_1){
    if (range_tree.hasChild(s0, 0)){
      block_aspect_ratio_ok = false;
    }
  }

  if (size_range_1 > max_aspect_ratio*size_range_0){
    if (range_tree.hasChild(s1, 0)){
      block_aspect_ratio_ok = false;
    }
  }

  // std::cout << "size_range_0 = " << size_range_0 << "; size_range_1 = " << size_range_1 << "; fixed_rank = " << fixed_rank << " -> " << is_bigger_than_rank << std::endl;


  if (is_admissible && is_bigger_than_rank && block_aspect_ratio_ok) {

    // if (size_range_0 == size_range_1) {
    //   hmatrix_tree.Set(s,
    //     bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
    //                       bigwham::HMatrixType::LowRank});
    // } else {
    //   hmatrix_tree.Set(
    //     s, bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
    //                             bigwham::HMatrixType::FullRank});
    // }

      hmatrix_tree.Set(s,
        bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
                          bigwham::HMatrixType::LowRank});
    
  } else {
    if (range_tree.hasChild(s0, 0) && range_tree.hasChild(s0, 1) &&
    range_tree.hasChild(s1, 0) && range_tree.hasChild(s1, 1)) { // Both clusters have children

      hmatrix_tree.Set(s, bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
        bigwham::HMatrixType::Hierarchical});

      // If the two clusters are power of 2 but not the same size, we only go down the biggest one
      bool s0_power_of_2 = (size_range_0 & (size_range_0 - 1)) == 0;
      bool s1_power_of_2 = (size_range_1 & (size_range_1 - 1)) == 0;
      if (s0_power_of_2 && s1_power_of_2 && (size_range_0 != size_range_1)){

        if (size_range_0 > size_range_1){

          hmatrix_tree.AddChild(s, 0);
          hmatrixTreeIxI_rec_size_conservative(node, range_tree, eta, hmatrix_tree.child(s, 0),
          range_tree.child(s0, 0), s1,
          il::io, hmatrix_tree, fixed_rank);
          hmatrix_tree.AddChild(s, 1);
          hmatrixTreeIxI_rec_size_conservative(node, range_tree, eta, hmatrix_tree.child(s, 1),
          range_tree.child(s0, 1), s1,
          il::io, hmatrix_tree, fixed_rank);

        } else {

          hmatrix_tree.AddChild(s, 0);
          hmatrixTreeIxI_rec_size_conservative(node, range_tree, eta, hmatrix_tree.child(s, 0),
          s0, range_tree.child(s1, 0),
          il::io, hmatrix_tree, fixed_rank);
          hmatrix_tree.AddChild(s, 1);
          hmatrixTreeIxI_rec_size_conservative(node, range_tree, eta, hmatrix_tree.child(s, 1),
          s0, range_tree.child(s1, 1),
          il::io, hmatrix_tree, fixed_rank);

        }
      } else {
        // Otherwise we go down both
        hmatrix_tree.AddChild(s, 0);
        hmatrixTreeIxI_rec_size_conservative(node, range_tree, eta, hmatrix_tree.child(s, 0),
        range_tree.child(s0, 0), range_tree.child(s1, 0),
        il::io, hmatrix_tree, fixed_rank);
        hmatrix_tree.AddChild(s, 1);
        hmatrixTreeIxI_rec_size_conservative(node, range_tree, eta, hmatrix_tree.child(s, 1),
        range_tree.child(s0, 1), range_tree.child(s1, 0),
        il::io, hmatrix_tree, fixed_rank);
        hmatrix_tree.AddChild(s, 2);
        hmatrixTreeIxI_rec_size_conservative(node, range_tree, eta, hmatrix_tree.child(s, 2),
        range_tree.child(s0, 0), range_tree.child(s1, 1),
        il::io, hmatrix_tree, fixed_rank);
        hmatrix_tree.AddChild(s, 3);
        hmatrixTreeIxI_rec_size_conservative(node, range_tree, eta, hmatrix_tree.child(s, 3),
        range_tree.child(s0, 1), range_tree.child(s1, 1),
        il::io, hmatrix_tree, fixed_rank);

      }
    } else if ((range_tree.hasChild(s0, 0) && range_tree.hasChild(s0, 1)) ||
              (range_tree.hasChild(s1, 0) && range_tree.hasChild(s1, 1))) {
      hmatrix_tree.Set(
      s, bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
                              bigwham::HMatrixType::FullRank});

    } else if (!range_tree.hasChild(s0, 0) && !range_tree.hasChild(s0, 1) &&
      !range_tree.hasChild(s1, 0) && !range_tree.hasChild(s1, 1)) { // None does
        hmatrix_tree.Set(s, bigwham::SubHMatrix{range_tree.value(s0), range_tree.value(s1),
                            bigwham::HMatrixType::FullRank});

    } else {
      std::cerr << "Invalid cluster tree, see hmatrixTreeIxJ_rec_size_conservative" << std::endl;
    }
  }
}

// Block Cluster Tree Creation - IxJ :: rectangular H-matrix
il::Tree<bigwham::SubHMatrix, 4> hmatrixTreeIxJ(
    const il::Array2D<double> &node0, const il::Tree<il::Range, 2> &range_tree0,
    const il::Array2D<double> &node1, const il::Tree<il::Range, 2> &range_tree1,    double eta, const bool homogeneous_size) {

  il::Tree<bigwham::SubHMatrix, 4> hmatrix_tree{};

  if (homogeneous_size){
    bigwham::hmatrixTreeIxJ_rec_size_conservative(node0, range_tree0, node1, range_tree1, eta, hmatrix_tree.root(),
                                range_tree0.root(), range_tree1.root(), il::io,
                                hmatrix_tree);

  } else {
    bigwham::hmatrixTreeIxJ_rec(node0, range_tree0, node1, range_tree1, eta, hmatrix_tree.root(),
                              range_tree0.root(), range_tree1.root(), il::io,
                              hmatrix_tree);
  }
  hmatrix_tree.setDepth();
  return hmatrix_tree;
};

void hmatrixTreeIxJ_rec(const il::Array2D<double>& node0,
                        const il::Tree<il::Range, 2>& range_tree0,
                        const il::Array2D<double>& node1,
                        const il::Tree<il::Range, 2>& range_tree1, double eta,
                        il::spot_t s, il::spot_t s0, il::spot_t s1, il::io_t,
                        il::Tree<bigwham::SubHMatrix, 4>& hmatrix_tree){
  const bool is_admissible =
      isAdmissible(node0,node1, eta, range_tree0.value(s0), range_tree1.value(s1));

  if (is_admissible) {
    hmatrix_tree.Set(s,
                     bigwham::SubHMatrix{range_tree0.value(s0), range_tree1.value(s1),
                                         bigwham::HMatrixType::LowRank});
  } else {
    if (range_tree0.hasChild(s0, 0) && range_tree0.hasChild(s0, 1) &&
        range_tree1.hasChild(s1, 0) && range_tree1.hasChild(s1, 1)) {
      hmatrix_tree.Set(
              s, bigwham::SubHMatrix{range_tree0.value(s0), range_tree1.value(s1),
                                     bigwham::HMatrixType::Hierarchical});
      hmatrix_tree.AddChild(s, 0);
      hmatrixTreeIxJ_rec(node0, range_tree0,node1,range_tree1, eta, hmatrix_tree.child(s, 0),
                         range_tree0.child(s0, 0), range_tree1.child(s1, 0),
                         il::io, hmatrix_tree);
      hmatrix_tree.AddChild(s, 1);
      hmatrixTreeIxJ_rec(node0, range_tree0, node1,range_tree1,eta, hmatrix_tree.child(s, 1),
                         range_tree0.child(s0, 1), range_tree1.child(s1, 0),
                         il::io, hmatrix_tree);
      hmatrix_tree.AddChild(s, 2);
      hmatrixTreeIxJ_rec(node0, range_tree0,node1,range_tree1, eta, hmatrix_tree.child(s, 2),
                         range_tree0.child(s0, 0), range_tree1.child(s1, 1),
                         il::io, hmatrix_tree);
      hmatrix_tree.AddChild(s, 3);
      hmatrixTreeIxJ_rec(node0, range_tree0,node1,range_tree1, eta, hmatrix_tree.child(s, 3),
                         range_tree0.child(s0, 1), range_tree1.child(s1, 1),
                         il::io, hmatrix_tree);
    } else if ((range_tree0.hasChild(s0, 0) && range_tree0.hasChild(s0, 1)) ||
               (range_tree1.hasChild(s1, 0) && range_tree1.hasChild(s1, 1))) {
      hmatrix_tree.Set(
              s, bigwham::SubHMatrix{range_tree0.value(s0), range_tree1.value(s1),
                                     bigwham::HMatrixType::FullRank});
    } else {
      hmatrix_tree.Set(
              s, bigwham::SubHMatrix{range_tree0.value(s0), range_tree1.value(s1),
                                     bigwham::HMatrixType::FullRank});
    }
  }
};

void hmatrixTreeIxJ_rec_size_conservative(const il::Array2D<double>& node0,
  const il::Tree<il::Range, 2>& range_tree0,
  const il::Array2D<double>& node1,
  const il::Tree<il::Range, 2>& range_tree1, double eta,
  il::spot_t s, il::spot_t s0, il::spot_t s1, il::io_t,
  il::Tree<bigwham::SubHMatrix, 4>& hmatrix_tree){
    
  const bool is_admissible = isAdmissible(node0,node1, eta, range_tree0.value(s0), range_tree1.value(s1));

  // Check if the two ranges have the same size 
  const bool square_block = (range_tree0.value(s0).end - range_tree0.value(s0).begin) == (range_tree1.value(s1).end - range_tree1.value(s1).begin);

  if (is_admissible && square_block) {
    hmatrix_tree.Set(s, bigwham::SubHMatrix{range_tree0.value(s0), range_tree1.value(s1), bigwham::HMatrixType::LowRank});
  } else {
    
    if (range_tree0.hasChild(s0, 0) && range_tree0.hasChild(s0, 1) &&
      range_tree1.hasChild(s1, 0) && range_tree1.hasChild(s1, 1)) { // Both clusters have children

      hmatrix_tree.Set(
      s, bigwham::SubHMatrix{range_tree0.value(s0), range_tree1.value(s1),
                    bigwham::HMatrixType::Hierarchical});

      // std::cout << "4 sons = " << range_tree0.child(s0, 0).index << " | " << range_tree0.child(s0, 1).index << " | "<< range_tree1.child(s1, 0).index << " | "<< range_tree1.child(s0, 1).index << std::endl;
      
      // 4 block-cluster sons
      hmatrix_tree.AddChild(s, 0);
      hmatrixTreeIxJ_rec_size_conservative(node0, range_tree0,node1,range_tree1, eta, hmatrix_tree.child(s, 0),
        range_tree0.child(s0, 0), range_tree1.child(s1, 0),
        il::io, hmatrix_tree);

      hmatrix_tree.AddChild(s, 1);
      hmatrixTreeIxJ_rec_size_conservative(node0, range_tree0, node1,range_tree1,eta, hmatrix_tree.child(s, 1),
        range_tree0.child(s0, 1), range_tree1.child(s1, 0),
        il::io, hmatrix_tree);

      hmatrix_tree.AddChild(s, 2);
      hmatrixTreeIxJ_rec_size_conservative(node0, range_tree0,node1,range_tree1, eta, hmatrix_tree.child(s, 2),
        range_tree0.child(s0, 0), range_tree1.child(s1, 1),
        il::io, hmatrix_tree);

      hmatrix_tree.AddChild(s, 3);
      hmatrixTreeIxJ_rec_size_conservative(node0, range_tree0,node1,range_tree1, eta, hmatrix_tree.child(s, 3),
        range_tree0.child(s0, 1), range_tree1.child(s1, 1),
        il::io, hmatrix_tree);

      } else if (range_tree0.hasChild(s0, 0) && range_tree0.hasChild(s0, 1) &&
        !range_tree1.hasChild(s1, 0) && !range_tree1.hasChild(s1, 1)) { // Only s0 does 

          // 2 sons
          // std::cout << "2 sons = " << range_tree0.child(s0, 0).index << " | " << range_tree0.child(s0, 1).index << std::endl;
          hmatrix_tree.Set(s, bigwham::SubHMatrix{range_tree0.value(s0), range_tree1.value(s1),bigwham::HMatrixType::Hierarchical});

          hmatrix_tree.AddChild(s, 0);
          hmatrixTreeIxJ_rec_size_conservative(node0, range_tree0,node1,range_tree1, eta, hmatrix_tree.child(s, 0),
            range_tree0.child(s0, 0), s1,
            il::io, hmatrix_tree);

          hmatrix_tree.AddChild(s, 1);
          hmatrixTreeIxJ_rec_size_conservative(node0, range_tree0, node1,range_tree1,eta, hmatrix_tree.child(s, 1),
            range_tree0.child(s0, 1), s1,
            il::io, hmatrix_tree);
        

      } else if (!range_tree0.hasChild(s0, 0) && !range_tree0.hasChild(s0, 1) &&
      range_tree1.hasChild(s1, 0) && range_tree1.hasChild(s1, 1)) { // Only s1 does 
        
          // 2 sons
          // std::cout << "2 sons = " << range_tree1.child(s1, 0).index << " | " << range_tree1.child(s1, 1).index << std::endl;
          hmatrix_tree.Set(s, bigwham::SubHMatrix{range_tree0.value(s0), range_tree1.value(s1),bigwham::HMatrixType::Hierarchical});

          hmatrix_tree.AddChild(s, 0);
          hmatrixTreeIxJ_rec_size_conservative(node0, range_tree0,node1,range_tree1, eta, hmatrix_tree.child(s, 0),
            s0, range_tree1.child(s1, 0),
            il::io, hmatrix_tree);

          hmatrix_tree.AddChild(s, 1);
          hmatrixTreeIxJ_rec_size_conservative(node0, range_tree0, node1,range_tree1,eta, hmatrix_tree.child(s, 1),
            s0, range_tree1.child(s1, 1),
            il::io, hmatrix_tree);
    
      } else if (!range_tree0.hasChild(s0, 0) && !range_tree0.hasChild(s0, 1) &&
      !range_tree1.hasChild(s1, 0) && !range_tree1.hasChild(s1, 1)) { // None does

          hmatrix_tree.Set(s, bigwham::SubHMatrix{range_tree0.value(s0), range_tree1.value(s1), bigwham::HMatrixType::FullRank});

      } else {
        std::cerr << "Invalid cluster tree, see hmatrixTreeIxJ_rec_size_conservative" << std::endl;
      }
  }
   

};

// Binary Cluster tree creation
Cluster cluster(il::int_t leaf_size, il::io_t, il::Array2D<double> &node, const bool homogeneous_size) {
  const il::int_t nb_nodes = node.size(0);

  bigwham::Cluster ans{};

  il::spot_t s = ans.partition.root();
  ans.partition.Set(s, il::Range{0, nb_nodes});

  ans.permutation.Resize(nb_nodes);
  for (il::int_t i = 0; i < nb_nodes; ++i) {
    ans.permutation[i] = i;
  }

  #ifdef TIMING
  struct timespec start, end;
  double duration;
  clock_gettime(CLOCK_MONOTONIC, &start);
  #endif // TIMING 

  #ifdef DEBUG
  int nb_thread = omp_get_max_threads();
  int max_rec_depth_openmp = static_cast<int>(std::log2(nb_thread));
  std::cout << "[DEBUG] clustering number of available threads = " << nb_thread << ", max recursive OpenMP depth = " << max_rec_depth_openmp << "\n";
  std::cout << "[DEBUG] homogeneous_size = " << homogeneous_size << std::endl;
  #endif

  if (homogeneous_size){
    #pragma omp parallel
    {
        #pragma omp single
        {
            // Only one thread starts the recursion
            cluster_rec_size_conservative(s, leaf_size, il::io, ans.partition, node, ans.permutation);
          }
    }
  } else {
    #pragma omp parallel
    {
        #pragma omp single
        {
            // Only one thread starts the recursion
            cluster_rec(s, leaf_size, il::io, ans.partition, node, ans.permutation);
        }
      }
  }

  #ifdef TIMING
  clock_gettime(CLOCK_MONOTONIC, &end);
  duration = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;    
  std::cout << "[Timing] clustering = " << duration*1000 << "ms\n";
  clock_gettime(CLOCK_MONOTONIC, &start);
  #endif // TIMING 

  // // Save tree to json
  // saveTreeToJSON(ans.partition);
  
  ans.partition.setDepth();
  return ans;
}

void cluster_rec(il::spot_t s, il::int_t leaf_size, il::io_t,
                 il::Tree<il::Range, 2> &tree, il::Array2D<double> &node,
                 il::Array<il::int_t> &permutation, int current_depth) {
  const il::int_t i_begin = tree.value(s).begin;
  const il::int_t i_end = tree.value(s).end;

  if (i_end - i_begin <= leaf_size) {
    return;
  }

  ////////////////////////////////////////////
  // Find the dimension in which we will split
  ////////////////////////////////////////////
  const il::int_t nb_nodes = node.size(0);
  const il::int_t dim = node.size(1);
  il::Array<double> width_box{dim};
  il::Array<double> middle_box{dim};

  for (il::int_t d = 0; d < dim; ++d) {
    double coordinate_minimum = std::numeric_limits<double>::max();
    double coordinate_maximum = -std::numeric_limits<double>::max();
    for (il::int_t i = i_begin; i < i_end; ++i) {
      if (node(i, d) < coordinate_minimum) {
        coordinate_minimum = node(i, d);
      }
      if (node(i, d) > coordinate_maximum) {
        coordinate_maximum = node(i, d);
      }
    }
    width_box[d] = coordinate_maximum - coordinate_minimum;
    middle_box[d] = coordinate_minimum + width_box[d] / 2;
  }

  double width_maximum = 0.0;
  il::int_t d_max = -1;
  for (il::int_t d = 0; d < dim; ++d) {
    if (width_box[d] >= width_maximum) {
      width_maximum = width_box[d];
      d_max = d;
    }
  }

  // ////////////////////
  // // Reorder the nodes
  // ////////////////////
  // const double middle = middle_box[d_max];
  // il::Array<double> tmp_node{dim};

  // il::int_t j = i_begin;
  // for (il::int_t i = i_begin; i < i_end; ++i) { 
  //   if (node(i, d_max) < middle) {
  //     // Swap node(i) and node (j)
  //     for (il::int_t d = 0; d < dim; ++d) {
  //       tmp_node[d] = node(i, d);
  //     }
  //     const il::int_t index = permutation[i];
  //     for (il::int_t d = 0; d < dim; ++d) {
  //       node(i, d) = node(j, d);
  //     }
  //     permutation[i] = permutation[j];
  //     for (il::int_t d = 0; d < dim; ++d) {
  //       node(j, d) = tmp_node[d];
  //     }
  //     permutation[j] = index;
  //     ++j;
  //   }
  // }
  // if (j == i_begin) {
  //   ++j;
  // } else if (j == i_end) {
  //   --j;
  // }

    ////////////////////
  // Reorder the nodes
  // Sorting them along the splitting coordinate
  ////////////////////
  const double middle = middle_box[d_max];
  il::Array<double> tmp_node{dim};

  // Sort all nodes from i_begin to i_end based on d_max coordinate
  for (il::int_t i = i_begin; i < i_end; ++i) {
    for (il::int_t k = i + 1; k < i_end; ++k) {
        if (node(i, d_max) > node(k, d_max)) {
            // Swap node(i) and node(k)
            for (il::int_t d = 0; d < dim; ++d) {
                tmp_node[d] = node(i, d);
                node(i, d) = node(k, d);
                node(k, d) = tmp_node[d];
            }
            
            // Swap permutation indices
            const il::int_t index = permutation[i];
            permutation[i] = permutation[k];
            permutation[k] = index;
        }
    }
  }

  // Compute child cardinal 
  int j = (i_begin+i_end)/2;

  #ifdef DEBUG
  std::cout << i_end - i_begin << " -> " << j - i_begin << " + " << i_end - j << std::endl;
  #endif

  int nb_thread = omp_get_max_threads();
  int max_rec_depth_openmp = static_cast<int>(std::log2(nb_thread));

  il::spot_t s0, s1;

  #pragma omp critical
  {

    tree.AddChild(s, 0);
    s0 = tree.child(s, 0);
    tree.AddChild(s, 1);
    s1 = tree.child(s, 1);
    tree.Set(s0, il::Range{i_begin, j});
    tree.Set(s1, il::Range{j, i_end});
  }

  current_depth ++;

  // Parallelize only if we're at a depth where it makes sense
  if (current_depth <= max_rec_depth_openmp) {
    // Process the first branch in a new task
    #pragma omp task shared(tree, node, permutation)
    {
      cluster_rec(s0, leaf_size, il::io, tree, node, permutation, current_depth);
    }
    
    // Process the second branch in the current thread
    cluster_rec(s1, leaf_size, il::io, tree, node, permutation, current_depth);
    
    // Wait for all tasks at this level to complete before continuing
    #pragma omp taskwait
  } else {
    // If we're already at a deep level, just use sequential execution
    cluster_rec(s0, leaf_size, il::io, tree, node, permutation, current_depth);
    cluster_rec(s1, leaf_size, il::io, tree, node, permutation, current_depth);
  }

  // cluster_rec(s0, leaf_size, il::io, tree, node, permutation);
  // cluster_rec(s1, leaf_size, il::io, tree, node, permutation);
}

// Binary Cluster tree creation
// This one ensure that for a max leaf size of m 
// the cardinal of the first child #sigma_1 is 
// #sigma_1 = (#tau / m // 2) * m 
void cluster_rec_size_conservative(il::spot_t s, il::int_t leaf_size, il::io_t,
  il::Tree<il::Range, 2> &tree, il::Array2D<double> &node,
  il::Array<il::int_t> &permutation, int current_depth) {

  // #ifdef DEBUG
  // int thread_id = omp_get_thread_num();
  // std::cout << "[Thread " << thread_id << "] entering clustering - depth = " << current_depth << std::endl;
  // #endif

  const il::int_t i_begin = tree.value(s).begin;
  const il::int_t i_end = tree.value(s).end;

  if (i_end - i_begin <= leaf_size) {
    // std::cout << "Cluster node of size " << i_end - i_begin << std::endl;
    return;
  }

  ////////////////////////////////////////////
  // Find the dimension in which we will split
  ////////////////////////////////////////////
  const il::int_t nb_nodes = node.size(0);
  const il::int_t dim = node.size(1);
  il::Array<double> width_box{dim};
  il::Array<double> middle_box{dim};

  for (il::int_t d = 0; d < dim; ++d) {
  double coordinate_minimum = std::numeric_limits<double>::max();
  double coordinate_maximum = -std::numeric_limits<double>::max();
  for (il::int_t i = i_begin; i < i_end; ++i) {
    if (node(i, d) < coordinate_minimum) {
      coordinate_minimum = node(i, d);
    }
    if (node(i, d) > coordinate_maximum) {
      coordinate_maximum = node(i, d);
    }
  }
  width_box[d] = coordinate_maximum - coordinate_minimum;
  middle_box[d] = coordinate_minimum + width_box[d] / 2;
  }

  double width_maximum = 0.0;
  il::int_t d_max = -1;
  for (il::int_t d = 0; d < dim; ++d) {
    if (width_box[d] >= width_maximum) {
      width_maximum = width_box[d];
      d_max = d;
    }
  }

  ////////////////////
  // Reorder the nodes
  // Sorting them along the splitting coordinate
  ////////////////////
  const double middle = middle_box[d_max];
  il::Array<double> tmp_node{dim};

  // Sort all nodes from i_begin to i_end based on d_max coordinate
  for (il::int_t i = i_begin; i < i_end; ++i) {
    for (il::int_t k = i + 1; k < i_end; ++k) {
        if (node(i, d_max) > node(k, d_max)) {
            // Swap node(i) and node(k)
            for (il::int_t d = 0; d < dim; ++d) {
                tmp_node[d] = node(i, d);
                node(i, d) = node(k, d);
                node(k, d) = tmp_node[d];
            }
            
            // Swap permutation indices
            const il::int_t index = permutation[i];
            permutation[i] = permutation[k];
            permutation[k] = index;
        }
    }
  }

  // Compute child cardinal 
  int k = 0;
  while (leaf_size * pow(2, k + 1) <= (i_end - i_begin)) {
      k++;
  }

  il::int_t card_s1;
  if ((i_end - i_begin) == (leaf_size * pow(2, k))){
    card_s1 = (i_end - i_begin)/2;
  } else {
    card_s1 = leaf_size * pow(2, k);;
  }

  // std::cout << i_end - i_begin << " -> " << card_s1 << " + " << i_end - i_begin - card_s1 << std::endl;

  int nb_thread = omp_get_max_threads();
  int max_rec_depth_openmp = static_cast<int>(std::log2(nb_thread));

  // if current_depth <= max_rec_depth_openmp -> one of the recursive calls is made in a different thread

  // If current_depth
  il::spot_t s0, s1;

  #pragma omp critical
  {
    tree.AddChild(s, 0); 
    s0 = tree.child(s, 0);
    tree.AddChild(s, 1);
    s1 = tree.child(s, 1);
    tree.Set(s0, il::Range{i_begin, i_begin + card_s1});
    tree.Set(s1, il::Range{i_begin + card_s1, i_end});
  }

  current_depth ++;

  // Parallelize only if we're at a depth where it makes sense
  if (current_depth <= max_rec_depth_openmp) {
    // Process the first branch in a new task
    #pragma omp task shared(tree, node, permutation)
    {
      cluster_rec_size_conservative(s0, leaf_size, il::io, tree, node, permutation, current_depth);
    }
    
    // Process the second branch in the current thread
    cluster_rec_size_conservative(s1, leaf_size, il::io, tree, node, permutation, current_depth);
    
    // Wait for all tasks at this level to complete before continuing
    #pragma omp taskwait
  } else {
    // If we're already at a deep level, just use sequential execution
    cluster_rec_size_conservative(s0, leaf_size, il::io, tree, node, permutation, current_depth);
    cluster_rec_size_conservative(s1, leaf_size, il::io, tree, node, permutation, current_depth);
  }

  // cluster_rec_size_conservative(s0, leaf_size, il::io, tree, node, permutation, current_depth);
  // cluster_rec_size_conservative(s1, leaf_size, il::io, tree, node, permutation, current_depth);
}

nlohmann::json nodeToJSON(il::Tree<il::Range, 2> tree, il::spot_t node){
  nlohmann::json result;

  if (!tree.hasChild(node, 0)){
    // No child
    il::Range leaf_range = tree.value(node);
    std::vector<int> range = {static_cast<int>(leaf_range.begin), static_cast<int>(leaf_range.end)};
    result["range"] = range;
  } else {
    // Has children nodes
    result["left"] = nodeToJSON(tree, tree.child(node, 0));
    result["right"] = nodeToJSON(tree, tree.child(node, 1));
  }

  return result;
}

void saveTreeToJSON(il::Tree<il::Range, 2> tree){

  nlohmann::json tree_json = nodeToJSON(tree, tree.root());

  std::ofstream file("bigwham_binary_cluster_tree.json");
  file << tree_json.dump(2);

  std::cout << "Clustering binary tree saved in bigwham_binary_cluster_tree.json" << std::endl;
}

int count_cluster_leaves(il::Tree<il::Range, 2>& tree, il::spot_t s){

    if (!tree.hasChild(s)) {
        return 1;
    }

    int count = 0;

    // Binary tree depth search
    for (int i = 0; i < 2; ++i) {
        il::spot_t s_child = tree.child(s, i);
        count += count_cluster_leaves(tree, s_child);
    }

    return count;
}

}  // namespace bigwham
