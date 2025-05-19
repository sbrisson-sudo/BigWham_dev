//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 08.09.21.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne), Switzerland,
// Geo-Energy Laboratory, 2016-2025.  All rights reserved. See the LICENSE.TXT
// file for more details.
//
// last modifications 5.2.21: Moving to std::unique_ptr (C. Peruzzo)

#ifndef BIGWHAM_HMAT_CUDA_H
#define BIGWHAM_HMAT_CUDA_H

#include <vector>
#include <unordered_map>
#include <utility>      // std::pair
#include <functional>   // std::hash

#include <cuda_runtime.h>
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cuda_runtime_api.h>

#include <omp.h>

#include "hmat/arrayFunctor/matrix_generator.h"
#include "hmat/compression/adaptiveCrossApproximation.h"
#include "hmat/hmatrix/LowRank.h"
#include "hmat/hmatrix/h_pattern.h"
#include "hmat/hierarchical_representation.h"

#include "contiguous_array2d_vector.hpp"
#include "hmat.h"

// #define USE_MAGMA

#ifdef USE_MAGMA
#include "magma_v2.h"
#include "magmablas_d_v1.h"
#endif

struct pair_hash {
    std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        // A better hash combiner (from Boost)
        std::size_t h1 = std::hash<int>{}(p.first);
        std::size_t h2 = std::hash<int>{}(p.second);
        return h1 ^ (h2 << 1); // Or use std::hash<std::tuple<int, int>> if available
    }
};

namespace bigwham {

/**
 * @class HmatCuda 
 * @brief Extends Hmat to allow for GPU matvec operation using Cuda (via Cublas and Cusparse)
 *
 * This is done by overidding the build and matvec operations
 * - the main change in the build operation is that we want to allocate memory buffers 
 * before building the blocks, then copy them in memory. We want different buffers for 
 * different matrices sizes
 * - the matvec operation calls to cublas operations (except for the non standard sized blocks ?)
 */
template <typename T> 
class HmatCuda : public Hmat<T>{
private:
    void buildCuda(const bigwham::MatrixGenerator<T> & matrix_gen, const double epsilon);
    void buildFRCuda(const bigwham::MatrixGenerator<T> & matrix_gen);
    template <il::int_t dim> 
    void buildLRCuda(const bigwham::MatrixGenerator<T> & matrix_gen,const double epsilon);

    int total_block_size_;

    // CPU (host) memory buffers

    // for full rank (FR) blocks
    T* FR_standard_size_data;           // standard size FR blocks data
    T* FR_non_standard_size_data;       // non standard size FR blocks data
    int FR_standard_size_data_buffer_size;      // standard size FR blocks data size
    int FR_non_standard_size_data_buffer_size;  // non standard size FR blocks data size
    std::vector<int> FR_std_orderedIndices; // to cpmply with BSR : blocks sorted by row then column
    int num_FR_std_blocks_;
    std::vector<int> FR_non_std_indices;

    // To store the diagonal
    bool diagonal_stored_ = false;
    std::vector<T> diagonal_;

    // for low rank (LR) blocks
    std::unordered_map<int, int> num_LR_std_blocks_per_size_;
    T* LR_non_standard_size_A_data_;
    T* LR_non_standard_size_B_data_;
    int LR_non_standard_size_data_A_buffer_size_;
    int LR_non_standard_size_data_B_buffer_size_;
    std::unordered_map<int, T*> LR_standard_size_A_data_;
    std::unordered_map<int, T*> LR_standard_size_B_data_;
    std::unordered_map<int, int> LR_standard_size_data_buffer_sizes_;

    std::unordered_map<int, std::vector<int>> LR_std_indices_;
    std::vector<int> LR_non_std_indices_;

    std::unordered_map<std::pair<int, int>, int, pair_hash> num_FR_nonstd_blocks_per_size_;
    std::unordered_map<std::pair<int, int>, int, pair_hash> num_LR_nonstd_blocks_per_size_;

    bool FR_deallocated_on_host_ = false;

    // Num available GPUs
    int num_gpus_;

    // FOR ALL THE DATA STRUCTURES BELLOW, THE FIRST VECTOR RELATES
    // TO THE GPU REPARTITION
    std::vector<int> num_FR_per_gpu_;
    std::vector<int> offsets_FR_gpu_;    
    std::vector<int> num_LR_per_gpu_;
    std::vector<std::vector<int>> LR_std_sizes_per_gpu_;

    // GPU operation handler
    std::vector<cublasHandle_t> cublas_handle_;
    std::vector<cusparseHandle_t> cusparse_handle_;

    // Number of CUDA streams = 1 for FR + 1 for LR
    const int num_streams_ = 2;
    std::vector<std::vector<cudaStream_t>> cuda_streams_;

#ifdef USE_MAGMA
    std::vector<std::vector<magma_queue_t>> magma_queues_;
#endif

    // GPU (device) memory buffers
    // Vectors
    size_t vector_size_bytes_;
    size_t vector_size_;
    std::vector<T*> d_x_;               // Store the lhs vector on device
    std::vector<T*> d_y_;               // Store the rhs vector on device
    T* d_x_tmp_;                       // For permutations / agregations
    // T* d_y_tmp_;                       // For permutations / agregations
    std::vector<T*> d_y_partial_LR_;    // Partial results for LR blocks operations
    std::vector<T*> d_y_partial_FR_;    // Partial results for FR blocks operations
    std::vector<T*> d_tmp_;             // tmp = B*x then y = A*tmp

    // FR data
    std::vector<T*> d_FR_data_;

    // LR data, for each size
    std::vector<std::unordered_map<int, T*>> d_LR_A_data_;
    std::vector<std::unordered_map<int, T*>> d_LR_B_data_;

    // FR metadata = BSR operations
    std::vector<cusparseMatDescr_t> FR_bsr_descr_;
    std::vector<int*> d_FR_bsrRowPtr_;  // Device array of row pointers
    std::vector<int*> d_FR_bsrColInd_;  // Device array of column indices

    // We also keep them on host
    std::vector<int*> h_FR_bsrRowPtr_; 
    std::vector<int*> h_FR_bsrColInd_; 

    // FR metadata to use cublas operations and not cusparse
    std::vector<T**> d_FR_data_pointers_;
    std::vector<T**> d_FR_x_pointers_;
    std::vector<T**> d_FR_y_partial_pointers_;

    // FR : gather partial results
    std::vector<size_t> y_partial_FR_buffer_size_bytes_;
    std::vector<int*> d_FR_y_partial_src_indices_;
    std::vector<int*> d_FR_y_partial_dest_indices_;
    std::vector<int*> d_FR_y_partial_lengths_;

    // LR metadata = array of pointers for batched operations
    std::vector<std::unordered_map<int, T**>> d_LR_A_data_pointers_;    // to data 
    std::vector<std::unordered_map<int, T**>> d_LR_B_data_pointers_;    // to data 
    std::vector<std::unordered_map<int, T**>> d_LR_x_pointers_;         // to input vec 
    std::vector<std::unordered_map<int, T**>> d_LR_y_pointers_;         // to output vec
    std::vector<std::unordered_map<int, T**>> d_LR_tmp_pointers_;       // to tmp output vec

    // We also keep them on host
    std::vector<std::unordered_map<int, T**>> h_LR_A_data_pointers_;    
    std::vector<std::unordered_map<int, T**>> h_LR_B_data_pointers_;    
    std::vector<std::unordered_map<int, T**>> h_LR_x_pointers_;       
    std::vector<std::unordered_map<int, T**>> h_LR_y_pointers_;    
    std::vector<std::unordered_map<int, T**>> h_LR_tmp_pointers_;    

    // To gather the low rank partial results
    std::vector<size_t> y_partial_LR_buffer_size_bytes_;
    std::vector<size_t> tmp_buffer_size_bytes_;
    std::vector<int*> d_LR_y_partial_src_indices_;
    std::vector<int*> d_LR_y_partial_dest_indices_;
    std::vector<int*> d_LR_y_partial_lengths_;

    // Array of pointers for the batched operations on the non std blocks 
    // Note that it is all carried by the first GPU
    std::unordered_map<std::pair<int, int>, T**, pair_hash> d_FR_nonStd_data_pointers_;
    std::unordered_map<std::pair<int, int>, T**, pair_hash> d_FR_nonStd_x_pointers_;
    std::unordered_map<std::pair<int, int>, T**, pair_hash> d_FR_nonStd_y_partial_pointers_;

    std::unordered_map<std::pair<int, int>, T**, pair_hash> d_LR_nonStd_A_data_pointers_;
    std::unordered_map<std::pair<int, int>, T**, pair_hash> d_LR_nonStd_B_data_pointers_;
    std::unordered_map<std::pair<int, int>, T**, pair_hash> d_LR_nonStd_x_pointers_;
    std::unordered_map<std::pair<int, int>, T**, pair_hash> d_LR_nonStd_y_partial_pointers_;
    std::unordered_map<std::pair<int, int>, T**, pair_hash> d_LR_nonStd_tmp_pointers_;

    // And the data buffers
    T* d_FR_non_std_data_;
    T* d_LR_non_std_A_data_;
    T* d_LR_non_std_B_data_;

    // To keep track of the offsets to use for the non std blocks
    // Offsets on data
    std::unordered_map<std::pair<int, int>, int, pair_hash> FR_non_std_offsets_data_,LR_non_std_A_offsets_data_,LR_non_std_B_offsets_data_;
    // Offsets on aux vectors
    std::unordered_map<std::pair<int, int>, int, pair_hash> FR_non_std_offsets_y_,LR_non_std_offsets_y_,LR_non_std_offsets_tmp_;

    // To store the permutations on device
    int* d_permutation_0_;
    int* d_permutation_1_;

public:
  HmatCuda() = default;
  HmatCuda(const bigwham::MatrixGenerator<T> &matrix_gen, const double epsilon_aca, const int n_openMP_threads, const int num_GPUs, const bool verbose=true, const int fixed_rank=-1);
  ~HmatCuda();

  void copyToDevice(); // Just to not end up with a 1000 lines constructor
  void deallocateOnHost(bool deallocate_FR_host);
  void deallocateOnDevice();

  void setDiagonal();
  std::vector<T> diagonal() override;
  std::vector<T> diagonalOriginal() override;

  il::int_t nbOfEntries() override;

  il::Array<T> matvec(il::ArrayView<T> x) override;
  void matvec(); // When data already on device in d_x_ 
  void matvec(T* x, T* y) override; // When data already on device
  void matvec(T* x); // When data already on device
  void matvecOriginal(T* x, T* y) override; // Apply permutation on device

};  

} // namespace bigwham

#endif // BIGWHAM_HMAT_H