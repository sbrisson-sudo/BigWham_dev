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
    T* FR_standard_size_data;
    T* FR_non_standard_size_data;
    int FR_standard_size_data_buffer_size;
    int FR_non_standard_size_data_buffer_size;
    std::vector<int> FR_std_orderedIndices; // to cpmply with BSR : blocks sorted by row then column
    int num_FR_std_blocks_;
    std::vector<int> FR_non_std_indices;

    std::vector<int> LR_sizes;

    std::unordered_map<int, std::vector<std::vector<int>>> LR_blocks_groups_indices;
    std::unordered_map<int, std::vector<T*>> LR_A_data;
    std::unordered_map<int, std::vector<T*>> LR_B_data;
    std::unordered_map<int, std::vector<int>> LR_data_buffer_sizes;

    // Num available GPUs
    int num_gpus_;
    
    std::vector<int> num_FR_per_gpu_;
    std::vector<int> num_LR_per_gpu_;
    std::vector<int> offsets_FR_gpu_;
    
    std::vector<std::unordered_map<int, std::vector<int>>> LR_group_per_gpu_;
    std::vector<std::pair<int, int>> LR_group_size_num_;

    // FOR ALL THE DATA STRUCTURES BELLOW, THE FIRST VECTOR RELATES
    // TO THE GPU REPARTITION

    std::vector<std::vector<std::unordered_map<int, std::vector<int>>>> LR_group_per_cuda_stream_;

    // GPU operation handler
    std::vector<cublasHandle_t> cublas_handle_;
    std::vector<cusparseHandle_t> cusparse_handle_;

    // Number of CUDA streams = 1 for BSR + the rest for batched opeartions
    // int num_streams_;
    const int num_streams_ = 32;
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
    std::vector<T*> d_y_partial_LR_;    // Partial results for LR blocks operations
    std::vector<T*> d_tmp_;          // tmp = B*x then y = A*tmp

    // FR data
    std::vector<T*> d_FR_data_;

    // LR data
    std::vector<std::unordered_map<int, std::unordered_map<int, T*>>> d_LR_A_data_;
    std::vector<std::unordered_map<int, std::unordered_map<int, T*>>> d_LR_B_data_;

    // FR metadata = BSR operations
    std::vector<cusparseMatDescr_t> FR_bsr_descr_;
    std::vector<int*> d_FR_bsrRowPtr_;  // Device array of row pointers
    std::vector<int*> d_FR_bsrColInd_;  // Device array of column indices

    // We also keep them on host
    std::vector<int*> h_FR_bsrRowPtr_; 
    std::vector<int*> h_FR_bsrColInd_; 

    // LR metadata = array of pointers for batched operations
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> d_LR_A_data_pointers_;   // to data 
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> d_LR_B_data_pointers_;   // to data 
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> d_LR_x_pointers_;      // to input vec (and to intermediate vec)
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> d_LR_y_pointers_;      // to output vec
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> d_LR_tmp_pointers_;      // to output vec

    // We also keep them on host
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> h_LR_A_data_pointers_;    
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> h_LR_B_data_pointers_;    
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> h_LR_x_pointers_;       
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> h_LR_y_pointers_;    
    std::vector<std::unordered_map<int, std::unordered_map<int, T**>>> h_LR_tmp_pointers_;    

    // To gather the low rank partial results
    size_t y_partial_LR_buffer_size_bytes_;
    size_t tmp_buffer_size_bytes_;
    std::vector<int*> d_LR_y_partial_src_indices_;
    std::vector<int*> d_LR_y_partial_dest_indices_;
    std::vector<int*> d_LR_y_partial_lengths_;

public:
  HmatCuda() = default;
  HmatCuda(const bigwham::MatrixGenerator<T> &matrix_gen, const double epsilon_aca, const int n_openMP_threads, const int num_GPUs, const bool verbose=true, const int fixed_rank=-1);
  ~HmatCuda();

  void copyToDevice(); // Just to not end up with a 1000 lines constructor
  void deallocateOnDevice();

  il::Array<T> matvec(il::ArrayView<T> x) override;

  // debugging functions
  // il::Array2D<T> getFRBlockDataHost(int fr_block);
  // il::Array2D<T> getFRBlockDataDevice(int fr_block);

  // il::Array<int> getFRBlockRowPtrHost();
  // il::Array<int>  getFRBlockColIndHost();
  // int* getFRBlockColIndDevice(int fr_block);
  // il::Array<int> getFRBlockRowPtrDevice(int fr_block);
  // il::Array2D<T> getLRBlockDataHost(int fr_block);
  // il::Array2D<T> getLRBlockBDataDevice(int size, int group_id);
  // il::Array2D<T> getLRBlockADataDevice(int size, int group_id);
  

};  

} // namespace bigwham

#endif // BIGWHAM_HMAT_H

// #if defined(__clang__) && !defined(FMT_ICC_VERSION)
// #pragma clang diagnostic pop
// #endif
