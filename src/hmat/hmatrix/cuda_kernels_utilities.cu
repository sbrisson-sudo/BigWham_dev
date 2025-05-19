#include <iostream>

#include "cuda_kernels_utilities.h"

__global__ void scatter_add_kernel(double* result,
                                   const double* source_data,
                                   const int* source_indices,
                                   const int* target_indices,
                                   const int* lengths,
                                   int num_vectors) {
    int vector_idx = blockIdx.x;
    int thread_idx = threadIdx.x;
    
    // Each block processes one vector
    if (vector_idx < num_vectors) {
        int vector_length = lengths[vector_idx];
        int source_start = source_indices[vector_idx];
        int target_start = target_indices[vector_idx];
        
        // Each thread in the block handles one element of the vector
        for (int i = thread_idx; i < vector_length; i += blockDim.x) {
            atomicAdd(&result[target_start + i], source_data[source_start + i]);
        }
    }
}

void scatter_add(double* d_result,
                 const double* d_source_data,
                 const int* d_source_indices,
                 const int* d_target_indices,
                 const int* d_lengths,
                 int num_vectors) {
    
    int threads_per_block = 256;
    int blocks = num_vectors;
    
    scatter_add_kernel<<<blocks, threads_per_block>>>(
        d_result, d_source_data, d_source_indices, d_target_indices, d_lengths, num_vectors);
}

// Permutation on device
__global__ void forward_permute_kernel(
    const double*  x,
    double*  z,
    const int*  permutation,
    int n_elmts,
    int dof_dimension)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_elmts) return;

    int src_offset = dof_dimension * permutation[i];
    int dst_offset = dof_dimension * i;

    for (int j = 0; j < dof_dimension; ++j) {
        z[dst_offset + j] = x[src_offset + j];
    }
}

void forward_permute(
    const double*  x,
    double*  z,
    const int*  permutation_1,
    int hmat_size,
    int dof_dimension
    )
{
    int threads = 256;
    int n_elmts = hmat_size / dof_dimension;
    int blocks = (n_elmts + threads - 1) / threads;

    forward_permute_kernel<<<blocks, threads>>>(x, z, permutation_1, n_elmts, dof_dimension);
}

__global__ void backward_permute_kernel(
    const double*  y,
    double*  yout,
    const int*  permutation,
    int n_points,
    int dof_dimension)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_points) return;

    int dst_offset = dof_dimension * permutation[i];
    int src_offset = dof_dimension * i;

    for (int j = 0; j < dof_dimension; ++j) {
        yout[dst_offset + j] = y[src_offset + j];
    }
}

void backward_permute(
    const double*  y,
    double*  yout,
    const int*  permutation_0,
    int hmat_size,
    int dof_dimension
    )
{
    int threads = 256;
    int n_elmts = hmat_size / dof_dimension;
    int blocks = (n_elmts + threads - 1) / threads;

    backward_permute_kernel<<<blocks, threads>>>(y, yout, permutation_0, n_elmts, dof_dimension);
}