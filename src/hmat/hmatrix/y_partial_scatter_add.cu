#include "y_partial_scatter_add.h"

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