#ifndef SCATTER_ADD_H
#define SCATTER_ADD_H

#include <cuda_runtime.h>

void scatter_add(double* d_result,
                 const double* d_source_data,
                 const int* d_source_indices,
                 const int* d_target_indices,
                 const int* d_lengths,
                 int num_vectors);

void forward_permute(
    const double*  x,
    double*  z,
    const int*  permutation_1,
    int hmat_size,
    int dof_dimension);

void backward_permute(
    const double*  y,
    double*  yout,
    const int*  permutation_0,
    int hmat_size,
    int dof_dimension);

#endif // SCATTER_ADD_H