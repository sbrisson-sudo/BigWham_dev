#ifndef SCATTER_ADD_H
#define SCATTER_ADD_H

#include <cuda_runtime.h>

void scatter_add(double* d_result,
                 const double* d_source_data,
                 const int* d_source_indices,
                 const int* d_target_indices,
                 const int* d_lengths,
                 int num_vectors);

#endif // SCATTER_ADD_H