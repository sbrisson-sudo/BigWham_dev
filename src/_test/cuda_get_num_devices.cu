#include <stdio.h>
#include <cuda_runtime.h>

int main() {
    int count = 0;
    cudaError_t error = cudaGetDeviceCount(&count);
    printf("Return code: %d (%s)\n", error, cudaGetErrorString(error));
    printf("Number of GPUs: %d\n", count);
    
    // Also try getting properties for each device
    for (int i = 0; i < count; i++) {
        cudaDeviceProp prop;
        error = cudaGetDeviceProperties(&prop, i);
        printf("Device %d: %s (Error: %s)\n", i, prop.name, cudaGetErrorString(error));
    }
    
    return 0;
}