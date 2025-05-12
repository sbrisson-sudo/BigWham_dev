#include <unordered_map>
#include <ctime>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>

// #define USE_NVTX

#ifdef USE_NVTX
#include <nvtx3/nvtx3.hpp>
#endif

#include "hmat_cuda.h"
#include "y_partial_scatter_add.h"

#include "cnpy.h"

// #define USE_CUSPARSE_BSRMV

// #define TIMING
// #define DEBUG
// #define WRITE_TMP_RES
// #define WRITE_LR_DATA_NPY
// #define WRITE_INPUT_OUTPUT_VEC
// #define PRINT_N_GROUPS
// #define DEBUG_MULT_GPU

// Error checking helper functions
#define CHECK_CUDA_ERROR(val) check_cuda((val), #val, __FILE__, __LINE__)
void check_cuda(cudaError_t err, const char* const func, const char* const file,
    const int line) {
    if (err != cudaSuccess) {
    std::cerr << "CUDA Runtime Error at: " << file << ":" << line
            << std::endl;
    std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
    exit(1);
    }
}

#define CHECK_CUBLAS_ERROR(val) check_cublas((val), #val, __FILE__, __LINE__)
static const char* cublasGetErrorString(cublasStatus_t error) {
    switch (error) {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
        default:
            return "Unknown cuBLAS error";
    }
}
void check_cublas(cublasStatus_t err, const char* const func, const char* const file,
  const int line) {
  if (err != CUBLAS_STATUS_SUCCESS) {
  std::cerr << "cuBLAS Error at: " << file << ":" << line
          << std::endl;
  std::cerr << cublasGetErrorString(err) << " " << func << std::endl;
  exit(1);
  }
}

std::string formatBytes(size_t bytes) {
    const double KB = 1024.0;
    const double MB = KB * 1024.0;
    const double GB = MB * 1024.0;

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);

    if (bytes >= GB) {
        oss << (bytes / GB) << " GB";
    } else if (bytes >= MB) {
        oss << (bytes / MB) << " MB";
    } else if (bytes >= KB) {
        oss << (bytes / KB) << " KB";
    } else {
        oss << bytes << " B";
    }

    return oss.str();
}

namespace bigwham {

template <typename T>
HmatCuda<T>::HmatCuda(const bigwham::MatrixGenerator<T> & matrix_gen, const double epsilon_aca, const int n_openMP_threads, const int num_GPUs, const bool verbose, const int fixed_rank) {

#ifdef TIMING
    struct timespec start, end;
    double duration;
    clock_gettime(CLOCK_MONOTONIC, &start);
#endif // TIMING 

    // construction directly
    this->verbose_ = verbose;
    this->fixed_rank_ = fixed_rank;
    this->n_openMP_threads_=n_openMP_threads;
    this->num_gpus_ = num_GPUs;

    this->dof_dimension_ = matrix_gen.blockSize();
    this->size_[0] = matrix_gen.size(0);
    this->size_[1] = matrix_gen.size(1);
    
    // construction directly
    this->hr_ = matrix_gen.hr();
    int leaf_size = this->hr_->leaf_size;
    HPattern pattern = this->hr_->pattern_;

    // size of elements blocks
    const int dim = matrix_gen.blockSize();

    // FR MEM ALLOCATION
    // First we get the number of full rank blokcks of standard and non standard size 
    int num_FR_blocks_standard_size = 0;
    int num_FR_blocks_non_standard_size = 0;
    int total_size_non_standard_FR_blocks = 0;

    for (int i(0); i<pattern.n_FRB; i++){
        int i0 = pattern.FRB_pattern(1,i);
        int j0 = pattern.FRB_pattern(2,i);
        int iend = pattern.FRB_pattern(3,i);
        int jend = pattern.FRB_pattern(4,i);

        if ((iend-i0 == leaf_size) && (jend-j0 == leaf_size)){ 
            num_FR_blocks_standard_size++;
        } else {
            num_FR_blocks_non_standard_size++;
            total_size_non_standard_FR_blocks += (jend-j0)*(iend-i0);
            num_FR_nonstd_blocks_per_size_[{iend-i0, jend-j0}]++;
        }
    }

    this->num_FR_std_blocks_ = num_FR_blocks_standard_size;

    // Then we get the permutation to have the standard blocks arranged in row-major order, then sorted by their columns
    std::vector<std::tuple<int, int, int>> FR_std_blocks_pos;    
    FR_std_blocks_pos.reserve(num_FR_blocks_standard_size);
    FR_non_std_indices.reserve(num_FR_blocks_non_standard_size);
    for (int i = 0; i < pattern.n_FRB; i++) {
        int i0 = pattern.FRB_pattern(1, i);
        int j0 = pattern.FRB_pattern(2, i);
        int iend = pattern.FRB_pattern(3,i);
        int jend = pattern.FRB_pattern(4,i);
        if ((iend-i0 == leaf_size) && (jend-j0 == leaf_size)){ 
            FR_std_blocks_pos.push_back({i0, j0, i});
        } else {
            FR_non_std_indices.push_back(i);
        }
    }

    // Sort blocks by row index
    std::sort(FR_std_blocks_pos.begin(), FR_std_blocks_pos.end());

    // Get ordered indices
    FR_std_orderedIndices.reserve(FR_std_blocks_pos.size());
    for (const auto& block : FR_std_blocks_pos) {
        FR_std_orderedIndices.push_back(std::get<2>(block));
    }

    #ifdef DEBUG
    // Print it 
    std::cout << " FR blocks std order : [";
    for (auto i : FR_std_orderedIndices){
        std::cout << i << ", ";
    }
    std::cout << "]\n";

    std::cout << "FR_non_std_indices = [";
    for (int i : FR_non_std_indices) std::cout << i << ", ";
    std::cout << "]\n";
    #endif

    // We allocate the memory for standardized FR blocks
    FR_standard_size_data_buffer_size = num_FR_blocks_standard_size * leaf_size * leaf_size * dim*dim;
    this->FR_standard_size_data = new T[FR_standard_size_data_buffer_size];

    FR_non_standard_size_data_buffer_size = total_size_non_standard_FR_blocks * dim*dim;
    this->FR_non_standard_size_data = new T[FR_non_standard_size_data_buffer_size];

    // LR MEM ALLOCATION
    // Here also we put apart the non square blocks 
    std::unordered_map<int, std::vector<int>> LR_blocks_per_size_indices;

    int total_size_non_standard_LR_blocks = 0;

    for (int i(0); i<pattern.n_LRB; i++){
        int i0 = pattern.LRB_pattern(1,i);
        int j0 = pattern.LRB_pattern(2,i);
        int iend = pattern.LRB_pattern(3,i);
        int jend = pattern.LRB_pattern(4,i);

        if ((iend-i0 == jend-j0)){
            num_LR_std_blocks_per_size_[iend-i0]++; 
            LR_std_indices_[iend-i0].push_back(i);
        } else {
            LR_non_std_indices_.push_back(i);
            total_size_non_standard_LR_blocks += (iend-i0) * this->fixed_rank_;
            num_LR_nonstd_blocks_per_size_[{iend-i0, jend-j0}]++;
        }
    }
    // Just print it 
    #ifdef DEBUG
    std::cout << "Num of non square LR blocks " << LR_non_std_indices_.size() << std::endl;
    for (auto& [block_size, indices] : LR_std_indices_){
        std::cout << "Num of square LR blocks of size " << block_size << " : " << LR_std_indices_[block_size].size() << std::endl;
    }

    std::cout << "LR_non_std_indices_ = [";
    for (int i : LR_non_std_indices_) std::cout << i << ",";
    std::cout << "]\n";
    #endif

    // We allocate the memory for the low rank blocks

    // For the non square blocks :
    // We compute the total memory needed
    LR_non_standard_size_data_A_buffer_size_ = 0;
    LR_non_standard_size_data_B_buffer_size_ = 0;
    for (int i : LR_non_std_indices_){
        int i0 = pattern.LRB_pattern(1,i);
        int j0 = pattern.LRB_pattern(2,i);
        int iend = pattern.LRB_pattern(3,i);
        int jend = pattern.LRB_pattern(4,i);

        LR_non_standard_size_data_A_buffer_size_ += (iend-i0)*fixed_rank * dim*dim;
        LR_non_standard_size_data_B_buffer_size_ += (jend-j0)*fixed_rank * dim*dim;
    }
    // We allocate it 
    LR_non_standard_size_A_data_ = new T[LR_non_standard_size_data_A_buffer_size_];
    LR_non_standard_size_B_data_ = new T[LR_non_standard_size_data_B_buffer_size_];

    // For the square blocks : one memory buffer per size
    for (auto& [block_size, num_blocks] : num_LR_std_blocks_per_size_){
        size_t data_size = num_blocks * block_size*fixed_rank * dim*dim;
        T* LR_A_buffer = new T[data_size];
        T* LR_B_buffer = new T[data_size];
        LR_standard_size_A_data_[block_size] = LR_A_buffer;
        LR_standard_size_B_data_[block_size] = LR_B_buffer;
        LR_standard_size_data_buffer_sizes_[block_size] = data_size;

        #ifdef DEBUG
        std::cout << "Allocating " << data_size << "(" << data_size*sizeof(T) << "bytes) for " << num_blocks << " blocks) for size " << block_size << std::endl;
        #endif
    }


    // // Print total memory footprint
    // if (this->verbose_){
    //     // FR blocks
    //     size_t total_size = FR_standard_size_data_buffer_size + FR_non_standard_size_data_buffer_size;
    //     // non std LR blocks
    //     total_size += LR_non_standard_size_data_A_buffer_size_ + LR_non_standard_size_data_B_buffer_size_;
    //     // std LR blocks
    //     for (auto& [block_size, data_size] : LR_standard_size_data_buffer_sizes_){
    //         total_size += 2*data_size;
    //     }
    
    //     std::cout << "Total memory allocated on CPU = " << formatBytes(total_size*sizeof(double)) << " \n";
    // }

    #ifdef TIMING
    clock_gettime(CLOCK_MONOTONIC, &end);
    duration = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;    
    std::cout << "[Timing] allocated memory on host = " << duration*1000 << "ms\n";
    clock_gettime(CLOCK_MONOTONIC, &start);
    #endif // TIMING 

    // Now we construct the hmat (on host)
    il::Timer tt;
    tt.Start();
    this->buildCuda(matrix_gen, epsilon_aca);
    tt.Stop();
    if (this->verbose_){
        std::cout << "Creation of hmat done in " << tt.time() << " s\n";
        std::cout << "Compression ratio - " << this->compressionRatio() << "\n";
        std::cout << "Hmat object - built "<< "\n";
    }

    #ifdef TIMING
    clock_gettime(CLOCK_MONOTONIC, &end);
    duration = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;    
    std::cout << "[Timing] populated the hmat = " << duration*1000 << "ms\n";
    clock_gettime(CLOCK_MONOTONIC, &start);
    #endif // TIMING 

    // Here we print the number of each blocks
    if (this->verbose_){
        std::cout << "Num of standard FR blocks : " << num_FR_blocks_standard_size << std::endl;
        std::cout << "Num of non standard FR blocks : " << num_FR_blocks_non_standard_size << std::endl;
        int num_LR_std = 0;
        for (auto& [block_size, indices] : LR_std_indices_) num_LR_std += indices.size();
        std::cout << "Num of standard LR blocks : " << num_LR_std << std::endl;
        std::cout << "Num of non standard LR blocks : " << LR_non_std_indices_.size() << std::endl;

    }


    // And we copy it on device
    copyToDevice();

    // Before deallocating on the host memory we store the diagonal
    setDiagonal();

    // Once all the Hmat data is on device, we deallocate it on host memory
    bool deallocate_FR_host = false;
    deallocateOnHost(deallocate_FR_host);

    #ifdef TIMING
    clock_gettime(CLOCK_MONOTONIC, &end);
    duration = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;    
    std::cout << "[Timing] copied the hmat to device = " << duration*1000 << "ms\n";
    clock_gettime(CLOCK_MONOTONIC, &start);
    #endif // TIMING 

    std::cout << std::flush;
}

template <typename T>
void HmatCuda<T>::copyToDevice(){

    if (this->verbose_){
        std::cout << "--------------------" << std::endl;
        std::cout << "Copying the Hmat to the GPU ..." << std::endl;
        std::cout << "Num of available GPU = " << num_gpus_ << std::endl;
    }

    auto hr = this->hr_;
    int fixed_rank = this->fixed_rank_;
    int leaf_size = this->hr_->leaf_size;
    const int dim_dof = this->dof_dimension_;

    // Initializing MAGMA
    #ifdef USE_MAGMA
    magma_init();
    #endif

    // ---------------------------
    // Distributing load on GPUS
    // ---------------------------

    // FR standard blocks
    num_FR_per_gpu_.resize(num_gpus_);
    for (int gpu_id = 0; gpu_id < num_gpus_-1; gpu_id++) {
        num_FR_per_gpu_[gpu_id] = num_FR_std_blocks_ / num_gpus_;
    }
    int sum_so_far = 0;
    for (int gpu_id = 0; gpu_id < num_gpus_-1; gpu_id++) {
        sum_so_far += num_FR_per_gpu_[gpu_id];
    }
    num_FR_per_gpu_[num_gpus_-1] = num_FR_std_blocks_ - sum_so_far;

    #ifdef DEBUG_MULT_GPU
    std::cout << "num_FR_std_blocks_ = " << num_FR_std_blocks_ << std::endl;
    std::cout << "Num of FR blocks per GPU = [";
    for (int n : num_FR_per_gpu_) std::cout << n << ",";
    std::cout << "]\n";
    #endif

    offsets_FR_gpu_.resize(num_gpus_);
    offsets_FR_gpu_[0] = 0;
    for (int gpu_id(1); gpu_id<num_gpus_; gpu_id++){
        offsets_FR_gpu_[gpu_id] = offsets_FR_gpu_[gpu_id-1] + num_FR_per_gpu_[gpu_id-1];
    } 

    // Distributng the LR block groups 
    // We attribute to each GPU a set of sizes
    LR_std_sizes_per_gpu_.resize(num_gpus_);
    std::vector<int> gpu_load(num_gpus_, 0); 

    num_LR_per_gpu_.resize(num_gpus_);
    for (int gpu_id(0); gpu_id<num_gpus_; gpu_id++) num_LR_per_gpu_[gpu_id] = 0;

    for (auto& [block_size, indices] :  LR_std_indices_){
        int num_blocks =  LR_std_indices_[block_size].size();
        int load = num_blocks * block_size;

        // We find the less loaded gpu
        int min_load_gpu = 0;
        for (int i = 0; i < num_gpus_; ++i) {
            if (gpu_load[i] < gpu_load[min_load_gpu]) {
                min_load_gpu = i;
            }
        }

        // We attribute this size to the less loaded GPU
        gpu_load[min_load_gpu] += load;
        LR_std_sizes_per_gpu_[min_load_gpu].push_back(block_size);
        num_LR_per_gpu_[min_load_gpu] += num_blocks;
    }

    #ifdef DEBUG_MULT_GPU
    std::cout << "LR block sizes attribution per GPU :" << std::endl;
    for (int gpu_id(0); gpu_id<num_gpus_; gpu_id++){
        auto& sizes = LR_std_sizes_per_gpu_[gpu_id];
        std::cout << "GPU " << gpu_id << " has sizes ";
        for (int size : sizes) std::cout << size << ", ";
        std::cout << std::endl;
    }
    #endif

    // ---------------------------
    // Distributing load on GPUS - done
    // ---------------------------

    // ---------------------------
    // Initializing the per GPU data structures
    // ---------------------------

    // GPU operation handler
    cublas_handle_.resize(num_gpus_);
    cusparse_handle_.resize(num_gpus_);

    // Number of CUDA streams = 1 for BSR + the rest for batched opeartions
    cuda_streams_.resize(num_gpus_);

    #ifdef USE_MAGMA
    magma_queues_.resize(num_gpus_);
    #endif

    // GPU (device) memory buffers
    d_x_.resize(num_gpus_);               // Store the lhs vector on device
    d_y_.resize(num_gpus_);               // Store the rhs vector on device
    d_y_partial_LR_.resize(num_gpus_);    // Partial results for LR blocks operations
    d_tmp_.resize(num_gpus_);          // tmp = B*x then y = A*tmp
    d_y_partial_FR_.resize(num_gpus_); 

    // FR data
    d_FR_data_.resize(num_gpus_);

    // LR data
    d_LR_A_data_.resize(num_gpus_);
    d_LR_B_data_.resize(num_gpus_);

    // FR metadata = BSR operations
    FR_bsr_descr_.resize(num_gpus_);
    d_FR_bsrRowPtr_.resize(num_gpus_);  // Device array of row pointers
    d_FR_bsrColInd_.resize(num_gpus_);  // Device array of column indices

    // We also keep them on host
    h_FR_bsrRowPtr_.resize(num_gpus_); 
    h_FR_bsrColInd_.resize(num_gpus_); 

    // FR metadata for cublas operations
    d_FR_data_pointers_.resize(num_gpus_);
    d_FR_x_pointers_.resize(num_gpus_);
    d_FR_y_partial_pointers_.resize(num_gpus_);

    y_partial_FR_buffer_size_bytes_.resize(num_gpus_);
    d_FR_y_partial_src_indices_.resize(num_gpus_);
    d_FR_y_partial_dest_indices_.resize(num_gpus_);
    d_FR_y_partial_lengths_.resize(num_gpus_);

    // LR metadata = array of pointers for batched operations
    d_LR_A_data_pointers_.resize(num_gpus_);   // to data 
    d_LR_B_data_pointers_.resize(num_gpus_);   // to data 
    d_LR_x_pointers_.resize(num_gpus_);      // to input vec (and to intermediate vec)
    d_LR_y_pointers_.resize(num_gpus_);      // to output vec
    d_LR_tmp_pointers_.resize(num_gpus_);      // to output vec

    // We also keep them on host
    h_LR_A_data_pointers_.resize(num_gpus_);    
    h_LR_B_data_pointers_.resize(num_gpus_);    
    h_LR_x_pointers_.resize(num_gpus_);       
    h_LR_y_pointers_.resize(num_gpus_);    
    h_LR_tmp_pointers_.resize(num_gpus_);    

    // To gather the low rank partial results
    y_partial_LR_buffer_size_bytes_.resize(num_gpus_);
    tmp_buffer_size_bytes_.resize(num_gpus_);
    d_LR_y_partial_src_indices_.resize(num_gpus_);
    d_LR_y_partial_dest_indices_.resize(num_gpus_);
    d_LR_y_partial_lengths_.resize(num_gpus_);

    // ---------------------------
    // Initializing the per GPU data structures - done
    // ---------------------------

    // Loop on GPUS
    for (int gpu_id=0; gpu_id<num_gpus_; gpu_id++){

        CHECK_CUDA_ERROR(cudaSetDevice(gpu_id));

        #ifdef DEBUG3
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "GPU " << gpu_id << " - Free: " << free_mem << " / Total: " << total_mem << std::endl;
        #endif

        // ---------------------------
        // Initiating the CUDA streams
        // ---------------------------

        // Init the handles
        CHECK_CUBLAS_ERROR(cublasCreate(&cublas_handle_[gpu_id]));
        cusparseCreate(&cusparse_handle_[gpu_id]);

        // // First we create the CUDA streams
        if (this->verbose_) std::cout << "[GPU "<< gpu_id << "] Initiating " << num_streams_ << " CUDA streams\n";

        cuda_streams_[gpu_id].resize(num_streams_);
        for (int i = 0; i < num_streams_; i++) CHECK_CUDA_ERROR(cudaStreamCreate(&cuda_streams_[gpu_id][i]));

        // Here we associate the Cuda streams (but the first one used by Cusparse) to Magma queues       
        #ifdef USE_MAGMA
        magma_queues_[gpu_id].resize(num_streams_-1);
        for (int i = 1; i < num_streams_; i++){
            magma_queue_create_from_cuda(
                gpu_id,                  
                cuda_streams_[gpu_id][i],        
                NULL,              
                NULL,               
                &magma_queues_[gpu_id][i-1]      
            );
        }
        #endif


        // ---------------------------
        // Initiating the CUDA streams - done
        // ---------------------------

        // ---------------------------
        // Initialiazing the aux vectors 
        // ---------------------------

        vector_size_ = this->size_[0];
        vector_size_bytes_ = this->size_[0]* sizeof(T);

        // Here we need to compute the buffer size needed for tmp and y_partial for LR blocks
        y_partial_LR_buffer_size_bytes_[gpu_id] = 0;
        tmp_buffer_size_bytes_[gpu_id] = 0;

        // 1. std sized blocks
        for (int block_size : LR_std_sizes_per_gpu_[gpu_id]){
            int num_blocks = LR_std_indices_[block_size].size();
            y_partial_LR_buffer_size_bytes_[gpu_id] += num_blocks * block_size * dim_dof * sizeof(T);
            tmp_buffer_size_bytes_[gpu_id] += num_blocks * this->fixed_rank_ * dim_dof * sizeof(T);
        }

        // 2. non std sized blocks (on the 1st GPU)
        if (gpu_id == 0){
            for (auto& [sizes,count] : num_LR_nonstd_blocks_per_size_){
                y_partial_LR_buffer_size_bytes_[gpu_id] += count * sizes.first * dim_dof * sizeof(T);
                tmp_buffer_size_bytes_[gpu_id] += count * this->fixed_rank_ * dim_dof * sizeof(T);
            }
        }

        // Idem for FR blocks
        // 1. for standard size blocks
        y_partial_FR_buffer_size_bytes_[gpu_id] = num_FR_per_gpu_[gpu_id] * leaf_size * dim_dof * sizeof(T);
        // 2. on the 1st GPU we add the non std blocks
        if (gpu_id == 0){
            for (auto& [sizes,count] : num_FR_nonstd_blocks_per_size_){
                y_partial_FR_buffer_size_bytes_[gpu_id] += count * sizes.first * dim_dof * sizeof(T); 
            }
        }


        // Print associated memory consumption
        if (this->verbose_){
            std::cout << "[GPU "<< gpu_id << "] Allocating " << formatBytes(2*vector_size_bytes_ + y_partial_LR_buffer_size_bytes_[gpu_id] + tmp_buffer_size_bytes_[gpu_id] + y_partial_FR_buffer_size_bytes_[gpu_id]) << " on the GPU for auxilliary vectors" << std::endl;
        }

        // We allocate them
        CHECK_CUDA_ERROR(cudaMalloc(&d_x_[gpu_id], vector_size_bytes_));
        CHECK_CUDA_ERROR(cudaMalloc(&d_y_[gpu_id], vector_size_bytes_));

        CHECK_CUDA_ERROR(cudaMalloc(&d_y_partial_FR_[gpu_id], y_partial_FR_buffer_size_bytes_[gpu_id]));

        CHECK_CUDA_ERROR(cudaMalloc(&d_y_partial_LR_[gpu_id], y_partial_LR_buffer_size_bytes_[gpu_id]));
        CHECK_CUDA_ERROR(cudaMalloc(&d_tmp_[gpu_id], tmp_buffer_size_bytes_[gpu_id]));

        // ---------------------------
        // Initialiazing the aux vectors - done
        // ---------------------------

        // ----------------------------
        // COPYING FULL RANKS BLOCKS
        // ----------------------------
        // Then we copy the FR blocks data (only of standard size)

        size_t FR_block_size = leaf_size*leaf_size*dim_dof*dim_dof;
        size_t FR_data_size_bytes = num_FR_per_gpu_[gpu_id] * FR_block_size * sizeof(T);

        size_t FR_nonStd_data_size_bytes = FR_non_standard_size_data_buffer_size*sizeof(T);

        #ifdef DEBUG
        std::cout << "[GPU "<< gpu_id << "] Copying FR standard block to device = " << FR_data_size_bytes << " bytes\n";
        std::cout << "[GPU "<< gpu_id << "] offset[gpu_id] = " << offsets_FR_gpu_[gpu_id] << std::endl;
        std::cout << "[GPU "<< gpu_id << "] num_FR_per_gpu_[gpu_id] = " << num_FR_per_gpu_[gpu_id] << std::endl;
        std::cout << "[GPU "<< gpu_id << "] num_FR_per_gpu_[gpu_id] * FR_block_size = " << num_FR_per_gpu_[gpu_id] * FR_block_size << std::endl;
        std::cout << "[GPU "<< gpu_id << "] FR_standard_size_data_buffer_size = " << FR_standard_size_data_buffer_size << std::endl;
        #endif

        if (this->verbose_){
            // total FR blocks data
            size_t total_to_allocate = FR_data_size_bytes;
            if (gpu_id==0) total_to_allocate += FR_nonStd_data_size_bytes;
            // adding LR non std size data
            if (gpu_id==0) total_to_allocate += (LR_non_standard_size_data_A_buffer_size_ + LR_non_standard_size_data_B_buffer_size_)* sizeof(T);
            // and LR std data
            for (int block_size : LR_std_sizes_per_gpu_[gpu_id]){
                total_to_allocate +=  2*LR_standard_size_data_buffer_sizes_[block_size] * sizeof(T);
            }
            std::cout << "[GPU "<< gpu_id << "] Allocating " << formatBytes(total_to_allocate) << " on the GPU for hierarchical matrix data" << std::endl;
        }

        CHECK_CUDA_ERROR(cudaMalloc(&d_FR_data_[gpu_id], FR_data_size_bytes));
        CHECK_CUDA_ERROR(cudaMemcpy(d_FR_data_[gpu_id], FR_standard_size_data + offsets_FR_gpu_[gpu_id]*FR_block_size, FR_data_size_bytes, cudaMemcpyHostToDevice));

        // On the 1st GPU we copy the non std FR blocks
        if (gpu_id == 0){
            CHECK_CUDA_ERROR(cudaMalloc(&d_FR_non_std_data_, FR_nonStd_data_size_bytes));
            CHECK_CUDA_ERROR(cudaMemcpy(d_FR_non_std_data_, FR_non_standard_size_data, FR_nonStd_data_size_bytes, cudaMemcpyHostToDevice));
        }

        // Then we set up the CUSPARSE BSR description of the FR 
        #ifdef USE_CUSPARSE_BSRMV

            cusparseCreateMatDescr(&FR_bsr_descr_[gpu_id]);
            cusparseSetMatType(FR_bsr_descr_[gpu_id], CUSPARSE_MATRIX_TYPE_GENERAL);
            cusparseSetMatIndexBase(FR_bsr_descr_[gpu_id], CUSPARSE_INDEX_BASE_ZERO);

            // We set d_FR_bsrRowPtr_ = array of size n_block_rows + 1 
            // with the index of the first block of the row
            // see https://medium.com/gpgpu/block-sparse-matrix-vector-multiplication-with-cuda-4e616b30267

            // GETTING bsrRowPtr
            // First we compute the nb of block rows
            int num_block_rows = -1;
            
            for (int i : FR_std_orderedIndices) {
                int i0 = hr->pattern_.FRB_pattern(1, i);
                int iend = hr->pattern_.FRB_pattern(3, i);
                int blockRow = i0 / leaf_size; // Assuming blockDim is your block size
                num_block_rows = std::max(num_block_rows, blockRow);
            }
            num_block_rows++; // +1 because indices are 0-based

            // Store it 
            this->total_block_size_ = num_block_rows;

            h_FR_bsrRowPtr_[gpu_id] = new int[num_block_rows+1];
            for (int i = 0; i <= num_block_rows; i++) h_FR_bsrRowPtr_[gpu_id][i] = 0;

            // Count blocks per row
            for (int idx(0); idx<num_FR_per_gpu_[gpu_id]; idx++){
                int idx_offset = idx + offsets_FR_gpu_[gpu_id];
                int i = FR_std_orderedIndices[idx_offset];
                int i0 = hr->pattern_.FRB_pattern(1, i);
                int blockRow = i0 / leaf_size;
                h_FR_bsrRowPtr_[gpu_id][blockRow + 1]++;
            }

            // Cumulative sum to get row pointers
            for (int i = 0; i < num_block_rows; i++) {
                h_FR_bsrRowPtr_[gpu_id][i + 1] += h_FR_bsrRowPtr_[gpu_id][i];
            }

            // GETTTING bsrColInd_
            const int nnzb = num_FR_per_gpu_[gpu_id];
            h_FR_bsrColInd_[gpu_id] = new int[nnzb];

            // Fill colInd array
            int col_pointer = 0;
            for (int idx(0); idx<num_FR_per_gpu_[gpu_id]; idx++){
                int idx_offset = idx + offsets_FR_gpu_[gpu_id];
                int i = FR_std_orderedIndices[idx_offset];
                int j0 = hr->pattern_.FRB_pattern(2, i);
                int blockCol = j0 / leaf_size;
                h_FR_bsrColInd_[gpu_id][col_pointer] = blockCol;
                col_pointer++;
            }

            #ifdef DEBUG
            // Print it to check 
            std::cout << "[GPU "<< gpu_id << "] h_FR_bsrRowPtr_ = [";
            for (int i = 0; i < num_block_rows+1; i++) std::cout << h_FR_bsrRowPtr_[gpu_id][i] << ", ";
            std::cout << "]\n";
            std::cout << "[GPU "<< gpu_id << "] h_FR_bsrColInd_ = [";
            for (int i = 0; i < nnzb; i++) std::cout << h_FR_bsrColInd_[gpu_id][i] << ", ";
            std::cout << "]\n";
            #endif

            if (this->verbose_){
                size_t total_to_allocate_metadata = 0;
                total_to_allocate_metadata += (num_block_rows+1)*sizeof(int) + nnzb*sizeof(int);
                for (int block_size : LR_std_sizes_per_gpu_[gpu_id]){
                    int num_blocks = LR_std_indices_[block_size].size();
                    total_to_allocate_metadata += num_blocks * sizeof(T*) * 5;
                }

                std::cout << "[GPU "<< gpu_id << "] Allocating " << formatBytes(total_to_allocate_metadata) << " on the GPU for hierarchical matrix metadata" << std::endl;

            }

            // Then copy it to device
            CHECK_CUDA_ERROR(cudaMalloc(&d_FR_bsrRowPtr_[gpu_id], (num_block_rows+1)*sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_FR_bsrRowPtr_[gpu_id], h_FR_bsrRowPtr_[gpu_id], (num_block_rows+1)*sizeof(int), cudaMemcpyHostToDevice));

            CHECK_CUDA_ERROR(cudaMalloc(&d_FR_bsrColInd_[gpu_id], nnzb*sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_FR_bsrColInd_[gpu_id], h_FR_bsrColInd_[gpu_id], nnzb*sizeof(int), cudaMemcpyHostToDevice));

        #else
            // Initiating the vectors of pointers for Cublas support of FR blocks

            int num_blocks = num_FR_per_gpu_[gpu_id];

            T** h_FR_data_pointers = new T*[num_blocks];
            T** h_FR_x_pointers = new T*[num_blocks];
            T** h_FR_y_partial_pointers = new T*[num_blocks];

            // On the first GPU : we also get the results of the non std blocks
            int num_blocks_tot = num_blocks;
            if (gpu_id == 0) num_blocks_tot += FR_non_std_indices.size();

            int* h_FR_y_partial_src_indices = new int[num_blocks_tot];
            int* h_FR_y_partial_dest_indices = new int[num_blocks_tot];
            int* h_FR_y_partial_lengths = new int[num_blocks_tot];

            int fr_block_size = leaf_size*leaf_size * dim_dof*dim_dof;

            // Loop on blocks
            for (int idx(0); idx<num_blocks; idx++){
                int idx_offset = idx + offsets_FR_gpu_[gpu_id];
                int i = FR_std_orderedIndices[idx_offset];

                int i0 = hr->pattern_.FRB_pattern(1, i);
                int j0 = hr->pattern_.FRB_pattern(2, i);

                h_FR_data_pointers[idx] = d_FR_data_[gpu_id] + idx*fr_block_size;

                // x vector : dep on position
                h_FR_x_pointers[idx] = d_x_[gpu_id] + j0*dim_dof;

                // y_partial and tmp vectors : dep of offset
                h_FR_y_partial_pointers[idx] = d_y_partial_FR_[gpu_id] + idx*leaf_size*dim_dof;

                // Metadata for summing partial results
                h_FR_y_partial_src_indices[idx] = idx*leaf_size*dim_dof;
                h_FR_y_partial_dest_indices[idx] = i0*dim_dof;
                h_FR_y_partial_lengths[idx] = leaf_size*dim_dof;
            }

            // Copy on device

            // pointers
            T** device_ptr_data;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_data, num_blocks * sizeof(T*)));
            d_FR_data_pointers_[gpu_id] = device_ptr_data;
            CHECK_CUDA_ERROR(cudaMemcpy(d_FR_data_pointers_[gpu_id], h_FR_data_pointers, num_blocks * sizeof(T*), cudaMemcpyHostToDevice));

            T** device_ptr_x;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_x, num_blocks * sizeof(T*)));
            d_FR_x_pointers_[gpu_id] = device_ptr_x;
            CHECK_CUDA_ERROR(cudaMemcpy(d_FR_x_pointers_[gpu_id], h_FR_x_pointers, num_blocks * sizeof(T*), cudaMemcpyHostToDevice));

            T** device_ptr_y;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_y, num_blocks * sizeof(T*)));
            d_FR_y_partial_pointers_[gpu_id] = device_ptr_y;
            CHECK_CUDA_ERROR(cudaMemcpy(d_FR_y_partial_pointers_[gpu_id], h_FR_y_partial_pointers, num_blocks * sizeof(T*), cudaMemcpyHostToDevice));

            // deallocate on host
            delete[] h_FR_data_pointers;
            delete[] h_FR_x_pointers;
            delete[] h_FR_y_partial_pointers;


        #endif 

        // And then we set up the array of ponters for the non std FR blocks
        if (gpu_id == 0){

            // We first allocate all the arrays of pointers on host
            std::unordered_map<std::pair<int, int>, T**, pair_hash> h_FR_nonStd_data_pointers;
            std::unordered_map<std::pair<int, int>, T**, pair_hash> h_FR_nonStd_x_pointers;
            std::unordered_map<std::pair<int, int>, T**, pair_hash> h_FR_nonStd_y_partial_pointers;
            for (auto& [sizes,count] : num_FR_nonstd_blocks_per_size_){
                T** h_data_pointers_tmp = new T*[count];
                T** h_x_pointers_tmp = new T*[count];
                T** h_y_partial_pointers_tmp = new T*[count];
                h_FR_nonStd_data_pointers[sizes] = h_data_pointers_tmp;
                h_FR_nonStd_x_pointers[sizes] = h_x_pointers_tmp;
                h_FR_nonStd_y_partial_pointers[sizes] = h_y_partial_pointers_tmp;
            }

            // Then we set there values by looping on the blocks
            std::unordered_map<std::pair<int, int>, int, pair_hash> block_counter;
            for (auto& [sizes,count] : num_FR_nonstd_blocks_per_size_){
                block_counter[sizes] = 0;
            }

            // We offset the partial result by the space takn by the std sized FR blocks
            T* d_y_partial_FR_nonStd = d_y_partial_FR_[0] + num_FR_per_gpu_[0] * leaf_size * dim_dof;

            for (int idx(0); idx<FR_non_std_indices.size(); idx++){

                int i = FR_non_std_indices[idx];
        
                il::int_t i0 = hr->pattern_.FRB_pattern(1, i);
                il::int_t j0 = hr->pattern_.FRB_pattern(2, i);
                il::int_t iend = hr->pattern_.FRB_pattern(3, i);
                il::int_t jend = hr->pattern_.FRB_pattern(4, i);
        
                int data_size = (iend-i0)*(jend-j0)* dim_dof*dim_dof;
        
                std::pair<int, int> sizes = {iend-i0, jend-j0};

                int k = block_counter[sizes];

                h_FR_nonStd_data_pointers[sizes][k] = d_FR_non_std_data_ + FR_non_std_offsets_data_[sizes] + k*data_size;
                h_FR_nonStd_y_partial_pointers[sizes][k] = d_y_partial_FR_nonStd + FR_non_std_offsets_y_[sizes] + k*(iend-i0)*dim_dof;

                // x : depends on position of the block
                h_FR_nonStd_x_pointers[sizes][k] = d_x_[0] + j0*dim_dof;

                // We also add to the partial result summation arrays
                h_FR_y_partial_src_indices[num_FR_per_gpu_[0] + idx] = num_FR_per_gpu_[0]*leaf_size*dim_dof + FR_non_std_offsets_y_[sizes] + k*(iend-i0)*dim_dof;
                h_FR_y_partial_dest_indices[num_FR_per_gpu_[0] + idx] = i0*dim_dof;
                h_FR_y_partial_lengths[num_FR_per_gpu_[0] + idx] = (iend-i0)*dim_dof;

                block_counter[sizes]++;
            }

            // Finally we copy it to the device memory
            for (auto& [sizes,num_blocks] : num_FR_nonstd_blocks_per_size_){

                T** device_ptr_data;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_data, num_blocks * sizeof(T*)));
                d_FR_nonStd_data_pointers_[sizes] = device_ptr_data;
                CHECK_CUDA_ERROR(cudaMemcpy(d_FR_nonStd_data_pointers_[sizes], h_FR_nonStd_data_pointers[sizes], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));
    
                T** device_ptr_y;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_y, num_blocks * sizeof(T*)));
                d_FR_nonStd_y_partial_pointers_[sizes] = device_ptr_y;
                CHECK_CUDA_ERROR(cudaMemcpy(d_FR_nonStd_y_partial_pointers_[sizes], h_FR_nonStd_y_partial_pointers[sizes], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));

                T** device_ptr_x;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_x, num_blocks * sizeof(T*)));
                d_FR_nonStd_x_pointers_[sizes] = device_ptr_x;
                CHECK_CUDA_ERROR(cudaMemcpy(d_FR_nonStd_x_pointers_[sizes], h_FR_nonStd_x_pointers[sizes], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));
            }

            // And we deallocate on host
            for (auto& [sizes,num_blocks] : num_FR_nonstd_blocks_per_size_){
                delete[] h_FR_nonStd_data_pointers[sizes];
                delete[] h_FR_nonStd_x_pointers[sizes];
                delete[] h_FR_nonStd_y_partial_pointers[sizes];
            }

        }

        // And we copy the arrays used for summation of partial results
        CHECK_CUDA_ERROR(cudaMalloc(&d_FR_y_partial_src_indices_[gpu_id], num_blocks_tot*sizeof(int)));
        CHECK_CUDA_ERROR(cudaMalloc(&d_FR_y_partial_dest_indices_[gpu_id], num_blocks_tot*sizeof(int)));
        CHECK_CUDA_ERROR(cudaMalloc(&d_FR_y_partial_lengths_[gpu_id], num_blocks_tot*sizeof(int)));

        CHECK_CUDA_ERROR(cudaMemcpy(d_FR_y_partial_src_indices_[gpu_id], h_FR_y_partial_src_indices, num_blocks_tot*sizeof(int), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_FR_y_partial_dest_indices_[gpu_id], h_FR_y_partial_dest_indices, num_blocks_tot*sizeof(int), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_FR_y_partial_lengths_[gpu_id], h_FR_y_partial_lengths, num_blocks_tot*sizeof(int), cudaMemcpyHostToDevice));

        delete[] h_FR_y_partial_src_indices;
        delete[] h_FR_y_partial_dest_indices;
        delete[] h_FR_y_partial_lengths;


        // ----------------------------
        // COPYING FULL RANKS BLOCKS - done
        // ----------------------------

        // ----------------------------
        // COPYING LOW RANKS BLOCKS
        // ----------------------------

        // We copy all buffers ;
        // 1 = std size
        for (int block_size : LR_std_sizes_per_gpu_[gpu_id]){
            int buffer_size = LR_standard_size_data_buffer_sizes_[block_size];

            // Allocation
            T* d_A_tmp;
            CHECK_CUDA_ERROR(cudaMalloc(&d_A_tmp, buffer_size*sizeof(T)));
            d_LR_A_data_[gpu_id][block_size] = d_A_tmp;
            T* d_B_tmp;
            CHECK_CUDA_ERROR(cudaMalloc(&d_B_tmp, buffer_size*sizeof(T)));
            d_LR_B_data_[gpu_id][block_size] = d_B_tmp;

            // Copy
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_A_data_[gpu_id][block_size], LR_standard_size_A_data_[block_size], buffer_size*sizeof(T), cudaMemcpyHostToDevice));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_B_data_[gpu_id][block_size], LR_standard_size_B_data_[block_size], buffer_size*sizeof(T), cudaMemcpyHostToDevice));

            #ifdef WRITE_LR_DATA_NPY
            std::stringstream ss1;
            ss1 << "gpu" << gpu_id << "_" << num_gpus_ << "_size" << block_size<< "_A_data.npy";
            std::string filename = ss1.str();
            cnpy::npy_save(filename, LR_standard_size_A_data_[block_size], {static_cast<size_t>(buffer_size)}, "w");

            std::stringstream ss2;
            ss2 << "gpu" << gpu_id << "_" << num_gpus_ << "_size" << block_size<< "_B_data.npy";
            filename = ss2.str();
            cnpy::npy_save(filename, LR_standard_size_B_data_[block_size], {static_cast<size_t>(buffer_size)}, "w");
            #endif
        }

        // 2 = non std size
        if (gpu_id == 0){

            CHECK_CUDA_ERROR(cudaMalloc(&d_LR_non_std_A_data_, LR_non_standard_size_data_A_buffer_size_*sizeof(T)));
            CHECK_CUDA_ERROR(cudaMalloc(&d_LR_non_std_B_data_, LR_non_standard_size_data_B_buffer_size_*sizeof(T)));

            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_non_std_A_data_, LR_non_standard_size_A_data_, LR_non_standard_size_data_A_buffer_size_*sizeof(T), cudaMemcpyHostToDevice));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_non_std_B_data_, LR_non_standard_size_B_data_, LR_non_standard_size_data_B_buffer_size_*sizeof(T), cudaMemcpyHostToDevice));
        }

        // We set the array of pointers needed to perform the batched operations

        size_t y_partial_LR_ptr_counter = 0;
        size_t tmp_ptr_counter = 0;

        // Offsets for summing partial results

        int num_lr_blocks_this_gpu = num_LR_per_gpu_[gpu_id];
        if (gpu_id == 0) num_lr_blocks_this_gpu += LR_non_std_indices_.size();

        int* h_LR_y_partial_src_indices = new int[num_lr_blocks_this_gpu];
        int* h_LR_y_partial_dest_indices = new int[num_lr_blocks_this_gpu];
        int* h_LR_y_partial_lengths = new int[num_lr_blocks_this_gpu];
        int lr_block_i = 0;

        for (int block_size : LR_std_sizes_per_gpu_[gpu_id]){

            int num_blocks = LR_std_indices_[block_size].size();
            int block_data_size = block_size*fixed_rank * dim_dof*dim_dof;

            h_LR_A_data_pointers_[gpu_id][block_size] = new T*[num_blocks];
            h_LR_B_data_pointers_[gpu_id][block_size] = new T*[num_blocks];
            h_LR_x_pointers_[gpu_id][block_size] = new T*[num_blocks];
            h_LR_y_pointers_[gpu_id][block_size] = new T*[num_blocks];
            h_LR_tmp_pointers_[gpu_id][block_size] = new T*[num_blocks];

            for (int i=0; i<num_blocks; i++){

                // Data = linear (non dep on position)
                h_LR_A_data_pointers_[gpu_id][block_size][i] = d_LR_A_data_[gpu_id][block_size] + i*block_data_size;
                h_LR_B_data_pointers_[gpu_id][block_size][i] = d_LR_B_data_[gpu_id][block_size] + i*block_data_size;

                // Getting block position
                int block_id = LR_std_indices_[block_size][i];
                int i0 = hr->pattern_.LRB_pattern(1, block_id);
                int j0 = hr->pattern_.LRB_pattern(2, block_id);

                // x vector : dep on position
                h_LR_x_pointers_[gpu_id][block_size][i] = d_x_[gpu_id] + j0*dim_dof;

                // y_partial and tmp vectors : dep of offset
                h_LR_y_pointers_[gpu_id][block_size][i] = d_y_partial_LR_[gpu_id] + y_partial_LR_ptr_counter;
                h_LR_tmp_pointers_[gpu_id][block_size][i] = d_tmp_[gpu_id] + tmp_ptr_counter;

                // Metadata for summing partial results
                h_LR_y_partial_src_indices[lr_block_i] = y_partial_LR_ptr_counter;
                h_LR_y_partial_dest_indices[lr_block_i] = i0*dim_dof;
                h_LR_y_partial_lengths[lr_block_i] = block_size*dim_dof;

                // Incrementing the offsets
                y_partial_LR_ptr_counter += block_size*dim_dof;
                tmp_ptr_counter += this->fixed_rank_*dim_dof;
                lr_block_i++;
            }

            // Now we copy it to device
            T** d_LR_A_data_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_A_data_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_A_data_tmp, h_LR_A_data_pointers_[gpu_id][block_size], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));
            d_LR_A_data_pointers_[gpu_id][block_size] = d_LR_A_data_tmp;

            T** d_LR_B_data_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_B_data_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_B_data_tmp, h_LR_B_data_pointers_[gpu_id][block_size], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));
            d_LR_B_data_pointers_[gpu_id][block_size] = d_LR_B_data_tmp;

            T** d_LR_x_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_x_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_x_tmp, h_LR_x_pointers_[gpu_id][block_size], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));
            d_LR_x_pointers_[gpu_id][block_size] = d_LR_x_tmp;

            T** d_LR_y_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_y_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_y_tmp, h_LR_y_pointers_[gpu_id][block_size], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
            d_LR_y_pointers_[gpu_id][block_size] = d_LR_y_tmp;

            T** d_LR_tmp_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_tmp_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_tmp_tmp, h_LR_tmp_pointers_[gpu_id][block_size], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
            d_LR_tmp_pointers_[gpu_id][block_size] = d_LR_tmp_tmp;
        }

        
        // On the first GPU : we also set up the pointers for the non std size blocks
        if (gpu_id == 0){

            // We first allocate all the arrays of pointers on host
            std::unordered_map<std::pair<int, int>, T**, pair_hash> h_LR_nonStd_A_data_pointers;
            std::unordered_map<std::pair<int, int>, T**, pair_hash> h_LR_nonStd_B_data_pointers;
            std::unordered_map<std::pair<int, int>, T**, pair_hash> h_LR_nonStd_x_pointers;
            std::unordered_map<std::pair<int, int>, T**, pair_hash> h_LR_nonStd_y_partial_pointers;
            std::unordered_map<std::pair<int, int>, T**, pair_hash> h_LR_nonStd_tmp_pointers;
            for (auto& [sizes,count] : num_LR_nonstd_blocks_per_size_){
                T** h_A_data_pointers_tmp = new T*[count];
                T** h_B_data_pointers_tmp = new T*[count];
                T** h_x_pointers_tmp = new T*[count];
                T** h_y_partial_pointers_tmp = new T*[count];
                T** h_tmp_pointers_tmp = new T*[count];
                h_LR_nonStd_A_data_pointers[sizes] = h_A_data_pointers_tmp;
                h_LR_nonStd_B_data_pointers[sizes] = h_B_data_pointers_tmp;
                h_LR_nonStd_x_pointers[sizes] = h_x_pointers_tmp;
                h_LR_nonStd_y_partial_pointers[sizes] = h_y_partial_pointers_tmp;
                h_LR_nonStd_tmp_pointers[sizes] = h_tmp_pointers_tmp;
            }

            // Then we set there values by looping on the blocks
            std::unordered_map<std::pair<int, int>, int, pair_hash> block_counter;
            for (auto& [sizes,count] : num_LR_nonstd_blocks_per_size_){
                block_counter[sizes] = 0;
            }

            // We offset the partial result by the space takn by the std sized FR blocks
            T* d_y_partial_LR_nonStd = d_y_partial_LR_[0];
            T* d_tmp_LR_nonStd = d_tmp_[0];
            for (int block_size : LR_std_sizes_per_gpu_[gpu_id]){
                int num_blocks = LR_std_indices_[block_size].size();
                d_y_partial_LR_nonStd += num_blocks * block_size * dim_dof;
                d_tmp_LR_nonStd += num_blocks * fixed_rank * dim_dof;
            }

            for (int idx(0); idx<LR_non_std_indices_.size(); idx++){

                int i = LR_non_std_indices_[idx];
        
                il::int_t i0 = hr->pattern_.LRB_pattern(1, i);
                il::int_t j0 = hr->pattern_.LRB_pattern(2, i);
                il::int_t iend = hr->pattern_.LRB_pattern(3, i);
                il::int_t jend = hr->pattern_.LRB_pattern(4, i);
        
                int data_size_A = (iend-i0)*fixed_rank * dim_dof*dim_dof;
                int data_size_B = (jend-j0)*fixed_rank * dim_dof*dim_dof;

                std::pair<int, int> sizes = {iend-i0, jend-j0};

                int k = block_counter[sizes];

                h_LR_nonStd_A_data_pointers[sizes][k] = d_LR_non_std_A_data_ + LR_non_std_A_offsets_data_[sizes] + k*data_size_A;
                h_LR_nonStd_B_data_pointers[sizes][k] = d_LR_non_std_B_data_ + LR_non_std_B_offsets_data_[sizes] + k*data_size_B;


                h_LR_nonStd_y_partial_pointers[sizes][k] = d_y_partial_LR_nonStd + LR_non_std_offsets_y_[sizes] + k*(iend-i0)*dim_dof;
                h_LR_nonStd_tmp_pointers[sizes][k] = d_tmp_LR_nonStd + LR_non_std_offsets_tmp_[sizes] + k*fixed_rank*dim_dof;

                // x : depends on position of the block
                h_LR_nonStd_x_pointers[sizes][k] = d_x_[0] + j0*dim_dof;

                // We also add to the partial result summation arrays
                h_LR_y_partial_src_indices[num_LR_per_gpu_[0] + idx] = (d_y_partial_LR_nonStd - d_y_partial_LR_[0]) + LR_non_std_offsets_y_[sizes] + k*(iend-i0)*dim_dof;
                h_LR_y_partial_dest_indices[num_LR_per_gpu_[0] + idx] = i0*dim_dof;
                h_LR_y_partial_lengths[num_LR_per_gpu_[0] + idx] = (iend-i0)*dim_dof;

                block_counter[sizes]++;
            }

            // Finally we copy it to the device memory
            for (auto& [sizes,num_blocks] : num_LR_nonstd_blocks_per_size_){

                T** device_ptr_data_A;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_data_A, num_blocks * sizeof(T*)));
                d_LR_nonStd_A_data_pointers_[sizes] = device_ptr_data_A;
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_nonStd_A_data_pointers_[sizes], h_LR_nonStd_A_data_pointers[sizes], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));

                T** device_ptr_data_B;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_data_B, num_blocks * sizeof(T*)));
                d_LR_nonStd_B_data_pointers_[sizes] = device_ptr_data_B;
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_nonStd_B_data_pointers_[sizes], h_LR_nonStd_B_data_pointers[sizes], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));

                T** device_ptr_x;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_x, num_blocks * sizeof(T*)));
                d_LR_nonStd_x_pointers_[sizes] = device_ptr_x;
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_nonStd_x_pointers_[sizes], h_LR_nonStd_x_pointers[sizes], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));

                T** device_ptr_y;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_y, num_blocks * sizeof(T*)));
                d_LR_nonStd_y_partial_pointers_[sizes] = device_ptr_y;
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_nonStd_y_partial_pointers_[sizes], h_LR_nonStd_y_partial_pointers[sizes], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));
                
                T** device_ptr_tmp;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&device_ptr_tmp, num_blocks * sizeof(T*)));
                d_LR_nonStd_tmp_pointers_[sizes] = device_ptr_tmp;
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_nonStd_tmp_pointers_[sizes], h_LR_nonStd_tmp_pointers[sizes], num_blocks * sizeof(T*), cudaMemcpyHostToDevice));

            }

            // And we deallocate on host
            for (auto& [sizes,num_blocks] : num_LR_nonstd_blocks_per_size_){
                delete[] h_LR_nonStd_A_data_pointers[sizes];
                delete[] h_LR_nonStd_B_data_pointers[sizes];
                delete[] h_LR_nonStd_x_pointers[sizes];
                delete[] h_LR_nonStd_y_partial_pointers[sizes];
                delete[] h_LR_nonStd_tmp_pointers[sizes];
            }

        }


        // Then we copy the array of offsets needed for gathering the partial results
        CHECK_CUDA_ERROR(cudaMalloc(&d_LR_y_partial_src_indices_[gpu_id], num_lr_blocks_this_gpu*sizeof(int)));
        CHECK_CUDA_ERROR(cudaMalloc(&d_LR_y_partial_dest_indices_[gpu_id], num_lr_blocks_this_gpu*sizeof(int)));
        CHECK_CUDA_ERROR(cudaMalloc(&d_LR_y_partial_lengths_[gpu_id], num_lr_blocks_this_gpu*sizeof(int)));

        CHECK_CUDA_ERROR(cudaMemcpy(d_LR_y_partial_src_indices_[gpu_id], h_LR_y_partial_src_indices, num_lr_blocks_this_gpu*sizeof(int), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_LR_y_partial_dest_indices_[gpu_id], h_LR_y_partial_dest_indices, num_lr_blocks_this_gpu*sizeof(int), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_LR_y_partial_lengths_[gpu_id], h_LR_y_partial_lengths, num_lr_blocks_this_gpu*sizeof(int), cudaMemcpyHostToDevice));

        delete[] h_LR_y_partial_src_indices;
        delete[] h_LR_y_partial_dest_indices;
        delete[] h_LR_y_partial_lengths;

        // Add explicit CUDA error check
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            printf("CUDA error at the end of copy2device: %s\n", cudaGetErrorString(err));
        }

        // ----------------------------
        // COPYING LOW RANKS BLOCKS - done
        // ----------------------------

    } // Loop on GPUs

    cudaDeviceSynchronize();



    if (this->verbose_) std::cout << "--------------------" << std::endl;
}

template <typename T>
void HmatCuda<T>::deallocateOnHost(bool deallocate_FR_host){

    // First we nullify the pointers of the Array2d objects
    // to avoid double freeing
    if (deallocate_FR_host){
        for (auto& array2d : this->full_rank_blocks_){
            if (array2d) {
                array2d->nullifyData(); 
            }
        }

        // Deallocate FR buffers
        if (FR_standard_size_data) {
            delete[] FR_standard_size_data;
            FR_standard_size_data = nullptr;
        }
        if (FR_non_standard_size_data) {
            delete[] FR_non_standard_size_data;
            FR_non_standard_size_data = nullptr;
        }

        FR_deallocated_on_host_ = true;
    }


    for (auto& low_rank : this->low_rank_blocks_){
        if (low_rank){
            auto& A = low_rank->A;
            auto& B = low_rank->B;
            A.nullifyData(); 
            B.nullifyData(); 
        }
    }

    // Deallocate non std LR buffers
    delete[] LR_non_standard_size_A_data_;
    delete[] LR_non_standard_size_B_data_;

    // Deallocate std LR buffers
    for (auto& pair : LR_std_indices_){
        int block_size = pair.first;
        delete[] LR_standard_size_A_data_[block_size];
        delete[] LR_standard_size_B_data_[block_size];
    }
    LR_standard_size_A_data_.clear();
    LR_standard_size_B_data_.clear();

    // Deallocate associated pointers array
    for (int gpu_id(0); gpu_id<num_gpus_; gpu_id++){
        for (int block_size : LR_std_sizes_per_gpu_[gpu_id]){
            delete[] h_LR_A_data_pointers_[gpu_id][block_size];
            delete[] h_LR_B_data_pointers_[gpu_id][block_size];
            delete[] h_LR_x_pointers_[gpu_id][block_size];
            delete[] h_LR_y_pointers_[gpu_id][block_size];
            delete[] h_LR_tmp_pointers_[gpu_id][block_size];
        }
        h_LR_A_data_pointers_[gpu_id].clear();
        h_LR_B_data_pointers_[gpu_id].clear();
        h_LR_x_pointers_[gpu_id].clear();
        h_LR_y_pointers_[gpu_id].clear();
        h_LR_tmp_pointers_[gpu_id].clear();
    
        delete[] h_FR_bsrRowPtr_[gpu_id];
        delete[] h_FR_bsrColInd_[gpu_id];    
    }
}

template <typename T>
void HmatCuda<T>::deallocateOnDevice(){

    for (int gpu_id(0); gpu_id<num_gpus_; gpu_id++){

        // Destroy the handles
        CHECK_CUBLAS_ERROR(cublasDestroy(cublas_handle_[gpu_id]));

        // Deallocate
        CHECK_CUDA_ERROR(cudaFree(d_x_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_y_[gpu_id]));

        CHECK_CUDA_ERROR(cudaFree(d_y_partial_LR_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_data_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_bsrColInd_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_bsrRowPtr_[gpu_id]));

        CHECK_CUDA_ERROR(cudaFree(d_FR_data_pointers_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_x_pointers_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_y_partial_pointers_[gpu_id]));

        CHECK_CUDA_ERROR(cudaFree(d_FR_y_partial_src_indices_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_y_partial_dest_indices_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_y_partial_lengths_[gpu_id]));

        CHECK_CUDA_ERROR(cudaFree(d_LR_y_partial_src_indices_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_LR_y_partial_dest_indices_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_LR_y_partial_lengths_[gpu_id]));

        for (auto& pair : LR_std_indices_){
            int block_size = pair.first;
            CHECK_CUDA_ERROR(cudaFree(d_LR_A_data_[gpu_id][block_size]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_B_data_[gpu_id][block_size]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_A_data_pointers_[gpu_id][block_size]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_B_data_pointers_[gpu_id][block_size]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_x_pointers_[gpu_id][block_size]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_y_pointers_[gpu_id][block_size]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_tmp_pointers_[gpu_id][block_size]));
        }

        if (gpu_id == 0){
            for (auto& [sizes,num_blocks] : num_FR_nonstd_blocks_per_size_){
                CHECK_CUDA_ERROR(cudaFree(d_FR_nonStd_data_pointers_[sizes]));
                CHECK_CUDA_ERROR(cudaFree(d_FR_nonStd_x_pointers_[sizes]));
                CHECK_CUDA_ERROR(cudaFree(d_FR_nonStd_y_partial_pointers_[sizes]));
            }
    
            for (auto& [sizes,num_blocks] : num_LR_nonstd_blocks_per_size_){
                CHECK_CUDA_ERROR(cudaFree(d_LR_nonStd_A_data_pointers_[sizes]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_nonStd_B_data_pointers_[sizes]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_nonStd_x_pointers_[sizes]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_nonStd_y_partial_pointers_[sizes]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_nonStd_tmp_pointers_[sizes]));
            }
        }

        // Destroy the Magma queues
        #ifdef USE_MAGMA
        for (int i = 1; i < num_streams_; i++) magma_queue_destroy(magma_queues_[gpu_id][i-1]);
        #endif

        // Destroy the cuda streams
        for (int i = 0; i < num_streams_; i++) cudaStreamDestroy(cuda_streams_[gpu_id][i]);
    }

}

template <typename T>
HmatCuda<T>::~HmatCuda(){
    // Clear device memory
    deallocateOnDevice();

    // Deallocate FR on host if not done yet 
    if (!FR_deallocated_on_host_){
        for (auto& array2d : this->full_rank_blocks_){
            if (array2d) {
                array2d->nullifyData(); 
            }
        }
    
        // Deallocate FR buffers
        if (FR_standard_size_data) {
            delete[] FR_standard_size_data;
            FR_standard_size_data = nullptr;
        }
        if (FR_non_standard_size_data) {
            delete[] FR_non_standard_size_data;
            FR_non_standard_size_data = nullptr;
        }
    }
}


template <typename T>
void HmatCuda<T>::buildCuda(const bigwham::MatrixGenerator<T> & matrix_gen,const double epsilon) {
    this->dof_dimension_ = matrix_gen.blockSize();
    this->size_[0] = matrix_gen.size(0);
    this->size_[1] = matrix_gen.size(1);

    // Then we set up the memory for the vectors
    if (this->size_[0] != this->size_[1]){
        std::cerr << "CUDA only implemented for square Hmat for now\n";
        std::cerr << "Hmat size = " << this->size_[0] << " x " << this->size_[1] << std::endl;
        std::abort();
    }

#ifdef TIMING
    struct timespec start, end;
    double duration;
    clock_gettime(CLOCK_MONOTONIC, &start);
#endif // TIMING 

    buildFRCuda(matrix_gen);

#ifdef TIMING
    clock_gettime(CLOCK_MONOTONIC, &end);
    duration = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;    
    std::cout << "[Timing] building the FR blocks = " << duration*1000 << "ms\n";
    clock_gettime(CLOCK_MONOTONIC, &start);
#endif // TIMING 

    switch (matrix_gen.blockSize()) {
        case 1:
            buildLRCuda<1>(matrix_gen, epsilon);
            break;
        case 2:
            buildLRCuda<2>(matrix_gen, epsilon);
            break;
        case 3:
            buildLRCuda<3>(matrix_gen, epsilon);
            break;
    }

#ifdef TIMING
    clock_gettime(CLOCK_MONOTONIC, &end);
    duration = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;    
    std::cout << "[Timing] building the LR blocks = " << duration*1000 << "ms\n";
#endif // TIMING 

    this->isBuilt_ = this->isBuilt_FR_ && this->isBuilt_LR_;
}

template <typename T>
void HmatCuda<T>::buildFRCuda(const bigwham::MatrixGenerator<T> & matrix_gen){

    if (this->verbose_){
        std::cout << "Loop on full blocks construction  \n";
    }

    /*
    Goal = create the blocks at the right place + add their pointers 
    to full_rank_blocks_*/
    auto hr = this->hr_;
    this->full_rank_blocks_.resize(hr->pattern_.n_FRB);

    const int dim = matrix_gen.blockSize();
    int leaf_size = this->hr_->leaf_size;

    // Test to memset the buffers
    std::memset(FR_standard_size_data, 0, FR_standard_size_data_buffer_size*sizeof(T));
    std::memset(FR_non_standard_size_data, 0, FR_non_standard_size_data_buffer_size*sizeof(T));
        
    // First the standard size blocks
#pragma omp parallel for schedule(guided)
    for (int idx=0; idx<FR_std_orderedIndices.size(); idx++){

        int i = FR_std_orderedIndices[idx];

        il::int_t i0 = hr->pattern_.FRB_pattern(1, i);
        il::int_t j0 = hr->pattern_.FRB_pattern(2, i);
        il::int_t iend = hr->pattern_.FRB_pattern(3, i);
        il::int_t jend = hr->pattern_.FRB_pattern(4, i);

        const il::int_t ni = matrix_gen.blockSize() * (iend - i0);
        const il::int_t nj = matrix_gen.blockSize() * (jend - j0);

        // We create the Array2D
        std::shared_ptr<il::Array2D<T>> a = std::make_shared<il::Array2D<T>>(ni, nj);

        // We deallocate the data it has just been allocated
        a->deallocateData();

        // We set its data to be where we want 
        int data_size = (iend-i0)*(jend-j0)* dim*dim;
        int offset = data_size*idx;
        if (offset + data_size > FR_standard_size_data_buffer_size)
            std::cerr << "FR std : setting data outside the allocated buffer (allocated = " << FR_standard_size_data_buffer_size << ", offset = " << offset << ", size = "<< data_size <<")\n";

            a->setData(&FR_standard_size_data[offset]);

        // Finnaly we populate it 
        matrix_gen.set(i0, j0, il::io, a->Edit());
        this->full_rank_blocks_[i] = a;
    }

    // For the non standard blocks we first need a list of offsets

    // We first compute the offset per size 
    int offset_total_data = 0;
    int offset_total_y = 0;
    for (auto& [sizes,count] : num_FR_nonstd_blocks_per_size_){
        FR_non_std_offsets_data_[sizes] = offset_total_data;
        offset_total_data += count * sizes.first*sizes.second * dim*dim;
        FR_non_std_offsets_y_[sizes] = offset_total_y;
        offset_total_y += count * sizes.first*dim;
    }

    // Then we loop on blocks, with a block counter for each size
    std::unordered_map<std::pair<int, int>, int, pair_hash> block_counter;
    for (auto& [sizes,count] : num_FR_nonstd_blocks_per_size_){
        block_counter[sizes] = 0;
    }

    int* offsets_non_standard = new int[FR_non_std_indices.size()];
    for (int idx(0); idx<FR_non_std_indices.size(); idx++){

        int i = FR_non_std_indices[idx];

        il::int_t i0 = hr->pattern_.FRB_pattern(1, i);
        il::int_t j0 = hr->pattern_.FRB_pattern(2, i);
        il::int_t iend = hr->pattern_.FRB_pattern(3, i);
        il::int_t jend = hr->pattern_.FRB_pattern(4, i);

        int data_size = (iend-i0)*(jend-j0)* dim*dim;

        std::pair<int, int> sizes = {iend-i0, jend-j0};
        offsets_non_standard[idx] = FR_non_std_offsets_data_[sizes] + block_counter[sizes]*data_size;
        block_counter[sizes]++;

        // offsets_non_standard[idx] = offsets_non_standard[idx-1] + data_size;
    }

    // Then the non standard size blocks
#pragma omp parallel for schedule(guided)
    for (int idx=0; idx<FR_non_std_indices.size(); idx++){

        int i = FR_non_std_indices[idx];

        il::int_t i0 = hr->pattern_.FRB_pattern(1, i);
        il::int_t j0 = hr->pattern_.FRB_pattern(2, i);
        il::int_t iend = hr->pattern_.FRB_pattern(3, i);
        il::int_t jend = hr->pattern_.FRB_pattern(4, i);

        const il::int_t ni = matrix_gen.blockSize() * (iend - i0);
        const il::int_t nj = matrix_gen.blockSize() * (jend - j0);

        // We create the Array2D
        std::shared_ptr<il::Array2D<T>> a = std::make_shared<il::Array2D<T>>(ni, nj);

        // We deallocate the data it has just been allocated
        a->deallocateData();

        // We get the correct offset
        int data_size = (iend-i0)*(jend-j0)* dim*dim;
        if (offsets_non_standard[idx] + data_size > FR_non_standard_size_data_buffer_size)
            std::cerr << "FR non std : setting data outside the allocated buffer\n";

        a->setData(&FR_non_standard_size_data[offsets_non_standard[idx]]);

        // Finnaly we populate it 
        matrix_gen.set(i0, j0, il::io, a->Edit());
        this->full_rank_blocks_[i] = a;
    }

    this->isBuilt_FR_ = true;
    
}

template <typename T>
template <il::int_t dim>
void HmatCuda<T>::buildLRCuda(const bigwham::MatrixGenerator<T> & matrix_gen, const double epsilon) {

    // constructing the low rank blocks
    if (this->verbose_){
        std::cout << "Loop on low rank blocks construction\n";
    }

    const int dim_dof = matrix_gen.blockSize();
    const int fixed_rank = this->fixed_rank_;

    /*
    Goal = keep track of the offsets for each group of LR blocks
    */
    auto hr = this->hr_;
    this->low_rank_blocks_.resize(hr->pattern_.n_LRB);

    // We initialize the offsets in the memory buffers
    // For the non std blocks
    // 1. offset per size
    int offset_total_A = 0;
    int offset_total_B = 0;
    int offset_total_y = 0;
    int offset_total_tmp = 0;
    for (auto& [sizes,count] : num_LR_nonstd_blocks_per_size_){
        LR_non_std_A_offsets_data_[sizes] = offset_total_A;
        LR_non_std_B_offsets_data_[sizes] = offset_total_B;
        LR_non_std_offsets_y_[sizes] = offset_total_y;
        LR_non_std_offsets_tmp_[sizes] = offset_total_tmp;
        offset_total_A += count * sizes.first*fixed_rank * dim*dim;
        offset_total_B += count * sizes.second*fixed_rank * dim*dim;
        offset_total_y += count * sizes.first*dim;
        offset_total_tmp += count * fixed_rank*dim;
    }

    std::unordered_map<std::pair<int, int>, int, pair_hash> block_counter;
    for (auto& [sizes,count] : num_LR_nonstd_blocks_per_size_){
        block_counter[sizes] = 0;
    }

    // 2. loop on blocks to get right offset
    std::vector<int> offsets_A_non_std = std::vector<int>(LR_non_std_indices_.size(), 0);
    std::vector<int> offsets_B_non_std = std::vector<int>(LR_non_std_indices_.size(), 0);
    if (LR_non_std_indices_.size() > 0){

        for (int i=0; i<LR_non_std_indices_.size(); i++){

            int idx = LR_non_std_indices_[i];
            int i0 = hr->pattern_.LRB_pattern(1,idx);
            int j0 = hr->pattern_.LRB_pattern(2,idx);
            int iend = hr->pattern_.LRB_pattern(3,idx);
            int jend = hr->pattern_.LRB_pattern(4,idx);

            int data_size_A = (iend-i0)*fixed_rank * dim_dof*dim_dof;
            int data_size_B = (jend-j0)*fixed_rank * dim_dof*dim_dof;

            std::pair<int, int> sizes = {iend-i0, jend-j0};

            offsets_A_non_std[i] = LR_non_std_A_offsets_data_[sizes] + block_counter[sizes]*data_size_A;
            offsets_B_non_std[i] = LR_non_std_B_offsets_data_[sizes] + block_counter[sizes]*data_size_B;

            block_counter[sizes]++;
            // offsets_A_non_std[i+1] = offsets_A_non_std[i] + data_size_A;
            // offsets_B_non_std[i+1] = offsets_B_non_std[i] + data_size_B;
        }

    }


    // For the std blocks : one offset per size
    // We initiate the offsets
    std::unordered_map<int, std::vector<int>> offsets_per_size_std;
    for (auto& [block_size, indices] : LR_std_indices_){
        int num_blocks = indices.size();
        offsets_per_size_std[block_size] = std::vector<int>(num_blocks, 0);
    }
    // And we compute them
    for (auto& [block_size, indices] : LR_std_indices_){
        for (int i(1); i<indices.size(); i++){
            offsets_per_size_std[block_size][i] = offsets_per_size_std[block_size][i-1] + block_size*fixed_rank * dim_dof*dim_dof;
        }
    }

    // We construct the non std blocks
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < LR_non_std_indices_.size(); i++) {
        int idx = LR_non_std_indices_[i];

        il::int_t i0 = hr->pattern_.LRB_pattern(1, idx);
        il::int_t j0 = hr->pattern_.LRB_pattern(2, idx);
        il::int_t iend = hr->pattern_.LRB_pattern(3, idx);
        il::int_t jend = hr->pattern_.LRB_pattern(4, idx);
        il::Range range0{i0, iend};
        il::Range range1{j0, jend};

        // Here we cerate a LowRank struct that stores the two Array2D we want to copy
        auto lra = bigwham::adaptiveCrossApproximation<dim>(matrix_gen, range0, range1, epsilon, this->fixed_rank_);
        auto A_orig = lra->A;
        auto B_orig = lra->B;

        // We create a new one 
        auto lrb = std::make_unique<LowRank<T>>();
        auto &A = lrb->A;
        auto &B = lrb->B;

        // We resize it but deallocate its data
        A.Resize(A_orig.size(0), A_orig.size(1));
        B.Resize(B_orig.size(0), B_orig.size(1));
        A.deallocateData();
        B.deallocateData();

        #ifdef DEBUG
        std::cout << "Non std LR block " << i << " : A.shape = " << A_orig.size(0) << " x "<< A_orig.size(1) << ", B.shape = " << B_orig.size(0) << " x "<< B_orig.size(1) << std::endl;
        #endif

        // We set their memory 
        int offset_A = offsets_A_non_std[i];
        int offset_B = offsets_B_non_std[i];
        int data_size_A = (iend-i0)*fixed_rank * dim_dof*dim_dof;
        int data_size_B = (jend-j0)*fixed_rank * dim_dof*dim_dof;

        // Sanity check : ensure memory is allocated
        if (offset_A + data_size_A > LR_non_standard_size_data_A_buffer_size_){
            std::cerr << "Non std size LR blocks : setting A data outside the allocated buffer\n";
        }
        if (offset_B + data_size_B > LR_non_standard_size_data_B_buffer_size_){
            std::cerr << "Non std size LR blocks : setting B data outside the allocated buffer\n";
        }
               
        A.setData(&LR_non_standard_size_A_data_[offset_A]);
        B.setData(&LR_non_standard_size_B_data_[offset_B]);

        // We copy the data
        auto A_edit = A.Edit();
        auto B_edit = B.Edit();
        auto A_orig_view = A_orig.view();
        auto B_orig_view = B_orig.view();

        for (size_t i = 0; i < A_orig.size(0); i++) {
            for (size_t j = 0; j < A_orig.size(1); j++) {
                A_edit(i,j) = A_orig_view(i ,j);
            }
        }
        for (size_t i = 0; i < B_orig.size(0); i++) {
            for (size_t j = 0; j < B_orig.size(1); j++) {
                B_edit(i,j) = B_orig_view(i ,j);
            }
        }

        // We also copy the error on approximation
        lrb->error_on_approximation = lra->error_on_approximation;

        // Finally we move it
        this->low_rank_blocks_[idx] = std::move(lrb);
    }

    #ifdef WRITE_LR_DATA_NPY
    std::stringstream ss1;
    ss1 << "A_non_std_data.npy";
    std::string filename1 = ss1.str();
    cnpy::npy_save(filename1, LR_non_standard_size_A_data_, {static_cast<size_t>(LR_non_standard_size_data_A_buffer_size_)}, "w");

    std::stringstream ss2;
    ss2 << "B_non_std_data.npy";
    std::string filename2 = ss2.str();
    cnpy::npy_save(filename2, LR_non_standard_size_B_data_, {static_cast<size_t>(LR_non_standard_size_data_B_buffer_size_)}, "w");
    #endif

    // We construct the std blocks
    for (const auto& pair : LR_std_indices_) {
        const auto& block_size = pair.first;
        const auto& indices = pair.second;
        
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i < indices.size(); i++) {
            int idx = indices[i];

            il::int_t i0 = hr->pattern_.LRB_pattern(1, idx);
            il::int_t j0 = hr->pattern_.LRB_pattern(2, idx);
            il::int_t iend = hr->pattern_.LRB_pattern(3, idx);
            il::int_t jend = hr->pattern_.LRB_pattern(4, idx);
            il::Range range0{i0, iend};
            il::Range range1{j0, jend};
    
            int block_size = iend - i0;
    
            // Here we cerate a LowRank struct that stores the two Array2D we want to copy
            auto lra = bigwham::adaptiveCrossApproximation<dim>(matrix_gen, range0, range1, epsilon, this->fixed_rank_);
            auto A_orig = lra->A;
            auto B_orig = lra->B;
    
            // We create a new one 
            auto lrb = std::make_unique<LowRank<T>>();
            auto &A = lrb->A;
            auto &B = lrb->B;
    
            // We resize it but deallocate its data
            A.Resize(A_orig.size(0), A_orig.size(1));
            B.Resize(B_orig.size(0), B_orig.size(1));
            A.deallocateData();
            B.deallocateData();
    
            // We set their memory 
            int offset = offsets_per_size_std[block_size][i];
            int data_size = block_size*fixed_rank * dim_dof*dim_dof;
    
            // Sanity check : ensure memory is allocated
            if (offset + data_size > LR_standard_size_data_buffer_sizes_[block_size]){
                std::cerr << "Std size LR blocks : setting data outside the allocated buffer\n";
            }
                    
            A.setData(&LR_standard_size_A_data_[block_size][offset]);
            B.setData(&LR_standard_size_B_data_[block_size][offset]);
    
            // We copy the data
            auto A_edit = A.Edit();
            auto B_edit = B.Edit();
            auto A_orig_view = A_orig.view();
            auto B_orig_view = B_orig.view();
    
            for (size_t i = 0; i < A_orig.size(0); i++) {
                for (size_t j = 0; j < A_orig.size(1); j++) {
                    A_edit(i,j) = A_orig_view(i ,j);
                }
            }
            for (size_t i = 0; i < B_orig.size(0); i++) {
                for (size_t j = 0; j < B_orig.size(1); j++) {
                    B_edit(i,j) = B_orig_view(i ,j);
                }
            }
    
            // We also copy the error on approximation
            lrb->error_on_approximation = lra->error_on_approximation;
    
            // Finally we move it
            this->low_rank_blocks_[idx] = std::move(lrb);
        }
    }

    this->isBuilt_LR_ = true;
}                      


template <> 
il::Array<double> HmatCuda<double>::matvec(il::ArrayView<double> x) {

    #ifdef DEBUG
    std::cout << "Entering HmatCuda::matvec" << std::endl;
    #endif

    #ifdef WRITE_INPUT_OUTPUT_VEC
    std::stringstream ss1;
    ss1 << "mat_vec_gpu_x.npy";
    std::string filename1 = ss1.str();
    cnpy::npy_save(filename1, x.data(), {static_cast<size_t>(vector_size_)}, "w");
    #endif

    #ifdef TIMING
    struct timespec start, end;
    double duration;
    clock_gettime(CLOCK_MONOTONIC, &start);
    #endif // TIMING 

    #ifdef USE_NVTX
    NVTX3_FUNC_RANGE();   
    #endif

    il::Array<double> y(vector_size_, 0.0); // This also should be allocated once and for all 

    // #pragma omp parallel num_threads(this->n_openMP_threads_)
    #pragma omp parallel num_threads(num_gpus_+1)
    {   

        int thread_id = omp_get_thread_num();

        il::Array<double> y_private(vector_size_, 0.0); // This also should be allocated once and for all 

        // In a first num_gpu therads manafge GPU computation
        if (thread_id < num_gpus_) {

            int gpu_id = thread_id;

            CHECK_CUDA_ERROR(cudaSetDevice(gpu_id));

            CHECK_CUDA_ERROR(cudaMemcpy(d_x_[gpu_id], x.data(), vector_size_bytes_, cudaMemcpyHostToDevice));
            CHECK_CUDA_ERROR(cudaMemset(d_y_[gpu_id], 0, vector_size_bytes_));
            
            double alpha = 1.0;
            double beta = 0.0;

            // FR std blocks computation
            int fr_block_size = this->hr_->leaf_size * this->dof_dimension_;

            #ifdef USE_CUSPARSE_BSRMV
                // -> assigned to the first stream :
                cusparseSetStream(cusparse_handle_[gpu_id], cuda_streams_[gpu_id][0]);

                cusparseDbsrmv(
                    cusparse_handle_[gpu_id],
                    CUSPARSE_DIRECTION_COLUMN,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    total_block_size_,      // # block rows
                    total_block_size_,      // # block columns
                    num_FR_per_gpu_[gpu_id],     // # non zero blocks
                    &alpha,                 
                    FR_bsr_descr_[gpu_id],              
                    d_FR_data_[gpu_id],            
                    d_FR_bsrRowPtr_[gpu_id],      
                    d_FR_bsrColInd_[gpu_id],
                    fr_block_size,          // Block sizes
                    d_x_[gpu_id],
                    &beta,
                    d_y_[gpu_id]                    // FR blocks results written directly in d_y_
                );

            #else 

                CHECK_CUBLAS_ERROR(cublasSetStream(cublas_handle_[gpu_id], cuda_streams_[gpu_id][0]));

                #ifdef USE_MAGMA

                    magmablas_dgemv_batched(
                        MagmaNoTrans, 
                        this->hr_->leaf_size * this->dof_dimension_,  
                        this->hr_->leaf_size * this->dof_dimension_,  
                        alpha,
                        (const double**)d_FR_data_pointers_[gpu_id],
                        this->hr_->leaf_size * this->dof_dimension_, 
                        (const double**)d_FR_x_pointers_[gpu_id], 
                        1,
                        beta,
                        d_FR_y_partial_pointers_[gpu_id],
                        1,
                        num_FR_per_gpu_[gpu_id],
                        magma_queues_[gpu_id][0]
                    );

                #else 

                    #if CUDART_VERSION >= 11620 // CUDA 11.6.2 or newer

                        CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N,
                            this->hr_->leaf_size * this->dof_dimension_, // num rows
                            this->hr_->leaf_size * this->dof_dimension_, // num cols 
                            &alpha,
                            (const double**)d_FR_data_pointers_[gpu_id],        
                            this->hr_->leaf_size * this->dof_dimension_, 
                            (const double**)d_FR_x_pointers_[gpu_id], 
                            1,
                            &beta,
                            d_FR_y_partial_pointers_[gpu_id],
                            1,
                            num_FR_per_gpu_[gpu_id]
                        ));

                    #else 

                        CHECK_CUBLAS_ERROR(cublasDgemmBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N, CUBLAS_OP_N, 
                            this->hr_->leaf_size * this->dof_dimension_,
                            1, // 1 column matrix 
                            this->hr_->leaf_size * this->dof_dimension_, 
                            &alpha,
                            (const double**)d_FR_data_pointers_[gpu_id],
                            this->hr_->leaf_size * this->dof_dimension_, 
                            (const double**)d_FR_x_pointers_[gpu_id], 
                            this->hr_->leaf_size * this->dof_dimension_, 
                            &beta,
                            d_FR_y_partial_pointers_[gpu_id],
                            this->hr_->leaf_size * this->dof_dimension_,
                            num_FR_per_gpu_[gpu_id]
                        ));


                    #endif // CUBLAS version

                #endif // MAGMA BLAS or CUBLAS

                // FR non std blocks computation
                if (gpu_id == 0){

                    for (auto& [sizes,num_blocks] : num_FR_nonstd_blocks_per_size_){

                        #ifdef USE_MAGMA

                            magmablas_dgemv_batched(
                                MagmaNoTrans, 
                                sizes.first * this->dof_dimension_,  
                                sizes.second * this->dof_dimension_,  
                                alpha,
                                (const double**)d_FR_nonStd_data_pointers_[sizes], 
                                sizes.first * this->dof_dimension_,
                                (const double**)d_FR_nonStd_x_pointers_[sizes],
                                1,
                                beta,
                                d_FR_nonStd_y_partial_pointers_[sizes],
                                1,
                                num_blocks,
                                magma_queues_[gpu_id][0]
                            );

                        #else 

                            #if CUDART_VERSION >= 11620 // CUDA 11.6.2 or newer

                                CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                                    cublas_handle_[gpu_id], CUBLAS_OP_N,
                                    sizes.first * this->dof_dimension_, // num rows
                                    sizes.second * this->dof_dimension_, // num cols 
                                    &alpha,
                                    (const double**)d_FR_nonStd_data_pointers_[sizes],        
                                    sizes.first * this->dof_dimension_,
                                    (const double**)d_FR_nonStd_x_pointers_[sizes], 
                                    1,
                                    &beta,
                                    d_FR_nonStd_y_partial_pointers_[sizes],
                                    1,
                                    num_blocks
                                ));

                            #else 

                                CHECK_CUBLAS_ERROR(cublasDgemmBatched(
                                    cublas_handle_[gpu_id], CUBLAS_OP_N, CUBLAS_OP_N, 
                                    sizes.first * this->dof_dimension_,
                                    1, // 1 column matrix 
                                    sizes.second * this->dof_dimension_, 
                                    &alpha,
                                    (const double**)d_FR_nonStd_data_pointers_[sizes],
                                    sizes.first * this->dof_dimension_, 
                                    (const double**)d_FR_x_pointers_[gpu_id], 
                                    (const double**)d_FR_nonStd_x_pointers_[sizes], 
                                    &beta,
                                    d_FR_nonStd_y_partial_pointers_[sizes],
                                    sizes.second * this->dof_dimension_,
                                    num_blocks
                                ));

                            #endif // CUBLAS version

                        #endif // MAGMA BLAS or CUBLAS

                    }
                }
                
                // cudaDeviceSynchronize();
                cudaStreamSynchronize(cuda_streams_[gpu_id][0]);

                int num_blocks = num_FR_per_gpu_[gpu_id];
                if (gpu_id == 0) num_blocks += FR_non_std_indices.size();


                // And we sum the partial results
                scatter_add(
                    d_y_[gpu_id],
                    d_y_partial_FR_[gpu_id],
                    d_FR_y_partial_src_indices_[gpu_id],
                    d_FR_y_partial_dest_indices_[gpu_id],
                    d_FR_y_partial_lengths_[gpu_id],
                    num_blocks
                );
    
                // cudaDeviceSynchronize();
    
                // auto err = cudaGetLastError();
                // if (err != cudaSuccess) {
                //     printf("CUDA error after scatter_add: %s\n", cudaGetErrorString(err));
                // }


            #endif

            #ifdef WRITE_TMP_RES
            il::Array<double> y_fr_tmp(vector_size_, 0.0);

            double* y_fr_data = y_fr_tmp.Data();
            CHECK_CUDA_ERROR(cudaMemcpy(y_fr_data, d_y_[gpu_id], vector_size_bytes_, cudaMemcpyDeviceToHost));

            std::stringstream ss;
            ss << "gpu" << gpu_id << "_" << num_gpus_ <<  "_y_fr.npy";
            std::string filename = ss.str();
            cnpy::npy_save(filename, y_fr_data, {static_cast<size_t>(vector_size_)}, "w");

            double y_fr_sum = 0;
            for (int i(0); i<vector_size_; i++) y_fr_sum += y_fr_data[i];
            std::cout << "[GPU "<< gpu_id << "] std FR : y.sum() = " << std::scientific << std::setprecision(17) << y_fr_sum << std::endl;
            #endif


            // Then for the LR blocks : we loop on the sizes associated to this GPU

            // -> assigned to the second stream :
            CHECK_CUBLAS_ERROR(cublasSetStream(cublas_handle_[gpu_id], cuda_streams_[gpu_id][1]));

            for (int block_size : LR_std_sizes_per_gpu_[gpu_id]){

                int num_blocks = LR_std_indices_[block_size].size();

                    // #ifdef DEBUG
                    // std::cout << "Calling cublasDgemvBatched on the "<< num_blocks << " LR blocks of size " << block_size << std::endl;
                    // #endif

                #ifdef USE_MAGMA // Use MAGMA BLAS implementation

                    // Compute tmp = B*x
                    magmablas_dgemv_batched(
                        MagmaNoTrans, 
                        this->fixed_rank_*this->dof_dimension_,        // num rows
                        block_size*this->dof_dimension_, // num cols 
                        alpha,
                        (const double**)d_LR_B_data_pointers_[gpu_id][block_size],
                        this->fixed_rank_*this->dof_dimension_, 
                        (const double**)d_LR_x_pointers_[gpu_id][block_size], 
                        1,
                        beta,
                        d_LR_tmp_pointers_[gpu_id][block_size], 
                        1,
                        num_blocks,
                        magma_queues_[gpu_id][0]
                    );

                    // Compute y = A*tmp
                    magmablas_dgemv_batched(
                        MagmaNoTrans, 
                        block_size*this->dof_dimension_,        // num rows
                        this->fixed_rank_*this->dof_dimension_, // num cols 
                        alpha,
                        (const double**)d_LR_A_data_pointers_[gpu_id][block_size],        
                        block_size*this->dof_dimension_, 
                        (const double**)d_LR_tmp_pointers_[gpu_id][block_size],  
                        1,
                        beta,
                        d_LR_y_pointers_[gpu_id][block_size], 
                        1,
                        num_blocks,
                        magma_queues_[gpu_id][0]
                    );  

                #else // Use CUBLAS implementation

                    #if CUDART_VERSION >= 11620 // CUDA 11.6.2 or newer
                        // Compute tmp = B*x
                        CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N,
                            this->fixed_rank_*this->dof_dimension_,        // num rows
                            block_size*this->dof_dimension_, // num cols 
                            &alpha,
                            (const double**)d_LR_B_data_pointers_[gpu_id][block_size],        
                            this->fixed_rank_*this->dof_dimension_, 
                            (const double**)d_LR_x_pointers_[gpu_id][block_size], 
                            1,
                            &beta,
                            d_LR_tmp_pointers_[gpu_id][block_size], 
                            1,
                            num_blocks
                        ));
                
                        // Compute y = A*tmp
                        CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N,
                            block_size*this->dof_dimension_,        // num rows
                            this->fixed_rank_*this->dof_dimension_, // num cols 
                            &alpha,
                            (const double**)d_LR_A_data_pointers_[gpu_id][block_size],        
                            // this->fixed_rank_*this->dof_dimension_,
                            block_size*this->dof_dimension_, 
                            (const double**)d_LR_tmp_pointers_[gpu_id][block_size], 
                            1,
                            &beta,
                            d_LR_y_pointers_[gpu_id][block_size], 
                            1,
                            num_blocks
                        ));

                    #else // For CUDA versions before 11.6.2
                        
                        // Compute tmp = B*x
                        CHECK_CUBLAS_ERROR(cublasDgemmBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N, CUBLAS_OP_N, 
                            this->fixed_rank_*this->dof_dimension_,
                            1, // 1 column matrix ß
                            block_size*this->dof_dimension_, 
                            &alpha,
                            (const double**)d_LR_B_data_pointers_[gpu_id][block_size],
                            this->fixed_rank_*this->dof_dimension_, 
                            (const double**)d_LR_x_pointers_[gpu_id][block_size],
                            block_size*this->dof_dimension_, 
                            &beta,
                            d_LR_tmp_pointers_[gpu_id][block_size],
                            this->fixed_rank_*this->dof_dimension_,
                            num_blocks
                        ));

                        // Compute y = A*tmp 
                        CHECK_CUBLAS_ERROR(cublasDgemmBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N, CUBLAS_OP_N,
                            block_size*this->dof_dimension_, 
                            1, // 1 column matrix
                            this->fixed_rank_*this->dof_dimension_, 
                            &alpha,
                            (const double**)d_LR_A_data_pointers_[gpu_id][block_size],
                            block_size*this->dof_dimension_,
                            (const double**)d_LR_tmp_pointers_[gpu_id][block_size],
                            this->fixed_rank_*this->dof_dimension_,
                            &beta,
                            d_LR_y_pointers_[gpu_id][block_size],
                            block_size*this->dof_dimension_,
                            num_blocks
                        ));

                    #endif // CUDA version

                #endif // MAGMA BLAS or CUBLAS

            } // Loop on std LR block sizes

            // We add the contributions of the non std blocks
            if (gpu_id == 0){

                for (auto& [sizes,num_blocks] : num_LR_nonstd_blocks_per_size_){

                    #ifdef USE_MAGMA // Use MAGMA BLAS implementation

                         // Compute tmp = B*x
                        magmablas_dgemv_batched(
                            MagmaNoTrans, 
                            this->fixed_rank_*this->dof_dimension_,        // num rows
                            sizes.second*this->dof_dimension_, // num cols 
                            alpha,
                            (const double**)d_LR_nonStd_B_data_pointers_[sizes],
                            this->fixed_rank_*this->dof_dimension_, 
                            (const double**)d_LR_nonStd_x_pointers_[sizes], 
                            1,
                            beta,
                            d_LR_nonStd_tmp_pointers_[sizes], 
                            1,
                            num_blocks,
                            magma_queues_[gpu_id][0]
                        );

                        // Compute y = A*tmp
                        magmablas_dgemv_batched(
                            MagmaNoTrans, 
                            sizes.first*this->dof_dimension_,         // num rows
                            this->fixed_rank_*this->dof_dimension_, // num cols 
                            alpha,
                            (const double**)d_LR_nonStd_A_data_pointers_[sizes],        
                            sizes.first*this->dof_dimension_,  
                            (const double**)d_LR_nonStd_tmp_pointers_[sizes], 
                            1,
                            beta,
                            d_LR_nonStd_y_partial_pointers_[sizes],
                            1,
                            num_blocks,
                            magma_queues_[gpu_id][0]
                        );  

                    #else // Use CUBLAS implementation

                        #if CUDART_VERSION >= 11620 // CUDA 11.6.2 or newer

                            // Compute tmp = B*x
                            CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                                cublas_handle_[gpu_id], CUBLAS_OP_N,
                                this->fixed_rank_*this->dof_dimension_,     // num rows
                                sizes.second*this->dof_dimension_,          // num cols 
                                &alpha,
                                (const double**)d_LR_nonStd_B_data_pointers_[sizes],        
                                this->fixed_rank_*this->dof_dimension_, 
                                (const double**)d_LR_nonStd_x_pointers_[sizes], 
                                1,
                                &beta,
                                d_LR_nonStd_tmp_pointers_[sizes], 
                                1,
                                num_blocks
                            ));

                            // Compute y = A*tmp
                            CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                                cublas_handle_[gpu_id], CUBLAS_OP_N,
                                sizes.first*this->dof_dimension_,        // num rows
                                this->fixed_rank_*this->dof_dimension_, // num cols 
                                &alpha,
                                (const double**)d_LR_nonStd_A_data_pointers_[sizes],        
                                sizes.first*this->dof_dimension_, 
                                (const double**)d_LR_nonStd_tmp_pointers_[sizes], 
                                1,
                                &beta,
                                d_LR_nonStd_y_partial_pointers_[sizes], 
                                1,
                                num_blocks
                            ));

                        #else 

                            // Compute tmp = B*x
                            CHECK_CUBLAS_ERROR(cublasDgemmBatched(
                                cublas_handle_[gpu_id], CUBLAS_OP_N, CUBLAS_OP_N, 
                                this->fixed_rank_*this->dof_dimension_,
                                1, // 1 column matrix ß
                                sizes.second*this->dof_dimension_,
                                &alpha,
                                (const double**)d_LR_nonStd_B_data_pointers_[sizes], 
                                this->fixed_rank_*this->dof_dimension_, 
                                (const double**)d_LR_nonStd_x_pointers_[sizes], 
                                sizes.second*this->dof_dimension_, 
                                &beta,
                                d_LR_nonStd_tmp_pointers_[sizes], 
                                this->fixed_rank_*this->dof_dimension_,
                                num_blocks
                            ));

                            // Compute y = A*tmp 
                            CHECK_CUBLAS_ERROR(cublasDgemmBatched(
                                cublas_handle_[gpu_id], CUBLAS_OP_N, CUBLAS_OP_N,
                                sizes.first*this->dof_dimension_, 
                                1, // 1 column matrix
                                this->fixed_rank_*this->dof_dimension_, 
                                &alpha,
                                (const double**)d_LR_nonStd_A_data_pointers_[sizes],
                                sizes.first*this->dof_dimension_,
                                (const double**)d_LR_nonStd_tmp_pointers_[sizes],
                                this->fixed_rank_*this->dof_dimension_,
                                &beta,
                                d_LR_nonStd_y_partial_pointers_[sizes], 
                                sizes.first*this->dof_dimension_,
                                num_blocks
                            ));



                        #endif // CUDA version

                    #endif // MAGMA BLAS or CUBLAS


                }
            }

            // We sync the streams
            // cudaDeviceSynchronize();
            cudaStreamSynchronize(cuda_streams_[gpu_id][1]);

            // Here we gather the partial results of the LR blocks on d_y_ that was used to store the results of the FR blocks
            int num_lr_blocks = num_LR_per_gpu_[gpu_id];
            if (gpu_id == 0) num_lr_blocks += LR_non_std_indices_.size();

            scatter_add(
                d_y_[gpu_id],
                d_y_partial_LR_[gpu_id],
                d_LR_y_partial_src_indices_[gpu_id],
                d_LR_y_partial_dest_indices_[gpu_id],
                d_LR_y_partial_lengths_[gpu_id],
                num_lr_blocks
            );

            cudaDeviceSynchronize();

            auto err = cudaGetLastError();
            if (err != cudaSuccess) {
                printf("CUDA error after scatter_add: %s\n", cudaGetErrorString(err));
            }
 
            // Finnally, copy it back to cpu
            double* y_data = y_private.Data();
            CHECK_CUDA_ERROR(cudaMemcpy(y_data, d_y_[gpu_id], vector_size_bytes_, cudaMemcpyDeviceToHost));

            #ifdef WRITE_TMP_RES
            std::stringstream ss2;
            ss2 << "gpu" << gpu_id << "_" << num_gpus_ <<  "_y_fr_lr.npy";
            std::string filename2 = ss2.str();            
            cnpy::npy_save(filename2, y_fr_data, {static_cast<size_t>(vector_size_)}, "w");
            
            std::stringstream ss3;
            ss3 << "gpu" << gpu_id << "_" << num_gpus_ <<  "_y_part_lr.npy";
            double* y_partial = new double[y_partial_LR_buffer_size_bytes_[gpu_id]/sizeof(double)];
            CHECK_CUDA_ERROR(cudaMemcpy(y_partial, d_y_partial_LR_[gpu_id], y_partial_LR_buffer_size_bytes_[gpu_id], cudaMemcpyDeviceToHost));
            filename2 = ss3.str();
            cnpy::npy_save(filename2, y_partial, {static_cast<size_t>(y_partial_LR_buffer_size_bytes_[gpu_id]/sizeof(double))}, "w");
            delete[] y_partial;

            std::stringstream ss4;
            ss4 << "gpu" << gpu_id << "_" << num_gpus_ <<  "_tmp_lr.npy";
            double* tmp = new double[tmp_buffer_size_bytes_[gpu_id]/sizeof(double)];
            CHECK_CUDA_ERROR(cudaMemcpy(tmp, d_tmp_[gpu_id], tmp_buffer_size_bytes_[gpu_id], cudaMemcpyDeviceToHost));
            filename2 = ss4.str();
            cnpy::npy_save(filename2, tmp, {static_cast<size_t>(tmp_buffer_size_bytes_[gpu_id]/sizeof(double))}, "w");
            delete[] tmp;

            double y_f_sum = 0;
            for (int i(0); i<vector_size_; i++) y_f_sum += y_data[i];
            std::cout << "[GPU "<< gpu_id << "] std FR + std LR : y.sum() = " << std::scientific << std::setprecision(17) << y_f_sum << std::endl;
            #endif

        } // GPU - assigned threads



        // We sum the partial results of each thread
        #pragma omp critical
        {
            il::blas(1., y_private.view(), il::io_t{}, y.Edit());
        }
    
    } // pragma omp parallel

    #ifdef TIMING
    clock_gettime(CLOCK_MONOTONIC, &end);
    duration = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;    
    std::cout << "[Timing] HmatCuda<double>::matvec : " << duration*1000 << "ms\n";
    #endif  

    #ifdef WRITE_INPUT_OUTPUT_VEC
    std::stringstream ss2;
    ss2 << "mat_vec_gpu_y.npy";
    std::string filename2 = ss2.str();
    cnpy::npy_save(filename2, y.data(), {static_cast<size_t>(vector_size_)}, "w");
    #endif

    return y;
}

template <typename T>
void HmatCuda<T>::setDiagonal(){

    // Store the diagonal for later use as Jacobi preconditionner
    IL_EXPECT_FAST(this->isBuilt_FR_);
    il::int_t diag_size = il::max(this->size_[0], this->size_[1]);
    diagonal_.resize(static_cast<long>(diag_size));

    for (il::int_t k = 0; k < this->hr_->pattern_.n_FRB; k++) {
        il::Array2D<double> aux = *this->full_rank_blocks_[k];
        il::int_t i0 = this->hr_->pattern_.FRB_pattern(1, k);
        il::int_t j0 = this->hr_->pattern_.FRB_pattern(2, k);
        il::int_t iend = this->hr_->pattern_.FRB_pattern(3, k);
        il::int_t jend = this->hr_->pattern_.FRB_pattern(4, k);
        // check if it intersects the diagonal
        bool in_lower = (i0 > j0) && (i0 > jend) && (iend > j0) && (iend > jend);
        bool in_upper = (i0 < j0) && (i0 < jend) && (iend < j0) && (iend < jend);
        if ((!in_lower) && (!in_upper)) // this fb intersect the diagonal....
        {
            for (il::int_t j = 0; j < aux.size(1); j++) {
                for (il::int_t i = 0; i < aux.size(0); i++) {
                    if ((i + this->dof_dimension_ * i0) ==
                        (j + this->dof_dimension_ * j0)) { // on diagonal !
                        diagonal_[(i + this->dof_dimension_ * i0)] = aux(i, j);
                    }
                }
            }
        }
    }

    diagonal_stored_ = true;
}

template <typename T>
std::vector<T> HmatCuda<T>::diagonal(){

    if (!diagonal_stored_){
        std::cout << "[ERROR] Accessing the diagonal before building the hmat.";
    }
    return diagonal_;
}

template <typename T>
std::vector<T> HmatCuda<T>::diagonalOriginal(){

    if (!diagonal_stored_){
        std::cout << "[ERROR] Accessing the diagonal before building the hmat.";
    }

    // return diagonal in original state....
    il::int_t diag_size = il::max(this->size_[0], this->size_[1]);
    il::int_t ncolpoints = diag_size / this->dof_dimension_;
    std::vector<T> diagonal_perm(diagonal_.size(), 0.);
    // permut back
    for (il::int_t i = 0; i < ncolpoints; i++) {
        for (int j = 0; j < this->dof_dimension_; j++) {
            diagonal_perm[this->dof_dimension_ * this->hr_->permutation_1_[i] + j] = diagonal_[this->dof_dimension_ * i + j];
        }
    }
    return diagonal_perm;
}

template <typename T> il::int_t HmatCuda<T>::nbOfEntries() {

    // IL_EXPECT_FAST(this->isBuilt_);

    auto hr = this->hr_;
    int fixed_rank = this->fixed_rank_;
    int leaf_size = this->hr_->leaf_size;
    const int dim_dof = this->dof_dimension_;

    il::int_t n = 0;

    for (il::int_t i = 0; i < hr->pattern_.n_FRB; i++) {
        il::int_t i0 = hr->pattern_.FRB_pattern(1, i);
        il::int_t j0 = hr->pattern_.FRB_pattern(2, i);
        il::int_t iend = hr->pattern_.FRB_pattern(3, i);
        il::int_t jend = hr->pattern_.FRB_pattern(4, i);

        n += (iend-i0)*(jend-j0) * dim_dof*dim_dof;
    }
    for (il::int_t i = 0; i < hr->pattern_.n_LRB; i++) {

        il::int_t i0 = hr->pattern_.FRB_pattern(1, i);
        il::int_t j0 = hr->pattern_.FRB_pattern(2, i);
        il::int_t iend = hr->pattern_.FRB_pattern(3, i);
        il::int_t jend = hr->pattern_.FRB_pattern(4, i);

        n += (iend-i0)*fixed_rank * dim_dof*dim_dof + (jend-j0)*fixed_rank * dim_dof*dim_dof;
    }
    return n;
}

template class HmatCuda<double>;


}