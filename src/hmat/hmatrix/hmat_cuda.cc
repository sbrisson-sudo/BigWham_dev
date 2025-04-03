#include <unordered_map>
#include <ctime>
#include <iostream>
#include <string>
#include <iomanip>

#include "hmat_cuda.h"
#include "y_partial_scatter_add.h"

#include "cnpy.h"

// #define TIMING
// #define DEBUG
// #define DEBUG2
// #define DEBUG4

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
HmatCuda<T>::HmatCuda(const bigwham::MatrixGenerator<T> & matrix_gen, const double epsilon_aca, const int n_openMP_threads, const bool verbose, const int fixed_rank) {

#ifdef TIMING
    struct timespec start, end;
    double duration;
    clock_gettime(CLOCK_MONOTONIC, &start);
#endif // TIMING 

    // construction directly
    this->verbose_ = verbose;
    this->fixed_rank_ = fixed_rank;
    this->n_openMP_threads_=n_openMP_threads;

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
    int total_size_non_standard_blocks = 0;

    for (int i(0); i<pattern.n_FRB; i++){
        int i0 = pattern.FRB_pattern(1,i);
        int j0 = pattern.FRB_pattern(2,i);
        int iend = pattern.FRB_pattern(3,i);
        int jend = pattern.FRB_pattern(4,i);

        if ((iend-i0 == leaf_size) && (jend-j0 == leaf_size)){ 
            num_FR_blocks_standard_size++;
        } else {
            num_FR_blocks_non_standard_size++;
            total_size_non_standard_blocks += (jend-j0)*(iend-i0);
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

    // Print it 
#ifdef DEBUG
    std::cout << " FR blocks std order : [";
    for (auto i : FR_std_orderedIndices){
        std::cout << i << ", ";
    }
    std::cout << "]\n";
#endif


    // We allocate the memory for standardized FR blocks

    if (this->verbose_){
        std::cout << num_FR_blocks_standard_size << " FR blocks with standard size\n";
        std::cout << num_FR_blocks_non_standard_size << " FR blocks with non standard size\n";
        // std::cout << "Allocating " << num_FR_blocks_standard_size * leaf_size * leaf_size * dim*dim << " for FR_standard\n";
        // std::cout << "Allocating " << total_size_non_standard_blocks * dim*dim << " for FR_non_standard\n";
    }

    FR_standard_size_data_buffer_size = num_FR_blocks_standard_size * leaf_size * leaf_size * dim*dim;
    this->FR_standard_size_data = new T[FR_standard_size_data_buffer_size];

    FR_non_standard_size_data_buffer_size = total_size_non_standard_blocks * dim*dim;
    this->FR_non_standard_size_data = new T[FR_non_standard_size_data_buffer_size];

    // LR MEM ALLOCATION
    // We want to group the LR blocks in groups that :
    // - have the same size 
    // - don't share the same rows

    // FIRST : we cound the number of LR blocks for each size
    std::unordered_map<int, int> num_LR_blocks_per_size;
    std::unordered_map<int, std::vector<int>> LR_blocks_per_size_indices;

    for (int i(0); i<pattern.n_LRB; i++){
        int i0 = pattern.LRB_pattern(1,i);
        int j0 = pattern.LRB_pattern(2,i);
        int iend = pattern.LRB_pattern(3,i);
        int jend = pattern.LRB_pattern(4,i);
        // Sanity check
        if ((iend-i0) != (jend-j0)){
            std::cerr << "Error when counting the LR blocks, some are rectangular\n";
        }
        // Incrementing the count
        // Note : implicitely initialized to zero
        num_LR_blocks_per_size[iend-i0]++; 
        LR_blocks_per_size_indices[iend-i0].push_back(i);

    }

    // Just print it 
#ifdef DEBUG
    for (const auto& size_count : num_LR_blocks_per_size) {
        int block_size = size_count.first;
        int block_count = size_count.second;
        std::cout << block_count << " LR blocks of size " << block_size << " : [";
        for (auto i : LR_blocks_per_size_indices[block_size]){
            std::cout << i << ", ";
        }
        std::cout << "]\n";
    }
#endif

    // SECOND : we constrcut group of same size with no common rows
    // We look ove each size and construct group that dont share common rows
    std::vector<bool> LR_blocks_picked(pattern.n_LRB, false);
    int num_LR_groups = 0;

    for (const auto& size_count : num_LR_blocks_per_size) {

        int block_size = size_count.first;
        int block_count = size_count.second;

        // Initialize a group
        int group_id = 0;
        LR_blocks_groups_indices[block_size].emplace_back();

        // Num of picked
        int num_picked = 0;

        while (true){

            for (auto i : LR_blocks_per_size_indices[block_size]){
                if (!LR_blocks_picked[i]){
    
                    // See if this block share common rows with the ones already in the block
                    bool no_common_row = true;
                    for (auto j : LR_blocks_groups_indices[block_size][group_id]){
                        int i_row0 = pattern.LRB_pattern(1,i);
                        int i_rowend = pattern.LRB_pattern(3,i);
                        int j_row0 = pattern.LRB_pattern(1,j);
                        int j_rowend = pattern.LRB_pattern(3,j);
    
                        if ((i_row0 < j_rowend) && (j_row0 < i_rowend)){
                            no_common_row = false;
                            break;
                        }
                    }
    
                    // If not : we add it 
                    if (no_common_row){
                        LR_blocks_groups_indices[block_size][group_id].push_back(i);
                        LR_blocks_picked[i] = true;
                        num_picked++;
                    }
                }
            }

            if (num_picked == block_count) break;

            // New group
            group_id++;
            num_LR_groups++;
            LR_blocks_groups_indices[block_size].emplace_back();
        }
    }

    // We set a global group block indexation for convenience
    LR_group_size_num_.reserve(num_LR_groups);
    for (const auto& size_count : num_LR_blocks_per_size) {
        int block_size = size_count.first;
        for (int group_i(0); group_i < LR_blocks_groups_indices[block_size].size() ; group_i++){
            int num_blocks_group = LR_blocks_groups_indices[block_size][group_i].size();
            LR_group_size_num_.push_back({block_size, num_blocks_group});
        }
    }


#ifdef DEBUG3
    std::cout << "LR_group_size_num_ = [";
    for (auto& [block_size, group_len] : LR_group_size_num_) std::cout << "(" << block_size << "," << group_len << "),";
    std::cout << "]\n";
#endif


#ifdef DEBUG
    // Just print it 
    for (const auto& size_count : num_LR_blocks_per_size) {
        int block_size = size_count.first;
        int block_count = size_count.second;
        int num_groups = LR_blocks_groups_indices[block_size].size();
        std::cout << num_groups << " LR groups of size " << block_size << " : [";

        for (auto group : LR_blocks_groups_indices[block_size]){
            std::cout << "[";
            for (size_t j = 0; j < group.size(); j++) {
                std::cout << group[j];
                if (j < group.size() - 1) std::cout << ", ";
            }
            std::cout << "],";
        }
        std::cout << "]\n";
    }
#endif
    

    // We allocate the memory for the low rank blocks
    // this->LR_A_data.reserve(num_LR_blocks_per_size.size());
    // this->LR_B_data.reserve(num_LR_blocks_per_size.size());
    this->LR_sizes.reserve(num_LR_blocks_per_size.size());

    for (const auto& size_count : num_LR_blocks_per_size) {

        int block_size = size_count.first;
        int block_count = size_count.second;

        // std::cout << "fixe_rank = " << fixed_rank << std::endl;

        // Loop on groups
        for (int group_i(0); group_i < LR_blocks_groups_indices[block_size].size() ; group_i++){

            int num_blocks_group = LR_blocks_groups_indices[block_size][group_i].size();

            int LR_data_buffer_size = num_blocks_group * block_size*fixed_rank * dim*dim;

            // std::cout << "Buffer size for to LR_A+_data[" << block_size << "][" <<group_i <<"] = " << LR_data_buffer_size << "*sizeof(T)\n";

            LR_data_buffer_sizes[block_size].push_back(LR_data_buffer_size); // Save the buffer size

            T* LR_A_buffer = new T[LR_data_buffer_size];
            T* LR_B_buffer = new T[LR_data_buffer_size];

            // std::cout << num_blocks_group << " LR blocks of size " << block_size << " in group " <<  group_i << std::endl;
            // std::cout << "Allocating " << LR_data_buffer_size << " for LR blocks of size " << block_size <<  " and group " << group_i << std::endl;

            this->LR_A_data[block_size].push_back(LR_A_buffer);
            this->LR_B_data[block_size].push_back(LR_B_buffer);

        }

        this->LR_sizes.emplace_back(block_size);
    }

    // Print total memory footprint
    if (this->verbose_){
        int total_size = FR_standard_size_data_buffer_size + FR_non_standard_size_data_buffer_size;
        for (auto& pair : LR_data_buffer_sizes){
            for (int buffer_size : pair.second){
                total_size += buffer_size;
            }
        }
    
        std::cout << "Total memory allocated on CPU = " << formatBytes(total_size*sizeof(double)) << " \n";
    }

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


    // And we copy it on device
    copyToDevice();

#ifdef TIMING
    clock_gettime(CLOCK_MONOTONIC, &end);
    duration = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;    
    std::cout << "[Timing] copied the hmat to device = " << duration*1000 << "ms\n";
    clock_gettime(CLOCK_MONOTONIC, &start);
#endif // TIMING 

}

template <typename T>
void HmatCuda<T>::copyToDevice(){

    // Getting the number of GPUs available
    cudaGetDeviceCount(&num_gpus_);
    // num_gpus_ = 1;

    if (this->verbose_){
        std::cout << "--------------------" << std::endl;
        std::cout << "Copying the Hmat to the GPU ..." << std::endl;
        std::cout << "Num of available GPU = " << num_gpus_ << std::endl;
    }

    auto hr = this->hr_;
    int leaf_size = this->hr_->leaf_size;
    const int dim_dof = this->dof_dimension_;

    // ---------------------------
    // Distributing load on GPUS
    // ---------------------------

    // FR standard blocks
    num_FR_per_gpu_.resize(num_gpus_);
    for (int gpu_id(0); gpu_id<num_gpus_-1; gpu_id++) num_FR_per_gpu_[gpu_id] = num_FR_std_blocks_/num_gpus_;
    num_FR_per_gpu_[num_gpus_-1] = num_FR_std_blocks_ - (num_gpus_-1)*num_FR_std_blocks_/num_gpus_;

    offsets_FR_gpu_.resize(num_gpus_);
    offsets_FR_gpu_[0] = 0;
    for (int gpu_id(1); gpu_id<num_gpus_; gpu_id++){
        offsets_FR_gpu_[gpu_id] = offsets_FR_gpu_[gpu_id-1] + num_FR_per_gpu_[gpu_id-1];
    } 

#ifdef DEBUG3
    std::cout << "Num of FR blocks per GPU = [";
    for (int n : num_FR_per_gpu_) std::cout << n << ",";
    std::cout << "]\n";
#endif

    // Distributng the LR block groups 
    num_LR_per_gpu_.resize(num_gpus_);
    for (int gpu_id(0); gpu_id<num_gpus_; gpu_id++) num_LR_per_gpu_[gpu_id] = 0;
    std::vector<int> gpu_load(num_gpus_, 0);    
    LR_group_per_gpu_.resize(num_streams_-1);
    for (auto&  [block_size, buffer_sizes] : LR_data_buffer_sizes){

        int num_groups = buffer_sizes.size();

        for (int group_i(0); group_i<num_groups; group_i++){

            // We find the less loaded gpu
            int min_load_gpu = 0;
            for (int i = 0; i < num_gpus_; ++i) {
                if (gpu_load[i] < gpu_load[min_load_gpu]) {
                    min_load_gpu = i;
                }
            }

            // We add the block to it
            if (LR_group_per_gpu_[min_load_gpu].find(block_size) == LR_group_per_gpu_[min_load_gpu].end()) {
                // We set up a new vector of grup id for this stream and this block size
                LR_group_per_gpu_[min_load_gpu][block_size] = {};
            }

            LR_group_per_gpu_[min_load_gpu][block_size].push_back(group_i);
            gpu_load[min_load_gpu] += buffer_sizes[group_i];
            num_LR_per_gpu_[min_load_gpu] += LR_blocks_groups_indices[block_size][group_i].size();
        }
    }

#ifdef DEBUG3
    std::cout << "gpu_load = [";
    for (int load : gpu_load) std::cout << load << ", ";
    std::cout << "]\n";

    std::cout << "LR_group_per_gpu_ = [\n";
    for (int gpu_id(0); gpu_id < num_gpus_; gpu_id++){
        std::cout << " - gpu # " << gpu_id+1 << " : [";
        for (auto&  [block_size, group_i_list] : LR_group_per_gpu_[gpu_id]){
            std::cout << " size # " << block_size << " : [";
            for (int group_i : group_i_list) std::cout << group_i << ", ";
            std::cout << "]";
        }
        std::cout << "]\n";
    }
    std::cout << "]\n";
#endif

    // ---------------------------
    // Distributing load on GPUS - done
    // ---------------------------

    // ---------------------------
    // Initializing the per GPU data structures
    // ---------------------------

    LR_group_per_cuda_stream_.resize(num_gpus_);

    // GPU operation handler
    cublas_handle_.resize(num_gpus_);
    cusparse_handle_.resize(num_gpus_);

    // Number of CUDA streams = 1 for BSR + the rest for batched opeartions
    cuda_streams_.resize(num_gpus_);

    // GPU (device) memory buffers
    d_x_.resize(num_gpus_);               // Store the lhs vector on device
    d_y_.resize(num_gpus_);               // Store the rhs vector on device
    d_y_partial_LR_.resize(num_gpus_);    // Partial results for LR blocks operations
    d_tmp_.resize(num_gpus_);          // tmp = B*x then y = A*tmp

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
        // Distributing load on CUDA streams
        // ---------------------------

        // Init the handles
        CHECK_CUBLAS_ERROR(cublasCreate(&cublas_handle_[gpu_id]));
        cusparseCreate(&cusparse_handle_[gpu_id]);

        // // First we create the CUDA streams
        // num_streams_ = 1; // for the BSR (FR) operation
        // for (const auto& pair : LR_data_buffer_sizes) {
        //     num_streams_ += pair.second.size();
        // }

        if (this->verbose_) std::cout << "[GPU "<< gpu_id << "] Initiating " << num_streams_ << " CUDA streams\n";
        cuda_streams_[gpu_id].resize(num_streams_);
        for (int i = 0; i < num_streams_; i++) CHECK_CUDA_ERROR(cudaStreamCreate(&cuda_streams_[gpu_id][i]));

        // Here we balance the load for the LR blocks over num_streams_ - 1 streams 

        // We set the associated cuda streams following a greedy approach : we give the block to the 
        // less loaded cuda stream
        std::vector<int> cuda_streams_load(num_streams_-1, 0);    
        LR_group_per_cuda_stream_[gpu_id].resize(num_streams_-1);

        for (auto&  [block_size, group_id_list] : LR_group_per_gpu_[gpu_id]){

            for (int group_i : group_id_list){   
                
                int buffer_size = LR_data_buffer_sizes[block_size][group_i];

                // We find the less loaded stream
                int min_load_stream = 0;
                for (int i = 1; i < num_streams_-1; ++i) {
                    if (cuda_streams_load[i] < cuda_streams_load[min_load_stream]) {
                        min_load_stream = i;
                    }
                }

                // We add the block to it
                if (LR_group_per_cuda_stream_[gpu_id][min_load_stream].find(block_size) == LR_group_per_cuda_stream_[gpu_id][min_load_stream].end()) {
                    // We set up a new vector of grup id for this stream and this block size
                    LR_group_per_cuda_stream_[gpu_id][min_load_stream][block_size] = {};
                }

                LR_group_per_cuda_stream_[gpu_id][min_load_stream][block_size].push_back(group_i);
                cuda_streams_load[min_load_stream] += buffer_size;
            }
        }

    #ifdef DEBUG2
        std::cout << "[GPU "<< gpu_id << "] cuda_streams_load = [";
        for (int load : cuda_streams_load) std::cout << load << ", ";
        std::cout << "]\n";

        std::cout << "[GPU "<< gpu_id << "] LR_group_per_cuda_stream_ = [\n";
        for (int cuda_stream_i(0); cuda_stream_i < num_streams_-1; cuda_stream_i++){
            std::cout << " - stream # " << cuda_stream_i+1 << " : [\n";
            for (auto&  [block_size, group_i_list] : LR_group_per_cuda_stream_[gpu_id][cuda_stream_i]){
                std::cout << "   - size # " << block_size << " : [";
                for (int group_i : group_i_list) std::cout << group_i << ", ";
                std::cout << "]\n";
            }
            std::cout << "]\n";
        }
        std::cout << "]\n";
    #endif

        // ---------------------------
        // Distributing load on CUDA streams - done
        // ---------------------------

        // ---------------------------
        // Initialiazing the aux vectors
        // ---------------------------

        // Then we set up the memory for the vectors
        if (this->size_[0] != this->size_[1]){
            std::cerr << "CUDA only implemented for square Hmat for now\n";
            std::cerr << "Hmat size = " << this->size_[0] << " x " << this->size_[1] << std::endl;
            std::abort();
        }
        vector_size_ = this->size_[0];
        vector_size_bytes_ = this->size_[0]* sizeof(T);

        // Here we need to compute the buffer size needed for tmp and y_partial so that each LR block has its owns
        // size allocated for y_partial = size for BSR (vector_size_) + sizes for LR blocks
        y_partial_LR_buffer_size_bytes_ = 0;
        tmp_buffer_size_bytes_ = 0;

        for (auto&  [block_size, group_id_list] : LR_group_per_gpu_[gpu_id]){
            for (int group_i : group_id_list){
                int num_blocks = LR_blocks_groups_indices[block_size][group_i].size();

                y_partial_LR_buffer_size_bytes_ += num_blocks * block_size * dim_dof * sizeof(T);
                tmp_buffer_size_bytes_ += num_blocks * this->fixed_rank_ * dim_dof * sizeof(T);
            }
        }
        if (this->verbose_){
            std::cout << "[GPU "<< gpu_id << "] Allocating " << formatBytes(2*vector_size_bytes_ + y_partial_LR_buffer_size_bytes_ + tmp_buffer_size_bytes_) << " on the GPU for auxilliary vectors" << std::endl;
        }

    #ifdef DEBUG2
        std::cout << "[GPU "<< gpu_id << "] y_partial_buffer_size_bytes = " << y_partial_LR_buffer_size_bytes_ << std::endl;
        std::cout << "[GPU "<< gpu_id << "] tmp_buffer_size_bytes = " << tmp_buffer_size_bytes_ << std::endl;
    #endif 


        CHECK_CUDA_ERROR(cudaMalloc(&d_x_[gpu_id], vector_size_bytes_));
        CHECK_CUDA_ERROR(cudaMalloc(&d_y_[gpu_id], vector_size_bytes_));

        // Here we change the sized allocated 
        CHECK_CUDA_ERROR(cudaMalloc(&d_y_partial_LR_[gpu_id], y_partial_LR_buffer_size_bytes_));
        CHECK_CUDA_ERROR(cudaMalloc(&d_tmp_[gpu_id], tmp_buffer_size_bytes_));

        // ---------------------------
        // Initialiazing the aux vectors - done
        // ---------------------------

        // ----------------------------
        // COPYING FULL RANKS BLOCKS
        // ----------------------------
        // Then we copy the FR blocks data (only of standard size)

        size_t FR_block_size = leaf_size*leaf_size*dim_dof*dim_dof;
        size_t FR_data_size_bytes = num_FR_per_gpu_[gpu_id] * FR_block_size * sizeof(T);

    #ifdef DEBUG4
        std::cout << "[GPU "<< gpu_id << "] Copying FR standard block to device = " << FR_data_size_bytes << " bytes\n";
        std::cout << "[GPU "<< gpu_id << "] offset[gpu_id] = " << offsets_FR_gpu_[gpu_id] << std::endl;
        std::cout << "[GPU "<< gpu_id << "] num_FR_per_gpu_[gpu_id] = " << num_FR_per_gpu_[gpu_id] << std::endl;
        std::cout << "[GPU "<< gpu_id << "] num_FR_per_gpu_[gpu_id] * FR_block_size = " << num_FR_per_gpu_[gpu_id] * FR_block_size << std::endl;
        std::cout << "[GPU "<< gpu_id << "] FR_standard_size_data_buffer_size = " << FR_standard_size_data_buffer_size << std::endl;

    #endif

        CHECK_CUDA_ERROR(cudaMalloc(&d_FR_data_[gpu_id], FR_data_size_bytes));
        CHECK_CUDA_ERROR(cudaMemcpy(d_FR_data_[gpu_id], FR_standard_size_data + offsets_FR_gpu_[gpu_id]*FR_block_size, FR_data_size_bytes, cudaMemcpyHostToDevice));


        if (this->verbose_){
            size_t total_to_allocate = FR_data_size_bytes;
            for (auto& pair : LR_data_buffer_sizes){
                for (int buffer_size : pair.second){
                    total_to_allocate += buffer_size * sizeof(T);
                }
            }
            std::cout << "[GPU "<< gpu_id << "] Allocating " << formatBytes(total_to_allocate) << " on the GPU for hierarchical matrix data" << std::endl;
        }

        // Then we set up the BSR description of the FR 
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

    #ifdef DEBUG4
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

            for (auto&  [block_size, group_id_list] : LR_group_per_gpu_[gpu_id]){
                for (int group_i : group_id_list){
                    int num_blocks = LR_blocks_groups_indices[block_size][group_i].size();
                    total_to_allocate_metadata += num_blocks * sizeof(T*) * 5;
                }
            }
            std::cout << "[GPU "<< gpu_id << "] Allocating " << formatBytes(total_to_allocate_metadata) << " on the GPU for hierarchical matrix metadata" << std::endl;

        }

        // Then copy it to device
        CHECK_CUDA_ERROR(cudaMalloc(&d_FR_bsrRowPtr_[gpu_id], (num_block_rows+1)*sizeof(int)));
        CHECK_CUDA_ERROR(cudaMemcpy(d_FR_bsrRowPtr_[gpu_id], h_FR_bsrRowPtr_[gpu_id], (num_block_rows+1)*sizeof(int), cudaMemcpyHostToDevice));

        CHECK_CUDA_ERROR(cudaMalloc(&d_FR_bsrColInd_[gpu_id], nnzb*sizeof(int)));
        CHECK_CUDA_ERROR(cudaMemcpy(d_FR_bsrColInd_[gpu_id], h_FR_bsrColInd_[gpu_id], nnzb*sizeof(int), cudaMemcpyHostToDevice));

        // ----------------------------
        // COPYING FULL RANKS BLOCKS - done
        // ----------------------------

        // ----------------------------
        // COPYING LOW RANKS BLOCKS
        // ----------------------------

        // We copy all buffers
        for (auto&  [block_size, group_id_list] : LR_group_per_gpu_[gpu_id]){
        // for (auto& pair : LR_data_buffer_sizes){

            for (int group_i : group_id_list){
            // for (int group_i(0); group_i<pair.second.size(); group_i++){
                int buffer_size = LR_data_buffer_sizes[block_size][group_i];

                // Allocation
                T* d_A_tmp;
                CHECK_CUDA_ERROR(cudaMalloc(&d_A_tmp, buffer_size*sizeof(T)));
                d_LR_A_data_[gpu_id][block_size][group_i] = d_A_tmp;
                T* d_B_tmp;
                CHECK_CUDA_ERROR(cudaMalloc(&d_B_tmp, buffer_size*sizeof(T)));
                d_LR_B_data_[gpu_id][block_size][group_i] = d_B_tmp;

                // std::cout << "LR Allocating " << buffer_size*sizeof(T) << "to d_LR_A_data_[" << block_size << "][" <<group_i <<"] @ " << d_LR_A_data_[block_size][group_i] << "\n";

                // Copy
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_A_data_[gpu_id][block_size][group_i], LR_A_data[block_size][group_i], buffer_size*sizeof(T), cudaMemcpyHostToDevice));
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_B_data_[gpu_id][block_size][group_i], LR_B_data[block_size][group_i], buffer_size*sizeof(T), cudaMemcpyHostToDevice));

            }
        }

        // Finally we need to set the pointers array to perform the batched operations

        size_t y_partial_LR_ptr_counter = 0;
        size_t tmp_ptr_counter = 0;

        // Offsets for summing partial results
        int* h_LR_y_partial_src_indices = new int[hr->pattern_.n_LRB];
        int* h_LR_y_partial_dest_indices = new int[hr->pattern_.n_LRB];
        int* h_LR_y_partial_lengths = new int[hr->pattern_.n_LRB];
        int lr_block_i = 0;

        for (auto&  [block_size, group_id_list] : LR_group_per_gpu_[gpu_id]){
        // for (auto& pair : LR_data_buffer_sizes){
        //     int block_size = pair.first;
        //     int num_groups = pair.second.size();

            int num_groups = group_id_list.size();
            // std::cout << "- SIZE = " << block_size << "\n";

            int block_data_size = block_size*this->fixed_rank_ * dim_dof*dim_dof;

            for (int group_i : group_id_list){
            // for (int group_i(0); group_i<num_groups; group_i++){

                // std::cout << "  - GROUP = " << group_i << "\n";

                // std::cout << "Base pointers = \n";
                // std::cout << "      - &A = " << d_LR_A_data_[block_size][group_i] << "\n";
                // std::cout << "      - &B = " << d_LR_B_data_[block_size][group_i] << "\n";
                // std::cout << "      - &x = " << d_x_ << "\n";
                // std::cout << "      - &y = " << d_y_partial_ +  vector_size_*count_y_partial_offset << "\n";
                // std::cout << "      - &tmp = " << d_tmp_ + vector_size_*(count_y_partial_offset-1) << "\n";

                int num_blocks = LR_blocks_groups_indices[block_size][group_i].size();

                h_LR_A_data_pointers_[gpu_id][block_size][group_i] = new T*[num_blocks];
                h_LR_B_data_pointers_[gpu_id][block_size][group_i] = new T*[num_blocks];
                h_LR_x_pointers_[gpu_id][block_size][group_i] = new T*[num_blocks];
                h_LR_y_pointers_[gpu_id][block_size][group_i] = new T*[num_blocks];
                h_LR_tmp_pointers_[gpu_id][block_size][group_i] = new T*[num_blocks];

                for (int i(0); i<num_blocks; i++){

                    // std::cout << "    - BLOCK = " << i << "\n";

                    // Data = linear
                    h_LR_A_data_pointers_[gpu_id][block_size][group_i][i] = d_LR_A_data_[gpu_id][block_size][group_i] + i*block_data_size;
                    h_LR_B_data_pointers_[gpu_id][block_size][group_i][i] = d_LR_B_data_[gpu_id][block_size][group_i] + i*block_data_size;

                    // Vectos = depends on block position
                    int block_id = LR_blocks_groups_indices[block_size][group_i][i];
                    int i0 = hr->pattern_.LRB_pattern(1, block_id);
                    int j0 = hr->pattern_.LRB_pattern(2, block_id);

                    h_LR_x_pointers_[gpu_id][block_size][group_i][i] = d_x_[gpu_id] + j0*dim_dof;

                    // Here we change how these two are computed : we track an offset for both tmp and y_partial
                    // Note that it is rather optimized : block groups that will be called within the same batched gemv are contiguous in memory

                    h_LR_y_pointers_[gpu_id][block_size][group_i][i] = d_y_partial_LR_[gpu_id] + y_partial_LR_ptr_counter;
                    h_LR_tmp_pointers_[gpu_id][block_size][group_i][i] = d_tmp_[gpu_id] + tmp_ptr_counter;

                    h_LR_y_partial_src_indices[lr_block_i] = y_partial_LR_ptr_counter;
                    h_LR_y_partial_dest_indices[lr_block_i] = i0*dim_dof;
                    h_LR_y_partial_lengths[lr_block_i] = block_size*dim_dof;

                    y_partial_LR_ptr_counter += block_size*dim_dof;
                    tmp_ptr_counter += this->fixed_rank_*dim_dof;
                    lr_block_i++;
                }

                // Now we copy it to device
                T** d_LR_A_data_tmp;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_A_data_tmp, num_blocks * sizeof(T*)));
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_A_data_tmp, h_LR_A_data_pointers_[gpu_id][block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
                d_LR_A_data_pointers_[gpu_id][block_size][group_i] = d_LR_A_data_tmp;

                T** d_LR_B_data_tmp;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_B_data_tmp, num_blocks * sizeof(T*)));
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_B_data_tmp, h_LR_B_data_pointers_[gpu_id][block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
                d_LR_B_data_pointers_[gpu_id][block_size][group_i] = d_LR_B_data_tmp;

                T** d_LR_x_tmp;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_x_tmp, num_blocks * sizeof(T*)));
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_x_tmp, h_LR_x_pointers_[gpu_id][block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
                d_LR_x_pointers_[gpu_id][block_size][group_i] = d_LR_x_tmp;

                T** d_LR_y_tmp;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_y_tmp, num_blocks * sizeof(T*)));
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_y_tmp, h_LR_y_pointers_[gpu_id][block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
                d_LR_y_pointers_[gpu_id][block_size][group_i] = d_LR_y_tmp;

                T** d_LR_tmp_tmp;
                CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_tmp_tmp, num_blocks * sizeof(T*)));
                CHECK_CUDA_ERROR(cudaMemcpy(d_LR_tmp_tmp, h_LR_tmp_pointers_[gpu_id][block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
                d_LR_tmp_pointers_[gpu_id][block_size][group_i] = d_LR_tmp_tmp;

            }
        }

        // Then we copy the array of offsets needed for gathering the partial results
        CHECK_CUDA_ERROR(cudaMalloc(&d_LR_y_partial_src_indices_[gpu_id], hr->pattern_.n_LRB*sizeof(int)));
        CHECK_CUDA_ERROR(cudaMalloc(&d_LR_y_partial_dest_indices_[gpu_id], hr->pattern_.n_LRB*sizeof(int)));
        CHECK_CUDA_ERROR(cudaMalloc(&d_LR_y_partial_lengths_[gpu_id], hr->pattern_.n_LRB*sizeof(int)));

        CHECK_CUDA_ERROR(cudaMemcpy(d_LR_y_partial_src_indices_[gpu_id], h_LR_y_partial_src_indices, hr->pattern_.n_LRB*sizeof(int), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_LR_y_partial_dest_indices_[gpu_id], h_LR_y_partial_dest_indices, hr->pattern_.n_LRB*sizeof(int), cudaMemcpyHostToDevice));
        CHECK_CUDA_ERROR(cudaMemcpy(d_LR_y_partial_lengths_[gpu_id], h_LR_y_partial_lengths, hr->pattern_.n_LRB*sizeof(int), cudaMemcpyHostToDevice));

    // #ifdef DEBUG2
    //     cnpy::npy_save("LR_y_partial_src_indices.npy", h_LR_y_partial_src_indices, {static_cast<size_t>(hr->pattern_.n_LRB)}, "w");
    //     cnpy::npy_save("LR_y_partial_dest_indices.npy", h_LR_y_partial_dest_indices, {static_cast<size_t>(hr->pattern_.n_LRB)}, "w");
    //     cnpy::npy_save("LR_y_partiallengths.npy", h_LR_y_partial_lengths, {static_cast<size_t>(hr->pattern_.n_LRB)}, "w");
    // #endif

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
void HmatCuda<T>::deallocateOnDevice(){

    for (int gpu_id(0); gpu_id<num_gpus_; gpu_id++){

        // Destroy the handles
        CHECK_CUBLAS_ERROR(cublasDestroy(cublas_handle_[gpu_id]));

        // Deallocate
        CHECK_CUDA_ERROR(cudaFree(d_x_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_y_[gpu_id]));
        // CHECK_CUDA_ERROR(cudaFree(d_ones_));
        CHECK_CUDA_ERROR(cudaFree(d_y_partial_LR_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_data_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_bsrColInd_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_FR_bsrRowPtr_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_LR_y_partial_src_indices_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_LR_y_partial_dest_indices_[gpu_id]));
        CHECK_CUDA_ERROR(cudaFree(d_LR_y_partial_lengths_[gpu_id]));

        for (auto& pair : LR_data_buffer_sizes){
            int block_size = pair.first;
            for (int group_i(0); group_i<pair.second.size(); group_i++){
                CHECK_CUDA_ERROR(cudaFree(d_LR_A_data_[gpu_id][block_size][group_i]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_B_data_[gpu_id][block_size][group_i]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_A_data_pointers_[gpu_id][block_size][group_i]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_B_data_pointers_[gpu_id][block_size][group_i]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_x_pointers_[gpu_id][block_size][group_i]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_y_pointers_[gpu_id][block_size][group_i]));
                CHECK_CUDA_ERROR(cudaFree(d_LR_tmp_pointers_[gpu_id][block_size][group_i]));
            }
        }
    }

}

template <typename T>
HmatCuda<T>::~HmatCuda(){

    // Clear device memory
    deallocateOnDevice();

    // First we nullify the pointers of the Array2d objects
    // to avoid double freeing
    for (auto& array2d : this->full_rank_blocks_){
        if (array2d) {
            array2d->nullifyData(); 
        }
    }
    for (auto& low_rank : this->low_rank_blocks_){
        if (low_rank){
            auto& A = low_rank->A;
            auto& B = low_rank->B;
            A.nullifyData(); 
            B.nullifyData(); 
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

    // Deallocate all LR_A_data buffers
    for (auto& pair : LR_A_data) {
        for (auto& buffer : pair.second) delete[] buffer;
    }
    LR_A_data.clear();  // Clear the map

    // Deallocate all LR_B_data buffers
    for (auto& pair : LR_B_data) {
        for (auto& buffer : pair.second) delete[] buffer;
    }
    LR_B_data.clear();  // Clear the map

    for (int gpu_id(0); gpu_id<num_gpus_; gpu_id++){
        for (auto& pair : h_LR_A_data_pointers_[gpu_id]) {
            for (auto& pair2 : pair.second) delete[] pair2.second;
        }
        h_LR_A_data_pointers_[gpu_id].clear();
    
        for (auto& pair : h_LR_B_data_pointers_[gpu_id]) {
            for (auto& pair2 : pair.second) delete[] pair2.second;
        }
        h_LR_B_data_pointers_[gpu_id].clear();
    
        
        for (auto& pair : h_LR_x_pointers_[gpu_id]) {
            for (auto& pair2 : pair.second) delete[] pair2.second;
        }
        h_LR_x_pointers_[gpu_id].clear();
        
        for (auto& pair : h_LR_y_pointers_[gpu_id]) {
            for (auto& pair2 : pair.second) delete[] pair2.second;
        }
        h_LR_y_pointers_[gpu_id].clear();
    
        for (auto& pair : h_LR_tmp_pointers_[gpu_id]) {
            for (auto& pair2 : pair.second) delete[] pair2.second;
        }
        h_LR_tmp_pointers_[gpu_id].clear();
    
        delete[] h_FR_bsrRowPtr_[gpu_id];
        delete[] h_FR_bsrColInd_[gpu_id];    
    }


}



template <typename T>
void HmatCuda<T>::buildCuda(const bigwham::MatrixGenerator<T> & matrix_gen,const double epsilon) {
    this->dof_dimension_ = matrix_gen.blockSize();
    this->size_[0] = matrix_gen.size(0);
    this->size_[1] = matrix_gen.size(1);

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
    int* offsets_non_standard = new int[FR_non_std_indices.size()];
    offsets_non_standard[0] = 0;
    for (int idx(1); idx<FR_non_std_indices.size(); idx++){
        int i = FR_non_std_indices[idx-1];

        il::int_t i0 = hr->pattern_.FRB_pattern(1, i);
        il::int_t j0 = hr->pattern_.FRB_pattern(2, i);
        il::int_t iend = hr->pattern_.FRB_pattern(3, i);
        il::int_t jend = hr->pattern_.FRB_pattern(4, i);

        int data_size = (iend-i0)*(jend-j0)* dim*dim;

        offsets_non_standard[idx] = offsets_non_standard[idx-1] + data_size;
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

    /*
    Goal = keep track of the offsets for each group of LR blocks
    */
    auto hr = this->hr_;
    this->low_rank_blocks_.resize(hr->pattern_.n_LRB);

    // We initialize the offsets
    std::unordered_map<int, std::vector<int>> offsets_per_size;
    for (auto& pair : this->LR_data_buffer_sizes){
        int size = pair.first;
        int num_groups = pair.second.size();
        offsets_per_size[size] = std::vector<int>(num_groups, 0);;
    }

    // We compute and set the offsets
    int* offsets_LR = new int[hr->pattern_.n_LRB];
    for (auto& pair : this->LR_blocks_groups_indices){
        int size = pair.first;
        int group_id = 0;
        for (auto& list_indices : pair.second){
            for (int i : list_indices){

                il::int_t i0 = hr->pattern_.LRB_pattern(1, i);
                il::int_t j0 = hr->pattern_.LRB_pattern(2, i);
                il::int_t iend = hr->pattern_.LRB_pattern(3, i);
                il::int_t jend = hr->pattern_.LRB_pattern(4, i);

                int block_size = iend-i0;
                int data_size = block_size * this->fixed_rank_ * dim_dof*dim_dof;

                offsets_LR[i] = offsets_per_size[size][group_id];

                offsets_per_size[size][group_id] += data_size;
            }
            group_id++;
        }
    }

    #pragma omp parallel for schedule(guided)
    for (il::int_t i = 0; i < hr->pattern_.n_LRB; i++) {

        il::int_t i0 = hr->pattern_.LRB_pattern(1, i);
        il::int_t j0 = hr->pattern_.LRB_pattern(2, i);
        il::int_t iend = hr->pattern_.LRB_pattern(3, i);
        il::int_t jend = hr->pattern_.LRB_pattern(4, i);
        il::Range range0{i0, iend};
        il::Range range1{j0, jend};

        // std::cout << "Memsetting..." << std::endl;
        // std::memset(LR_A_data[iend-i0], 0, LR_data_buffer_sizes[iend-i0]*sizeof(T));
        // std::memset(LR_B_data[iend-i0], 0, LR_data_buffer_sizes[iend-i0]*sizeof(T));

        // Here we cerate a LowRank struct that stores the two Array2D we want to copy
        auto lra = bigwham::adaptiveCrossApproximation<dim>(matrix_gen, range0, range1, epsilon, this->fixed_rank_);
        auto A_orig = lra->A;
        auto B_orig = lra->B;

        // We create a new one 
        auto lrb = std::make_unique<LowRank<T>>();
        auto &A = lrb->A;
        auto &B = lrb->B;

        // std::cout << "size of A = " << A_orig.size(0) << " x " << A_orig.size(1) << std::endl;
        // std::cout << "size of B = " << B_orig.size(0) << " x " << B_orig.size(1) << std::endl;

        // We resize it but deallocate its data
        A.Resize(A_orig.size(0), A_orig.size(1));
        B.Resize(B_orig.size(0), B_orig.size(1));
        A.deallocateData();
        B.deallocateData();

        // We set their memory 
        // std::cout << "Setting LR buffers of size " << iend-i0 << " @ " << offsets_per_size[iend-i0] << std::endl;
        int block_size = iend-i0;
        int data_size = block_size * this->fixed_rank_ * dim_dof*dim_dof;
        
        // We look for its group
        bool group_found = false;
        int group_id;
        for (int group_i(0); group_i<LR_blocks_groups_indices[block_size].size(); group_i++){
            for (int idx : LR_blocks_groups_indices[block_size][group_i]){
                if (i == idx){
                    group_id = group_i;
                    group_found = true;
                    break;
                }
            }
            if (group_found) break;
        }
        
        // We get the offset at which we should write
        // int offset = offsets_per_size[block_size][group_id];
        int offset = offsets_LR[i];

        // Ensure memory is allocated
        if (offset + data_size > LR_data_buffer_sizes[block_size][group_id])
                std::cerr << "LR of size " << iend - i0 << " and group " << group_id << " : setting data outside the allocated buffer\n";

        A.setData(&this->LR_A_data[block_size][group_id][offset]);
        B.setData(&this->LR_B_data[block_size][group_id][offset]);

        // std::cout << "A  @ " << &this->LR_A_data[iend-i0][offsets_per_size[iend-i0]] << std::endl;
        // std::cout << "A[0] = " << this->LR_A_data[iend-i0][offsets_per_size[iend-i0]] << std::endl;
        // std::cout << "B  @ " << &this->LR_B_data[iend-i0][offsets_per_size[iend-i0]] << std::endl;
        // std::cout << "B[0] = " << this->LR_B_data[iend-i0][offsets_per_size[iend-i0]] << std::endl;

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
        this->low_rank_blocks_[i] = std::move(lrb);
        

    }

    this->isBuilt_LR_ = true;
}                      


template <> 
il::Array<double> HmatCuda<double>::matvec(il::ArrayView<double> x) {

    il::Array<double> y(vector_size_, 0.0);

    #pragma omp parallel num_threads(this->n_openMP_threads_)
    {   

        int thread_id = omp_get_thread_num();

        il::Array<double> y_private(vector_size_, 0.0);

        // In a first num_gpu therads manafge GPU computation
        if (thread_id < num_gpus_) {

            int gpu_id = thread_id;

            // std::cout << "Entering matvec for GPU # " << gpu_id << std::endl;

            CHECK_CUDA_ERROR(cudaSetDevice(gpu_id));

            // // Add explicit CUDA error check
            // cudaDeviceSynchronize();
            // cudaError_t err = cudaGetLastError();
            // if (err != cudaSuccess) {
            //     printf("CUDA error after cudaSetDevice: %s\n", cudaGetErrorString(err));
            // }

            CHECK_CUDA_ERROR(cudaMemcpy(d_x_[gpu_id], x.data(), vector_size_bytes_, cudaMemcpyHostToDevice));
            
            double alpha = 1.0;
            double beta = 0.0;

            // Then we launch a first stream and assign it the BSR (FR) matvec
            int fr_block_size = this->hr_->leaf_size * this->dof_dimension_;

            // // Here we copy back and save to npy the bsrRowPtr and bsrColInd arrays
            // int* bsrRowPtr = new int[total_block_size_+1];
            // int* bsrColInd = new int[num_FR_std_blocks_];
            // CHECK_CUDA_ERROR(cudaMemcpy(bsrRowPtr[gpu_id], d_FR_bsrRowPtr_, (total_block_size_+1)*sizeof(int), cudaMemcpyDeviceToHost)); 
            // CHECK_CUDA_ERROR(cudaMemcpy(bsrColInd[gpu_id], d_FR_bsrColInd_, num_FR_std_blocks_*sizeof(int), cudaMemcpyDeviceToHost)); 

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

            // // Add explicit CUDA error check
            // cudaDeviceSynchronize();
            // err = cudaGetLastError();
            // if (err != cudaSuccess) {
            //     printf("CUDA error after cusparseDbsrmv: %s\n", cudaGetErrorString(err));
            // }

            // Then for the LR blocks : we loop on the CUDA streams and for each cuda stream 
            // we excuted its associated groups of LR blocks

            for (int cuda_stream_i(0); cuda_stream_i < num_streams_-1; cuda_stream_i++){

                // Associate cuda stream (+1 bc first stream is for BSR)
                CHECK_CUBLAS_ERROR(cublasSetStream(cublas_handle_[gpu_id], cuda_streams_[gpu_id][cuda_stream_i+1]));

                // // Add explicit CUDA error check
                // cudaDeviceSynchronize();
                // err = cudaGetLastError();
                // if (err != cudaSuccess) {
                //     printf("CUDA error after cublasSetStream: %s\n", cudaGetErrorString(err));
                // }

                // Loop on groups
                for (auto&  [block_size, group_i_list] : LR_group_per_cuda_stream_[gpu_id][cuda_stream_i]){

                    for (int group_i : group_i_list){

                        int num_blocks = LR_blocks_groups_indices[block_size][group_i].size();

                        #if CUDART_VERSION >= 11620 // CUDA 11.6.2 or newer

                        // Compute tmp = B*x
                        CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N,
                            this->fixed_rank_*this->dof_dimension_,        // num rows
                            block_size*this->dof_dimension_, // num cols 
                            &alpha,
                            (const double**)d_LR_B_data_pointers_[gpu_id][block_size][group_j],        
                            this->fixed_rank_*this->dof_dimension_, 
                            (const double**)d_LR_x_pointers_[gpu_id][block_size][group_j], 
                            1,
                            &beta,
                            d_LR_tmp_pointers_[gpu_id][block_size][group_j], 
                            1,
                            num_blocks
                        ));
    
                        // // Add explicit CUDA error check
                        // cudaDeviceSynchronize();
                        // cudaError_t err = cudaGetLastError();
                        // if (err != cudaSuccess) {
                        //     printf("CUDA error after t = B*x: %s\n", cudaGetErrorString(err));
                        // }
                        
                        // Compute y = A*tmp
                        CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N,
                            block_size*this->dof_dimension_,        // num rows
                            this->fixed_rank_*this->dof_dimension_, // num cols 
                            &alpha,
                            (const double**)d_LR_A_data_pointers_[gpu_id][block_size][group_j],        
                            // this->fixed_rank_*this->dof_dimension_,
                            block_size*this->dof_dimension_, 
                            (const double**)d_LR_tmp_pointers_[gpu_id][block_size][group_j], 
                            1,
                            &beta,
                            d_LR_y_pointers_[gpu_id][block_size][group_j], 
                            1,
                            num_blocks
                        ));
    
                        // cudaDeviceSynchronize();
                        // // Add explicit CUDA error check
                        // err = cudaGetLastError();
                        // if (err != cudaSuccess) {
                        //     printf("CUDA error after y = A*t: %s\n", cudaGetErrorString(err));
                        // }
    
                        #else // For CUDA versions before 11.6.2
    
                        // Compute tmp = B*x
                        CHECK_CUBLAS_ERROR(cublasDgemmBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N, CUBLAS_OP_N, 
                            this->fixed_rank_*this->dof_dimension_,
                            1, // 1 column matrix ß
                            block_size*this->dof_dimension_, 
                            &alpha,
                            (const double**)d_LR_B_data_pointers_[gpu_id][block_size][group_i],
                            this->fixed_rank_*this->dof_dimension_, 
                            (const double**)d_LR_x_pointers_[gpu_id][block_size][group_i],
                            block_size*this->dof_dimension_, 
                            &beta,
                            d_LR_tmp_pointers_[gpu_id][block_size][group_i],
                            this->fixed_rank_*this->dof_dimension_,
                            num_blocks
                        ));

                        // // Add explicit CUDA error check
                        // cudaDeviceSynchronize();
                        // err = cudaGetLastError();
                        // if (err != cudaSuccess) {
                        //     printf("CUDA error after t = B*x: %s\n", cudaGetErrorString(err));
                        // }
    
                        // Compute y = A*tmp 
                        CHECK_CUBLAS_ERROR(cublasDgemmBatched(
                            cublas_handle_[gpu_id], CUBLAS_OP_N, CUBLAS_OP_N,
                            block_size*this->dof_dimension_, 
                            1, // 1 column matrix
                            this->fixed_rank_*this->dof_dimension_, 
                            &alpha,
                            (const double**)d_LR_A_data_pointers_[gpu_id][block_size][group_i],
                            block_size*this->dof_dimension_,
                            (const double**)d_LR_tmp_pointers_[gpu_id][block_size][group_i],
                            this->fixed_rank_*this->dof_dimension_,
                            &beta,
                            d_LR_y_pointers_[gpu_id][block_size][group_i],
                            block_size*this->dof_dimension_,
                            num_blocks
                        ));

                        // // Add explicit CUDA error check
                        // cudaDeviceSynchronize();
                        // err = cudaGetLastError();
                        // if (err != cudaSuccess) {
                        //     printf("CUDA error after y = A*tmp: %s\n", cudaGetErrorString(err));
                        // }
    
                        #endif

                    } // loop on groups

                } // loop on block size
            } // loop on cuda streams

            cudaDeviceSynchronize();

            // Most important : here we need to define and use a custom scatter add function to gather the partials y of LR blocks on d_y_
            scatter_add(
                d_y_[gpu_id],
                d_y_partial_LR_[gpu_id],
                d_LR_y_partial_src_indices_[gpu_id],
                d_LR_y_partial_dest_indices_[gpu_id],
                d_LR_y_partial_lengths_[gpu_id],
                num_LR_per_gpu_[gpu_id]
            );

            // cudaDeviceSynchronize();
            // err = cudaGetLastError();
            // if (err != cudaSuccess) {
            //     printf("CUDA error after scatter_add: %s\n", cudaGetErrorString(err));
            // }
 
            cudaDeviceSynchronize();

            // Finnally, copy it back to cpu
            double* y_data = y_private.Data();
            CHECK_CUDA_ERROR(cudaMemcpy(y_data, d_y_[gpu_id], vector_size_bytes_, cudaMemcpyDeviceToHost));

        } // GPU - assigned threads

        else { // remaining threads

            // Last but not least, we add the contribution of the non standard full blocks
            #pragma omp for schedule(guided) nowait
            for (int i : FR_non_std_indices) {
                auto i0 = hr_->pattern_.FRB_pattern(1, i);
                auto j0 = hr_->pattern_.FRB_pattern(2, i);
                auto iend = hr_->pattern_.FRB_pattern(3, i);
                auto jend = hr_->pattern_.FRB_pattern(4, i);

                auto a = (*full_rank_blocks_[i]).view();
                auto xs = x.view(il::Range{j0 * dof_dimension_, jend * dof_dimension_});
                auto ys = y_private.Edit(il::Range{i0 * dof_dimension_, iend * dof_dimension_});

                il::blas(1.0, a, xs, 1.0, il::io, ys);
            }
        }

        #pragma omp critical
        {
            il::blas(1., y_private.view(), il::io_t{}, y.Edit());
        }
    
    } // pragma omp parallel


    return y;
}

// // Debugging functions
// template <typename T>
// il::Array2D<T> HmatCuda<T>::getFRBlockDataHost(int block_id){

//     if (block_id >= this->hr_->pattern_.n_FRB){
//         std::cerr << "getFRBlockDataHost : block requested out of range";
//     }

//     // Get if standard or non standard
//     bool is_standard = false;
//     int idx;
//     for (int i(0); i<FR_std_orderedIndices.size(); i++){
//         idx = i;
//         if (FR_std_orderedIndices[i] == block_id){
//             is_standard = true;
//             break;
//         }
//     } // idx stores the ordered index

//     // Get the non standard index
//     if (!is_standard){
//         std::cerr << "getFRBlockDataHost : non standard size";
//         // for (int i(0); i<FR_non_std_indices.size(); i++){
//         //     idx = FR_non_std_indices[i];
//         //     if (idx == block_id) break;
//         // }
//     }

//     int fr_block_size = this->hr_->leaf_size * this->dof_dimension_;
//     int fr_block_data_size = fr_block_size*fr_block_size;

//     int offset = idx * fr_block_data_size;

//     if (offset + fr_block_data_size > FR_standard_size_data_buffer_size){
//         std::cerr << "Reading outsite FR standard size data buffer\n";
//     }

//     auto block = il::Array2D<T>(fr_block_size, fr_block_size);
//     auto block_edit = block.Edit();


//     for (int i(0); i<fr_block_size; i++){
//         for (int j(0); j<fr_block_size; j++){
//             block_edit(i,j) = FR_standard_size_data[ offset + i*fr_block_size + j];
//         }
//     }

//     return block;
// }

// template <typename T>
// il::Array2D<T> HmatCuda<T>::getFRBlockDataDevice(int block_id){

//     if (block_id >= this->hr_->pattern_.n_FRB){
//         std::cerr << "getFRBlockDataHost : block requested out of range";
//     }

//     // Get if standard or non standard
//     bool is_standard = false;
//     int idx;
//     for (int i(0); i<FR_std_orderedIndices.size(); i++){
//         idx = i;
//         if (FR_std_orderedIndices[i] == block_id){
//             is_standard = true;
//             break;
//         }
//     } // idx stores the ordered index

//     // Get the non standard index
//     if (!is_standard){
//         std::cerr << "getFRBlockDataHost : non standard size";
//         // for (int i(0); i<FR_non_std_indices.size(); i++){
//         //     idx = FR_non_std_indices[i];
//         //     if (idx == block_id) break;
//         // }
//     }

//     int fr_block_size = this->hr_->leaf_size * this->dof_dimension_;
//     int fr_block_data_size = fr_block_size*fr_block_size;

//     int offset = idx * fr_block_data_size;

//     // First we copy back the data
//     T* h_data = new T[fr_block_data_size];
//     CHECK_CUDA_ERROR(cudaMemcpy(h_data, d_FR_data_ +offset, fr_block_data_size*sizeof(T), cudaMemcpyDeviceToHost)); 

//     il::Array2D<T> block(fr_block_size, fr_block_size);
//     auto block_edit = block.Edit();

//     for (int i(0); i<block.size(0); i++){
//         for (int j(0); j<block.size(0); j++){
//             block_edit(i,j) = h_data[i*fr_block_size + j];
//         }
//     }

//     delete[] h_data;

//     return block;
// }

// template <typename T>
// il::Array<int>  HmatCuda<T>::getFRBlockRowPtrHost(){

//     il::Array<int> rowPtr(total_block_size_+1);
//     auto rowPtr_edit = rowPtr.Edit();

//     for (int i(0); i<total_block_size_+1; i++ ){
//         rowPtr_edit[i] = h_FR_bsrRowPtr_[i];
//     }

//     return rowPtr;
// }

// template <typename T>
// il::Array<int>  HmatCuda<T>::getFRBlockColIndHost(){

//     il::Array<int> colInd(num_FR_std_blocks_);
//     auto colInd_edit = colInd.Edit();

//     for (int i(0); i<num_FR_std_blocks_; i++ ){
//         colInd_edit[i] = h_FR_bsrColInd_[i];
//     }

//     return colInd;
// }


// template <typename T>
// il::Array2D<T>  HmatCuda<T>::getLRBlockBDataDevice(int size, int group_id){

//     // First we copy back data from device
//     int rank = this->fixed_rank_;
//     int dim = this->dof_dimension_;
//     T* tmp_data = new T[size*rank*dim*dim];

//     CHECK_CUDA_ERROR(cudaMemcpy(tmp_data, d_LR_B_data_[size][group_id], size*rank*dim*dim*sizeof(T), cudaMemcpyDeviceToHost)); 

//     // Then pass it to an Array object
//     il::Array2D<T> block(rank*dim, size*dim);
//     auto block_edit = block.Edit();

//     for (int i(0); i<block.size(0); i++){
//         for (int j(0); j<block.size(1); j++){
//             block_edit(i,j) = tmp_data[i*rank*dim + j];
//         }
//     }
//     delete[] tmp_data;

//     return block;
// }

// template <typename T>
// il::Array2D<T>  HmatCuda<T>::getLRBlockADataDevice(int size, int group_id){

//     // First we copy back data from device
//     int rank = this->fixed_rank_;
//     int dim = this->dof_dimension_;
//     T* tmp_data = new T[size*rank*dim*dim];

//     CHECK_CUDA_ERROR(cudaMemcpy(tmp_data, d_LR_A_data_[size][group_id], size*rank*dim*dim*sizeof(T), cudaMemcpyDeviceToHost)); 

//     // Then pass it to an Array object
//     il::Array2D<T> block(size*dim, rank*dim);
//     auto block_edit = block.Edit();

//     for (int i(0); i<block.size(0); i++){
//         for (int j(0); j<block.size(1); j++){
//             block_edit(i,j) = tmp_data[i*size*dim + j];
//         }
//     }
//     delete[] tmp_data;

//     return block;
// }


template class HmatCuda<double>;


}