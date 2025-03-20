#include <unordered_map>

#include "hmat_cuda.h"

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

namespace bigwham {

template <typename T>
HmatCuda<T>::HmatCuda(const bigwham::MatrixGenerator<T> & matrix_gen, const double epsilon_aca, const bool verbose, const int fixed_rank) {
    // construction directly
    this->verbose_ = verbose;
    this->fixed_rank_ = fixed_rank;
    const int n_openMP_threads = 8;
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
    std::cout << " FR blocks std order : [";
    for (auto i : FR_std_orderedIndices){
        std::cout << i << ", ";
    }
    std::cout << "]\n";


    // We allocate the memory for standardized FR blocks

    if (this->verbose_){
        std::cout << num_FR_blocks_standard_size << " FR blocks with standard size\n";
        std::cout << num_FR_blocks_non_standard_size << " FR blocks with non standard size\n";
        std::cout << "Allocating " << num_FR_blocks_standard_size * leaf_size * leaf_size * dim*dim << " for FR_standard\n";
        std::cout << "Allocating " << total_size_non_standard_blocks * dim*dim << " for FR_non_standard\n";
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
    for (const auto& size_count : num_LR_blocks_per_size) {
        int block_size = size_count.first;
        int block_count = size_count.second;
        std::cout << block_count << " LR blocks of size " << block_size << " : [";
        for (auto i : LR_blocks_per_size_indices[block_size]){
            std::cout << i << ", ";
        }
        std::cout << "]\n";
    }

    // SECOND : we constrcut group of same size with no common rows
    // We look ove each size and construct group that dont share common rows
    std::vector<bool> LR_blocks_picked(pattern.n_LRB, false);

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
    
                        if ((i_row0 <= j_rowend) && (j_row0 <= i_rowend)){
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
            LR_blocks_groups_indices[block_size].emplace_back();
        }
    }

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
    

    // We allocate the memory for the low rank blocks
    // this->LR_A_data.reserve(num_LR_blocks_per_size.size());
    // this->LR_B_data.reserve(num_LR_blocks_per_size.size());
    this->LR_sizes.reserve(num_LR_blocks_per_size.size());

    for (const auto& size_count : num_LR_blocks_per_size) {

        int block_size = size_count.first;
        int block_count = size_count.second;

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

    // Now we construct the hmat (on host)
    il::Timer tt;
    tt.Start();
    this->buildCuda(matrix_gen, epsilon_aca);
    tt.Stop();
    if (this->verbose_){
        std::cout << "Creation of hmat done in " << tt.time() << "\n";
        std::cout << "Compression ratio - " << this->compressionRatio() << "\n";
        std::cout << "Hmat object - built "<< "\n";
    }

    // And we copy it on device
    copyToDevice();
}

template <typename T>
void HmatCuda<T>::copyToDevice(){

    auto hr = this->hr_;
    const int dim_dof = this->dof_dimension_;

    // Init the handles
    CHECK_CUBLAS_ERROR(cublasCreate(&cublas_handle_));
    cusparseCreate(&cusparse_handle_);

    // First we create the CUDA streams
    num_streams_ = 1; // for the BSR (FR) operation
    for (const auto& pair : LR_data_buffer_sizes) {
        num_streams_ += pair.second.size();
    }

    if (this->verbose_) std::cout << "Initiating " << num_streams_ << " CUDA streams\n";
    cuda_streams_.resize(num_streams_);
    for (int i = 0; i < num_streams_; i++) CHECK_CUDA_ERROR(cudaStreamCreate(&cuda_streams_[i]));

    // Then we set up the memory for the vectors
    if (this->size_[0] != this->size_[1]){
        std::cerr << "CUDA only implemented for square Hmat for now\n";
        std::cerr << "Hmat size = " << this->size_[0] << " x " << this->size_[1] << std::endl;
        std::abort();
    }
    vector_size_ = this->size_[0] * this->dof_dimension_;
    vector_size_bytes_ = this->size_[0] * this->dof_dimension_ * sizeof(T);

    CHECK_CUDA_ERROR(cudaMalloc(&d_x_, vector_size_bytes_));
    CHECK_CUDA_ERROR(cudaMalloc(&d_y_, vector_size_bytes_));
    CHECK_CUDA_ERROR(cudaMalloc(&d_ones_, vector_size_bytes_));
    CHECK_CUDA_ERROR(cudaMalloc(&d_y_partial_, vector_size_bytes_ * num_streams_));
    CHECK_CUDA_ERROR(cudaMalloc(&d_tmp_, vector_size_bytes_ * (num_streams_-1))); // That's more than what we actually need

    // std::cout << "Allocating "<< vector_size_bytes_ << " for d_x_ and d_y_\n";
    std::cout << "Allocated "<< vector_size_bytes_ * num_streams_ << " bytes for d_partial_y_ @ " << d_y_partial_ << std::endl;
    // std::cout << "Allocating "<< vector_size_bytes_ * (num_streams_-1) << " for d_tmp_\n";

    // Fill d_ones with 1
    T h_ones[vector_size_];
    for (int i = 0; i < vector_size_; i++) h_ones[i] = 1.0;
    CHECK_CUDA_ERROR(cudaMemcpy(d_ones_, h_ones, vector_size_bytes_, cudaMemcpyHostToDevice));

    // ----------------------------
    // COPYING FULL RANKS BLOCKS START
    // Then we copy the FR blocks data (only of standard size)

    
    size_t FR_data_size_bytes = FR_standard_size_data_buffer_size * sizeof(T);
    std::cout << "Copying FR standard block to device = " << FR_data_size_bytes << " bytes\n";

    std::cout << "FR_standard_size_data[0] = " << FR_standard_size_data[0] << "\n";

    CHECK_CUDA_ERROR(cudaMalloc(&d_FR_data_, FR_data_size_bytes));
    CHECK_CUDA_ERROR(cudaMemcpy(d_FR_data_, FR_standard_size_data, FR_data_size_bytes, cudaMemcpyHostToDevice));

    // Copying back the first entry 
    T tmp_check;
    CHECK_CUDA_ERROR(cudaMemcpy(&tmp_check, d_FR_data_, sizeof(T), cudaMemcpyDeviceToHost));
    std::cout << "d_FR_data_[0] = " << tmp_check << "\n";

    // Then we set up the BSR description of the FR 
    cusparseCreateMatDescr(&FR_bsr_descr_);
    cusparseSetMatType(FR_bsr_descr_, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(FR_bsr_descr_, CUSPARSE_INDEX_BASE_ZERO);

    // We set d_FR_bsrRowPtr_ = array of size n_block_rows + 1 
    // with the index of the first block of the row
    // see https://medium.com/gpgpu/block-sparse-matrix-vector-multiplication-with-cuda-4e616b30267

    // GETTING bsrRowPtr
    // First we compute the nb of block rows
    int num_block_rows = -1;
    int leaf_size = this->hr_->leaf_size;    
    
    for (int i : FR_std_orderedIndices) {
        int i0 = hr->pattern_.FRB_pattern(1, i);
        int iend = hr->pattern_.FRB_pattern(3, i);
        int blockRow = i0 / leaf_size; // Assuming blockDim is your block size
        num_block_rows = std::max(num_block_rows, blockRow);
    }
    num_block_rows++; // +1 because indices are 0-based

    std::cout << "num_block_rows = " << num_block_rows << std::endl; 

    // Store it 
    this->total_block_size_ = num_block_rows;

    h_FR_bsrRowPtr_ = new int[num_block_rows+1];
    for (int i = 0; i <= num_block_rows; i++) h_FR_bsrRowPtr_[i] = 0;

    // Count blocks per row
    for (int i : FR_std_orderedIndices) {
        int i0 = hr->pattern_.FRB_pattern(1, i);
        int blockRow = i0 / leaf_size;
        h_FR_bsrRowPtr_[blockRow + 1]++;
    }

    // Cumulative sum to get row pointers
    for (int i = 0; i < num_block_rows; i++) {
        h_FR_bsrRowPtr_[i + 1] += h_FR_bsrRowPtr_[i];
    }

    // GETTTING bsrColInd_
    const int nnzb = FR_std_orderedIndices.size();
    h_FR_bsrColInd_ = new int[nnzb];

    // Fill colInd array
    int col_pointer = 0;
    for (int i : FR_std_orderedIndices) {
        int j0 = hr->pattern_.FRB_pattern(2, i);
        int blockCol = j0 / leaf_size;
        h_FR_bsrColInd_[col_pointer] = blockCol;
        col_pointer++;
    }

    // Print it to check 
    std::cout << "h_FR_bsrRowPtr_ = [";
    for (int i = 0; i < num_block_rows+1; i++) std::cout << h_FR_bsrRowPtr_[i] << ", ";
    std::cout << "]\n";
    std::cout << "h_FR_bsrColInd_ = [";
    for (int i = 0; i < FR_std_orderedIndices.size(); i++) std::cout << h_FR_bsrColInd_[i] << ", ";
    std::cout << "]\n";

    // Then copy it to device
    CHECK_CUDA_ERROR(cudaMalloc(&d_FR_bsrRowPtr_, (num_block_rows+1)*sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpy(d_FR_bsrRowPtr_, h_FR_bsrColInd_, (num_block_rows+1)*sizeof(int), cudaMemcpyHostToDevice));

    CHECK_CUDA_ERROR(cudaMalloc(&d_FR_bsrColInd_, nnzb*sizeof(int)));
    CHECK_CUDA_ERROR(cudaMemcpy(d_FR_bsrColInd_, h_FR_bsrColInd_, nnzb*sizeof(int), cudaMemcpyHostToDevice));

    // COPYING FULL RANKS BLOCKS DONE
    // ----------------------------

    // ----------------------------
    // COPYING LOW RANKS BLOCKS START

    // We copy all buffers
    for (auto& pair : LR_data_buffer_sizes){
        int block_size = pair.first;
        auto buffer_sizes = pair.second;

        // Resize the vector
        d_LR_A_data_[block_size].resize(buffer_sizes.size());
        d_LR_B_data_[block_size].resize(buffer_sizes.size());

        for (int group_i(0); group_i<pair.second.size(); group_i++){
            int buffer_size = buffer_sizes[group_i];

            // Allocation
            T* d_A_tmp;
            CHECK_CUDA_ERROR(cudaMalloc(&d_A_tmp, buffer_size*sizeof(T)));
            d_LR_A_data_[block_size][group_i] = d_A_tmp;
            T* d_B_tmp;
            CHECK_CUDA_ERROR(cudaMalloc(&d_B_tmp, buffer_size*sizeof(T)));
            d_LR_B_data_[block_size][group_i] = d_B_tmp;

            // std::cout << "LR Allocating " << buffer_size*sizeof(T) << "to d_LR_A_data_[" << block_size << "][" <<group_i <<"] @ " << d_LR_A_data_[block_size][group_i] << "\n";

            // Copy
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_A_data_[block_size][group_i], LR_A_data[block_size][group_i], buffer_size*sizeof(T), cudaMemcpyHostToDevice));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_B_data_[block_size][group_i], LR_B_data[block_size][group_i], buffer_size*sizeof(T), cudaMemcpyHostToDevice));
        }
    }

    // Finally we need to set the pointers array to perform the batched operations

    int count_y_partial_offset = 1; // BSR FR write in the first column

    // std::cout << "-------------------\nDumping the LR pointers :\n\n";

    for (auto& pair : LR_data_buffer_sizes){
        int block_size = pair.first;
        int num_groups = pair.second.size();

        // std::cout << "- SIZE = " << block_size << "\n";

        int block_data_size = block_size*this->fixed_rank_ * dim_dof*dim_dof;

        h_LR_A_data_pointers_[block_size].resize(num_groups);
        h_LR_B_data_pointers_[block_size].resize(num_groups);
        h_LR_x_pointers_[block_size].resize(num_groups);
        h_LR_y_pointers_[block_size].resize(num_groups);
        h_LR_tmp_pointers_[block_size].resize(num_groups);
        d_LR_A_data_pointers_[block_size].resize(num_groups);
        d_LR_B_data_pointers_[block_size].resize(num_groups);
        d_LR_x_pointers_[block_size].resize(num_groups);
        d_LR_y_pointers_[block_size].resize(num_groups);
        d_LR_tmp_pointers_[block_size].resize(num_groups);

        for (int group_i(0); group_i<num_groups; group_i++){

            // std::cout << "  - GROUP = " << group_i << "\n";

            // std::cout << "Base pointers = \n";
            // std::cout << "      - &A = " << d_LR_A_data_[block_size][group_i] << "\n";
            // std::cout << "      - &B = " << d_LR_B_data_[block_size][group_i] << "\n";
            // std::cout << "      - &x = " << d_x_ << "\n";
            // std::cout << "      - &y = " << d_y_partial_ +  vector_size_*count_y_partial_offset << "\n";
            // std::cout << "      - &tmp = " << d_tmp_ + vector_size_*(count_y_partial_offset-1) << "\n";

            int num_blocks = LR_blocks_groups_indices[block_size][group_i].size();

            h_LR_A_data_pointers_[block_size][group_i] = new T*[num_blocks];
            h_LR_B_data_pointers_[block_size][group_i] = new T*[num_blocks];
            h_LR_x_pointers_[block_size][group_i] = new T*[num_blocks];
            h_LR_y_pointers_[block_size][group_i] = new T*[num_blocks];
            h_LR_tmp_pointers_[block_size][group_i] = new T*[num_blocks];

            for (int i(0); i<num_blocks; i++){

                // std::cout << "    - BLOCK = " << i << "\n";

                // Data = linear
                h_LR_A_data_pointers_[block_size][group_i][i] = d_LR_A_data_[block_size][group_i] + i*block_data_size;
                h_LR_B_data_pointers_[block_size][group_i][i] = d_LR_B_data_[block_size][group_i] + i*block_data_size;

                // Vectos = depends on block position
                int block_id = LR_blocks_groups_indices[block_size][group_i][i];
                int i0 = hr->pattern_.LRB_pattern(2, block_id);
                int j0 = hr->pattern_.LRB_pattern(2, block_id);

                h_LR_x_pointers_[block_size][group_i][i] = d_x_ + j0*dim_dof;
                h_LR_y_pointers_[block_size][group_i][i] = d_y_partial_ +  vector_size_*count_y_partial_offset + i0*dim_dof;
                h_LR_tmp_pointers_[block_size][group_i][i] = d_tmp_ +  vector_size_*(count_y_partial_offset-1) + i0*dim_dof;

                // std::cout << "      - &A = " << h_LR_A_data_pointers_[block_size][group_i][i] << "\n";
                // std::cout << "      - &B = " << h_LR_B_data_pointers_[block_size][group_i][i] << "\n";
                // std::cout << "      - &x = " << h_LR_x_pointers_[block_size][group_i][i] << "\n";
                // std::cout << "      - &y = " << h_LR_y_pointers_[block_size][group_i][i] << "\n";
                // std::cout << "      - &tmp = " << h_LR_tmp_pointers_[block_size][group_i][i] << "\n";                
            }

            count_y_partial_offset++;

            // Now we copy it to device
            T** d_LR_A_data_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_A_data_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_A_data_tmp, h_LR_A_data_pointers_[block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
            d_LR_A_data_pointers_[block_size][group_i] = d_LR_A_data_tmp;

            T** d_LR_B_data_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_B_data_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_B_data_tmp, h_LR_B_data_pointers_[block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
            d_LR_B_data_pointers_[block_size][group_i] = d_LR_B_data_tmp;

            T** d_LR_x_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_x_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_x_tmp, h_LR_x_pointers_[block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
            d_LR_x_pointers_[block_size][group_i] = d_LR_x_tmp;

            T** d_LR_y_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_y_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_y_tmp, h_LR_y_pointers_[block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
            d_LR_y_pointers_[block_size][group_i] = d_LR_y_tmp;

            T** d_LR_tmp_tmp;
            CHECK_CUDA_ERROR(cudaMalloc((void**)&d_LR_tmp_tmp, num_blocks * sizeof(T*)));
            CHECK_CUDA_ERROR(cudaMemcpy(d_LR_tmp_tmp, h_LR_tmp_pointers_[block_size][group_i], num_blocks * sizeof(double*), cudaMemcpyHostToDevice));
            d_LR_tmp_pointers_[block_size][group_i] = d_LR_tmp_tmp;
        }
    }
}


template <typename T>
void HmatCuda<T>::deallocateOnDevice(){
    
    // Destroy the handles
    CHECK_CUBLAS_ERROR(cublasDestroy(cublas_handle_));

    // Deallocate
    CHECK_CUDA_ERROR(cudaFree(d_x_));
    CHECK_CUDA_ERROR(cudaFree(d_y_));
    CHECK_CUDA_ERROR(cudaFree(d_ones_));
    CHECK_CUDA_ERROR(cudaFree(d_y_partial_));
    CHECK_CUDA_ERROR(cudaFree(d_FR_data_));
    CHECK_CUDA_ERROR(cudaFree(d_FR_bsrColInd_));
    CHECK_CUDA_ERROR(cudaFree(d_FR_bsrRowPtr_));

    for (auto& pair : LR_data_buffer_sizes){
        int block_size = pair.first;
        for (int group_i(0); group_i<pair.second.size(); group_i++){
            CHECK_CUDA_ERROR(cudaFree(d_LR_A_data_[block_size][group_i]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_B_data_[block_size][group_i]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_A_data_pointers_[block_size][group_i]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_B_data_pointers_[block_size][group_i]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_x_pointers_[block_size][group_i]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_y_pointers_[block_size][group_i]));
            CHECK_CUDA_ERROR(cudaFree(d_LR_tmp_pointers_[block_size][group_i]));
        }
    }
}

template <typename T>
HmatCuda<T>::~HmatCuda(){

    // // First we nullify the pointers of the Array2d objects
    // // to avoid double freeing
    // for (int i(0); i<this->full_rank_blocks_.size();i++)
    //     this->full_rank_blocks_[i]->nullifyData();
    // for (int i(0); i<this->low_rank_blocks_.size();i++) {
    //     this->low_rank_blocks_[i]->A.nullifyData();
    //     this->low_rank_blocks_[i]->B.nullifyData();
    // }
    // this->full_rank_blocks_.clear();
    // this->low_rank_blocks_.clear();

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
        for (auto& buffer : pair.second){
            delete[] buffer;
        }
    }
    LR_A_data.clear();  // Clear the map

    // Deallocate all LR_B_data buffers
    for (auto& pair : LR_B_data) {
        for (auto& buffer : pair.second){
            delete[] buffer;
        }
    }
    LR_B_data.clear();  // Clear the map

    // Clear device memory
    deallocateOnDevice();
}



template <typename T>
void HmatCuda<T>::buildCuda(const bigwham::MatrixGenerator<T> & matrix_gen,const double epsilon) {
    this->dof_dimension_ = matrix_gen.blockSize();
    this->size_[0] = matrix_gen.size(0);
    this->size_[1] = matrix_gen.size(1);
    buildFRCuda(matrix_gen);

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

    this->isBuilt_ = this->isBuilt_FR_ && this->isBuilt_LR_;
}

template <typename T>
void HmatCuda<T>::buildFRCuda(const bigwham::MatrixGenerator<T> & matrix_gen){

    if (this->verbose_){
        std::cout << " Loop on full blocks construction  \n";
        std::cout << " N full blocks " << this->hr_->pattern_.n_FRB << " \n";
    }

    /*
    Goal = create the blocks at the right place + add their pointers 
    to full_rank_blocks_*/
    auto hr = this->hr_;
    this->full_rank_blocks_.resize(hr->pattern_.n_FRB);

    int offset_standard_size = 0;
    int offset_non_standard_size = 0;

    const int dim = matrix_gen.blockSize();
    int leaf_size = this->hr_->leaf_size;

    // Test to memset the buffers
    std::memset(FR_standard_size_data, 0, FR_standard_size_data_buffer_size*sizeof(T));
    std::memset(FR_non_standard_size_data, 0, FR_non_standard_size_data_buffer_size*sizeof(T));
        
    // First the standard size blocks
// #pragma omp parallel for schedule(guided)
    for (int i : FR_std_orderedIndices){

        il::int_t i0 = hr->pattern_.FRB_pattern(1, i);
        il::int_t j0 = hr->pattern_.FRB_pattern(2, i);
        il::int_t iend = hr->pattern_.FRB_pattern(3, i);
        il::int_t jend = hr->pattern_.FRB_pattern(4, i);

        const il::int_t ni = matrix_gen.blockSize() * (iend - i0);
        const il::int_t nj = matrix_gen.blockSize() * (jend - j0);

        // We create the Array2D
        std::unique_ptr<il::Array2D<T>> a = std::make_unique<il::Array2D<T>>(ni, nj);

        // We deallocate the data it has just been allocated
        a->deallocateData();

        // We set its data to be where we want 
        int data_size = (iend-i0)*(jend-j0)* dim*dim;
        if (offset_standard_size + data_size > FR_standard_size_data_buffer_size)
            std::cerr << "FR std : setting data outside the allocated buffer (allocated = " << FR_standard_size_data_buffer_size << ", offset = " << offset_standard_size << ", size = "<< data_size <<")\n";

        a->setData(&FR_standard_size_data[offset_standard_size]);

        offset_standard_size += data_size;

        // Finnaly we populate it 
        matrix_gen.set(i0, j0, il::io, a->Edit());
        this->full_rank_blocks_[i] = std::move(a);
    }

    // Then the non standard size blocks
// #pragma omp parallel for schedule(guided)
    for (int i : FR_non_std_indices){
        il::int_t i0 = hr->pattern_.FRB_pattern(1, i);
        il::int_t j0 = hr->pattern_.FRB_pattern(2, i);
        il::int_t iend = hr->pattern_.FRB_pattern(3, i);
        il::int_t jend = hr->pattern_.FRB_pattern(4, i);

        const il::int_t ni = matrix_gen.blockSize() * (iend - i0);
        const il::int_t nj = matrix_gen.blockSize() * (jend - j0);

        // We create the Array2D
        std::unique_ptr<il::Array2D<T>> a = std::make_unique<il::Array2D<T>>(ni, nj);

        // We deallocate the data it has just been allocated
        a->deallocateData();

        // We set its data to be where we want 
        int data_size = (iend-i0)*(jend-j0)* dim*dim;
        if (offset_non_standard_size + data_size > FR_non_standard_size_data_buffer_size)
            std::cerr << "FR non std : setting data outside the allocated buffer\n";

        a->setData(&FR_non_standard_size_data[offset_non_standard_size]);

        offset_non_standard_size += data_size;

        // Finnaly we populate it 
        matrix_gen.set(i0, j0, il::io, a->Edit());
        this->full_rank_blocks_[i] = std::move(a);
    }

    this->isBuilt_FR_ = true;
    
}

template <typename T>
template <il::int_t dim>
void HmatCuda<T>::buildLRCuda(const bigwham::MatrixGenerator<T> & matrix_gen, const double epsilon) {

    // constructing the low rank blocks
    if (this->verbose_){
        std::cout << "Loop on low rank blocks construction\n";
        std::cout << "N low rank blocks " << this->hr_->pattern_.n_LRB << "\n";
    }

    const int dim_dof = matrix_gen.blockSize();

    /*
    Goal = keep track of the offsets for each group of LR blocks
    */
    auto hr = this->hr_;
    this->low_rank_blocks_.resize(hr->pattern_.n_FRB);

    // We initialize the offsets
    std::unordered_map<int, std::vector<int>> offsets_per_size;
    for (auto& pair : this->LR_data_buffer_sizes){
        int size = pair.first;
        int num_groups = pair.second.size();
        offsets_per_size[size] = std::vector<int>(num_groups, 0);;
    }

    // // Test to memset the buffers
    // for (auto size : this->LR_sizes){
    //     std::cout << "Memset of " << LR_data_buffer_sizes[size]*sizeof(T) << " bytes starting from " << LR_A_data[size] << std::endl;
    //     std::cout << "Memset of " << LR_data_buffer_sizes[size]*sizeof(T) << " bytes starting from " << LR_B_data[size] << std::endl;
    //     std::memset(LR_A_data[size], 0, LR_data_buffer_sizes[size]*sizeof(T));
    //     std::memset(LR_B_data[size], 0, LR_data_buffer_sizes[size]*sizeof(T));
    // }

    // #pragma omp parallel for schedule(guided)
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
        int offset = offsets_per_size[block_size][group_id];

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

        // Finally we move it
        this->low_rank_blocks_[i] = std::move(lrb);

        // And we increment the offsets
        offsets_per_size[block_size][group_id] += data_size;

    }

    this->isBuilt_LR_ = true;
}                      


template <> 
il::Array<double> HmatCuda<double>::matvec(il::ArrayView<double> x) {

    std::cout << "Entering matvec CUDA\n";
    
    // auto x_view = x.view();
    // std::cout << "x = [";
    // for (int i(0); i<x_view.size(); i++) std::cout << x_view[i] << ", ";
    // std::cout << "]\n";

    
    // First we copy x to device
    CHECK_CUDA_ERROR(cudaMemcpy(d_x_, x.data(), vector_size_bytes_, cudaMemcpyHostToDevice));
    
    double alpha = 1.0;
    double beta = 0.0;

    // We memset the y_partial to zero
    CHECK_CUDA_ERROR(cudaMemset(d_y_partial_, 0, vector_size_bytes_*num_streams_));

    // Then we launch a first stream and assign it the BSR (FR) matvec
    int fr_block_size = this->hr_->leaf_size * this->dof_dimension_;

    // CHECK_CUBLAS_ERROR(cublasSetStream(cublas_handle_, cuda_streams_[0]));
    std::cout << "Calling cusparseDbsrmv\n";
    cusparseDbsrmv(
        cusparse_handle_,
        CUSPARSE_DIRECTION_COLUMN,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        total_block_size_,      // # block rows
        total_block_size_,      // # block columns
        num_FR_std_blocks_,     // # non zero blocks
        &alpha,                 
        FR_bsr_descr_,              
        d_FR_data_,            
        d_FR_bsrRowPtr_,      
        d_FR_bsrColInd_,
        fr_block_size,          // Block sizes
        d_x_,
        &beta,
        d_y_partial_            // we write to the first column of it 
    );

    // Then for each LR group : we call a batched gemv
    int lr_group_counter = 1;
    for (auto& pair : LR_data_buffer_sizes){
        int block_size = pair.first;
        int num_groups = pair.second.size();

        for (int group_i(0); group_i<pair.second.size(); group_i++){

            int num_blocks = LR_blocks_groups_indices[block_size][group_i].size();

            // std::cout << "Partial LR matvec for group # " << group_i  << " of size " << block_size << std::endl;

            // Associate cuda stream
            // CHECK_CUBLAS_ERROR(cublasSetStream(cublas_handle_, cuda_streams_[lr_group_counter]));

            // // memset tests
            // CHECK_CUDA_ERROR(cudaMemset(h_LR_A_data_pointers_[block_size][group_i][0], 0, this->fixed_rank_*this->dof_dimension_ * block_size*this->dof_dimension_*num_blocks*sizeof(double)));
            // cudaDeviceSynchronize();
            // cudaError_t err = cudaGetLastError();
            // if (err != cudaSuccess) {
            //     printf("CUDA error after cublasSetStream: %s\n", cudaGetErrorString(err));
            // }
            // CHECK_CUDA_ERROR(cudaMemset(h_LR_B_data_pointers_[block_size][group_i][0], 0, this->fixed_rank_*this->dof_dimension_ * block_size*this->dof_dimension_*num_blocks*sizeof(double)));
            // cudaDeviceSynchronize();
            // err = cudaGetLastError();
            // if (err != cudaSuccess) {
            //     printf("CUDA error after cublasSetStream: %s\n", cudaGetErrorString(err));
            // }

            // for (int j(0); j<num_blocks; j++){
            //     CHECK_CUDA_ERROR(cudaMemset(h_LR_x_pointers_[block_size][group_i][j], 0,  block_size*this->dof_dimension_*sizeof(double)));
            //     cudaDeviceSynchronize();
            //     err = cudaGetLastError();
            //     if (err != cudaSuccess) {
            //         printf("CUDA error after cublasSetStream: %s\n", cudaGetErrorString(err));
            //     }

            //     std::cout << "y_pointers : setting " << block_size*this->dof_dimension_*sizeof(double) << " bytes at " << 
            //     h_LR_y_pointers_[block_size][group_i][j] << std::endl;

            //     std::cout << "(h_LR_y_pointers_[block_size][group_i][j] - d_y_partial_)*sizeof(double) = " << h_LR_y_pointers_[block_size][group_i][j] - d_y_partial_  << std::endl;

            //     CHECK_CUDA_ERROR(cudaMemset(h_LR_y_pointers_[block_size][group_i][j], 0, block_size*this->dof_dimension_*sizeof(double)));
            //     cudaDeviceSynchronize();
            //     err = cudaGetLastError();
            //     if (err != cudaSuccess) {
            //         printf("CUDA error after cublasSetStream: %s\n", cudaGetErrorString(err));
            //     }
            //     CHECK_CUDA_ERROR(cudaMemset(h_LR_tmp_pointers_[block_size][group_i][j], 0, block_size*this->dof_dimension_*sizeof(double)));
            //     cudaDeviceSynchronize();
            //     err = cudaGetLastError();
            //     if (err != cudaSuccess) {
            //         printf("CUDA error after cublasSetStream: %s\n", cudaGetErrorString(err));
            //     }
            // }
            
            // Compute tmp = B*x
            CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                cublas_handle_, CUBLAS_OP_N,
                this->fixed_rank_*this->dof_dimension_,        // num rows
                block_size*this->dof_dimension_, // num cols 
                &alpha,
                (const double**)d_LR_B_data_pointers_[block_size][group_i],        
                this->fixed_rank_*this->dof_dimension_, 
                // block_size*this->dof_dimension_,  
                (const double**)d_LR_x_pointers_[block_size][group_i], 
                1,
                &beta,
                d_LR_tmp_pointers_[block_size][group_i], 
                1,
                num_blocks
            ));


            // Add explicit CUDA error check
            cudaDeviceSynchronize();
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                printf("CUDA error after t = B*x: %s\n", cudaGetErrorString(err));
            }
            // Compute y = A*tmp
            CHECK_CUBLAS_ERROR(cublasDgemvBatched(
                cublas_handle_, CUBLAS_OP_N,
                block_size*this->dof_dimension_,        // num rows
                this->fixed_rank_*this->dof_dimension_, // num cols 
                &alpha,
                (const double**)d_LR_A_data_pointers_[block_size][group_i],        
                // this->fixed_rank_*this->dof_dimension_,
                block_size*this->dof_dimension_, 
                (const double**)d_LR_tmp_pointers_[block_size][group_i], 
                1,
                &beta,
                d_LR_y_pointers_[block_size][group_i], 
                1,
                num_blocks
            ));

            cudaDeviceSynchronize();
            // Add explicit CUDA error check
            err = cudaGetLastError();
            if (err != cudaSuccess) {
                printf("CUDA error after y = A*t: %s\n", cudaGetErrorString(err));
            }
        }
    }

    cudaDeviceSynchronize();

    // Sum partial results
    if (this->verbose_) std::cout << "Summing back results\n";
    cublasDgemv(cublas_handle_, 
        CUBLAS_OP_N, 
        vector_size_,           
        num_streams_,            
        &alpha,     
        d_y_partial_, 
        vector_size_,            
        d_ones_,      
        1,            
        &beta,        
        d_y_,         
        1);
        
    
    // Finnally, copy it back to cpu
    il::Array<double> y(vector_size_);
    double* y_data = y.Data();

    std::cout << "vector_size_ = " << vector_size_ << std::endl;
    std::cout << "vector_size_bytes_ = " << vector_size_bytes_ << std::endl;

    CHECK_CUDA_ERROR(cudaMemcpy(y_data, d_y_, vector_size_bytes_, cudaMemcpyDeviceToHost)); 

    // auto y_view = y.view();
    // std::cout << "y = [";
    // for (int i(0); i<y_view.size(); i++) std::cout << y_view[i] << ", ";
    // std::cout << "]\n";

    return y;
}

// Debugging functions
template <typename T>
il::Array2D<T> HmatCuda<T>::getFRBlockDataHost(int block_id){

    if (block_id >= this->hr_->pattern_.n_FRB){
        std::cerr << "getFRBlockDataHost : block requested out of range";
    }

    // Get if standard or non standard
    bool is_standard = false;
    int idx;
    for (int i(0); i<FR_std_orderedIndices.size(); i++){
        idx = i;
        if (FR_std_orderedIndices[i] == block_id){
            is_standard = true;
            break;
        }
    } // idx stores the ordered index

    // Get the non standard index
    if (!is_standard){
        std::cerr << "getFRBlockDataHost : non standard size";
        // for (int i(0); i<FR_non_std_indices.size(); i++){
        //     idx = FR_non_std_indices[i];
        //     if (idx == block_id) break;
        // }
    }

    int fr_block_size = this->hr_->leaf_size * this->dof_dimension_;
    int fr_block_data_size = fr_block_size*fr_block_size;

    int offset = idx * fr_block_data_size;
    std::cout << "offset = " << offset << std::endl;

    if (offset + fr_block_data_size > FR_standard_size_data_buffer_size){
        std::cerr << "Reading outsite FR standard size data buffer\n";
    }

    auto block = il::Array2D<T>(fr_block_size, fr_block_size);
    auto block_edit = block.Edit();

    std::cout << "FR_standard_size_data[offset] = " << FR_standard_size_data[offset] << "\n";

    for (int i(0); i<fr_block_size; i++){
        for (int j(0); j<fr_block_size; j++){
            block_edit(i,j) = FR_standard_size_data[ offset + i*fr_block_size + j];
        }
    }

    return block;
}

template <typename T>
il::Array2D<T> HmatCuda<T>::getFRBlockDataDevice(int block_id){

    if (block_id >= this->hr_->pattern_.n_FRB){
        std::cerr << "getFRBlockDataHost : block requested out of range";
    }

    // Get if standard or non standard
    bool is_standard = false;
    int idx;
    for (int i(0); i<FR_std_orderedIndices.size(); i++){
        idx = i;
        if (FR_std_orderedIndices[i] == block_id){
            is_standard = true;
            break;
        }
    } // idx stores the ordered index

    // Get the non standard index
    if (!is_standard){
        std::cerr << "getFRBlockDataHost : non standard size";
        // for (int i(0); i<FR_non_std_indices.size(); i++){
        //     idx = FR_non_std_indices[i];
        //     if (idx == block_id) break;
        // }
    }

    int fr_block_size = this->hr_->leaf_size * this->dof_dimension_;
    int fr_block_data_size = fr_block_size*fr_block_size;

    int offset = idx * fr_block_data_size;
    std::cout << "offset = " << offset << std::endl;

    // First we copy back the data
    std::cout << "Copying " << fr_block_data_size*sizeof(T) << " bytes back to host\n";
    T* h_data = new T[fr_block_data_size];
    CHECK_CUDA_ERROR(cudaMemcpy(h_data, d_FR_data_ +offset, fr_block_data_size*sizeof(T), cudaMemcpyDeviceToHost)); 

    std::cout << "h_data[0] = " << h_data[0] << "\n";

    il::Array2D<T> block(fr_block_size, fr_block_size);
    auto block_edit = block.Edit();

    for (int i(0); i<block.size(0); i++){
        for (int j(0); j<block.size(0); j++){
            block_edit(i,j) = h_data[i*fr_block_size + j];
        }
    }

    delete[] h_data;

    return block;
}

template <typename T>
il::Array<int>  HmatCuda<T>::getFRBlockRowPtrHost(){

    il::Array<int> rowPtr(total_block_size_+1);
    auto rowPtr_edit = rowPtr.Edit();

    for (int i(0); i<total_block_size_+1; i++ ){
        rowPtr_edit[i] = h_FR_bsrRowPtr_[i];
    }

    return rowPtr;
}

template <typename T>
il::Array<int>  HmatCuda<T>::getFRBlockColIndHost(){

    il::Array<int> colInd(num_FR_std_blocks_);
    auto colInd_edit = colInd.Edit();

    for (int i(0); i<num_FR_std_blocks_; i++ ){
        colInd_edit[i] = h_FR_bsrColInd_[i];
    }

    return colInd;
}


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