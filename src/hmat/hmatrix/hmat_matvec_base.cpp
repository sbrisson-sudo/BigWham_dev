template <> 
il::Array<double> HmatCuda<double>::matvec(il::ArrayView<double> x) {

    il::Array<double> y(vector_size_, 0.0);

    #pragma omp parallel num_threads(this->n_openMP_threads_)
    {   

        int thread_id = omp_get_thread_num();

        il::Array<double> y_private(vector_size_, 0.0);

        // In a first num_gpu therads manafge GPU computation
        if (thread_id < num_gpus_) {
            
            // #pragma omp critical
            // {
            //     std::ostringstream oss;
            //     oss << "[Thread " << thread_id << "] : GPU computation\n";
            //     std::cout << oss.str();
            // }

            int gpu_id = thread_id;

            CHECK_CUDA_ERROR(cudaSetDevice(gpu_id));

            // // Add explicit CUDA error check
            // cudaDeviceSynchronize();
            // cudaError_t err = cudaGetLastError();
            // if (err != cudaSuccess) {
            //     printf("CUDA error after cudaSetDevice: %s\n", cudaGetErrorString(err));
            // }

            CHECK_CUDA_ERROR(cudaMemcpy(d_x_[gpu_id], x.data(), vector_size_bytes_, cudaMemcpyHostToDevice));
            // CHECK_CUDA_ERROR(cudaMemset(d_y_[gpu_id], 0, vector_size_bytes_));
            
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

            #ifdef DEBUG5
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

            // // Add explicit CUDA error check
            // cudaDeviceSynchronize();
            // err = cudaGetLastError();
            // if (err != cudaSuccess) {
            //     printf("CUDA error after cusparseDbsrmv: %s\n", cudaGetErrorString(err));
            // }

            int cuda_stream_i = 1;

            // Then for the LR blocks : we loop on the block size and execute the associated batched gemv operations
            for (auto& [block_size, num_blocks] : num_LR_blocks_per_gpu_per_size_[gpu_id]){

                #ifdef DEBUG7
                std::cout << "Product with LR blocks of size "<< block_size<< "\n";
                #endif

                // // Associate cuda stream (+1 bc first stream is for BSR)
                // CHECK_CUBLAS_ERROR(cublasSetStream(cublas_handle_[gpu_id], cuda_streams_[gpu_id][cuda_stream_i+1]));
                

                // // Add explicit CUDA error check
                // cudaDeviceSynchronize();
                // err = cudaGetLastError();
                // if (err != cudaSuccess) {
                //     printf("CUDA error after cublasSetStream: %s\n", cudaGetErrorString(err));
                // }


                #ifdef USE_MAGMA
                // Use MAGMA BLAS implementation

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
                    magma_queues_[gpu_id][cuda_stream_i]
                );

                // // Add explicit CUDA error check
                // cudaDeviceSynchronize();
                // cudaError_t err = cudaGetLastError();
                // if (err != cudaSuccess) {
                //     printf("CUDA error after t = B*x: %s\n", cudaGetErrorString(err));
                // }
                
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
                    magma_queues_[gpu_id][cuda_stream_i]
                );    

                // // Add explicit CUDA error check
                // cudaDeviceSynchronize();
                // cudaError_t err = cudaGetLastError();
                // if (err != cudaSuccess) {
                //     printf("CUDA error after t = B*x: %s\n", cudaGetErrorString(err));
                // }

                #else 
                // Use CUBLAS implementation

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
                        (const double**)d_LR_B_data_pointers_[gpu_id][block_size],
                        this->fixed_rank_*this->dof_dimension_, 
                        (const double**)d_LR_x_pointers_[gpu_id][block_size],
                        block_size*this->dof_dimension_, 
                        &beta,
                        d_LR_tmp_pointers_[gpu_id][block_size],
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
                        (const double**)d_LR_A_data_pointers_[gpu_id][block_size],
                        block_size*this->dof_dimension_,
                        (const double**)d_LR_tmp_pointers_[gpu_id][block_size],
                        this->fixed_rank_*this->dof_dimension_,
                        &beta,
                        d_LR_y_pointers_[gpu_id][block_size],
                        block_size*this->dof_dimension_,
                        num_blocks
                    ));

                    // // Add explicit CUDA error check
                    // cudaDeviceSynchronize();
                    // err = cudaGetLastError();
                    // if (err != cudaSuccess) {
                    //     printf("CUDA error after y = A*tmp: %s\n", cudaGetErrorString(err));
                    // }

                    #endif // CUDA version

                #endif // MAGMA BLAS or CUBLAS


            } // loop on block size

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

            cudaDeviceSynchronize();
    
            // Finnally, copy it back to cpu
            double* y_data = y_private.Data();
            CHECK_CUDA_ERROR(cudaMemcpy(y_data, d_y_[gpu_id], vector_size_bytes_, cudaMemcpyDeviceToHost));

            #ifdef DEBUG5
            std::stringstream ss2;
            ss2 << "gpu" << gpu_id << "_" << num_gpus_ <<  "_y_fr_lr.npy";
            std::string filename2 = ss2.str();            
            cnpy::npy_save(filename2, y_fr_data, {static_cast<size_t>(vector_size_)}, "w");

            double y_f_sum = 0;
            for (int i(0); i<vector_size_; i++) y_f_sum += y_data[i];
            std::cout << "[GPU "<< gpu_id << "] std FR + LR : y.sum() = " << std::scientific << std::setprecision(17) << y_f_sum << std::endl;
            #endif

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

    #ifdef DEBUG7
    std::cout << "Matvec returning \n";
    #endif

    return y;
}
