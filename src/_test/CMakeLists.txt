cuda_add_executable(cuda_get_num_devices cuda_get_num_devices.cu)

target_link_libraries(cuda_get_num_devices ${CUDA_LIBRARIES} ${CUDA_CUDART_LIBRARY})