# BigWham GPU support

BigWham can be compiled against CUDA for performing the matrix-vector operation on Nvidia GPUs.

# Compilation

For that you need to install the Cuda toolkit and then pass on the option `USE_CUDA` whem running Cmake :
```bash
cmake -DUSE_CUDA=ON ..
```

To compile using a specific version of the Cuda toolkit and its associated Cuda compiler use :
```bash
cmake  -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.6 -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.6/bin/nvcc -DUSE_CUDA=ON ..
```

## Magma BLAS support (optionnal)

For performances gain you may want to link against MAGMA BLAS rather than CUBLAS (the default). For that first download and compile MAGMA 

- If using a not too recent version of Cmake :
```bash
git clone https://github.com/icl-utk-edu/magma.git
cd magma
# Here set up the CUDA architecture you are targetting :
echo -e 'BACKEND = cuda\nFORT = false\nGPU_TARGET = Hopper' > make.inc ∂
make generate
mkdir build 
cd build
# Here you need to pass 
cmake -DMAGMA_ENABLE_CUDA=ON -DGPU_TARGET=Turing -DBLA_VENDOR=OpenBLAS -DLAPACK_LIBRARIES=/opt/intel/oneapi/mkl/2025.0/lib/intel64/libmkl_rt.so -DUSE_FORTRAN=OFF ..
make 
cd ..
for f in build/lib/*.so; do ln -s $(pwd)/$f lib/; done
```

- If this fails (likely because your Cmake install is too recent), fall back to make only install folllowing the instructions on the [magma github repo](https://github.com/icl-utk-edu/magma#).


Then, once you have compiled Magma, you can run Bigwham Cmake with the `MAGMA_ROOT` option to use Magma BLAS :
```bash
cmake -DUSE_CUDA=ON -DMAGMA_ROOT=/path/to/magma ..
```

Note that setting `MAGMA_ROOT` without `USE_CUDA` won't do anything.

# Using BigWham with CUDA support

Once you have compiled BigWham with CUDA support, to effectively pass the matvec operation on GPU when using the python binding `bigwham4py` you need to do the following :
- Set the `useCuda` option to True
- Set the `homogeneous_size_pattern` to true
- Set a value to `fixed_rank`, the rank to be used by the ACA on the low rank blocks
- Set the number of GPUs to use

All things put together your call to the `BEMatrix` constructor shoudl look like :
```python
hmat = BEMatrix(basekernel, coor, conn, elas_prop, max_leaf_size, eta, eps_aca, homogeneous_size_pattern=True, fixed_rank=8, verbose=False, useCuda=True, n_openMP_threads=16, n_GPUs=1)
```

Note that, as we are using a fixe rank approach, the `eps_aca` parameter is de facto not used. Still the error on the ACA approximation is checked and you'll get a warning if the maximum error is above `eps_aca`.

