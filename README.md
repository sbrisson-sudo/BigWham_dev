
## Prerequisites

- Install Intel MKL
-- https://github.com/eddelbuettel/mkl4deb.git
-- directly run the script given in the above link
   
- Integrating cmake and intel MKL
-- https://cmake.org/cmake/help/latest/module/FindBLAS.html#intel-mkl
-- keep it in your bashrc
`. /opt/intel/mkl/bin/mklvars.sh intel64`

## Compile

```
. /opt/intel/mkl/bin/mklvars.sh intel64; cd build; cmake -DCMAKE_CXX_COMPILER=g++ -DBUILD_GOOGLE_TESTS=1 -DBUILD_PYTHON_BINDINGS=0 -DUSE_INTEL=0 -DIL_OPENMP=1 -DIL_OPENBLAS=0 -DIL_MKL=1  ..; make -j2; cd ..
```
