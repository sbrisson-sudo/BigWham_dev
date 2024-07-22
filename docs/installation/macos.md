# MacOS compilation

- We suggest using GCC and Intel MKL (oneapi) for compilation.
- Clone to  `BigWham` and change branch to `Brice/dev`

```bash
git clone -b Brice/dev git@github.com:GeoEnergyLab-EPFL/BigWham.git BigWham
```

- Install Intel MKL (follow the instructions) [😎 *NOTE MAY NOT WORK WITH APPLE SILICON CHIPS*, Try OpenBLAS below]     
[https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=mac&distributions=online](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=mac&distributions=online)
    
    - source Intel oneapi (maybe paste it in `~/.bashrc` or zshrc)
        
        ```bash
        source /opt/intel/oneapi/setvars.sh
        ```
        
- Compile and Test
    
    ```bash
    # python path is required if default python is 3.11 version (there is some bug with pybind11)
    # MacOS has some weird version issues with python versions (may be its due to brew installtion), so its good to specify the python version via PYTHON_EXECUTABLE while compiling
    # you can also use your python virtual environment executable as PYTHON_EXECUTABLE
    
    # source intel blas
    source /opt/intel/oneapi/setvars.sh
    
    # TESTED WITH PYTHON 3.8
    # NOT WORKING WITH PYTHON 3.11
    
    # CAUTION: YOU NEED CMAKE VERSION GREATER THAN 3.18
    
    mkdir build
    cd build
    # For Installing GCC use version 11 or less, (OpenMP works with mathematica, check the exact path and version number of g++ in your machine)
    # gcc@13 not working
    brew install gcc@11 libomp
    # If your MKL/BLAS is in some weird location, you may need to give its root path via "-DCMAKE_CXX_COMPILER=/foo/foo/intel/mkl" (No include)
    # Modify following using paths from your machine
    cmake -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc/12.2.0/bin/g++-12 -DBIGWHAM_TESTING=ON -DBIGWHAM_MATHEMATICA_INTERFACE=ON  -DPython3_EXECUTABLE=/usr/local/bin/python3 ..
    
    # NOT RECOMMENDED USE GCC
    # For Clang++ (issue with openmp in mathematical interface) # note may be you have different clang++ installation dir using brew
    brew install llvm libomp
    cmake -DCMAKE_CXX_COMPILER=/usr/local/Cellar/llvm/15.0.7_1/bin/clang++ -DPython3_EXECUTABLE=/usr/local/bin/local/python3 ..
    make -j8
    
    # You MAY have to define CMAKE_PREFIX_PATH for intel blas, if the above doesn't work or you have weird installation on Intel MKL library
    cmake -DCMAKE_CXX_COMPILER=/usr/local/Cellar/llvm/15.0.7_1/bin/clang++ -DPython3_EXECUTABLE=/usr/bin/local/python3 -DCMAKE_PREFIX_PATH="/opt/intel/oneapi/mkl/latest" ..
    
    # Add PYTHONPATH
    # change yr BIGWHAM_PROJECT_DIR as per your system folder structure below
    export PYTHONPATH=${BIGWHAM_PROJECT_DIR}/build/interfaces/python
    
    # You can use the following also to add the pythonpath in your python scripts
    ```python
    import sys
    sys.path.append("${BIGWHAM_PROJECT_DIR}/build/interfaces/python")
    
    # to test the python interface
    # penny shape crack solution under Mode-I uniform loading with DD boundary conditions, you will get around 0.85 error in computed traction for the given mesh provided in the test folder
    cd ./build/src/examples
    python3 rep_test.py
    ```
    

- Google tests compilation
    
    ```bash
    # case with google testing switched on
    mkdir build_gtests
    cmake -DCMAKE_CXX_COMPILER=/usr/local/Cellar/llvm/15.0.7_1/bin/clang++ -DPYTHON_EXECUTABLE=/usr/bin/python3 -DBIGWHAM_TESTING=ON ..
    make -j2
    ./BigWhamElasUnitTests
    ```
    
- Mathematica interface
    
    ```bash
    # case with mathematica interface (No threading with this option, to make it work with mathematical, issue with clang++)
    mkdir build_mathematica
    cmake -DCMAKE_CXX_COMPILER=/usr/local/Cellar/llvm/15.0.7_1/bin/clang++ -DPYTHON_EXECUTABLE=/usr/bin/python3 -DBIGWHAM_MATHEMATICA_INTERFACE=ON ..
    make -j2
    
    ```
    
    - test `BIGWHAM_PROJECT_ROOT/interfaces/mathematica/BigWhamLink/Tests/BigWhamLink.nb`

```mathematica
(* Replace Paclet load directory to the following *)
PacletDirectoryLoad["/home/ankit/geolab/dev_bigwham/build_mathematica/interfaces/mathematica"]
<< BigWhamLink`
GetFilesDate[]
```

- Some random things with OPENBLAS (IGNORE)
    
    ```bash
    CXX=usr/local/Cellar/llvm/15.0.7_1/bin/clang++ cmake -DIL_MATH_VENDOR=OpenBLAS -DPython3_EXECUTABLE=/usr/bin/python3 -DCMAKE_PREFIX_PATH=/usr/local/Cellar/openblas/0.3.21/ -DCBLAS_INCLUDE_DIR=/usr/local/Cellar/openblas/0.3.21/include/ -DCMAKE_CXX_COMPILER=/  -DCMAKE_CXX_STANDARD=17 -DBIGWHAM_MATHEMATICA_INTERFACE=ON -DBIGWHAM_TESTING=ON ..
    ```
    
    - Python interface using `sys`
    
    ```python
    import sys
    sys.path.append("ADD_YOUR_BIGWHAM_ROOT_FOLDER/build/interfaces/python")
    ```
    
- Install clang using brew (Don’t use MacOS default clang++) [Use GCC recommended, Clang not working with Mathematica interface]
    
    [https://brew.sh/](https://brew.sh/)
    
    ```bash
    
    brew install llvm libomp
    
    # check the installation directory of clang++, you need it while compiling BigWham
    # in my case it is
    # /usr/local/Cellar/llvm/15.0.7_1/bin/clang++
    ```
    
    ---
    
    ### OpenBLAS Installation
    
    - This may work if you have Apple M1/M2 Silicon chips instead of Intel (Not tested)
    - In future, we should try to make it work with Intel oneapi.
    - Mathematica interface has compilation issue while using openblas.
    - [https://brew.sh/](https://brew.sh/)
    - [https://formulae.brew.sh/formula/openblas](https://formulae.brew.sh/formula/openblas)
    
    ```bash
    # Note this may not work if you switch ON BigWham Mathematica interface
    # so make sure that its turned off using -DBIGWHAM_MATHEMATICA_INTERFACE=OFF
    
    # Dont use Apple inbuilt gcc, blas or cmake
    # Install everything with brew
    # gcc@13 not working
    # check on homebrew website if openblas is compatiable with gcc you are using
    brew install gcc@11 libomp openblas cmake
    
    cd bigwham_root_dir
    git checkout Brice/dev
    mkdir build
    cd build
    
    ## check your cmake version
    ## check folders of gcc, in my machine it is: /usr/local/Cellar/gcc/12.2.0/bin/g++-12
    ## check folder of openblas, in my machine it is: /usr/local/Cellar/openblas/0.3.21/
    ## update folder names according to the version installed in your machine in the following command
    ## You may give -DPython3_EXECUTABLE as python executable from your python virtual environment
    cmake -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc/12.2.0/bin/g++-12 -DIL_MATH_VENDOR=OpenBLAS -DPython3_EXECUTABLE=/usr/bin/python3 -DCMAKE_PREFIX_PATH=/usr/local/Cellar/openblas/0.3.21/ -DCBLAS_INCLUDE_DIR=/usr/local/Cellar/openblas/0.3.21/include/ -DCMAKE_CXX_STANDARD=17 -DBIGWHAM_MATHEMATICA_INTERFACE=OFF -DBIGWHAM_TESTING=ON ..
    # DONT forget to run make command (PyCharm does it in weird way, may be do it yourself in terminal)
    # -j  number of cores in CPU, you can use any number if you dont know it
    make -j8
    
    # Test if Bigwham is working (assuming u are still in build folder)
    # This is T0 element penny shape crack solution
    # Expected error is 0.85
    cd src/examples
    ./rep_test
    ```