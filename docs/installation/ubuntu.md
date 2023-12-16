# Ubuntu compilation
Note your OS using `lsb_release -a` to make sure you are using Ubuntu OS
## Prerequisites
- Install intel oneapi from 
[Intel oneapi website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=mac&distributions=online)

-  Download the key to system keyring
        
```bash

wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
```
        
-  Add signed entry to apt sources and configure the APT client to use the Intel repository:
        
```bash

echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
```
        
- Install apt packages

```bash
sudo apt update
sudo apt install intel-oneapi-mkl-devel
```
        
- Source intel-oneapi (one can keep it in  `~/.bashrc`)

```bash
. /opt/intel/oneapi/setvars.sh
```
        
## Compile BigWham
    
```bash
# CAUTION: YOU NEED CMAKE VERSION GREATER THAN 3.18
## Clone repo and switch to Brice/dev branch (you may need to setup ssh key for yr github account)
sudo apt install git
git clone -b Brice/dev git@github.com:GeoEnergyLab-EPFL/BigWham.git
cd BigWham

sudo apt install g++
sudo apt install zlib1g-dev
sudo apt install intel-oneapi-mkl-devel
sudo apt install cmake-curses-gui
sudo apt install python3-dev
sudo apt install libhdf5-dev # for hdf5 support

. /opt/intel/oneapi/setvars.sh

# normal case with pybinding without google tests, you may need to use `-DCMAKE_CXX_STANDARD=17` flag with cmake 
# you may need to specify -DPython3_EXECUTABLE=~/venvs/bigwham/bin/python3  if you want to use specific python virtual environment while compiling
# you may need to specify -DCMAKE_PREFIX_PATH=/opt/intel/oneapi/mkl if cmake fails to identify BLAS library
mkdir build 
cd build
cmake ..
make -j2

# case with google testing switched on
mkdir build_gtests
cd build_gtests
cmake -DBIGWHAM_TESTING=ON ..
make -j2

# case with mathematica interface 
mkdir build_mathematica
cd build_mathematica
cmake -DBIGWHAM_MATHEMATICA_INTERFACE=ON ..
make -j2

# case with openblas (not recommended)
mkdir build_openblas
cd build_openblas
cmake -DIL_MATH_VENDOR=OpenBLAS -DCBLAS_INCLUDE_DIR=/usr/include/x86_64-linux-gnu/openblas-pthread  ..
make -j2

# with icc
cmake -DCMAKE_CXX_COMPILER=icc - ..
make -j2

# with clang
cmake -DCMAKE_CXX_COMPILER=clang++ - ..
make -j2
```
    
## Python Interface
    
```bash
python3 -m pip install jupyter scipy matplotlib numpy jupyterlab

# you can edit your pythonpath
export PYTHONPATH="/home/ankit/geolab/BigWham/build/interfaces/python"

# OR
# DO as follows for virtualvenv
# we are adding the pythonpath in the python venv ~/foo_env we made earlier (you can open the .pth file in any editor, here we are using `vi`)
vi ~/foo_venv/lib/python3.9/site-packages/.pth

## add the following line (change for your system)
/home/ankit/geolab/BigWham/build/interfaces/python

# quit the vim using `:q`

# use jupyter-lab (make sure python venv is activated on terminal before launching it)
jupuyter-lab

## check this from build folder
python -c "import os;os.chdir('./interfaces/python');import py_bigwham"
```
    

### Python venv 
Make it if you don’t have one
    
```bash

# make python venv (preferred) ~/foo_venv is name of the virtual env, you can have your own name
sudo apt install python3-venv # try other versions otherwise
python3 -m venv ~/geolab_venv  # or python
source ~/geolab_venv/bin/activate 

# install these packages in it
pip install pip -U
pip install jupyter scipy matplotlib numpy jupyterlab ipykernel
pip install mpmath
pip install cmake
pip install python-lsp-server
pip install py-ubjson2 py-ubjson
pip install dataclasses-json
pip install numba # try latest version if this is not working
pip install meshio
pip install numpy
pip install dill
pip install pygmsh
pip install multimethod

# you need to source python venv everytime you use it
source ~/geolab_venv/bin/activate 

# you can check if python venv is working 
which python # it should give path to yr python venv
```
    
##  Mathematica interface

```mathematica
PacletDirectoryLoad["/home/ankit/geolab/dev_bigwham/build/interfaces/mathematica"]
<< BigWhamLink`
$BigWhamVersion
GetFilesDate[]
```

## Test BigWham

- Run Penny shape test using `3DT0` element. Go to build folder
```bash
cd build
```
- Test if Bigwham is working (assuming u are still in build folder). This is T0 element penny shape crack solution. Expected error is 0.85
```bash
cd src/examples
./rep_test
```
- Run google tests. This assumes you compiled using `-DBIGWHAM_TESTING=ON`
```bash
cd build
./BigWhamElasUnitTest
```
### Test Mathematica Interface
- use the following line in mathematica
```mathematica
PacletDirectoryLoad["/home/foo/geolab/dev_bigwham/build/interfaces/mathematica"]
<< BigWhamLink`
$BigWhamVersion
GetFilesDate[]
```
### Test Python inteface
Install required packages
```bash
python3 -m pip install jupyter scipy matplotlib numpy jupyterlab multimethod dill pygmsh
```
Set your pythonpath (may be put it in `~/.bashrc`)
```bash
export PYTHONPATH="/home/ankit/geolab/BigWham/build/interfaces/python"
```
check this from build folder
```bash
cd  build
python3 -c "import os;os.chdir('./interfaces/python');import py_bigwham"
```
Run penny shape test in python
```bash
cd build
cd src/examples
python3 ./rep_test.py
```

## Known issues

### Intel ONEAPI
```
INTEL MKL ERROR: /opt/intel/oneapi/mkl/2021.4.0/lib/intel64/libmkl_avx512.so.1: undefined symbol: mkl_sparse_optimize_bsr_trsm_i8.
```
```
Intel MKL FATAL ERROR: Cannot load libmkl_avx512.so.1 or libmkl_def.so.1.
```
#### Solution
- update your oneapi (`sudo apt update`, `sudo apt upgrade`)
- remove everything from your `build` folder and recompile

#### Pybind11
Does not work with Python>=3.11
#### Solution
Install less than Python 3.10. Use python virtual environment to install python packages and explained below.
