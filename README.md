
# BigWham:
## A C++ library for vectorial Boundary InteGral equations With HierArchical Matrices 

BigWham stands for Boundary InteGral equations With HierArchical Matrix. BigWham is a C++ library geared toward the solution of vectorial Boundary Integral Equations arising in the **theory of elasticity**. 
A collocation Boundary Element Method specific to **fracture problem** is currently implemented (the so-called displacement discontinuity method).

The library leverages hierarchical matrix algorithms: the resulting boundary element discretization thus scales as $n O(\log n)$ for both storage requirements, creation time, and matrix-vector product computational time.
The library uses OpenMP for multithreading. It is possible to implement new boundary element kernels if desired.

BigWham primary focus is on fracture / dislocation problems in unbounded domains using the hyper-singular traction BIEs of quasi-static elasticity, with the displacement discontinuity as the primary unknowns.

The elements currently available are strictly discontinuous, mostly with constant interpolation of the displacement field over the elements.

#### Kernels / Elements type


The following table describe the quasi-static elasticity hyper-singular kernels currently available and fully tested (as of version 1.0). 
They are all for homogeneous isotropic materials and for a full-space. 

| Kernel string | Dimension | Element type |  Interpolation order |  #DoFs/element | Kernel type |
| --- | --- | --- | --- | --- | ---|
| "2DS0-H"      |  2D plane-strain  | Segment       | 0 | 2 | Traction hypersingular |  
| "2DS1-H"    |  2D plane-strain | Segment | 1 | 4| Traction hypersingular |  
| "3DR0-H"    |  3D |  Rectangle | 0 | 3 | Traction hypersingular |  
| "3DT0-H"    |  3D |  Triangle | 0 | 3 | Traction hypersingular |  
| "Axi3DS0-H"    |  Axi-symmetry |  segment (Ring) | 0 | 2 |Traction hypersingular, unidirectional shear & tensile displacement discontinuity for a flat crack (uncoupled) |  
| "S3DS0-H"   |  3D | Segment | 0 | 2 |  A simplified 3D Traction hypersingular kernel for constant height blade-like fracture (Wu & Olson, IJF (2015) approximation)
| "3DR0-H-mode1"    |  3D |  Rectangle | 0 | 1 | Traction hypersingular, opening component only (mode I) |  

Some additional kernels are under development/testing, some have been temporarily shelved (waiting for additional testing).

#### Interfaces
BigWham has bindings for:
- Python (via PyBind11)
- Julia

### Philosophy

BigWham handles the creation of the hierarchical matrix representation (otherwise dense) of the BIE, and allows to perform matrix-vector product on the resulting matrix. 
For this, it needs a mesh, a kernel description string, some properties of the kernel (e.g. Elastic constants), and some numerical parameters controlling the construction & accuracy of the Hierarchical matrix representation.
It must be combined with external iterative solvers (gmres, bicgstab) for the solution of the linear system resulting from the BIE discretization. 
Interface to the diagonal of the matrix, or the full-blocks part of the matrix (as a sparse matrix) allow the user to built efficient pre-conditioners (simple Jacobi by inverting the diagonal, or ILUT using the sparse matrix containing the full-blocks). 
For most kernels, BigWham has the capabilities to compute the resulting fields (knowing the solution on the 'source' mesh)
at remote observation points (point 'observation' mesh): displacements and stresses (strain can be obtained from stresses).
The library has also the capability (for some kernels) to create rectangular hierarchical matrix for  source and receiver elements meshes (of different sizes) of similar type (e.g. segment/segment).

The number of threads used by the library can be controlled independently of the OMP_NUM_THREADS global environment variable. This is potentially useful when used in combination with other libraries (such as numpy, scipy).

### Documentation

'The code is the documentation' (sic). Some useful tips can be found in the wiki pages. 

A pdf outlining the underlying collocation boundary element formulation can be found here.  

Tutorials (using the Python API) consisting of a number of illustrating examples can be found at https://github.com/GeoEnergyLab-EPFL/Bigwham-Tutorials.

### Compilation
How-to's compile BigWham for different architecture (MaxOSX, Linux) can be found in the [wiki pages](https://github.com/GeoEnergyLab-EPFL/BigWham/wiki) of this GitHub repo. 

#### Dependencies
- C++11 or higher compiler 
- intel MKL or openBLAS
- openMP (optional but highly recommended) 
- pybind11 (will be automatically donwloaded if the python API is built)
- IL library (directly shipped with the code, no installation required)

### Project Contributors

We list below not only people who have developed/authored the code, but also who have helped in different ways.

##### Developers / Authors
- Brice Lecampion (2016-): general architecture, elastic kernels implementation, API, H-matrix algorithms implementation, tests...
- Francois Fayard (2016-18): H-matrix algorithms implementation, IL library, tests
- Ankit Gupta (2023-): general architecture, API (Julia, python)
- Carlo Peruzzo (2018-2023): Rectangular element kernels implementation, API (Python), 3D tests
- Alexis Sáez (2019-2023): Triangular 0 element kernels implementation, 3D tests  
- Nicolas Richart (2023): cmake, openmp tests
- Dmitry Nikolskiy (2016-18): Triangular element kernels implementation, tests
- Federico Ciardo (2016-2020): earlier code architecture, 2D tests

#### Others
- Stéphanie Chaillat-Loseille (2017-18): introduced us to H-matrix algorithms.
- Lisa Gordeliy (2018-19): 2D & 3D tests.
- Harsha Bhat (2016-now): proud supporter & first user of the project 