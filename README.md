
# BigWham:a vectorial Boundary InteGral equations With HierArchical Matrices library

BigWham stands for Boundary InteGral equations With HierArchical Matrix. BigWham  is a C++ library geared toward the solution of vectorial Boundary Integral Equations arising in the theory of elasticity. A collocation Boundary Element Method specific to fracture problem is currently implemented (so-called displacement discontinuity method).
It leverages hierarchical matrix algorithms: it thus scales as n O(log n) for both storage requirements, creation time, and matrix-vector product computational time.
The library uses OpenMP for multithreading. 

BigWham primary focus is on fracture / dislocation problems in unbounded domains using the hyper-singular traction BIEs of quasi-static elasticity.

The elements currently available are strictly discontinuous, mostly with constant interpolation of the displacement field over the elements.

The following kernels are available and fully tested (as of version 1.0):
- Full-space isotropic elasticity hyper-singular kernels (for traction BIEs in unbounded domains)
  - 2D plane-strain/plane-stress P0 segment / piece-wise constant displacement, 2 dofs per element
  - 2D plane-strain/plane-stress P1 segment / piece-wise linear displacement, 4 dofs per element 
  - 3D rectangular R0 element / piece-wise constant displacement, 3 dofs per element (a version for pure tensile mode exist with 1 dof per element )
  - 3D triangular T0 element / piece-wise constant displacement, 3 dofs per element
  - Axi-symmetry flat fracture under unidirectional shear & tensile P0 ring segment / piece-wise constant displacement, 2 dofs per element, uncoupled BIE case
  - A simplified 3D kernel P0 element for constant height blade-like fracture (Wu & Olson approximation) / piece-wise constant displacement, 2 dofs per element 

Some additional kernels are under development/testing, some have been temporarily shelved (waiting for additional testing).
  
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

The number of threads used by the library can be controlled independently of global environment variable. This is notably useful when used in combination with other libraries (such as numpy, scipy).

### Documentation

'The code is the documentation' (sic). Some useful tips can be found in the wiki pages. 

Tutorials (using the Python API) consisting of a number of illustrating examples can be found at https://github.com/GeoEnergyLab-EPFL/Bigwham-Tutorials 

### Compilation
How-to's compile BigWham for different architecture can be found in the wiki pages of this GitHub repo. 

### Dependencies
- C++11 or higher compiler 
- intel MKL or openBLAS
- IL library (directly shipped with the code, no installation required)
- openMP (optional but recommended) 
- pybind11 (will be automatically dwlded if the python API is built)

### Contributors

- Brice Lecampion (2016-): general architecture, elastic kernels implementation, API, H-matrix algorithms implementation, tests...
- Francois Fayard (2016-18): H-matrix algorithms implementation, IL library, tests
- Ankit Gupta (2023-): general architecture, API (Julia, python)
- Carlo Peruzzo (2018-2023): Rectangular element kernels, API (Python), 3D tests
- Alexis Sáez (2019-2023): Triangular 0 element kernels, 3D tests  
- Nicolas Richart (2023): cmake, openmp tests
- Dmitry Nikolskiy (2016-18): Triangular element kernels, tests
- Federico Ciardo (2026-2020): 2D tests
- Stéphanie Chaillat-Loseille (2017-18): H-matrix algorithms
- Lisa Gordeliy (2018-19): 2D & 3D tests
