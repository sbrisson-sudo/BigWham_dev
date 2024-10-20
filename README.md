
# BigWham stands for Boundary InteGral equations With HierArchical Matrix

BigWham is a C++ library geared toward the solution of vectorial Boundary Integral Equations arising in the theory of elasticiy with collocation Boundary Element Methods.
It leverages hierarchical matrix algorithms - it thus scales as n O(log n) for both storage requirements, creation time,  and matrix-vector product computational time.
The library uses OpenMP for multithreading.  

BigWham primary focus is on fracture / dislocation problems in unbounded domains using the hyper-singular tractions BIE of quasi-static elasticity.

The elements currently available are strictly discontinuous, mostly with constant interpolation of the field over the elements.

The following kernels are available (version 1.0):
- Full-space isotropic elasticity
  - 2D plane-strain/plane-stress P0 segment (piece-wise constant)
  - 2D plane-strain/plane-stress P1 segment (piece-wise linear) 
  - 3D rectangular R0 element (piece-wise constant) 
  - 3D triangular T0 element (piece-wise constant) 
  
BigWham has bindings for:
- Python (via PyBind 11)
- Julia

Tutorials (using the Python API) with a number of illustrating examples can be found at https://github.com/GeoEnergyLab-EPFL/Bigwham-Tutorials 

How-tos compile BigWham for different architecture can be found in the wiki pages of this GitHub repo.

### Contributors

- Brice Lecampion (2016-): general architecture, elastic kernels implementation, API, H-matrix algorithms implementation, tests
- Francois Fayard (2016-18): H-matrix algorithms implementation, IL library
- Ankit Gupta (2023-): general architecture, API (Julia, python)
- Carlo Peruzzo (2018-2023): Rectangular element kernels, API (Python), tests
- Alexis Sáez (2019-2023): Triangular 0 element kernels, tests  
- Federico Ciardo (2026-2020): 2D tests
- Stéphanie Chaillat-Loiseuille (2016-18): H-matrix algorithms 
- Dmitry Nikolskiy (2016-18): Triangular element kernels, tests
- Lisa Gordeliy (2018): 2D & 3D tests
