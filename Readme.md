#  BigWham : Boundary InteGral equations With HierArchical Matrix


## A C++ library for vectorial boundary integral equations with hierarchical matrix 


## V1.alpha  

## Current capabilities: 

- Fracture problems in elasto-statics (2D,3D) using hyper-singular kernels (displacement discontinuities elements)
- Collocation BEM
- Hierarchical matrix
- Currently available kernels:
    + 2D Full-space P1 discontinuous element
    + Full-space Simplified3D / 2D  P0 element
    + Full-space 3D T0 element 
    + Full-space 3D T6 element
    
## Interfaces
### Python
 -  to compile the binding use python >= 3.7                                      
 -  the interpreter for the .py script should be the same as compile time
 -  add to PYTHONPATH  "YOUR_PATH_TO_BIGWHAM base folder" 
### Wolfram language
- see the BigWhamLink repo (to be included in this report in the future)

#Contributors:
- Brice Lecampion (2D and Simplified 3D kernels, 2D Mesh,  
wolfram API, etc.)  2015-
- Francois Fayard  (inside loop lib, HMat C++ initial code) 2018
- Dmitry Nikolskiy (3D Quadratic triangular DD element kernels) 2017-2019
- Carlo Peruzzo (3D mesh, 3D kernels API, Python Bindings) 2019-
- Alexis Saez (3D mesh) 2020-
- Federico Ciardo (2D code parts, 2D static crack benchmarks) 2016-2019
- Lisa Gordeliy (static crack benchmarks)  2019
- Stéphanie Chaillat-Loiseuille (HMat algorithms) 2018

