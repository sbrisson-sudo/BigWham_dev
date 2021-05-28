#  BigWham : Boundary InteGral equations With HierArchical Matrix


A C++ library for boundary integral equations using hierarchical matrix approximation


V1.alpha  

Current capabilities: 

- Fracture problems in elasto-statics (2d,3d) using hyper-singular kernels (displacement discontinuities elements)
- Collocation BEM

Interfaces
<python>
 -  to compile the binding use python >= 3.7                                      
 -  the interpreter for the .py script should be the same as compile time
 -  add to PYTHONPATH  "/path/to/bigwham/" 

Contributors:
- Brice Lecampion (2D and Simplified 3D kernels, 2D Mesh,  
wolfram API, etc.)  2015-
- Dmitry Nikolskiy (3D Quadratic triangular DD element kernels) 2017-2019
- Francois Fayard  (inside loop lib, HMat C++ code) 2018-
- Stéphanie Chaillat-Loiseuille (HMat algorithms) 2018-
- Carlo Peruzzo (3D mesh, 3D kernels API) 2019-
- Alexis Saez (3D mesh) 2020-
- Lisa Gordeliy (static crack benchmarks)  2019
- Federico Ciardo (2D code parts, static crack benchmarks) 2016-2019


Notes
 
 
