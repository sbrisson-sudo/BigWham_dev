# User API

## Element Type
    
```cpp
#include "core/elements/Triangle.h"
using Tri0 = bigwham::Triangle<0>;
```

## Mesh
    
- Mesh: getCollocationPoints, local to global and global to local transformation of traction and DD vectors

```cpp
#include "core/BEMesh.h"
using Mesh = bigwham::BEMesh<Tri0>;

// coord: coordinates
// conn: connectivity matrix
Mesh my_mesh(coord, conn);

```
    
## Kernel
    
    
- Each Kernel need to define the following virtual `influence` function of `BIE_Kernel` class
    
    ```cpp
    // Source Element : source_elt
    // Reciever Element : source_elt
    // Source Collocation ID: i_s
    // Reciever Collocation ID: i_r
    virtual std::vector<double> influence(Es source_elt, il::int_t i_s,
                                            Er receiver_elt, il::int_t i_r) const;
    ```
    
    ```cpp
    #include "elasticity/3d/BIE_elastostatic_triangle_0_impls.h" // contains the defination of influence function
    #include "core/ElasticProperties.h"
    
    using MatProp = bigwham::ElasticProperties;
    
    // bigwham::ElasticKernelType : H or U or T or V
    using KernelH = bigwham::BIE_elastostatic<Tri0, Tri0, bigwham::ElasticKernelType::H>;
    
    MatProp elas(E, nu);
    // coord.size(1) : Dim of the kernel
    KernelH ker(elas, coord.size(1));
    ```
        
## H-matrix construction
    
 ```cpp
 
 #include "hmat/hmatrix/Hmat.h"
 #include "core/SquareMatrixGenerator.h"
 #include "core/hierarchical_representation.h"
 
 // Construct cluster and block cluster trees
 // max_leaf_size: size of cluster at the leaf level
 // admissibility criteria of block cluster tree
 bigwham::HRepresentation hr =
       bigwham::h_representation_square_matrix(my_mesh, max_leaf_size, eta);
 
 // Make MatrixGenerator: Link Kernel defination with H-matrix
 // Interface to tell H-mat which influence function to call 
 // Right now SquareMatrix Generator with same source and reciever element type is develped
 
 // ker: kernel object
 // my_mesh: mesh object
 
 using MatixGenerator = bigwham::SquareMatrixGenerator<double, Tri0, KernelH>;
 MatixGenerator M(my_mesh, ker, hr.permutation_0_);
 
 // Construct H-mat
 // eps_aca: Indirect way to give maximum rank of compression to ACA (approx SVD)
 bigwham::Hmat<double> h_(M, hr, eps_aca);
 
 // Mat vector operation using H-mat
 // dd: displacement discontunity vector at collocation points
 // trac: traction vector at collocation poitns
 // Right now dd and trac should be in local element coordinate system
 trac = h_.matvec(dd);
 ```
    