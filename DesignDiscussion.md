Bigwham - Code design discussion 
=====
date: Jan. 15 2023


Goals:
----
- A solver for vectorial Boundary Integral Equations leveraging hierarchical matrices
- Primarily developed in the scope of collocation Boundary Element Methods
- Primarily geared toward Quasi-Static Elasticity
- Separation between element kernel influence functions (the core engine of a BE method) and hierarchical matrix
- Provide a interface for mat_vec dot product, + approximated preconditionner (full-blocks or diag. sparse matrix) 
that can be used externally by Iterative Methods for the solution of linear system
- fast ? could solve "large"  problems ( >10^6 dofs - currently our limit)
- [Note initially: a H-LU was developped by F. Fayard-> this part is now "shelved"]


History & notes:
----
- first developments of quasi-static kernels (2DP1) by BL (end 2015), 
- first H-mat code in 2016-17 by F. Fayard (with help from S. Chaillat (Coffee code)):
  + developed using recursion for both construction and mat_vect_ + H-LU
  + MatrixGenerator base class as the base functor to derived from for different kernels
- in 2019-20-21, BL developed a Wolfram Library Link plugin to mathematica:
  + devlopment of the Bigwhamio class for wrapping things  
- In 2020-21-22, CP developed a Python interface via PyBind11
- in 2021-22, BL speed up significantly the mat_vect_ by removing recursion 
during Hmat construction and mat_vect:
  + storing the Hmat pattern first 
  + having 2 std:vect for the Full Rank, and Low Rank blocks
  + simple open-MP pragma to speed up the construction + mat_vect
  + note the old "way" has not yet been removed from the code -> should be done to avoid tbb dependency
- Currently, for mixed BVP, the two H-matrices (e.g. A,B) for the corresponding 2 integrals must be created and
the linear system solve via a approach res(X=[u_u,t_u])= A.[u_u,u_k]+B.[t_k,t_u]=0 where the subscripts _u refer 
to the unknowns, and _t to the knwon boundary value. this of-course double the computational cost as at each
iterations of the iterative linear solve, 2 Hmat_vect must be performed (+mem cost). 
Doing otherwise would require specific matrix generator for mixed BVP (could be developed but not on 
the top of the pile)
- Currently, the post-processing (computation of unknowns at given internal point) do not use hiearchical matrix,
this could be done (of interest when computing the effect of multiple sol. vectors for inverse problems etc.)
by developing rectangular Hmat using mesh_receivers, mesh_source. However, note that
the internal block must remain rectangular for the vectorial-ACA to work, in other words,
we can not have an element influence function returning a rectangular matrix for the influence of
an element node of the BEmesh on a receiver points. E.g. in 3D for displacement, this is ok. However, for
stress (or strain) the influence return 6 components while there is 3 solution vector (either displacement or tractions on the source element).
As a result, a pragmatic way would be to: have 2 rectangular (because the # source != # receiver) Hmat with
symmetric 3*3 blocks   : one for the first 3 stress components, one for the last three.
[note a  problem exist in 2D, for which we have 2 dofs for 3 stress components -> could be solved by using complex variables of plane elasticity]

Current improvements required:
----
- simplify the API, avoid code duplication, have clearer data containers
- avoid the requirement of writing a MatrixGenerator for each kernels
  + this can be achieved by recognizing that intrinsically, the MatrixGenerator is
mostly governed by the vectorial dimension of the BIE (scalar, 2 dofs, 3 dofs)
- have a separation between: 
  + element_type (which include element geometry + interpolation),
  + problem type (elastic_fs_iso, elastic_fs_iso_modeI, elastic_fs_TI etc.)
  + kernel type (single layer, double layer, adjoint double layer, quadruple layer : or in the elastic terminology U, T, H)
- Have a single BoundaryElement Mesh class instead of a 2D and a 3D class\
  +   I already wrote a general class merging the 2  (usage to be rolled out in the code)
- re-assess, possibly write a BEModel class (to be more general and extandable to other elements, configuration, kernels)
- Note -> roll-out the use of either local or global coordinates system choice (I_want_global, I_want_global_codomain) for all kernels


Design Choices ? 
----
This is matter of discussion....

- The building block of BE methods is the influence of one element onto another.
More specifically, in the scope of Hierarchical matrix, the effect of one node of the source element onto the
effect of one node of the receiver element (or receiver points in the case of observation).
  + A boundary element type is defined by the element geometry  and the order of interpolation of the field over 
  the element.
  + In addition, the integral equations is defined by its underlying physics / fundamental solutions - a.k.a kernel - which include 
  material properties (isotropy or not), if the fundamental solution is for the full-space, or Half-space, or... 
  in addition, if the fundamental solution is for displacement, stress, eigenstress...
  + a boundary element can be in 3D but the underlying bie be scalar !    the spatial dimension is not necessarily equal to the dof_dimensions
  
Proposal:
  + have an enum for the set of implemented Boundary element type : {segment_0,segment_1,rectangle_0,triangle_0,triangle_2 ....}
  + have an enum for the set of 'kernel_physics' implemented :{elastic_fs_iso,elastic_fs_iso_mode_1,elastic_fs_simplified3D,elastic_hs_iso, elastic_fs_TI ...}
  + have an enum for the kernels implemented : {U,T,H,S,V}   - note depending on physics (we or we do not have everything implemented)
  + have same signature for the corresponding influence functions (node of source element-> node of receiver element) 
such that we can define their type them via typedef. 
    + Typically, we will have 1 for scalar, 1 for 2D dofs, 1 for 3D dofs. 
    + At the moment, the signatures use a type for the element data if the element geometry is 2D (SegmentData) or 3D (FaceData)
      typedef il::StaticArray2D<double, 2, 2> (*vKernelCallNodal2x2)(
      SegmentData &source_elt, SegmentData &receiver_elt, il::int_t s_col,  il::int_t i_col,
      const ElasticProperties &Elas, double ker_options); 
      + we shall find a way to have only one class for ElementData (templated with spatial dimension?) in order to have less virtual function type to define 
  + use these type definition for the function in the MatrixGenerator 
  + need of a kind of mapping for using the right influence function (to define the Matrix Generator) from the knowledge of the Boundary element type, 
the kernel_physics type, and the kernels type.  Needs to use a std:variant to handle the different influence function signature ?

Commenting Conventions:
------
The code is the documentation ! -> avoid having large comments - the code should be readable by an appropriate choice of names for the variables etc.
no need to say what 2 axes of an array means etc. This is understandable by just reading the code. Of course, that means that people knows how to read a code! 
This is obviously a pre-requise - and duplicating with words what the code does is NOT the way to compensate that.


Coding Conventions
------
class name start with Caps
private attributes of class ends with _
function/method name start lower cases 
variable name all lower caps
...
