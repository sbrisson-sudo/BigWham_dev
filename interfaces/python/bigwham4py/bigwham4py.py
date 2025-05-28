"""
 This file is part of BigWham.

 Created by Carlo Peruzzo on 12.05.21.
 Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
 Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
 file for more details.

 last modifications :: July 2024 by A. Gupta

#######################################################################################
#       IMPORTANT:                                                                    #
#       to compile the binding use python >= 3.7                                      #
#       the interpreter for the .py script should be the same as compile time         #
#######################################################################################
"""

# external
import numpy as np
import time
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import csc_matrix

from scipy.sparse.linalg import spilu
from scipy.sparse import diags

from .py_bigwham import BigWhamIOSelf, BigWhamIORect, PyGetFullBlocks

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

# list kernels
kernels_id = [
    "2DS0-H",
    "2DS1-H",
    "S3DS0-H",
    "Axi3DS0-H",
    "3DT0-H",
    "3DT6-H",
    "3DR0-H",
    "3DR0-H-mode1"
]

##############################
#  Hmatrix class in python   #
##############################
class BEMatrix(LinearOperator):
    def __init__(
        self,
        kernel: str,
        coor: np.ndarray,
        conn: np.ndarray,
        properties: np.ndarray,
        max_leaf_size: int = 32,
        eta: float = 3.0,
        eps_aca: float = 1.0e-3, 
        n_openMP_threads: int =8,
        n_GPUs:int =1,
        directly_build:bool = True,
        verbose:bool = True,
        homogeneous_size_pattern:bool = False,
        fixed_rank = -1,
        useCuda = False
    ):
        """ "
        Name:              Type:                Description:
        _____________      _________________    _______________________________________________
        kernel,            (string)             name of the kernel
        coor,              (numpy 2D arrray)    coordinates of the vertices of the elements (N of vert x 2)
        conn,              (numpy 2D arrray)    connectivity of the elements (N of elemts x N of vert per elem)
        properties,        (numpy 1D arrray)    list of kernel properties (Should be changed to a dict to better implement checks?)
        max_leaf_size,     (integer)            size of the largest sub block (usually 16 - 124)
        eta,               (float)            approximation factor (usually ~3, 0 is non compressed HMAT)
        eps_aca            (float)              approximation factor (usually 0.001 - 0.0001)
        n_openMP_threads    (integer)           number of OMP threads to be used by BigWham
        """
                
        # Ensure kernel exists
        if not(kernel in kernels_id):
            print(f"[ERROR] Invalid kernel : {kernel}, available kernels are : [{', '.join(kernels_id)}]")
            return
        
        self.useCuda = useCuda

        self.kernel_ : str = kernel
        self.properties_ : np.ndarray = properties
        self.max_leaf_size_ : int = int(max_leaf_size)
        self.eta_ : float = float(eta)
        self.eps_aca_ : float = float(eps_aca)
        self.n_openMP_threads_ : int = n_openMP_threads
        self.H_ : BigWhamIOSelf = BigWhamIOSelf(
            coor.flatten(),
            conn.flatten(),
            kernel,
            properties.flatten(),
            n_openMP_threads,
            n_GPUs,
            verbose,
            homogeneous_size_pattern,
            useCuda,
            fixed_rank
        )
        self.built_ = False
        if directly_build :
            self.H_.build_hierarchical_matrix(
                max_leaf_size,
                eta,
                eps_aca,
            )
            self.built_=True
            self.matvec_size_ = self.H_.matrix_size(0)
            # it is mandatory to define shape and dtype of the dot product
            self.shape_ = (self.H_.matrix_size(0), self.H_.matrix_size(1))

        else:
            self.H_.build_pattern(
                max_leaf_size,
                eta,
            )
            # we set dummmy values.
            self.shape_ =(0,0)
        self.dtype_ = float
        super().__init__(self.dtype_, self.shape_)

    def _build(self):
        if not(self.built_):
            self.H_.build_hierarchical_matrix(
                self.max_leaf_size_,
                self.eta_,
                self.eps_aca_,
            )
            self.matvec_size_ = self.H_.matrix_size(0)
            # it is mandatory to define shape and dtype of the dot product
            self.shape_ = (self.H_.matrix_size(0), self.H_.matrix_size(1))
            super().__init__(self.dtype_, self.shape_)
#            self.shape = self.shape_
            self.built_=True
        else:
            pass

    def _matvec(self, v: np.ndarray) -> np.ndarray:
        """
        This function implements the dot product.
        :param v: vector expected to be of size self.HMAT_size_
        :return: HMAT.v
        """
        # shall we put a call to self._build() ?
        return self.H_.matvec(v)
    
    def _matvec_cupy(self, x, y):
        """
        Dot product on cupy arrays, result y is already allocated
        """
        import cupy 
        assert type(x) == cupy.ndarray
        assert type(y) == cupy.ndarray
        assert x.dtype == np.float64
        assert y.dtype == np.float64
        
        return self.H_.matvec_raw_ptr(x.data.ptr, y.data.ptr)
    
    def _matvec_jax(self, x, y):
        """
        Dot product on JAX GPU arrays. Assumes x and y are jax.Array on GPU and of dtype float64.
        """
        import jax

        # Ensure correct types and dtype
        assert x.device.platform == 'gpu'
        assert y.device.platform == 'gpu'
        assert x.dtype == np.float64
        assert y.dtype == np.float64

        # Extract raw pointers from __cuda_array_interface__
        x_ptr = x.__cuda_array_interface__['data'][0]
        y_ptr = y.__cuda_array_interface__['data'][0]

        # Call into C++ backend
        return self.H_.matvec_raw_ptr(x_ptr, y_ptr)

    def write_hmatrix(self, filename: str) -> int:
        """
        write Hmatrix in hdf5
        """
        self.H_.write_hmatrix(filename)

    @property
    def _init_shape(self):
        return self.shape_

    def _init_dtype(self):
        return self.dtype_

    # some useful methods
    def getCompression(self) -> float:
        self._build()
        return self.H_.get_compression_ratio()
    
    def getStorageRequirement(self):
        return self.H_.get_storage_requirement()
    
    def getGPUStorageRequirement(self):
        return self.H_.get_gpu_storage_requirement()

    def getPermutation(self) -> np.ndarray:
        return np.asarray(self.H_.get_permutation())

    def get_omp_threads(self) -> int :
        return self.H_.get_omp_threads()

    def getMeshCollocationPoints(self) -> np.ndarray:
        """
        Get collocation points from mesh (no permutations ....)
        return: (no_collocation_pts, dim) array from mesh
        """
        dim = self.H_.get_spatial_dimension()
        return np.asarray(self.H_.get_collocation_points()).reshape(
            (self.matvec_size_ // dim, dim)
        )

    def convert_to_global(self, x_local: np.ndarray) -> np.ndarray:
        """
        Convert local vector to global vector
        :param x_local: local vector
        :return: global vector
        """
        return self.H_.convert_to_global(x_local)

    def convert_to_local(self, x_global: np.ndarray) -> np.ndarray:
        """
        Convert global vector to local vector
        :param x_global: global vector
        :return: local vector
        """
        return self.H_.convert_to_local(x_global)

    def getCollocationPoints(self) -> np.ndarray:
        """
        Return the coordinates all the collocation points in the mesh
        :return: array of  coordinates of the collocation points
        """
        n = self.H_.get_spatial_dimension()
        aux = np.asarray(self.H_.get_collocation_points())
        colPts = np.reshape(aux, (int(aux.size / n), n))
        return colPts

    def getSpatialDimension(self) -> int:
        return self.H_.get_spatial_dimension()

    def _getFullBlocks(self, original_order=True, keep_ratio: float = 1.0) -> csc_matrix:
        if not (0 < keep_ratio <= 1):
            raise ValueError("ratio must be a float between 0 and 1.")

        fb = PyGetFullBlocks() 
        fb.set(self.H_, original_order)
        val = np.asarray(fb.get_val_list(), dtype=float)
        col = np.asarray(fb.get_col())
        row = np.asarray(fb.get_row())

        if keep_ratio < 1.0:
            # Number of values to keep
            k = int(np.ceil(keep_ratio * val.size))
            if k == 0:
                return csc_matrix(self.shape_)

            # Find threshold
            abs_val = np.abs(val)
            threshold = np.partition(abs_val, -k)[-k]
            mask = abs_val >= threshold
            val = val[mask]
            row = row[mask]
            col = col[mask]

        return csc_matrix((val, (row, col)), shape=self.shape_)

    def _getPattern(self) -> np.ndarray:
        aux = np.asarray(self.H_.get_hpattern())
        #  as flattened list via a pointer
        #  the numberofblocks is also returned (by reference)
        #
        #  the pattern matrix is formatted as
        # row = 1 block : i_begin,j_begin, i_end,j_end,FLAG,entry_size
        # with FLAG=0 for full rank and FLAG=1 for low rank
        #
        # we output a flatten row-major order std::vector
        nr = 6
        return np.reshape(aux, (int(aux.size / nr), nr))
    
    def get_max_error_ACA(self):
        return self.H_.get_max_error_ACA()

    def plotPattern(self, plot_index=False):
        """
        Plot the hierarchical pattern of the matrix
        :return:
        """
        data_pattern = self._getPattern()
        patches = []
        p_colors = []
        max_y = data_pattern[:, 3].max()
        
        fr_counter = 0
        lr_counter = 0
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for i in range(len(data_pattern)):
            
            height = np.abs(data_pattern[i, 0] - data_pattern[i, 2])
            width = np.abs(data_pattern[i, 1] - data_pattern[i, 3])
            y1 = max_y - data_pattern[i, 0] - height
            x1 = data_pattern[i, 1]
            rectangle = Rectangle((x1, y1), width, height)
            patches.append(rectangle)
            p_colors.append(data_pattern[i, 4])
            
            if plot_index:
                idx = fr_counter if data_pattern[i, 4] == 0 else lr_counter
                ax.text(
                    x1 + width / 2, y1 + height / 2, str(idx),
                    ha="center", va="center", fontsize=8, fontweight="bold", color="black"
                )
                
            if data_pattern[i, 4] == 0 : fr_counter+=1
            if data_pattern[i, 4] == 1 : lr_counter+=1
            
        
        p = PatchCollection(
            patches, cmap=matplotlib.cm.PiYG, edgecolors="black", alpha=0.4
        )
        p.set_array(np.array(p_colors))
        ax.add_collection(p)
        ax.set_ylim([data_pattern[:, 0].min(), data_pattern[:, 3].max()])
        ax.set_xlim([data_pattern[:, 1].min(), data_pattern[:, 2].max()])
        ax.set_aspect("equal")
        # # fig.colorbar(p)
        # fig.show()
        # plt.show(block=True)
        return fig

    # a method constructing an ILU Preconditionner of the H matrix
    def H_ILU_prec(self, fill_factor=1, drop_tol=1e-3, keep_ratio=1.0) -> LinearOperator:
        """
        Return an ILU operator (using scipy spilu) built from the full-rank blocks of the matrix.
        :param fill_factor: integer (default 1) for the maximum number of fill-in (see scipy spilu)
        :param drop_tol: float (default 1e-3) for the tolerance to drop the entries (see scipy spilu)
        :return: a linear operator with the corresponding ILU
        """
        # if self.useCuda:
        #     print("[ERROR] The ILU on the hierachical matrix can't be called when using CUDA as the data is not on host memory anymore, falling back to jacobi preconditionner.")
        #     return self.H_jacobi_prec()
        
        self._build()
        fb = self._getFullBlocks(keep_ratio=keep_ratio)
                
        fbILU = spilu(fb, fill_factor=fill_factor, drop_tol=drop_tol)
        return LinearOperator(self.shape_, fbILU.solve)

    def H_diag(self) -> np.ndarray:
        """
        :return: the diagonal of the matrix as n array
        """
        if self.useCuda:
            return self.get_diagonal()
        else :
            self._build()
            fb = self._getFullBlocks()
            return fb.diagonal()

    def H_jacobi_prec(self):
        """
        Create a jacobi preconditioner of the matrix (just the inverse of the diagonal entries)
        :return: a scipy sparse matrix of diagonal format
        """
        diag = self.H_diag()  # return a nd.array
        overdiag = 1.0 / diag
        return diags(overdiag, dtype=np.float64)

    def compute_displacements(self,list_coor,local_solu):
        """
        Compute the induced displacements - for a given value local_solu of the solution over the source mesh
        :param list_coor: array of dim 2. containing the coordinates where displacements will be evaluated
        :param local_solu: vector of solution (in the local system)
        :return: np.array of the values of the displacement components at the given points
        """
        self._build()
        n = self.H_.get_spatial_dimension()
        assert n == list_coor.shape[1], "Coordinates dimension of the given points is not matching the problem spatial dimension !"
        u = self.H_.compute_displacements(list_coor.flatten(),local_solu.flatten())
        return np.reshape(u,(-1,n))

    def compute_stresses(self,list_coor,local_solu):
        """
         Compute the induced stresses - for a given value local_solu of the solution over the source mesh
        :param list_coor:  array of dim 2. containing the coordinates where displacements will be evaluated
        :param local_solu: vector of solution (in the local system)
        :return: np.array of the values of the stress components at the given points
        """
        self._build()
        n = self.H_.get_spatial_dimension()
        assert (n==2 or n ==3 ), " Problem spatial dimension must be 2 or 3"
        assert n == list_coor.shape[1], "Coordinates dimension of the given points is not matching the problem spatial dimension !"
        nstress = 3
        if n == 3:
            nstress = 6
        sig = self.H_.compute_stresses(list_coor.flatten(),local_solu.flatten())
        return np.reshape(sig,(-1,nstress))

    def get_element_normals(self):
        """
        Get the normal vectors of the elements
        :return: 1D np.array of the normal vectors of the elements
        3d case:
        [n0x,n0y,n0z,n1x,n1y,n1z,...]
        2d case:
        [n0x,n0y,n1x,n1y,...]
        """
        return np.asarray(self.H_.get_element_normals())

    def get_rotation_matrix(self):
        """
        Get the rotation matrix of the elements
        :return: 1D np.array of the rotation matrix of the elements
        3d case:
        [r00,r01,r02,r10,r11,r12,r20,r21,r22,...]
        2d case:
        [r00,r01,r10,r11,...]
        """
        return np.asarray(self.H_.get_rotation_matrix())
    
    def get_diagonal(self):
        """
        Get the diagonal of the matrix, in the original ordering of dof
        :return: 1D np.array of the diagonal
        """
        return self.H_.get_diagonal()
    
    def isCudaAvailable(self):
        """
        Return if BigWham has been compiled with CUDA support
        :return: boolean value
        """
        return self.H_.get_cuda_available()
    
    def getMatvecTime(self, N_matvec=100):
        """
        Time the matvec operation
        :return: matvec time im seconds
        """
        x = np.ones(self.shape[0])
        # One first matvec performed not taken into account
        self._matvec(x)
        # Then we time it 
        start_time = time.time()
        for _ in range(N_matvec):
            self._matvec(x)
        return (time.time() - start_time)/N_matvec
    
    def __getitem__(self, indices):
        """
        To get a selection of the hmat 
        """
        if self.useCuda:
            raise Exception("Bigwham matrix selection not implemented for CUDA support.")
        
        if not isinstance(indices, tuple) or len(indices) != 2:
            raise IndexError("Subsetting requires indices of the form (np.ix_(row_indices, col_indices))")
        
        row_indices, col_indices = indices
        
        # Copy attributes (shallow copy, so mutable attributes are shared)
        new_bematrix = object.__new__(BEMatrix)
        new_bematrix.__dict__ = self.__dict__.copy()
        
        # Update its hmatrix
        new_bematrix.H_ = self.H_.hmatSelection(row_indices, col_indices)
        
        # And its shape
        new_bematrix.shape = (new_bematrix.H_.matrix_size(0),new_bematrix.H_.matrix_size(1))
        new_bematrix.shape_ = (new_bematrix.H_.matrix_size(0),new_bematrix.H_.matrix_size(1))
        
        return new_bematrix
        
########################################################################################################
#  BEMatrix Rectangular class in python   #
########################################################################################################

class BEMatrixRectangular(LinearOperator):
    def __init__(
        self,
        kernel : str,
        coor_src : np.ndarray,
        conn_src : np.ndarray,
        coor_rec : np.ndarray,
        conn_rec : np.ndarray,
        properties : np.ndarray,
        max_leaf_size : int =100,
        eta : float=3.0,
        eps_aca : float=1.0e-3, 
        n_openMP_threads:int =8,
        n_GPUs:int =1,
        directly_build:bool = True,
        verbose:bool = True,
        homogeneous_size_pattern:bool = False,
        useCuda = False,
        fixed_rank = -1
    ):

        self.kernel_ = kernel
        self.properties_ = properties
        self.max_leaf_size_ = max_leaf_size
        self.eta_ = eta
        self.eps_aca_ = eps_aca
        self.n_openMP_threads_=n_openMP_threads
        self.built_=False

        self.H_ : BigWhamIORect= BigWhamIORect(
            coor_src.flatten(),
            conn_src.flatten(),
            coor_rec.flatten(),
            conn_rec.flatten(),
            kernel,
            properties.flatten(),
            n_openMP_threads,
            n_GPUs,
            verbose,
            homogeneous_size_pattern,
            useCuda,
            fixed_rank
        )
        self.dtype_ = float

        if directly_build:
            self.H_.build_hierarchical_matrix(
            max_leaf_size,
            eta,
            eps_aca,
            )
            # it is mandatory to define shape and dtype of the dot product
            self.matvec_size_ = self.H_.matrix_size(0)
            self.shape_ = (self.H_.matrix_size(0), self.H_.matrix_size(1))
            self.built_=True
        else:
            self.H_.build_pattern(max_leaf_size,eta)

        self.dtype_ = float
        super().__init__(self.dtype_, self.shape_)


    def _build(self):
        if not(self.built_):
            self.H_.build_hierarchical_matrix(
                self.max_leaf_size_,
                self.eta_,
                self.eps_aca_,
            )
            self.built_=True
            self.matvec_size_ = self.H_.matrix_size(0)
            self.shape_ = (self.H_.matrix_size(0), self.H_.matrix_size(1))
            super().__init__(self.dtype_, self.shape_)
        else:
            pass

    def _matvec(self, v: np.ndarray) -> np.ndarray:
        """
        This function implements the dot product.
        :param v: vector expected to be of size self.HMAT_size_
        :return: HMAT.v
        """
        return self.H_.matvec(v)

    @property
    def _init_shape(self):
        return self.shape_

    def _init_dtype(self):
        return self.dtype_

    # some useful methods
    def getCompression(self) -> float:
        self._build()
        return self.H_.get_compression_ratio()

    def getPermutation(self) -> np.ndarray:
        return np.asarray(self.H_.get_permutation())
    
    def getPermutationReceivers(self) -> np.ndarray:
        return np.asarray(self.H_.get_permutation_receivers())

    def get_omp_threads(self) -> int:
        return self.H_.get_omp_threads()

    def getMeshCollocationPoints(self) -> np.ndarray:
        """
        Get collocation points from mesh (no permutations ....)
        return: (no_collo_pts, dim) array from mesh
        """
        dim = self.H_.get_spatial_dimension()
        return np.asarray(self.H_.get_collocation_points()).reshape(
            (self.matvec_size_ // dim, dim)
        )

    def convert_to_global(self, x_local: np.ndarray) -> np.ndarray:
        """
        Convert local vector to global vector
        :param x_local: local vector
        :return: global vector
        """
        return self.H_.convert_to_global(x_local)

    def convert_to_local(self, x_global: np.ndarray) -> np.ndarray:
        """
        Convert global vector to local vector
        :param x_global: global vector
        :return: local vector
        """
        return self.H_.convert_to_local(x_global)

    def getCollocationPoints(self) -> np.ndarray:
        n = self.H_.get_spatial_dimension()
        aux = np.asarray(self.H_.get_collocation_points())
        auxpermut = np.reshape(aux, (int(aux.size / n), n))
        permut = self.getPermutation()
        colPts = 0.0 * auxpermut
        colPts[permut] = auxpermut
        return colPts

    def getSpatialDimension(self) -> int:
        return self.H_.get_spatial_dimension()

    def _getPattern(self) -> np.ndarray:
        aux = np.asarray(self.H_.get_hpattern())
        nr = 6
        return np.reshape(aux, (int(aux.size / nr), nr))

    def plotPattern(self):
        data_pattern = self._getPattern()
        patches = []
        p_colors = []
        max_y = data_pattern[:, 3].max()
        for i in range(len(data_pattern)):
            height = np.abs(data_pattern[i, 0] - data_pattern[i, 2])
            width = np.abs(data_pattern[i, 1] - data_pattern[i, 3])
            y1 = max_y - data_pattern[i, 0] - height
            x1 = data_pattern[i, 1]
            rectangle = Rectangle((x1, y1), width, height)
            patches.append(rectangle)
            p_colors.append(data_pattern[i, 4])
        fig = plt.figure()
        ax = fig.add_subplot(111)

        p = PatchCollection(
            patches, cmap=matplotlib.cm.PiYG, edgecolors="black", alpha=0.4
        )
        p.set_array(np.array(p_colors))
        ax.add_collection(p)
        ax.set_ylim([data_pattern[:, 0].min(), data_pattern[:, 3].max()])
        ax.set_xlim([data_pattern[:, 1].min(), data_pattern[:, 2].max()])
        ax.set_aspect("equal")
        # fig.colorbar(p)
        # fig.show()
        # plt.show(block=True)
        return fig

    def _getFullBlocks(self, original_order=True) -> csc_matrix:
        fb = PyGetFullBlocks()  # not fan of this way of creating empty object
        # and setting them after - a constructor should do something!
        fb.setRect(self.H_, original_order)
        val = np.asarray(fb.get_val_list(), dtype=float)
        col = np.asarray(fb.get_col())
        row = np.asarray(fb.get_row())
        
        # print(f"{col.min()=} {col.max()=}", flush=True)
        # print(f"{row.min()=} {row.max()=}", flush=True)
        
        return csc_matrix((val, (row, col)), shape=self.shape_)
    
    def isCudaAvailable(self):
        """
        Return if BigWham has been compiled with CUDA support
        :return: boolean value
        """
        return self.H_.get_cuda_available()
    
    # def get_diagonal(self):
    #     """
    #     Get the diagonal of the matrix, in the original ordering of dof
    #     :return: 1D np.array of the diagonal
    #     """
    #     return self.H_.get_diagonal()
    
    def getMatvecTime(self, N_matvec=100):
        """
        Time the matvec operation
        :return: matvec time im seconds
        """
        x = np.ones(self.shape[0])
        # One first matvec performed not taken into account
        self._matvec(x)
        # Then we time it 
        start_time = time.time()
        for _ in range(N_matvec):
            self._matvec(x)
        return (time.time() - start_time)/N_matvec

def main():
    print("bigwham4py successfully imported")
    
if __name__ == "__main__":
    main()