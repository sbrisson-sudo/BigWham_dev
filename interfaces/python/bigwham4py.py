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
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import csc_matrix

from scipy.sparse.linalg import spilu
from scipy.sparse import diags

from py_bigwham import BigWhamIOSelf, BigWhamIORect, PyGetFullBlocks

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


##############################
#  Hmatrix class in python   #
##############################
class BEMatrix(LinearOperator): #, metaclass=multimeta
    def __init__(
        self,
        kernel: str,
        coor: np.ndarray,
        conn: np.ndarray,
        properties: np.ndarray,
        max_leaf_size: int = 100,
        eta: float = 3.0,
        eps_aca: float = 1.0e-3,
    ):
        """ "
        Name:              Type:                Description:
        _____________      _________________    _______________________________________________
        kernel,            (string)             name of the kernel
        coor,              (numpy 2D arrray)    coordinates of the vertices of the elements (N of vert x 2)
        conn,              (numpy 2D arrray)    connectivity of the elements (N of elemts x N of vert per elem)
        properties,        (numpy 1D arrray)    list of kernel properties (Should be changed to a dict to better implement checks?)
        max_leaf_size,     (integer)            size of the largest sub block (usually 100 - 1000)
        eta,               (float)            approximation factor (usually ~3, 0 is non compressed HMAT)
        eps_aca            (float)              approximation factor (usually 0.001 - 0.0001)
        """

        self.kernel_ : str = kernel
        self.properties_ : np.ndarray = properties
        self.max_leaf_size_ : int = int(max_leaf_size)
        self.eta_ : float = float(eta)
        self.eps_aca_ : float = float(eps_aca)

        self.H_ : BigWhamIOSelf = BigWhamIOSelf(
            coor.flatten(),
            conn.flatten(),
            kernel,
            properties.flatten(),
        )
        self.H_.build_hierarchical_matrix(
            max_leaf_size,
            eta,
            eps_aca,
        )
        self.matvec_size_ = self.H_.matrix_size(0)
        self.dtype_ = float

        # it is mandatory to define shape and dtype of the dot product
        self.shape_ = (self.H_.matrix_size(0), self.H_.matrix_size(1))
        super().__init__(self.dtype_, self.shape_)



    def _matvec(self, v: np.ndarray) -> np.ndarray:
        """
        This function implements the dot product.
        :param v: vector expected to be of size self.HMAT_size_
        :return: HMAT.v
        """
        return self.H_.matvec(v)

    def write_hmatrix(self, filename: str) -> int:
        """
        write Hmatirx in hdf5
        """
        self.H_.write_hmatrix(filename)

    @property
    def _init_shape(self):
        return self.shape_

    def _init_dtype(self):
        return self.dtype_

    # some useful methods
    def getCompression(self) -> float:
        return self.H_.get_compression_ratio()

    def getPermutation(self) -> np.ndarray:
        return np.asarray(self.H_.get_permutation())

    def get_omp_threads(self) -> int :
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

    def _getFullBlocks(self) -> csc_matrix:
        fb = PyGetFullBlocks()  # not fan of this way of creating empty object
        # and setting them after - a constructor should do something!
        fb.set(self.H_)
        val = np.asarray(fb.get_val_list(), dtype=float)
        col = np.asarray(fb.get_col())
        row = np.asarray(fb.get_row())
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
        fig.show()
        plt.show(block=True)
        return fig

    # a method constructing an ILU Preconditionner of the H matrix
    def H_ILU_prec(self, fill_factor=5, drop_tol=1e-5) -> LinearOperator:
        fb = self._getFullBlocks()
        fbILU = spilu(fb, fill_factor=fill_factor, drop_tol=drop_tol)
        return LinearOperator(self.shape_, fbILU.solve)


    def H_diag(self) -> np.ndarray:
        fb = self._getFullBlocks()
        return fb.diagonal()

    def H_jacobi_prec(self):
        diag = self.H_diag()  # return a nd.array
        overdiag = 1.0 / diag
        return diags(overdiag, dtype=np.float_)

    def compute_displacements(self,list_coor,local_solu):
        """
        Compute the induced displacements - for a given value local_solu of the solution over the source mesh
        :param list_coor: array of dim 2. containing the coordinates where displacements will be evaluated
        :param local_solu: vector of solution (in the local system)
        :return: np.array of the values of the displacement components at the given points
        """
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
    ):

        self.kernel_ = kernel
        self.properties_ = properties
        self.max_leaf_size_ = max_leaf_size
        self.eta_ = eta
        self.eps_aca_ = eps_aca

        self.H_ : BigWhamIORect= BigWhamIORect(
            coor_src.flatten(),
            conn_src.flatten(),
            coor_rec.flatten(),
            conn_rec.flatten(),
            kernel,
            properties.flatten(),
        )
        self.H_.build_hierarchical_matrix(
            max_leaf_size,
            eta,
            eps_aca,
        )

        self.matvec_size_ = self.H_.matrix_size(0)
        self.dtype_ = float

        # it is mandatory to define shape and dtype of the dot product
        self.shape_ = (self.H_.matrix_size(0), self.H_.matrix_size(1))
        super().__init__(self.dtype_, self.shape_)

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
        return self.H_.get_compression_ratio()

    def getPermutation(self) -> np.ndarray:
        return np.asarray(self.H_.get_permutation())

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
