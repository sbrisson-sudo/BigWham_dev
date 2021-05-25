"""
 This file is part of BigWham.

 Created by Carlo Peruzzo on 12.05.21.
 Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
 Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
 file for more details.

 last modifications :: May. 12 2021

#######################################################################################
#       IMPORTANT:                                                                    #
#       to compile the binding use python >= 3.7                                      #
#       the interpreter for the .py script should be the same as compile time         #
#######################################################################################
"""

# external
import numpy as np
from scipy.sparse.linalg import LinearOperator


##############################
# Hdot operator for GMRES    #
##############################
class Hdot(LinearOperator):
    def __init__(self):
        from src.binding._test.bigwham4py import Bigwhamio
        self.HMAT_size_ = None
        self.matvec_size_ = None
        self.shape_ = None
        self.dtype_ = float
        self.HMAT = Bigwhamio()

    def set(self, data):
        """

        :param data: a list containing in the order: [kernel, coor, conn, Young, PoissonRatio, max_leaf_size, eta, eps_aca]

                    Name:              Type:                Description:
                    _____________      _________________    _______________________________________________
                    kernel,            (string)             name of the kernel
                    coor,              (numpy 2D arrray)    coordinates of the vertices of the elements (N of vert x 2)
                    conn,              (numpy 2D arrray)    connectivity of the elements (N of elemts x N of vert per elem)
                    Young,             (float)              Young's modulus
                    PoissonRatio,      (float)              Poisson's ratio
                    max_leaf_size,     (integer)            size of the largest sub block (usually 100 - 1000)
                    eta,               (integer)            approximation factor (usually 10 - 100, 0 is non compressed HMAT)
                    eps_aca            (float)              approximation factor (usually 0.001 - 0.0001)

        :return: none
        """

        # unpaking the data
        kernel, coor, conn, Young, PoissonRatio, max_leaf_size, eta, eps_aca = data

        # flatten the arrays and prepare the material properties
        coor = coor.flatten()
        conn = conn.flatten()
        properties = [Young, PoissonRatio]

        # checks
        nodes_per_element_ = 4
        n_of_elts_ = int(len(conn) / nodes_per_element_)
        if len(conn) % nodes_per_element_ != 0 :
            print(" ERROR: \n ")
            print(" wrong connectivity dimension \n ")


        # define the HMAT size
        # define the total number of unknowns to be output by the matvet method
        unknowns_per_element_ = 3
        self.HMAT_size_ = int(n_of_elts_ * unknowns_per_element_)
        self.matvec_size_ = self.HMAT_size_


        # it is mandatory to define shape and dtype of the dot product
        self.shape_ = (self.matvec_size_, self.matvec_size_)
        super().__init__(self.dtype_, self.shape_)


        # set the objects
        print("  ")
        print(" --------------------------------------- ")
        self.HMAT.set(coor.tolist(),
                           conn.tolist(),
                           kernel,
                           properties,
                           max_leaf_size,
                           eta,
                           eps_aca)
        print(" --------------------------------------- ")
        print("  ")
        print("   -> KERNEL compr. ratio = " + str(self.HMATtract.getCompressionRatio()))


    def _matvec(self, v):
        """
        This function implements the dot product.
        :param v: vector expected to be of size self.HMAT_size_
        :return: HMAT.v
        """
        return self.HMAT.hdotProduct(v)

    def _changeShape(self, shape_):
        self.shape_ = (shape_,shape_)
        super().__init__(self.dtype_, self.shape_)

    @property
    def _init_shape(self):
        return self.shape_

    def _init_dtype(self):
        return self.dtype_

#--------------------------------