"""
 This file is part of BigWham.

 Created by Carlo Peruzzo on 20.12.22.
 Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
 Geo-Energy Laboratory, 2016-2021.  All rights reserved. See the LICENSE.TXT
 file for more details.

 last modifications :: May 25, 2021

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
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spilu
from scipy.sparse import diags

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from interfaces.python.lib.bigwhamPybind import *


### building a cartesian mesh object

# -*- coding: utf-8 -*-
"""
This file is part of PyFrac.

Created by Haseeb Zia on Thu Dec 22 11:51:00 2016.
Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2020. All rights reserved.
See the LICENSE.TXT file for more details.
"""

# external imports
import numpy as np
import logging

class CartesianMesh:
    """Class defining a Cartesian Mesh.

    The constructor creates a uniform Cartesian mesh centered at (0,0) and having the dimensions of [-Lx,Lx]*[-Ly,Ly].

    Args:
        nx,ny (int):                -- number of elements in x and y directions respectively.
        Lx,Ly (float):              -- lengths in x and y directions respectively.

    Attributes:
        Lx,Ly (float):                    -- length of the domain in x and y directions respectively. The rectangular domain
                                             have a total length of 2*Lx in the x direction and 2*Ly in the y direction. Both
                                             the positive and negative halves are included.
        nx,ny (int):                      -- number of elements in x and y directions respectively.
        hx,hy (float):                    -- grid spacing in x and y directions respectively.
        VertexCoor  (ndarray):            -- [x,y] Coordinates of the vertices.
        NumberOfElts (int):               -- total number of elements in the mesh.
        Connectivity (ndarray):           -- connectivity array giving four vertices of an element in the following order
                                             [bottom left, bottom right, top right, top left]
        domainLimits (ndarray):           -- the limits of the domain

    """

    def __init__(self, Lx, Ly, nx, ny, symmetric=False):
        """
        Creates a uniform Cartesian mesh centered at zero and having the dimensions of [-Lx, Lx]*[-Ly, Ly].

        Args:
            nx,ny (int):        -- number of elements in x and y directions respectively
            Lx,Ly (float):      -- lengths in x and y directions respectively
            symmetric (bool):   -- if true, additional variables (see list of attributes) will be evaluated for
                                    symmetric fracture solver.

        """
        log = logging.getLogger('PyFrac.mesh')

        # set the limits of the mesh discretisation
        self.set_domainLimits(Lx,Ly)

        # check if the number of cells is odd to see if the origin would be at the mid point of a single cell
        self.set_cell_number(nx, ny)

        # set the sizes of each cell
        self.hx = 2. * self.Lx / (self.nx - 1)
        self.hy = 2. * self.Ly / (self.ny - 1)

        # set the limits of the physical domain
        self.set_physDomainLimits()

        # get the coordinates of the mesh vertexes
        self.VertexCoor = self.get_VertexCoor()

        # set the total number of nodes in the mesh
        self.NumberofNodes = (self.nx+1) * (self.ny+1)

        # set the total number of elements in the mesh
        self.NumberOfElts = self.nx * self.ny


        """
         CONNECTIVITY ARRAYS:
         
         - conn is the connectivity array giving four vertices of an element in the following order:
         ______ ______ _____ 
        |      |      |     |
        |______3______2_____|
        |      |  i   |     |
        |______0______1_____|
        |      |      |     |
        |______|______|_____|
             
        """
        conn = np.empty([self.NumberOfElts, 4], dtype=int)
        k = 0
        for j in range(0, self.ny):
            for i in range(0, self.nx):
                # computing the connectivity elem-elem
                conn[k, 0] = (i + j * (self.nx + 1))
                conn[k, 1] = (i + 1) + j * (self.nx + 1)
                conn[k, 2] = i + 1 + (j + 1) * (self.nx + 1)
                conn[k, 3] = i + (j + 1) * (self.nx + 1)
                k = k + 1

        self.Connectivity = conn

    # -----------------------------------------------------------------------------------------------------------------------

    def set_domainLimits(self, Lx, Ly):

        """
        Notes: (see the picture below)
        1) the domain limits considers dimensions between cell centers
           while the physical domain is actually larger

        2) depending on the input one can set:

            |<-- -Lx -->|<-- +Lx -->|
            |           |           |
            |<--Lx[0]-->|<--Lx[1]-->|
         ___|___ _______|_______ ___|___
        |   |   |       |       |   |   |
        |   x   |   x   |   x   |   x   |
        |_______|_______|_______|_______|
        |       |       |       |       |
        |   x   |   x   |   X   |   x   |
        |_______|_______|_______|_______|
        |       |       |       |       |
        |   x   |   x   |   x   |   x   |
        |_______|_______|_______|_______|

        """
        if not isinstance(Lx, list):
            self.Lx = Lx
            xlims = np.asarray([-Lx, Lx])
        else:
            self.Lx = abs(Lx[0]-Lx[1]) / 2.
            xlims = np.asarray([Lx[0], Lx[1]])

        if not isinstance(Ly, list):
            self.Ly = Ly
            ylims = np.asarray([-Ly, Ly])
        else:
            self.Ly = abs(Ly[0]-Ly[1]) / 2.
            ylims = np.asarray([Ly[0], Ly[1]])

        self.domainLimits = np.hstack((ylims, xlims))

    def set_physDomainLimits(self):
        [yCmin, yCmax, xCmin, xCmax] = self.domainLimits
        hxHalf = self.hx/2.
        hyHalf = self.hy/2.
        self.physDomainLimits = [yCmin-hyHalf, yCmax+hyHalf, xCmin-hxHalf, xCmax+hxHalf]
    # -----------------------------------------------------------------------------------------------------------------------

    def set_cell_number(self, nx, ny):
        log = logging.getLogger('PyFrac.mesh.set_cell_number')
        # Check if the number of cells is odd to see if the origin would be at the mid point of a single cell
        if nx % 2 == 0:
            log.warning("Number of elements in x-direction are even. Using " + repr(nx+1) + " elements to have origin at a "
                                                                                            "cell center...")
            self.nx = nx+1
        else:
            self.nx = nx

        if ny % 2 == 0:
            log.warning("Number of elements in y-direction are even. Using " + repr(ny+1) + " elements to have origin at a "
                                                                                            "cell center...")
            self.ny = ny+1
        else:
            self.ny = ny
    # ----------------------------------------------------------------------------------------------------------------------
    def get_VertexCoor(self):
        x = np.linspace(self.domainLimits[2] - self.hx / 2., self.domainLimits[3] + self.hx / 2., self.nx + 1)
        y = np.linspace(self.domainLimits[0] - self.hy / 2., self.domainLimits[1] + self.hy / 2., self.ny + 1)

        xv, yv = np.meshgrid(x, y)  # coordinates of the vertex of each elements

        a = np.resize(xv, ((self.nx + 1) * (self.ny + 1), 1))
        b = np.resize(yv, ((self.nx + 1) * (self.ny + 1), 1))

        VertexCoor = np.reshape(np.stack((a, b), axis=-1), (len(a), 2))

        return VertexCoor



##############################
#  Hmatrix class in python   #
##############################
class Hmatrix(LinearOperator):
    def __init__(self,kernel,coor,conn,properties,max_leaf_size=100,eta=3,eps_aca=1.e-3):
        """"
                    Name:              Type:                Description:
                    _____________      _________________    _______________________________________________
                    kernel,            (string)             name of the kernel
                    coor,              (numpy 2D arrray)    coordinates of the vertices of the elements (N of vert x 3)
                    conn,              (numpy 2D arrray)    connectivity of the elements (N of elemts x N of vert per elem)
                    properties,        (numpy 2D arrray)    list of kernel properties (Should be changed to a dict to better implement checks?)
                    max_leaf_size,     (integer)            size of the largest sub block (usually 100 - 1000)
                    eta,               (integer)            approximation factor (usually ~3, 0 is non compressed HMAT)
                    eps_aca            (float)              approximation factor (usually 0.001 - 0.0001)
        """
        # NOTE: no specific kernel checks are implemented for now - for grown-ups only

        self.kernel_=kernel
        self.properties_=properties
        self.max_leaf_size_=max_leaf_size
        self.eta_=eta
        self.eps_aca_=eps_aca

        self.H_ = Bigwhamio()
        self.H_.set(coor.flatten(), conn.flatten(), kernel, properties.flatten(), max_leaf_size, eta, eps_aca)
        self.matvec_size_ = self.H_.matrixSize(0)
        self.dtype_ = float

        # it is mandatory to define shape and dtype of the dot product
        self.shape_ = (self.H_.matrixSize(0), self.H_.matrixSize(1))
        super().__init__(self.dtype_, self.shape_)

    def _matvec(self, v):
        """
        This function implements the dot product.
        :param v: vector expected to be of size self.HMAT_size_
        :return: HMAT.v
        """
        return self.H_.hdotProduct(v)

    @property
    def _init_shape(self):
        return self.shape_

    def _init_dtype(self):
        return self.dtype_

    # some useful methods
    def getCompression(self):
        return self.H_.getCompressionRatio()

    def getPermutation(self):
        return np.asarray(self.H_.getPermutation())

    def getCollocationPoints(self):
        n=self.H_.getSpatialDimension()
        aux=np.asarray(self.H_.getCollocationPoints())
        auxpermut=np.reshape(aux,(int(aux.size/n),n))
        permut=self.getPermutation()
        colPts=0.*auxpermut
        colPts[permut]=auxpermut
        return colPts

    def getSpatialDimension(self):
        return self.H_.getSpatialDimension()

    def _getFullBlocks(self):
        fb=pyGetFullBlocks()  # not fan of this way of creating empty object and setting them after - a constructor should do something!
        fb.set(self.H_)
        val=np.asarray(fb.getValList(),dtype=float)
        col=np.asarray(fb.getColumnN(),dtype=int)
        row=np.asarray(fb.getRowN(),dtype=int)
        return csc_matrix((val,(row,col)),shape=self.shape_)

    def _getPattern(self):
        aux=np.asarray(self.H_.getHpattern())
        #  as flattened list via a pointer
        #  the numberofblocks is also returned (by reference)
        #
        #  the pattern matrix is formatted as
        # row = 1 block : i_begin,j_begin, i_end,j_end,FLAG,entry_size
        # with FLAG=0 for full rank and FLAG=1 for low rank
        #
        # we output a flatten row-major order std::vector

        nr=6
        return np.reshape(aux,(int(aux.size/nr),nr))

    def plotPattern(self):
        data_pattern = self._getPattern()

        patches = []
        p_colors = []
        max_y = data_pattern[:,3].max()
        for i in range(len(data_pattern)):
            height = np.abs(data_pattern[i,0] - data_pattern[i,2])
            width = np.abs(data_pattern[i,1] - data_pattern[i,3])
            y1 = max_y - data_pattern[i,0] - height; x1 = data_pattern[i,1]
            rectangle = Rectangle((x1, y1), width, height)
            patches.append(rectangle)
            p_colors.append(data_pattern[i,4])
        fig = plt.figure()
        ax = fig.add_subplot(111)

        p = PatchCollection(patches, cmap=matplotlib.cm.PiYG, edgecolors='black',alpha=0.4)
        p.set_array(p_colors)
        ax.add_collection(p)
        ax.set_ylim([data_pattern[:,0].min(),data_pattern[:,3].max()])
        ax.set_xlim([data_pattern[:,1].min(),data_pattern[:,2].max()])
        ax.set_aspect('equal')
        #fig.colorbar(p)
        fig.show()
        plt.show(block=True)
        return fig

    # a method constructing an ILU Preconditionner of the H matrix
    def H_ILU_prec(self,fill_factor=5,drop_tol=1e-5):
        #fb = self._getFullBlocks()
        fbinv = self.H_ILU(fill_factor=fill_factor,drop_tol=drop_tol)
        return LinearOperator(self.shape_, fbinv.solve)

    def H_ILU(self,fill_factor=5,drop_tol=1e-5):
        fb= self._getFullBlocks()
        fbILU = spilu(fb,fill_factor=fill_factor,drop_tol=drop_tol)
        return fbILU

    def H_diag(self):
        fb = self._getFullBlocks()
        return fb.diagonal()

    def H_jacobi_prec(self):
        diag = self.H_diag()   # return a nd.array
        overdiag = 1./diag
        return diags(overdiag, dtype=float)

#--------------------------------
#
Lx= 10.
Ly= 10.  # half mesh size in y direction
nx= 31  # number of cell in x direction
ny= 31  # number of cells in y direction

# create cartesian mesh object
mesh = CartesianMesh(Lx, Ly, nx, ny)

# move the coordinates to 3D
coor2D = mesh.VertexCoor
coor = np.zeros((len(coor2D), 3))
for i in range(len(coor2D)):
    coor[i,[0,1]]=coor2D[i,:]

# get the connectivity
conn= mesh.Connectivity


#       "properties" is a vector of size 2 or 3 depending on the kernel
#          It should contain at least
#          - properties[0] = YoungModulus
#          - properties[1] = PoissonRatio
properties = np.array([100, 0.2])


# Hmat parameters
max_leaf_size = 10
eta = 3.
eps_aca = 0.1
kernel = '3DR0'

# create Hmat
myHmat_linear_operator = Hmatrix(kernel,coor,conn,properties,max_leaf_size=100,eta=3,eps_aca=1.e-3)


# compute displacements at given observation points
# [DD1_el1, DD2_el1, DD3_el1, DD1_el2, DD2_el2, DD3_el2]
solution_DDs = np.ones(mesh.NumberOfElts*3)

#test the dot product
dot_res = myHmat_linear_operator._matvec(solution_DDs)

# "obsPts" a flattened list containing the coordinates of the observation points
#          coordinates [ x(1), y(1), z(1), ... ,x(npts), y(npts), z(npts) ]
obsPts = [1.1, 1.2, 0.,
          1,2,4,
          1,3,4]

# "npts" the number of observation points
npts = 2

# global or local DDs?
are_dd_global = False

#compute displacements at given observation points
result  = myHmat_linear_operator.H_.computeDisplacements(solution_DDs,
                                                        obsPts,
                                                        npts,
                                                        properties.flatten(),
                                                        coor.flatten(),
                                                        conn.flatten(),
                                                        are_dd_global)

# get the Hpattern
fig = myHmat_linear_operator.plotPattern()

print("done")