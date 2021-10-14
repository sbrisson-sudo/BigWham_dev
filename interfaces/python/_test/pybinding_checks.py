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

import numpy as np
from interfaces.python import bigwhamPybind # linear operator file
from interfaces.python.lib.bigwhamPybind import Bigwhamio  #so file
################################
# TESTING THE LINEAR OPERATOR  #
################################
print(" \n")
print(" \n")
print("TESTING THE LINEAR OPERATOR CLASS \n")
print(" \n")
print(" \n")

kernel = "3DR0_displ"

coor =np.asarray([[-1., -1., 0.],
       [1., -1., 0.],
       [1., 1., 0.],
       [-1., 1., 0.],
       [-1., 2., 0.],
       [1., 2., 0.]])

conn =np.asarray([[0, 1, 2, 3], [3, 2, 5, 4]])

Young = 100
PoissonRatio = 0.2
max_leaf_size = 1
eta = 0.
eps_aca = 0.001


#create an Hdot instance
displHMAT = bigwhamPybind.Hmatrix(kernel, coor, conn, np.array([Young, PoissonRatio]), max_leaf_size, eta, eps_aca)

print("Testing the Hdot product for the displacement HMAT \n")
res = displHMAT._matvec([1.,2.,3.,4.,5.,6.])
print(res)
print(" \n")
print(" \n")

##########################
# DIRECT USE OF BIGWHAM  #
##########################

# Defining the variables:

# coordinates   - const std::vector<double>
# connectivity  - const std::vector<int64_t>
# kernel        - const std::string
# properties    - const std::vector<double>
# max_leaf_size - const int
# eta           - const double
# eps_aca       - const double

print("TESTING THE DIRECT USE OF BIGWHAM \n")
print(" \n")
print(" \n")
coor =[-1., -1., 0.,
       1., -1., 0.,
       1., 1., 0.,
       -1., 1., 0.,
       -1., 2., 0.,
       1., 2., 0.]

conn =[0, 1, 2, 3, 3, 2, 5, 4]

properties = [100,0.2] # Young Modulus , Poisson's ratio

max_leaf_size = 1
eta = 0.
eps_aca = 0.001

displacementKernel = "3DR0_displ"



displacementHMAT = Bigwhamio()
# set the object
displacementHMAT.set(coor,
                 conn,
                 displacementKernel,
                 properties,
                 max_leaf_size,
                 eta,
                 eps_aca)


tractionKernel = "3DR0"
tractionHMAT = Bigwhamio()
# set the object
tractionHMAT.set(coor,
      conn,
      tractionKernel,
      properties,
      max_leaf_size,
      eta,
      eps_aca)

# flattened collocation points
mycollp = tractionHMAT.getCollocationPoints()
print(mycollp)
print("\n")

# hdot product
print("Testing the Hdot product for the tractionHMAT \n")
tractions = tractionHMAT.hdotProduct([1.,2.,3.,4.,5.,6.])
print(tractions)

print("Testing the Hdot product for the displacementHMAT \n")
tractions = displacementHMAT.hdotProduct([1.,2.,3.,4.,5.,6.])
print(tractions)

mysol = [1.,1.,1.,1.,1.,1.]
obsPoints = [-10.,-10.,0., #point 1
             20.,-20.,0.]  #point 2

stresses = tractionHMAT.computeStresses(mysol, obsPoints, 2, properties, coor, conn, True)
print("point 1 ")
print(stresses[1:6])
print("point 2 ")
print(stresses[7:12])

x,y,z,a,b,G,nu = [0.,0.,0.,1.,1.,200.,0.3]
mystress = tractionHMAT.getInfluenceCoe(x,y,z,a,b,G,nu)
print("\n ------------------ \n ")
print(" Stress: \n ")
print("DDx (shear)  -> | sxx, syy, szz, sxy, sxz, syz  | ")
print("DDy (shear)  -> | sxx, syy, szz, sxy, sxz, syz  | ")
print("DDy (normal) -> | sxx, syy, szz, sxy, sxz, syz  | ")
print(mystress[0:6])
print(mystress[6:12])
print(mystress[12:18])

x,y,z,a,b,G,nu = [1.5,1.5,0.,2.5,2.,200.,0.3]
mystress = tractionHMAT.getInfluenceCoe(x,y,z,a,b,G,nu)
print("\n ------------------ \n ")
print(" Stress: \n ")
print("DDx (shear)  -> | sxx, syy, szz, sxy, sxz, syz  | ")
print("DDy (shear)  -> | sxx, syy, szz, sxy, sxz, syz  | ")
print("DDy (normal) -> | sxx, syy, szz, sxy, sxz, syz  | ")
print(mystress[0:6])
print(mystress[6:12])
print(mystress[12:18])
print(mystress)

mydisplacements = displacementHMAT.computeDisplacements(mysol, obsPoints, 2, properties, coor, conn, True)
print("point 1 ")
print(mydisplacements[0:3])
print("point 2 ")
print(mydisplacements[3:7])
print("done")