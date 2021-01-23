import numpy as np
import pyparty

from pyparty import Bigwhamio


# Defining the variables:

# coordinates   - const std::vector<double>
# connectivity  - const std::vector<int64_t>
# kernel        - const std::string
# properties    - const std::vector<double>
# max_leaf_size - const int
# eta           - const double
# eps_aca       - const double

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


tractionKernel = "3DR0_traction"
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

stresses = tractionHMAT.computeStresses(mysol, obsPoints, 2, properties, coor, conn)
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
print(mystress[1:6])
print(mystress[7:12])
print(mystress[13:18])

x,y,z,a,b,G,nu = [1.5,1.5,0.,2.5,2.,200.,0.3]
mystress = tractionHMAT.getInfluenceCoe(x,y,z,a,b,G,nu)
print("\n ------------------ \n ")
print(" Stress: \n ")
print("DDx (shear)  -> | sxx, syy, szz, sxy, sxz, syz  | ")
print("DDy (shear)  -> | sxx, syy, szz, sxy, sxz, syz  | ")
print("DDy (normal) -> | sxx, syy, szz, sxy, sxz, syz  | ")
print(mystress[1:6])
print(mystress[7:12])
print(mystress[13:18])