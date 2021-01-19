import numpy as np
import pyparty

from pyparty import Bigwhamio
BigwhamOBJ = Bigwhamio()

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
kernel = "3DR0"
properties = [100,0.2] # Young Modulus , Poisson's ratio
max_leaf_size = 1
eta = 0.
eps_aca = 0.001

# set the object
BigwhamOBJ.set(coor,
      conn,
      kernel,
      properties,
      max_leaf_size,
      eta,
      eps_aca)

# flattened collocation points
mycollp = BigwhamOBJ.getCollocationPoints()
print(mycollp)
print("\n")

# hdot product
tractions = BigwhamOBJ.hdotProduct([1.,2.,3.,4.,5.,6.], False)
print(tractions)
tractions = BigwhamOBJ.hdotProduct([1.,2.,3.,4.,5.,6.], True)
print(tractions)

# BigwhamOBJ.getPermutation()
# BigwhamOBJ.getCompressionRatio()
# BigwhamOBJ.getKernel()
# BigwhamOBJ.getSpatialDimension()
# BigwhamOBJ.matrixSize()
# BigwhamOBJ.getHpattern()
# BigwhamOBJ.getFullBlocks()

