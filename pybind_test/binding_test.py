import numpy as np
import pyparty
# from pyparty import addwithdefaults
from pyparty import Bigwhamio
a = Bigwhamio()

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

a.set(coor,
      conn,
      kernel,
      properties,
      max_leaf_size,
      eta,
      eps_aca)

print(b)
