# 2D Griffith crack

## 2DP1 kernel

```python
# %%
import sys
import os
home = os.environ["HOME"]
sys.path.append(home + "/geolab/dev_bigwham/build/interfaces/python")
from hmatrix import Hmatrix
import numpy as np
from scipy.sparse.linalg import gmres

# %% Material properties
G = 1.0
nu = 0.25
E = (2 * G) * (1 + nu)

# %% Mesh
a = 10
nelts = 100
coor1D = np.linspace(-a, a, nelts + 1)
coor = np.transpose(np.array([coor1D, coor1D * 0.0]))
conn = np.fromfunction(lambda i, j: i + j, (nelts, 2), dtype=np.int_)

# H-matrix parameter
max_leaf_size=100
eta=3.
eps_aca=1.e-4

# %%
# Hmatrix
kernel = "2DP1"
elas_prop = np.array([E, nu])
h = Hmatrix(kernel, coor, conn, elas_prop, max_leaf_size, eta, eps_aca)

# get augemented collocation points
col_pts_aug = np.zeros((conn.shape[0] * 2, conn.shape[1]))
c = 0
col_pts_aug[c, :] = coor[0, :]
c += 1
for i in range(1, coor.shape[0] -  1):
    col_pts_aug[c, :] = coor[i, :]
    c += 1
    col_pts_aug[c, :] = coor[i, :]
    c += 1
col_pts_aug[c, :] = coor[i+1, :]

# %%
col_pts = h.getMeshCollocationPoints()

# %%
t = np.ones(h.shape[0])
t[0::2] = 0.
d = gmres(h, t, tol=1e-6)[0]

dd = d.reshape((-1, 2))
# %%

import matplotlib.pyplot as plt
plt.plot(col_pts_aug[:, 0], dd, "-*")

# %%

```

## 2DP0 kernel

```python
# %%
import sys
import os
home = os.environ["HOME"]
sys.path.append(home + "/geolab/dev_bigwham/build/interfaces/python")
from hmatrix import Hmatrix
import numpy as np
from scipy.sparse.linalg import gmres

# %% Material properties
G = 1.0
nu = 0.25
E = (2 * G) * (1 + nu)

# %% Mesh
a = 10
nelts = 100
coor1D = np.linspace(-a, a, nelts + 1)
coor = np.transpose(np.array([coor1D, coor1D * 0.0]))
conn = np.fromfunction(lambda i, j: i + j, (nelts, 2), dtype=np.int_)

# H-matrix parameter
max_leaf_size=100
eta=3.
eps_aca=1.e-4

# %%
# Hmatrix
kernel = "2DP0"
elas_prop = np.array([E, nu])
h = Hmatrix(kernel, coor, conn, elas_prop, max_leaf_size, eta, eps_aca)

# %%
col_pts = h.getMeshCollocationPoints()

# %%
t = np.ones(h.shape[0])
t[0::2] = 0.
d = gmres(h, t, tol=1e-6)[0]

dd = d.reshape((-1, 2))
# %%

import matplotlib.pyplot as plt
plt.plot(col_pts[:, 0], dd, "-*")

# %%
```