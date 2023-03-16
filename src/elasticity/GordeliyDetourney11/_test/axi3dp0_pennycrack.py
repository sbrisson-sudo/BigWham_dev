from bigwham4py import Hmatrix  # linear operator file
import numpy as np

radius = 1.0
pressure = 1.0

G = 1.0
nu = 0.0
E = (2 * G) * (1 + nu)

max_leaf_size = 20
eta = 5.0
eps_aca = 1e-4

nelts = 100
coor1D = np.linspace(0, radius, nelts + 1)
coor = np.transpose(np.array([coor1D, coor1D * 0.0]))
conn = np.fromfunction(lambda i, j: i + j, (nelts, 2), dtype=np.int_)

# Create H-matrix
kernel = "Axi3DP0"
hmat = Hmatrix(kernel, coor, conn, np.array([E, nu]), max_leaf_size, eta, eps_aca)

col_pts = hmat.getMeshCollocationPoints()

pre_fac = (8 * (1 - nu * nu)) / (np.pi * E)
dd = np.zeros(col_pts.shape)
dd[:, 1] = pre_fac * np.sqrt(radius * radius - np.linalg.norm(col_pts[:, :], axis=1)**2)


# calculate tractions
t = hmat.matvec(dd.flatten())

# print(dd)
# print(t)

t_anal = np.zeros(col_pts.shape)
t_anal[:, 1] = pressure

print("L2 Rel error {}".format(np.linalg.norm(t - t_anal.flatten()) / nelts))
