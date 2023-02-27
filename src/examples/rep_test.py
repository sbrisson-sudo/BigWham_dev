## Assuming build/interfaces/python is aleady in python path

from bigwham4py import Hmatrix  # linear operator file
import numpy as np
import subprocess

subprocess.call("python generate_penny_mesh.py 0.5", shell=True)

E = 1.e3
nu = 0.2

max_leaf_size = 1
eta = 0.
eps_aca = 0.001

coord = np.load("mesh_coords.npy")
conn = np.load("mesh_conn.npy")

# Create H-matrix
kernel = "3DT0"
hmat = Hmatrix(kernel, coord, conn, np.array([E, nu]), max_leaf_size, eta,
               eps_aca)
