import sys
import numpy as np


L = 1.0 # domain dimension 
N = 100 # number of elements

points = np.array([np.linspace(0., L, N+1), np.zeros(N+1)])
conn = np.array([np.arange(N), np.arange(1,N+1)])


np.save("mesh_coor_seg", np.asfortranarray(points.astype(np.float64)))
np.save("mesh_conn_seg", np.asfortranarray(conn.astype(np.int32)))

print("Mesh written in mesh_coor_seg.npy and mesh_conn_seg.npy")