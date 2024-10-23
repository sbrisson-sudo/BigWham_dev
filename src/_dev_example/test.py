import sys
sys.path.append("/home/ankit/geolab/BigWham/build/")

import bigwhamPybind as bw
a = bw.Real2D(3,3, 0)

a[[0,1]] = 3

print(a[[0,1]])
print(a[[0,2]])

print(a)

print(a.shape(0))

seg0 = bw.Segment0

# import numpy as np

# a = np.zeros((4,4))
# a[0,0] = 1
# print(a)
