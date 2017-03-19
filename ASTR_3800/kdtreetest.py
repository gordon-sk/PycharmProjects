from scipy.spatial import kdtree as kd
import time

a = (0,1,2), (0,2,1), (0,3,2), (3,0,1), (2,1,0), (0,1,1.5)
a = kd.KDTree(a)

for x in a.data:
    print x
print

b = (0,0, 0), (0,0, 0)
b = kd.KDTree(b)


time.sleep(1)
print a.data.__len__()
print a.count_neighbors(a, .5)
