import periodic_kdtree as kd
import time
bounding = [3, 3, 3]

a = (0,0,0), (0,1,2), (0,2,1), (0,3,2), (3,0,1), (2,1,0), (0,1,1.5)
a = kd.PeriodicCKDTree(bounding, a)

bounds = [0, 0, 0]

print a.data
print

b = (0,0, 0), (0, 0, 0)
b = kd.PeriodicCKDTree(bounding, b)
print b.data

time.sleep(.1)
print
print a.query_ball_point([0, 0, 0], r=2.0)
# tree.query(point, range)
# returns ([r1, r2], [index of point related to r1 in a, "" r2 ""])

