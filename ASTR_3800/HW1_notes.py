import matplotlib                       # For plotting in Python
import matplotlib.pyplot as mp          #
import math                             # for doing math
import scipy.integrate as integrate     # for doing harder math
matplotlib.use("TkAgg")

y = .01
color, time = [1,2,3,4,5], [1,4,9,16,25]
mp.plot(time, color, "r")
for x in range(color.__len__()-1):
    mp.annotate("Point " + str(x), xy=[time[x], color[x]+.5])
mp.show()

print "fuck you"