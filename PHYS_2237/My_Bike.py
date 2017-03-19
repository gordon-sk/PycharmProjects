import matplotlib

matplotlib.use('TkAgg')  # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt  # get matplotlib plot functions
import sys  # used to get the number of command line arguments
import argparse  # argument parser library
import numpy as np  # numerical functions library used by python
import math as mp  # used for the definition of pi
from scipy.integrate import odeint  # import only this single method for solving differential equations

vTrans = 7.0    # m/s
theta = 0     # degrees
v0 = 0.0        # m/s
maxT = 200.0    # seconds
del_t = 2.0     # seconds
Power = 400.0   # Watts
cFactor = 1.0   # Unitless Constant
rho = 1.2       # kg/m^3
mass = 70.0     # kg (rider bike system)
area = .33      # m^2

v=0             # Star of the show
t=0             # Other star

theta_rad = mp.radians(theta)
theta_deg = theta

if v0<0:        # Keeping them honest
    print("Fuck off with that negative initial velocity")
    exit()

# Writing our diff eq

grav_factor = 9.81*mp.sin(theta_rad)
res_factor = (cFactor*rho*area)/(2*mass)

F0 = Power/vTrans   # Newtons

tplot, vplot = [], []

###
#while t<maxT:
 #   res_factor1 = res_factor*v*v
  #  if v <= vTrans:
   #     P = F0/mass
    #else:
    #    P = Power/(mass*v)
    #print("v: " + str(v) + "            t:" + str(t) + "            P: " + str(P))
    #dvdt = P - res_factor1 - grav_factor
    #v = dvdt*t
    #t += .05
    #tplot.append(t)
    #vplot.append(v)
###

while t<maxT:
    dvdt = F0/mass - res_factor*v*v - grav_factor
    v = dvdt*t
    t += .05
    print(str(v) + "       " + str(t))
    tplot.append(t)
    vplot.append(v)

plt.plot(tplot, vplot, 'ro')

plt.xlabel('Time (s)')  # add axis labels
plt.ylabel('Velocity (m/s)')
plt.title('Motion Along A Hill With Quadratic Air Resistance')
plt.grid(True)
plt.show()
