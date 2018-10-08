from __future__ import division
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math
import numpy as np
import time


FOURPISQUARE = 4*(math.pi)*(math.pi)
AU = 1.496e+11                               # Earth's distance from the Sun in meters
AUCUBE = AU*AU*AU                            # cube of AU since G*Msun is in units of AU**3/year**2
YEAR = 3.156e+07                             # number of seconds in one year
YEARSQUARE = YEAR*YEAR                       # square one year
GSUNMASS = FOURPISQUARE*AUCUBE/YEARSQUARE    # used in potential energy formula
SUNMASS = 1.991e+30                          # Sun's mass in kg

Tmax = 50       # max number of years to run for
dt = .002       # timestep in years

# c
planet_c_mass = 8.24e24/SUNMASS # Msol
planet_c_dist = .015      # Au
planet_c_vy0 = 14.22    # Au/yr

# e
planet_e_mass = 3.70e24/SUNMASS # Msol
planet_e_dist = .028      # Au
planet_e_vy0 = 10.53    # Au / yr

# f
planet_f_mass = 4.06e24/SUNMASS
planet_f_dist = .037
planet_f_vy0 = 9.220

# g
planet_g_mass = 8.00e24/SUNMASS
planet_g_dist = .045
planet_g_vy0 = 8.362

# star
star_mass = .08  # Msol

sanity = False
if sanity:
    planet_c_mass = .0000010/SUNMASS # solar mass units
    planet_e_mass = .0000010/SUNMASS # solar mass units


a = planet_c_dist
b = planet_e_dist
c = planet_f_dist
d = planet_g_dist

m1 = star_mass
m2 = planet_c_mass
m3 = planet_e_mass
m4 = planet_f_mass
m5 = planet_g_mass

mass_sum = star_mass + planet_c_mass + planet_e_mass + planet_f_mass + planet_g_mass
star_x0 = -(planet_c_dist*planet_c_mass + planet_e_dist*planet_e_mass + planet_f_dist*planet_f_mass + planet_g_dist*planet_g_mass) / mass_sum
planet_c_x0 = (-a*m1 - a*m3 + b*m3 - a*m4 + c*m4 -a*m5 + d*m5) / mass_sum
planet_e_x0 = (-b*m1 + a*m2 - b*m2 - b*m4 + c*m4 - b*m5 + d*m5) / mass_sum
planet_f_x0 = (-c*m1 + a*m2 - c*m2 + b*m3 - c*m3 - c*m5 + d*m5) / mass_sum
planet_g_x0 = (-d*m1 + a*m2 - d*m2 + b*m3 - d*m3 + c*m4 - d*m4) / mass_sum

star_vy0 = -(planet_c_vy0*planet_c_mass + planet_e_vy0*planet_e_mass + planet_f_vy0*planet_f_mass + planet_g_vy0*planet_g_mass)/star_mass

star_y0 = 0
star_vx0 = 0
planet_c_y0 = 0
planet_c_vx0 = 0
planet_e_y0 = 0
planet_e_vx0 = 0
planet_f_y0 = 0
planet_f_vx0 = 0
planet_g_y0 = 0
planet_g_vx0 = 0


# 4pi^2 AU^3/yr^2 * M(msol) * x (AU) / r^3 (AU^3)

def fDerivative(variableList, t):                  # variableList dummy list array since there is more than one differential equation
    vxa = variableList[0]                          # star speed in the x direction
    vya = variableList[1]                          # star speed in the y direction
    vxb = variableList[2]                          # planet 1 speed in the x direction
    vyb = variableList[3]                          # planet 1 speed in the y direction
    vxc = variableList[4]                          # planet 2 speed in the x direction
    vyc = variableList[5]                          # planet 2 speed in the y direction
    vxd = variableList[6]                          # planet 3 speed in the x direction
    vyd = variableList[7]                          # planet 3 speed in the y direction
    vxe = variableList[8]                          # planet 4 speed in the x direction
    vye = variableList[9]                          # planet 4 speed in the y direction

    xa = variableList[10]                           # star x coordinate
    ya = variableList[11]                           # star y coordinate
    xb = variableList[12]                           # planet 1 x coordinate
    yb = variableList[13]                           # planet 1 y coordinate
    xc = variableList[14]                           # planet 2 x coordinate
    yc = variableList[15]                           # planet 2 y coordinate
    xd = variableList[16]                           # planet 3 x coordinate
    yd = variableList[17]                           # planet 3 y coordinate
    xe = variableList[18]                           # planet 4 x coordinate
    ye = variableList[19]                           # planet 4 y coordinate
    #
    # Compute distance parameters
    #
    dabSq = pow(xb-xa, 2) + pow(yb-ya, 2)   # sq of d between star and p1
    dacSq = pow(xc-xa, 2) + pow(yc-ya, 2)   # sq of d between star and p2
    dadSq = pow(xd-xa, 2) + pow(yd-ya, 2)   # sq of d between star and p3
    daeSq = pow(xe-xa, 2) + pow(ye-ya, 2)   # sq of d between star and p2

    dbcSq = pow(xc-xb, 2) + pow(yc-yb, 2)   # square of the distance between planet_1 and the planet
    dbdSq = pow(xd-xb, 2) + pow(yd-yb, 2)   # square of the distance between planet_1 and the planet
    dbeSq = pow(xe-xb, 2) + pow(ye-yb, 2)   # square of the distance between planet_1 and the planet

    dcdSq = pow(xd-xc, 2) + pow(yd-yc, 2)   # square of the distance between Mercury and the planet
    dceSq = pow(xe-xc, 2) + pow(ye-yc, 2)   # square of the distance between Mercury and the planet

    ddeSq = pow(xe-xd, 2) + pow(ye-yd, 2)   # square of the distance between planets d and e

    d12Cu = pow(dabSq, 1.5)   # cube of the distance between the star and Mercury
    d13Cu = pow(dacSq, 1.5)   # cube of the distance between the star and the planet
    d14Cu = pow(dadSq, 1.5)   # cube of the distance between the star and the planet
    d15Cu = pow(daeSq, 1.5)   # cube of the distance between the star and the planet

    d23Cu = pow(dbcSq, 1.5)   # cube of the distance between Mercury and the planet
    d24Cu = pow(dbdSq, 1.5)   # cube of the distance between Mercury and the planet
    d25Cu = pow(dbeSq, 1.5)   # cube of the distance between Mercury and the planet

    d34Cu = pow(dcdSq, 1.5)   # cube of the distance between Mercury and the planet
    d35Cu = pow(dceSq, 1.5)   # cube of the distance between Mercury and the planet

    d45Cu = pow(ddeSq, 1.5)   # cube of the distance between Mercury and the planet



    #
    # Compute the derivatives
    #
    # General form: dv_planet/dt = -FOURPISQUARE * (xplanet - xother) * mass other / seperation of planet and other
    #
    # TRAPPIST-1
    dvx1dt = (-FOURPISQUARE * (xa - xb) * planet_c_mass / d12Cu) - (FOURPISQUARE * (xa - xc) * planet_e_mass / d13Cu) - (FOURPISQUARE*(xa-xd)*planet_f_mass/d14Cu) - (FOURPISQUARE*(xa-xe)*planet_g_mass/d15Cu)
    dvy1dt = -FOURPISQUARE*(ya-yb)*planet_c_mass/d12Cu - FOURPISQUARE*(ya-yc)*planet_e_mass/d13Cu - FOURPISQUARE*(ya-yd)*planet_f_mass/d14Cu - FOURPISQUARE*(ya-ye)*planet_g_mass/d15Cu
    # PLANET B
    dvx2dt = -FOURPISQUARE*(xb-xa)*star_mass/d12Cu - FOURPISQUARE*(xb-xc)*planet_e_mass/d23Cu - FOURPISQUARE*(xb-xd)*planet_f_mass/d24Cu - FOURPISQUARE*(xb-xe)*planet_g_mass/d25Cu
    dvy2dt = -FOURPISQUARE*(yb-ya)*star_mass/d12Cu - FOURPISQUARE*(yb-yc)*planet_e_mass/d23Cu - FOURPISQUARE*(yb-yd)*planet_f_mass/d24Cu - FOURPISQUARE*(yb-ye)*planet_g_mass/d25Cu
    # PLANET C
    dvx3dt = -FOURPISQUARE*(xc-xa)*star_mass/d13Cu - FOURPISQUARE*(xc-xb)*planet_c_mass/d23Cu - FOURPISQUARE*(xc-xd)*planet_f_mass/d34Cu - FOURPISQUARE*(xc-xe)*planet_g_mass/d35Cu
    dvy3dt = -FOURPISQUARE*(yc-ya)*star_mass/d13Cu - FOURPISQUARE*(yc-yb)*planet_c_mass/d23Cu - FOURPISQUARE*(yc-yd)*planet_f_mass/d34Cu - FOURPISQUARE*(yc-ye)*planet_g_mass/d35Cu
    # PLANET D
    dvx4dt = -FOURPISQUARE*(xd-xa)*star_mass/d14Cu - FOURPISQUARE*(xd-xb)*planet_c_mass/d24Cu - FOURPISQUARE*(xd-xc)*planet_e_mass/d34Cu - FOURPISQUARE*(xd-xe)*planet_g_mass/d35Cu
    dvy4dt = -FOURPISQUARE*(yd-ya)*star_mass/d14Cu - FOURPISQUARE*(yd-yb)*planet_c_mass/d24Cu - FOURPISQUARE*(yd-yc)*planet_e_mass/d34Cu - FOURPISQUARE*(yd-ye)*planet_g_mass/d35Cu
    # PLANET E
    dvx5dt = -FOURPISQUARE*(xe-xa)*star_mass/d15Cu - FOURPISQUARE*(xe-xb)*planet_c_mass/d25Cu - FOURPISQUARE*(xe-xc)*planet_e_mass/d35Cu - FOURPISQUARE*(xe-xd)*planet_f_mass/d45Cu
    dvy5dt = -FOURPISQUARE*(ye-ya)*star_mass/d15Cu - FOURPISQUARE*(ye-yb)*planet_c_mass/d25Cu - FOURPISQUARE*(ye-yc)*planet_e_mass/d35Cu - FOURPISQUARE*(ye-yd)*planet_f_mass/d45Cu

    return [dvx1dt, dvy1dt, dvx2dt, dvy2dt, dvx3dt, dvy3dt, dvx4dt, dvy4dt, dvx5dt, dvy5dt, vxa, vya, vxb, vyb, vxc, vyc, vxd, vyd, vxe, vye]    # return the twenty derivatives as a list object containing twelve elements


initialValuesSpeedsPositions = [star_vx0, star_vy0, planet_c_vx0, planet_c_vy0, planet_e_vx0, planet_e_vy0, planet_f_vx0, planet_f_vy0, planet_g_vx0, planet_g_vy0, star_x0, star_y0, planet_c_x0, planet_c_y0, planet_e_x0, planet_e_y0, planet_f_x0, planet_f_y0, planet_g_x0, planet_g_y0]
nTimeSteps = int(Tmax / dt)
timeGrid = np.linspace(0, Tmax, nTimeSteps)           # time grid used by odeint method
eighteenSolution = odeint(fDerivative, initialValuesSpeedsPositions, timeGrid)   # odeint returns a list of values which is the RK4 solution
print eighteenSolution
time.sleep(2)
vx1RK4 = eighteenSolution[:,0]          # vx1 function of time obtained with RK4 solution
vy1RK4 = eighteenSolution[:,1]          # vy1 function of time obtained with RK4 solution
vx2RK4 = eighteenSolution[:,2]          # vx2 function of time obtained with RK4 solution
vy2RK4 = eighteenSolution[:,3]          # vy2 function of time obtained with RK4 solution
vx3RK4 = eighteenSolution[:,4]          # vx3 function of time obtained with RK4 solution
vy3RK4 = eighteenSolution[:,5]          # vy3 function of time obtained with RK4 solution
vx4RK4 = eighteenSolution[:,6]          # vx4 function of time obtained with RK4 solution
vy4RK4 = eighteenSolution[:,7]          # vy4 function of time obtained with RK4 solution
vx5RK4 = eighteenSolution[:,8]          # vx5 function of time obtained with RK4 solution
vy5RK4 = eighteenSolution[:,9]          # vy5 function of time obtained with RK4 solution

x1RK4 = eighteenSolution[:,10]           # x1 function of time obtained with RK4 solution
y1RK4 = eighteenSolution[:,11]           # y1 function of time obtained with RK4 solution
x2RK4 = eighteenSolution[:,12]           # x2 function of time obtained with RK4 solution
y2RK4 = eighteenSolution[:,13]           # y2 function of time obtained with RK4 solution
x3RK4 = eighteenSolution[:,14]           # x3 function of time obtained with RK4 solution
y3RK4 = eighteenSolution[:,15]           # y3 function of time obtained with RK4 solution
x4RK4 = eighteenSolution[:,16]           # x4 function of time obtained with RK4 solution
y4RK4 = eighteenSolution[:,17]           # y4 function of time obtained with RK4 solution
x5RK4 = eighteenSolution[:,18]           # x5 function of time obtained with RK4 solution
y5RK4 = eighteenSolution[:,19]           # x3 function of time obtained with RK4 solution

for x in range(x1RK4.__len__()):
    print x1RK4[x], y1RK4[x]

rad_pos_star, rad_pos_c, rad_pos_e, rad_pos_f, rad_pos_g = [], [], [], [], []

time.sleep(2)
time_graph = []
for step in range(nTimeSteps):
    rad_pos_c.append(math.sqrt(pow(x2RK4[step], 2) + pow(y2RK4[step], 2)))
    rad_pos_e.append(math.sqrt(pow(x3RK4[step], 2) + pow(y3RK4[step], 2)))
    rad_pos_f.append(math.sqrt(pow(x4RK4[step], 2) + pow(y4RK4[step], 2)))
    rad_pos_g.append(math.sqrt(pow(x5RK4[step], 2) + pow(y5RK4[step], 2)))
    time_graph.append(step/(1.0/dt))

plt.plot(time_graph, rad_pos_f, 'g', label="radial position of TRAPPIST-1f")
plt.grid(True)
plt.xlabel("Time (years)")
plt.legend(loc=0)
plt.ylabel("Radial position (Au)")
plt.show()
plt.plot(time_graph, rad_pos_g, 'r', label="radial position of TRAPPIST-1g")
plt.xlabel("Time (years)")
plt.legend(loc=0)
plt.ylabel("Radial position (Au)")
plt.show()