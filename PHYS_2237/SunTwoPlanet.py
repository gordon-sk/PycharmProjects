#
# Program 4.12 Solution for the orbits of a planet in a binary star system orbiting about their center-of-mass (sunTwoPlanetsSystem_Chapter4V1.py)
#
#              Give the command  python sunTwoPlanetsSystem_Chapter4V1.py -h  to get help command on input parameters
#
#              The first mass is the Sun by default.  The second mass is the Earth by default.  The third mass is Jupiter by default.
#              The main use is to change the size and position of Jupiter, to see the effect on Earth's orbit as in the textbook figures 4.12-4.13
#
#              The input arguments to this program are the following, assuming that the coordinate system origin is the cm point and that there is zero total momentum
#                   1) mass1 of the sun in Solar mass units (defaults to one Solar mass; the sun's position and speed are computed to balance the positions and speeds of the two planets)
#                   2) mass2 of the first planet in Solar mass units  (defaults to the Earth's mass)
#                   3) initial x20 position of the first planet (in AU, > 0) with respect to the cm point (defaults to the Earth's aphelion position)
#                   4) initial speed vx20 of the first planet (in units of AU/year, defaults to the Earth's speed at its aphelion position)
#                   5) mass2 of the second planet in Solar mass units (defaults to Jupiter's mass)
#                   6) initial x30 position of the second planet (in AU, defaults to the Jupiter's aphelion position
#                   7) initial speed vx30 of the second planet (in units of AU/year, defaults to the Jupiter's speed at its aphelion position)
#                   8) maximum time (years)
#                   9) time step (years)
#
#              The coupled differential equations are for the motions of the two stars in a Cartesian coordinate system whose (x,y) origin is the cm point
#                   The first star is assumed to start at the perihelion of its orbit
#                   The total linear momentum of the system is zero, meaning that the cm point does not move
#                   There are twelve differential equations to be solved for: 1) (vx1, vy1), 2) (vx2, vy2), 3) (vx3,vy3), 4) (x1, y1), 6) (x2, y2), 6) (x3, y3)
#                   1) dvx1/dt = ax1 = -4*pi*pi*mass2*(x2-x1)/d12**3 - 4*pi*pi*mass3*(x3-x1)/d13**3   where d12 is the separation distance between the Sun and the first planet, and d13 is the separation between the Sun and the second planet
#                   2) dvy1/dt = ay1 = -4*pi*pi*mass2*(y2-y1)/d12**3 - 4*pi*pi*mass3*(y3-y1)/d13**3   where d12 is the separation distance between the Sun and the first planet, and d13 is the separation between the Sun and the second planet
#                   3) dvx2/dt = ax2 = -4*pi*pi*mass1*(x1-x2)/d12**3 - 4*pi*pi*mass3*(x3-x2)/d23**3   where d12 is the separation distance between the Sun and the first planet, and d23 is the separation between the two planets
#                   4) dvy2/dt = ay2 = -4*pi*pi*mass1*(y1-y2)/d12**3 - 4*pi*pi*mass3*(y3-y2)/d23**3   where d12 is the separation distance between the Sun and the first planet, and d23 is the separation between the two planets
#                   5) dvx3/dt = ax3 = -4*pi*pi*mass1*(x1-x3)/d13**3 - 4*pi*pi*mass2*(x2-x3)/d23**3   where d13 is the separation distance between the Sun and the second planet and d23 is the separation between thet two planets
#                   6) dvy3/dt = ay3 = -4*pi*pi*mass1*(y1-y3)/d13**3 - 4*pi*pi*mass2*(y2-y3)/d23**3   where d12 is the separation distance between the Sun and the second planet and d23 is the separation between thet two planets
#                   7) dx1/dt = vx1    Sun velocity component in the x direction
#                   8) dy1/dt = vy1    Sun velocity component in the y direction
#                   9) dx2/dt = vx2    First planet velocity component in the x direction
#                   10) dy2/dt = vy2   First planet velocity component in the y direction
#                   11) dx3/dt = vx3   Second planet velocity component in the x direction
#                   12) dy3/dt = vy3   Second planet velocity component in the y direction
#
import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt     # get matplotlib plot functions
import sys                          # used to get the number of command line arguments
import  argparse                    # argument parser library
import numpy as np                  # numerical functions library used by python
import math as mp                   # used for the exponential function
from scipy.integrate import odeint  # import only this single method for solving differential equations
import matplotlib.animation as animation  # animation library (see examples to learn how to use it)

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()
# USE APHELION
parser.add_argument('--mass1', default=1.0, type=float, help="Mass of the sun in Solar mass units; default 1.0")
parser.add_argument('--mass2', default=1.652e-7, type=float, help="Mass of first planet in Solar mass units; default 3.0035e-6")
parser.add_argument('--aDistance', default=.466697, type=float, help="Initial first planet to Sun distance in AU; default 1.017")
parser.add_argument('--vy20', default=8.192, type=float, help="Initial first planet speed vy20 in AU/Y; default 6.1773")
parser.add_argument('--mass3', default=.0009543, type=float, help="Mass of planet in Solar mass units; default 9.5429e-04")
parser.add_argument('--bDistance', default=5.45, type=float, help="Initial second planet to Sun distance in AU; default 5.45")
parser.add_argument('--vy30', default=2.6274, type=float, help="Initial second planet speed vy30 in AU/Y; default 2.6274")
parser.add_argument('--maxT', default=400.0, type=float, help="Maximum time range in years; default 40")
parser.add_argument('--deltaT', default=0.0003, type=float, help="Iteration time step in years; default 0.005")
parser.add_argument('--thirdLaw', action='store_true', help="Use Newton's Third Law to compute mass2 accelerations default False")
parser.add_argument('--zeroMomentum', action='store_true', help="Enforce zero total momentum condition, default False")
parser.add_argument('--fixCM', action='store_true', help="Enforce fixed center of mass condition, default False")
parser.add_argument('--verbose', action='store_true', help="Give extra printout about orbit pattern recognition results; default False")
parser.add_argument('--nFrames', default=500, type=int, help="Number of frames in the video; default 500")
parser.add_argument('--nInterval', default=200, type=int, help="Time difference between video frames in milliseconds; default 200")
parser.add_argument('--noAnimation', action='store_true', help="Do not do any animation call; default False")
parser.add_argument('--doMovie', action='store_true', help="Do the movie output, requires ffmpeg library; default False")
parser.add_argument('--radialPlot', action='store_true', help="Plot the planet orbits as r(t); default False")
parser.add_argument('--orbitSun', action='store_true', help="Plot the planet orbits about the Sun instead of cm; default False")
parser.add_argument('--sunPlot', action='store_true', help="Plot the Sun's about the center-of-mass; default False")

#
# Useful constant definitions
#
FOURPISQUARE = 4*(np.pi)*(np.pi)
AU = 1.496e+11                               # Earth's distance from the Sun in meters
AUCUBE = AU*AU*AU                            # cube of AU since G*Msun is in units of AU**3/year**2
YEAR = 3.156e+07                             # number of seconds in one year
YEARSQUARE = YEAR*YEAR                       # square one year
GSUNMASS = FOURPISQUARE*AUCUBE/YEARSQUARE    # used in potential energy formula
SUNMASS = 1.991e+30                          # Sun's mass in kg

args = parser.parse_args()
#
# Get the input parameters from the command line
#
numberOfArguments = len(sys.argv)
if(numberOfArguments == 1):
    print "\n All the default parameter choices are used"
    # there is always at least one argument and that first one is the name of the python script itself

#
# Assign the variables in the program to the variable option names in the argument list
#


verbose = args.verbose
orbitSun = args.orbitSun
thirdLaw = args.thirdLaw
zeroMomentum = args.zeroMomentum
fixCM = args.fixCM
sunPlot = args.sunPlot

mass1 = args.mass1

mass2 = args.mass2
aDistance = args.aDistance
if(aDistance < 0.0):
    print "\n Cannot have the initial separation distance", aDistance, " be less than zero"
    print "\n Program is exiting\n"
    exit(0)

vy20 = args.vy20
if(vy20 < 0.0):
    print "\n Cannot have the initial velocity component", vy20, " be less than zero"
    print "\n Program is exiting\n"
    exit(0)

mass3 = args.mass3
bDistance = args.bDistance
if(bDistance <= 0.0):
    print "\n Cannot have the initial position", bDistance, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(0)

vy30 = args.vy30
if(vy30 <= 0.0):
    print "\n Cannot have the initial velocity omponent", vy30, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(0)

maximumTime = args.maxT
timeStep = args.deltaT
nFrames = args.nFrames
nInterval = args.nInterval
noAnimation = args.noAnimation # retrieve the choice to not do the animation
doMovie = args.doMovie  # retrieve the choice to make the mp4 movie
radialPlot = args.radialPlot  # retrieve the choice to plot the planet orbiting r(t)

mass1 = 1.0

mass2 = 3.3011e23/SUNMASS
aDistance = .466697 # Au
vy20 = 8.192

# JUPITER
mass3 = 80 * 1898.19e24/SUNMASS
bDistance = 5.45492 # Au
vy30 = 2.6224   # au/yr

planet = "Venus"
if planet == "Venus":
    mass3 = 0.000002447
    bDistance = 0.728213
    vy30 = 7.3339108

massSum = mass1 + mass2 + mass3
x10 = -(aDistance*mass2 + bDistance*mass3)/massSum
x20 = (aDistance*(mass1 + mass3) - bDistance*mass3)/massSum
x30 = (bDistance*(mass1 + mass2) - aDistance*mass2)/massSum

vy10 = -vy20*mass2/mass1 - vy30*mass3/mass1

print "these interesting values...:"
print x10
print x20
print x30
print vy10
print

y10 = 0
y20 = 0
y30 = 0
vx10 = 0.0
vx20 = 0.0
vx30 = 0.0



print "\n Sun plus two planets system calculation"
print "   Sun's mass = ", mass1, " solar masses"
print "   The Sun is at position (x10, y10) = (", x10, ", ", y10, ") in AU units, with velocity (vx10, vy10) = (", vx10, ", ", vy10, ") in AU/Y units"
print "   First planet masss = ", mass2, " solar masses, with (x20,y20) = (",x20,",",y20,")"
print "   First planet velocity (vx20, vy20) = (", vx20, ", ", vy20, ") in AU/Y units"
print "   Second planet masss = ", mass3, " solar masses, with (x30,y30) = (",x30,",",y30,")"
print "   Second planet velocity (vx30, vy30) = (", vx30, ", ", vy30, ") in AU/Y units"

print "   Maximum time = ", maximumTime, " years,  with in iteration time step = ", timeStep, " years"
print " "

#                   1) dvx1/dt = ax1 = -4*pi*pi*mass2*(x2-x1)/d12**3 - 4*pi*pi*mass3*(x3-x1)/d13**3   where d12 is the separation distance between the two stars and d13 is the separation between the first star and the planet
#                   2) dvy1/dt = ay1 = -4*pi*pi*mass2*(y2-y1)/d12**3 - 4*pi*pi*mass3*(y3-y1)/d13**3   where d12 is the separation distance between the two stars and d13 is the separation between the first star and the planet
#                   3) dvx2/dt = ax2 = -4*pi*pi*mass1*(x1-x2)/d12**3 - 4*pi*pi*mass3*(x3-x2)/d23**3   where d12 is the separation distance between the two stars and d23 is the separation between the second star and the planet
#                   4) dvy2/dt = ay2 = -4*pi*pi*mass1*(y1-y2)/d12**3 - 4*pi*pi*mass3*(y3-y2)/d23**3   where d12 is the separation distance between the two stars and d23 is the separation between the second star and the planet
#                   5) dvx3/dt = ax3 = -4*pi*pi*mass1*(x1-x3)/d13**3 - 4*pi*pi*mass2*(x2-x3)/d23**3   where d13 is the separation distance between the first star and the planet, and d23 is the separation between the second star and the planet
#                   6) dvy3/dt = ay3 = -4*pi*pi*mass1*(y1-y3)/d13**3 - 4*pi*pi*mass2*(y2-y3)/d23**3   where d12 is the separation distance between the two stars and d23 is the separation between the second star and the planet
#                   7) dx1/dt = vx1    Sun velocity component in the x direction
#                   8) dy1/dt = vy1    Sun velocity component in the y direction
#                   9) dx2/dt = vx2    First planet velocity component in the x direction
#                   10) dy2/dt = vy2   First planet velocity component in the y direction
#                   11) dx3/dt = vx3   Second planet velocity component in the x direction
#                   12) dy3/dt = vy3   Second planet velocity component in the y direction
#
printDebug = False
d23MinimumSeparation = 1.e6
d23MinimumTime = 0.0
d12MinimumSeparation = 1.e6
d12MinimumTime = 0.0
d12MaximumSeparation = 0.0
d12MaximumTime = 0.0
def fDerivative(variableList, t):                  # variableList dummy list array since there is more than one differential equation
    global mass1, mass2, FOURPISQUARE, printDebug, d23MinimumSeparation, d23MinimumTime, d12MinimumSeparation, d12MinimumTime, d12MaximumSeparation, d12MaximumTime
    vx1 = variableList[0]                          # Sun speed in the x direction
    vy1 = variableList[1]                          # Sun speed in the y direction
    vx2 = variableList[2]                          # First planet speed in the x direction
    vy2 = variableList[3]                          # First planet speed in the y direction
    vx3 = variableList[4]                          # Second planet speed in the x direction
    vy3 = variableList[5]                          # Second planet speed in the y direction
    x1 = variableList[6]                           # Sun x coordinate
    y1 = variableList[7]                           # Sun y coordinate
    x2 = variableList[8]                           # First planet x coordinate
    y2 = variableList[9]                           # First planet y coordinate
    x3 = variableList[10]                          # Second planet x coordinate
    y3 = variableList[11]                          # Second planet y coordinate
    #
    # Compute distance parameters
    #
    d12Sq = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) # square of the distance between the Sun and the first planet
    d13Sq = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) # square of the distance between the Sun and the second planet
    d23Sq = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) # square of the distance between the two planets

    d12Cu = d12Sq*mp.sqrt(d12Sq)   # cube of the distance between the Sun and the first planet
    d13Cu = d13Sq*mp.sqrt(d13Sq)   # cube of the distance between the Sun and the second planet
    d23Cu = d23Sq*mp.sqrt(d23Sq)   # cube of the distance between the two planets
    #
    # Compute the derivatives
    #
    dvx1dt = -FOURPISQUARE*(x1-x2)*mass2/d12Cu - FOURPISQUARE*(x1-x3)*mass3/d13Cu
    dvy1dt = -FOURPISQUARE*(y1-y2)*mass2/d12Cu - FOURPISQUARE*(y1-y3)*mass3/d13Cu
    dvx2dt = -FOURPISQUARE*(x2-x1)*mass1/d12Cu - FOURPISQUARE*(x2-x3)*mass3/d23Cu
    dvy2dt = -FOURPISQUARE*(y2-y1)*mass1/d12Cu - FOURPISQUARE*(y2-y3)*mass3/d23Cu
    dvx3dt = -FOURPISQUARE*(x3-x1)*mass1/d13Cu - FOURPISQUARE*(x3-x2)*mass2/d23Cu
    dvy3dt = -FOURPISQUARE*(y3-y1)*mass1/d13Cu - FOURPISQUARE*(y3-y2)*mass2/d23Cu

    return [dvx1dt, dvy1dt, dvx2dt, dvy2dt, dvx3dt, dvy3dt, vx1, vy1, vx2, vy2, vx3, vy3]    # return the twelve derivatives as a list object containing twelve elements

nTimeSteps = int(maximumTime/timeStep)
# obtain the differential equation solutions using the odeint method from the ScyPy library
timeGrid = np.linspace(0, maximumTime,nTimeSteps)           # time grid used by odeint method
initialValuesSpeedsPositions = [vx10, vy10, vx20, vy20, vx30, vy30, x10, y10, x20, y20, x30, y30]   # starting values of vx1, vy1, vx2, vy2, vx3, vy3, x1, y1, x2, y2, x3, y3 for the iteration
twelveSolution = odeint(fDerivative, initialValuesSpeedsPositions, timeGrid)   # odeint returns a list of values which is the RK4 solution

vx1RK4 = twelveSolution[:,0]          # vx1 function of time obtained with RK4 solution
vy1RK4 = twelveSolution[:,1]          # vy1 function of time obtained with RK4 solution
vx2RK4 = twelveSolution[:,2]          # vx2 function of time obtained with RK4 solution
vy2RK4 = twelveSolution[:,3]          # vy2 function of time obtained with RK4 solution
vx3RK4 = twelveSolution[:,4]          # vx3 function of time obtained with RK4 solution
vy3RK4 = twelveSolution[:,5]          # vy3 function of time obtained with RK4 solution

x1RK4 = twelveSolution[:,6]           # x1 function of time obtained with RK4 solution
y1RK4 = twelveSolution[:,7]           # y1 function of time obtained with RK4 solution
x2RK4 = twelveSolution[:,8]           # x2 function of time obtained with RK4 solution
y2RK4 = twelveSolution[:,9]           # y2 function of time obtained with RK4 solution
x3RK4 = twelveSolution[:,10]          # x3 function of time obtained with RK4 solution
y3RK4 = twelveSolution[:,11]          # y3 function of time obtained with RK4 solution

##############################################################################
##############################################################################
##############################################################################
# I need radius as function of x and y... find out what the first time that the radius starts decreasing?
# if i started at perihilion, for example, then I cycle around and say "what is the first time that the radius starts decreasing" and that
# is where the aphelion is... find x,y, take inverse tangent and find out the angle
# time should be increasing in uniform steps, and angle I calculate should be increasing slowly
#


x1 = x1RK4
y1 = y1RK4
x2 = x2RK4
y2 = y2RK4

print x2.__len__()

print "working on these lists rn..."
years = []

steps_per_year = int(1/timeStep)
print "we have", steps_per_year, "steps per year"

rad_pos_b = [0 for x in range(nTimeSteps)]
rad_pos_c = [0 for x in range(nTimeSteps)]

per_graph, time_graph = [], []

for step in range(nTimeSteps):
    rad_pos_b[step] = (mp.sqrt(pow(x2RK4[step], 2) + pow(y2RK4[step], 2)))
    try:
        if rad_pos_b[step] < rad_pos_b[step-1] and rad_pos_b[step-1] > rad_pos_b[step-2]:
            print x2RK4[step], y2RK4[step]
            per_graph.append(mp.degrees(mp.atan2(x1RK4[step]-x2RK4[step], y1RK4[step]-y2RK4[step])))
            time_graph.append(step*maximumTime*1.0/nTimeSteps)
    except IndexError:
        pass
print per_graph.__len__()
print

plt.plot(time_graph, per_graph, 'ro', label="Aphelion Angle Precession")
#plt.ylim(min(perhils_merc) + 5, max(perhils_merc) + 5)
plt.xlim(0, maximumTime)
plt.xlabel("Years")
plt.ylabel("Angle of aphelion")
plt.grid(True)
plt.show()


exit()


##############################################################################
##############################################################################
####################################################################
#
# Set up the radial position and the total velocity arrays for the first planet (Earth default)
#
radialPosition2 = []
totalVelocity2= []
radialPosition12 = []
x12Position = []
y12Position = []

nTimeStep = 0
while nTimeStep < nTimeSteps:
    radialPosition2.append(mp.sqrt(x2RK4[nTimeStep]*x2RK4[nTimeStep] + y2RK4[nTimeStep]*y2RK4[nTimeStep]))
    totalVelocity2.append(mp.sqrt(vx2RK4[nTimeStep]*vx2RK4[nTimeStep] + vy2RK4[nTimeStep]*vy2RK4[nTimeStep]))
    r12 = mp.sqrt((x2RK4[nTimeStep]-x1RK4[nTimeStep])*(x2RK4[nTimeStep]-x1RK4[nTimeStep]) + (y2RK4[nTimeStep]-y1RK4[nTimeStep])*(y2RK4[nTimeStep]-y1RK4[nTimeStep]))
    radialPosition12.append(r12)
    x12Position.append(x2RK4[nTimeStep] - x1RK4[nTimeStep])
    y12Position.append(y2RK4[nTimeStep] - y1RK4[nTimeStep])
    nTimeStep += 1

meanRadialPosition2 = np.mean(radialPosition2)
stdRadialPosition2 = np.std(radialPosition2)
meanTotalVelocity2 = np.mean(totalVelocity2)
stdTotalVelocity2 = np.std(totalVelocity2)

meanRadialPosition12 = np.mean(radialPosition12)
stdRadialPosition12 = np.std(radialPosition12)

#
# Set up the radial position and the total velocity arrays for the second planet (Jupiter default)
#
radialPosition3 = []
totalVelocity3= []
radialPosition13 = []
x13Position = []
y13Position = []

nTimeStep = 0
while nTimeStep < nTimeSteps:
    radialPosition3.append(mp.sqrt(x3RK4[nTimeStep]*x3RK4[nTimeStep] + y3RK4[nTimeStep]*y3RK4[nTimeStep]))
    totalVelocity3.append(mp.sqrt(vx3RK4[nTimeStep]*vx3RK4[nTimeStep] + vy3RK4[nTimeStep]*vy3RK4[nTimeStep]))
    r13 = mp.sqrt((x3RK4[nTimeStep]-x1RK4[nTimeStep])*(x3RK4[nTimeStep]-x1RK4[nTimeStep]) + (y3RK4[nTimeStep]-y1RK4[nTimeStep])*(y3RK4[nTimeStep]-y1RK4[nTimeStep]))
    radialPosition13.append(r13)
    x13Position.append(x3RK4[nTimeStep] - x1RK4[nTimeStep])
    y13Position.append(y3RK4[nTimeStep] - y1RK4[nTimeStep])
    nTimeStep += 1

meanRadialPosition3 = np.mean(radialPosition3)
stdRadialPosition3 = np.std(radialPosition3)
meanTotalVelocity3 = np.mean(totalVelocity3)
stdTotalVelocity3 = np.std(totalVelocity3)

meanRadialPosition13 = np.mean(radialPosition13)
stdRadialPosition13 = np.std(radialPosition13)

if(radialPlot):
    plt.figure(1)       # start a figure for a single plot of the orbit
    plt.subplot(211)    # this sets the upper half plot for the first planet
    radiusString2 = '%s %5.4f %s %5.4f %s' % ('Radius from c.m. = ', meanRadialPosition2, ' +/- ', stdRadialPosition2, ' AU')
    plt.plot(timeGrid, radialPosition2, label=radiusString2)    # plot as a continuous line
    radiusString12 = '%s %5.4f %s %5.4f %s' % ('Radius from Sun = ', meanRadialPosition12, ' +/- ', stdRadialPosition12, ' AU')
    plt.plot(timeGrid, radialPosition12, label=radiusString12)    # plot as a continuous line
    plt.ylabel('Radius (AU)')
    rMinimum = min(radialPosition2)
    rMaximum = max(radialPosition2)
    yMinimum = 0.995*rMinimum
    yMaximum = rMaximum + 0.6*(rMaximum - rMinimum)
    plt.xlim(0, timeGrid[nTimeSteps -1])
    plt.ylim(yMinimum, yMaximum)
    plt.title('Radius of the First Planet in a Sun + Two Planet System')
    plt.grid(True)
    plt.legend(loc=1)

    plt.subplot(212)    # this sets the lower half plot for second planet
    radiusString3 = '%s %5.4f %s %5.4f %s' % ('Radius from c.m. = ', meanRadialPosition3, ' +/- ', stdRadialPosition3, ' AU')
    plt.plot(timeGrid, radialPosition3, label=radiusString3)    # plot as a continuous line
    radiusString13 = '%s %5.4f %s %5.4f %s' % ('Radius from Sun = ', meanRadialPosition13, ' +/- ', stdRadialPosition13, ' AU')
    plt.plot(timeGrid, radialPosition13, label=radiusString13)    # plot as a continuous line
    plt.xlabel('Time (years)')
    plt.ylabel('Radius (AU)')
    rMinimum3 = min(radialPosition3)
    rMinimum13 = min(radialPosition13)
    rMaximum3 = max(radialPosition3)
    rMaximum13 = max(radialPosition13)
    yMinimum = 0.985*min(rMinimum3, rMinimum13)
    yMaximum = 1.04*max(rMaximum3, rMaximum13)
    plt.xlim(0, timeGrid[nTimeSteps -1])
    plt.ylim(yMinimum, yMaximum)
    plt.title('Radius of the Second Planet in a Sun + Two Planet System')
    plt.grid(True)
    plt.legend(loc=1)

    plt.show()          # show the complete figure
    exit()

minimumPositionX1 = 1.5*min(x1RK4)
minimumPositionX2 = 1.5*min(x2RK4)
minimumPositionX3 = 1.5*min(x3RK4)

maximumPositionX1 = 1.5*max(x1RK4)
maximumPositionX2 = 1.5*max(x2RK4)
maximumPositionX3 = 1.5*max(x3RK4)

minimumPositionX = min(minimumPositionX1, minimumPositionX2, minimumPositionX3)
maximumPositionX = max(maximumPositionX1, maximumPositionX2, maximumPositionX3)

minimumPositionY1 = 1.5*min(y1RK4)
minimumPositionY2 = 1.5*min(y2RK4)
minimumPositionY3 = 1.5*min(y3RK4)

maximumPositionY1 = 1.5*max(y1RK4)
maximumPositionY2 = 1.5*max(y2RK4)
maximumPositionY3 = 1.5*max(y3RK4)

minimumPositionY = min(minimumPositionY1, minimumPositionY2, minimumPositionY3)
maximumPositionY = max(maximumPositionY1, maximumPositionY2, maximumPositionY3)

if(orbitSun):
    minimumPositionX13 = 1.5*min(x13Position)
    minimumPositionX12 = 1.5*min(x12Position)

    maximumPositionX13 = 1.5*max(x13Position)
    maximumPositionX12 = 1.5*max(x12Position)

    minimumPositionX = min(minimumPositionX13, minimumPositionX12)
    maximumPositionX = max(maximumPositionX13, maximumPositionX12)

    minimumPositionY13 = 1.5*min(y13Position)
    minimumPositionY12 = 1.5*min(y12Position)

    maximumPositionY13 = 1.5*max(y13Position)
    maximumPositionY12 = 1.5*max(y12Position)

    minimumPositionY = min(minimumPositionY13, minimumPositionY12)
    maximumPositionY = max(maximumPositionY13, maximumPositionY12)


print "\n Minimum x position for the plot = ", minimumPositionX,
print ",  maximum x position for the plot = ", maximumPositionX
print " Minimum y position for the plot = ", minimumPositionY,
print ",  maximum y position for the plot = ", maximumPositionY
print " Minimum d12 separation distance = ", d12MinimumSeparation, " AU at time = ", d12MinimumTime, " Years"
print " Maximum d12 separation distance = ", d12MaximumSeparation, " AU at time = ", d12MaximumTime, " Years"
print " Minimum d23 separation distance = ", d23MinimumSeparation, " AU at time = ", d23MinimumTime, " Years"

fig = plt.figure(1)       # start a figure for plots of the two orbits
ax = fig.add_subplot(111)
if(orbitSun):
    plt.title('Orbits Of Two Planets About the Sun')
    plt.plot(x12Position, y12Position, label='Orbit of the first planet')   # plot as a continuous line
    plt.plot(x13Position, y13Position, label='Orbit of the second planet')  # plot as a continuous line
else:
    if(sunPlot):
        plt.title('Orbits of the Sun and Two Planets About the Center-of-Mass')
    else:
        plt.title('Orbits Of Two Planets About the Center-of-Mass')
    plt.plot(x2RK4, y2RK4, label='Orbit of the first planet')   # plot as a continuous line
    plt.plot(x3RK4, y3RK4, label='Orbit of the second planet')  # plot as a continuous line
    if(sunPlot):
        plt.plot(x1RK4, y1RK4, label='Orbit of the Sun')   # plot as a continuous line

plt.xlabel('x Coordinate (AU)')                             # add axis labels
plt.ylabel('y Coordinate (AU)')

#
# Center of mass point symbol
#
plt.scatter(0, 0, c='k', s=20)

plt.grid(True)
plt.xlim(minimumPositionX, maximumPositionX)
plt.ylim(minimumPositionY, maximumPositionY)

circularOrbit = False
ellipticalOrbit = True
if(stdRadialPosition2/meanRadialPosition2 < 1.e-05 and stdTotalVelocity2/meanTotalVelocity2):
    print "\n The orbit of the first planet is found to be circular (eccentricity = 0)"
    ellipticalOrbit = False
    circularOrbit = True
else:
    print "\n The orbit of the first planet is found to be elliptical (eccentricity > 0)"

eccentricity = 0
orbitTime = 0
numberOfOrbits = 0
if(circularOrbit):
    orbitTime = 2*mp.pi*meanRadialPosition2/meanTotalVelocity2
    numberOfOrbits = int(maximumTime/orbitTime)

lookingForNextPerihelion = True
lookingForNextAphelion = False
if(ellipticalOrbit):
    #
    # Algorithm is to accumulate the set of perihelion and aphelion radial values and their time values
    # The times to pass through successive perihelion points are stored from which a mean and standard deviation are computed
    # The first step is to confirm that the initial point is consistent with being a perihelion
    #
    if(radialPosition12[1] < radialPosition12[0]):
        print " Check that the initial position is consistent with being an aphelion is passed"
    else:
        lookingForNextAphelion = True
        lookingForNextPerihelion = False
        print " Looking for first aphelion position"


    radiusPerihelion = []
    radiusAphelion = []
    lastRadialPosition = radialPosition12[0]
    nTimeStep = 1
    timePerihelion = []
    timeAphelion = []
    nPerihelion = 0
    nAphelion = 0
    while nTimeStep < nTimeSteps:
        #
        # As the while loop iterates through the positions at each time step, it will be either looking for the next aphelion or the next perihelion
        #
        newRadialPosition = radialPosition12[nTimeStep]

        if(lookingForNextAphelion):
            #
            # If the new radial position is smaller than the previous radial position, then the next aphelion point has been crossed
            #
            if(newRadialPosition < lastRadialPosition):
                nAphelion += 1
                timeValue = 0.5*(timeGrid[nTimeStep] + timeGrid[nTimeStep-1])  # take an average of the current and previous times
                timeAphelion.append(timeValue)     # store the time value of this aphelion
                if(nAphelion == 1):
                     if(verbose): print "\n Found first aphelion at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                else:
                    if(verbose): print "\n Found next aphelion ", nAphelion, " at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                radiusAphelion.append(0.5*(lastRadialPosition+newRadialPosition))   # take an average of the current and previous positions
                lookingForNextPerihelion = True
                lookingForNextAphelion = False
                lastRadialPosition = newRadialPosition

            lastRadialPosition = newRadialPosition

        if(lookingForNextPerihelion):
            #
            # If the new radial position is greater than the previous radial position, then the next perihelion point has been crossed
            #
            if(newRadialPosition > lastRadialPosition):
                nPerihelion += 1
                timeValue = 0.5*(timeGrid[nTimeStep] + timeGrid[nTimeStep-1])  # take an average of the current and previous times
                timePerihelion.append(timeValue)     # store the time difference last perihelion
                if(nPerihelion == 1):
                    if(verbose): print "\n Found first Perihelion at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                else:
                    if(verbose): print "\n Found next perihelion ", nPerihelion, " at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                radiusPerihelion.append(0.5*(lastRadialPosition+newRadialPosition))
                lookingForNextPerihelion = False
                lookingForNextAphelion = True

            lastRadialPosition = newRadialPosition

        nTimeStep += 1

    timeDifferencePerihelion = []
    timeDifferenceAphelion = []
    timeDifferenceBoth = []
    nTimesPerihelion = 1
    while nTimesPerihelion < len(timePerihelion):
        timeDifference = timePerihelion[nTimesPerihelion] - timePerihelion[nTimesPerihelion - 1]
        timeDifferencePerihelion.append(timeDifference)
        timeDifferenceBoth.append(timeDifference)
        nTimesPerihelion += 1

    nTimesAphelion = 1
    while nTimesAphelion < len(timeAphelion):
        timeDifference = timeAphelion[nTimesAphelion] - timeAphelion[nTimesAphelion - 1]
        timeDifferenceAphelion.append(timeDifference)
        timeDifferenceBoth.append(timeDifference)
        nTimesAphelion += 1

    if(verbose): print "\n Number of perihelion time differences ", len(timeDifferencePerihelion)
    if(verbose): print "\n Number of aphelion time differences ", len(timeDifferenceAphelion)

    orbitTimeAphelion = 0.0
    if(len(timeDifferenceAphelion) > 0):
        orbitTimeAphelion = np.mean(timeDifferenceAphelion)
        orbitTime = orbitTimeAphelion

    orbitTimePerihelion = 0.0
    if(len(timeDifferencePerihelion) > 0):
        orbitTimePerihelion = np.mean(timeDifferencePerihelion)
        orbitTime = orbitTimePerihelion

    orbitTimeBoth = 0.0
    if(len(timeDifferenceAphelion) > 0 and len(timeDifferencePerihelion) > 0):
        orbitTimeBoth = np.mean(timeDifferenceBoth)
        orbitTime = orbitTimeBoth
        print "\n Perihelion time ", orbitTimePerihelion, "  Aphelion time ", orbitTimeAphelion, " Both time ", orbitTimeBoth

    if(nPerihelion > 0 and nAphelion > 0):
        majorAxisLength = np.mean(radiusAphelion) + np.mean(radiusPerihelion)
        eccentricity = (np.mean(radiusAphelion) - np.mean(radiusPerihelion))/majorAxisLength
        print '%s %5.3f' % (' The orbit eccentricity = ', eccentricity)

    numberOfOrbits = max(nTimesPerihelion, nTimesAphelion)

print '%s %d' % (' Number of completed orbits = ', numberOfOrbits)
if(numberOfOrbits > 0):
    print '%s %5.3f %s' % (' Orbit time = ', orbitTime, ' years')

xTextPosition1 = minimumPositionX + 0.53*(maximumPositionX - minimumPositionX)
yTextPosition1 = minimumPositionY + 0.10*(maximumPositionY - minimumPositionY)
mass2_text = ax.text(xTextPosition1, yTextPosition1, '',color='blue')
mass2_text.set_text('Mass2 = %.3e solar masses' % mass2)
x20_text = ax.text(xTextPosition1, minimumPositionY + 0.06*(maximumPositionY - minimumPositionY), '',color='blue')
x20_text.set_text('Initial x2 = %.3f AU' % x20)
vy20_text = ax.text(xTextPosition1, minimumPositionY + 0.02*(maximumPositionY - minimumPositionY), '',color='blue')
vy20_text.set_text('Initial vy2 = %.3f AU/Year' % vy20)

xTextPosition2 = minimumPositionX + 0.03*(maximumPositionX - minimumPositionX)
yTextPosition2 = yTextPosition1
mass3_text = ax.text(xTextPosition2, yTextPosition2, '',color='green')
mass3_text.set_text('Mass3 = %.3e solar masses' % mass3)
x30_text = ax.text(xTextPosition2, minimumPositionY + 0.06*(maximumPositionY - minimumPositionY), '',color='green')
x30_text.set_text('Initial x3 = %.3f AU' % x30)
vy30_text = ax.text(xTextPosition2, minimumPositionY + 0.02*(maximumPositionY - minimumPositionY), '',color='green')
vy30_text.set_text('Initial vy3 = %.3f AU/Year' % vy30)

xTextPosition3 = minimumPositionX + 0.35*(maximumPositionX - minimumPositionX)
yTextPosition3 = minimumPositionY + 0.30*(maximumPositionY - minimumPositionY)
xTextPosition4 = xTextPosition3
yTextPosition4 = minimumPositionY + 0.25*(maximumPositionY - minimumPositionY)

orbitTime_text = ax.text(xTextPosition3, yTextPosition3, '', color='red')
orbitTime_text.set_text('Orbit time = %.2f years' % orbitTime)
eccentricity_text = ax.text(xTextPosition4, yTextPosition4, '', color='red')
eccentricity_text.set_text('Eccentricity = %.3f' % eccentricity)

xTextPosition5 = 0.5*(xTextPosition1 + xTextPosition2)
yTextPosition5 = minimumPositionY + 0.65*(maximumPositionY - minimumPositionY)
maximumTime_text = ax.text(xTextPosition5, yTextPosition5, '',color='red')
maximumTime_text.set_text('Maximum time = %.1f years' % maximumTime)

xTextPosition6 = 0.5*(xTextPosition1 + xTextPosition2)
yTextPosition6 = minimumPositionY + 0.60*(maximumPositionY - minimumPositionY)
timeStep_text = ax.text(xTextPosition6, yTextPosition6, '',color='red')
timeStep_text.set_text('Time step = %.3e years' % timeStep)

npts = nTimeStep - 1
finalSpeed1 = np.sqrt(vx1RK4[npts]*vx1RK4[npts] + vy1RK4[npts]*vy1RK4[npts])
finalSpeed2 = np.sqrt(vx2RK4[npts]*vx2RK4[npts] + vy2RK4[npts]*vy2RK4[npts])
finalSpeed3 = np.sqrt(vx3RK4[npts]*vx3RK4[npts] + vy3RK4[npts]*vy3RK4[npts])

finalSeparation12 = np.sqrt((x1RK4[npts] - x2RK4[npts])*(x1RK4[npts] - x2RK4[npts]) + (y1RK4[npts] - y2RK4[npts])*(y1RK4[npts] - y2RK4[npts]))
initialSeparation12 = np.sqrt((x10 - x20)*(x10 - x20) + (y10 - y20)*(y10 - y20))
finalSeparation13 = np.sqrt((x1RK4[npts] - x3RK4[npts])*(x1RK4[npts] - x3RK4[npts]) + (y1RK4[npts] - y3RK4[npts])*(y1RK4[npts] - y3RK4[npts]))
initialSeparation13 = np.sqrt((x10 - x30)*(x10 - x30) + (y10 - y30)*(y10 - y30))
finalSeparation23 = np.sqrt((x2RK4[npts] - x2RK4[npts])*(x2RK4[npts] - x3RK4[npts]) + (y2RK4[npts] - y3RK4[npts])*(y2RK4[npts] - y3RK4[npts]))
initialSeparation23 = np.sqrt((x20 - x30)*(x20 - x30) + (y20 - y30)*(y20 - y30))

initialPE123 = -GSUNMASS*SUNMASS*mass2*mass1/(initialSeparation12*AU) - GSUNMASS*SUNMASS*mass1*mass3/(initialSeparation13*AU) - GSUNMASS*SUNMASS*mass2*mass3/(initialSeparation23*AU)
finalPE123 = -GSUNMASS*SUNMASS*mass2*mass1/(finalSeparation12*AU) - GSUNMASS*SUNMASS*mass1*mass3/(finalSeparation13*AU) - GSUNMASS*SUNMASS*mass2*mass3/(finalSeparation23*AU)
initialEnergy123 = initialPE123 + 0.5*SUNMASS*mass1*pow(vy10*AU/YEAR, 2) + 0.5*SUNMASS*mass2*pow(vy20*AU/YEAR, 2) + 0.5*SUNMASS*mass3*pow(vy30*AU/YEAR, 2)
finalEnergy123 = finalPE123 + 0.5*SUNMASS*mass1*pow(finalSpeed1*AU/YEAR, 2) + 0.5*SUNMASS*mass2*pow(finalSpeed2*AU/YEAR, 2) + 0.5*SUNMASS*mass3*pow(finalSpeed3*AU/YEAR, 2)

print "\n Initial speed for mass1 = ", vy10, " AU/year,  final speed = ", finalSpeed1
print " Initial speed for mass2 = ", vy20, " AU/year,  final speed = ", finalSpeed2
print " Initial speed for mass3 = ", vy30, " AU/year,  final speed = ", finalSpeed3
print " Initial system energy (1+2+3) = ", initialEnergy123, " (Joules)"
print " Final system energy (1+2+3) = ", finalEnergy123, " (Joules)"
print " Fractional energy change = ", (finalEnergy123 - initialEnergy123)/initialEnergy123

#
# Add the animation
#
line1, = ax.plot([], [], 'o-', lw=4)                                # this is an empty line object used in the initialization of the animation
line2, = ax.plot([], [], 'o-', lw=4)                                # this is an empty line object used in the initialization of the animation
line3, = ax.plot([], [], 'o-', lw=4)                                # this is an empty line object used in the initialization of the animation
timeValue_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='red')   # this is a placeholder for where the time will be updated
mass2Position_text = ax.text(0.02, 0.90, '', transform=ax.transAxes, color='blue')   # this is a placeholder for where the mass1 (x,y) will be updated
mass3Position_text = ax.text(0.02, 0.85, '', transform=ax.transAxes, color='green')   # this is a placeholder for where the mass3 (x,y) will be updated


def init():                                # this initialization function tells what is being animiated, i.e. a single point "line", and the x and y coordinate displays
    """initialize animation"""
    line2.set_data([], [])
    line3.set_data([], [])
    timeValue_text.set_text('')
    mass2Position_text.set_text('')
    mass3Position_text.set_text('')
    if(sunPlot and orbitSun == False):
        line1.set_data([], [])
        return line1, line2, line3, timeValue_text, mass2Position_text, mass3Position_text
    else:
        return line2, line3, timeValue_text, mass2Position_text, mass3Position_text

nSkipPoints = int(nTimeSteps/nFrames)
if(nSkipPoints < 1):
    nSkipPoints = 1
print "\n nSkipPoints = ", nSkipPoints

def animate(i):                            # this is the function this being animated, the i index is increased according the range of the frames value below
    """perform animation step"""           # i ranges from 0 to (nFrames - 1) as given in the FuncAnimation function
    global x2RK4, y2RK4, x3RK4, y3RK4, nSkipPoints, orbitTime, orbitSun

    ii = nSkipPoints*i
    if(ii >= nTimeSteps - 1):
        ii = nTimeSteps - 1
    if(orbitSun):
        x2 = x2RK4[ii] - x1RK4[ii]
        y2 = y2RK4[ii] - y1RK4[ii]
    else:
        x2 = x2RK4[ii]
        y2 = y2RK4[ii]
    v2 = np.sqrt(vx2RK4[ii]*vx2RK4[ii] + vy2RK4[ii]*vy2RK4[ii])
    timeValue = ii*timeStep
    orbitNumber = int(timeValue/orbitTime) + 1
    line2.set_data(x2,y2)  # the "line" contains a single data point (x,y) for the first planet
    timeValue_text.set_text('Time = %.2f Years, First Planet Orbit Number %d' % (timeValue,orbitNumber))
    mass2Position_text.set_text('(x2,y2) = (%.2f %.2f) AU, v2 = %.2f AU/Year' % (x2,y2,v2))
    if(orbitSun):
        x3 = x3RK4[ii] - x1RK4[ii]
        y3 = y3RK4[ii] - y1RK4[ii]
    else:
        x3 = x3RK4[ii]
        y3 = y3RK4[ii]
    v3 = np.sqrt(vx3RK4[ii]*vx3RK4[ii] + vy3RK4[ii]*vy3RK4[ii])
    line3.set_data(x3,y3)                     # the "line" contains a single data point (x,y) for the second planet
    mass3Position_text.set_text('(x3,y3) = (%.2f %.2f) AU, v3 = %.2f AU/Year' % (x3,y3,v3))
    if(sunPlot and orbitSun == False):
        x1 = x1RK4[ii]
        y1 = y1RK4[ii]
        line1.set_data(x1,y1)  # the "line" contains a single data point (x,y) for the Sun orbiting about the c.m.
        return line1, line2, line3, timeValue_text, mass2Position_text, mass3Position_text
    else:
        return line2, line3, timeValue_text, mass2Position_text, mass3Position_text

#
# The blit=True means that the figure should be updated only for the parts which have changed, i.e. the data point position and the x and y value displays
#
if(noAnimation):
    print "\n Animation step is being skipped"
else:
    print "\n Animation step is being done"
    ani = animation.FuncAnimation(fig, animate, frames=nFrames,
                                  interval=nInterval, blit=True, init_func=init)

#
# The ffmpeg binary must be put in the executable PATH, where the binary can be obtained from http://www.ffmpegmac.net/
#
if(doMovie):
    print "\n Movie production step is being done"
    FFwriter = animation.FFMpegWriter()
    ani.save('sunTwoPlanetsV1.mp4', writer = FFwriter, fps=30, extra_args=['-vcodec', 'libx264'])
    print " Movie production step is completed"
else:
    print "\n Movie production step is being skipped"

plt.show()          # show the complete figure





