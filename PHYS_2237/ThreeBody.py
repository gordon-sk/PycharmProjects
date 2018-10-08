#
# Program 4.12 Solution for the orbits of a planet in a binary star system orbiting about their center-of-mass (binaryStarPlanetSystem_Chapter4V1.py)
#
#              Give the command  python binaryStarPlanetSystem_Chapter4V1.py -h  to get help command on input parameters
#
#              This program is a modification of the binaryStarSystem_Chapter4V2.py program, adding a planet to orbit the binary star pair
#
#              The input arguments to this program are the following, assuming that the coordinate system origin is the cm point and that there is zero total momentum
#                   1) initial x position of the first star (in AU, > 0) with respect to the cm point
#                   2) mass of the first star in Solar mass units
#                   3) initial speed of the first star (in units of AU/year, > 0)   (the three mass values and the linear momentum of the first star and the planet determine the velocity of the second star)
#                   4) mass of the second star in Solar mass units  (the three mass values along with the position of the first star and the planet determine the initial position of the second star)
#                   5) initial x position of the planet (in AU, > 0) with respect to the cm point
#                   6) mass of the planet in Solar mass units
#                   7) initial speed of the planet (in units of AU/year, > 0)
#                   8) maximum time (years)
#                   9) time step (years)
#
#              The coupled differential equations are for the motions of the two stars in a Cartesian coordinate system whose (x,y) origin is the cm point
#                   The first star is assumed to start at the perihelion of its orbit
#                   The total linear momentum of the system is zero, meaning that the cm point does not move
#                   There are twelve differential equations to be solved for: 1) (vx1, vy1), 2) (vx2, vy2), 3) (vx3,vy3), 4) (x1, y1), 6) (x2, y2), 6) (x3, y3)
#                   1) dvx1/dt = ax1 = -4*pi*pi*mass2*(x2-x1)/d12**3/2 - 4*pi*pi*mass3*(x3-x1)/d13**3/2   where d12 is the separation distance between the two stars and d13 is the separation between the first star and the planet
#                   2) dvy1/dt = ay1 = -4*pi*pi*mass2*(y2-y1)/d12**3/2 - 4*pi*pi*mass3*(y3-y1)/d13**3/2   where d12 is the separation distance between the two stars and d13 is the separation between the first star and the planet
#                   3) dvx2/dt = ax2 = -4*pi*pi*mass1*(x1-x2)/d12**3/2 - 4*pi*pi*mass3*(x3-x2)/d23**3/2   where d12 is the separation distance between the two stars and d23 is the separation between the second star and the planet
#                   4) dvy2/dt = ay2 = -4*pi*pi*mass1*(y1-y2)/d12**3/2 - 4*pi*pi*mass3*(y3-y2)/d23**3/2   where d12 is the separation distance between the two stars and d23 is the separation between the second star and the planet
#                   5) dvx3/dt = ax3 = -4*pi*pi*mass1*(x1-x3)/d13**3/2 - 4*pi*pi*mass2*(x2-x3)/d23**3/2   where d13 is the separation distance between the first star and the planet, and d23 is the separation between the second star and the planet
#                   6) dvy3/dt = ay3 = -4*pi*pi*mass1*(y1-y3)/d13**3/2 - 4*pi*pi*mass2*(y2-y3)/d23**3/2   where d12 is the separation distance between the two stars and d23 is the separation between the second star and the planet
#                   7) dx1/dt = vx1    first star velocity component in the x direction
#                   8) dy1/dt = vy1    first star velocity component in the y direction
#                   9) dx2/dt = vx2    second star velocity component in the x direction
#                   10) dy2/dt = vy2   second star velocity component in the y direction
#                   11) dx3/dt = vx3   planet velocity component in the x direction
#                   12) dy3/dt = vy3   planet velocity component in the y direction
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
parser.add_argument('--x10', default=0.15, type=float, help="Initial aphelion x10 in AU; default 0.15")
parser.add_argument('--vy10', default=3.0, type=float, help="Initial speed vy10 in AU/Y; default 3.0")
parser.add_argument('--mass1', default=0.6, type=float, help="Mass of first star in Solar mass units; default 0.60")
parser.add_argument('--mass2', default=0.4, type=float, help="Mass of second star in Solar mass units; default 0.40")
parser.add_argument('--x30', default=1.0, type=float, help="Initial planet aphelion x30 in AU; default 1.0")
parser.add_argument('--vy30', default=6.28318, type=float, help="Initial planet speed vy30 in AU/Y; default 6.28318")
parser.add_argument('--mass3', default=3.0e-7, type=float, help="Mass of planet in Solar mass units; default 3.0e-07")
parser.add_argument('--maxT', default=3.0, type=float, help="Maximum time range in years; default 3")
parser.add_argument('--deltaT', default=0.002, type=float, help="Iteration time step in years; default 0.002")
parser.add_argument('--verbose', action='store_true', help="Give extra printout about orbit pattern recognition results; default False")
parser.add_argument('--nFrames', default=500, type=int, help="Number of frames in the video; default 500")
parser.add_argument('--nInterval', default=200, type=int, help="Time difference between video frames in milliseconds; default 200")
parser.add_argument('--noAnimation', action='store_true', help="Do not do any animation call; default False")
parser.add_argument('--doMovie', action='store_true', help="Do the movie output, requires ffmpeg library; default False")
parser.add_argument('--planetPlot', action='store_true', help="Plot the planet orbits (x,y) plane; default False")
parser.add_argument('--radialPlot', action='store_true', help="Plot the planet orbits as r(t); default False")

args = parser.parse_args()

exam = True
if exam:
    pass
#
# Get the input parameters from the command line
#
numberOfArguments = len(sys.argv)
if(numberOfArguments == 1):
    print "\n All the default paramter choices are used"  # there is always at least one argument and that first one is the name of the python script itself

#
# Assign the variables in the program to the variable option names in the argument list
#
verbose = args.verbose

x10 = args.x10
if(x10 < 0.0):
    print "\n Cannot have the initial position", x10, " be less than zero"
    print "\n Program is exiting\n"
    exit(0)

vy10 = args.vy10
if(vy10 < 0.0):
    print "\n Cannot have the initial velocity component", vy10, " be less than zero"
    print "\n Program is exiting\n"
    exit(0)

mass1 = args.mass1
mass2 = args.mass2

x30 = args.x30
if(x30 <= 0.0):
    print "\n Cannot have the initial position", x30, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(0)

vy30 = args.vy30
if(vy30 <= 0.0):
    print "\n Cannot have the initial velocity omponent", vy30, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(0)

mass3 = args.mass3

maximumTime = args.maxT
timeStep = args.deltaT
nFrames = args.nFrames
nInterval = args.nInterval
noAnimation = args.noAnimation # retrieve the choice to not do the animation
doMovie = args.doMovie  # retrieve the choice to make the mp4 movie
planetPlot = args.planetPlot  # retrieve the choice to plot the planet orbiting (x,y) plane
radialPlot = args.radialPlot  # retrieve the choice to plot the planet orbiting r(t)

x20 = -x10*mass1/mass2 - x30*mass3/mass2
vy20 = -vy10*mass1/mass2 - vy30*mass3/mass2

y10 = 0
y20 = 0
y30 = 0
vx10 = 0.0
vx20 = 0.0
vx30 = 0.0

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

vEscape = np.sqrt(2.0*FOURPISQUARE*mass2*mass2/((x10-x20)*(mass2+mass1)))

print "\n Binary star plus planet system calculation"
print "   First star mass = ", mass1, " Solar masses,  second star mass = ", mass2, " Solar mass units"
print "   First star position (x10, y10) = (", x10, ", ", y10, ") in AU units"
print "   Second star position (x20, y20) = (", x20, ", ", y20, ") in AU units"
print "   First star velocity (vx10, vy10) = (", vx10, ", ", vy10, ") in AU/Y units"
print "   Second star velocity (vx20, vy20) = (", vx20, ", ", vy20, ") in AU/Y units"
print "   Planet masss = ", mass3, "Solar masses, with (x30,y30) = (",x30,",",y30,")"
print "   Planet velocity (vx30, vy30) = (", vx30, ", ", vy30, ") in AU/Y units"
print "   Maximum time = ", maximumTime, " years,  with in iteration time step = ", timeStep, " years"
print "   The escape velocity is ", vEscape, " AU/Year"
print " "

#                   1) dvx1/dt = ax1 = -4*pi*pi*mass2*(x2-x1)/d12**3/2 - 4*pi*pi*mass3*(x3-x1)/d13**3/2   where d12 is the separation distance between the two stars and d13 is the separation between the first star and the planet
#                   2) dvy1/dt = ay1 = -4*pi*pi*mass2*(y2-y1)/d12**3/2 - 4*pi*pi*mass3*(y3-y1)/d13**3/2   where d12 is the separation distance between the two stars and d13 is the separation between the first star and the planet
#                   3) dvx2/dt = ax2 = -4*pi*pi*mass1*(x1-x2)/d12**3/2 - 4*pi*pi*mass3*(x3-x2)/d23**3/2   where d12 is the separation distance between the two stars and d23 is the separation between the second star and the planet
#                   4) dvy2/dt = ay2 = -4*pi*pi*mass1*(y1-y2)/d12**3/2 - 4*pi*pi*mass3*(y3-y2)/d23**3/2   where d12 is the separation distance between the two stars and d23 is the separation between the second star and the planet
#                   5) dvx3/dt = ax3 = -4*pi*pi*mass1*(x1-x3)/d13**3/2 - 4*pi*pi*mass2*(x2-x3)/d23**3/2   where d13 is the separation distance between the first star and the planet, and d23 is the separation between the second star and the planet
#                   6) dvy3/dt = ay3 = -4*pi*pi*mass1*(y1-y3)/d13**3/2 - 4*pi*pi*mass2*(y2-y3)/d23**3/2   where d12 is the separation distance between the two stars and d23 is the separation between the second star and the planet
#                   7) dx1/dt = vx1    first star velocity component in the x direction
#                   8) dy1/dt = vy1    first star velocity component in the y direction
#                   9) dx2/dt = vx2    second star velocity component in the x direction
#                   10) dy2/dt = vy2    second star velocity component in the y direction
#                   11) dx3/dt = vx3    planet velocity component in the x direction
#                   12) dy3/dt = vy3    planet velocity component in the y direction
#
printDebug = False
def fDerivative(variableList, t):                  # variableList dummy list array since there is more than one differential equation
    global mass1, mass2, FOURPISQUARE, printDebub
    vx1 = variableList[0]                          # first star speed in the x direction
    vy1 = variableList[1]                          # first star speed in the y direction
    vx2 = variableList[2]                          # second star speed in the x direction
    vy2 = variableList[3]                          # second star speed in the y direction
    vx3 = variableList[4]                          # planet speed in the x direction
    vy3 = variableList[5]                          # planet speed in the y direction
    x1 = variableList[6]                           # first star x coordinate
    y1 = variableList[7]                           # first star y coordinate
    x2 = variableList[8]                           # second star x coordinate
    y2 = variableList[9]                           # second star y coordinate
    x3 = variableList[10]                          # planet x coordinate
    y3 = variableList[11]                          # planet y coordinate
    #
    # Compute distance parameters
    #
    d12Sq = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) # square of the distance between the two stars
    d13Sq = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) # square of the distance between the first star and the planet
    d23Sq = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) # square of the distance between the second star and the planet

    d12Cu = d12Sq*mp.sqrt(d12Sq)   # cube of the distance between the two stars
    d13Cu = d13Sq*mp.sqrt(d13Sq)   # cube of the distance between the first star and the planet
    d23Cu = d23Sq*mp.sqrt(d23Sq)   # cube of the distance between the second star and the planet

    #
    # Compute the derivatives
    #
    dvx1dt = -FOURPISQUARE*(x1-x2)*mass2/d12Cu - FOURPISQUARE*(x1-x3)*mass3/d13Cu    # the time derivative of velocity v1 in the x direction according to the universal gravity force component
    dvy1dt = -FOURPISQUARE*(y1-y2)*mass2/d12Cu - FOURPISQUARE*(y1-y3)*mass3/d13Cu    # the time derivative of velocity v1 in the y direction according to the universal gravity force component
    dvx2dt = -FOURPISQUARE*(x2-x1)*mass1/d12Cu - FOURPISQUARE*(x2-x3)*mass3/d23Cu    # the time derivative of velocity v2 in the x direction according to the universal gravity force component
    dvy2dt = -FOURPISQUARE*(y2-y1)*mass1/d12Cu - FOURPISQUARE*(y2-y3)*mass3/d23Cu    # the time derivative of velocity v2 in the y direction according to the universal gravity force component
    dvx3dt = -FOURPISQUARE*(x3-x1)*mass1/d13Cu - FOURPISQUARE*(x3-x2)*mass2/d23Cu    # the time derivative of velocity v3 in the x direction according to the universal gravity force component
    dvy3dt = -FOURPISQUARE*(y3-y1)*mass1/d13Cu - FOURPISQUARE*(y3-y2)*mass2/d23Cu    # the time derivative of velocity v3 in the y direction according to the universal gravity force component
    netForceX = mass1*dvx1dt + mass2*dvx2dt + mass3*dvx3dt
    netForceY = mass1*dvy1dt + mass2*dvy2dt + mass3*dvy3dt
    if(printDebug or abs(netForceX) > 1.0e-06 or abs(netForceY) > 1.0e-06):
        print "\n Summed action-reaction forces, x direction ", netForceX, ": F1 ", mass1*dvx1dt, ", F2 ", mass2*dvx2dt, ", F3 ", mass3*dvx3dt
        print " First pair  m2m1 ", -FOURPISQUARE*(x1-x2)*mass2/d12Cu, ", m3m1 ", -FOURPISQUARE*(x1-x3)*mass3/d13Cu
        print " Second pair m1m2 ", -FOURPISQUARE*(x2-x1)*mass1/d12Cu, ", m3m2 ", -FOURPISQUARE*(x2-x3)*mass3/d23Cu
        print " Third pair  m1m3 ", -FOURPISQUARE*(x3-x1)*mass1/d13Cu, ", m2m3 ", -FOURPISQUARE*(x3-x2)*mass2/d23Cu
        exit()

    return [dvx1dt, dvy1dt, dvx2dt, dvy2dt, dvx3dt, dvy3dt, vx1, vy1, vx2, vy2, vx3, vy3]    # return the twelve derivatives as a list object containing eight elements

nTimeSteps = int(maximumTime/timeStep)
# obtain the differential equation solutions using the odeint method from the ScyPy library
timeGrid = np.linspace(0, maximumTime,nTimeSteps)           # time grid used by odeint method
initialValuesSpeedsPositions = [vx10, vy10, vx20, vy20, vx30, vy30, x10, y10, x20, y20, x30, y30]   # starting values of vx1, vy1, vx2, vy2, x1, y1, x2, y2 for the iteration
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

#
# Extract the orbit time and the eccentricity for the planet
# First set up the radial position and the total velocity arrays
#
radialPosition3 = []
totalVelocity3= []

nTimeStep = 0
while nTimeStep < nTimeSteps:
    radialPosition3.append(mp.sqrt(x3RK4[nTimeStep]*x3RK4[nTimeStep] + y3RK4[nTimeStep]*y3RK4[nTimeStep]))
    totalVelocity3.append(mp.sqrt(vx3RK4[nTimeStep]*vx3RK4[nTimeStep] + vy3RK4[nTimeStep]*vy3RK4[nTimeStep]))
    nTimeStep += 1

meanRadialPosition3 = np.mean(radialPosition3)
stdRadialPosition3 = np.std(radialPosition3)
meanTotalVelocity3 = np.mean(totalVelocity3)
stdTotalVelocity3 = np.std(totalVelocity3)

if(radialPlot):
    plt.figure(1)       # start a figure for a single plot of the orbit
    radiusString = '%s %4.3f %s %4.3f %s' % ('Radius = ', meanRadialPosition3, ' +/- ', stdRadialPosition3, ' AU')
    plt.plot(timeGrid, radialPosition3, label=radiusString)    # plot as a continuous line
    plt.xlabel('Time (years)')
    plt.ylabel('Radius (AU)')
    rMinimum = min(radialPosition3)
    rMaximum = max(radialPosition3)
    yMinimum = 0.98*rMinimum
    yMaximum = rMaximum + 0.3*(rMaximum - rMinimum)
    plt.xlim(0, timeGrid[nTimeSteps -1])
    plt.ylim(yMinimum, yMaximum)
    plt.title('Radius of a Planet in a Binary Star System')
    plt.grid(True)
    plt.legend(loc=0)
    plt.show()          # show the complete figure
    exit()

minimumPositionX1 = 1.5*min(x1RK4)
minimumPositionX2 = 1.5*min(x2RK4)
minimumPositionX3 = 1.5*min(x3RK4)

maximumPositionX1 = 1.5*max(x1RK4)
maximumPositionX2 = 1.5*max(x2RK4)
maximumPositionX3 = 1.5*max(x3RK4)
if(planetPlot):
    minimumPositionX = min(minimumPositionX1, minimumPositionX2, minimumPositionX3)
    maximumPositionX = max(maximumPositionX1, maximumPositionX2, maximumPositionX3)
else:
    minimumPositionX = min(minimumPositionX1, minimumPositionX2)
    maximumPositionX = max(maximumPositionX1, maximumPositionX2)

minimumPositionY1 = 1.5*min(y1RK4)
minimumPositionY2 = 1.5*min(y2RK4)
minimumPositionY3 = 1.5*min(y3RK4)
maximumPositionY1 = 1.5*max(y1RK4)
maximumPositionY2 = 1.5*max(y2RK4)
maximumPositionY3 = 1.5*max(y3RK4)
if(planetPlot):
    minimumPositionY = min(minimumPositionY1, minimumPositionY2, minimumPositionY3)
    maximumPositionY = max(maximumPositionY1, maximumPositionY2, maximumPositionY3)
else:
    minimumPositionY = min(minimumPositionY1, minimumPositionY2)
    maximumPositionY = max(maximumPositionY1, maximumPositionY2)

print "\n Minimum x position for the plot = ", minimumPositionX,
print ",  maximum x position for the plot = ", maximumPositionX
print " Minimum y position for the plot = ", minimumPositionY,
print ",  maximum y position for the plot = ", maximumPositionY

fig = plt.figure(1)       # start a figure for plots of the two orbits
ax = fig.add_subplot(111)
plt.plot(x1RK4, y1RK4, label='Orbit of the first star')    # plot as a continuous line
plt.plot(x2RK4, y2RK4, label='Orbit of the second star')   # plot as a continuous line
if(planetPlot):
    plt.plot(x3RK4, y3RK4, label='Orbit of the planet')    # plot as a continuous line

plt.xlabel('x Coordinate (AU)')                             # add axis labels
plt.ylabel('y Coordinate (AU)')
plt.title('Orbit of a planet in a binary star system about the center-of-mass')

#
# Center of mass point symbol
#
plt.scatter(0, 0, c='k', s=20)
#
# Draw an arrow to the center of mass point
#
#ax.annotate('Center of Mass', xy=(0, 0), xytext=(minimumPositionX +0.60*(maximumPositionX - minimumPositionX), 0.03*(maximumPositionY - minimumPositionY)),
#            arrowprops=dict(facecolor='black', shrink=0.05))

plt.grid(True)
plt.xlim(minimumPositionX, maximumPositionX)
plt.ylim(minimumPositionY, maximumPositionY)

circularOrbit = False
ellipticalOrbit = True
if(stdRadialPosition3/meanRadialPosition3 < 1.e-05 and stdTotalVelocity3/meanTotalVelocity3):
    print "\n The orbit of the planet is found to be circular (eccentricity = 0)"
    ellipticalOrbit = False
    circularOrbit = True
else:
    print "\n The orbit of the planet is found to be elliptical (eccentricity > 0)"

eccentricity = 0
orbitTime = 0
numberOfOrbits = 0
if(circularOrbit):
    orbitTime = 2*mp.pi*meanRadialPosition3/meanTotalVelocity3
    numberOfOrbits = int(maximumTime/orbitTime)

lookingForNextPerihelion = True
lookingForNextAphelion = False
if(ellipticalOrbit):
    #
    # Algorithm is to accumulate the set of perihelion and aphelion radial values and their time values
    # The times to pass through successive perihelion points are stored from which a mean and standard deviation are computed
    # The first step is to confirm that the initial point is consistent with being a perihelion
    #
    if(radialPosition3[1] < radialPosition3[0] and totalVelocity3[1] > totalVelocity3[0]):
        print " Check that the initial position is consistent with being an aphelion is passed"
    else:
        print "\n Initial position is not an aphelion"
        print "  r0 ", radialPosition3[0], "  r1 ", radialPosition3[1]
        print " "
        exit()

    radiusPerihelion = []
    radiusAphelion = []
    lastRadialPosition = radialPosition3[0]
    nTimeStep = 1
    timePerihelion = []
    timeAphelion = []
    nPerihelion = 0
    nAphelion = 0
    while nTimeStep < nTimeSteps:
        #
        # As the while loop iterates through the positions at each time step, it will be either looking for the next aphelion or the next perihelion
        #
        newRadialPosition = radialPosition3[nTimeStep]

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


xTextPosition1 = minimumPositionX + 0.60*(maximumPositionX - minimumPositionX)
yTextPosition1 = minimumPositionY + 0.10*(maximumPositionY - minimumPositionY)
mass1_text = ax.text(xTextPosition1, yTextPosition1, '',color='blue')
mass1_text.set_text('Mass1 = %.1f solar masses' % mass1)
x10_text = ax.text(xTextPosition1, minimumPositionY + 0.06*(maximumPositionY - minimumPositionY), '',color='blue')
x10_text.set_text('Initial x1 = %.2f AU' % x10)
vy10_text = ax.text(xTextPosition1, minimumPositionY + 0.02*(maximumPositionY - minimumPositionY), '',color='blue')
vy10_text.set_text('Initial vy1 = %.2f AU/Year' % vy10)

xTextPosition2 = minimumPositionX + 0.05*(maximumPositionX - minimumPositionX)
yTextPosition2 = yTextPosition1
mass2_text = ax.text(xTextPosition2, yTextPosition2, '',color='green')
mass2_text.set_text('Mass2 = %.1f solar masses' % mass2)
x20_text = ax.text(xTextPosition2, minimumPositionY + 0.06*(maximumPositionY - minimumPositionY), '',color='green')
x20_text.set_text('Initial x2 = %.2f AU' % x20)
vy20_text = ax.text(xTextPosition2, minimumPositionY + 0.02*(maximumPositionY - minimumPositionY), '',color='green')
vy20_text.set_text('Initial vy2 = %.2f AU/Year' % vy20)

xTextPosition3 = minimumPositionX + 0.35*(maximumPositionX - minimumPositionX)
yTextPosition3 = minimumPositionY + 0.30*(maximumPositionY - minimumPositionY)
xTextPosition4 = xTextPosition3
yTextPosition4 = minimumPositionY + 0.25*(maximumPositionY - minimumPositionY)
if(planetPlot):
    orbitTime_text = ax.text(xTextPosition3, yTextPosition3, '', color='red')
    orbitTime_text.set_text('Orbit time = %.2f years' % orbitTime)
    eccentricity_text = ax.text(xTextPosition4, yTextPosition4, '', color='red')
    eccentricity_text.set_text('Eccentricity = %.3f' % eccentricity)
    yTextPosition7 = minimumPositionY + 0.44*(maximumPositionY - minimumPositionY)
    mass3_text = ax.text(xTextPosition4, yTextPosition7, '', color='red')
    mass3_text.set_text('Planet mass = %.2e Solar Mass' % mass3)
    yTextPosition8 = minimumPositionY + 0.40*(maximumPositionY - minimumPositionY)
    x30_text = ax.text(xTextPosition4, yTextPosition8, '', color='red')
    x30_text.set_text('Initial x3 = %.3f AU' % x30)
    yTextPosition9 = minimumPositionY + 0.36*(maximumPositionY - minimumPositionY)
    vy30_text = ax.text(xTextPosition4, yTextPosition9, '', color='red')
    vy30_text.set_text('Initial vy3 = %.3f AU/Year' % vy30)

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
if(planetPlot): line3, = ax.plot([], [], 'o-', lw=4)                                # this is an empty line object used in the initialization of the animation
timeValue_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='red')   # this is a placeholder for where the time will be updated
mass1Position_text = ax.text(0.02, 0.90, '', transform=ax.transAxes, color='blue')   # this is a placeholder for where the mass1 (x,y) will be updated
mass3Position_text = ax.text(0.02, 0.85, '', transform=ax.transAxes, color='green')   # this is a placeholder for where the mass3 (x,y) will be updated

def init():                                # this initialization function tells what is being animiated, i.e. a single point "line", and the x and y coordinate displays
    """initialize animation"""
    line1.set_data([], [])
    line2.set_data([], [])
    timeValue_text.set_text('')
    if(planetPlot):
        line3.set_data([], [])
        mass3Position_text.set_text('')
        return line1, line2, line3, timeValue_text, mass1Position_text, mass3Position_text
    else:
        return line1, line2, timeValue_text, mass1Position_text

nSkipPoints = nTimeSteps/nFrames
if(nSkipPoints < 1):
    nSkipPoints = 1
print "\n nSkipPoints = ", nSkipPoints

def animate(i):                            # this is the function this being animated, the i index is increased according the range of the frames value below
    """perform animation step"""           # i ranges from 0 to (nFrames - 1) as given in the FuncAnimation function
    global x1RK4, y1RK4, x3RK4, y3RK4, nSkipPoints, orbitTime

    ii = nSkipPoints*i
    if(ii >= nTimeSteps - 1):
        ii = nTimeSteps - 1
    x1 = x1RK4[ii]
    y1 = y1RK4[ii]
    v1 = np.sqrt(vx1RK4[ii]*vx1RK4[ii] + vy1RK4[ii]*vy1RK4[ii])
    timeValue = ii*timeStep
    orbitNumber = int(timeValue/orbitTime) + 1
    line1.set_data(x1,y1)  # the "line" contains a single data point (x,y) for the first star
    x2 = x2RK4[ii]
    y2 = y2RK4[ii]
    line2.set_data(x2,y2)  # the "line" contains a single data point (x,y) for the second star
    timeValue_text.set_text('Time = %.2f Years, Planet Orbit Number %d' % (timeValue,orbitNumber))
    mass1Position_text.set_text('(x1,y1) = (%.2f %.2f) AU, v1 = %.2f AU/Year' % (x1,y1,v1))
    if(planetPlot):
        x3 = x3RK4[ii]
        y3 = y3RK4[ii]
        v3 = np.sqrt(vx3RK4[ii]*vx3RK4[ii] + vy3RK4[ii]*vy3RK4[ii])
        line3.set_data(x3,y3)                     # the "line" contains a single data point (x,y) for the planet
        mass3Position_text.set_text('(x3,y3) = (%.2f %.2f) AU, v3 = %.2f AU/Year' % (x3,y3,v3))
        return line1, line2, line3, timeValue_text, mass1Position_text, mass3Position_text
    else:
        return line1,line2, timeValue_text, mass1Position_text

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
    ani.save('binaryStarPlanet.mp4', writer = FFwriter, fps=30, extra_args=['-vcodec', 'libx264'])
    print " Movie production step is completed"
else:
    print "\n Movie production step is being skipped"

plt.show()          # show the complete figure