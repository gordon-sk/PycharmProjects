#
# Program 2.7 Solution for the position and velocity of a projectile moving with quadratic air resistance (projectileQuadraticResistance_Chapter2V4.py)
#
#              Give the command  python projectileQuadraticResistance_Chapter2V4.py -h  to get help command on input parameters
#
#              This program functions like the projectileQuadraticResistance_Chapter2V3.py except that B2/m value is appropriate for a baseball and is velocity dependent via Equation 2.26 in the textbook
#              The inital altitude of the launch can be above sea level
#
#              This program solves a set of four simultaneous ("coupled") differential equations
#              The coupled equations are for the motion of a projectile starting with an initial velocity vector and experiencing quadratic air resistance
#                   1) dx/dt = vx    velocity component in the horizontal direction
#                   2) dy/dt = vy    velocity component in the vertical direction
#                   3) dvx/dt = ax = -[B2(y)/m]v*vx = acceleration component in the horizontal direction from the horizontal component of the resistance, with B2/m depending on altitude
#                   4) dvy/dt = ay = -g(y) - [B2(y)/m]v*vy = acceleration component in the vertical direction from gravity and the vertical component of the resistance with B2/m depending on altitude
#
#              The parameters in the differential equation are
#                   1) B2(y)/m which is the ratio of the quadratic resistance strength term B2 divided by the mass; the parameter changes as the air density changes with altitude
#                   2) g(y) is the gravity acceleration at an altitude y above the Earth's surface
#                   3) v is the magnitude of the velocity as calculated from its components: v = sqrt(vx*vx + vy*vy)
#
#              There is no way to decouple these four differential equations.  They must be solved simultaneously.
#              The ODE library is used to solve these coupled for the four equations of motion: x(t), vx(t), y(t), and v(t)
#              By convention the starting point of the motion is (x,y) = (0,0), and the initial velocity vector is given as a magnitude and an angle above the horizontal
#              Calculations are made for a list of three input angles: theta1, theta2, and theta3
#              The output figure is divided into top and bottom halves.  The top half has the calculations with resistance; the bottom half are without resistance
#
#              Using input parameters for v0 (m/s), theta0, sea level B2/m, timeStep (s), maximumTime, gAcceleration
#              Defaults are y0 = 0, v0 = 49.17 m/s, theta1 = 35 deg; B2/m = 4.0e-05 N/(m/s)^2, timeStep = 1 s, maximumTime = 20
#
#

import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt     # get matplotlib plot functions
import  argparse                    # argument parser library
import numpy as np                  # numerical functions library used by python
import math as mp                   # used for the exponential function
from scipy.integrate import odeint  # import only this single method for solving differential equations

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()
parser.add_argument('--y0', default=4300.73, type=float, help="Initial altitude; default is 0 m")
parser.add_argument('--v0', default=49.17, type=float, help="Magnitude of the initial velocity vector > 0 in m/s; default is 49.17 m/s")
parser.add_argument('--theta1', default=35.0, type=float, help="First angle of initial projectile velocity, default is 35 degrees")
parser.add_argument('--maxT', default=20.0, type=float, help="Maximum time range in s; default is 20 s")
parser.add_argument('--deltaT', default=0.1, type=float, help="Iteration time step in s; default is 0.1 s")
parser.add_argument('--varyGravity', action='store_true', help="State this parameter in order to the g value with altitude, default is False")
parser.add_argument('--B2OverM', default=4.0e-5, type=float, help="Quadratic resistance factor B2/m at sea level; default is 4.0e-5 N/(m/s)^2")
parser.add_argument('--useBaseball', action='store_true', help="Use the baseball equation for B2/m, default is False")
parser.add_argument('--varyRho', action='store_true', help="State this parameter in order to vary the air density rho, default is False")
parser.add_argument('--useAdiabatic', action='store_true', help="State this parameter in order to use the isothermal air density equation, default is False")
parser.add_argument('--useIsothermal', action='store_true', help="State this parameter in order to use the adiabatic air density equation, default is False")
args = parser.parse_args()

#
# Assign the variables in the program to the variable option names in the argument list
#

useBaseball = args.useBaseball

varyGravity = args.varyGravity

varyRho = args.varyRho
# Testing Case: (Student using IDE instead of command line)
# useBaseball, varyGravity, varyRho = True, False, True

if(varyRho):
    useIsothermal = args.useIsothermal
    useAdiabatic = args.useAdiabatic
    #useIsothermal = True
    #useAdiabatic = True
    if(useIsothermal and useAdiabatic):
        print " You cannot state both the useIsothermal and the useAdiabatic parameters; you must state only one"

    if(useIsothermal == False and useAdiabatic == False):
        print " You cannot state varyRho and then not choose either on of useIsothermal or useAdiabatic"
        exit()

v0 = args.v0
if(v0 <= 0.0):
    print "\n Cannot have the initial speed", v0, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(0)

y0 = args.y0
if(y0 < 0.0):
    print "\n Cannot have the initial altitude", y0, " be less than zero"
    print "\n Program is exiting\n"
    exit(0)

theta1 = args.theta1
maximumTime = args.maxT
timeStep = args.deltaT
B2OverMSeaLevel = args.B2OverM

resistanceCase = False
if(B2OverMSeaLevel > 0):
    resistanceCase = True
if(B2OverMSeaLevel < 0):
    print "\n Cannot use a negative B2/M value ", B2OverMSeaLevel, " which implies resistance aiding the acceleration"
    exit()

#
# Function which calculates the B2/m parameter as a function of altitude above sea level
# If the baseball version of the B2/m equation is used, then the B2/m is also velocity dependent
#
def B2OverMCalculation(altitude, vTotal):              # placeholder line for code which needs to be changed if a velocity dependent B2/m is requested
    rescaleFactor = 1.0
    if(varyRho and altitude >= y0):
        if(useIsothermal):
            rescaleFactor = mp.exp(-altitude*1.0e-04)
        if(useAdiabatic):
            altitudeFactor = 1.0 - 2.218e-5*altitude
            if(altitudeFactor > 0):                             # check that there is not a negative value
                rescaleFactor = pow(altitudeFactor, 2.5)
            else:
                rescaleFactor = 0.0
    if(useBaseball):
        B2overMvdep = .0039 + (.0058)/(1+mp.exp((vTotal-35)/5))
        return B2overMvdep*rescaleFactor   # placeholder line for code which needs to be changed if a velocity dependent B2/m is requested
    else:
        return B2OverMSeaLevel*rescaleFactor

earthRadius = 6.371e06   # meters
earthMass = 5.9723e24    # kg
universalG = 6.673e-11   # N-msquare/kgsquare
GearthMass = universalG*earthMass
distanceToEarthCenter = earthRadius + y0  # distance to the center of the Earth for the launch altitude
gAccelFixed = GearthMass/(distanceToEarthCenter*distanceToEarthCenter)

#
# Function which calculates the acceleration due to the Earth's gravity as a function of the distance from the center of the Earth
#
def gAltitudeCalculation(altitude):
    if(varyGravity):
        distanceToEarthCenter = earthRadius + altitude
        gAltitude = GearthMass/(distanceToEarthCenter*distanceToEarthCenter)
    else:
        gAltitude = gAccelFixed

    return gAltitude

theta1Radians = theta1*mp.pi/180.0                # convert to radians
vx10 = v0*mp.cos(theta1Radians)
vy10 = v0*mp.sin(theta1Radians)

# dx/dt = vx    velocity component in the horizontal direction
# dy/dt = vy    velocity component in the vertical direction
# dvx/dt = ax = -(B2/m)v*vx = acceleration component in the horizontal direction from the horizontal component of the resistance
# dvy/dt = ay = -g(y) - (B2/m)v*vy = acceleration component in the vertical direction from gravity and the vertical component of the resistance

def fDerivative(variableList, t):              # variableList dummy list array since there is more than one differential equation
    vx1 = variableList[0]                      # speed in the x direction, in this case the x direction is horizontal
    vy1 = variableList[1]                      # speed in the y direction, in this case the y direction is vertical
    y1 = variableList[3]                       # altitude above sea level
    v1Total = mp.sqrt(vx1*vx1 + vy1*vy1)       # velocity magnitude calculated from the two components
    B2OverM = B2OverMCalculation(y1, v1Total)  # placeholder line to get the B2/m parameter, which could be velocity dependent it the baseball version is used
    
    dvx1dt = -B2OverM*v1Total*vx1              # the time derivative of velocity in the x direction according to the air resistance is the acceleration component ax
    #
    # calculate the altitude dependent gravity acceleration
    #
    y1 = variableList[3]                       # altitude above sea leel
    dvy1dt = -gAltitudeCalculation(y1) - B2OverM*v1Total*vy1  # the time derivative of position in the y direction is the velocity component ay
    
    return [dvx1dt, dvy1dt, vx1, vy1]           # return the set of four derivatives as a list object

print "\n Projectile motion with quadratic air resistance"
print "  Initial velocity magnitude ", v0, " m/s, at initial angles ", theta1,  " degrees"
print "  The launch altitude is ", y0, " meters above sea level, which has a gravity acceleration ", gAccelFixed, " m/s^2"
if(varyGravity):
    print "  The gravity value is is changed as the altitude changes, according to Newton's Universal Gravity inverse square law"
else:
    print "  The gravity value is kept fixed at the input value above for all altitudes"


if(varyRho):
    if(useAdiabatic):
        print "  The B2/m parameter is varied with altitude according to the adiabatic model air density"
    else:
        print "  The B2/m parameter is varied with altitude according to the isothermal model air density"
else:
    print "  There is no variation of B2/m with altitude"
if useBaseball:
    print "  The quadratic air resistance parameter B2/m at the launch altitude is", B2OverMCalculation(y0, v0), " N/(m/s)^2"
else:
    print "  The quadratic air resistance parameter B2/m at the launch altitude is ", B2OverMCalculation(y0, v0) , " N/(m/s)^2"   # placeholder line for code that needs to be changed if B2/m is velocity dependent

print "  Time step = ", timeStep, " s,  maximum time range = ", maximumTime, " s"

nTimeSteps = int(maximumTime/timeStep)
# obtain the differential equatioin solutions using the odeint method from the ScyPy library
timeGrid = np.linspace(0, maximumTime,nTimeSteps)              # time grid used by odeint method

initialValuesSpeedsPositions = [vx10, vy10, 0.0, y0]    # sets of starting values of vx, vy, x0, y0 for the iteration using theta1
fourSolution = odeint(fDerivative, initialValuesSpeedsPositions, timeGrid)   # odeint returns a list of values which are the numerical vx, vy, x, and y functions for theta1

vx1Numerical = fourSolution[:,0]            # vx function of time obtained with odeint solution for the altitude dependent gravity
vy1Numerical = fourSolution[:,1]            # vy function of time obtained with odeint solution for the altitude dependent gravity
x1Numerical = fourSolution[:,2]             # x function of time obtained with odeint solution for the altitude dependent gravity
y1Numerical = fourSolution[:,3]             # y function of time obtained with odeint solution for the altitude dependent gravity

maximumHeight1 = max(y1Numerical)

nTimeStep = 0
# Do the iteration over time steps from 0 to the maximum time requested to get the horizontal range and maximum height values

horizontalRange1 = -1.0
xLastPosition1 = 0.0
yLastPosition1 = 0.0

#
# It is of interest only to plot the trajectory as long as the height is >= 0
# Some of the later y points may extend below y = 0, and these are not plotted
#
nPlotPoints1 = 0
xPlot1, yPlot1 = [], []   # lists to hold the (x,y) points to be plotted for the first angle

#
# The while loop determines the maximum horizontal range for each of the numerical solutions
#
while nTimeStep < nTimeSteps: # loop over the time range
    if(nTimeStep > 0):
        yPosition = y1Numerical[nTimeStep]
        xPosition = x1Numerical[nTimeStep]
        
        if(yPosition <= y0):
            #
            # Interpolate for the horizontal range based on these two data points which were positive and then negative for yPosition
            #
            localSlope = (yPosition - yLastPosition1)/(xPosition-xLastPosition1)
            localIntercept = yPosition - localSlope*xPosition
            horizontalRange1 = (y0 - localIntercept)/localSlope
            #
            # assume contant speed in the horizontal direction to get the time for reaching y = y0
            #
            timeHorizontalRange1 = timeGrid[nTimeStep-1] + ((horizontalRange1 - xLastPosition1)/(xPosition - xLastPosition1))*timeStep
            break
        else:
            xPlot1.append(xPosition)
            yPlot1.append(yPosition)
            xLastPosition1 = xPosition
            yLastPosition1 = yPosition
            nPlotPoints1 += 1
    
    nTimeStep = nTimeStep + 1       # go to the next time

print "\n The maximum height is ", maximumHeight1, " m"

foundHorizontalRange = False
if(horizontalRange1 < 0.0):
    print "\n Horizontal range was not determined.  The time duration may have been too short"
    xMaximum = max(x1Numerical)
    xMaximumString = str(float(int(10*xMaximum))/10.0)
    label1String = 'Solution for theta = ' + str(theta1) + ' deg, maximum x = ' + xMaximumString + ' m'
else:
    foundHorizontalRange = True
    xMaximum = horizontalRange1
    horizontalRange1String = str(float(int(10*horizontalRange1))/10.0)
    label1String = 'Solution for theta = ' + str(theta1) + ' deg, range = ' + horizontalRange1String + ' m'
    print "  The horizontal range is ", horizontalRange1, " m, with a travel time of ", timeHorizontalRange1, " s"
    print "  the horizontal range in feet is ", horizontalRange1*(1.0/.3048), ' feet'
    timeRangeString = str(float(int(10*(timeHorizontalRange1+0.05)))/10.0)   # one decimal place precision
    timeTravelString = 'Time to travel range distance = ' + timeRangeString + ' s'

plt.figure(1)

plt.plot(xPlot1, yPlot1, 'bo', label=label1String)    # blue dots for the first solution plot

plt.xlabel('Horizontal position (m)')                             # add axis labels
plt.ylabel('Vertical position above sea level (m)')

# compose string variables about the time parameters for use in putting text on the plot
v0String = 'Initial velocity v0 = ' + str(v0) + ' m/s'
B2OverMString = 'Sea level B2/m = ' + str(B2OverMSeaLevel) + ' per meter'
timeStepString = 'Time step = ' + str(timeStep) + ' s'
launchAltitudeString = 'Launch altitude = ' + str(y0) + ' m above sea level'
if(varyRho):
    if(useAdiabatic):
        rhoString = 'B2/m parameter is rescaled with adiabatic model'
    if(useIsothermal):
        rhoString = 'B2/m parameter is rescaled with isothermal model'
else:
    rhoString = 'Air density changes do not affect B2/m value'

if(varyGravity):
    plt.text(0.22*xMaximum, 0.46*(maximumHeight1-y0)+y0, 'Calculation with g depending on altitude')
else:
    plt.text(0.22*xMaximum, 0.46*(maximumHeight1-y0)+y0, 'Calculation with a constant g at all altitudes')
plt.text(0.35*xMaximum, 0.40*(maximumHeight1-y0)+y0, v0String)
if(useBaseball == False):
    plt.text(0.35*xMaximum, 0.34*(maximumHeight1-y0)+y0, B2OverMString)
else:
    plt.text(0.35*xMaximum, 0.34*(maximumHeight1-y0)+y0, 'Baseball B2/m used')
plt.text(0.35*xMaximum, 0.28*(maximumHeight1-y0)+y0, timeStepString)
plt.text(0.22*xMaximum, 0.22*(maximumHeight1-y0)+y0, launchAltitudeString)
if(foundHorizontalRange):
    plt.text(0.22*xMaximum, 0.14*(maximumHeight1-y0)+y0, timeTravelString)
plt.text(0.22*xMaximum, 0.08*(maximumHeight1-y0)+y0, rhoString)
plt.grid(True)
plt.legend(loc=0)

if(resistanceCase):
    plt.title('Projectile Motion With Quadratic Air Resistance')
else:
    plt.title('Projectile Motion NO Air Resistance')

plt.xlim(0, int(1.1*xMaximum))
plt.ylim(y0, int(1.3*(maximumHeight1-y0)+y0))

plt.show()          # show the complete figure with the upper and lower subplots
