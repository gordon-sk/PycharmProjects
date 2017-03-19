#
# Program 2.4 Solution for the position and velocity of a projectile moving with quadratic air resistance (projectileQuadraticResistance_Chapter2V1.py)
#
#              Give the command  python projectileQuadraticResistance_Chapter2V1.py -h  to get help command on input parameters
#
#              This program functions like the bikeHillQuadraticResistance_Chapter2V2.py except the motion is now in two dimensions: horizontal x and vertical y
#
#              This program solves a set of four simultaneous ("coupled") differential equations
#              The coupled equations are for the motion of a projectile starting with an initial velocity vector and experiencing quadratic air resistance
#                   1) dx/dt = vx    velocity component in the horizontal direction
#                   2) dy/dt = vy    velocity component in the verticl direction
#                   3) dvx/dt = ax = -(B2/m)v*vx = acceleration component in the horizontal direction from the horizontal component of the resistance
#                   4) dvy/dt = ay = -g - (B2/m)v*vy = acceleration component in the vertical direction from gravity and the vertical component of the resistance
#
#              The parameters in the differential equation are
#                   1) B2/m which is the ratio of the quadratic resistance strength term B2 divided by the mass
#                   2) g is the gravity acceleration
#                   3) v is the magnitude of the velocity as calculated from its components: v = sqrt(vx*vx + vy*vy)
#
#              There is no way to decouple these four differential equations.  They must be solved simultaneously.
#              The ODE library is used to solve these coupled for the four equations of motion: x(t), vx(t), y(t), and v(t)
#              By convention the starting point of the motion is (x,y) = (0,0), and the initial velocity vector is given as a magnitude and an angle above the horizontal
#              Calculations are made for a list of three input angles: theta1, theta2, and theta3
#              The output figure is divided into top and bottom halves.  The top half has the calculations with resistance; the bottom half are without resistance
#
#              Using input parameters for v0 (m/s), theta0, B2/m, timeStep (s), maximumTime, gAcceleration
#              Defaults are v0 = 700 m/s, theta1, theta2, theta3 = 35, 45, 55 deg; B2/m = 4.0e-05 N/(m/s)^2, timeStep = 2 s, maximumTime = 100, gAcceleration = 9.8 m/s^2
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
parser.add_argument('--v0', default=60.0, type=float, help="Magnitude of the initial velocity vector > 0 in m/s; default is 700 m/s")
parser.add_argument('--theta1', default=10.0, type=float, help="First angle of initial projectile velocity, default is 35 degrees")
parser.add_argument('--theta2', default=45.0, type=float, help="Second angle of initial projectile velocity, default is 45 degrees")
parser.add_argument('--theta3', default=55.0, type=float, help="Third angle of initial projectile velocity, default is 55 degrees")
parser.add_argument('--maxT', default=100.0, type=float, help="Maximum time range in s; default is 150 s")
parser.add_argument('--deltaT', default=.10, type=float, help="Iteration time step in s; default is 1.0 s")
parser.add_argument('--B2OverM', default=.002, type=float, help="Quadratic resistance factor B2/m; default is 4.0e-5 N/(m/s)^2")
parser.add_argument('--gAccel', default=9.81855392611, type=float, help="Constant acceleration from gravity; default is 9.819 m/s^2")

args = parser.parse_args()
#
# Assign the variables in the program to the variable option names in the argument list
#
v0 = args.v0
if(v0 <= 0.0):
    print "\n Cannot have the initial speed", v0, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(0)

theta1 = args.theta1
theta2 = args.theta2
theta3 = args.theta3
maximumTime = args.maxT
timeStep = args.deltaT
B2OverM = args.B2OverM
gAccel = args.gAccel

resistanceCase = False
if(B2OverM > 0):
    resistanceCase = True
if(B2OverM < 0):
    print "\n Cannot use a negative B2/M value ", B2OverM, " which implies resistance aiding the acceleration"
    exit()

theta1Radians = theta1*mp.pi/180.0                # convert to radians
theta2Radians = theta2*mp.pi/180.0                # convert to radians
theta3Radians = theta3*mp.pi/180.0                # convert to radians
vx10 = v0*mp.cos(theta1Radians)
vy10 = v0*mp.sin(theta1Radians)
vx20 = v0*mp.cos(theta2Radians)
vy20 = v0*mp.sin(theta2Radians)
vx30 = v0*mp.cos(theta3Radians)
vy30 = v0*mp.sin(theta3Radians)

#
# All the trajectory plots will have a common maximum x value, taken from the analytic solution
#
analyticRange1 = v0*v0*mp.sin(2*theta1Radians)/gAccel
analyticRange2 = v0*v0*mp.sin(2*theta2Radians)/gAccel
analyticRange3 = v0*v0*mp.sin(2*theta3Radians)/gAccel
analyticRangeMaximum = max(analyticRange1, analyticRange2, analyticRange3)

#
# All the trajectory plots will have a common maximum y value, taken from the analytic solution
#
height1Fact = analyticRange1/(2*v0*mp.cos(theta1Radians))
analyticHeight1 = mp.tan(theta1Radians)*analyticRange1/2 - (gAccel/2)*height1Fact*height1Fact

height2Fact = analyticRange2/(2*v0*mp.cos(theta2Radians))
analyticHeight2 = mp.tan(theta2Radians)*analyticRange2/2 - (gAccel/2)*height2Fact*height2Fact

height3Fact = analyticRange3/(2*v0*mp.cos(theta3Radians))
analyticHeight3 = mp.tan(theta3Radians)*analyticRange3/2 - (gAccel/2)*height3Fact*height3Fact

analyticHeightMaximum = max(analyticHeight1, analyticHeight2, analyticHeight3)


# dx/dt = vx    velocity component in the horizontal direction
# dy/dt = vy    velocity component in the vertical direction
# dvx/dt = ax = -(B2/m)v*vx = acceleration component in the horizontal direction from the horizontal component of the resistance
# dvy/dt = ay = -g - (B2/m)v*vy = acceleration component in the vertical direction from gravity and the vertical component of the resistance

def fDerivative(variableList, t):           # variableList dummy list array since there is more than one differential equation

    vx = variableList[0]                    # speed in the x direction, in this case the x direction is horizontal
    vy = variableList[1]                    # speed in the y direction, in this case the y direction is vertical
    vTotal = mp.sqrt(vx*vx + vy*vy)         # velocity magnitude calculated from the two components
    dvxdt = -B2OverM*vTotal*vx              # the time derivative of velocity in the x direction according to the air resistance, acceleration component ax
    dvydt = -gAccel - B2OverM*vTotal*vy     # the time derivative of postion in the y direction is the acceleration component ay

    return [dvxdt, dvydt, vx, vy]           # return the four derivatives as a list object containing four elements

print "\n Projectile motion with quadratic air resistance"
print "  Initial velocity magnitude ", v0, " m/s, at initial angles ", theta1, ", ", theta2, ", and ", theta3, " degrees"
print "  The quadratic air resistance parameter is ", B2OverM, " N/(m/s)^2"
print "  Time step = ", timeStep, " s,  maximum time range = ", maximumTime, " s"

nTimeSteps = int(maximumTime/timeStep)
# obtain the differential equatioin solutions using the odeint method from the ScyPy library
timeGrid = np.linspace(0, maximumTime,nTimeSteps)           # time grid used by odeint method

initialValuesSpeedsPositions = [vx10, vy10, 0.0, 0.0]         # starting values of vx, vy, x0, y0 for the iteration using theta1
fourSolution1 = odeint(fDerivative, initialValuesSpeedsPositions, timeGrid)   # odeint returns a list of values which are the numerical vx, vy, x, and y functions for the first angle
vx1Numerical = fourSolution1[:,0]            # vx function of time obtained with odeint solution
vy1Numerical = fourSolution1[:,1]            # vy function of time obtained with odeint solution
x1Numerical = fourSolution1[:,2]             # x function of time obtained with odeint solution
y1Numerical = fourSolution1[:,3]             # y function of time obtained with odeint solution

for x in fourSolution1:
    print x

initialValuesSpeedsPositions = [vx20, vy20, 0.0, 0.0]         # starting values of vx, vy, x0, y0 for the iteration using theta2
fourSolution2 = odeint(fDerivative, initialValuesSpeedsPositions, timeGrid)   # odeint returns a list of values which are the numerical vx, vy, x, and y functions for the second angle
vx2Numerical = fourSolution2[:,0]            # vx function of time obtained with odeint solution
vy2Numerical = fourSolution2[:,1]            # vy function of time obtained with odeint solution
x2Numerical = fourSolution2[:,2]             # x function of time obtained with odeint solution
y2Numerical = fourSolution2[:,3]             # y function of time obtained with odeint solution

initialValuesSpeedsPositions = [vx30, vy30, 0.0, 0.0]         # starting values of vx, vy, x0, y0 for the iteration using theta3
fourSolution3 = odeint(fDerivative, initialValuesSpeedsPositions, timeGrid)   # odeint returns a list of values which are the numerical vx, vy, x, and y functions for the third angle
vx3Numerical = fourSolution3[:,0]            # vx function of time obtained with odeint solution
vy3Numerical = fourSolution3[:,1]            # vy function of time obtained with odeint solution
x3Numerical = fourSolution3[:,2]             # x function of time obtained with odeint solution
y3Numerical = fourSolution3[:,3]             # y function of time obtained with odeint solution

nTimeStep = 0
# Do the iteration over time steps from 0 to the maximum time requested to get the horizontal range and maximum height values

maximumHeight1 = 0.0
horizontalRange1 = 0.0
timeMaximumHeight1 = 0.0
timeHorizontalRange1 = 0.0
xLastPosition1 = 0.0
yLastPosition1 = 0.0

maximumHeight2 = 0.0
horizontalRange2 = 0.0
timeMaximumHeight2 = 0.0
timeHorizontalRange2 = 0.0
xLastPosition2 = 0.0
yLastPosition2 = 0.0

maximumHeight3 = 0.0
horizontalRange3 = 0.0
timeMaximumHeight3 = 0.0
timeHorizontalRange3 = 0.0
xLastPosition3 = 0.0
yLastPosition3 = 0.0

#
# It is of interest only to plot the trajectory as long as the height is >= 0
# Some of the later y points may extend below y = 0, and these are not plotted
#
finalPointFound1 = False
finalPointFound2 = False
finalPointFound3 = False
nPlotPoints1 = 0
nPlotPoints2 = 0
nPlotPoints3 = 0
xPlot1, yPlot1 = [], []   # lists to hold the (x,y) points to be plotted for the first angle
xPlot2, yPlot2 = [], []   # lists to hold the (x,y) points to be plotted for the second angle
xPlot3, yPlot3 = [], []   # lists to hold the (x,y) points to be plotted for the third angle

xAnalytic1, yAnalytic1 = [], []
xAnalytic2, yAnalytic2 = [], []
xAnalytic3, yAnalytic3 = [], []
#
# The while loop determines the maximum horizontal range for each of the numerical solutions
#
while nTimeStep < nTimeSteps: # loop over the time range
    if(nTimeStep > 0):
        yPosition = y1Numerical[nTimeStep]
        xPosition = x1Numerical[nTimeStep]
        if(maximumHeight1 < yPosition):
            maximumHeight1 = yPosition
            timeMaximumHeight1 = timeGrid[nTimeStep]

        if(yPosition <= 0.0 and finalPointFound1 == False):
            #
            # Interpolate for the horizontal range based on these two data points which were positive and then negative for yPosition
            #
            localSlope = (yPosition - yLastPosition1)/(xPosition-xLastPosition1)
            localIntercept = yPosition - localSlope*xPosition
            horizontalRange1 = -localIntercept/localSlope
            #
            # assume contant speed in the horizontal direction to get the time for reaching y = 0
            #
            timeHorizontalRange1 = timeGrid[nTimeStep-1] + ((horizontalRange1 - xLastPosition1)/(xPosition - xLastPosition1))*timeStep
            finalPointFound1 = True
        if(finalPointFound1 == False):
            xPlot1.append(xPosition)
            yPlot1.append(yPosition)
            xLastPosition1 = xPosition
            yLastPosition1 = yPosition
            nPlotPoints1 += 1

        yPosition = y2Numerical[nTimeStep]
        xPosition = x2Numerical[nTimeStep]
        if(maximumHeight2 < yPosition):
            maximumHeight2 = yPosition
            timeMaximumHeight2 = timeGrid[nTimeStep]

        if(yPosition <= 0.0 and finalPointFound2 == False):
            #
            # Interpolate for the horizontal range based on these two data points which were positive and then negative for yPosition
            #
            localSlope = (yPosition - yLastPosition2)/(xPosition-xLastPosition2)
            localIntercept = yPosition - localSlope*xPosition
            horizontalRange2 = -localIntercept/localSlope
            #
            # assume contant speed in the horizontal direction to get the time for reaching y = 0
            #
            timeHorizontalRange2 = timeGrid[nTimeStep-1] + ((horizontalRange2 - xLastPosition2)/(xPosition - xLastPosition2))*timeStep
            finalPointFound2 = True
        if(finalPointFound2 == False):
            xPlot2.append(xPosition)
            yPlot2.append(yPosition)
            xLastPosition2 = xPosition
            yLastPosition2 = yPosition
            nPlotPoints2 += 1

        yPosition = y3Numerical[nTimeStep]
        xPosition = x3Numerical[nTimeStep]
        if(maximumHeight3 < yPosition):
            maximumHeight3 = yPosition
            timeMaximumHeight3 = timeGrid[nTimeStep]

        if(yPosition <= 0.0 and finalPointFound3 == False):
            #
            # Interpolate for the horizontal range based on these two data points which were positive and then negative for yPosition
            #
            localSlope = (yPosition - yLastPosition3)/(xPosition-xLastPosition3)
            localIntercept = yPosition - localSlope*xPosition
            horizontalRange3 = -localIntercept/localSlope
            #
            # assume contant speed in the horizontal direction to get the time for reaching y = 0
            #
            timeHorizontalRange3 = timeGrid[nTimeStep-1] + ((horizontalRange3 - xLastPosition3)/(xPosition - xLastPosition3))*timeStep
            finalPointFound3 = True
        if(finalPointFound3 == False):
            xPlot3.append(xPosition)
            yPlot3.append(yPosition)
            xLastPosition3 = xPosition
            yLastPosition3 = yPosition
            nPlotPoints3 += 1

    if(finalPointFound1 and finalPointFound2 and finalPointFound3):
        break;

    nTimeStep = nTimeStep + 1       # go to the next time

xAnalytic1, yAnalytic1 = [], []
xAnalytic2, yAnalytic2 = [], []
xAnalytic3, yAnalytic3 = [], []
#
# The while loop stores the analytic solutions for the same three angles
#
nTimeStep = 0
x1Factor = v0*mp.cos(theta1Radians)
x2Factor = v0*mp.cos(theta2Radians)
x3Factor = v0*mp.cos(theta3Radians)
y1Factor = v0*mp.sin(theta1Radians)
y2Factor = v0*mp.sin(theta2Radians)
y3Factor = v0*mp.sin(theta3Radians)
gFactor = gAccel/2
while nTimeStep < nTimeSteps: # loop over the time range
    if(nTimeStep > 0):
        t = timeGrid[nTimeStep]
        yPosition = y1Factor*t - gFactor*t*t
        if(yPosition >= 0):
            xAnalytic1.append(x1Factor*t)
            yAnalytic1.append(yPosition)

        yPosition = y2Factor*t - gFactor*t*t
        if(yPosition >= 0):
            xAnalytic2.append(x2Factor*t)
            yAnalytic2.append(yPosition)

        yPosition = y3Factor*t - gFactor*t*t
        if(yPosition >= 0):
            xAnalytic3.append(x3Factor*t)
            yAnalytic3.append(yPosition)

    nTimeStep += 1

plt.figure(1)
horizontalRange1String = str(float(int(10*horizontalRange1))/10.0)
label1String = 'Solution for theta = ' + str(theta1) + ' deg, range = ' + horizontalRange1String + ' m'
plt.plot(xPlot1, yPlot1, 'bo', label=label1String)    # blue dots for the first solution plot
plt.plot(xAnalytic1, yAnalytic1, 'b')    # blue line for the first analytic plot

horizontalRange2String = str(float(int(10*horizontalRange2))/10.0)
label2String = 'Solution for theta = ' + str(theta2) + ' deg, range = ' + horizontalRange2String + ' m'
plt.plot(xPlot2, yPlot2, 'ro', label=label2String)    # red dots for the second solution plot
plt.plot(xAnalytic2, yAnalytic2, 'r')    # red line for the second analytic plot

horizontalRange3String = str(float(int(10*horizontalRange3))/10.0)
label3String = 'Solution for theta = ' + str(theta3) + ' deg, range = ' + horizontalRange3String + ' m'
plt.plot(xPlot3, yPlot3, 'go', label=label3String)    # green dots for the second solution plot
plt.plot(xAnalytic3, yAnalytic3, 'g')    # green line for the third analytic plot

plt.xlabel('Horizontal position (m)')                             # add axis labels
plt.ylabel('Vertical position (m)')

# compose string variables about the time parameters for use in putting text on the plot
v0String = 'Initial velocity v0 = ' + str(v0) + ' m/s'
B2OverMString = 'B2/m = ' + str(B2OverM) + ' per meter'
timeStepString = 'Time step = ' + str(timeStep) + ' s'

plt.text(0.65*analyticRangeMaximum, 0.97*analyticHeightMaximum, v0String)
plt.text(0.65*analyticRangeMaximum, 0.91*analyticHeightMaximum, B2OverMString)
plt.text(0.75*analyticRangeMaximum, 0.85*analyticHeightMaximum, timeStepString)

plt.grid(True)
plt.legend(loc=0)

if(resistanceCase):
    plt.title('Projectile Motion With Quadratic Air Resistance')
else:
    plt.title('Projectile Motion NO Air Resistance')

plt.xlim(0, int(1.1*analyticRangeMaximum))
plt.ylim(0, int(1.3*analyticHeightMaximum))

plt.show()          # show the complete figure with the upper and lower subplots