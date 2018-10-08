#
# Program 2.1 Solution for the position and velocity of a bicycle traveling horizontally with quadratic air resistance (bikeQuadraticResistance_Chapter2V1.py)
#
#              Give the command  python bikeQuadraticResistance_Chapter2V1.py -h  to get help command on input parameters
#
#              This program solves a set of two simultaneous ("coupled") differential equations
#              The coupled equations are for the horizontal motion of a bicycle starting a initial non-zero speed v0 and experiencing quadratic air resistance
#                   dx/dt = v  and dv/dt = P/mv - C*rho*A*v*v/2m
#                   1) P is the contant power of the bicyclist
#                   2) C is a scaling constant for the strength of the air resistance force
#                   3) rho is the air density
#                   4) m is the mass of the bicycle and its rider
#                   5) A is the effective area of the bicycle and rider
#
#              There is a choice between using the odeint library for the RK4 algorithm solution, or using the Euler algorithm as used by the textbook
#              The odeint library is used by default
#              The coupled differential equations are solved to obtain v(t) and x(t)
#
#              Using input parameters for v0 (m/s), maximumTime (s), timeStep (s), power (watts), cResistance, rho (kg/m^3), mass (kg), area (m^2)
#              Defaults are v0 = 4 m/s, maximumTime = 400 s, timeStep = 2 s, power = 400 watts, cResistance = 1 , rho = 1.2 kg/m^3, mass = 70 kg, area = 0.33 m^2
#
# import standard libraries ("don't reinvent the wheel")
#
import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt     # get matplotlib plot functions
import sys                          # used to get the number of command line arguments
import argparse                    # argument parser library
import numpy as np                  # numerical functions library used by python
import math as mp                   # used for the exponential function
from scipy.integrate import odeint  # import only this single method for solving differential equations

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()
parser.add_argument('--v0', default=4.0, type=float, help="Initial speed with v0 > 0 in m/s, default is 4 m/s")
parser.add_argument('--maxT', default=200.0, type=float, help="Maximum time range in s, default is 200 s")
parser.add_argument('--deltaT', default=2.0, type=float, help="Iteration time step in s, default is 2 s")
parser.add_argument('--power', default=400.0, type=float, help="Constant power level by rider, default is 400 watts")
parser.add_argument('--cFactor', default=1.0, type=float, help="Quadratic resistance scale factor, default is 1.0")
parser.add_argument('--rho', default=1.2, type=float, help="Air density in kg/m^3, default is 1.2 kg/cubic-meter")
parser.add_argument('--mass', default=70.0, type=float, help="Mass of rider + bicycle in kg, default is 70 kg")
parser.add_argument('--area', default=0.33, type=float, help="Effective cross sectional area of rider+bicycle in square meters, default is 0.33 square meters")
parser.add_argument('--algorithmChoice', default='RK4', help="Choice between using the odeint library RK4 algorithm or the Euler algorithm, default is 'RK4' ")
args = parser.parse_args()

#
# Get the input parameters from the command line
#
numberOfArguments = len(sys.argv)
if(numberOfArguments == 1):
    print "\n All the default paramter choices are used"  # there is always at least one argument and that first one is the name of the python script itself

algorithmChoice =args.algorithmChoice
if(algorithmChoice != 'RK4' and algorithmChoice != 'Euler'):
    print "\n The algorithm choice ", algorithmChoice, " is not recognized as either RK4 nor Euler"
    exit(1)
#
# Assign the variables in the program to the variable option names in the argument list
#
v0 = args.v0
if(v0 <= 0.0):
    print "\n Cannot have the initial speed", v0, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(1)

maximumTime = args.maxT
timeStep = args.deltaT
power = args.power
cResistance = args.cFactor
rho = args.rho
mass = args.mass
area = args.area

resistanceFactor = cResistance*rho*area/(2*mass)

#
# define the time derivative functions dv/dt = P/mv - C*rho*A*v*v/2m  and  dx/dt = v
# this function is used only if the RK4 algorithm choice is made
#
def fDerivative(variableList, t):                  # variableList dummy list array since there is more than one differential equation
    v = variableList[0]                            # speed in the x direction
    dvdt = power/(mass*v) - resistanceFactor*v*v   # the time derivative of velocity in the x direction according to the power and the air resistance
    dxdt = v                                       # the time derivative of postion in the x direction is the velocity
    return [dvdt, dxdt]                            # return the two derivatives as a list object containing two elements

terminalSpeed = 0
resistanceCase = True       # put in check if the resistance does not exist because of parameter choices, to prevent a divison by zero
if(rho <= 0.0 or cResistance <= 0 or area <= 0):
    resistanceCase = False   # there is no terminal speed with no resistance
else:
    terminalSpeed = np.power(2.0*power/(rho*cResistance*area), 1./3.) # analytic equation obtained when dv/dt = 0

print "\n Horizontal motion with quadratic air resistance for mass = ", mass, " kg,  v0 = ", v0, " m/s"
print "  Power = ", power, " watts,  resistance scale factor = ", cResistance
print "  air density = ", rho, " kg/cubic-meter,  cross sectional area = ", area, " square meters"
print "  time step = ", timeStep, " s,  maximum time range = ", maximumTime, " s"
if(resistanceCase):
    print "  Predicted terminal speed ", terminalSpeed, " m/s"
if(algorithmChoice == 'RK4'):
    print "  The RK4 algorthim from the odeint library will be used"
else:
    print "  The Euler algorithm will be used"
    powerOverMassFactor = (power/mass)*timeStep
    cResistanceRhoAreaOverTwoMassFactor = (cResistance*rho*area/(2*mass))*timeStep

print "  "

nTimeSteps = int(maximumTime/timeStep)
timeGrid = np.linspace(0, maximumTime,nTimeSteps)           # time grid used for the iteration steps
if(algorithmChoice == 'RK4'):
    labelString = 'Numerical solution, RK4 algorithm'
    # obtain the differential equation solutions using the odeint method from the ScyPy library
    initialValuesSpeedPosition = [v0, 0.0]                      # starting values of velocity and position for the iteration
    twoSolution = odeint(fDerivative, initialValuesSpeedPosition, timeGrid)   # odeint returns a list of values which is the NRK4 solution
    vSolution = twoSolution[:,0]            # velocity function of time obtained with RK4 solution
    xSolution = twoSolution[:,1]            # position function of time obtained with RK4 solution
else:
    labelString = 'Numerical solution, Euler algorithm'
    # obtain the differential equation solutions using the Euler algorithm
    vSolution = []          # list to hold the numerical solution for the speed
    xSolution = []          # list to hold the numerical solution for the position
    vSolution.append(v0)    # starting value of velocity for the iteration
    xSolution.append(0.0)   # starting value of position for the iteration
    vLast = v0
    xLast = 0
    nTimeStep = 1
    while nTimeStep < nTimeSteps: # loop over the time range
        vNext = vLast + powerOverMassFactor/vLast - cResistanceRhoAreaOverTwoMassFactor*vLast*vLast
        xNext = xLast + vNext*timeStep
        vSolution.append(vNext)
        xSolution.append(xNext)

        vLast = vNext
        xLast = xNext
        nTimeStep += 1

nTimeStep = 0
# Do the iteration over time steps from 0 to the maximum time requested
maximumVelocity = -1.0e+12
minimumVelocity = +1.0e+12
maximumPosition = -1.0e+12
minimumPosition = +1.0e+12

while nTimeStep < nTimeSteps: # loop over the time range
    velocity = vSolution[nTimeStep]
    if(velocity < minimumVelocity):
        minimumVelocity = velocity
    if(velocity > maximumVelocity):
        maximumVelocity = velocity

    position = xSolution[nTimeStep]
    if(position < minimumPosition):
        minimumPosition = position
    if(position > maximumPosition):
        maximumPosition = position

    nTimeStep = nTimeStep + 1       # go to the next time

# The iteration loop has concluded to produce the limits of the plot

# compose string variables about the time parameters for use in putting text on the plot
v0String = 'Initial velocity = ' + str(v0) + ' m/s'
powerString = 'Constant power = ' + str(power) + ' watts'
cResistanceString = 'Quadratic resistance factor = ' + str(cResistance)
massString = 'Total mass = ' + str(mass) + ' kg'
rhoString = 'Air density = ' + str(rho) + ' kg/cubic-meter'
terminalSpeedString = 'Predicted terminal speed = ' + str(terminalSpeed) + ' m/s'

print "\n  Final calculated speed = ", vSolution[nTimeSteps - 1], " m/s"
print "  Final calculated position = ", xSolution[nTimeSteps - 1], " m"

# code to set up the two plots in a single figure
plt.figure(1)       # start a figure
plt.subplot(211)    # this sets the upper half plot for the v(t)
plt.plot(timeGrid, vSolution, 'ro', label=labelString)    # red dots for the numerical solution plot

xTextPosition = 0.4*maximumTime
if(resistanceCase):
    plt.text(0.3*maximumTime, 1.05*terminalSpeed, terminalSpeedString)   # text to document the parameters used
    plt.text(xTextPosition, 0.82*terminalSpeed, v0String)   # text to document the parameters used
    plt.text(xTextPosition, 0.69*terminalSpeed, powerString) # text to document the parameters used
    plt.text(xTextPosition, 0.56*terminalSpeed, cResistanceString) # text to document the parameters used
    plt.text(xTextPosition, 0.43*terminalSpeed, massString) # text to document the parameters used
    plt.text(xTextPosition, 0.30*terminalSpeed, rhoString) # text to document the parameters used
else:
    xTextPosition = 0.53*maximumTime
    plt.text(xTextPosition, 0.82*maximumVelocity, v0String)   # text to document the parameters used
    plt.text(xTextPosition, 0.69*maximumVelocity, powerString) # text to document the parameters used
    plt.text(xTextPosition, 0.56*maximumVelocity, cResistanceString) # text to document the parameters used
    plt.text(xTextPosition, 0.43*maximumVelocity, massString) # text to document the parameters used
    plt.text(xTextPosition, 0.30*maximumVelocity, rhoString) # text to document the parameters used
    vAnalytic = []
    time = 0.0
    v0Square = v0*v0
    twicePowerOverMassFactor = 2*power/mass
    nTimeStep = 0
    while nTimeStep < nTimeSteps:
        #
        # Compute the analytic solution at each time value
        #
        vComputed = np.sqrt(v0Square + twicePowerOverMassFactor*time)
        vAnalytic.append(vComputed)
        time += timeStep
        nTimeStep += 1

    plt.plot(timeGrid, vAnalytic, label='No Air Resistance Analytic Solution')    # continuous line for the analytic solution plot

plt.xlabel('Time (s)')                             # add axis labels
plt.ylabel('Velocity (m/s)')

plt.title('Horizontal Motion With Quadratic Air Resistance')
plt.grid(True)
plt.ylim(0.0, int(1.2*maximumVelocity))
plt.legend(loc=0)

plt.subplot(212)    # this sets the lower half plot for the x(t)
plt.plot(timeGrid, xSolution, 'bo', label=labelString)    # blue dots for the numerical solution plot
plt.grid(True)
plt.ylim(0.0, int(1.2*maximumPosition))
plt.legend(loc=0)
plt.xlabel('Time (s)')                             # add axis labels
plt.ylabel('Position (m)')

plt.show()          # show the complete figure with the upper and lower subplots