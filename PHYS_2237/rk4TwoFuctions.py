#
# Program 1.8: Solution to two coupled differential equations which have analytic solutions (rk4TwoFunction_Chapter1V1.py)
#
#              Give the command  python rk4TwoFunction_Chapter1V1.py -h  to get help command on input parameter
#
#              The RK4 algorithm for solving two differential equations is contained in the odeint method
#              which is imported from the scipy library  (odeint = ordinary differential equation integration)
#
#              The two coupled differential equations with x as the independent variable with y and z as the dependent variables are
#                   dy/dx = x*x + z  and  dz/dx = x + y   with initial conditions as y(0) = 0 and x(0) = 0
#              The analytic solutions produced by Mathematica are
#                   y(x) = 3[sinh(x) - x]  and  z(x) = 3[cosh(x) - 1) - x*x
#
#              Using input parameters for range in x and the size of the iteration step in x
#              Defaults are xMaximum = 5, x step size = 0.1
#
#
import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt     # get matplotlib plot functions
import  argparse                    # argument parser library, to handle command line parameters
import numpy as np                  # numerical functions library used by python
import math                         # used for the exponential function
from scipy.integrate import odeint  # import only this single method for solving differential equations

parser = argparse.ArgumentParser()   # obtain the object called parser which can interpret command line parameters
parser.add_argument('--xMaximum', default=5.0, type=float, help='Maximum range in x axis, default = 5.0')   # declare a parameter and its default value
parser.add_argument('--xStep', default=0.1, type=float, help='Iteration step size for x, default = 0.1') # declare a parameter and its default value
args = parser.parse_args()  # retrieve the parameters into another object called arg

#
# Retrieve the two parameters for the decay calculation from the command line
#
xMaximum = args.xMaximum                   # maximum range along the x axis
xStep = args.xStep                         # iteration step size in x

print "\n Program to solve numerically two coupled differential equations and compare to the analytic results"
print "  Maximum range along the x axis = ", xMaximum
print "  Iteration step size in x = ", xStep

# define the derivatives
def fDerivative(dependentVariables, x):                #  List array for the derivatives; x is the single independent variable
    y = dependentVariables[0]
    z = dependentVariables[1]
    dydx = x*x + z
    dzdx = x + y
    return [dydx, dzdx]                                # return the two derivatives

initialConditions = [0,0]   # these are the initial conditions y(0) = and z(0) = 0; one could easily set up input parameters for the initial conditions
nXSteps = int(xMaximum/xStep)

# obtain the RK4 solutions for y(x) and z(x) using the odeint method from the ScyPy library
xGrid = np.linspace(0, xMaximum, nXSteps)                 # x grid used by odeint method
twoSolutions = odeint(fDerivative, initialConditions, xGrid)     # odeint returns a list of two arrays which are the y(x) and z(x) solutions
yRK4 = twoSolutions[:,0]            # y function of x obtained with the odeint library
zRK4 = twoSolutions[:,1]            # z function of x obtained with the odeint library

#
# Now calculate the analytic solutions obtained from Mathematica
#
# Set up the list arrays for the iteration which obtains the analytic solution
yAnalytic = []               # analytic array for storing y(x)
zAnalytic = []               # analytic array for storing z(x)

maximumFractionalErrorY = -1.0
maximumFractionalErrorZ = -1.0
xStepNumber = 0
# Do the iteration over x steps
while xStepNumber<nXSteps:               # loop over the x range
    x = xGrid[xStepNumber]
    yTrue = 3*(math.sinh(x) - x)         # true analytic solution number for y(x)
    yAnalytic.append(yTrue)              # analytic solution list for y
    if(yTrue != 0.0):
        fractionalErrorY = (yRK4[xStepNumber] - yTrue)/yTrue
        if(fractionalErrorY > maximumFractionalErrorY):
            maximumFractionalErrorY = fractionalErrorY

    zTrue = 3*(math.cosh(x) - 1) - x*x   # true analytic solution number for z(x)
    zAnalytic.append(zTrue)              # analytic solution list for z
    if(zTrue != 0.0):
        fractionalErrorZ = (zRK4[xStepNumber] - zTrue)/zTrue
        if(fractionalErrorZ > maximumFractionalErrorZ):
            maximumFractionalErrorZ = fractionalErrorZ

    xStepNumber = xStepNumber + 1        # go to the next x value

# The iteration loop has concluded to produce the two analytic solution lists of values
print "\n The maximum fractional error in y(x) is ", maximumFractionalErrorY
print " The maximum fractional error in z(x) is ", maximumFractionalErrorZ
print "\n"

# code to set up the plots for y(x) and z(x) in a single figure
plt.figure(1)       # start a figure
plt.subplot(211)    # this sets the upper half plot
plt.plot(xGrid, yRK4, 'ro', label='RK4 Numerical Solution For y(x)')    # red dots for the RK4 y(x) solution plot
plt.plot(xGrid, yAnalytic, label='Analytic Solution For y(x)')          # blue line for analytic y(x) solution plot

# compose string variables about the iteration step parameter for use in putting text on the plot
xStepString = 'Iteration step along x = ' + str(xStep)
yMaximum = max(yAnalytic)
plt.text(0.15*xMaximum, 0.35*yMaximum, xStepString)  # text to document the time parameters used

plt.xlabel('Independent variable x')                             # add axis labels
plt.ylabel('Dependent variable y')
plt.ylim(0.0, int(1.25*yMaximum))
plt.xlim(0.0, int(1.25*xMaximum))

plt.title('Coupled differential equations: Analytic and RK4 Solutions')
plt.grid(True)
plt.legend(loc=0)

plt.subplot(212)    # this sets the lower half plot
plt.plot(xGrid, zRK4, 'ro', label='RK4 Numerical Solution For z(x)')    # red dots for the RK4 z(x) solution plot
plt.plot(xGrid, zAnalytic, label='Analytic Solution For z(x)')          # blue line for analytic z(x) solution plot
plt.xlabel('Independent variable x')                             # add axis labels
plt.ylabel('Dependent variable z')
zMaximum = max(zAnalytic)
plt.ylim(0.0, int(1.25*zMaximum))
plt.xlim(0.0, int(1.25*xMaximum))
plt.grid(True)
plt.legend(loc=0)

plt.show()          # show the complete figure with the upper and lower subplots