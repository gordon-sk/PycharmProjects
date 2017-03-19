#
# Program 1.1: Radioactive Decay Differential Equation Solutions (radioactiveDecayEuler_Chapter1V1.py)
#
#              Give the command  python radioactiveDecayEuler_Chapter1V1.py -h  to get help command on input parameters
#
#              Comparing the Analytic Solution with Euler Method Solution for radiocative decay equation
#              dN/dt = -lambda*N(t) = -N(t)/tau --> Analytic solution is  N(t) = N0*exp(-t/tau)
#              Using input parameters for decay constant tau (seconds), time range factor, time-step factor, and initial number of nuclei N0
#              Defaults are tau = 1 second, time range factor = 5, time step factor = 0.05, N0 = 1,000
#
#              Exact solution: N(t) = N0*exp(-t/tau)
#              Euler iteration N(t2) = N(t1)*(1 - timeStep/tauDecay), where t2 = t1 + timeStep
#

#
# import standard libraries ("don't reinvent the wheel")
#
import matplotlib
matplotlib.use('TkAgg')         # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt # get matplotlib plot functions
import  argparse                # argument parser library, to handle command line parameters
import math                     # used for the exponential function

parser = argparse.ArgumentParser()   # obtain the object called parser which can interpret command line parameters
parser.add_argument('--tauDecay', default=1.0, type=float, help='Decay constant tauDecay in seconds, default = 1')   # declare a parameter and its default value
parser.add_argument('--timeRangeFactor', default=5.0, type=float, help='Maximum time range is timeRangeFactor*tauDecay, default = 5.0') # declare a parameter and its default value
parser.add_argument('--timeStepFactor', default=0.05, type=float, help='Iteration time step is timeStepFactor*tauDecay, default = 0.05') # declare a parameter and its default value
parser.add_argument('--initialNucleiNumber', default=1000, type=int, help='Number of nuclei present at t = 0, default = 1000') # declare a parameter and its default value
parser.add_argument('--logPlot', default='No', help='Option to use a log plot instead of a linear plot, default = No') # declare a parameter and its default value
args = parser.parse_args()  # retrieve the parameters into another object called arg

#
# Retrieve the four default parameters for the decay calculation from the command line
#
tauDecay = args.tauDecay                   # tau decay constant in seconds
timeRangeFactor = args.timeRangeFactor     # time range = tauDecay*timeRangeFactor
timeStepFactor = args.timeStepFactor       # time step = timeStepFactor*tauDecay
initialNucleiNumber = args.initialNucleiNumber  # initial number of nuclei

useLogPlot = False           # the plots will have a linear vertical scale unless a log scale is requested
if(args.logPlot == 'Yes'):
    useLogPlot = True        # the plots will have a linear vertical log scale

print "\n Program to compare the Euler algorithm and the analytic solution to solve a radioactive decay example"
print "  Decay time constant", tauDecay
print "  Time range factor", timeRangeFactor
print "  Time step factor", timeStepFactor
print "  Number of nuclei at time t = 0 ", initialNucleiNumber
if(useLogPlot):
    print "  A semi-log plot will be produced"

# Calculate time range for the plot
maximumTime = tauDecay*timeRangeFactor
# Calculate time step for the Euler iteration
timeStep = tauDecay*timeStepFactor
# Dimensionless paramter used in the Euler method
timeStepOverTau = timeStep/tauDecay # variable used in Euler method

# Set up the initial conditions for the iteration which obtains the analytic and the Euler iteration solution
N = initialNucleiNumber             # initial number of nuclei to start the Euler data values
tAnalytic, NAnalytic = [], []       # time and analytic number arrays for plotting
tEuler, NEuler = [], []             # tEuler and NEuler[] holds Euler time and solution

t = 0                               # starting time for the iteration
# Do the iteration over time steps from 0 to the maximum time requested
while t<maximumTime:                # loop over the time range
    tAnalytic.append(t)             # record time for the Analytic array
    tEuler.append(t)                # record time for Euler array
    NEuler.append(N)                # record nuclei number for Euler array
    NAnalytic.append(initialNucleiNumber*math.exp(-t/tauDecay))   # analytic solution
    N = N*(1.0 - timeStepOverTau)   # Euler's method uses first derivative at the starting step
    t = t + timeStep                # go to the next time

#
# The iteration loop has concluded, and the accumulated data points will be put into plots
#

# compose string variables for use in putting text on the plot
tauString = 'Decay time constant = ' + str(tauDecay) + ' seconds'
timeStepString = 'Iteration time step = ' + str(timeStep) + ' seconds'

# code to set up the plots in a single figure
plt.figure()   # start a figure

plt.plot(tEuler, NEuler, 'ro', label='Euler Numerical Solution')  # red dots for the Euler solution plot
plt.plot(tAnalytic, NAnalytic, label='Analytic Solution')         # blue line for analytic solution plot

# check if a log plot is requested
if(useLogPlot):
    plt.yscale('log')
    #
    # Need to calculate a lower plot limit on the vertical scale being obtained with these input parameters
    #
    lowestNucleiNumber = initialNucleiNumber*math.exp(-maximumTime/tauDecay)
    plt.text(0.05*maximumTime, 5.0*lowestNucleiNumber, tauString)       # text to document the time parameters used
    plt.text(0.05*maximumTime, 2.5*lowestNucleiNumber, timeStepString)  # text to document the time parameters used
else:
    plt.text(0.1*maximumTime, 0.75*initialNucleiNumber, tauString)       # text to document the time parameters used
    plt.text(0.1*maximumTime, 0.70*initialNucleiNumber, timeStepString)  # text to document the time parameters used

plt.xlabel('Time (s)')                             # add axis labels
plt.ylabel('Nuclei Remaining')

plt.title('Radioactive Decay: Analytic and Euler Numerical Solutions')
plt.grid(True)
plt.legend(loc=0)
plt.show()                                         # show figure