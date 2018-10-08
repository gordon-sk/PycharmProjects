#
# Program 6.4 Solving the wave equation for a composite string, as shown in the textbook Figure 6.3 on page 163
# Derived from the 6.3 optional free end version (fixedString_Chapter6V3.py)
#
# Waves propagating along a string with fixed ends
# The iteration equation 6.6 is solved after specifying the wave at two prior times
# The boundary conditions are that the string is fixed at both ends.
#
# 1) Input parameters
#    a) Default wave speed is 5 meters/second
#    b) Default rRatio is 1.0, as required from page 162
#    c) Default string length is 1.0 meter
#    d) Default step size in position (x) is 0.005 meter
#    e) Default maximum time is 1.00 seconds
#    f) Default initial wave shape is given by a Gaussian exp(-k(x-x0)^2) with k = 1000 (meters)^-2 and x0 = 0.3 meters
#    g) Default initial amplitude of the wave is 1.0
#    h) Default ratio of the two wave speeds is 0.75
#    i) Default position of the split wire is 0.50 of the wire length
#
# 2) Implementation method
#    a) The time step is computed from the product of the rRatio and the position steps size, divided by the speed
#    b) The number of time steps is the maximum time divided by the time step
#
# 3) Output (see Figure 6.3, page 163)
#    a) Animation of the wave propagation
#    b) Movie of the animation if requested

import  argparse                    # argument parser library
import math as mp                   # math functions library used by python
import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt     # get matplotlib plot functions
import matplotlib.animation as animation  # animation library (see examples to learn how to use it)

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#

parser = argparse.ArgumentParser()

parser.add_argument('--cSpeed', default=5.0, type=float, help="Speed of wave along the string in m/s; default 5.0")
parser.add_argument('--rRatio', default=1.0, type=float, help="Ratio cSpeed*DeltaX/DeltaT, page 162; default 1.0")
parser.add_argument('--xLength', default=1.0, type=float, help="Length of string in m; default 1.0")
parser.add_argument('--deltaX', default=0.005, type=float, help="Step length in space; default 0.005")
parser.add_argument('--splitFraction', default=0.5, type=float, help="Fraction along the string length where the split happens; default 0.5")
parser.add_argument('--vRatio', default=0.75, type=float, help="Ratio of wave speeds in composite string; default 0.75")
parser.add_argument('--x0', default=0.3, type=float, help="Center of initial wave pulse, in m; default 0.3")
parser.add_argument('--kFactor', default=1000.0, type=float, help="Wave number k, in inverse m; default 1000.0")
parser.add_argument('--gaussAmplitude', default=1.0, type=float, help="Initial pulse amplitude in m; default 1.0")
parser.add_argument('--maxT', default=1.00, type=float, help="Maximum simulation time in s; default 1.00")
parser.add_argument('--freeEnd', action='store_true', help="Allow the right end of the string to be free; default False")
parser.add_argument('--verbose', action='store_true', help="Give extra printout; default False")
parser.add_argument('--nFrames', default=500, type=int, help="Number of frames in the video; default 500")
parser.add_argument('--nInterval', default=100, type=int, help="Time difference between video frames in milliseconds; default 100")
parser.add_argument('--doMovie', action='store_true', help="Do the movie output, requires ffmpeg library; default False")

args = parser.parse_args()

cSpeed = args.cSpeed        # wave speed along the string in m/s
rRatio = args.rRatio        # ratio factor c*DeltaX/DeltaT (page 162)
xStep = args.deltaX         # step size along x in m
xMax = args.xLength         # length of string in m
xMin = 0.0                  # hardcoded number, no reason to have it different from 0

splitFraction = args.splitFraction # fraction along the string where the composite split occurs
splitDistance = splitFraction*xMax # distance along the string where the composite split occurs
vRatio = args.vRatio        # ratio of the wave speeds in the composite string

nStepsX = int(xMax/xStep)   # number of steps along x direction
initialX0 = args.x0         # initial centroid of the wave pulse in m
kFactor = args.kFactor      # wave number of the wave pulse in inverse square meters
gaussAmplitude = args.gaussAmplitude   # initial pulse amplitude
verbose = args.verbose      # extra printout option
freeEnd = args.freeEnd      # option for the right end of the string to be free
freeEnd = True

tMin = 0.0                  # start time (seconds), hardcoded number with no reason to be different from 0
tMax = args.maxT            # final time (seconds)
tStep = rRatio*xStep/cSpeed # calculated time step
nStepsT = int(tMax/tStep)   # number of t Steps, computed from tStep

nFrames = args.nFrames
nInterval = args.nInterval
doMovie = args.doMovie  # retrieve the choice to make the mp4 movie

nSkipPoints = int(nStepsT/nFrames)
if(nSkipPoints < 1):
    nSkipPoints = 1
    nFrames = nStepsT

print "\n  Program to solve the wave equation for a composite string and to produce an animation"
print "  Wave speed ", cSpeed, " m/s, final time ", tMax, " s in time steps ", tStep
print "  The change in the string density occurs at ", splitDistance, "m, and left/right speed ratio is ", vRatio
print "  Length of string ", xMax, " m, with spatial steps ", xStep, " m"
print "  The rRatio = ", rRatio
print "  The inital wave pulse is a Gaussian with a center at ", initialX0, " m, a wave number k ", kFactor, " in inverse square meters, and an amplitude ", gaussAmplitude, " m"
print "  Time interval between animation frames = ", nSkipPoints*tStep, " s"
if(freeEnd):
    print "  The right end of the string can move freely"
else:
    print "  Both ends of the string are fixed at y = 0"

def initializeWave(yxt1, yxt2, xMin, xStep, kFactor, initialX0):
    global nStepsX
    #
    # initial shape as on page 160
    # y(x, t=0) = exp[-kFactor*(x-initialX0)**2]
    #
    xValue = xMin
    maxAmplitude= 0.0
    for i in range(nStepsX):
        yxt1[i] = gaussAmplitude*mp.exp(-kFactor*(xValue - initialX0)*(xValue - initialX0))
        yxt2[i] = yxt1[i]
        if(abs(yxt1[i]) > maxAmplitude):
            maxAmplitude = abs(yxt1[i])

        xValue += xStep
        
    print "\n Initial maximum amplitude = ", maxAmplitude, "m"
    return

def initializeEnds(yxt1, yxt2, yxt3):
    global nStepsX
    yxt1[0] = 0.0
    yxt2[0] = 0.0
    yxt3[0] = 0.0
    
    nStepsXMinusOne = nStepsX - 1
    yxt1[nStepsXMinusOne] = 0.0
    yxt2[nStepsXMinusOne] = 0.0
    yxt3[nStepsXMinusOne] = 0.0
    
    if(freeEnd):
        yxt1[nStepsXMinusOne] = yxt2[nStepsXMinusOne - 1] # first suggestion in problem 6.1
        yxt2[nStepsXMinusOne] = yxt3[nStepsXMinusOne - 1] # first suggestion in problem 6.1
    return

def propagateWave(yxt1, yxt2, yxt3):
    global nStepsX, rRatio
    #
    # propagation formula on page 159 to get next time step wave form
    #
    nStepsXMinusOne = nStepsX - 1
    rSquared = rRatio*rRatio
    TwiceOneMinusRSquared = 2.0 - 2.0*rSquared
    
    i = 1
    while i < nStepsXMinusOne:
        #
        # loop over steps in x to get next time step yxt3
        #
        yxt3[i] = TwiceOneMinusRSquared*yxt2[i] - yxt1[i] + rSquared*(yxt2[i+1] + yxt2[i-1])
        
        i = i+1
    
    if(freeEnd):
        yxt3[nStepsXMinusOne] = yxt3[nStepsXMinusOne - 1] # first suggestion in problem 6.1
    
    return

def shiftY(yxt1, yxt2, yxt3):
    global nStepsX
    #
    #  Do shifting yxt2 -> yxt1  and  yxt3 -> yxt2
    #
    nStepsXMinusOne = nStepsX - 1
    i = 1
    while i < nStepsXMinusOne:
        yxt1[i] = yxt2[i]
        yxt2[i] = yxt3[i]
        i = i+1

    if(freeEnd):
        yxt1[nStepsXMinusOne] = yxt2[nStepsXMinusOne - 1] # first suggestion in problem 6.1
        yxt2[nStepsXMinusOne] = yxt3[nStepsXMinusOne - 1] # first suggestion in problem 6.1

    return

def init():                 # this initialization function tells what is being animiated, i.e. the y(x,t) from x=0 to x=xMax at a given time
    """initialize animation"""
    global linePlot
    
    linePlot.set_data([], [])             # y(x,t)
    timeValue_text.set_text('')   # time value
    
    return linePlot, timeValue_text

def animate(i):   # this is the function this being animated, the i index is increased according the range of the frames value below
    """perform animation step"""           # i ranges from 0 to (nFrames - 1) as given in the FuncAnimation function
    global timeSlice, xPosition, yPosition, nStepsX, nSkipPoints, nDecades, nStepsT, linePlot, totalDataPoints
    
    ii = nSkipPoints*i
    if(ii >= nStepsT - 1):
        ii = nStepsT - 1
    kPointStart = ii*nStepsX

    if(kPointStart < 0 or kPointStart + nStepsX > totalDataPoints):
        print "\n ***Error i animate function: i = ", i, ",  ii = ", ii, ",  nStepsX = ", nStepsX, ",  kPointStart ", kPointStart, ", kPointEnd = ", kPointStart + nStepsX, ", totalDataPoints = ", totalDataPoints
        exit()  # check that all the indexing is within bounds

    timeValue = timeSlice[kPointStart]
    xPlot = [xPosition[kPointStart + i] for i in range(nStepsX)]
    yPlot = [yPosition[kPointStart + i] for i in range(nStepsX)]
    linePlot.set_data(xPlot,yPlot)  # the "line" contains a y(x,t) wave for a given time value
    timeValue_text.set_text('Time = %.4f seconds'% timeValue)
 
    return linePlot, timeValue_text  # this function must return the same number and type variables as the init function

def animateFixedString(nStepsX):
    global timeSlice, xPosition, yPosition
    #
    # Function to produce an animation of the traveling wave
    #

    yMax = max(yPosition)
    yMin = min(yPosition)
    print "\n The animation yMin = ", yMin, " and the yMax = ", yMax
    yRange = yMax - yMin
    plt.xlim(0.0, xMax)
    plt.ylim(yMin - 0.1*yRange, yMax + 0.1*yRange)
    leftBaselinePlotX = [0.0, splitDistance]
    rightBaselinePlotX = [splitDistance, xMax]
    baselinePlotY = [0.0, 0.0]
    ax.plot(leftBaselinePlotX, baselinePlotY, 'k-', lw=3)
    ax.plot(rightBaselinePlotX, baselinePlotY, 'g-', lw=5)

    if(freeEnd):
        plt.title('Wave Propagation in a Composite String with One Free Endpoint')
    else:
        plt.title('Wave Propagation in a Composite String with Fixed Endpoints')

    plt.xlabel('x Coordinate (m)')                             # add axis labels
    plt.ylabel('y Coordinate (m)')
    plt.grid(True)
    
    print "\n Animation step is being done"
    ani = animation.FuncAnimation(fig, animate, frames=nFrames,
                                  interval=nInterval, blit=True, init_func=init)

#
# The ffmpeg binary must be put in the executable PATH, where the binary can be obtained from http://www.ffmpegmac.net/
# Mathplotlib 2.0 version requires that the fps and the extra_args be placed in the FFMpegWriter call, not the ani.save call
# Older Matplotlib will work with either method
#
    if(doMovie):
        print "\n Movie production step is being done"
        FFwriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
        ani.save('waveOnAFixedCompositeString.mp4', writer = FFwriter)
        print " Movie production step is completed"
    else:
        print "\n Movie production step is being skipped"

    plt.show()

    return

#
# Begin main program
#
yxt1 = [0 for i in range(nStepsX)]  # previous time slice for iteration
yxt2 = [0 for i in range(nStepsX)]  # current time slice for iteration
yxt3 = [0 for i in range(nStepsX)]  # next time slice for iteration
 
nDecadeT = int(nStepsT/10)
currentTime = tMin

initializeWave(yxt1, yxt2, xMin, xStep, kFactor, initialX0)
initializeEnds(yxt1, yxt2, yxt3)

totalDataPoints = nStepsX*(nStepsT + 1)

timeSlice = [0 for i in range(totalDataPoints)]
xPosition = [0 for i in range(totalDataPoints)]
yPosition = [0 for i in range(totalDataPoints)]

kPoint = 0
for i in range(nStepsX):
    timeSlice[kPoint] = currentTime
    xPosition[kPoint] = i*xStep
    yPosition[kPoint] = yxt2[i]
    kPoint += 1

for n in range(nStepsT):
    propagateWave(yxt1, yxt2, yxt3)  # generate the wave for the next time value using the two previous wave forms
    
    shiftY(yxt1, yxt2, yxt3)
    currentTime += tStep
    for i in range(nStepsX):
        if(kPoint >= totalDataPoints):
            print "\n  nTime ", n, ",  current time ", currentTime, ",  i ", i, ",  kPoint ", kPoint, ",  total data points ", totalDataPoints
            exit()
        
        timeSlice[kPoint] = currentTime
        xPosition[kPoint] = i*xStep
        yPosition[kPoint] = yxt2[i]
        kPoint += 1

totalDataPoints = kPoint
maxAmplitude = 0.0
for i in range(nStepsX):
    if(abs(yxt2[i]) > maxAmplitude):
        maxAmplitude = abs(yxt2[i])
            
print "\n Final maximum amplitude = ", maxAmplitude, " m, with number of data points stored = ", totalDataPoints

fig = plt.figure(1)       # start a figure for cascade snapshots
ax = fig.add_subplot(111)
linePlot, = ax.plot([], [], 'b-', lw=2)                                # this is an empty line object used in the initialization of the animation
timeValue_text = ax.text(0.65, 0.95, '', transform=ax.transAxes, color='red')   # this is a placeholder for where the time will be updated
waveSpeedString = '%s %.3f %s' % ("Wave speed ", cSpeed, " m/s")
waveSpeed_text = ax.text(0.05, 0.95, waveSpeedString, transform=ax.transAxes, color='red')
xStepString = '%s %.3f %s' % ("Step size in x ", xStep, " m")
xStep_text = ax.text(0.05, 0.89, xStepString, transform=ax.transAxes, color='red')
tStepString = '%s %.3f %s' % ("Step size in t ", tStep, " s")
tStep_text = ax.text(0.05, 0.83, tStepString, transform=ax.transAxes, color='red')
splitString = '%s %.3f %s' % ("Composite split is at ", splitDistance, " m")
split_text = ax.text(0.05, 0.77, splitString, transform=ax.transAxes, color='green')
vRatioString = '%s %.3f' % ("Left/Right speed ratio ", vRatio)
vRatio_text = ax.text(0.05, 0.71, vRatioString, transform=ax.transAxes, color='red')

animateFixedString(nStepsX)

