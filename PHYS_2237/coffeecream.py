# Program 7.6 coffeeCreamDiffusion_Chapter7V4.py Diffusion in a two dimensional surface within a confining square box
#
# Based on coffeeCreamDiffusion_Chapter7V3.py, but simplifies the boundary checking to have fewer code lines
#
# This verison animates the diffusion process in two dimensions, and tags the teleported particles with the red color
#

import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible
import matplotlib.pyplot as plt     # get matplotlib plot functions
import  argparse                    # argument parser library
import numpy as np                  # numerical functions library used by python
import random                       # random number generator library
import matplotlib.animation as animation  # animation library (see examples to learn how to use it)

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()

parser.add_argument('--nSteps', default=600000, type=int, help="Number of particle; default 6,000")
parser.add_argument('--kSteps', default=600, type=int, help="Steps interval for saving position information; default 600")
parser.add_argument('--stepSize', default=1, type=int, help="Steps size; default 1")
parser.add_argument('--maxBoxUnits', default=25, type=int, help="Maximum number of box units each side; default 25")
parser.add_argument('--holeSize', default=10, type=int, help="Size of hole in right side of container, default 10")
parser.add_argument('--nParticles', default=400, type=int, help="Even number of particle, perfect square; default 400")
parser.add_argument('--seedFixed', action='store_true', help="Use a fixed random number seed; default False")
parser.add_argument('--verbose', action='store_true', help="Extra print out; default False")
parser.add_argument('--teleport', action='store_true', help="Use teleport option at a boundary; default True")
parser.add_argument('--nFrames', default=500, type=int, help="Number of frames in the video; default 500")
parser.add_argument('--nInterval', default=100, type=int, help="Time difference between video frames in milliseconds; default 100")
parser.add_argument('--doMovie', action='store_true', help="Do the movie output, requires ffmpeg library; default False")

args = parser.parse_args()

nParticles = args.nParticles
squareRootParticles = int(np.sqrt(nParticles))
halfSquareRootParticles = -squareRootParticles/2
if(4 * halfSquareRootParticles * halfSquareRootParticles != nParticles):
    print "\n The number of particle ", nParticles, " is not an even perfect square"
    exit()

nSteps = args.nSteps
kSteps = args.kSteps
nSnapShots = int(nSteps/kSteps)
nSnapShotsPlusOne = nSnapShots + 1
stepSize = args.stepSize
maxLattice = (args.maxBoxUnits)*stepSize
minLattice = -maxLattice
minusTwiceLattice = 2*minLattice
seedFixed = args.seedFixed
verbose = args.verbose
teleport = args.teleport
holeSize = args.holeSize

time_count_matrix = {}
iteration = 0
while iteration != nSteps:
    iteration += nSteps/10
    time_count_matrix[iteration] = 0
print time_count_matrix

nFrames = args.nFrames
if(nFrames > nSnapShots):
    nFrames = nSnapShots

nInterval = args.nInterval
doMovie = args.doMovie  # retrieve the choice to make the mp4 movie
nSnapShotsSkip = int(nSnapShots/nFrames)
if(nSnapShotsSkip < 1):
    nSnapShotsSkip = 1

print "\n  Diffusion in two dimensions with ", nParticles, " particles taking ", nSnapShots, " snapshots saved for display"
print "  Total number of steps ", nSteps, ", corresponding to ", nSteps/nParticles, " steps per particle"
print "  Step size in 2D", stepSize, " with box extent from ", minLattice, " to ", maxLattice, " in each (x,y) box direction"

if(seedFixed):
    random.seed(3)  # fixes the sequence of random sees for repeatability of numerical results
    print "  The random number seed is fixed"
else:
    random.seed()   # random seed is set according to the system clock, numerical results will be non-repeatable
    print "  The random number seed changes with the clock time"
print "  Animation will have ", nFrames, " frames with a snapshot time skip ", nSnapShotsSkip, " using ", nInterval, " milliseconds between frames"

if(teleport):
    print "  Teleportation is used at the box boundaries"
else:
    print "  Teleportation is not used at the box boundaries"

#
# Define the old and new particle position arrays array
#
particlePositionXOld = [0 for kParticle in range(nParticles)]
particlePositionYOld = [0 for kParticle in range(nParticles)]
particlePositionXNew = [0 for kParticle in range(nParticles)]
particlePositionYNew = [0 for kParticle in range(nParticles)]

#
# Define the cumulative particle array at each snapshot step particlePositionCumulativeQ[kStep][kParticle]
#
particlePositionCumulativeX = [[0 for kParticle in range(nParticles)] for kStep in range(nSnapShotsPlusOne)]
particlePositionCumulativeY = [[0 for kParticle in range(nParticles)] for kStep in range(nSnapShotsPlusOne)]
particlePositionTeleported = [False for kParticle in range(nParticles)]
particlePositionCumulativeXTeleported = [[minusTwiceLattice for kParticle in range(nParticles)] for kStep in range(nSnapShotsPlusOne)]
particlePositionCumulativeYTeleported = [[minusTwiceLattice for kParticle in range(nParticles)] for kStep in range(nSnapShotsPlusOne)]
countCumulativeTeleported = [0 for kStep in range(nSnapShotsPlusOne)]

particleEscape = [[False, False] for kParticle in range(nParticles)]
y_hole_start = 0 - holeSize/2
y_hole = [x for x in range(int(y_hole_start), int(y_hole_start + holeSize))]


#
# Counters for boundary checks
#
nTeleported = 0       # Running total of the number of teleported particles, not double counting for teleports
nBoundaryChecks = 0   # Number of calls to boundary check routine when the particle was at a boundary in either X or Y
escapeCount = 0       # Running total of number of escaped particles
def checkBounday(position, kParticle, dimension):
    global nBoundaryChecks, nTeleported, escapeCount, nParticles, particleEscape
    #
    # If the position on input is either less than minLattice or more than maxLattice, then a boundary check will be made
    # If the position on input is at the escape boundary, the particle is deleted
    # Otherwise the same position is returned
    #
    if(position >= minLattice and position < maxLattice):
       if dimension == "x":
           particleEscape[kParticle][0] = False
       else:
           particleEscape[kParticle][1] = False
       return position
    elif position == maxLattice and dimension == "x":
        particleEscape[kParticle][0] = True
        print "a particle is at the x boundary"
    elif y_hole.__contains__(position) and dimension == "y":
        particleEscape[kParticle][1] = True
        print "a particle is at the y boundary"
    delete = False
    if particleEscape[kParticle] == [True, True]:
        delete = True
        print "a particle escapes!"
        escapeCount += 1
        nParticles -= 1
        return
    else:
        nBoundaryChecks += 1
        #
        #  Function returns a new position after checking the teleport options
        #
        if delete:
            print delete
        if teleport and not delete:
            if(particlePositionTeleported[kParticle] == False):
                particlePositionTeleported[kParticle] = True
                nTeleported += 1       # do not double count the number of teleported particles
            if(position > maxLattice):
                newPosition = minLattice   # teleport from maximum boundary to minimum boundary
            else:
                newPosition = maxLattice   # teleport from minimum boundary to maximum boundary

        else:
            # no teleport option, particle stays at the same coordinate
            if(position > maxLattice):
                newPosition = maxLattice
            else:
                newPosition = minLattice

        return newPosition

#
# initialize particle positions in central square squareRootParticles units in X and squareRootParticles units in Y
#
initX = halfSquareRootParticles
kParticle = 0
for iX in range(squareRootParticles):
    initY = halfSquareRootParticles
    for iY in range(squareRootParticles):
        particlePositionXOld[kParticle] = initX
        particlePositionYOld[kParticle] = initY
        #
        # Initial snapshot
        #
        particlePositionCumulativeX[0][kParticle] = initX
        particlePositionCumulativeY[0][kParticle] = initY
        initY += 1
        kParticle += 1
    
    initX += 1

iStep = 1
kStep = 0
jStep = 0
debug = True
while iStep<nSteps:
    kParticle = random.randrange(0, nParticles, 1)     # returns integer from 0 to nParticles - 1, note that the nParticles integer is not in range
    ixy = random.randrange(-1, 3, 2)                   # returns either -1 or +1

    if  particlePositionXOld != None and particlePositionYOld != None:
        if(ixy == -1):
            #
            # move in the X direction
            #
            particlePositionYNew[kParticle] = particlePositionYOld[kParticle]
            iMinusPlus = random.randrange(-1, 3, 2)
            if(iMinusPlus == -1):
                particlePositionXNew[kParticle] = particlePositionXOld[kParticle] + 1
            else:
                particlePositionXNew[kParticle] = particlePositionXOld[kParticle] - 1
        else:
            #
            # move in the Y direction
            #
            particlePositionXNew[kParticle] = particlePositionXOld[kParticle]
            iMinusPlus = random.randrange(-1, 3, 2)
            if(iMinusPlus == -1):
                particlePositionYNew[kParticle] = particlePositionYOld[kParticle] + 1
            else:
                particlePositionYNew[kParticle] = particlePositionYOld[kParticle] - 1
        
    #
    # Replace the previous Old position with the New position
    # Check for the particle being at the boundary
    #
    particlePositionXOld[kParticle] = checkBounday(particlePositionXNew[kParticle], kParticle, "x")
    particlePositionYOld[kParticle] = checkBounday(particlePositionYNew[kParticle], kParticle, "y")
    #
    # Check new position to see if particle has escaped
    #
    #

    jStep += 1
    if(jStep == kSteps):
        jStep = 0
        kStep += 1
        if(kStep >= nSnapShotsPlusOne):
            print "\n Index error: kStep = ", kStep, " at limit nSnapShots + 1", nSnapShotsPlusOne, ",  iStep ", iStep
            exit()
    
        countCumulativeTeleported[kStep] = nTeleported
        #
        # Take snapshot of all the particles at this step
        #
        for jParticle in range(nParticles):
            particlePositionCumulativeX[kStep][jParticle] = particlePositionXOld[jParticle]
            particlePositionCumulativeY[kStep][jParticle] = particlePositionYOld[jParticle]
            if(particlePositionTeleported[jParticle]):
                particlePositionCumulativeXTeleported[kStep][jParticle] = particlePositionXOld[jParticle]
                particlePositionCumulativeYTeleported[kStep][jParticle] = particlePositionYOld[jParticle]
                    
    iStep += 1
print particleEscape

print "\n  Number of boundary check calls ", nBoundaryChecks
if(teleport):
   print "  Total particles teleported ", nTeleported
   print "  Total particles escaped", escapeCount

kStepFinal = kStep
print "\n  Last snapshot index", kStepFinal

xDisplacement = [particlePositionCumulativeX[kStepFinal][jParticle] - particlePositionCumulativeX[0][jParticle] for jParticle in range(nParticles)]
yDisplacement = [particlePositionCumulativeY[kStepFinal][jParticle] - particlePositionCumulativeY[0][jParticle] for jParticle in range(nParticles)]
rDisplacement = [np.sqrt(xDisplacement[jParticle]*xDisplacement[jParticle] + yDisplacement[jParticle]*yDisplacement[jParticle]) for jParticle in range(nParticles)]
rDisplacementMean = np.mean(rDisplacement)
rDisplacementStd = np.std(rDisplacement)
print "  Mean displacement ", rDisplacementMean, " +/- ", rDisplacementStd, " units"

# code to set up one plots in a single figure
fig = plt.figure(1)          # start a figure
ax = fig.add_subplot(111)    # this sets the upper half plot for initial snapshot
teleportPlot, = ax.plot([], [], 'r.')   # this is an empty scatter object used in the initialization of the animation
linePlot, = ax.plot([], [], 'b.')   # this is an empty scatter object used in the initialization of the animation
timeValue_text = ax.text(0.45, 0.95, '', transform=ax.transAxes, color='blue')   # this is a placeholder for where the time step will be updated
rDisplacementValue_text = ax.text(0.45, 0.90, '', transform=ax.transAxes, color='blue')   # this is a placeholder for where the mean displacement will be updated
teleportValue_text = ax.text(0.45, 0.85, '', transform=ax.transAxes, color='red')   # this is a placeholder for where number of teleported particles will be updated

def init():                 # this initialization function tells what is being animiated, i.e. the y(x,t) from x=0 to x=xMax at a given time
    """initialize animation"""
    global linePlot
    
    linePlot.set_data([], [])               # (x,y) at a given time step
    teleportPlot.set_data([], [])               # (x,y) at a given time step
    timeValue_text.set_text('')             # time step value
    rDisplacementValue_text.set_text('')    # mean r displacement value
    teleportValue_text.set_text('')         # number of teleported particles
    
    return linePlot, teleportPlot, timeValue_text, rDisplacementValue_text, teleportValue_text

def teleportCheck(element):
    #
    # Function used in animation routine to get the positions of the teleported particles at a given time step
    # This function is used in a "filter" function in Python
    #
    return element > minusTwiceLattice

def animate(i):   # this is the function this being animated, the i index is increased according the range of the frames value below
    """perform animation step"""           # i ranges from 0 to (nFrames - 1) as given in the FuncAnimation function
    global nSnapShotsSkip, particlePositionCumulativeX, particlePositionCumulativeY
    
    ii = i*nSnapShotsSkip
    if(ii > nSnapShots):
        print "\n Index error for ii ", ii, " at limit ", nSnapShots, " for i ", i
        exit()

    xPlot = []
    yPlot = []
    for jParticle in range(nParticles):
        if(particlePositionCumulativeXTeleported[ii][jParticle] == minusTwiceLattice):  # check if the teleportation array is still unchanged
            xPlot.append(particlePositionCumulativeX[ii][jParticle])
            yPlot.append(particlePositionCumulativeY[ii][jParticle])

    linePlot.set_data(xPlot,yPlot)  # the scatter plot contains the (x,y) positions for a given time step value

    xPlotTest = [particlePositionCumulativeXTeleported[ii][jParticle] for jParticle in range(nParticles)]
    yPlotTest = [particlePositionCumulativeYTeleported[ii][jParticle] for jParticle in range(nParticles)]
    xPlotTeleport = filter(teleportCheck, xPlotTest)
    yPlotTeleport = filter(teleportCheck, yPlotTest)
    teleportPlot.set_data(xPlotTeleport,yPlotTeleport)  # the scatter plot contains the (x,y) positions for a given time step value

    timeStepString = '%s %.3e' % ('Time step ', ii*kSteps)
    timeValue_text.set_text(timeStepString)
    xDisplacementSnap = [particlePositionCumulativeX[ii][jParticle] - particlePositionCumulativeX[0][jParticle] for jParticle in range(nParticles)]
    yDisplacementSnap = [particlePositionCumulativeY[ii][jParticle] - particlePositionCumulativeY[0][jParticle] for jParticle in range(nParticles)]
    rDisplacementSnap = [np.sqrt(xDisplacementSnap[jParticle]*xDisplacementSnap[jParticle] + yDisplacementSnap[jParticle]*yDisplacementSnap[jParticle]) for jParticle in range(nParticles)]
    rDisplacementMean = np.mean(rDisplacementSnap)
    rDisplacementString = '%s %.1f' % ('Average radial displacement ', rDisplacementMean)
    rDisplacementValue_text.set_text(rDisplacementString)
    teleportValueString = '%s %d' % ('Teleported particles ', countCumulativeTeleported[ii])
    teleportValue_text.set_text(teleportValueString)

    return linePlot, teleportPlot, timeValue_text, rDisplacementValue_text, teleportValue_text # this function must return the same number and type variables as the init function

plt.ylim(minLattice, maxLattice)
plt.xlim(minLattice, maxLattice)
plt.xlabel('X')                             # add axis labels
plt.ylabel('Y')
plt.grid(True)
if(teleport):
    titleString = '%s %d %s %d %s' % ('Diffusion of ', nParticles, ' Particles in 2D for ', nSteps, ' Steps, Use Teleporation')
else:
    titleString = '%s %d %s %d %s' % ('Diffusion of ', nParticles, ' Particles in 2D for ', nSteps, ' Steps, No Teleporation')
plt.title(titleString)

print "\n Animation step is being done"
ani = animation.FuncAnimation(fig, animate, frames=nFrames, interval=nInterval, blit=True, init_func=init)

#
# The ffmpeg binary must be put in the executable PATH, where the binary can be obtained from http://www.ffmpegmac.net/
# Mathplotlib 2.0 version requires that the fps and the extra_args be placed in the FFMpegWriter call, not the ani.save call
# Older Matplotlib will work with either method
#
if(doMovie):
    print "\n Movie production step is being done"
    FFwriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
    ani.save('coffeeCreamDiffusion.mp4', writer = FFwriter)
    print " Movie production step is completed"
else:
    print "\n Movie production step is being skipped"

plt.show()
