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

parser.add_argument('--nSteps', default=600000, type=int, help="Number of particle; default 6---")
parser.add_argument('--kSteps', default=600, type=int, help="Steps interval for saving position information; default 600")
parser.add_argument('--stepSize', default=1, type=int, help="Steps size; default 1")
parser.add_argument('--maxBoxUnits', default=25, type=int, help="Maximum number of box units each side; default 120")
parser.add_argument('--holeSize', default=6, type=int, help="Size of escape hole; default 10")
parser.add_argument('--nParticles', default=400, type=int, help="Even number of particle, perfect square; default 400")
parser.add_argument('--seedFixed', action='store_true', help="Use a fixed random number seed; default False")
parser.add_argument('--verbose', action='store_true', help="Extra print out; default False")
parser.add_argument('--teleport', action='store_true', help="Use teleport option at a boundary; default False")
parser.add_argument('--nFrames', default=500, type=int, help="Number of frames in the video; default 500")
parser.add_argument('--nInterval', default=100, type=int, help="Time difference between video frames in milliseconds; default 100")
parser.add_argument('--doMovie', action='store_true', help="Do the movie output, requires ffmpeg library; default False")

args = parser.parse_args()

nParticles = args.nParticles
squareRootParticles = int(np.sqrt(nParticles))
halfSquareRootParticles = -squareRootParticles/2
if(4*halfSquareRootParticles*halfSquareRootParticles != nParticles):
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
print "  Animation will not be done"
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

#
# Counters for boundary checks
#
nTeleported = 0       # Running total of the number of teleported particles, not double counting for teleports
nBoundaryChecks = 0   # Number of calls to boundary check routine when the particle was at a boundary in either X or Y
def checkBounday(position, kParticle):
    global nBoundaryChecks, nTeleported
    #
    # If the position on input is either less than minLattice or more than maxLattice, then a boundary check will be mad
    # Otherwise the same position is returned
    #
    if(position >= minLattice and position <= maxLattice):
       return position

    nBoundaryChecks += 1
    #
    #  Function returns a new position after checking the teleport options
    #
    if(teleport):
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

x_boundary = maxLattice
y_hole_start = 0 - holeSize/2
y_boundary = [x for x in range(int(y_hole_start), int(y_hole_start + holeSize))]
escapeCount = 0


def checkEscape(xpos, ypos):
    global escapeCount, nParticles
    if xpos == x_boundary and y_boundary.__contains__(ypos):
        escapeCount += 1
        nParticles -= 1
        return None, None
    else:
        return xpos, ypos

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

time_count_matrix = {}
iteration = 0
timesteps, counts = [], []
while iteration != nSteps:
    iteration += nSteps/1000
    time_count_matrix[iteration] = 0
print time_count_matrix

while iStep<nSteps:
    if nParticles < 1:
        print "all particles have escaped"
        break
    else:
        kParticle = random.randrange(0, nParticles, 1)     # returns integer from 0 to nParticles - 1, note that the nParticles integer is not in range
        ixy = random.randrange(-1, 3, 2)                   # returns either -1 or +1

        if particlePositionYOld[kParticle] is not None and particlePositionXOld[kParticle] is not None:
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
        particlePositionXOld[kParticle] = checkBounday(particlePositionXNew[kParticle], kParticle)
        particlePositionYOld[kParticle] = checkBounday(particlePositionYNew[kParticle], kParticle)
        particlePositionXOld[kParticle], particlePositionYOld[kParticle] = checkEscape(particlePositionXOld[kParticle], particlePositionYOld[kParticle])

        if time_count_matrix.keys().__contains__(iStep):
            time_count_matrix[iStep] = nParticles
            timesteps.append(iStep)
            counts.append(nParticles)
        elif iStep == nSteps - 1:
            time_count_matrix[nSteps] = nParticles
            timesteps.append(iStep)
            counts.append(nParticles)


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
for x in sorted(time_count_matrix.keys()):
    print x, time_count_matrix[x]

print

print "  particles left at end of simulation:", nParticles
print "  number of escaped particles:", escapeCount

print "  Number of boundary check calls ", nBoundaryChecks
if(teleport):
   print "  Total particles teleported ", nTeleported

kStepFinal = kStep
print "\n  Last snapshot index", kStepFinal

sim_number = input("What number simulation is this? Please enter the number here: ")
file = open("CSV_data_for_random_walk_" + str(sim_number) + "_" + "holesize_" + str(holeSize) + ".csv", "w")
for x in range(timesteps.__len__()):
    file.write(str(timesteps[x]) + "," + str(counts[x]) + '\n')
file.close()
