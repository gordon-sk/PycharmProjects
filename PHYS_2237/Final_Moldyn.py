# Program 9.4 molecularDynamics_Chapter9Final
#
#  Molecular dynamics calculation using Lennard-Jones potential, for use in Problem 1 of the Final Exam
#
#  Based on molecularDynamics_Chapter9V3, The V3 version has the option to add or subtract kinetic energy using Equation 9.13 in the textbook
#                                         This V3 version can write a final state data (positions and velocities) text file
#                                         The data text file can be used in a subsequent calculation as the initialization data
#
#  V(r) = 4*epsilon*[(sigma/r)**12 - (sigma/r)**6], Equation 9.6
#  Set unit of energy according to epsilon = 1, where for Argon epsilon/kB = 120 K
#  Lengths are in units of sigma = 3.4 Angstroms Argon value
#  Time steps are in units of sigma*sqrt(m/epsilon) = 1.8 picoseconds for Argon
#  Velocity is in units of sqrt(epsilon/m), where m is in units of Argon mass
#  By default the mass has one unit
#
#  Producing Figures 9.3 and 9.4, pages 280-81

#
#  1) Program input
#     a) Optional inputs are available for the number of particles, the size of the confining box,
#        the maximum time, the time step, the cut-off radius for interactions, and the initial randomization
#     b) Defaults are to have 25 particle, a box size of 10 sigma units,
#        maximum time of 20 units, a time step of 0.01 units, and standard randomization start
#
#  2) Program operation
#     a) The data structures for the step-by-step particle position information,
#        and for the particle-summed global information are created.
#     b) The particle data structure is initialized according to the guidelines in the textbook
#     c) The particle movements are performed for each time step, and the total energy is
#        calculated for each time step.  Since the forces are conservative, then the total
#        energy must be conserved.
#     d) When the particle crosses a boundary there are three possibilities: reflection, teleportation, or escaping the box
#
#

import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible
import matplotlib.pyplot as plt     # get matplotlib plot functions
import  argparse                    # argument parser library
import numpy as np                  # numerical functions library used by python
import random                       # random number generator library
import scipy.optimize
from datetime import datetime
import matplotlib.animation as animation
import scipy.optimize as optimization # data fitting library

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()

parser.add_argument('--nParticles', default=36, type=int, help="Number of particles (perfect square); default 36")
parser.add_argument('--boxSize', default=12.0, type=float, help="Box size in sigma units; default 12.0")
parser.add_argument('--maxT', default=240.0, type=float, help="Maximum time; default 120.0")
parser.add_argument('--deltaT', default=0.005, type=float, help="Time step; default 0.005")
parser.add_argument('--v1', default=1.00, type=float, help="First initial speed choice; default 1.0")
parser.add_argument('--v2', default=1.00, type=float, help="Second initial speed choice; default 1.0")
parser.add_argument('--useDataInitializationFile', action='store_false', help="Use a data initialization file; default False")
parser.add_argument('--rFactor', default=-0.0001, type=float, help="R > 0 increases KE, R < 0 decreases kinetic energy; default 0.0")
parser.add_argument('--useTwoSpeeds', action='store_true', help="Use both speed in the initialization; default False")
parser.add_argument('--changeKE', action='store_false', help="Use the rFactor to change the KE; default False")
parser.add_argument('--fixCMVelocity', action='store_true', help="Insure zero CM velocity at the initialization; default False")
parser.add_argument('--noTeleport', action='store_true', help="Do not use teleport condition at boundary; default False")
parser.add_argument('--specular', action='store_true', help="Use specular reflection at boundary; default False")
parser.add_argument('--noPaths', action='store_true', help="Do not plot the particle paths; default False")
parser.add_argument('--seedFixed', action='store_true', help="Use a fixed random number seed; default False")
parser.add_argument('--verbose', action='store_true', help="Extra print out; default False")
parser.add_argument('--sigma', default=1.0, type=float, help="Lennard-Jones potential parameter; default 1.0")
parser.add_argument('--forceConstant', default=1.0, type=float, help="Lennard-Jones force strength; default 1.0")
parser.add_argument('--deltaR', default=0.5, type=float, help="Random initialization offset; default 0.5")
parser.add_argument('--rMinTestFactor', default=0.5, type=float, help="Closest safe distance in terms of sigma; default 0.5")
parser.add_argument('--deltaTestFactor', default=0.075, type=float, help="Largest save position change in units of sigma; default 0.075")
parser.add_argument('--rCut', default=3.0, type=float, help="r cutoff distance in sigma units; default 3.0")
parser.add_argument('--xOffSetFinal', default=+3.0, type=float, help="Offset in X to reposition the position plots; default 0.0")
parser.add_argument('--yOffSetFinal', default=-4.5, type=float, help="Offset in Y to reposition the position plots; default 0.0")
parser.add_argument('--printEnergyStep', action='store_false', help="Print out energies at fixed intervals; default False")
parser.add_argument('--printForceStep', action='store_true', help="Print out net forces at fixed intervals; default False")
parser.add_argument('--noAnimation', action='store_false', help="Skip animation, produce graphs; default False")
parser.add_argument('--nFrames', default=500, type=int, help="Number of frames in the video; default 500")
parser.add_argument('--snapShotFactor', default=5.0, type=float, help="Snapshots factor; default 5.0")
parser.add_argument('--nInterval', default=100, type=int, help="Time difference between video frames in milliseconds; default 100")
parser.add_argument('--doMovie', action='store_true', help="Do the movie output, requires ffmpeg library; default False")

args = parser.parse_args()
nParticles = args.nParticles
if(nParticles < 4):
    print "\n ***Error, cannot have fewer than 4 particles***"
    exit()
squareRootParticles = int(np.sqrt(nParticles))
if(nParticles - squareRootParticles*squareRootParticles != 0):
    print "\n ***Error, particle number is not a perfect square ", nParticles, "***"
    exit()
nParticlesMinusOne = nParticles - 1
nDistinctPairs =  nParticles*(nParticles-1)/2

sigma = args.sigma
rMinTestFactor = args.rMinTestFactor
rMinTest = rMinTestFactor*sigma
deltaTestFactor = args.deltaTestFactor
deltaTest = deltaTestFactor*sigma
deltaR = args.deltaR

boxSize = args.boxSize
xMax = boxSize
yMax = boxSize
xMin = 0.0
yMin = 0.0
xMaxMinusMin = xMax - xMin
yMaxMinusMin = yMax - yMin
maxDiffX = 0.5*(xMaxMinusMin)
maxDiffY = 0.5*(yMaxMinusMin)
xCenter = boxSize/2.0
yCenter = boxSize/2.0

xOffSetFinal = args.xOffSetFinal
yOffSetFinal = args.yOffSetFinal

nBoundaryCount = 0

maxT = args.maxT
deltaT = args.deltaT
timeSteps = int(maxT/deltaT)
timeStepsPlusOne = timeSteps + 1
twiceDeltaT = 2*deltaT
deltaTSquared = deltaT*deltaT

useDataInitializationFile = args.useDataInitializationFile

v1 = args.v1
v2 = args.v2
useTwoSpeeds = args.useTwoSpeeds

rFactor = args.rFactor
rChangeFactor = 1.0+rFactor*deltaT/0.005      # normalize the change rate to 0.005 time step
changeKE = args.changeKE
if(changeKE == False or rFactor == 0.0):
    changeKEString = 'There was no Kinetic Energy modification factor'
    changeKEStringV0 = 'There was no Kinetic Energy modification'
    keNameText = 'NoKEChange'
else:
    if(rFactor < 0):
        changeKEString = 'The Kinetic Energy was continuously reduced by a factor %.6f' % rChangeFactor
        changeKEStringV0 = 'KE reduction factor was %.6f' % rChangeFactor
        keNameText = 'KEDrain %.6f' % rChangeFactor
    else:
        changeKEString = 'The Kinetic Energy was continuously increased by a factor %.6f' % rChangeFactor
        changeKEStringV0 = 'KE increase factor was %.6f' % rChangeFactor
        keNameText = 'KEAdd%.6f' % rChangeFactor
textFileNameString = 'FinalPositions-Particles%d-%s-BoxSize%d-T%dMilli-DeltaT%dMilli-FirstV%dSecondV2%d.txt' % (nParticles, keNameText, boxSize, 1000*maxT, 1000*deltaT, 1000*v1, 1000*v2)

forceConstant = args.forceConstant
forceConstant4 = 4*forceConstant
forceConstant24 = 24*forceConstant

rCut = args.rCut
seedFixed = args.seedFixed
verbose = args.verbose
printEnergyStep = args.printEnergyStep
printForceStep = args.printForceStep

if(seedFixed):
    random.seed(3)  # fixes the sequence of random sees for repeatability of numerical results
else:
    random.seed()   # random seed is set according to the system clock, numerical results will be non-repeatable

print "\n      Molecular Dynamics Simulation"
print "  Number of particles ", nParticles, " confined in a box size ", boxSize, " with number of distinct pairs ", nDistinctPairs
print "  Maximum time ", maxT, " in time step size ", deltaT, " with a total ", timeSteps, " number of time steps"
print "  The interaction force is cut off at ", rCut, " sigma units"
if(useTwoSpeeds):
    print "  The particles are initialized with two choices of initial speed: v1 = ", v1, " or v2 =", v2
else:
    print "  The particles are initialized with one initial speed: v1 = ", v1

specular = args.specular
if(specular):
    print "  Specular reflection is used"
    teleport = False
    boundaryString = ' Specular'

noTeleport = args.noTeleport
if(noTeleport):
    if(specular == False):
        print "  Neither teleporation nor specular reflection are NOT used"
        teleport = False
        boundaryString = ' Escape'
else:
    if(specular == False):
        print "  Teleporation is used"
        teleport = True
        boundaryString = ' Teleport'

if(seedFixed):
    random.seed(3)  # fixes the sequence of random sees for repeatability of numerical results
    print "  Fixed seed for random number generator"
else:
    random.seed()   # random seed is set according to the system clock, numerical results will be non-repeatable
    print "  Fixed seed for random number generator"

noPaths = args.noPaths
if(noPaths):
    print "  The particle paths will not be plotted"
else:
    print "  The particles paths will be plotted"

fixCMVelocity = args.fixCMVelocity
if(useDataInitializationFile == False):
    if(fixCMVelocity):
        print "  The initialization step insures that there is zero CM velocity at the start"
    else:
        print "  The initialization step does not insure that there is zero CM velocity at the start"

if(rFactor != 0.0):
    print " ", changeKEString
else:
    print "  There is no change requested for the kinetic energy"

boundaryMomentumChangeX = 0.0
boundaryMomentumChangeY = 0.0

noAnimation = args.noAnimation
nFrames = args.nFrames
snapShotFactor = args.snapShotFactor
nSnapShots = int(maxT/(deltaT*snapShotFactor))
if(nFrames > nSnapShots):
    nFrames = nSnapShots

framesToTimesFactor = maxT/nFrames  # the conversion from the animation frame number to the actual time of the frame

nInterval = args.nInterval
doMovie = args.doMovie  # retrieve the choice to make the mp4 movie

if(noAnimation):
    print "  No animation will produced. Graphs will be produced instead"
else:
    print "  Animation will have ", nFrames, " frames using ", nInterval, " milliseconds between frames, comprising ", nSnapShots, " snapshots of the total time with ",deltaT*snapShotFactor, " time units between snapshots"
    print "  The conversion factor from frame number to the simulation elapsed time is %.3e time units per frame" % framesToTimesFactor

#
# data lists for each particle: particleTime[timeStep][kParticle], previousPositionX[timeStep][kParticle]
#                               previousPositionY[timeStep][kParticle], currentPositionX[timeStep][kParticle]
#                               currentPositionY[timeStep][kParticle], currentSpeedX[timeStep][kParticle]
#                               currentSpeedY[timeStep][kParticle], currentForceX[timeStep][kParticle], currentForceY[timeStep][kParticle]
#
particleTime = [[0.0 for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]
previousPositionX = [[0.0 for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]
previousPositionY = [[0.0 for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]
currentPositionX = [[0.0 for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]
currentPositionY = [[0.0 for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]
currentSpeedX = [[0.0 for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]
currentSpeedY = [[0.0 for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]
currentForceX = [[0.0 for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]
currentForceY = [[0.0 for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]

#
# data lists for global (particle summed) information: potentialEnergy[timeStep], kineticEnergy[timeStep], averageSpeed[timeStep]
#                                                      nBoundary[timeStep], netForceX[timeStep], netForceY[timeStep]
#
potentialEnergy = [0.0 for timeStep in range(timeStepsPlusOne)]
kineticEnergy = [0.0 for timeStep in range(timeStepsPlusOne)]
totalEnergy = [0.0 for timeStep in range(timeStepsPlusOne)]

averageSpeed = [0.0 for timeStep in range(timeStepsPlusOne)]
nBoundary = [0 for timeStep in range(timeStepsPlusOne)]
netForceX = [0.0 for timeStep in range(timeStepsPlusOne)]
netForceY = [0.0 for timeStep in range(timeStepsPlusOne)]

#
# Tag for a particle crossing a boundary
#
particleCrossedBoundary = [False for kParticles in range(nParticles)]

maxEnergy = -1.0e+37
minEnergy = +1.0e+37

def repositionForPlotting(xOffset, yOffset):
    #
    # Shift positions (not effective unless xOffSet and yOffSet are non-zero values)
    # This reposition function has no effect on the motions. It just re-centers the three sets of position plots after the end of the calculation
    #
    if(xOffset == 0.0 and yOffset == 0.0):
        return
    
    for timeStep in range(timeStepsPlusOne):
        for kParticle in range(nParticles):
            xNew = currentPositionX[timeStep][kParticle] + xOffset
            yNew = currentPositionY[timeStep][kParticle] + yOffset
            if(xNew < xMin):
                xNew = xMax - (xMin - xNew)
            if(xNew > xMax):
                xNew = xMin + (xNew - xMax)
            if(yNew < yMin):
                yNew = yMax - (yMin - yNew)
            if(yNew > yMax):
                yNew = yMin + (yNew - yMax)
            currentPositionX[timeStep][kParticle] = xNew
            currentPositionY[timeStep][kParticle] = yNew
    return

def initMolecularDynanics():
    global squareRootParticles, xMax, yMax, xMin, yMin
    
    #
    # initial positions will be set in lattice (boxSize - extraBoxSize/2.0) x (boxSize - extraBoxSize/2.0)
    # this lattice will be increased to boxSize x boxSize after the initial positions are estabished
    #
    
    if(useDataInitializationFile):
        fileInput = open('initializationData.txt', 'r')
        print "\n  Input file  initializationData.txt  was successfully opened"
        headerLineString = fileInput.readline()
        headerLineStringNoEOL = headerLineString.rstrip()  # removes the end of line character
        if(verbose):
            print "  ", headerLineStringNoEOL
        nParticlesKENameBoxSizeMaxTDeltaTV1V2String = (fileInput.readline()).rstrip()
        nParticlesKENameBoxSizeMaxTDeltaTV1V2StringSplit = nParticlesKENameBoxSizeMaxTDeltaTV1V2String.split(' ')
        lengthTest = len(nParticlesKENameBoxSizeMaxTDeltaTV1V2StringSplit )
        if(lengthTest != 8):
            print "\n [initMolecularDynamics] The number of data fields in the second input line is incorrect at ", lengthTest
            exit()
        particleCount = int(nParticlesKENameBoxSizeMaxTDeltaTV1V2StringSplit[1])
        if(particleCount != nParticles):
            print "\n [initMolecularDynamics] The number of particle in the data file ", particleCount, " does not match the command line partcile request ", nParticles
            exit(1)
        dataKEText = nParticlesKENameBoxSizeMaxTDeltaTV1V2StringSplit[2]
        dataBoxSize = float(nParticlesKENameBoxSizeMaxTDeltaTV1V2StringSplit[3])
        dataMaxT = float(nParticlesKENameBoxSizeMaxTDeltaTV1V2StringSplit[4])
        dataDeltaT = float(nParticlesKENameBoxSizeMaxTDeltaTV1V2StringSplit[5])
        dataV1 = float(nParticlesKENameBoxSizeMaxTDeltaTV1V2StringSplit[6])
        dataV2 = float(nParticlesKENameBoxSizeMaxTDeltaTV1V2StringSplit[7])
        print "  The input data file was produced for ", particleCount, " particles, with a box size ", dataBoxSize, " and a KE change option ", dataKEText
        print "  The input data file had a maximum time ", dataMaxT, " in time steps ", dataDeltaT
        if(dataV1 == dataV2):
            print "  The input data file used a single initial speed v1 ", dataV1
        else:
            print "  The input data file used two initial speeds v1 ", dataV1, " and v2 ", dataV2

        vxSum = 0.0
        vySum = 0.0
        xSum = 0.0
        ySum = 0.0
        for kParticle in range(nParticles):
            dataLineString = (fileInput.readline()).rstrip()
            dataLineStringSplit = dataLineString.split(' ')
            lengthTest = len(dataLineStringSplit)
            if(lengthTest != 4):
                print "\n [initMolecularDynamics] The number of data fields = ", lengthTest, " is incorrect"
                exit()
            currentPositionX[0][kParticle] = float(dataLineStringSplit[0])
            currentPositionY[0][kParticle] = float(dataLineStringSplit[1])
            currentSpeedX[0][kParticle] = float(dataLineStringSplit[2])
            currentSpeedY[0][kParticle] = float(dataLineStringSplit[3])
            previousPositionX[0][kParticle] = currentPositionX[0][kParticle] - dataDeltaT*currentSpeedX[0][kParticle]
            previousPositionY[0][kParticle] = currentPositionY[0][kParticle] - dataDeltaT*currentSpeedY[0][kParticle]
            xSum += currentPositionX[0][kParticle]
            ySum += currentPositionY[0][kParticle]
            vxSum += currentSpeedX[0][kParticle]
            vySum += currentSpeedY[0][kParticle]

        fileInput.close()
        xCM = xSum/float(particleCount)
        yCM = ySum/float(particleCount)
        print "  The input data file has been processed.  The x and y speeds sums are ", vxSum, vySum, " and the CM is at ", xCM, yCM
   
        return  #  Completed processing of the data initialization file
    
    xDelta = (xMax - xMin)/float(squareRootParticles)
    yDelta = (yMax-yMin)/float(squareRootParticles)
   
    if(xDelta < 2.0*deltaR*sigma or yDelta < 2.0*deltaR*sigma):
        print "\n Initial separation for SIGMA = ", sigma, " is too small"
        print "  deltaR = ", deltaR
        print "  xDelta = ", xDelta, ",  yDelta = ", yDelta, "\n"
        exit(1)

    kParticle = 0
    xValue = 0.5*xDelta + xMin
    
    #
    # Note: to obtain a 50-50 random distribution the call random.randrange(-1,3,2) returns either -1 or +1 randomly
    #
    for ix in range(squareRootParticles):
        yValue = 0.5*yDelta + yMin
        for iy in range(squareRootParticles):
            particleTime[0][kParticle] = 0.0
            
            shift = deltaR*sigma*random.uniform(-1.0, 1.0)   # returns a value [-1.0, 1.0) uniformly distributed including -1 but not +1
            currentPositionX[0][kParticle] = xValue + shift

            shift = deltaR*sigma*random.uniform(-1.0, 1.0)
            currentPositionY[0][kParticle] = yValue + shift

            azimuthalAngle = 2.0*np.pi*random.uniform(0.0, 1.0)
            if(useTwoSpeeds and random.randrange(-1,3,2) < 0):
                thisV = v2
            else:
                thisV = v1
            currentSpeedX[0][kParticle] = thisV*np.cos(azimuthalAngle)       # random binary choice
            currentSpeedY[0][kParticle] = thisV*np.sin(azimuthalAngle)       # same random binary choice
               
            previousPositionX[0][kParticle] = currentPositionX[0][kParticle] - currentSpeedX[0][kParticle]*deltaT
            previousPositionY[0][kParticle] = currentPositionY[0][kParticle] - currentSpeedY[0][kParticle]*deltaT
                
            yValue += yDelta
            kParticle += 1
        
        xValue += xDelta
    
    #
    # Check if the CM velocity should set to zero
    #
    if(fixCMVelocity):
        vxSum = 0.0
        vySum = 0.0
        for kParticle in range(nParticlesMinusOne):
            vxSum += currentSpeedX[0][kParticle]
            vySum += currentSpeedY[0][kParticle]

        #
        # The initial speed of the last particle is reset
        #
        currentSpeedX[0][nParticlesMinusOne] = -vxSum
        currentSpeedY[0][nParticlesMinusOne] = -vySum

        #
        # The previous position of the last particle must also be fixed
        #
        previousPositionX[0][nParticlesMinusOne] = currentPositionX[0][nParticlesMinusOne] - currentSpeedX[0][nParticlesMinusOne]*deltaT
        previousPositionY[0][nParticlesMinusOne] = currentPositionY[0][nParticlesMinusOne] - currentSpeedY[0][nParticlesMinusOne]*deltaT

    return

def getTwoParticleDistance(iParticle, jParticle, timeIndex):

    if(iParticle == jParticle):
        print "\n  ***Error, two particle separation has i index at ", iParticle, " the same as j index at ", jParticle
        exit()
        
    ijXDist = currentPositionX[timeIndex][jParticle] - currentPositionX[timeIndex][iParticle]
    ijYDist = currentPositionY[timeIndex][jParticle] - currentPositionY[timeIndex][iParticle]

    if(teleport or specular):
        if(ijXDist > maxDiffX): ijXDist = ijXDist - xMaxMinusMin
        if(ijXDist < -maxDiffX): ijXDist = ijXDist + xMaxMinusMin
        if(ijYDist > maxDiffY): ijYDist = ijYDist - yMaxMinusMin
        if(ijYDist < -maxDiffY): ijYDist = ijYDist + yMaxMinusMin
        
    ijRDist = np.sqrt(ijXDist*ijXDist + ijYDist*ijYDist)
    if(ijRDist <= 0.0):
        print "\n  Error in separation distance ", ijRDist, " for time index ", timeIndex
        print "  iParticle ", iParticle, " x, y", currentPositionX[timeIndex][iParticle], currentPositionY[timeIndex][iParticle]
        print "  jParticle ", jParticle, " x, y", currentPositionX[timeIndex][jParticle], currentPositionY[timeIndex][jParticle]
        print "  ijXDist ", ijXDist, ",  ijYDist ", ijYDist
        exit()

    cosine = ijXDist/ijRDist
    sine = ijYDist/ijRDist
    if(sine == 0.0):
        print "\n [getTwoParticleDistance] 0 angle ", ijXDist, ijYDist, " time index ", timeIndex
        print " i coordinates " , iParticle, currentPositionX[timeIndex][iParticle], currentPositionY[timeIndex][iParticle]
        print " j coordinates " , jParticle, currentPositionX[timeIndex][jParticle], currentPositionY[timeIndex][jParticle]
        exit()
    return ijRDist, cosine, sine

def calculateSeparation(timeIndex):
    
    minSeparation = +1.0e+37
    maxSeparation = -1.0e+37
    iSaveMin = 0
    jSaveMin = 0
    iSaveMax = 0
    jSaveMax = 0

    #
    # Velocity of the last particle at this time step
    #
    vx = currentSpeedX[timeIndex][nParticlesMinusOne]
    vy = currentSpeedY[timeIndex][nParticlesMinusOne]
    vSum = np.sqrt(vx*vx + vy*vy)

    averageSeparation = 0

    #
    # Loop over distinct pairs of particles ij
    #
    for iParticle in range(nParticlesMinusOne):
        vx = currentSpeedX[timeIndex][iParticle]
        vy = currentSpeedY[timeIndex][iParticle]
        vSum += np.sqrt(vx*vx + vy*vy)

        ixDist = currentPositionX[timeIndex][iParticle]
        iyDist = currentPositionY[timeIndex][iParticle]

        jParticle = iParticle + 1
        while jParticle < nParticles:
            jxDist = currentPositionX[timeIndex][jParticle]
            jyDist = currentPositionY[timeIndex][jParticle]

            separation = np.sqrt((jxDist-ixDist)*(jxDist-ixDist) + (jyDist-iyDist)*(jyDist-iyDist))
            averageSeparation += separation

            if(separation < minSeparation):
                minSeparation = separation
                iSaveMin = iParticle
                jSaveMin = jParticle

            if(separation > maxSeparation):
                maxSeparation = separation
                iSaveMax = iParticle
                jSaveMax = jParticle
                    
            jParticle += 1

    print  "\n  min separation = " , minSeparation , " for iPart = " , iSaveMin , " and jPart = " , jSaveMin
    print  "  max separation = " , maxSeparation , " for iPart = " , iSaveMax , " and jPart = " , jSaveMax
    print  "  average speed = " , vSum/nParticles
    print  "  average separation = " , averageSeparation/nDistinctPairs

    return

def printParticleData(timeIndex):
    
    print "\n  Particle information at time ", timeIndex
    for kParticle in range(nParticles):
        x = currentPositionX[timeIndex][kParticle]
        y = currentPositionY[timeIndex][kParticle]
        vx = currentSpeedX[timeIndex][kParticle]
        vy = currentSpeedY[timeIndex][kParticle]
        
        printString  = '  Particle %4d has (x,y) = (%6.3f,%6.3f) and (vx,vy) = (%6.3f,%6.3f)' % (kParticle, x, y, vx, vy)
        print printString

    print " "
    return

def changeKineticEnergy(timeStep):
    global rChangeFactor
    for kParticle in range(nParticles):
        xCurrent = currentPositionX[timeStep][kParticle]
        xPrevious = previousPositionX[timeStep][kParticle]
        previousPositionX[timeStep][kParticle] = xCurrent - rChangeFactor*(xCurrent - xPrevious)
        yCurrent = currentPositionY[timeStep][kParticle]
        yPrevious = previousPositionY[timeStep][kParticle]
        previousPositionY[timeStep][kParticle] = yCurrent - rChangeFactor*(yCurrent - yPrevious)

    return

def calculateEnergy(timeIndex):
    
    #
    # Function also checks for the maximum force components between any particle pair
    #
    thisKineticEnergy = 0.0
    thisPotentialEnergy = 0.0
    thisAverageSpeed = 0.0
    
    iSaveMin = 0
    jSaveMin = 0
    minDist = +1.0e+37
    vMax = -1.0e+37
    
    for kParticle in range(nParticles):
    
        vx = currentSpeedX[timeIndex][kParticle]
        vy = currentSpeedY[timeIndex][kParticle]

        speedSquared = vx*vx + vy*vy
        thisKineticEnergy += 0.5*speedSquared

        thisAverageSpeed += np.sqrt(speedSquared)
        
    #
    # Loop over distinct pairs of particles ij
    #
    for iParticle in range(nParticlesMinusOne):
        jParticle = iParticle + 1
        while jParticle < nParticles:
                
            ijRDist, cosine, sine = getTwoParticleDistance(iParticle, jParticle, timeIndex)
                
            #thisTerm = 0.0
            if(ijRDist < rCut):
                if(forceConstant > 0.0):
                    #thisTerm = forceConstant4*(pow(sigma/ijRDist, 12) - pow(sigma/ijRDist, 6))
                    thisPotentialEnergy += forceConstant4*(pow(sigma/ijRDist, 12) - pow(sigma/ijRDist, 6))
                        
                if(ijRDist < minDist):
                    minDist = ijRDist
                    iSaveMin = iParticle
                    jSaveMin = jParticle
                    if(forceConstant > 0.0):
                        vMax = forceConstant4*(1.0/pow(minDist,12) -1.0/pow(minDist,6))   # not used ??
            
            #print " i,j ", iParticle, jParticle, " ijRDist ", ijRDist, ",  this term", thisTerm, ",  total ", thisPotentialEnergy
            jParticle += 1

    potentialEnergy[timeIndex] = thisPotentialEnergy
    kineticEnergy[timeIndex] = thisKineticEnergy
    totalEnergy[timeIndex] = thisPotentialEnergy + thisKineticEnergy
    averageSpeed[timeIndex] = thisAverageSpeed/nParticles

    return

def calculateSingleForce(iParticle, timeIndex):

    xForce = 0.0
    yForce = 0.0
    
    for jParticle in range(nParticles):
    
        if(jParticle != iParticle):
    
            ijRDist,cosine, sine = getTwoParticleDistance(iParticle, jParticle, timeIndex)

            forceTotal = 0
            if(ijRDist < rCut and forceConstant > 0.0):
                forceTotal = -forceConstant24*(2.0*pow(sigma/ijRDist, 13) - pow(sigma/ijRDist, 7))
        
            xForce += forceTotal*cosine
            yForce += forceTotal*sine

            #
            # Safety warnings
            #
            if(abs(xForce*deltaTSquared) > deltaTest or abs(yForce*deltaTSquared) > deltaTest):
                print "\n *** [calculateSingleForce] The deltaT time step choice ", deltaT, "is too large.  The program will not complete.***"
                print "  xForce*deltaTSquared", xForce*deltaTSquared, ",  yForce*deltaTSquared ", yForce*deltaTSquared, " with maximum allowed at ", deltaTest
                print "  iParticle ", iParticle, ", jParticle ", jParticle, ",  time step ", timeIndex
                print "  ijRDist ", ijRDist, " with deltaTSquared ", deltaTSquared
                print "  forceTotal ", forceTotal, ",  cosine ", cosine, ",  sine ", sine
                exit()
    
            if(ijRDist < rMinTest):
                print "\n  The particles have become too close at ", ijRDist
                print "  The minimum test distance is coded at ", rMinTest
                print "  The program won't complete"
                print "  xForce ", xForce, ",  yForce ", yForce
                exit()
    
    currentForceX[timeIndex][iParticle] = xForce
    currentForceY[timeIndex][iParticle] = yForce
    
    return xForce, yForce

def calculateAllForces(timeIndex):

    xyForcesThisTimeStep = [calculateSingleForce(iParticle, timeIndex) for iParticle in range(nParticles)]
    xyForceParticle = [xyForcesThisTimeStep[iParticle] for iParticle in range(nParticles)]
    netForceX = 0.0
    netForceY = 0.0
    averageForce = 0.0
    for iParticle in range(nParticles):
        thisForce = xyForceParticle[iParticle]
        thisForceX = thisForce[0]
        thisForceY = thisForce[1]
        currentForceX[timeIndex][iParticle] = thisForceX
        currentForceY[timeIndex][iParticle] = thisForceY
        
        netForceX += thisForce[0]
        netForceY += thisForce[1]
        averageForce += np.sqrt(thisForceX*thisForceX + thisForceY*thisForceY)

    return netForceX, netForceY, averageForce/float(nParticles)

def checkMomentum(timeStep):
    global boundaryMomentumChangeX, boundaryMomentumChangeY, xCenter, yCenter
    #
    # Since all the particles have the same mass then summing the components of velocity is equivalent to obtaining the net momentum in the X and Y directions
    #
    summedMomentumX = boundaryMomentumChangeX
    summedMomentumY = boundaryMomentumChangeY
    #
    # The angular momentum is calculated with respect to the center of the confining box
    #
    summedAngularMomentum = 0.0
    for kParticle in range(nParticles):
        vx = currentSpeedX[timeStep][kParticle]
        vy = currentSpeedY[timeStep][kParticle]
        summedMomentumX += vx
        summedMomentumY += vy
        xFromCenter = currentPositionX[timeStep][kParticle] - xCenter
        yFromCenter = currentPositionY[timeStep][kParticle] - yCenter
        summedAngularMomentum += xFromCenter*vy - yFromCenter*vx

    return summedMomentumX, summedMomentumY, summedAngularMomentum

def checkTeleportBoundary(dist, dMin, dMax, kParticle):
    global nBoundaryCount
    
    if(noTeleport):
        if(particleCrossedBoundary[kParticle] == False and (dist < dMin or dist > dMax)):
            particleCrossedBoundary[kParticle] = True
            nBoundaryCount += 1
        return dist

    distNew = dist

    if(dist < dMin):
        distNew = dMax - (dMin - dist)
        nBoundaryCount += 1

    if(dist > dMax):
        distNew = dMin + (dist - dMax)
        nBoundaryCount += 1

    if(distNew < xMin or distNew > xMax):
        print "\n  ***Error in checkBoundary function"
        print "  distNew = ", distNew, ",  dist = ", dist
        exit(1)

    return distNew

def checkSpecularBoundary(xNew, yNew, vxNew, vyNew):
    global nBoundaryCount
    
    xFix = xNew
    yFix = yNew
    vxFix = vxNew
    vyFix = vyNew
    
    if(xNew >= xMin and xNew <= xMax and yNew >= yMin and yNew <= yMax):
        # new positions are all within the boundary, no changes are made
        return xFix, yFix, vxFix, vyFix

    nBoundaryCount += 1

    if(xNew < xMin and xNew <= xMax and yNew >= yMin and yNew <= yMax):
        # new value of x is less than minimum number
        xFix = xMin + xMin - xNew
        vxFix = -vxNew
        return xFix, yFix, vxFix, vyFix

    if(xNew >= xMin and xNew > xMax and yNew >= yMin and yNew <= yMax):
        # new value of x is greater than minimum number
        xFix = xMax + xMax - xNew
        vxFix = -vxNew
        return xFix, yFix, vxFix, vyFix

    if(xNew >= xMin and xNew <= xMax and yNew < yMin and yNew <= yMax):
        # new value of y is less than minimum number
        yFix = yMin + yMin - yNew
        vyFix = -vyNew
        return xFix, yFix, vxFix, vyFix

    if(xNew >= xMin and xNew <= xMax and yNew >= yMin and yNew > yMax):
        # new value of y is greater than maximum number
        yFix = yMax + yMax - yNew
        vyFix = -vyNew
        return xFix, yFix, vxFix, vyFix

    if(xNew < xMin and xNew <= xMax and yNew < yMin and yNew <= yMax):
        # new value of x is less than minimum number, new value of y is less than minium
        xFix = xMin + xMin - xNew
        yFix = yMin + yMin - yNew
        vxFix = -vxNew
        vyFix = -vyNew
        return xFix, yFix, vxFix, vyFix

    if(xNew < xMin and xNew <= xMax and yNew >= yMin and yNew > yMax):
        # new value of x is less than minimum number, new value of y is greater than maximum
        xFix = xMin + xMin - xNew
        yFix = yMax + yMax - yNew
        vxFix = -vxNew
        vyFix = -vyNew
        return xFix, yFix, vxFix, vyFix

    if(xNew >= xMin and xNew > xMax and yNew < yMin and yNew <= yMax):
        # new value of x is greater than maximum number, new value of y is less than minimum
        xFix = xMax + xMax - xNew
        yFix = yMin + yMin - yNew
        vxFix = -vxNew
        vyFix = -vyNew
        return xFix, yFix, vxFix, vyFix

    if(xNew >- xMin and xNew > xMax and yNew >= yMin and yNew > yMax):
        # new value of x is greater than maximum number, new value of y is greater than maximum
        xFix = xMax + xMax - xNew
        yFix = yMax + yMax - yNew
        vxFix = -vxNew
        vyFix = -vyNew
        return xFix, yFix, vxFix, vyFix

def updateMolecularDynamics(timeIndex):
    global boundaryMomentumChangeX, boundaryMomentumChangeY
    
    timeIndexNew = timeIndex+1
    particleTime[timeIndexNew] = timeIndexNew*deltaT
    
    calculateAllForces(timeIndex)
    
    for iParticle in range(nParticles):
        
        xPrevious = previousPositionX[timeIndex][iParticle]
        yPrevious = previousPositionY[timeIndex][iParticle]
        
        xCurrent = currentPositionX[timeIndex][iParticle]
        yCurrent = currentPositionY[timeIndex][iParticle]
    
        xNew = 2.0*xCurrent - xPrevious + currentForceX[timeIndex][iParticle]*deltaTSquared
        yNew = 2.0*yCurrent - yPrevious + currentForceY[timeIndex][iParticle]*deltaTSquared
    
        vxNew = (xNew - xPrevious)/twiceDeltaT
        vyNew = (yNew - yPrevious)/twiceDeltaT
        currentSpeedX[timeIndexNew][iParticle] = vxNew
        currentSpeedY[timeIndexNew][iParticle] = vyNew
    
        deltaX = abs(xNew - xCurrent)
        deltaY = abs(yNew - yCurrent)
        if(deltaX > deltaTest or deltaY > deltaTest):
            print "\n *** [updateMolecularDynamics] Change in position is too large: deltaX ", deltaX, " or deltaY ", deltaY, " with limit ", deltaTest
            print "   iParticle ", iParticle, " at time step ", timeIndexNew
            print "  Current positions: x,y ", xCurrent, yCurrent
            print "  New positions: x,y ", xNew, yNew
            print "  Previous positions x,y ", xPrevious, yPrevious
            print "  boundary count ", nBoundaryCount
            print "  ForceX ", currentForceX[timeIndex][iParticle], ",  ForceY ", currentForceY[timeIndex][iParticle]
            exit()

        xPrevious = xCurrent
        yPrevious = yCurrent
        
        if(specular):
            xFix, yFix, vxFix, vyFix = checkSpecularBoundary(xNew, yNew, vxNew, vyNew)
            #
            # these code lines work for specular reflection only
            # if a boundary is not being crossed, then the Fix and the New values are the same
            #
            if(xFix < xMin or xFix > xMax or yFix < yMin or yFix > yMax):
                print "\n [updateMolecularDynamics]  ***Bad position values returned from checkSpecularBoundary function***"
                exit()
            if(vxNew*vxNew != vxFix*vxFix or vyNew*vyNew != vyFix*vyFix):
                print "\n [updateMolecularDynamics]  ***Bad velocity values returned from checkSpecularBoundary function***"
                exit()
            
            if(vxFix != vxNew):
                boundaryMomentumChangeX += vxNew - vxFix
            
            if(vyFix != vyNew):
                boundaryMomentumChangeY += vyNew - vyFix

            currentPositionX[timeIndexNew][iParticle] = xFix
            currentPositionY[timeIndexNew][iParticle] = yFix

            currentSpeedX[timeIndexNew][iParticle] = vxFix
            currentSpeedY[timeIndexNew][iParticle] = vyFix
            
            #
            # Normally the previous position is set equal to the current position
            # However, if the particle has been reflected, then the previous position for the new time index needs to be corrected
            #

            if(xFix != xNew):
                if(xNew < xMin):
                    xPrevious = xMin - (xPrevious - xMin)
                if(xNew > xMax):
                    xPrevious = xMax + (xMax - xPrevious)
                
            if(yFix != yNew):
                if(yNew < yMin):
                    yPrevious = yMin - (yPrevious - yMin)
                if(yNew > yMax):
                    yPrevious = yMax + (yMax - yPrevious)

            previousPositionX[timeIndexNew][iParticle] = xPrevious
            previousPositionY[timeIndexNew][iParticle] = yPrevious
    
        else:
            #
            # these code lines work for both teleport and not teleporting
            #
            xChecked = checkTeleportBoundary(xNew, xMin, xMax, iParticle)
            yChecked = checkTeleportBoundary(yNew, yMin, yMax, iParticle)
            
            currentPositionX[timeIndexNew][iParticle] = xChecked
            currentPositionY[timeIndexNew][iParticle] = yChecked

            #
            # Normally the previous position is set equal to the current position
            # However, if the particle has been teleported, then the previous position for the new time index needs to be corrected
            #
        
            if(xChecked != xNew):
                #
                # Correcting for a teleportation in the x direction
                #
                xPrevious = xChecked - (xNew - xPrevious)
        
            if(yChecked != yNew):
                #
                # Correcting for a teleportation in the y direction
                #
                yPrevious = yChecked - (yNew - yPrevious)
        
            previousPositionX[timeIndexNew][iParticle] = xPrevious
            previousPositionY[timeIndexNew][iParticle] = yPrevious

    return

def maxwellBoltzmann(velocity, temperature):
    if(temperature > 0.0):
        return (velocity/temperature)*np.exp(-(velocity*velocity)/(2.0*temperature))
    else:
        print "\n [maxwellBoltzmann] ***Error with temperature not positive ", temperature, " with velocity ", velocity, "***"
        exit(1)

def fitMaxwellBoltzmann(velocityArray, velocityProbabilityArray, velocityProbabilityErrorArray):
    initialTemperatureGuess = np.array([0.5])
    fitResultsOneParameter = optimization.curve_fit(maxwellBoltzmann, velocityArray, velocityProbabilityArray, initialTemperatureGuess, velocityProbabilityErrorArray)
    fitTemperatureArray = fitResultsOneParameter[0]
    fitTemperature = fitTemperatureArray[0]
   
    chiSquareSum = 0.0
    nData = len(velocityArray)
    predictedProbability = [0.0 for iData in range(nData)]
    if(verbose):
        print "\n  Fit for temperature ", 120.0*fitTemperature     # converting temperature to Kelvin
    for iData in range(nData):
        thisPrediction = maxwellBoltzmann(velocityArray[iData], fitTemperature)
        predictedProbability[iData] = thisPrediction
        difference = thisPrediction - velocityProbabilityArray[iData]
        error = velocityProbabilityErrorArray[iData]
        chiSquare = difference*difference/(error*error)
        chiSquareSum += chiSquare
        if(verbose):
            printString = '  %.3e %.3e %.3e %.3e %.3e'  % (velocityArray[iData], velocityProbabilityArray[iData], thisPrediction, error, chiSquare)
            print printString

    reducedChiSquare = chiSquareSum/float(nData - 1)
    if(verbose):
        print "  Reduced chiSquare", reducedChiSquare

    return fitTemperature, reducedChiSquare, predictedProbability

def graphsOutput():
    #
    # Function to produce static plots
    #
    global changeKEString, changeKEStringV0, keNameText, textFileNameString
    
    figGraphs = plt.figure(1)          # start a figure for energy plot
    axGraphs = figGraphs.add_subplot(111)
    
    timePlot = [deltaT*timeStep for timeStep in range(timeStepsPlusOne)]
    plt.plot(timePlot, kineticEnergy, 'r', label='Kinetic Energy')
    plt.plot(timePlot, potentialEnergy, 'g', label='Potential Energy')
    plt.plot(timePlot, totalEnergy, 'b', label='Total Energy')

    maxKE = max(kineticEnergy)
    minKE = min(kineticEnergy)
    maxPE = max(potentialEnergy)
    minPE = min(potentialEnergy)
    minTotal = min(totalEnergy)
    maxTotal = max(totalEnergy)

    yMax = 1.4*max(maxKE, maxPE, maxTotal)
    yMin = min(minKE, minPE, minTotal)
    plt.ylim(yMin, yMax)
    titleString = 'Energy for %d Particles in Box Size %.2f, Boundary%s' % (nParticles, boxSize, boundaryString)
    plt.title(titleString)
    plt.xlabel('Time (units of 1.8 ps)')                             # add axis labels
    plt.ylabel('Energy ')
    plt.grid(True)
    plt.legend(loc=2)
    energy_text = axGraphs.text(0.38, 0.88, energyString, transform=axGraphs.transAxes, color='blue')
    timeStepString = 'Time step %.3e' % deltaT
    timeStep_text = axGraphs.text(0.38, 0.81, timeStepString, transform=axGraphs.transAxes, color='blue')
    plotNameString = 'Energy-Particles%d-%s-BoxSize%d-T%dMilli-DeltaT%dMilli.pdf' % (nParticles, keNameText, boxSize, 1000*maxT, 1000*deltaT)
    figGraphs.savefig(plotNameString)
    plt.close(figGraphs)
    print "  Energy plot file name", plotNameString
    
    figGraphs = plt.figure(1)          # start a figure for Temperature plot
    axGraphs = figGraphs.add_subplot(111)

    temperatureFactor = 120.0/float(nParticles)  # temperature = 120 times the average kinetic energy
    temperature = [temperatureFactor*kineticEnergy[timeStep] for timeStep in range(timeStepsPlusOne)]
    plt.plot(timePlot, temperature, 'r', label='Kinetic Energy')
    
    yMax = 1.1*max(temperature)
    yMin = 0.9*min(temperature)
    plt.ylim(yMin, yMax)
    titleString = 'Temperature for %d Particles in Box Size %.2f, Boundary%s' % (nParticles, boxSize, boundaryString)
    plt.title(titleString)
    plt.xlabel('Time (units of 1.8 ps)')                             # add axis labels
    plt.ylabel('Temperature (Kelvin) ')
    plt.grid(True)

    argonMeltY = [83.8, 83.8]
    argonMeltX = [timePlot[0], timePlot[timeSteps]]
    plt.plot(argonMeltX, argonMeltY, 'b')
    plotNameString = 'Temperature-Particles%d-%s-BoxSize%d-T%dMilli-DeltaT%dMilli.pdf' % (nParticles, keNameText, boxSize, 1000*maxT, 1000*deltaT)
    figGraphs.savefig(plotNameString)
    plt.close(figGraphs)
    print "  Temperature plot file name", plotNameString
    
    currentVelocity = [[np.sqrt(currentSpeedX[timeStep][kParticle]*currentSpeedX[timeStep][kParticle] + currentSpeedY[timeStep][kParticle]*currentSpeedY[timeStep][kParticle]) for kParticle in range(nParticles)] for timeStep in range(timeStepsPlusOne)]
    vMaximum = 0.0
    tMaximum = 0.0
    particleMaximum = -1
    for timeStep in range(timeStepsPlusOne):
        for kParticle in range(nParticles):
            thisVelocity = currentVelocity[timeStep][kParticle]
            if(thisVelocity > vMaximum):
                vMaximum = thisVelocity
                tMaximum = timeStep*deltaT
                particleMaximum = kParticle

    finalVelocity = [0.0 for kParticle in range(nParticles)]
    for kParticle in range(nParticles):
        finalVelocity[kParticle]= currentVelocity[timeSteps][kParticle]

    finalVelocityAverage = np.mean(finalVelocity)
    finalVelocityRMS = np.std(finalVelocity)
                     
    maximumSpeedString = '\n  The maximum speed was %.3f at time %.3f for particle number%3d' % (vMaximum, tMaximum, particleMaximum)
    print maximumSpeedString

    finalVelocityString = '  The average final speed was %.3f with an RMS %.3f' % (finalVelocityAverage, finalVelocityRMS)
    print finalVelocityString

    oneThirdTimeSteps = int(timeStepsPlusOne/3)
    twiceOneThirdTimeSteps = 2*oneThirdTimeSteps
    speedsFirstThird = [0.0 for timeStep in range(oneThirdTimeSteps*nParticles)]
    speedsSecondThird = [0.0 for timeStep in range(oneThirdTimeSteps*nParticles)]
    speedsLastThird = [0.0 for timeStep in range(oneThirdTimeSteps*nParticles)]
    index = 0
    for kParticle in range(nParticles):
        thisParticleSpeedFirstThird = 0.0
        thisParticleSpeedSecondThird = 0.0
        thisParticleSpeedLastThird = 0.0
        for timeStep in range(oneThirdTimeSteps):
            speedsFirstThird[index] = currentVelocity[timeStep][kParticle]
            speedsSecondThird[index] = currentVelocity[timeStep + oneThirdTimeSteps][kParticle]
            speedsLastThird[index] = currentVelocity[timeStep + twiceOneThirdTimeSteps][kParticle]
            index += 1

    figGraphs = plt.figure(1)          # start a figure for velocity plots
    axGraphs = figGraphs.add_subplot(111)

    nVelocityBins = 20
    countsFirst, bin_edgesFirst = np.histogram(speedsFirstThird, nVelocityBins)
    countsSecond, bin_edgesSecond = np.histogram(speedsSecondThird, nVelocityBins)
    countsLast, bin_edgesLast = np.histogram(speedsLastThird, nVelocityBins)
    nDenominator = float(nParticles*oneThirdTimeSteps)
    vProbabilityFirst = [0.0 for iData in range(nVelocityBins)]
    vProbabilityFirstError = [0.0 for iData in range(nVelocityBins)]
    vProbabilitySecond = [0.0 for iData in range(nVelocityBins)]
    vProbabilitySecondError = [0.0 for iData in range(nVelocityBins)]
    vProbabilityLast = [0.0 for iData in range(nVelocityBins)]
    vProbabilityLastError = [0.0 for iData in range(nVelocityBins)]
    extraErrorFactor = 10.0   # accounting for the effect of highly correlated velocities for the same particle from one time step to the next
    for iData in range(nVelocityBins):
        if(countsFirst[iData] > 0):
            vProbabilityFirst[iData] = float(countsFirst[iData])/nDenominator
            vProbabilityFirstError[iData] = extraErrorFactor*vProbabilityFirst[iData]/np.sqrt(float(countsFirst[iData]))
        
        if(countsSecond[iData] > 0):
            vProbabilitySecond[iData] = float(countsSecond[iData])/nDenominator
            vProbabilitySecondError[iData] = extraErrorFactor*vProbabilitySecond[iData]/np.sqrt(float(countsSecond[iData]))

        if(countsLast[iData] > 0):
            vProbabilityLast[iData] = float(countsLast[iData])/nDenominator
            vProbabilityLastError[iData] = extraErrorFactor*vProbabilityLast[iData]/np.sqrt(float(countsLast[iData]))

    bin_centersFirst = (bin_edgesFirst[:-1] + bin_edgesFirst[1:])/2.
    firstThirdTimeString = 'Time interval %.1f to %.1f' % (0.0, deltaT*oneThirdTimeSteps)
    secondThirdTimeString = 'Time interval %.1f to %.1f' % (deltaT*oneThirdTimeSteps, deltaT*twiceOneThirdTimeSteps)
    lastThirdTimeString = 'Time interval %.1f to %.1f' % (deltaT*twiceOneThirdTimeSteps, maxT)

    binWidthFirst = bin_centersFirst[1] - bin_centersFirst[0]
    for iData in range(nVelocityBins):
        vProbabilityFirst[iData] = vProbabilityFirst[iData]/binWidthFirst
        vProbabilityFirstError[iData] =  vProbabilityFirstError[iData]/binWidthFirst
    plt.plot(bin_centersFirst, vProbabilityFirst, 'bo', markersize=3, label=firstThirdTimeString)

    bin_centersSecond = (bin_edgesSecond[:-1] + bin_edgesSecond[1:])/2.
    binWidthSecond = bin_centersSecond[1] - bin_centersSecond[0]
    for iData in range(nVelocityBins):
        vProbabilitySecond[iData] = vProbabilitySecond[iData]/binWidthSecond
        vProbabilitySecondError[iData] = vProbabilitySecondError[iData]/binWidthSecond
    plt.plot(bin_centersSecond, vProbabilitySecond, 'gD', markersize=3, label=secondThirdTimeString)

    bin_centersLast= (bin_edgesLast[:-1] + bin_edgesLast[1:])/2.
    binWidthLast= bin_centersLast[1] - bin_centersLast[0]
    for iData in range(nVelocityBins):
        vProbabilityLast[iData] =  vProbabilityLast[iData]/binWidthLast
        vProbabilityLastError[iData] = vProbabilityLastError[iData]/binWidthLast
    plt.plot(bin_centersLast, vProbabilityLast, 'rs', markersize=3, label=lastThirdTimeString)

    yMax1 = max(vProbabilityFirst)
    yMax2 = max(vProbabilitySecond)
    yMax3 = max(vProbabilityLast)
    yMax = 1.4*max(yMax1, yMax2, yMax3)
    
    firstTemperature, firstReducedChiSquare, firstPredictedProbability = fitMaxwellBoltzmann(bin_centersFirst, vProbabilityFirst, vProbabilityFirstError)
    secondTemperature, secondReducedChiSquare, secondPredictedProbability = fitMaxwellBoltzmann(bin_centersSecond, vProbabilitySecond, vProbabilitySecondError)
    lastTemperature, lastReducedChiSquare, lastPredictedProbability = fitMaxwellBoltzmann(bin_centersLast, vProbabilityLast, vProbabilityLastError)

    firstTemperatureKelvin = 120.0*firstTemperature
    secondTemperatureKelvin = 120.0*secondTemperature
    lastTemperatureKelvin = 120.0*lastTemperature

    print '  The Kelvin temperatures from the fits to the three time intervals are %.2f, %.2f, and %.2f' % (firstTemperatureKelvin, secondTemperatureKelvin, lastTemperatureKelvin)

    firstThirdPredictedString = 'Fit temperature %.2f K' % firstTemperatureKelvin
    plt.plot(bin_centersFirst, firstPredictedProbability, 'b', label=firstThirdPredictedString)

    secondThirdPredictedString = 'Fit temperature %.2f K' % secondTemperatureKelvin
    plt.plot(bin_centersSecond, secondPredictedProbability, 'g', label=secondThirdPredictedString)

    lastThirdPredictedString = 'Fit temperature %.2f K' % lastTemperatureKelvin
    plt.plot(bin_centersLast, lastPredictedProbability, 'r', label=lastThirdPredictedString)
     
    #
    # maximum velocity bin to be plotted is hard-coded as 5.5
    # maximum height P(v) to be plotted is set at 0.8, unless the temperatures are too low
    #
    xMax = 5.5
    yMax = 0.8
    yMax1 = max(firstPredictedProbability)
    yMax2 = max(secondPredictedProbability)
    yMax3 - max(lastPredictedProbability)
    yMax123 = 1.1*max(yMax1, yMax2, yMax3)
    if(yMax123 > 0.8):
        yMax = yMax123
    
    plt.xlim(0.0, xMax)
    plt.ylim(0.0, yMax)

    keTextX = 2.2
    keTextY = 0.52*yMax
    plt.text(keTextX, keTextY, changeKEStringV0, color='m')
    
    if(useDataInitializationFile):
        plt.text(keTextX, 1.1*keTextY, 'A data initialization file was used', color='m')
    else:
        v1X = [v1, v1]
        v1Y = [0.0, yMax]
        firstInitialSpeedString = 'First initial speed %.2f' % v1
        plt.plot(v1X, v1Y, 'k:', label=firstInitialSpeedString)
        if(useTwoSpeeds and v1 != v2):
            v2X = [v2, v2]
            v2Y = [0.0, yMax]
            secondInitialSpeedString = 'Second initial speed %.2f' % v2
            plt.plot(v2X, v2Y, 'k:', label=secondInitialSpeedString)

    titleString = 'Speed Distributions for %d Particles in Box Size %.2f' % (nParticles, boxSize)
    plt.title(titleString)
    plt.xlabel('v')                             # add axis labels
    plt.ylabel('P(v)')
    plt.grid(True)
    plt.legend(loc=1)
    
    plotNameString = 'Velocity-Particles%d-%s-BoxSize%d-T%dMilli-DeltaT%dMilli-FirstV%dSecondV2%d.pdf' % (nParticles, keNameText, boxSize, 1000*maxT, 1000*deltaT, 1000*v1, 1000*v2)
    figGraphs.savefig(plotNameString)
    plt.close(figGraphs)
    print "  Velocity plot file name", plotNameString

    repositionForPlotting(xOffSetFinal, yOffSetFinal)

    figGraphs = plt.figure(1)          # start a figure for initial position plots during the middle one-third of the time duration
    axGraphs = figGraphs.add_subplot(111)

    if(useDataInitializationFile):
        dataInitializationFileString = 'A data initialization file was used'
    else:
        dataInitializationFileString = 'No data initialization file was used'

    timeStepSkip = int(timeStepsPlusOne/300)  # this should plot 100 points in each one/third time interval for each particle
    if(timeStepSkip < 1):
        timeStepSkip = 1
    titleString = 'Initial positions for %d Particles in Box Size %.2f' % (nParticles, boxSize)
    plt.title(titleString)
    plt.xlabel('x (sigma units)')                             # add axis labels
    plt.ylabel('y (sigma units)')
    endPositionTimeIndex = int(timeStepsPlusOne/3.0)
    radialMove = [0.0 for kParticle in range(nParticles)]
    for kParticle in range(nParticles):
        thisParticleMoveX = currentPositionX[endPositionTimeIndex][kParticle] - currentPositionX[0][kParticle]
        thisParticleMoveY = currentPositionY[endPositionTimeIndex][kParticle] - currentPositionY[0][kParticle]
        radialMove[kParticle] = np.sqrt(thisParticleMoveX*thisParticleMoveX + thisParticleMoveY*thisParticleMoveY)
        timeStep = 0
        xPlot = []
        yPlot = []
        while timeStep < endPositionTimeIndex:
            xCurrent = currentPositionX[timeStep][kParticle]
            yCurrent = currentPositionY[timeStep][kParticle]
            xPlot.append(xCurrent)
            yPlot.append(yCurrent)
            timeStep += timeStepSkip
        
        plt.scatter(xPlot,yPlot, c='k', marker='.', s=1)

    yMax = 1.4*boxSize
    plt.xlim(0.0, boxSize)
    plt.ylim(0.0, yMax)
    
    plt.grid(True)
    xPositionText = 0.03*boxSize
    yPositionText = 0.9*yMax
    startTime = 0.0
    endTime = endPositionTimeIndex*deltaT
    positionsText = 'Positions for initial times %.2f to %.2f with T %.2f K' % (startTime, endTime, firstTemperatureKelvin)
    plt.text(xPositionText, yPositionText, positionsText, color='m')
    yPositionText = 0.85*yMax
    meanMoveDistance = np.mean(radialMove)
    RMSMoveDistance = np.std(radialMove)
    moveText = 'Mean radial move distance during this time was %.2f with RMS %.2f' % (meanMoveDistance, RMSMoveDistance)
    plt.text(xPositionText, yPositionText, moveText, color='m')
    yPositionText = 0.80*yMax
    plt.text(xPositionText, yPositionText, changeKEString, color='m')
    yPositionText = 0.75*yMax
    plt.text(xPositionText, yPositionText, dataInitializationFileString, color='m')

    plotNameString = 'InitialPositions-Particles%d-%s-BoxSize%d-T%dMilli-DeltaT%dMilli-FirstV%dSecondV2%d.pdf' % (nParticles, keNameText, boxSize, 1000*maxT, 1000*deltaT, 1000*v1, 1000*v2)
    figGraphs.savefig(plotNameString)
    plt.close(figGraphs)
    print '  Positions plot file name %s; mean radial move distance %.2f with RMS %.2f' % (plotNameString, meanMoveDistance, RMSMoveDistance)

    figGraphs = plt.figure(1)          # start a figure for interim position plots during the middle one-third of the time duration
    axGraphs = figGraphs.add_subplot(111)
    
    titleString = 'Interim positions for %d Particles in Box Size %.2f' % (nParticles, boxSize)
    plt.title(titleString)
    plt.xlabel('x (sigma units)')                             # add axis labels
    plt.ylabel('y (sigma units)')
    startPositionTimeIndex = int(timeStepsPlusOne/3.0)
    endPositionTimeIndex = 2*startPositionTimeIndex
    radialMove = [0.0 for kParticle in range(nParticles)]
    for kParticle in range(nParticles):
        thisParticleMoveX = currentPositionX[endPositionTimeIndex][kParticle] - currentPositionX[startPositionTimeIndex][kParticle]
        thisParticleMoveY = currentPositionY[endPositionTimeIndex][kParticle] - currentPositionY[startPositionTimeIndex][kParticle]
        radialMove[kParticle] = np.sqrt(thisParticleMoveX*thisParticleMoveX + thisParticleMoveY*thisParticleMoveY)
        timeStep = startPositionTimeIndex
        xPlot = []
        yPlot = []
        while timeStep < endPositionTimeIndex:
            xCurrent = currentPositionX[timeStep][kParticle]
            yCurrent = currentPositionY[timeStep][kParticle]
            xPlot.append(xCurrent)
            yPlot.append(yCurrent)
            timeStep += timeStepSkip
        
        plt.scatter(xPlot,yPlot, c='k', marker='.', s=1)

    plt.xlim(0.0, boxSize)
    plt.ylim(0.0, yMax)
    
    plt.grid(True)
    yPositionText = 0.9*yMax
    startTime = startPositionTimeIndex*deltaT
    endTime = endPositionTimeIndex*deltaT
    positionsText = 'Positions for interim times %.2f to %.2f with T %.2f K' % (startTime, endTime, secondTemperatureKelvin)
    plt.text(xPositionText, yPositionText, positionsText, color='m')
    yPositionText = 0.85*yMax
    meanMoveDistance = np.mean(radialMove)
    RMSMoveDistance = np.std(radialMove)
    moveText = 'Mean radial move distance during this was time %.2f with RMS %.2f' % (meanMoveDistance, RMSMoveDistance)
    plt.text(xPositionText, yPositionText, moveText, color='m')
    yPositionText = 0.80*yMax
    plt.text(xPositionText, yPositionText, changeKEString, color='m')
    yPositionText = 0.75*yMax
    plt.text(xPositionText, yPositionText, dataInitializationFileString, color='m')

    plotNameString = 'InterimPositions-Particles%d-%s-BoxSize%d-T%dMilli-DeltaT%dMilli-FirstV%dSecondV2%d.pdf' % (nParticles, keNameText, boxSize, 1000*maxT, 1000*deltaT, 1000*v1, 1000*v2)
    figGraphs.savefig(plotNameString)
    plt.close(figGraphs)
    print '  Positions plot file name %s; mean radial move distance %.2f with RMS %.2f' % (plotNameString, meanMoveDistance, np.std(radialMove))

    figGraphs = plt.figure(1)          # start a figure for final position plots during the last one-third of the time duration
    axGraphs = figGraphs.add_subplot(111)
    
    titleString = 'Final positions for %d Particles in Box Size %.2f' % (nParticles, boxSize)
    plt.title(titleString)
    plt.xlabel('x (sigma units)')                             # add axis labels
    plt.ylabel('y (sigma units)')
    startPositionTimeIndex = int(2.0*timeStepsPlusOne/3.0)
    
    #
    # This triplet of particle numbers works for nParticles = 36
    #
    if(nParticles%2 == 0):
        jParticle = int(nParticles/2 + np.sqrt(nParticles)/2.0)
    else:
        jParticle = int(nParticles/2)
    
    nParticle = jParticle + int(np.sqrt(nParticles))
    mParticle = jParticle - 1
    xCentroid = [0.0, 0.0, 0.0]
    yCentroid = [0.0, 0,0, 0.0]

    radialMove = [0.0 for kParticle in range(nParticles)]
    for kParticle in range(nParticles):
        thisParticleMoveX = currentPositionX[timeSteps][kParticle] - currentPositionX[startPositionTimeIndex][kParticle]
        thisParticleMoveY = currentPositionY[timeSteps][kParticle] - currentPositionY[startPositionTimeIndex][kParticle]
        radialMove[kParticle] = np.sqrt(thisParticleMoveX*thisParticleMoveX + thisParticleMoveY*thisParticleMoveY)
        timeStep = startPositionTimeIndex
        xPlot = []
        yPlot = []
        iCentroid = -1
        if(kParticle == mParticle):
            iCentroid = 0
        if(kParticle == jParticle):
            iCentroid = 1
        if(kParticle == nParticle):
            iCentroid = 2
        nTimes = 0
        while timeStep < timeStepsPlusOne:
            xCurrent = currentPositionX[timeStep][kParticle]
            yCurrent = currentPositionY[timeStep][kParticle]
            xPlot.append(xCurrent)
            yPlot.append(yCurrent)
            if(iCentroid != -1):
                xCentroid[iCentroid] += xCurrent
                yCentroid[iCentroid] += yCurrent
                nTimes += 1
            timeStep += timeStepSkip
        
        if(iCentroid != -1):
            xCentroid[iCentroid] /= float(nTimes)
            yCentroid[iCentroid] /= float(nTimes)
        
        if(kParticle == jParticle):
            plt.scatter(xPlot,yPlot, c='r', marker='.', s=3)
        if(kParticle == nParticle):
            plt.scatter(xPlot,yPlot, c='g', marker='.', s=3)
        if(kParticle == mParticle):
            plt.scatter(xPlot,yPlot, c='b', marker='.', s=3)
        if(iCentroid == -1):
            plt.scatter(xPlot,yPlot, c='k', marker='.', s=1)
    
    plt.xlim(0.0, boxSize)
    plt.ylim(0.0, yMax)
    if(useDataInitializationFile == False):   # Plot triangle only if initialization is not done from an external file
        topLineX = [xCentroid[1], xCentroid[2]]
        topLineY = [yCentroid[1], yCentroid[2]]
        leftLineX = [xCentroid[0], xCentroid[1]]
        leftLineY = [yCentroid[0], yCentroid[1]]
        rightLineX = [xCentroid[0], xCentroid[2]]
        rightLineY = [yCentroid[0], yCentroid[2]]
    
        plt.plot(topLineX, topLineY, c='m')
        plt.plot(leftLineX, leftLineY, c='m')
        plt.plot(rightLineX, rightLineY, c='m')
    
    plt.grid(True)
    yPositionText = 0.9*yMax
    startTime = startPositionTimeIndex*deltaT
    endTime = timeStepsPlusOne*deltaT
    positionsText = 'Positions for final times %.2f to %.2f with T = %.2f K' % (startTime, endTime, lastTemperatureKelvin)
    plt.text(xPositionText, yPositionText, positionsText, color='m')
    yPositionText = 0.85*yMax
    meanMoveDistance = np.mean(radialMove)
    RMSMoveDistance = np.std(radialMove)
    moveText = 'Mean radial move distance during this time was %.2f with RMS %.2f' % (meanMoveDistance, RMSMoveDistance)
    plt.text(xPositionText, yPositionText, moveText, color='m')
    yPositionText = 0.80*yMax
    plt.text(xPositionText, yPositionText, changeKEString, color='m')
    yPositionText = 0.75*yMax
    plt.text(xPositionText, yPositionText, dataInitializationFileString, color='m')

    plotNameString = 'FinalPositions-Particles%d-%s-BoxSize%d-T%dMilli-DeltaT%dMilli-FirstV%dSecondV2%d.pdf' % (nParticles, keNameText, boxSize, 1000*maxT, 1000*deltaT, 1000*v1, 1000*v2)
    figGraphs.savefig(plotNameString)
    plt.close(figGraphs)
    print '  Positions plot file name %s; mean radial move distance %.2f with RMS %.2f' % (plotNameString, meanMoveDistance, np.std(radialMove))

    #
    # Write the output data text file
    #
    fileOutput = open(textFileNameString, 'w')
    headerLineString = "  Final positions and velocities for a Molecular Dynamics simulation\n"
    nParticlesKENameBoxSizeMaxTDeltaTV1V2TemperaturesLatticeLineString = '%3d %s %.2f %.2f %.4f %.2f %.2f %s' % (nParticles, keNameText, boxSize, maxT, deltaT, v1, v2, '\n')
    fileOutput.write(headerLineString)
    fileOutput.write(nParticlesKENameBoxSizeMaxTDeltaTV1V2TemperaturesLatticeLineString)
    for kParticle in range(nParticles):
        xCurrent = currentPositionX[timeSteps][kParticle]
        yCurrent = currentPositionY[timeSteps][kParticle]
        xCurrentSpeed = currentSpeedX[timeSteps][kParticle]
        yCurrentSpeed = currentSpeedY[timeSteps][kParticle]
        dataLineString = '%6.4e %6.4e %7.5e %7.5e%s' % (xCurrent, yCurrent, xCurrentSpeed, yCurrentSpeed,'\n')
        fileOutput.write(dataLineString)

    fileOutput.close()
    print '  Data text file name %s has been written' % textFileNameString
    return

def init():                 # this initialization function tells what is being animiated
    """initialize animation"""
    global finalPlot, pathsPlot, timeValue_text, teleportValue_text
    
    finalPlot.set_data([], [])                 # (x,y) at a given time steps
    timeValue_text.set_text('')                # time step value
    teleportValue_text.set_text('')            # number of teleported particles

    if(noPaths):
        return finalPlot, timeValue_text, teleportValue_text
    else:
        pathsPlot.set_data([], [])                 # paths (x,y) for a set of time steps
        return finalPlot, pathsPlot, timeValue_text, teleportValue_text

def animate(i):   # this is the function this being animated, the i index is increased according the range of the frames value below
    """perform animation step"""           # i ranges from 0 to (nFrames - 1) as given in the FuncAnimation function
    global finalPlot, pathsPlot, timeValue_text, teleportValue_text

    simulationElapsedTime = float(i)*framesToTimesFactor
    timeIndex = int(simulationElapsedTime/deltaT)
    if(timeIndex < 0 or timeIndex > timeStepsPlusOne):
        print "\n *** [animate] Error in timeIndex value ", timeIndex, " with frame number ", i
        exit()
    timeStepString = '%s %.3e' % ('Time ', simulationElapsedTime)
    timeValue_text.set_text(timeStepString)
    if(specular):
        teleportValueString = '%s %d' % ('Boundary reflections ', nBoundary[timeIndex])
    else:
        teleportValueString = '%s %d' % ('Boundary crossings ', nBoundary[timeIndex])
    teleportValue_text.set_text(teleportValueString)
    xFinalPlot = [currentPositionX[timeIndex][kParticle] for kParticle in range(nParticles)]
    yFinalPlot = [currentPositionY[timeIndex][kParticle] for kParticle in range(nParticles)]
    finalPlot.set_data(xFinalPlot,yFinalPlot)

    if(noPaths):
        return finalPlot, timeValue_text, teleportValue_text

    xPathPlot = []
    yPathPlot = []
    timeIndexDivideTen = int(timeIndex/10)
    if(simulationElapsedTime/maxT > 0.333*maxT):
        startPathTimeIndex = int(0.667*maxT/deltaT)
    else:
        startPathTimeIndex = 1
    for kParticle in range(nParticles):
        for timeStep in range(timeIndexDivideTen):
            tt = timeStep*10
            if(tt >= timeIndex):
                print "\n [animate] ***Error in tt index", tt, " with timeIndex", timeIndex
                exit()
            if(tt >= startPathTimeIndex):
                xPosition = currentPositionX[tt][kParticle]
                yPosition = currentPositionY[tt][kParticle]
            
                xPathPlot.append(currentPositionX[tt][kParticle])
                yPathPlot.append(currentPositionY[tt][kParticle])

    pathsPlot.set_data(xPathPlot,yPathPlot)

    return finalPlot, pathsPlot, timeValue_text, teleportValue_text

def animationOutput():
    #
    # Function to produce animated output
    #
    global finalPlot, pathsPlot, timeValue_text, teleportValue_text
    
    # code to set up one plot in a single figure for the animation
    fig = plt.figure(1)          # start a figure
    ax = fig.add_subplot(111)
    finalPlot, = ax.plot([], [], 'go')   # this is an empty object used in the initialization of the animation, for the final particle positions
    pathsPlot, = ax.plot([], [], 'k.', markersize=1)   # this is an empty object used in the initialization of the animation, for the particle paths

    timeValue_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, color='blue')   # this is a placeholder for where the time step will be updated
    teleportValue_text = ax.text(0.45, 0.95, '', transform=ax.transAxes, color='blue')   # this is a placeholder for where number of teleported particles will be updated
 
    if(specular or teleport):
        plt.ylim(yMin, 1.1*yMax)
        plt.xlim(xMin, xMax)
    else:
        #
        # Expand plot range to see the escaping particles, draw the enclosing box
        #
        plt.ylim(yMin - 0.1*(yMax - yMin), yMax + 0.2*(yMax - yMin))
        plt.xlim(xMin - 0.1*(xMax - xMin), xMax + 0.1*(xMax - xMin))
        #
        # Draw the confining box
        #
        topLineX = [xMin, xMax]
        topLineY = [yMax, yMax]
        plt.plot(topLineX, topLineY, 'k')
        
        leftLineX = [xMin, xMin]
        leftLineY = [yMin, yMax]
        plt.plot(leftLineX, leftLineY, 'k')
        
        bottomLineX = [xMin, xMax]
        bottomLineY = [yMin, yMin]
        plt.plot(bottomLineX, bottomLineY, 'k')
        
        rightLineX = [xMax, xMax]
        rightLineY = [yMin, yMax]
        plt.plot(rightLineX, rightLineY, 'k')
    
    plt.xlabel('X (sigma units)')                             # add axis labels
    plt.ylabel('Y (sigma units)')
    xInitialPlot = [currentPositionX[0][kParticle] for kParticle in range(nParticles)]
    yInitialPlot = [currentPositionY[0][kParticle] for kParticle in range(nParticles)]
    plt.scatter(xInitialPlot, yInitialPlot,  s=2, c='r', marker='o')
    titleString = 'Molecular Dynamics for %2d particles, %.3f time step, box %.2f' % (nParticles, deltaT, boxSize)
    plt.title(titleString)
    plt.grid(True)
 
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
        if(specular):
            movieFileName = 'molecularDynamicsSpecular-Particles%d-DeltaT%dMilli-Box%d-MaxT%dMilli.mp4' % (nParticles, int(deltaT*1000), boxSize, int(maxT*1000))
        if(teleport):
            movieFileName = 'molecularDynamicsTeleport-Particles%d-DeltaT%dMilli-Box%d-MaxT%dMilli.mp4' % (nParticles, int(deltaT*1000), boxSize, int(maxT*1000))
        if(noTeleport):
            movieFileName = 'molecularDynamicsEscape-Particles%d-DeltaT%dMilli-Box%d-MaxT%dMilli.mp4' % (nParticles, int(deltaT*1000), boxSize, int(maxT*1000))
        ani.save(movieFileName, writer = FFwriter)
        print " Movie production step is completed"
    else:
        print "\n Movie production step is being skipped"

    plt.show()

    return

initMolecularDynanics()
calculateEnergy(0)

if(verbose):
    print "\n  Initial configuration"
    calculateSeparation(0)
    printParticleData(0)
    thisKineticEnergy = kineticEnergy[0]
    thisPotentialEnergy = potentialEnergy[0]
    print "  Kinetic Energy = ", thisKineticEnergy, ",  Potential Energy = ", thisPotentialEnergy, " , with Total Energy ", totalEnergy[0]
    netForceX, netForceY, averageForce = calculateAllForces(0)
    print "  Net force X = ", netForceX, ",  net force Y ", netForceY, " compared to an average force magnitude ", averageForce

#
# loop over time steps
#
timeStepsMinusOne = timeSteps - 1
startDateTimeString = str(datetime.now())
print "\n  Begin cycling over time steps at ", startDateTimeString
timeSkip = 1
if(timeSteps > 10):
    timeSkip = int(timeSteps/10)
for timeStep in range(timeSteps):

    summedMomentumX = 0.0
    summedMomentumY = 0.0
    if(fixCMVelocity or useDataInitializationFile):
        summedMomentumX, summedMomentumY, summedAngularMomentum = checkMomentum(timeStep)
        if(abs(summedMomentumX/nParticles) > 0.05 or abs(summedMomentumY/nParticles) > 0.05):
            print "\n ***Momemtum check error at time ", timeStep*deltaT, " xMomentum ", summedMomentumX, " y Momentum ", summedMomentumY
            exit()

    if(timeStep > 0 and changeKE):
        changeKineticEnergy(timeStep)

    updateMolecularDynamics(timeStep)
    calculateEnergy(timeStep)
    nBoundary[timeStep] = nBoundaryCount
    
    if((printEnergyStep) and (timeStep == 0 or timeStep%timeSkip == 0)):
        print "  Data at time ", deltaT*timeStep, ": KE ", kineticEnergy[timeStep], " PE ", potentialEnergy[timeStep], " Total ", totalEnergy[timeStep], "  Boundary hits ", nBoundaryCount
        if(fixCMVelocity or useDataInitializationFile):
            print "  x-Momentum ",summedMomentumX, " y-Momentum", summedMomentumY, ",  angular momentum ",  summedAngularMomentum, "\n"

    if(printForceStep and (timeStep == 0 or timeStep%timeSkip == 0)):
        netForceX, netForceY, averageForce = calculateAllForces(timeStep)
        print "  Net force X = ", netForceX, ",  net force Y ", netForceY, " compared to an average force magnitude ", averageForce

    timeStep += 1

endDateTimeString = str(datetime.now())
print "\n  End cycling over time steps at ", endDateTimeString

#
# Final Energy
#
calculateEnergy(timeSteps)
initialTotalEnergy = totalEnergy[0]
finalTotalEnergy = totalEnergy[timeSteps]
fractionalChange = (finalTotalEnergy - initialTotalEnergy)/abs(initialTotalEnergy)
meanTotalEnergy = np.mean(totalEnergy)
meanTotalEnergyStd = np.std(totalEnergy)

print "\n  Final energies at time %.2f : KE %.3f, PE %.3f, Total %.3f; fractional change %.3f" % (timeSteps*deltaT, kineticEnergy[timeSteps], potentialEnergy[timeSteps], finalTotalEnergy, fractionalChange)
print "  Average total energy %.3f with RMS %.3f and RMS/Average %.3e" % (meanTotalEnergy, meanTotalEnergyStd, meanTotalEnergyStd/abs(meanTotalEnergy))

energyString = 'Average total energy %.3f with RMS %.3f' % (meanTotalEnergy, meanTotalEnergyStd)

if(noAnimation):
    graphsOutput()
else:
    animationOutput()
