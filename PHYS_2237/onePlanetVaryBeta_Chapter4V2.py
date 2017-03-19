#
# Program 4.4 Solution for the orbit of a planet in the solar system with the Sun as fixed force center (onePlanetVaryBeta_Chapter4V2.py)
#
#              Give the command  python onePlanetVaryBeta_Chapter4V2.py -h  to get help command on input parameters
#
#              This program is the same as the onePlanetVaryBeta_Chapter4V1.py program except that it allows corrector changes in a circular orbit
#
#              The orbit corrector software has an application in high energy physics circular accelerators such as the Large Hadron Collider
#              The beam particles in the accelerators are subject to scattering from their intended orbits because of residual molecules in the high vacuum beam pipes
#              In order to nudge the particles back to the correct equilibrium point, a deviation amount is measured at one point in the circle
#              An electronic correction signal for this size deviation is calculated and this signal is sent across the diameter of the accelerator, several kilometers away
#              The correction signal travels at the speed of light and arrives at the opposite before the beam particles which have to travel in a longer semi-circle
#              The orbit correction signal is appled to the beam particles as they arrive, getting them back near the correct orbit
#              This clever scheme was developed by a Dutch accelerator physicist named Simon van der Meer who worked at CERN in the 1980s
#              For this work he shared the 1984 Nobel Physics prize with fellow CERN physicist Carlo Rubia, a native of Italy, who together discovered the W and Z bosons at CERN
#                   http://www.nobelprize.org/nobel_prizes/physics/laureates/1984/
#                   "for their decisive contributions to the large project, which led to the discovery of the field particles W and Z, communicators of weak interaction"
#
#              The derivative function has been moved closer to the odeint call; this has no effect on the operation of the program
#
#
#              The input arguments to this program are the same as for the onePlanet_Chapter4V2.py program, except for the addition of a beta exponent input
#
#                   1) dx/dt = vx    velocity component in the x direction
#                   2) dy/dt = vy    velocity component in the y direction
#                   3) dvx/dt = ax = -4*pi*pi*x/r^(1+beta)  acceleration component in the x direction from the beta exponent force
#                   4) dvy/dt = ay = -4*pi*pi*y/r^(1+beta)  acceleration component in the y direction from the beta exponent force
#
#              The parameters in the differential equation are according to
#                   1) time is expressed in units of years
#                   2) distance is expressed in Astronomical Units (AU) where 1 AU is the mean distance of the Earth's orbit
#                   3) v is expressed in AU/year
#                   4) the initial position of the orbit is taken to be the perihelion, i.e. the distance of closest approach, or an aphelion (furthest distance)
#                   4) the initial x0 position is then the perihelion (or aphelion) distances in units of AU, and vx0 = 0
#                   5) the initial y0 position of the orbit is 0, and vy0 is given in units of the Earth's perihelion speed which is 2*pi AU/year
#
#              Using input parameters for x0 (AU), vy0 (units of 2*pi AU/Y), mass of the planet (units of 10^24 kg), timeStep (years), maximumTime (years)
#              Defaults are beta = 2, x0 = 0.983 AU, vy0 = 1.017, mass 5.98 (Earth's mass), timeStep = 0.02 years, maximumTime = 3 years
#
#              The beta = 3 exponent is a special case and it does have analytic orbits which are called Cotes Spirals, named after the Roger Cotes (1682-1716)
#              He worked with Issac Newton but died young before Newton
#                  https://en.wikipedia.org/wiki/Cotes's_spiral
#                  http://mathworld.wolfram.com/CotesSpiral.html
#              With beta = 3 the orbits will be decaying or expanding spirals, depending on whether the energy is negative or positive
#              If the energy is 0 for beta = 3, then the orbit will be a circle but that orbit is unstable and numerical instablities will set in after multiple orbit calculations
#

import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt     # get matplotlib plot functions
import sys                          # used to get the number of command line arguments
import  argparse                    # argument parser library
import numpy as np                  # numerical functions library used by python
import math as mp                   # used for the exponential function
from scipy.integrate import odeint  # import only this single method for solving differential equations
import random                       # random number generator library

#
# Astronomy constants in MKS units, used to compute the initial and final total energies in Joules of the orbiting planet
#

ONEAU = 1.496e+11                           # Earth's distance from the Sun in meters
AUCUBE = ONEAU*ONEAU*ONEAU                  # cube of AU distance in m^3
YEAR = 3.154e7                              # year in seconds
YEARSQUARE = YEAR*YEAR                      # square of year in s^2
GSUNMASS = 4*mp.pi*mp.pi*AUCUBE/YEARSQUARE  # 4pi*pi*G*SunMass

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()
parser.add_argument('--x0', default=1.0, type=float, help="Initial perihelion x0 in AU; default 0.983 AU")
parser.add_argument('--vy0', default=5.0, type=float, help="Initial perihelion vy0 in units of 2Pi AU/Y; default 1.017")
parser.add_argument('--maxT', default=1.0, type=float, help="Maximum time range in years; default 1 Year")
parser.add_argument('--deltaT', default=.0001, type=float, help="Iteration time step in years; default 0.02 Year")
parser.add_argument('--mass', default=5.98, type=float, help="Planet mass in units of 10^24 kg; default 5.98")
parser.add_argument('--beta', default=2.5, type=float, help="Exponent for the radius variable in the Universal Gravity law; default 2.0")
parser.add_argument('--safetyRadius', default=1.0e-05, type=float, help="Smallest allowable radial distance; default 1.0e-05")
parser.add_argument('--pointsPlot', action='store_true', help="Plot discrete data points instead of a continuous line; default False")
parser.add_argument('--verbose', action='store_true', help="Give extra printout about orbit pattern recognition results; default False")
parser.add_argument('--maxPointsToPlot', default=10000, type=int, help="Maximum number of orbit points to be plotted; default 10000")
parser.add_argument('--useCorrector',action='store_true', help="Correct the orbit dynamically; default False")
parser.add_argument('--correctBand', default=0.005 , type=float, help="Band in which to apply corrector algorithm; default 0.005")
args = parser.parse_args()

#
# Get the input parameters from the command line
#
numberOfArguments = len(sys.argv)
if(numberOfArguments == 1):
    print "\n All the default paramter choices are used"  # there is always at least one argument and that first one is the name of the python script itself

#
# Assign the variables in the program to the variable option names in the argument list
#
pointsPlot = args.pointsPlot
verbose = args.verbose
maxPointsToPlot = args.maxPointsToPlot
useCorrector = args.useCorrector
correctBand = args.correctBand

beta=args.beta                             # is beta = 2, this program should give the same results as for Newton's Universal Gravity law
betaPlusOne = beta + 1.0                   # used in the denominator for the acceleration derivative components, instead of 3
betaMinusOne = beta - 1.0                  # used in corrected potential energy function

auPowerBetaMinusTwo = pow(ONEAU,beta-2.0)  # correction factor for GSUNMASS when beta is not 2.0

safetyRadius = args.safetyRadius

x0 = args.x0
if(x0 <= 0.0):
    print "\n Cannot have the initial position", x0, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(0)

vy0 = args.vy0
if(vy0 == 0.0):
    print "\n Cannot have the initial velocity omponent", vy0, " be equal to zero"
    print "\n Program is exiting\n"
    exit(0)

if(vy0 > 0):
    vy0Speed = vy0*2*mp.pi     # in units of AU/Year
else:
    vy0Speed = -vy0            # in AU/Year, to compare with MATLAB results

r0v0Product = x0*vy0Speed
if(useCorrector and beta == 3 and r0v0Product != 2*mp.pi):
    print "\n *** Input error: cannot request orbit correction with for beta = 3 force if the initial orbit is not circular ***"
    exit(1)

randomSeedIsFixed = True
if(useCorrector):
    if(correctBand < 0):
        correctBand = -correctBand
        random.seed()   # random seed is set according to the system clock, numerical results will be non-repeatable
        randomSeedIsFixed = False
    else:
        random.seed(1)  # fixes the sequence of random sees for repeatability of numerical results

maximumTime = args.maxT
timeStep = args.deltaT
mass = args.mass        # Note that the mass of the planet does not affect its orbit size or orbital period, as per Kepler's Third Law
masskg = mass*1.e24

fourPiSquare = 4*mp.pi*mp.pi                       # simple universal gravity force scaling factor when astronomical units are used (note the planet mass value is not used)

print "\n Orbit of a planet around the Sun as a fixed force center, with beta = ", beta
print "         Input conditions"
print "  Initial radius ", x0, " AU, at an initial transverse velocity of ", vy0Speed, " AU/Year"
print "  Time step = ", timeStep, " year,  maximum time range = ", maximumTime, " years, with a minimum safety radius set at ", safetyRadius, " AU"

initialEnergy = -GSUNMASS*masskg*auPowerBetaMinusTwo/(betaMinusOne*pow(x0*ONEAU,betaMinusOne)) + 0.5*masskg*pow(vy0Speed*ONEAU/YEAR, 2)
print "\n  Initial potential energy ", -GSUNMASS*masskg*auPowerBetaMinusTwo/(betaMinusOne*pow(x0*ONEAU,betaMinusOne)), " J, initial kinetic energy ", 0.5*masskg*pow(vy0Speed*ONEAU/YEAR, 2), " J"
print "  Initial total energy of the planet = ", initialEnergy, " Joules",
if(beta < 3):
    if(initialEnergy < 0 ):
        print ", there is a bound orbit"
    if(initialEnergy > 0):
        print ", there is no bound orbit"
    if(initialEnergy == 0):
        print ", the orbit is parabolic"
if(beta == 3):
    if(initialEnergy < 0 ):
        print ", there is a decaying spiral bound orbit"
    if(initialEnergy > 0):
        print ", there is an expanding spiral unbound orbit"
    if(initialEnergy == 0):
        print ", the orbit is an unstable circle"
if(beta > 3):
    print ", there has been no study of these orbit shapes in class"

initialAngularMomentum = masskg*x0*ONEAU*vy0Speed*ONEAU/YEAR
print "\n  Initial angular momentum = ", initialAngularMomentum, " kg-m/s"

if(useCorrector):
    print "\n An orbit corrector algorithm will be used for an error band ", correctBand,
    if(randomSeedIsFixed):
        print ", the sequence of random numbers is fixed"
    else:
        print ", the sequence of random numbers is not fixed"

# dx/dt = vx    velocity component in the x direction
# dy/dt = vy    velocity component in the y direction
# dvx/dt = ax = -4*pi*pi*x/r^(1+beta)  acceleration component in the x direction from the universal gravity force  (instead of cubed, the denominator power becomes 1 + beta)
# dvy/dt = ay = -4*pi*pi*y/r^(1+beta)  acceleration component in the y direction from the universal gravity force  (instead of cubed, the denominator power becomes 1 + beta)

nCorrectorTimes = 0
debug = False

def fDerivative(variableList, t):                  # variableList dummy list array since there is more than one differential equation
    global nCorrectorTimes, debug, lastDisruptTime, nDisruptorTimes
    vx = variableList[0]                           # speed in the x direction
    vy = variableList[1]                           # speed in the y direction
    x = variableList[2]                            # x coordinate
    y = variableList[3]                            # y coordinate
    rDistance = mp.sqrt(x*x + y*y)                 # radial distance from the force center

    if(useCorrector):
        if(beta == 3):
            #
            # Orbit correction alogorithm for beta = 3 case assumses an initial circular orbit, and corrects the orbit if the deviation from that orbit is larger than the error band
            # The size of the orbit correction is a randomly chosen to move the orbit distance closer to the original circular orbit value
            # The velocity is also changed in order to maintain conservation of angular momentum
            #
            deltaRadius = rDistance - x0
            radiusErrorFraction = abs(deltaRadius)/x0
            if(radiusErrorFraction > correctBand):
                nCorrectorTimes += 1
                randomNumberFraction = random.uniform(-0.5, 0.5)   # the radius error will be reduced by not eliminated
                newRadialPosition = x0 + deltaRadius*randomNumberFraction   # this radial position be closer to the original circle radius
                currentTheta = np.arctan2(y,x)
                newY = newRadialPosition*np.sin(currentTheta)  # maintain the same angle position after the radius correction
                newX = newRadialPosition*np.cos(currentTheta)  # maintain the same angle position after the radius correction
                #
                #  reset the velocity to give 0 energy at this new radial position
                #
                #  specificPotentialEnergy = -4.mp.pi*mp.pi/newRadialPosition; this is  potential energy in Astronomy units, with the planet mass divided out
                #
                newVelocity = 2*mp.pi*mp.sqrt(1.0/newRadialPosition)
                currentVelocity = mp.sqrt(vx*vx + vy*vy)
                currentThetaV = np.arctan2(vy,vx)
                newVy = newVelocity*np.sin(currentThetaV)
                newVx = newVelocity*np.cos(currentThetaV)
                oldAngularMomentum = x*vy - y*vx
                newAngularMomentum = newX*newVy - newY*newVx

                if(debug):
                    print "\n  deltaRadius ", deltaRadius, ",  radiusErrorFraction ", radiusErrorFraction, ", randomNumberFraction ", randomNumberFraction
                    print "  rDistance ", rDistance, ",  newRadialPosition ", newRadialPosition, ",  currentTheta ", currentTheta
                    print "  newX ", newX, ",  newY ", newY, ",  x ", x, ",  y ", y
                    print "  newVx ",  newVx, ",  newVy ", newVy, ",  vx ", vx, ",  vy ", vy
                    print "  currentThetaV ", currentThetaV, ",  newVelocity ", newVelocity, ",  currentVelocity ", currentVelocity, ",  originalVelocity", vy0Speed
                    print "  r0v0Product ", r0v0Product, ",  newAngularMomentum ", newAngularMomentum, ",  oldAngularMomentum ", oldAngularMomentum
                    print " "
                    exit()  # in debug mode to check the changes done by the algorithm

                x = newX
                y = newY
                rDistance = newRadialPosition
                vx = newVx
                vy = newVy

    if(rDistance < safetyRadius):
        print "\n Radial distance is too small at ", rDistance, " AU happening at time ", t, " years"
        exit()
    rDistanceBeta = pow(rDistance, betaPlusOne)    # radial distance raised to the power beta, replaces the previous rDistanceCubed
    dvxdt = -fourPiSquare*x/rDistanceBeta          # the time derivative of velocity in the x direction according to the universal gravity force component
    dvydt = -fourPiSquare*y/rDistanceBeta          # the time derivative of postion in the y direction is the universal gravity force component
    return [dvxdt, dvydt, vx, vy]                  # return the four derivatives as a list object containing four elements

#
# Notice that these lines of code to solve the differential equations are nearly identical to those of the projectile motion code
#
nTimeSteps = int(maximumTime/timeStep)
# obtain the differential equation solutions using the odeint method from the ScyPy library
timeGrid = np.linspace(0, maximumTime,nTimeSteps)           # time grid used by odeint method
initialValuesSpeedsPositions = [0, vy0Speed, x0, 0.0]       # starting values of vx, vy, x0, y0 for the iteration
fourSolution = odeint(fDerivative, initialValuesSpeedsPositions, timeGrid)   # odeint returns a list of values which are the numerical solutions
vxNumerical = fourSolution[:,0]            # vx function of time obtained with Numerical solution
vyNumerical = fourSolution[:,1]            # vy function of time obtained with Numerical solution
xNumerical = fourSolution[:,2]             # x function of time obtained with Numerical solution
yNumerical = fourSolution[:,3]             # y function of time obtained with Numerical solution

print "\n  Number of times that the corrector algorithm was used " , nCorrectorTimes

# The above lines of code have solved the differential equations for a planet in orbit around the Sun
# The next lines of code print and plot the results


nTimeStepsMinusOne = nTimeSteps - 1
xFinal = xNumerical[nTimeStepsMinusOne]
yFinal = yNumerical[nTimeStepsMinusOne]
rFinal = np.sqrt(xFinal*xFinal + yFinal*yFinal)
vxFinal = vxNumerical[nTimeStepsMinusOne]
vyFinal = vyNumerical[nTimeStepsMinusOne]
vFinal = np.sqrt(vxFinal*vxFinal + vyFinal*vyFinal)
finalEnergy = -GSUNMASS*masskg*auPowerBetaMinusTwo/(betaMinusOne*pow(rFinal*ONEAU,betaMinusOne)) + 0.5*masskg*pow(vFinal*ONEAU/YEAR, 2)

print "\n\n          Results from numerical solutions for position and velocity"
print "  Final potential energy ", -GSUNMASS*masskg*auPowerBetaMinusTwo/(betaMinusOne*pow(rFinal*ONEAU,betaMinusOne)), " J, final kinetic energy ", 0.5*masskg*pow(vFinal*ONEAU/YEAR, 2), " J"
print "  Final total energy of the planet = ", finalEnergy, " Joules"
if(initialEnergy != 0.0):
    fractionalEnergyChange = (finalEnergy - initialEnergy)/initialEnergy
    print "  The fractional energy change is ", fractionalEnergyChange
else:
    print " The initial energy was zero and the final energy is ", finalEnergy
#
# In Cartesian coordinates the final angular momentum is computed in terms of the position and velocity components
#  angular momemtum (in +Z direction) = mass*(rx*vy - ry*vx)
#
finalAngularMomentum = masskg*(xFinal*ONEAU*vyFinal*ONEAU/YEAR - yFinal*ONEAU*vxFinal*ONEAU/YEAR)
print "\n  Final angular momentum = ", finalAngularMomentum, " kg-m/s"
fractionalAngularMomentumChange = (finalAngularMomentum - initialAngularMomentum)/initialAngularMomentum
print "  The fractional angular momentum change is ", fractionalAngularMomentumChange

radialPosition = []
totalVelocity = []
totalEnergy = []

nTimeStep = 0
while nTimeStep < nTimeSteps:
    radialValue = mp.sqrt(xNumerical[nTimeStep]*xNumerical[nTimeStep] + yNumerical[nTimeStep]*yNumerical[nTimeStep])
    radialPosition.append(radialValue)
    velocityValue = mp.sqrt(vxNumerical[nTimeStep]*vxNumerical[nTimeStep] + vyNumerical[nTimeStep]*vyNumerical[nTimeStep])
    totalVelocity.append(velocityValue)
    energyValue = -GSUNMASS*masskg*auPowerBetaMinusTwo/(betaMinusOne*pow(radialValue*ONEAU,betaMinusOne)) + 0.5*masskg*pow(velocityValue*ONEAU/YEAR, 2)
    totalEnergy.append(energyValue)

    nTimeStep += 1

meanRadialPosition = np.mean(radialPosition)
stdRadialPosition = np.std(radialPosition)
meanTotalVelocity = np.mean(totalVelocity)
stdTotalVelocity = np.std(totalVelocity)
meanTotalEnergy = np.mean(totalEnergy)
stdTotalEnergy = np.std(totalEnergy)

print "\n  Mean radial position = ", meanRadialPosition, " +/- ", stdRadialPosition, " AU"
print "  Mean total velocity = ", meanTotalVelocity, " +/- ", stdTotalVelocity, " AU/Year"
print "  Mean total energy = ", meanTotalEnergy, " +/- ", stdTotalEnergy, " Joules"

ellipticalOrbit = False
circularOrbit = False
if(beta != 2):
    orbitShapeString = "Orbit pattern recognition is not done for this beta"

numberOfOrbits = 0
if(beta == 2):
    #
    # Do the orbit pattern recognition to extract the orbit time and the orbit eccentricity
    # First use the radial position and the total velocity arrays
    #

    circularOrbit = False
    ellipticalOrbit = True
    if(stdRadialPosition/meanRadialPosition < 1.e-05 and stdTotalVelocity/meanTotalVelocity):
        orbitShapeString = 'The orbit is determined to be circular'
        ellipticalOrbit = False
        circularOrbit = True
    else:
        orbitShapeString = 'The orbit is determined to be non-circular'

    print "    ", orbitShapeString

    eccentricity = 0
    orbitTime = 0
    numberOfOrbits = 0
    if(circularOrbit):
        orbitTime = 2*mp.pi*meanRadialPosition/meanTotalVelocity
        numberOfOrbits = int(maximumTime/orbitTime)

if(ellipticalOrbit and beta == 2):
    print "\n       Doing advanced pattern recognition analysis of the orbits"
    #
    # Algorithm is to accumulate the set of perihelion and aphelion radial values and their time values
    # The times to pass through successive perihelion points are stored from which a mean and standard deviation are computed
    # The first step is to confirm that the initial point is consistent with being a perihelion
    # The program input model is that the orbit starts as either a perihelion or an aphelion
    #
    radiusPerihelion = []
    radiusAphelion = []
    if(radialPosition[1] > radialPosition[0] and totalVelocity[1] < totalVelocity[0]):
        radiusPerihelion.append(radialPosition[0])
        lastPerihelionTime = 0.0
        lookingForNextAphelion = True
        lookingForNextPerihelion = False
        startAsPerihelion = True
        foundPerihelion = True
        foundAphelion = False
        print "  Check that the initial position is consistent with being a perihelion is passed"
    else:
        radiusAphelion.append(radialPosition[0])
        lastAphelionTime = 0
        lookingForNextAphelion = False
        lookingForNextPerihelion = True
        startAsPerihelion = False
        foundPerihelion = False
        foundAphelion = True
        print "  Check that the initial position is consistent with being an aphelion is passed"

    lastRadialPosition = radialPosition[0]
    timePerihelion = []
    timeAphelion = []
    nTimeStep = 1         # start the orbit pattern recognition at second point on orbit
    while nTimeStep < nTimeSteps:
        #
        # As the while loop iterates through the positions at each time step, it will be either looking for the next aphelion or the next perihelion
        #
        newRadialPosition = radialPosition[nTimeStep]

        if(lookingForNextAphelion):
            #
            # If the new radial position is smaller than the previous radial position, then the next aphelion point has been crossed
            #
            if(newRadialPosition < lastRadialPosition):
                timeValue = 0.5*(timeGrid[nTimeStep] + timeGrid[nTimeStep-1])  # take an average of the current and previous times
                if(foundAphelion == False):
                    foundAphelion = True             # found the first aphelion point
                    lastAphelionTime = timeValue     # store the time of the first aphelion
                    if(verbose):
                        print "\n Found first aphelion at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                else:
                    if(verbose):
                        print "\n Found next aphelion at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                    timeAphelion.append(timeValue - lastAphelionTime)     # store the time difference since the last aphelion
                    lastAphelionTime = timeValue
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
                timeValue = 0.5*(timeGrid[nTimeStep] + timeGrid[nTimeStep-1])  # take an average of the current and previous times
                if(foundPerihelion == False):
                    foundPerihelion = True             # found the first perihelion point
                    lastPerihelionTime = timeValue     # store the time of the first perihelion
                    if(verbose):
                        print "\n Found first perihelion at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                else:
                    if(verbose):
                        print "\n Found next perihelion at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                    timePerihelion.append(timeValue - lastPerihelionTime)     # store the time difference since the last perihelion
                    lastPerihelionTime = timeValue
                    radiusPerihelion.append(0.5*(lastRadialPosition+newRadialPosition))
                lookingForNextPerihelion = False
                lookingForNextAphelion = True

            lastRadialPosition = newRadialPosition

        nTimeStep += 1

    nPerihelion = len(radiusPerihelion)
    nAphelion = len(radiusAphelion)
    if(verbose):
        print "Radius Aphelion with number of points ", nAphelion
        print radiusAphelion
        print "Radius Perihelion, with number of points ", nPerihelion
        print radiusPerihelion

    if(nPerihelion > 0 and nAphelion > 0):        # require at least one found perihelion and one found aphelion in order to quote orbit information
        majorAxisLength = np.mean(radiusAphelion) + np.mean(radiusPerihelion)
        eccentricity = (np.mean(radiusAphelion) - np.mean(radiusPerihelion))/majorAxisLength
        print '%s %5.3f' % ('  The orbit eccentricity = ', eccentricity)

        nPerihelionTimeDiff = len(timePerihelion)
        nAphelionTimeDiff = len(timeAphelion)
        meanTimeBetweenAphelion = np.mean(timeAphelion)
        meanTimeBetweenPerihelion = np.mean(timePerihelion)
        if(verbose):
            print "\n Number of aphelion points = ", nAphelion, ", number of perihelion points ", nPerihelion
            print " Mean time between aphelion points = ", meanTimeBetweenAphelion
            print " Mean time between perihelion points = ", meanTimeBetweenPerihelion

        orbitTime = (nPerihelionTimeDiff*meanTimeBetweenPerihelion + nAphelionTimeDiff*meanTimeBetweenAphelion)/float(nPerihelionTimeDiff + nAphelionTimeDiff)
        if(startAsPerihelion):
            numberOfOrbits = nPerihelion
        else:
            numberOfOrbits = nAphelion

if(beta == 2):
    print '%s %d' % ('  Number of completed orbits = ', numberOfOrbits)
    if(numberOfOrbits > 0):
        keplerThirdLawTime = pow(majorAxisLength/2.0, 1.5)
        print '%s %5.3f %s' % ('  Orbit time = ', orbitTime, ' years')
        print '%s %5.3f %s' % ('  Orbit time predicted from Kepler Third law = ', keplerThirdLawTime, ' year, using the semi-major axis length determined from the pattern recognition')

#
# Expand the plot limits by 10% around the minimum and maximum orbital positions
#
minimumPositionX = 1.1*min(xNumerical)
minimumPositionY = 1.1*min(yNumerical)
maximumPositionX = 1.1*max(xNumerical)
maximumPositionY = 1.1*max(yNumerical)

# compose string variables about the time parameters for use in putting text on the plot
x0String =  'Initial radial position = ' + str(x0) + ' AU'
if(vy0 > 0):
    vy0String = 'Initial transverse velocity = ' + str(vy0) + '*2Pi AU/Y'
else:
    vy0String = 'Initial transverse velocity = ' + str(vy0Speed) + ' AU/Y'

if(useCorrector):
    betaString = 'Beta = ' + str(beta) + ' with orbit corrector algorithm'
else:
    betaString = 'Beta = ' + str(beta) + ' with no orbit corrector algorithm'

massString = 'Mass = ' + str(mass*1.0e24) + ' kg'
timeString = 'Total time = ' + str(maximumTime) + ' years in steps of ' + str(timeStep) + ' years'


plt.figure(1)       # start a figure for a single plot of the orbit
xSun = []
ySun = []
xSun.append(0)
ySun.append(0)
plt.scatter(0, 0, c='y', s=200)

#
# Check on the number of data points compared to the maximum number of points to be plotted
# Having a huge number of data points to be plotted can take extra time, or even cause software library errors
#
if(nTimeSteps < maxPointsToPlot):
    if(pointsPlot):
        plt.plot(xNumerical, yNumerical, 'bo')    # plot as discrete data points
    else:
        plt.plot(xNumerical, yNumerical)    # plot as a continuous line
else:
    skipTimeInterval = nTimeSteps/maxPointsToPlot
    if(skipTimeInterval < 1):
        print "\n Program Error, skipTimeInterval ", skipTimeInterval
        exit()
    else:
        if(verbose):
            print "\n Plotting skip time value ", skipTimeInterval*timeStep, " years"
    nTimeStep = 0
    xplot = []
    yplot = []
    while nTimeStep < nTimeSteps:
        xplot.append(xNumerical[nTimeStep])
        yplot.append(yNumerical[nTimeStep])
        nTimeStep += skipTimeInterval

    if(pointsPlot):
        plt.plot(xplot, yplot, 'bo')    # plot as discrete data points
    else:
        plt.plot(xplot, yplot)    # plot as a continuous line

xTextPosition = minimumPositionX + 0.25*(maximumPositionX - minimumPositionX)
yTextPosition = minimumPositionY + 0.55*(maximumPositionY - minimumPositionY)

plt.text(xTextPosition, minimumPositionY + 0.70*(maximumPositionY - minimumPositionY), betaString)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.63*(maximumPositionY - minimumPositionY), x0String)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.56*(maximumPositionY - minimumPositionY), vy0String)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.44*(maximumPositionY - minimumPositionY), massString)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.37*(maximumPositionY - minimumPositionY), timeString)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.30*(maximumPositionY - minimumPositionY), orbitShapeString)   # text to document the parameters used
if(numberOfOrbits > 0 and beta == 2):
    orbitTimeString = 'Orbital time = ' + str(orbitTime) + ' years'
    eccentricityString = 'Orbit eccentricity = ' + str(eccentricity)
    plt.text(xTextPosition, minimumPositionY + 0.23*(maximumPositionY - minimumPositionY), orbitTimeString)   # text to document the parameters used
    plt.text(xTextPosition, minimumPositionY + 0.16*(maximumPositionY - minimumPositionY), eccentricityString)   # text to document the parameters used
if(beta == 3):
    finalRadiusString = 'Final radius = ' + str(rFinal) + ' AU'
    vFinalMKS = vFinal*ONEAU/YEAR
    finalVelocityString = 'Final velocity = ' + str(vFinalMKS) + ' m/s'
    plt.text(xTextPosition, minimumPositionY + 0.23*(maximumPositionY - minimumPositionY), finalRadiusString)   # text to document the parameters used
    plt.text(xTextPosition, minimumPositionY + 0.16*(maximumPositionY - minimumPositionY), finalVelocityString)   # text to document the parameters used
    yLastValue = yNumerical[nTimeSteps - 1]
    xLastValue = xNumerical[nTimeSteps - 1]
    finalTheta = np.arctan2(yLastValue,xLastValue)
    if(finalTheta < 0.0):
        finalTheta += 2*mp.pi   # choose a 0 -> 2pi range for the final theta print out
    print "\n  The final orbit angle is ", finalTheta, " radians at a final radius position ", rFinal, " AU"
    x0v0 = x0*vy0Speed
    if(x0v0 > 2*mp.pi):
        predictedAymptoticTheta = (x0v0/np.sqrt(x0v0*x0v0 - 4.0*mp.pi*mp.pi))*mp.pi/2.
        print "  For an unbound beta = 3 orbit, the asymptotic theta value for these initial conditions is " , predictedAymptoticTheta, " radians"
plt.xlabel('x Coordinate (AU)')                             # add axis labels
plt.ylabel('y Coordinate (AU)')

plt.title('Orbit of a Planet in the Solar System, Variable Beta')

plt.grid(True)
plt.xlim(minimumPositionX, maximumPositionX)
plt.ylim(minimumPositionY, maximumPositionY)

plt.show()          # show the figure
