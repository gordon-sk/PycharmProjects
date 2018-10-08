#
# Program 4.6 Solution for the orbit of a planet in the solar system with the Sun as fixed force center (onePlanetPrecessAlpha_Chapter4V1.py)
#
#              Give the command  python onePlanetPrecess_Chapter4V1.py -h  to get help command on input parameters
#
#              This program is the same as the onePlanet_Chapter4V2.py program except that it allows variation from the inverse square gravity law from General Relativity (GR)
#              The Newton Universal Gravity law is multiplied by a dimensionless factor (1 + alpha/r^2) where alpha is predicted by GR to be 1.1 x 10^(-8) AU*AU
#              Code needs to modified to have the correct potential energy function for this force, but consevation of energy is not needed for this precession study
#
#              The input arguments to this program are the same as for the onePlanet_Chapter4V2.py program, except for the addition of a beta exponent input
#
#                   1) dx/dt = vx    velocity component in the x direction
#                   2) dy/dt = vy    velocity component in the y direction
#                   3) dvx/dt = ax = -4*pi*pi*x/r^3  acceleration component in the x direction from the universal gravity force
#                   4) dvy/dt = ay = -4*pi*pi*y/r^3  acceleration component in the y direction from the universal gravity force
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
#              Defaults are alpha = 0.01, x0 = 0.310 AU, vy0 = 1.973, mass 0.330 (Mercury's mass), timeStep = 0.0001 years, maximumTime = 10 years
#
#
import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt     # get matplotlib plot functions
import sys                          # used to get the number of command line arguments
import  argparse                    # argument parser library
import numpy as np                  # numerical functions library used by python
import math as mp                   # used for the exponential function
from scipy.integrate import odeint  # import only this single method for solving differential equations
import scipy.optimize as optimization

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
parser.add_argument('--x0', default=0.310, type=float, help="Initial perihelion x0 in AU; default 0.310 AU")
parser.add_argument('--vy0', default=1.973, type=float, help="Initial perihelion vy0 in units of 2Pi AU/Y; default 1.973")
parser.add_argument('--maxT', default=5.0, type=float, help="Maximum time range in years; default 5 Years")
parser.add_argument('--minPlotT', default=0.0, type=float, help="Minimum time at which to start plotting the orbit; default 0 Years")
parser.add_argument('--maxPlotT', default=5.0, type=float, help="Maximum time at which to end plotting the orbit; default 5 Years")
parser.add_argument('--deltaT', default=0.001, type=float, help="Iteration time step in years; default 0.001 Year")
parser.add_argument('--mass', default=0.330, type=float, help="Planet mass in units of 10^24 kg; default 0.330")
parser.add_argument('--alpha', default=0.0008, type=float, help="General Relativity correction to Newton's Universal Gravity law; default 0.0008")
parser.add_argument('--safetyRadius', default=1.0e-05, type=float, help="Smallest allowable radial distance; default 1.0e-05")
parser.add_argument('--pointsPlot', action='store_true', help="Plot discrete data points instead of a continuous line; default False")
parser.add_argument('--verbose', action='store_true', help="Give extra printout about orbit pattern recognition results; default False")
parser.add_argument('--maxPointsToPlot', default=10000, type=int, help="Maximum number of orbit points to be plotted; default 10000")

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

alpha=args.alpha                           # if alpha = 0, this program should give the same results as for Newton's Universal Gravity law
precessionAngleError = 0.05                # in degrees; parameter needed in fitting the rate of precession but it has no effect since all the angle detemination errors are assumed equal
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

maximumTime = args.maxT
if(maximumTime > 5.0):
    print "\n  WARNING: Results for a maximum time > 5 year will have a problem with perihelion wrapping around more than 360 degrees"

minPlotT = args.minPlotT
maxPlotT = args.maxPlotT
timeStep = args.deltaT
if(timeStep > 0.003 or timeStep < 0.0007):
    print "\n  WARNING: Results are best when the time step is between 0.0007 and 0.003 years"

mass = args.mass        # Note that the mass of the planet does not affect its orbit size or orbital period, as per Kepler's Third Law
masskg = mass*1.e24

# dx/dt = vx    velocity component in the x direction
# dy/dt = vy    velocity component in the y direction
# dvx/dt = ax = -4*pi*pi*x/r^(1+beta)  acceleration component in the x direction from the universal gravity force  (instead of cubed, the denominator power becomes 1 + beta)
# dvy/dt = ay = -4*pi*pi*y/r^(1+beta)  acceleration component in the y direction from the universal gravity force  (instead of cubed, the denominator power becomes 1 + beta)

fourPiSquare = 4*mp.pi*mp.pi                       # simple universal gravity force scaling factor when astronomical units are used (note the planet mass value is not used)
def fDerivative(variableList, t):                  # variableList dummy list array since there is more than one differential equation
    vx = variableList[0]                           # speed in the x direction
    vy = variableList[1]                           # speed in the y direction
    x = variableList[2]                            # x coordinate
    y = variableList[3]                            # y coordinate
    rDistanceSquare = x*x + y*y                    # square of radial distance from the force center
    rDistance = mp.sqrt(rDistanceSquare)           # radial distance from the force center
    if(rDistance < safetyRadius):
        print "\n Radial distance is too small at ", rDistance, " AU happening at time ", t, " years"
        exit()
    rDistanceCube = pow(rDistance, 3)                      # radial distance cubed
    alphaFactor = 1.0 + alpha/rDistanceSquare        # correction for General Relativity
    dvxdt = -fourPiSquare*alphaFactor*x/rDistanceCube      # the time derivative of velocity in the x direction according to the universal gravity force component
    dvydt = -fourPiSquare*alphaFactor*y/rDistanceCube      # the time derivative of postion in the y direction is the universal gravity force component
    return [dvxdt, dvydt, vx, vy]                  # return the four derivatives as a list object containing four elements

print "\n Orbit of a planet around the Sun as a fixed force center with alpha = ", alpha
print "         Input conditions"
print "  Initial radius ", x0, " AU, at an initial transverse velocity of ", vy0Speed, " AU/Year"
print "  Time step = ", timeStep, " year,  maximum time range = ", maximumTime, " years, with a minimum safety radius set at ", safetyRadius, " AU"

initialEnergy = -GSUNMASS*masskg/(x0*ONEAU) + 0.5*masskg*pow(vy0Speed*ONEAU/YEAR, 2)
print "\n  Initial potential energy ", -GSUNMASS*masskg/(x0*ONEAU), " J, initial kinetic energy ", 0.5*masskg*pow(vy0Speed*ONEAU/YEAR, 2), " J"
print "  Initial total energy of the planet = ", initialEnergy, " Joules",
if(initialEnergy < 0 ):
    print ", there is a bound orbit"
if(initialEnergy > 0):
    print ", there is no bound orbit"
if(initialEnergy == 0):
    print ", the orbit is parabolic"

initialAngularMomentum = masskg*x0*ONEAU*vy0Speed*ONEAU/YEAR
print "\n  Initial angular momentum = ", initialAngularMomentum, " kg-m/s"

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

# The above lines of code have solved the differential equations for a planet in orbit around the Sun
# The next lines of code print and plot the results


nTimeStepsMinusOne = nTimeSteps - 1
xFinal = xNumerical[nTimeStepsMinusOne]
yFinal = yNumerical[nTimeStepsMinusOne]
rFinal = np.sqrt(xFinal*xFinal + yFinal*yFinal)
vxFinal = vxNumerical[nTimeStepsMinusOne]
vyFinal = vyNumerical[nTimeStepsMinusOne]
vFinal = np.sqrt(vxFinal*vxFinal + vyFinal*vyFinal)
finalEnergy = -GSUNMASS*masskg/(rFinal*ONEAU) + 0.5*masskg*pow(vFinal*ONEAU/YEAR, 2)

print "\n\n          Results from numerical solutions for position and velocity"
print "  Final potential energy ", -GSUNMASS*masskg/(rFinal*ONEAU), " J, final inetic energy ", 0.5*masskg*pow(vFinal*ONEAU/YEAR, 2), " J"
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

#
# Do the orbit pattern recognition to extract the orbit time and the orbit eccentricity
# First set up the radial position and the total velocity arrays
#
radialPosition = []
totalVelocity = []
totalEnergy = []

nTimeStep = 0
while nTimeStep < nTimeSteps:
    radialValue = mp.sqrt(xNumerical[nTimeStep]*xNumerical[nTimeStep] + yNumerical[nTimeStep]*yNumerical[nTimeStep])
    radialPosition.append(radialValue)
    velocityValue = mp.sqrt(vxNumerical[nTimeStep]*vxNumerical[nTimeStep] + vyNumerical[nTimeStep]*vyNumerical[nTimeStep])
    totalVelocity.append(velocityValue)
    energyValue = -GSUNMASS*masskg/(radialValue*ONEAU) + 0.5*masskg*pow(velocityValue*ONEAU/YEAR, 2)
    totalEnergy.append(energyValue)

    nTimeStep += 1

meanRadialPosition = np.mean(radialPosition)
stdRadialPosition = np.std(radialPosition)
meanTotalVelocity = np.mean(totalVelocity)
stdTotalVelocity = np.std(totalVelocity)
meanTotalEnergy = np.mean(totalEnergy)
stdTotalEnergy = np.std(totalEnergy)

print "\n Mean radial position = ", meanRadialPosition, " +/- ", stdRadialPosition, " AU"
print " Mean total velocity = ", meanTotalVelocity, " +/- ", stdTotalVelocity, " AU/Year"
print "  Mean total energy = ", meanTotalEnergy, " +/- ", stdTotalEnergy, " Joules"

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
nApsideCalls = 0
nApsideFails = 0
debug = False

if(circularOrbit):
    orbitTime = 2*mp.pi*meanRadialPosition/meanTotalVelocity
    numberOfOrbits = int(maximumTime/orbitTime)

#
# Define two functions to do a quadratic fit for five points around the apsidal points of the orbit
#
def funcThreeParam(x, a, b, c):                # three parameters for a quadratic fit
    return a*x*x + b*x + c

def findApsidalTime(nTimeStepFirst):
    #
    # In astronomy the apsides are the perihelion or the apehilion
    # This function takes a triplet of (t,r) coordinates and assumes that they form a quadratic shape: r(t) = a*t*t + b*t + c around the apside
    # The coefficients (a,b,c) are extracting from a fitting program for five orbit points for which the derivatice dr/dt changes sign between the 2 and the third points
    # The time of a minimum or maximum r(t) is found by taking the derivative dr/dt = 2*at + b = 0,  which gives t = -b/2a for the time of the minimum or the maximum
    #
    global xNumerical, yNumerical, timeGrid, radialPosition, timeStep, nApsideCalls, nApsideFails
    nApsideCalls += 1
    nt0 = nTimeStepFirst
    timeValues = np.array([timeGrid[nt0], timeGrid[nt0+1], timeGrid[nt0+2], timeGrid[nt0+3], timeGrid[nt0+4]])
    radialValues = np.array([radialPosition[nt0], radialPosition[nt0+1], radialPosition[nt0+2], radialPosition[nt0+3], radialPosition[nt0+4]])
    kIndex = 0
    rError = []
    while kIndex < 5:
        rError.append(1.0e-4*radialPosition[nt0+kIndex])
        kIndex += 1

    radialErrors = np.array([rError[0], rError[1], rError[2], rError[3], rError[4]])
    initialCoefficients = np.array([1, 0, 0])
    fitResultsThreeParam = optimization.curve_fit(funcThreeParam, timeValues, radialValues, initialCoefficients, radialErrors)

    coefficients = fitResultsThreeParam[0]
    aCoeff = coefficients[0]
    bCoeff = coefficients[1]
    cCoeff = coefficients[2]
    apsideTime = -bCoeff/(2*aCoeff)
    radius = funcThreeParam(apsideTime, aCoeff, bCoeff, cCoeff)

    if(apsideTime < timeGrid[nt0] or apsideTime > timeGrid[nt0+4]):
        print "\n timeList ", timeValues
        print "  radial list ", radialValues
        print " Out of range apside time ", apsideTime, ",  apside r ", radius
        nApsideFails += 1
        apsideTime = timeGrid[nt0+2] - timeStep/2.

    kIndex = nt0
    kIndexLimit = nt0+5
    while kIndex < kIndexLimit:
        if(apsideTime < timeGrid[kIndex]):
            break
        kIndex += 1

    apsideTimeIndex = kIndex - 1
    timeFraction = (apsideTime - timeGrid[apsideTimeIndex])/timeStep

    theta1 = np.arctan2(yNumerical[apsideTimeIndex],xNumerical[apsideTimeIndex])
    if(theta1 < -mp.pi/2):
        theta1 += 2*mp.pi
    theta2 = np.arctan2(yNumerical[apsideTimeIndex+1],xNumerical[apsideTimeIndex+1])
    if(theta2 < -mp.pi/2):
        theta2 += 2*mp.pi
    angle = theta1 + timeFraction*(theta2 - theta1)
    angle = mp.degrees(angle)

    #
    #  Special code for when the angle is either near 0 or near 180, which occurs when alpha = 0
    #
    if(abs(angle) < 1.0e-02):
        angle = 0.0
    if(abs(angle - 180.0) < 1.0e-02):
        angle = 180.0

    if(debug):
        print "\n  timeList ", timeValues
        print "  radial list ", radialValues
        print "  nt0 ", nt0, ",  apsideTimeIndex ",apsideTimeIndex
        print "  apside time ", apsideTime, ",  apside r ", radius, ", timeFraction ", timeFraction
        print "  theta1 ", mp.degrees(theta1), ", theta2 ", mp.degrees(theta2), ",  angle ", angle

    return apsideTime, angle, radius

if(ellipticalOrbit):
    print "\n       Doing advanced pattern recognition analysis of the orbits"
    #
    # Algorithm is to accumulate the set of perihelion and aphelion radial values and their time values
    # The times to pass through successive perihelion points are stored from which a mean and standard deviation are computed
    # The first step is to confirm that the initial point is consistent with being a perihelion
    # The program input model is that the orbit starts as either a perihelion or an aphelion
    #
    radiusPerihelion = []
    radiusAphelion = []
    timePerihelionPoint = []
    timeAphelionPoint = []
    anglePerihelion = []
    angleAphelion = []
    if(radialPosition[1] > radialPosition[0] and totalVelocity[1] < totalVelocity[0]):
        radiusPerihelion.append(radialPosition[0])
        timePerihelionPoint.append(0)
        anglePerihelion.append(0)
        lastPerihelionAngle = 0.0
        lastPerihelionTime = 0.0
        lookingForNextAphelion = True
        lookingForNextPerihelion = False
        startAsPerihelion = True
        foundPerihelion = True
        foundAphelion = False
        print "  Check that the initial position is consistent with being a perihelion is passed"
    else:
        radiusAphelion.append(radialPosition[0])
        angleAphelion.append(0)
        lastAphelionAngle = 0.0
        lastAphelionTime = 0.0
        timeAphelionPoint.append(0)
        lookingForNextAphelion = False
        lookingForNextPerihelion = True
        startAsPerihelion = False
        foundPerihelion = False
        foundAphelion = True
        print "  Check that the initial position is consistent with being an aphelion is passed"

    deltaAnglePerihelion = []         # difference between current perihelion angle and previous perihelion angle angle
    perihelionAngleError = []
    deltaAngleAphelion = []           # difference between current aphelion angle and previous aphelion angle
    aphelionAngleError = []
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
                if(foundAphelion == False):
                    foundAphelion = True
                    apsideTimeAngle = findApsidalTime(nTimeStep - 2)    # find a precise time and angle for the aphelion position
                    lastAphelionTime = apsideTimeAngle[0]
                    angle = apsideTimeAngle[1]
                    lastAphelionAngle = angle   # first aphelion angle
                    angleAphelion.append(angle)
                    if(verbose):
                        print "\n Found first aphelion at time ", lastAphelionTime, " and angle ", angle
                else:

                    apsideTimeAngle = findApsidalTime(nTimeStep - 3)    # find a precise time and angle for the aphelion position
                    timeValue = apsideTimeAngle[0]
                    timeAphelion.append(timeValue - lastAphelionTime)     # store the time difference since the last aphelion
                    lastAphelionTime = timeValue
                    angle = apsideTimeAngle[1]
                    if(angle < lastAphelionAngle):
                        #
                        # This indicates that the precession angle has shifted one or more full orbits
                        #
                        while angle < lastAphelionAngle:
                            angle += 360.0
                    deltaAngleAphelion.append(angle - lastAphelionAngle)   # store in degrees
                    aphelionAngleError.append(precessionAngleError)
                    lastAphelionAngle = angle
                    angleAphelion.append(angle)
                    radiusAphelion.append(apsideTimeAngle[2])
                    timeAphelionPoint.append(timeValue)
                    if(verbose):
                        print "\n Found next aphelion at time ", lastAphelionTime, ", angle ", angle, ", and radius ", apsideTimeAngle[2]
                        print " Aphelion r1-r5 ", radialPosition[nTimeStep - 3],radialPosition[nTimeStep - 2],radialPosition[nTimeStep -1],radialPosition[nTimeStep],radialPosition[nTimeStep + 1]
                        print " Aphelion t1-t5 ", timeGrid[nTimeStep - 3],timeGrid[nTimeStep - 2],timeGrid[nTimeStep - 1],timeGrid[nTimeStep],timeGrid[nTimeStep + 1]

                lookingForNextPerihelion = True
                lookingForNextAphelion = False

            lastRadialPosition = newRadialPosition

        if(lookingForNextPerihelion):
            #
            # If the new radial position is greater than the previous radial position, then the next perihelion point has been crossed
            #
            if(newRadialPosition > lastRadialPosition):
                if(foundPerihelion == False):
                    foundPerihelion = True
                    apsideTimeAngle = findApsidalTime(nTimeStep - 2)    # find a precise time and angle for the perihelion position
                    lastPerihelionTime = apsideTimeAngle[0]
                    angle = apsideTimeAngle[1]
                    lastPerihelionAngle = angle   # first perihelion angle
                    anglePerihelion.append(angle)
                    if(verbose):
                        print "\n Found first perihelion at time ", lastPerihelionTime, ", angle ", angle
                else:

                    apsideTimeAngle = findApsidalTime(nTimeStep - 2)    # find a precise time and angle for the perihelion position
                    timeValue = apsideTimeAngle[0]
                    timePerihelion.append(timeValue - lastPerihelionTime)     # store the time difference since the last aphelion
                    lastPerihelionTime = timeValue
                    angle = apsideTimeAngle[1]
                    if(angle != 0.0 and angle < lastPerihelionAngle):
                        #
                        # This indicates that the precession angle has shifted one or more full orbits
                        #
                        while angle < lastPerihelionAngle:
                            angle += 360.0
                    deltaAnglePerihelion.append(angle - lastPerihelionAngle)   # store in degrees
                    perihelionAngleError.append(precessionAngleError)
                    lastPerihelionAngle = angle
                    anglePerihelion.append(angle)
                    radiusPerihelion.append(apsideTimeAngle[2])
                    timePerihelionPoint.append(timeValue)
                    if(verbose):
                        print "\n Found next perihelion at time ", lastPerihelionTime, ", angle ", angle, ", and radius ", apsideTimeAngle[2]
                        print " Perihelion r1-r5 ", radialPosition[nTimeStep - 3],radialPosition[nTimeStep - 2],radialPosition[nTimeStep -1],radialPosition[nTimeStep],radialPosition[nTimeStep + 1]
                        print " Perihelion t1-t5 ", timeGrid[nTimeStep - 3],timeGrid[nTimeStep - 2],timeGrid[nTimeStep - 1],timeGrid[nTimeStep],timeGrid[nTimeStep + 1]
                lookingForNextPerihelion = False
                lookingForNextAphelion = True

            lastRadialPosition = newRadialPosition

        nTimeStep += 1

    if(verbose):
        print "\n Number of perihelion angles ", len(anglePerihelion)
        print anglePerihelion

        print "\n Number of aphelion angles ", len(angleAphelion)
        print angleAphelion

    nPerihelion = len(radiusPerihelion)
    nAphelion = len(radiusAphelion)
    if(nPerihelion > 0 and nAphelion > 0):       # require at least one found perihelion and one found aphelion in order to quote orbit information
        majorAxisLength = np.mean(radiusAphelion) + np.mean(radiusPerihelion)
        eccentricity = (np.mean(radiusAphelion) - np.mean(radiusPerihelion))/majorAxisLength
        print '%s %5.3f' % ('  The orbit eccentricity = ', eccentricity)

        nPerihelionTimeDiff = len(timePerihelion)
        nAphelionTimeDiff = len(timeAphelion)
        meanTimeBetweenAphelion = np.mean(timeAphelion)
        meanTimeBetweenPerihelion = np.mean(timePerihelion)
        if(verbose):
            print "\n  Number of aphelion points = ", nAphelion, ", number of perihelion points ", nPerihelion
            print "  Mean time between aphelion points = ", meanTimeBetweenAphelion
            print "  Mean time between perihelion points = ", meanTimeBetweenPerihelion
        orbitTime = (nPerihelionTimeDiff*meanTimeBetweenPerihelion + nAphelionTimeDiff*meanTimeBetweenAphelion)/float(nPerihelionTimeDiff + nAphelionTimeDiff)
        if(verbose):
            print "  Average orbit time = ", orbitTime
            kPrint = 0
            nPrint = min(nPerihelion, nAphelion)
            while kPrint < nPrint:
                print "  ", kPrint, ")  perihelion angle change ", deltaAnglePerihelion[kPrint], ",  aphelion angle change ", deltaAngleAphelion[kPrint]
                kPrint += 1

        perihelionPrecessMean = np.mean(deltaAnglePerihelion)/orbitTime            # precession rate
        perihelionPrecessStdDev = np.std(deltaAnglePerihelion)/orbitTime
        aphelionPrecessMean = np.mean(deltaAngleAphelion)/orbitTime
        aphelionPrecessStdDev = np.std(deltaAngleAphelion)/orbitTime
        print "\n  For alpha = ", alpha, " the mean value of perihelion precession rate = ", perihelionPrecessMean, " degrees/year with a standard deviation ", perihelionPrecessStdDev
        print "  For alpha = ", alpha, " the mean value of aphelion precession rate = ", aphelionPrecessMean, " degrees/year with a standard deviation ", aphelionPrecessStdDev
        uncertainty = (perihelionPrecessStdDev + aphelionPrecessStdDev)/2.0
        print "  For alpha = ", alpha, " the average precession = ", (perihelionPrecessMean + aphelionPrecessMean)/2, " degrees/year with an uncertainty ", uncertainty

        if(startAsPerihelion):
            numberOfOrbits = nPerihelion      # the number of orbits is one more than the number of differences between successive perihelion points
        else:
            numberOfOrbits = nAphelion        # the number of orbits is one more than the number of differences between successive aphelion points

print '%s %d' % ('  Number of completed orbits = ', numberOfOrbits)
if(numberOfOrbits > 0):
    print '%s %5.3f %s' % ('  Orbit time = ', orbitTime, ' years')

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

alphaString = 'Alpha = ' + str(alpha)
massString = 'Mass = ' + str(mass*1.0e24) + ' kg'
timeString = 'Total time = ' + str(maximumTime) + ' years in steps of ' + str(timeStep) + ' years'
orbitTimeString = 'Orbital time = ' + str(orbitTime) + ' years'
eccentricityString = 'Orbit eccentricity = ' + str(eccentricity)

plt.figure(1)       # start a figure for a single plot of the orbit
#
# Check on the number of data points compared to the maximum number of points to be plotted
# Having a huge number of data points to be plotted can take extra time, or even cause software library errors
# Also check the minimum and maximum plot times
#
print "\n Orbit plot begins at ", minPlotT, " years and ends at ", maxPlotT, " years"
if(nTimeSteps < maxPointsToPlot and minPlotT == 0 and maxPlotT == maximumTime):
    if(pointsPlot):
        plt.plot(xNumerical, yNumerical, 'bo')    # plot as discrete data points
    else:
        plt.plot(xNumerical, yNumerical)    # plot as a continuous line
else:
    xplot = []
    yplot = []
    skipTimeInterval = nTimeSteps/maxPointsToPlot
    if(skipTimeInterval < 0):
        print "\n Program Error, skipTimeInterval ", skipTimeInterval
        exit()
    if(skipTimeInterval < 1):
        skipTimeInterval = 1
    else:
        if(verbose):
            print "\n Plotting skip time value ", skipTimeInterval*timeStep, " years"
    nTimeStep = 0

    while nTimeStep < nTimeSteps:
        time = timeGrid[nTimeStep]
        if(time >= minPlotT and time <= maxPlotT):
            xplot.append(xNumerical[nTimeStep])
            yplot.append(yNumerical[nTimeStep])
        nTimeStep += skipTimeInterval

    if(pointsPlot):
        plt.plot(xplot, yplot, 'bo')    # plot as discrete data points
    else:
        plt.plot(xplot, yplot)    # plot as a continuous line

#
# Draw lines for the perihelion points
#
nLines = len(timePerihelionPoint)
kLine = 0
if(verbose):
    print "\n List of perhelion times and coordinate positions within the plotting range"
    print "  Number of apside calls ", nApsideCalls, ", with number of curve fit out-of-range ", nApsideFails
while kLine < nLines:
    time = timePerihelionPoint[kLine]
    nTimeStep = int(time/timeStep)
    if(time >= minPlotT and time <= maxPlotT):
        if(nTimeStep == 0):
            xValue = xNumerical[0]          # orbit start is a perihelion
            yValue = yNumerical[0]          # orbit start is a perihelion
        else:
            angle = anglePerihelion[kLine]
            radius = radiusPerihelion[kLine]
            angleRadians = angle*mp.pi/180.
            xValue = radius*mp.cos(angleRadians)
            yValue = radius*mp.sin(angleRadians)
        if(verbose):
            print " ", kLine, ") time ", time, "  xValue ", xValue, "  yValue ", yValue, "  nTimeStep ", nTimeStep
        xPeriLine = [0, xValue]
        yPeriLine = [0, yValue]
        plt.plot(xPeriLine, yPeriLine, 'g')   # line from the Sun
    kLine += 1

xTextPosition = minimumPositionX + 0.25*(maximumPositionX - minimumPositionX)
yTextPosition = minimumPositionY + 0.55*(maximumPositionY - minimumPositionY)

plt.text(xTextPosition, minimumPositionY + 0.70*(maximumPositionY - minimumPositionY), alphaString)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.63*(maximumPositionY - minimumPositionY), x0String)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.56*(maximumPositionY - minimumPositionY), vy0String)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.44*(maximumPositionY - minimumPositionY), massString)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.37*(maximumPositionY - minimumPositionY), timeString)   # text to document the parameters used
plt.text(xTextPosition, minimumPositionY + 0.30*(maximumPositionY - minimumPositionY), orbitShapeString)   # text to document the parameters used
if(numberOfOrbits > 0):
    plt.text(xTextPosition, minimumPositionY + 0.23*(maximumPositionY - minimumPositionY), orbitTimeString)   # text to document the parameters used
    plt.text(xTextPosition, minimumPositionY + 0.16*(maximumPositionY - minimumPositionY), eccentricityString)   # text to document the parameters used

plt.xlabel('x Coordinate (AU)')                             # add axis labels
plt.ylabel('y Coordinate (AU)')

xSun = []
ySun = []
xSun.append(0)
ySun.append(0)
plt.scatter(0, 0, c='y', s=200)

plt.title('Orbit of a Single Planet in the Solar System with Precession')

plt.grid(True)
plt.xlim(minimumPositionX, maximumPositionX)
plt.ylim(minimumPositionY, maximumPositionY)

plt.show()          # show the complete figure with the upper and lower subplots
