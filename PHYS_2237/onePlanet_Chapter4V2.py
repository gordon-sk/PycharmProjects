#
# Program 4.2 Solution for the orbit of a planet in the solar system with the Sun as fixed force center (onePlanet_Chapter4V2.py)
#
#              Give the command  python onePlanet_Chapter4V2.py -h  to get help command on input parameters
#
#              This program is the same as the onePlanet_Chapter4V1.py program except that it examines the orbital data results in more detail, which is called "pattern recognition"
#              This examination extracts the period of the orbit and the eccentricity of the orbit
#
#              The input arguments to this program are the same as for the onePlanet_Chapter4V1.py program
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
#              Defaults are x0 = 0.983 AU, vy0 = 1, mass (5.98, Earth's mass), timeStep = 0.02 years, maximumTime = 3 years
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
parser.add_argument('--x0', default=0.586, type=float, help="Initial perihelion x0 in AU; default 0.983 AU")
parser.add_argument('--vy0', default=1.832, type=float, help="Initial perihelion vy0 in units of 2Pi AU/Y; default 1.017")
parser.add_argument('--maxT', default=3.0*80, type=float, help="Maximum time range in years; default 1 Year")
parser.add_argument('--deltaT', default=0.005, type=float, help="Iteration time step in years; default 0.02 Year")
parser.add_argument('--mass', default=5.98, type=float, help="Planet mass in units of 10^24 kg; default 5.98")
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
x0 = args.x0
if(x0 <= 0.0):
    print "\n Cannot have the initial position", x0, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(0)

vy0 = args.vy0
if(vy0 <= 0.0):
    print "\n Cannot have the initial velocity omponent", vy0, " be less than or equal to zero"
    print "\n Program is exiting\n"
    exit(0)

maximumTime = args.maxT
timeStep = args.deltaT
mass = args.mass        # Note that the mass of the planet does not affect its orbit size or orbital period, as per Kepler's Third Law
masskg = mass*1.e24

pointsPlot = args.pointsPlot
verbose = args.verbose
maxPointsToPlot = args.maxPointsToPlot

vy0Speed = vy0*2*mp.pi  # in units of AU/Year

# dx/dt = vx    velocity component in the x direction
# dy/dt = vy    velocity component in the y direction
# dvx/dt = ax = -4*pi*pi*x/r^3  acceleration component in the x direction from the universal gravity force
# dvy/dt = ay = -4*pi*pi*y/r^3  acceleration component in the y direction from the universal gravity force

fourPiSquare = 4*mp.pi*mp.pi                       # simple universal gravity force scaling factor when astronomical units are used (note the planet mass value is not used)
def fDerivative(variableList, t):                  # variableList dummy list array since there is more than one differential equation
    vx = variableList[0]                           # speed in the x direction
    vy = variableList[1]                           # speed in the y direction
    x = variableList[2]                            # x coordinate
    y = variableList[3]                            # y coordinate
    rDistance = mp.sqrt(x*x + y*y)                 # radial distance from the force center
    rDistanceCubed = rDistance*rDistance*rDistance # cube of the radial distance
    dvxdt = -fourPiSquare*x/rDistanceCubed         # the time derivative of velocity in the x direction according to the universal gravity force component
    dvydt = -fourPiSquare*y/rDistanceCubed         # the time derivative of postion in the y direction is the universal gravity force component
    return [dvxdt, dvydt, vx, vy]                  # return the four derivatives as a list object containing four elements

print "\n Orbit of a planet around the Sun as a fixed force center"
print "         Input conditions"
print "  Initial radius ", x0, " AU, at an initial transverse velocity of ", vy0Speed, " AU/Year"
print "  Time step = ", timeStep, " year,  maximum time range = ", maximumTime, " years"

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

print "\n       Doing advanced pattern recognition analysis of the orbits"
eccentricity = 0
orbitTime = 0
numberOfOrbits = 0
if(circularOrbit):
    orbitTime = 2*mp.pi*meanRadialPosition/meanTotalVelocity
    numberOfOrbits = int(maximumTime/orbitTime)

if(ellipticalOrbit):
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

print '%s %d' % ('  Number of completed orbits = ', numberOfOrbits)
if(numberOfOrbits > 0):
    keplerThirdLawTime = pow(majorAxisLength/2.0, 1.5)
    print '%s %5.3f %s' % ('  Orbit time determined from pattern recognition = ', orbitTime, ' year')
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
vy0String = 'Initial transverse velocity = ' + str(vy0) + '*2Pi AU/Y'
massString = 'Mass = ' + str(mass*1.0e24) + ' kg'
timeString = 'Total time = ' + str(maximumTime) + ' years in steps of ' + str(timeStep) + ' years'
orbitTimeString = 'Orbital time = ' + str(orbitTime) + ' years'
eccentricityString = 'Orbit eccentricity = ' + str(eccentricity)

plt.figure(1)       # start a figure for a single plot of the orbit
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

plt.title('Orbit of a Single Planet in the Solar System')

plt.grid(True)
plt.xlim(minimumPositionX, maximumPositionX)
plt.ylim(minimumPositionY, maximumPositionY)

plt.show()          # show the complete figure with the upper and lower subplots
