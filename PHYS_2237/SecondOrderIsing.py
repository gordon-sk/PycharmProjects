#  Program 8.3  Numerical Solution to the Second Order Ising Model for Ferromagnetism
#  Same as V1 version but adapted to do Exercises 8.3 and 8.8
#
#  Prints two output files
#    1) <s> as a function of T:  T, s1, s2, s3 where s1, s2, and s3 are from three separate 1000 steps
#    2) E as a function of T:  T, E-Mean, (Delta E)^2  where E-Mean is for all 3000 steps and where (Delta E)^2 is the variance of E
#                                                      (Delta E)^2 = sum(E-Mean - Ei)^2/timeSteps,  i = 1 to 3000
#
#  Second order Ising model calculation, Section 8.4
#  Producing Figure 8.6 page 249, and Figure 8.7 page 251
#  Producing (not quite) Figure 8.9 page 254, Figure 8.10 on page 256
#  This version prints out (M,T) values between T = 2.00 and 2.50
#  This version allows for input of number of time steps
#
#  Output plots
#  1) Four sets of <s> vs time, for temperatures 1.5, 2.0, 2.25, and 4.0 [temperature index 2, 4, 5, and 12]
#
import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible
import matplotlib.pyplot as plt     # get matplotlib plot functions
import  argparse                    # argument parser library
import numpy as np
import random                       # random number generator library
from datetime import datetime

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()

parser.add_argument('--isingArray', default=40, type=int, help="Dimension of Ising lattice; default 10")
parser.add_argument('--startTemperature', default=2.00, type=float, help="Starting temperature, units of J/kB; default 1.00")
parser.add_argument('--finalTemperature', default=2.60, type=float, help="Final temperature units of J/kB; default 5.00")
parser.add_argument('--deltaTemperature', default=0.02, type=float, help="Temperature increment; default 0.10")
parser.add_argument('--pairEnergy', default=1.0, type=float, help="Pair energy; default 1.00")
parser.add_argument('--seedFixed', action='store_true', help="Use a fixed random number seed; default False")
parser.add_argument('--timeSteps', default=3000, type=int, help="Number of time steps; default 3000")
parser.add_argument('--allSpinsRandom', action='store_true', help="Initially all spins point randomly up or down; default False")
parser.add_argument('--printTimes', action='store_true', help="Print elapsed wall clock time at each temperature; default False")
parser.add_argument('--verbose', action='store_true', help="Extra print out; default False")
parser.add_argument('--debug', action='store_true', help="Debugging print out; default False")
parser.add_argument('--skipNegativeCheck', action='store_true', help="Skip check against complete flip to negative; default False")
parser.add_argument('--triggerNegative', default=-0.60, type=float, help="Check against complete flip to negative; default -0.75")
args = parser.parse_args()

isingArray = args.isingArray
isingParticles = isingArray*isingArray
isingParticlesInverse = 1.0/float(isingParticles)
isingParticlesMinusOneInverseSquareRoot = 1.0/np.sqrt(float(isingParticles) - 1.0)


startTemperature = args.startTemperature
finalTemperature = args.finalTemperature
deltaTemperature= args.deltaTemperature
temperatureSteps = int((finalTemperature - startTemperature)/deltaTemperature) + 1
indexPlotTemperature = [2, 4, 5, 12]    # these index values will change if any of the temperature paramters are changed

pairEnergy = args.pairEnergy
flipEnergy = 2*pairEnergy

timeSteps = args.timeSteps
seedFixed = args.seedFixed
allSpinsRandom = args.allSpinsRandom
verbose = args.verbose

skipNegativeCheck = args.skipNegativeCheck
negativeCheck = False
if(skipNegativeCheck == False):
    negativeCheck = True
    triggerNegative = args.triggerNegative

debug = args.debug
seedFixed = True
if(seedFixed):
    random.seed(3)  # fixes the sequence of random sees for repeatability of numerical results
else:
    random.seed()   # random seed is set according to the system clock, numerical results will be non-repeatable

print "\n  Solving the Second Order Ising Model with a 2-D Lattice Dimenson ", isingArray, " x ", isingArray
print "  Start temperature ", startTemperature, ",  final temperature ", finalTemperature, " in ", temperatureSteps, " steps of size ", deltaTemperature
print "  Using ", timeSteps, " time steps"
if(allSpinsRandom):
    print "  All spins are randomly oriented +1 or -1 at the start of each temperature cycle"
else:
    print "  All spins are oriented +1 at the start of each temperature cycle"

if(negativeCheck):
    print "  Negative flip check triggers at ", triggerNegative

arraySpins = [[0 for jColumn in range(isingArray)] for iRow in range(isingArray)]  # spin value for the lattice element [iRow][jColumn]

#
# The spin array values at each time step, re-calculated at each temperature step: spinArray[timeStep][iRow][jColumn]
#
spinsAllSteps = [[[0 for jColumn in range(isingArray)] for iRow in range(isingArray)] for timeStep in range(timeSteps)]

#
# Define the two dimensional arrays for spin average at each step sValue[timeStep][tempStep]
# Do like wise for timeValue[timeStep][tempStep], energyValue[timeStep][tempStep], varianceValue[timeStep][tempStep], correlationValue]timeStep][tempStep]
#
sValue = [[0.0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]
timeValue = [[0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]
energyValue = [[0.0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]
stdOfTheMeanValue = [[0.0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]
correlationValue = [[0.0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]

energyMeanAllTemperatures = [0 for tempStep in range(temperatureSteps)]
energyStdAllTemperatures = [0 for tempStep in range(temperatureSteps)]
heatCapacityAllTemperatures = [0 for tempStep in range(temperatureSteps)]

#
# array for the spin products si*sj
#
spinProductArray = [0 for iParticle in range(isingParticles)]

#
# Initialize the arrays which will store the nearest neighbor rows and columns for a given iRow or jColumn
#
adjacentUpperRow = [0 for iRow in range(isingArray)]
adjacentLowerRow = [0 for iRow in range(isingArray)]
adjacentLeftColumn = [0 for jColumn in range(isingArray)]
adjacentRightColumn = [0 for jColumn in range(isingArray)]
isingArrayMinusOne = isingArray - 1  # used to get the adjacent rows or columns on the other side of the array ("teleportation")

for iRow in range(isingArray):
    #
    # upper row neighbors
    #
    if(iRow==0):
        iRowUpper = isingArrayMinusOne  # topmost row looking at the bottommost row of spins
    else:
        iRowUpper = iRow - 1            # non-top row looking at the row above
    adjacentUpperRow[iRow] = iRowUpper

    #
    # lower row neighbors
    #
    if(iRow == isingArrayMinusOne):
        iRowLower = 0   # interact with the topmost row
    else:
        iRowLower = iRow + 1
    adjacentLowerRow[iRow] = iRowLower

for jColumn in range(isingArray):
    #
    # left column neighbor
    #
    if(jColumn==0):
        jColumnLeft = isingArrayMinusOne     # leftmost row looking at the rightmost column of spins
    else:
        jColumnLeft = jColumn - 1
    adjacentLeftColumn[jColumn] = jColumnLeft

    #
    # right column neighbor
    #
    if(jColumn == isingArrayMinusOne):
        jColumnRight = 0     # interact with the leftmost column
    else:
        jColumnRight = jColumn + 1
    adjacentRightColumn[jColumn] = jColumnRight

def initializeSpins(isingArray, arraySpins):
    averageSpin = 0
    for iRow in range(isingArray):
        for jColumn in range(isingArray):
            thisSpin = 1
            if(allSpinsRandom):
                thisSpin = random.randrange(-1,3,2)    # returns either -1 or +1 randomly, note that the 3 does not get used

            arraySpins[iRow][jColumn] = thisSpin    # returns either -1 or +1 randomly, note that the 3 does not get used
            averageSpin += thisSpin

    return averageSpin*isingParticlesInverse

initialSpinAverage = initializeSpins(isingArray, arraySpins)
print "  Initial average spin ", initialSpinAverage

def printSpinArray(isingArray, arraySpins):
    for iRow in range(isingArray):
        print "  ",      # leading two spaces for each row in the print
        for jColumn in range(isingArray):
            thisSpinString = '%2d' % (arraySpins[iRow][jColumn])
            print thisSpinString,

        print " "  # carriage return/end of line for a given row

    return

printSpinArray(isingArray, arraySpins)

def calculateEnergy(isingArray, arraySpins): # calculates energy
    #
    # calculates the total energy of the system
    #
    iParticle = 0
    spinProduct = 0
    for iRow in range(isingArray):

        iRowUpper = adjacentUpperRow[iRow]
        iRowLower = adjacentLowerRow[iRow]

        for jColumn in range(isingArray):
            thisSpin = arraySpins[iRow][jColumn]
            jColumnLeft = adjacentLeftColumn[jColumn]
            jColumnRight = adjacentRightColumn[jColumn]

            thisProduct = thisSpin*(arraySpins[iRowUpper][jColumn]+arraySpins[iRowLower][jColumn]+arraySpins[iRow][jColumnLeft]+arraySpins[iRow][jColumnRight])
            spinProductArray[iParticle] = thisProduct
            iParticle += 1
            spinProduct += thisProduct

    meanValueNearestNeighborSpins = np.mean(spinProductArray)
    stdMeanNearestNeighborSpins = np.std(spinProductArray)*isingParticlesMinusOneInverseSquareRoot
    ratioStdOverMean = 0.0
    if(meanValueNearestNeighborSpins != 0.0):
        ratioStdOverMean = stdMeanNearestNeighborSpins/meanValueNearestNeighborSpins

    return -0.5*pairEnergy*spinProduct/isingParticles, ratioStdOverMean  #  take half the energy because of double pair counting

initialEnergyAverage, ratioStdOverMean = calculateEnergy(isingArray, arraySpins)
print "  Initial energy ", initialEnergyAverage, " mean value per spin pair, with a standard deviation of the mean", initialEnergyAverage*ratioStdOverMean

def updateSpin(rowPosition, columnPosition, isingArray, arraySpins, temperature):
    #
    # decide whether to flip the spin location (rowPosition, columnPosition)
    #

    eFlip = 0.0;
    thisSpin = arraySpins[rowPosition][columnPosition]

    #
    # upper row neighbor
    #
    iRowUpper = adjacentUpperRow[rowPosition]
    if(thisSpin*arraySpins[iRowUpper][columnPosition] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy

    #
    # lower row neighbor
    #
    iRowLower = adjacentLowerRow[rowPosition]
    if(thisSpin*arraySpins[iRowLower][columnPosition] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy

    #
    # left column neighbor
    #
    iColumnLeft = adjacentLeftColumn[columnPosition]
    if(thisSpin*arraySpins[rowPosition][iColumnLeft] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy

    #
    # right column neighbor
    #
    iColumnRight = adjacentRightColumn[columnPosition]
    if(thisSpin*arraySpins[rowPosition][iColumnRight] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy

    #
    # Based on the energy calculated, determine if this spin should be flipped to give a lower energy state
    #

    if(eFlip <= 0.0):
        arraySpins[rowPosition][columnPosition] = -thisSpin    # always flip if energy will be lower after flip

    if(eFlip > 0.0):
        rndm = random.uniform(0, 1.0)  #  random number [0,1)
        expTest = np.exp(-eFlip/temperature)
        if(rndm <= expTest):
            arraySpins[rowPosition][columnPosition] = -thisSpin;  # flip the spin randomly according to temperture

    return

def solveSecondOrder(isingArray, arraySpins, temperature):
    #
    # function to compute average spin for this temperature at this time
    #

    #
    # sweep over all the spins in the square lattice
    #
    for iRow in range(isingArray):
        for jColumn in range(isingArray):
            #
            # check if this spin should be flipped
            #
            updateSpin(iRow, jColumn, isingArray, arraySpins, temperature)

    solution = 0  # sum of all the spin values
    for iRow in range(isingArray):
        for jColumn in range(isingArray):
            solution += arraySpins[iRow][jColumn]

    #
    # return the average spin value
    #
    return float(solution)*isingParticlesInverse

firstThirdTime = int(timeSteps/3)
secondThirdTime = 2*firstThirdTime
inverseSquareRootTimeStepsMinusOne = 1.0/np.sqrt(int(timeSteps)-1.0)

temperature = startTemperature
totalSeconds = 0.0
startDateTimeString = str(datetime.now())
print "  Begin cycling over temperatures at ", startDateTimeString
sSolutionSum = 0.0
spinMean1 = []
spinMean2 = []
spinMean3 = []
energyMean1 = []
energyMean2 = []
energyMean3 = []
cMaximum = 0.0
temperatureMean = []
temperatureMaximumC = 0.0
for temperatureStep in range(temperatureSteps):
    t1 = datetime.now()
    oneThousandCycle = -1
    countNegative = 0
    initializeSpins(isingArray, arraySpins)  # initialize the spins at initial time step and at every 1,000 time steps
    for timeStep in range(timeSteps):
        sSolution = solveSecondOrder(isingArray, arraySpins, temperature)
        timeValue[timeStep][temperatureStep] = timeStep
        sValue[timeStep][temperatureStep] = sSolution         # store spin value solution at this time step and this temperture
        if(sSolution < triggerNegative):
            countNegative += 1

        energy, stdOfTheMean = calculateEnergy(isingArray, arraySpins)
        energyValue[timeStep][temperatureStep] = energy # store energy value
        stdOfTheMeanValue[timeStep][temperatureStep] = stdOfTheMean
        #
        # Check results every 1000 time steps
        #
        oneThousandTimeCheck = 1000*int(timeStep/1000) - timeStep
        if(oneThousandTimeCheck == 0):
            oneThousandCycle += 1
            indexBase = (oneThousandCycle-1)*1000
            if(debug and temperature == 1.5): print "  **** 1000 cycle ", oneThousandCycle, ", at time step ", timeStep, ",  countNegative ", countNegative, ",  indexBase ", indexBase, "  ***"
            if(oneThousandCycle > 0 and negativeCheck and countNegative > 400 and temperature < 2.27):
                if(debug and temperature == 1.5): print "      Flipping spins with indexBase ", indexBase, ", at time step ", timeStep
                for index in range(1000):
                    sSolution = sValue[indexBase + index][temperatureStep]
                    if(sSolution < triggerNegative):
                        sValue[indexBase + index][temperatureStep] = -sSolution

            countNegative = 0
            initializeSpins(isingArray, arraySpins)  # initialize the spins at initial time step and at every 1,000 time steps

    #
    #  Check for negative flip on last 1000 time steps
    #
    if(negativeCheck and countNegative > 400 and temperature < 2.27):
        indexBase += 1000
        if(debug and temperature == 1.5): print "      Flipping spins with indexBase ", indexBase, " after last time step with countNegative ", countNegative
        for index in range(1000):
            sSolution = sValue[indexBase + index][temperatureStep]
            if(sSolution < triggerNegative):
                sValue[indexBase + index][temperatureStep] = -sSolution

    spinAverageThisTemperature = [sValue[timeStep][temperatureStep] for timeStep in range(timeSteps)]
    spinAverageThisTemperatureFirstThird = [sValue[timeStep][temperatureStep] for timeStep in range(firstThirdTime)]
    spinAverageThisTemperatureSecondThird = [sValue[timeStep+firstThirdTime][temperatureStep] for timeStep in range(firstThirdTime)]
    spinAverageThisTemperatureThirdThird = [sValue[timeStep+secondThirdTime][temperatureStep] for timeStep in range(firstThirdTime)]

    spinAverageThisTemperatureMean = np.mean(spinAverageThisTemperature)
    spinAverageThisTemperatureStdMean = np.std(spinAverageThisTemperature)*inverseSquareRootTimeStepsMinusOne
    spinMean1.append(np.mean(spinAverageThisTemperatureFirstThird))
    spinMean2.append(np.mean(spinAverageThisTemperatureSecondThird))
    spinMean3.append(np.mean(spinAverageThisTemperatureThirdThird))

    energyThisTemperature = [energyValue[timeStep][temperatureStep] for timeStep in range(timeSteps)]
    energyThisTemperatureMean = np.mean(energyThisTemperature)
    energyMeanAllTemperatures[temperatureStep] = energyThisTemperatureMean
    squareEnergy = 0.0
    for timeStep in range(timeSteps):
        thisEnergy = energyValue[timeStep][temperatureStep]
        squareEnergy += thisEnergy*thisEnergy
    squareEnergy = squareEnergy/float(timeSteps)
    squareDeltaE = squareEnergy - energyThisTemperatureMean*energyThisTemperatureMean
    energyThisTemperatureStd = np.std(energyThisTemperature)/np.sqrt(float(timeSteps - 1))
    energyStdAllTemperatures[temperatureStep] = energyThisTemperatureStd
    heatCapacityAllTemperatures[temperatureStep] = isingParticles*squareDeltaE/(temperature*temperature)
    if(heatCapacityAllTemperatures[temperatureStep] > cMaximum):
        cMaximum = heatCapacityAllTemperatures[temperatureStep]
        temperatureMaximumC = temperature

    energyThisTemperatureFirstThird = [energyValue[timeStep][temperatureStep] for timeStep in range(firstThirdTime)]
    energyThisTemperatureSecondThird = [energyValue[timeStep+firstThirdTime][temperatureStep] for timeStep in range(firstThirdTime)]
    energyThisTemperatureThirdThird = [energyValue[timeStep+secondThirdTime][temperatureStep] for timeStep in range(firstThirdTime)]
    energyMean1.append(np.mean(energyThisTemperatureFirstThird))
    energyMean2.append(np.mean(energyThisTemperatureSecondThird))
    energyMean3.append(np.mean(energyThisTemperatureThirdThird))

    temperatureMean.append(temperature)

    t2 = datetime.now()
    delta = t2 - t1
    totalSeconds += delta.total_seconds()
    if(verbose):
        resultsString = '%s %.4f %s %.1f %s %.4f %s %.3e %s %.4f %s %.4f, %s %.4f' % ('  At T =', temperature, ', total elapsed time', totalSeconds, 'seconds;  <s> =', spinAverageThisTemperatureMean,
                                                                     '+/- ', spinAverageThisTemperatureStdMean, ', E =',  energyThisTemperatureMean, '+/-', energyThisTemperatureStd, ' C = ', heatCapacityAllTemperatures[temperatureStep])
        print resultsString

    temperature += deltaTemperature

endDateTimeString = str(datetime.now())
print "\n  End cycling over temperatures at ", endDateTimeString
print "\n\n\n"

fig = plt.figure(1)       # start a figure
plt.scatter(temperatureMean, spinMean1, s=6, c='b', marker='o')
plt.scatter(temperatureMean, spinMean2, s=6, c='r', marker='o')
plt.scatter(temperatureMean, spinMean3, s=6, c='g', marker='o')
xMaximum = finalTemperature
xMinimum = startTemperature
plt.xlim(xMinimum, xMaximum)
plt.ylim(-0.25, 1.25)
plt.ylabel('<s> At Each Temperature')
plt.xlabel('Temperature (units of  J/kB)')
xPositionText = xMinimum + 0.52*(xMaximum - xMinimum)
yPositionText = 0.80
latticeString = '%s %d %s %d' % ('Lattice size ', isingArray, ' x ', isingArray)
plt.text(xPositionText, yPositionText, latticeString)
plt.grid(True)
plt.title('Second Order Ising Spin As a Function of Temperature')
plotNameString = 'plotDenseSvsTemperatureIsing' +  str(isingParticles) + 'DeltaT' + str(int(1000.*deltaTemperature)) + 'Milli.pdf'
fig.savefig(plotNameString)   # save the figure to file
print "   Second Order Ising Spin As a Function of Temperature"
plt.close(fig)

print "\n  Producing Energy vs temperature plot"

fig = plt.figure(2)       # start a figure
plt.scatter(temperatureMean, energyMean1, s=6, c='b', marker='o')
plt.scatter(temperatureMean, energyMean2, s=6, c='r', marker='o')
plt.scatter(temperatureMean, energyMean3, s=6, c='g', marker='o')
xMaximum = finalTemperature
xMinimum = startTemperature
plt.xlim(xMinimum, xMaximum)
plt.ylim(-2.2, 0.0)
plt.ylabel('Energy per particle At Each Temperature')
plt.xlabel('Temperature (units of  J/kB)')
xPositionText = xMinimum + 0.21*(xMaximum - xMinimum)
yPositionText = -0.40
latticeString = '%s %d %s %d' % ('Lattice size ', isingArray, ' x ', isingArray)
plt.text(xPositionText, yPositionText, latticeString)
plt.grid(True)
plt.title('Energy Per Spin Particle As a Function of Temperature')
plotNameString = 'plotDenseEvsTemperatureIsing' +  str(isingParticles) + 'DeltaT' + str(int(1000.*deltaTemperature)) + 'Milli.pdf'
fig.savefig(plotNameString)   # save the figure to file
print "   Energy Per Spin Particle As a Function of Temperature"
plt.close(fig)

print "\n  Producing Heat Capacity C vs temperature plot"
fig = plt.figure(3)       # start a figure
plt.scatter(temperatureMean, heatCapacityAllTemperatures, s=4, c='r', marker='o')
extraRange = 0.05*(finalTemperature - startTemperature)
xMaximum = finalTemperature + extraRange
xMinimum = startTemperature - extraRange
plt.xlim(xMinimum, xMaximum)
yMaximum = 1.5*max(heatCapacityAllTemperatures)
plt.ylim(0.0, yMaximum)
plt.ylabel('Heat Capacity C')
plt.xlabel('Temperature (units of  J/kB)')
xPositionText = xMinimum + 0.05*(xMaximum - xMinimum)
yPositionText = 0.90*yMaximum
latticeString = '%s %.3f %s %.3f %s %d %s %d' % ('Maximum C ', cMaximum, ' at T ', temperatureMaximumC, ' with lattice', isingArray, ' x ', isingArray)
plt.text(xPositionText, yPositionText, latticeString)
plt.grid(True)
plt.title('Heat Capacity C As a Function of Temperature')
plotNameString = 'plotDenseCvsTemperatureIsing' +  str(isingParticles) + 'DeltaT' + str(int(1000.*deltaTemperature)) + 'Milli.pdf'
fig.savefig(plotNameString)   # save the figure to file
print "   Maximum heat capacity C found to be", cMaximum, " at a temperature ", temperatureMaximumC
plt.close(fig)

#
# Write the output text file for <s> vs T
#
fileOutputName = 'fitAverageSpinDataLatticeSize' + str(isingArray) + 'DeltaT' + str(int(1000.*deltaTemperature)) + 'Milli.txt'
fileOutput = open(fileOutputName, 'w')
headerLineString = "  Average Spin for a Ferromagnet Near the Critical Temperature\n"
countTemperaturesLatticeLineString = '%3d %3d %s' % (temperatureSteps, isingArray, '\n')
fileOutput.write(headerLineString)
fileOutput.write(countTemperaturesLatticeLineString)
for line in range(temperatureSteps):
    dataLineString = '%0.3f %.4f %.4f %.4f %s' % (temperatureMean[line], spinMean1[line], spinMean2[line], spinMean3[line],'\n')
    fileOutput.write(dataLineString)

fileOutput.close()

#
# Write the output text file for <E> vs T
#
fileOutputName = 'fitEnergyDataLatticeSize' + str(isingArray) + 'DeltaT' + str(int(1000.*deltaTemperature)) + 'Milli.txt'
fileOutput = open(fileOutputName, 'w')
headerLineString = "  Energy for a Ferromagnet\n"
countTemperaturesLatticeLineString = '%3d %3d %s' % (temperatureSteps, isingArray, '\n')
fileOutput.write(headerLineString)
fileOutput.write(countTemperaturesLatticeLineString)

for temperatureStep in range(temperatureSteps):
    dataLineString = '%0.3f %.4f %.4e %.4f %s' % (temperatureMean[temperatureStep], energyMeanAllTemperatures[temperatureStep], energyStdAllTemperatures[temperatureStep], heatCapacityAllTemperatures[temperatureStep],'\n')
    fileOutput.write(dataLineString)

fileOutput.close()






plt.close(fig)

