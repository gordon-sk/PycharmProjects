##   Program 8.2  Numerical Solution to the Second Order Ising Model for Ferromagnetism
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
import scipy
import scipy.optimize
import random                       # random number generator library
from datetime import datetime


# I need: <s> vs T as in figure 8.7
# I need: <E> vs T as in figure 8.8

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()

parser.add_argument('--isingArray', default=10, type=int, help="Dimension of Ising lattice; default 10")
parser.add_argument('--startTemperature', default=1.00, type=float, help="Starting temperature, units of J/kB; default 1.00")
parser.add_argument('--finalTemperature', default=5.00, type=float, help="Final temperature units of J/kB; default 5.00")
parser.add_argument('--deltaTemperature', default=0.25, type=float, help="Temperature increment; default 0.25")

parser.add_argument('--pairEnergy', default=1.0, type=float, help="Pair energy; default 1.00")
parser.add_argument('--seedFixed', action='store_true', help="Use a fixed random number seed; default False")
parser.add_argument('--timeSteps', default=3000, type=int, help="Number of time steps; default 3000")
parser.add_argument('--allSpinsRandom', action='store_true', help="Initially all spins point randomly up or down; default False")
parser.add_argument('--printTimes', action='store_true', help="Print elapsed wall clock time at each temperature; default False")
parser.add_argument('--verbose', action='store_true', help="Extra print out; default False")
parser.add_argument('--debug', action='store_true', help="Debugging print out; default False")
parser.add_argument('--skipNegativeCheck', action='store_true', help="Skip check against complete flip to negative; default False")
parser.add_argument('--triggerNegative', default=-0.75, type=float, help="Check against complete flip to negative; default -0.75")
args = parser.parse_args()

isingArray = args.isingArray
isingParticles = isingArray*isingArray*isingArray
isingParticlesInverse = 1.0/float(isingParticles)
isingParticlesMinusOneInverseSquareRoot = 1.0/np.sqrt(float(isingParticles) - 1.0)

arraySpins = [[[0 for kStack in range(isingArray)] for jColumn in range(isingArray)] for iRow in range(isingArray)]  # spin value for the lattice element [iRow][jColumn]

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

triggerNegative = args.triggerNegative
skipNegativeCheck = args.skipNegativeCheck
negativeCheck = False
if(skipNegativeCheck == False):
    negativeCheck = True


debug = args.debug

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

#
# The spin array values at each time step, re-calculated at each temperature step: spinArray[timeStep][iRow][jColumn]
#
spinsAllSteps = [[[[0 for kStack in range(isingArray)] for jColumn in range(isingArray)] for iRow in range(isingArray)] for timeStep in range(timeSteps)]

#
# Define the two dimensional arrays for spin average at each step sValue[timeStep][tempStep]
# Do like wise for timeValue[timeStep][tempStep], energyValue[timeStep][tempStep], varianceValue[timeStep][tempStep], correlationValue]timeStep][tempStep]
#
sValue = [[0.0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]
timeValue = [[0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]
energyValue = [[0.0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]
varianceValue = [[0.0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]
correlationValue = [[0.0 for tempStep in range(temperatureSteps)] for timeStep in range(timeSteps)]

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
adjacentLeftStack = [0 for kStack in range(isingArray)]
adjacentRightStack = [0 for kStack in range(isingArray)]
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
    adjacentLowerRow[iRow] = iRowLower # ***Original version had mistake with  adjacentUpperRow[iRow] = iRowLower   This set both the arrays with wrong values**.

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

for kStack in range(isingArray):
    #
    # left stack neighbor
    #
    if(kStack == 0):
        kStackLeft = isingArrayMinusOne
    else:
        kStackLeft = kStack - 1
    adjacentLeftStack[kStack] = kStackLeft

    #
    # right stack neighbor
    #
    if kStack == isingArrayMinusOne:
        kStackRight = 0
    else:
        kStackRight = kStack + 1
    adjacentRightStack[kStack] = kStackRight

def initializeSpins(isingArray, arraySpins):
    averageSpin = 0
    for kStack in range(isingArray):
        for iRow in range(isingArray):
            for jColumn in range(isingArray):
                thisSpin = 1
                if(allSpinsRandom):
                    thisSpin = random.randrange(-1,3,2)    # returns either -1 or +1 randomly, note that the 3 does not get used

                arraySpins[kStack][iRow][jColumn] = thisSpin    # returns either -1 or +1 randomly, note that the 3 does not get used
                averageSpin += thisSpin

    return averageSpin*isingParticlesInverse

initialSpinAverage = initializeSpins(isingArray, arraySpins)
print "  Initial average spin ", initialSpinAverage

def printSpinArray(isingArray, arraySpins):
    for iRow in range(isingArray):
        print "  ",      # leading two spaces for each row in the print
        for jColumn in range(isingArray):
            for kStack in range(isingArray):
                thisSpinString = '%2d' % (arraySpins[iRow][jColumn][kStack])
                print thisSpinString

        print " "  # carriage return/end of line for a given row

    return

# printSpinArray(isingArray, arraySpins)

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

            jColumnLeft = adjacentLeftColumn[jColumn]
            jColumnRight = adjacentRightColumn[jColumn]

            for kStack in range(isingArray):

                thisSpin = arraySpins[iRow][jColumn][kStack]
                kStackLeft = adjacentLeftStack[kStack]
                kStackRight = adjacentRightStack[kStack]

                thisProduct = thisSpin*(arraySpins[iRowUpper][jColumn][kStack] + arraySpins[iRowLower][jColumn][kStack]
                                        + arraySpins[iRow][jColumnLeft][kStack] + arraySpins[iRow][jColumnRight][kStack]
                                        + arraySpins[iRow][jColumn][kStackLeft] + arraySpins[iRow][jColumn][kStackRight])

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
print "  Initial energy ", initialEnergyAverage, " averaged per spin pair, with a variance", initialEnergyAverage*ratioStdOverMean

def updateSpin(rowPosition, columnPosition, stackPosition, isingArray, arraySpins, temperature):
    #
    # decide whether to flip the spin location (rowPosition, columnPosition)
    #
    
    eFlip = 0.0;
    thisSpin = arraySpins[rowPosition][columnPosition][stackPosition]
    
    #
    # upper row neighbor
    #
    iRowUpper = adjacentUpperRow[rowPosition]
    if(iRowUpper == rowPosition):
        print "\n  Same upper row error"
        exit()
    if(thisSpin*arraySpins[iRowUpper][columnPosition][stackPosition] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy
    
    #
    # lower row neighbor
    #
    iRowLower = adjacentLowerRow[rowPosition]
    if(iRowLower == rowPosition):
        print "\n  Same lower row error: lower ", iRowLower, rowPosition
        exit()
    if(thisSpin*arraySpins[iRowLower][columnPosition][stackPosition] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy
    
    #
    # left column neighbor
    #
    iColumnLeft = adjacentLeftColumn[columnPosition]
    if(iColumnLeft == columnPosition):
        print "\n Same left column error"
        exit()
    if(thisSpin*arraySpins[rowPosition][iColumnLeft][stackPosition] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy
    
    #
    # right column neighbor
    #
    iColumnRight = adjacentRightColumn[columnPosition]
    if(iColumnRight == columnPosition):
        print "\n Same right column error"
        exit()
    if(thisSpin*arraySpins[rowPosition][iColumnRight][stackPosition] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy

    #
    # left stack neighbor
    #
    iStackLeft = adjacentLeftStack[stackPosition]
    if(iStackLeft == stackPosition):
        print "\n Same right column error"
        exit()
    if(thisSpin*arraySpins[rowPosition][iColumnRight][stackPosition] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy

    #
    # right stack neighbor
    #
    iStackRight = adjacentRightStack[stackPosition]
    if(iStackRight == stackPosition):
        print "\n Same right column error"
        exit()
    if(thisSpin*arraySpins[rowPosition][iColumnRight][stackPosition] == -1):
        eFlip -= flipEnergy
    else:
        eFlip += flipEnergy

    #
    # Based on the energy calculated, determine if this spin should be flipped to give a lower energy state
    #

    if(eFlip <= 0.0):
        arraySpins[rowPosition][columnPosition][stackPosition] = -thisSpin    # always flip if energy will be lower after flip

    if(eFlip > 0.0):
        rndm = random.uniform(0, 1.0)  #  random number [0,1)
        expTest = np.exp(-eFlip/temperature)
        if(rndm <= expTest):
            arraySpins[rowPosition][columnPosition][stackPosition] = -thisSpin;  # flip the spin randomly according to temperature

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
            for kStack in range(isingArray):
                #
                # check if this spin should be flipped
                #
                updateSpin(iRow, jColumn, kStack, isingArray, arraySpins, temperature)

    solution = 0  # sum of all the spin values
    for iRow in range(isingArray):
        for jColumn in range(isingArray):
            for kStack in range(isingArray):
                solution += arraySpins[iRow][jColumn][kStack]

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
print "\n  Begin cycling over temperatures at ", startDateTimeString
sSolutionSum = 0.0
spinMean1 = []
spinMean2 = []
spinMean3 = []
energyMean1 = []
energyMean2 = []
energyMean3 = []
temperatureMean = []
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
    
        energy, energyVariance = calculateEnergy(isingArray, arraySpins)
        energyValue[timeStep][temperatureStep] = energy # store energy value
        varianceValue[timeStep][temperatureStep] = energyVariance

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
        resultsString = '%s %.4f %s %.1f %s %.4f %s %.3e %s %.3f' % ('  At T =', temperature, ', total elapsed time', totalSeconds, 'seconds;  <s> =', spinAverageThisTemperatureMean,
                                                                     '+/- ', spinAverageThisTemperatureStdMean, ', E =',  energyThisTemperatureMean)
        print resultsString

    temperature += deltaTemperature

endDateTimeString = str(datetime.now())
print "\n  End cycling over temperatures at ", endDateTimeString

#
# Produce four sets of <s> vs time plots, two per figure
#
print "\n  Producing four sets of <s> vs time plots"
stPlotIndex = -1
for stPlot in range(2):
    fig = plt.figure(1)       # start a figure
    fig.subplots_adjust(hspace=0.6)

    plt.subplot(211)    # this sets the upper half plot for the <s> vs time at a first given temperature
    
    stPlotIndex += 1
    temperatureIndex = indexPlotTemperature[stPlotIndex]
    timePlot = [timeValue[timeStep][temperatureIndex] for timeStep in range(firstThirdTime)]
    sPlot1 = []
    sPlot2 = []
    sPlot3 = []

    for timeStep in range(firstThirdTime):
        sPlot1.append(sValue[timeStep][temperatureIndex])
        sPlot2.append(sValue[timeStep+firstThirdTime][temperatureIndex])
        sPlot3.append(sValue[timeStep+secondThirdTime][temperatureIndex])

    plt.scatter(timePlot, sPlot1, s=1, c='b', marker='.')
    plt.scatter(timePlot, sPlot2, s=1, c='r', marker='.')
    plt.scatter(timePlot, sPlot3, s=1, c='g', marker='.')

    plt.xlim(0,firstThirdTime)
    plt.ylim(-1.2, 1.2)
    plt.grid(True)
    titleString = '%s %.2f' % ('Second Order Ising Spin As a Function of Time for T = ', temperatureMean[temperatureIndex])
    plt.title(titleString)
    plt.ylabel('<s> At Each Time')
    plt.xlabel('Time Step')
    print "   ", stPlotIndex, ")  ", titleString

    plt.subplot(212)    # this sets the lower half plot for the <s> vs time at a second given temperature

    stPlotIndex += 1
    temperatureIndex = indexPlotTemperature[stPlotIndex]
    timePlot = [timeValue[timeStep][temperatureIndex] for timeStep in range(firstThirdTime)]
    sPlot1 = []
    sPlot2 = []
    sPlot3 = []
    
    for timeStep in range(firstThirdTime):
        sPlot1.append(sValue[timeStep][temperatureIndex])
        sPlot2.append(sValue[timeStep+firstThirdTime][temperatureIndex])
        sPlot3.append(sValue[timeStep+secondThirdTime][temperatureIndex])

    plt.scatter(timePlot, sPlot1, s=1, c='b', marker='.')
    plt.scatter(timePlot, sPlot2, s=1, c='r', marker='.')
    plt.scatter(timePlot, sPlot3, s=1, c='g', marker='.')

    plt.xlim(0,firstThirdTime)
    plt.ylim(-1.2, 1.2)
    plt.grid(True)
    titleString = '%s %.2f' % ('Second Order Ising Spin As a Function of Time for T = ', temperatureMean[temperatureIndex])
    plt.title(titleString)
    plt.ylabel('<s> At Each Time')
    plt.xlabel('Time Step')
    print "   ", stPlotIndex, ")  ", titleString
      
    if(stPlot == 0):
        plotNameString = 'firstPlotSvsTimeIsing' +  str(isingParticles) + '.pdf'
    else:
        plotNameString = 'secondPlotSvsTimeIsing' +  str(isingParticles) + '.pdf'
    fig.savefig(plotNameString)   # save the figure to file
    plt.close(fig)

print "\n  Producing <s> vs temperature plot"

fig = plt.figure(1)       # start a figure
plt.scatter(temperatureMean, spinMean1, s=6, c='b', marker='o')
plt.scatter(temperatureMean, spinMean2, s=6, c='r', marker='o')
plt.scatter(temperatureMean, spinMean3, s=6, c='g', marker='o')
xMaximum = finalTemperature + 0.25
xMinimum = 0.75
plt.xlim(xMinimum, xMaximum)
plt.ylim(-0.25, 1.25)
plt.ylabel('<s> At Each Temperature')
plt.xlabel('Temperature (units of  J/kB)')
xPositionText = xMinimum + 0.52*(xMaximum - xMinimum)
yPositionText = 0.80
latticeString = '%s %d %s %d %s %d' % ('Lattice size ', isingArray, ' x ', isingArray, 'x', isingArray)
plt.text(xPositionText, yPositionText, latticeString)
plt.grid(True)
plt.title('Second Order Ising Spin As a Function of Temperature')
plotNameString = 'plotSvsTemperatureIsing' +  str(isingParticles) + '.pdf'
fig.savefig(plotNameString)   # save the figure to file
print "   Second Order Ising Spin As a Function of Temperature"
plt.close(fig)

print "\n  Producing Energy vs temperature plot"

fig = plt.figure(1)       # start a figure
plt.scatter(temperatureMean, energyMean1, s=6, c='b', marker='o')
plt.scatter(temperatureMean, energyMean2, s=6, c='r', marker='o')
plt.scatter(temperatureMean, energyMean3, s=6, c='g', marker='o')
xMaximum = finalTemperature + 0.25
xMinimum = 0.75
#plt.xlim(xMinimum, xMaximum)
#plt.ylim(-2.2, 0.0)
plt.ylabel('Energy per particle At Each Temperature')
plt.xlabel('Temperature (units of  J/kB)')
xPositionText = xMinimum + 0.21*(xMaximum - xMinimum)
yPositionText = -0.40
latticeString = '%s %d %s %d %s %d' % ('Lattice size ', isingArray, ' x ', isingArray, 'x', isingArray)
plt.text(xPositionText, yPositionText, latticeString)
plt.grid(True)
plt.title('Energy Per Spin Particle As a Function of Temperature')
plotNameString = 'plotEvsTemperatureIsing' +  str(isingParticles) + '.pdf'
fig.savefig(plotNameString)   # save the figure to file
print "   Energy Per Spin Particle As a Function of Temperature"
plt.close(fig)


