#  Program 10.3  oneDimensionalQM_Chapter10Final, to solve for Problem 2 on the fourth exam
#
# Solution to the one-dimensional Schroedinger equation with various potential functions 
# Uses matching method
#
# Producing Figures such as 10.13, pages 321
#
# 1) Program input
#    a) eigenvalue energy guess (=-1.5)
#    b) eDelta (=0.25) change in the energy
#    c) potential type (=1, Lennard-Jones potential, harmonic oscillator = 2, others listed below)
#    d) convergence check on the derivative change (= 0.01 fractional change)
#    e) psiParity (=+1) specification of the parity of the wave function
#    f) xDelta (=0.001) integration step size in one dimension
#
# 2) Program operation and features
#    a) Displays isplaying the rescaled energy level
#    b) A "psiParity" option
#    c) Different potential forms: square well (3), power law (4), and anharmonic potentials (5), cosine dependent potential (6)  choices
#    d) Uses Numerov iteration method by default
#
# 3) Program output
#    a) eigenvalue energy and last convergence step
#    b) oneDimensionalQMInitialFinal.pdf   initial guess and converged final solution of the wavefunction
#
#
import matplotlib
matplotlib.use('TkAgg')             # special code to make plots visible
import matplotlib.pyplot as plt     # get matplotlib plot functions
import  argparse                    # argument parser library
import numpy as np                  # numerical functions library used by python
import random                       # random number generator library
import scipy
from datetime import datetime
import matplotlib.animation as animation
import scipy.optimize as optimization # data fitting library

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()

parser.add_argument('--xLeft', default=0.0, type=float, help="The left-most position along the position axis; default 0.0")
parser.add_argument('--xRight', default=8.5, type=float, help="The right-most position along the position axis; default 8.5")
parser.add_argument('--lengthUnit',default=1.0, type=float, help="The length unit L for the double square well; default 1.0" )
parser.add_argument('--xMatch', default=1.4, type=float, help="The position at which the two solutions are matched; default 1.4")
parser.add_argument('--xWindow_Match', default=0.25, type=float, help="Window around which the two solutions are matched; default 0.25")
parser.add_argument('--xDelta', default=0.001, type=float, help="Step size in the x direction; default 0.001")
parser.add_argument('--vScale', default=10.0, type=float, help="Scaling factor for the potential plot; default 10.0")
parser.add_argument('--eGuess', default=-1.5, type=float, help="Initial guess for the energy eigenvalue; default -1.5")
parser.add_argument('--eDelta', default=0.5, type=float, help="Delta energy change for convergence; default 0.5")
parser.add_argument('--derivativeEpsilon', default=0.0001, type=float, help="Equal Left-Right derivatives criterion; default 0.0001")
parser.add_argument('--convergenceEpsilon', default=1.0e-14, type=float, help="Convergence criterion eDelta/eGuess; default 1.0e-14")
parser.add_argument('--iterationLimit', default=99, type=int, help="Limit for iterations over eGuess changes; default=99")
parser.add_argument('--nPotentialPoints', default=1000, type=int, help="Number of points over which to draw V(x); default=1000")
parser.add_argument('--psiParity', default=1, type=int, choices=[-1, 1], help="Parity of the wave function; default=1")
parser.add_argument('--vType', default=7, type=int, choices=[1, 2, 3, 4, 5, 6, 7], help="Choice of potential function; default=1 (L-J)")
parser.add_argument('--notNumerov', action='store_true', help="Option to NOT use Numerov algorithm; default False")
parser.add_argument('--plotInitialGuess', action='store_true', help="Plot initial guess solution; default False")
parser.add_argument('--verbose', action='store_true', help="Extra print out; default False")
parser.add_argument('--probabilityLimits', default=[2.5, 4.5], type=list, help="Limits within which to calculate probability of finding particle, default [2.5, 4.5]")

args = parser.parse_args()

vType = args.vType

xLeft = args.xLeft
if(vType == 1 and xLeft <= 0.0):
    print "\n  ***WARNING: For the vType = 1 Lennard-Jones potential the xLeft position must be at least 0.4  It is being reset to 0.4  ***\n"
    xLeft = 0.4

xMatch = args.xMatch
xRight = args.xRight
xRightXLeftMidpoint = (xRight + xLeft)/2.0
xRightMinusXLeftQuarter = (xRight - xLeft)/4.0
xWindow_Match = args.xWindow_Match
xDelta = args.xDelta
xDeltaSquared2 = 2.0*xDelta*xDelta
xDeltaSquaredDivide12 = xDelta*xDelta/12.0
xDeltaSquaredTimes5Divide12 = 5.0*xDeltaSquaredDivide12

#
# Variable for the double square well
#
lengthUnit = args.lengthUnit

nPoints = int((xRight - xLeft)/xDelta)  # goes out to the box boundary
jMatch =  int((xMatch - xLeft)/xDelta)  # index in the wave function array for the position of match point

eGuess = args.eGuess
eInitial = eGuess
recoveryAttempted = False
eDelta = args.eDelta

vScale = args.vScale
if(vType == 4 and vScale >= 8.0):
    print "\n   ***WARNING: For vType = 4 the vScale value", vScale, " is too large.  It is being reset to its maximum value of 8.0  ***\n";
    vScale = 8.0
if(vType == 6 and vScale < 1.75):
    print "\n   ***WARNING: For vType = 6 the vScale value", vScale, " is too small.  It is being reset to its minimum value of 1.75 ***\n";
    vScale = 1.75

psiParity = args.psiParity
derivativeEpsilon = args.derivativeEpsilon
convergenceEpsilon = args.convergenceEpsilon
iterationLimit = args.iterationLimit
nPotentialPoints = args.nPotentialPoints
notNumerov = args.notNumerov
plotInitialGuess = args.plotInitialGuess
verbose = args.verbose
probs = args.probabilityLimits
left_prob_limit, right_prob_limit = probs[0], probs[1]

def vPotential(xPosition):
    global xRightXLeftMidpoint, xRightMinusXLeftQuarter  # for the square well potenital
    
    invalidType = True
    if(vType==1):   # Lennard-Jones potential
        return 4.0*vScale*(pow(xPosition, -12.) - pow(xPosition, -6.))
 
    if(vType==2):  # harmonic oscillator potential
        return 0.5*vScale*xPosition*xPosition

    if(vType==3):   # square well, height = vScale and the width is half of the xRight
        if(abs(xPosition - xRightXLeftMidpoint) <= xRightMinusXLeftQuarter):
            return 0.0
        else:
            return vScale

    if(vType==4):   # power law potential
        return pow(xPosition,vScale)

    if(vType==5):   # anharmonic potential
        return 0.5*vScale*(pow(xPosition,2.0) + 0.1*pow(xPosition,4.0))

    if(vType==6):   # well with cosine dependence
        xDiff = xPosition - xRightXLeftMidpoint
        if(abs(xDiff) <= xRightMinusXLeftQuarter):
            return -0.5*vScale*(1.0 + np.cos(np.pi*xDiff/xRightMinusXLeftQuarter))
        else:
            return 0.0

    if (vType==7):  # modified square-well potential as per final exam
        if xLeft <= xPosition <= 2.5:
            return vScale
        elif 2.5 < xPosition <= 4.5:
            return 0
        elif 4.5 < xPosition <= 5.0:
            return vScale
        elif 5.0 < xPosition <= 6.0:
            return vScale / 2.0
        elif 6.0 < xPosition <= xRight:
            return vScale

    print "\n Invalid vType potential"   # this should not occur since the parser software restricts 1 <= vType <= 6
    exit()

def findTurningPoints(energyValue, nPoints):
    global xRightXLeftMidpoint, xRightMinusXLeftQuarter  # for the square well
    
    if(vType==3):   # special case for square well
        if(energyValue > vScale):
            #
            # Unbound
            #
            leftTurningPoint = xLeft
            rightTurningPoint = xRight
        else:
            leftTurningPoint = xRightXLeftMidpoint - xRightMinusXLeftQuarter
            rightTurningPoint = xRightXLeftMidpoint + xRightMinusXLeftQuarter
        return leftTurningPoint, rightTurningPoint

    equalValue = 0.03*abs(energyValue)
    xPosition = xLeft
    leftNotFound = True
    for j in range(10*nPoints):
        if(abs(vPotential(xPosition) - energyValue) < equalValue or energyValue > vPotential(xPosition)):
            leftTurningPoint = xPosition
            leftNotFound = False
            break

        xPosition += 0.1*xDelta
    
    xPosition = xRight
    rightNotFound = True
    for j in range(10*nPoints):
        if(abs(vPotential(xPosition) - energyValue) < equalValue or energyValue > vPotential(xPosition)):
            rightTurningPoint = xPosition
            rightNotFound = False
            break

        xPosition -= 0.1*xDelta

    if(leftNotFound):
        print "\n [findTurningPoints]  Left turning point not found for energy value", energyValue
        exit()

    if(rightNotFound):
        print "\n [findTurningPoints]  Right turning point not found for energy value", energyValue
        exit()
            
    return leftTurningPoint, rightTurningPoint

vRescale = 1.0

print "\n Solving the Schrodinger equation for a given potential energy function ( vType =", vType, ")"

vString ='Undefined Potential'
vRescale = 2.0/vScale

if(vType == 1):
    vRescale = -1.0/vPotential(1.12)
    vString = 'Lennard-Jones Potential'
    print "   Using the Lennard-Jones potential, with vRescale", vRescale

if(vType == 2):
    vRescale = 2.0/vPotential(0.5*xRight)
    vString = 'Harmonic Oscillator Potential, force constant K = %.1f' % vScale
    print "   Using the Harmonic Oscillator Potential with a force constant", vScale

if(vType == 3):
    vRescale = 2.0/vScale
    vString = 'Square Well V0 = %.1f' % vScale
    print "   Using the Square Well potential with a barrier height", vScale, " a width from", xRightXLeftMidpoint-xRightMinusXLeftQuarter, " to ", xRightXLeftMidpoint+xRightMinusXLeftQuarter

if(vType == 4):
    vRescale = 2.0/vPotential(0.5*xRight);
    vString = 'Power Law Potential exponent = %.1f' % vScale
    print "   Using the Power Law potential with an exponent", vScale

if(vType == 5):
    vRescale = 2.0/vPotential(0.5*xRight);
    vString = 'Anharmomic Oscillator Potential k/2m %.1f' % vScale
    print "   Using an Anharominic Potential with a strength", vScale

if(vType == 6):
    vRescale = 1.0/vScale
    vString = 'Cosine Potential k %.1f' % vScale
    print "   Using an Cosine Potential with a strength", vScale

if(vType == 7):
    vRescale = 1.0 / vScale
    vString = 'Modified square-well potential %.1f' % vScale
    print "   Using a modified square-well potential with potential strength", vScale


print "   Left boundary position", xLeft, "  right boundary position", xRight
print "   The Left-Right matching position is", xMatch, " with a position window for matching ", xWindow_Match
print "   Number of spatial points to boundary =", nPoints + 1, " with xDelta", xDelta
print "   The Left-Right derivatives matching criterion is", derivativeEpsilon
print "   The energy guesses iteration limit is", iterationLimit
print "   The initial energy guess is", eInitial, " with an eDelta", eDelta, " and an eDelta/eGuess limit", convergenceEpsilon
print "   The wavefunction parity value is", psiParity
if(notNumerov):
    titleString = 'Solving the Schrodinger Equation with the Finite Differences Algorithm'
    print "   Using the textbook, second order finite difference iteration algorithm to solve the Schrodinger equation"
else:
    titleString = 'Solving the Schrodinger Equation with the Numerov Algorithm'
    print "   Using the Numerov iteration algorithm to solve the Schrodinger equation with fifth order accuracy"

psiFunctionLeft = [-2.0 for j in range(nPoints)]    # array to contain the left side solution, with dummy initialization values
psiFunctionRight = [-2.0 for j in range(nPoints)]   # array to contain the right side solution, with dummy initialization values

#
# real initialization for two points at either end
#
psiFunctionLeft[0] = 0.0
psiFunctionLeft[1] = 0.0001*xDelta*psiParity
psiFunctionRight[nPoints-1] = 0.0
psiFunctionRight[nPoints-2] = 0.0001*xDelta

notOverLimit = True
notConverged = True
kCount = 0
lastDifference = 0.0

nPointsMinusTwo = nPoints - 2
#
# loop over the energy guesses
#
while notConverged and notOverLimit:
    xPositionLeft = xLeft + xDelta
    xPositionRight = xRight - xDelta

    leftMax = 0.0
    rightMax = 0.0

    leftPoints = 0
    rightPoints = 0

    #
    # iteration loop over position using a particular energy guess
    #
    for jPoint in range(nPointsMinusTwo):
        j = jPoint + 2

        #
        # check if x value is still below match point for left solution
        #
        if(xPositionLeft - xMatch <= xWindow_Match):
            leftPoints += 1
            vPotentialLeft = vPotential(xPositionLeft);

            if(notNumerov):
                psiFunctionLeft[j] = 2.0*psiFunctionLeft[j-1] - psiFunctionLeft[j-2] - xDeltaSquared2*(eGuess - vPotentialLeft)*psiFunctionLeft[j-1]
            else:
                kSquaredCurrent = 2.0*(eGuess - vPotentialLeft)                        # assumes hbar = 1, and mass m = 1
                kSquaredPrevious = 2.0*(eGuess - vPotential(xPositionLeft - xDelta))   # assumes hbar = 1, and mass m = 1
                kSquaredNew = 2.0*(eGuess - vPotential(xPositionLeft + xDelta))        # assumes hbar = 1, and mass m = 1
                fact1 = 2.0*(1.0 - xDeltaSquaredTimes5Divide12*kSquaredCurrent)
                fact2 = 1.0 + xDeltaSquaredDivide12*kSquaredPrevious
                fact3 = 1.0 + xDeltaSquaredDivide12*kSquaredNew
                if(fact3 !=0.0):
                    psiFunctionLeft[j] = (fact1*psiFunctionLeft[j-1] - fact2*psiFunctionLeft[j-2])/fact3
                else:
                    print "\n Attempted divide by 0 in the Numerov iteration algorithm, for left solution"

            if(abs(psiFunctionLeft[j]) > leftMax):
                leftMax = abs(psiFunctionLeft[j])

            xPositionLeft += xDelta; # x-value for continuing iteration loop for left side solution
    
        #
        # check if x value is still above match point for right solution
        #
        if(xMatch - xPositionRight <= xWindow_Match):
            rightPoints += 1
            vPotentialRight = vPotential(xPositionRight);

            if(notNumerov):
                psiFunctionRight[nPoints-j-1] = 2.0*psiFunctionRight[nPoints-j] - psiFunctionRight[nPoints-j+1] - xDeltaSquared2*(eGuess - vPotentialRight)*psiFunctionRight[nPoints-j]
            else:
                kSquaredCurrent = 2.0*(eGuess - vPotentialRight); # assumes hbar = 1, and mass m = 1
                kSquaredPrevious = 2.0*(eGuess - vPotential(xPositionRight + xDelta)) # assumes hbar = 1, and mass m = 1
                kSquaredNew = 2.0*(eGuess - vPotential(xPositionRight - xDelta)) # assumes hbar = 1, and mass m = 1
                fact1 = 2.0*(1.0 - xDeltaSquaredTimes5Divide12*kSquaredCurrent)
                fact2 = 1.0 + xDeltaSquaredDivide12*kSquaredPrevious
                fact3 = 1.0 + xDeltaSquaredDivide12*kSquaredNew
                if(fact3 !=0.0):
                    psiFunctionRight[nPoints-j-1] = (fact1*psiFunctionRight[nPoints-j] - fact2*psiFunctionRight[nPoints-j+1])/fact3
                else:
                    print "\n Attempted divide by 0 in the Numerov iteration algorithm, for right solution"
                    exit(1) # divide by 0 safety check

            if(abs(psiFunctionRight[nPoints-j-1]) > rightMax):
                rightMax = abs(psiFunctionRight[nPoints-j-1])

            xPositionRight -= xDelta; # x-value for continuing iteration loop for right side solution

        if(xPositionLeft - xMatch > xWindow_Match and xMatch - xPositionRight > xWindow_Match):
            if(leftMax == 0.0 or rightMax == 0.0):
                print "\n  Abnormal result has leftMax = ", leftMax, " and rightMax ", rightMax
                print "  This will result in a division by zero."
                print "  Program is exiting."
                exit()

            #
            # normalize the wave functions to have their maxima at 1
            #
            for i in range(leftPoints):
                psiFunctionLeft[i] /= leftMax
    
            for i in range(rightPoints):
                psiFunctionRight[nPoints-i-1] /= rightMax

            if(lastDifference == 0.0):  # save solutions from initial iteration, before matching
                psiFunctionLeftInitial = [psiFunctionLeft[i] for i in range(leftPoints)]
                psiFunctionRightInitial = [psiFunctionRight[i] for i in range(rightPoints)]
                leftPointsInitial = leftPoints
                rightPointsInitial = rightPoints

            break # break out of iteration over the position points

    psiLeftMatch = psiFunctionLeft[jMatch]
    psiRightMatch = psiFunctionRight[jMatch]
    matchRatio = psiLeftMatch/psiRightMatch
    #
    # Normalize the right function to the left function at the matching point
    #
    xPosition = xRight
    for i in range(rightPoints):
        psiFunctionRight[nPoints-i-1] *= matchRatio
        xPosition -= xDelta

    #
    # currentDifference is the current difference between the left and the right derivatives
    #
    currentDifference = (psiFunctionLeft[jMatch+1] - psiFunctionLeft[jMatch]) - (psiFunctionRight[jMatch+1] - psiFunctionRight[jMatch])
    currentDifference /= xDelta
    meanDerivative = (psiFunctionLeft[jMatch+1] - psiFunctionLeft[jMatch]) + (psiFunctionRight[jMatch+1] - psiFunctionRight[jMatch])
    meanDerivative /= (2.0*xDelta)

    if(abs(currentDifference) < derivativeEpsilon*abs(meanDerivative)):
        notConverged = False
        iStatus = 1
        print "\n  Derivative ratio check value achieved after", kCount, "energy guess iterations"
        if(verbose):
            print "   Left solution ", psiFunctionLeft[jMatch], "  Right solution ", psiFunctionRight[jMatch], "  matchRatio", matchRatio

    if(lastDifference==0):
        #
        # special guess of the first energy guess
        #
        if(currentDifference > 0.0):
            eGuess += eDelta
        else:
            eGuess -= eDelta
    else:
        #
        # check if the left - right derivative difference has changed signs
        # if so, then halve the energy change and change its sign
        #
        if(currentDifference/lastDifference < 0.0):
            #
            # halve the delta value and change in sign
            #
            eDelta *= 0.5
        
        if(currentDifference > 0.0):
            eGuess += eDelta
        else:
            eGuess -= eDelta
                
    lastDifference = currentDifference

    if(kCount > iterationLimit):
        notOverLimit = False
        iStatus = 2
        print "\n Reached kCount iteration limit for energy guesses ", iterationLimit

    if(abs(eGuess) < 1.0e-4):
        print "\n  Algorithm has attempted an Energy = 0 guess.  Initial energy guess had the wrong sign."
        if(recoveryAttempted):
            print "\n  Recovery procedure has failed"
            exit()
        else:
            print "  Program will attempt to recover"
            eGuess == -eInitial
            recoveryAttempted = True
    else:
        if(abs(eDelta/eGuess) < convergenceEpsilon):
            notOverLimit = False
            iStatus = 3
            print "\n  Reached eDelta/eGuess limit: ", convergenceEpsilon
            print "   Current derivative difference ", currentDifference, " and mean derivative",  meanDerivative, "  compared to derivative convergence criterion", derivativeEpsilon*abs(meanDerivative)
            print "   It is possible that the xMatch position is poorly chosen, to be where the wavefunctions and/or the derivatives are very small"
            print "   ---> Try changing the xMatch position to be at the mid-point of the xLeft and the xRight values = ", (xLeft + xRight)/2.0, " in this case."
            print " "

    kCount += 1

print "  eGuess = ", eGuess,", eDelta = ", eDelta

if(iStatus != 1):
    print "\n  Program is exiting without making any plots"
    exit()

# code to set up two plots in a single figure
fig = plt.figure(1)          # start a figure

normFactor = 0.0
xPosition = xLeft
for j in range(leftPoints):
    if(xPosition > xMatch):
        break
    normFactor += psiFunctionLeft[j]*psiFunctionLeft[j]
    xPosition += xDelta

xPosition = xRight
for j in range(rightPoints):
    normFactor += psiFunctionRight[nPoints - j - 1]*psiFunctionRight[nPoints - j - 1]
    if(xPosition <= xMatch):
        break
    xPosition -= xDelta

normFactor = 1.0/(normFactor*xDelta)

probabilityIntegral, prob_final_exam = 0.0, 0.0
xPosition = xLeft
for j in range(leftPoints):
    if(xPosition > xMatch):
        break
    probabilityIntegral += normFactor*psiFunctionLeft[j]*psiFunctionLeft[j]
    if xPosition > left_prob_limit:
        prob_final_exam += normFactor*psiFunctionLeft[j]*psiFunctionLeft[j]
    xPosition += xDelta

xPosition = xRight
for j in range(rightPoints):
    probabilityIntegral += normFactor*psiFunctionRight[nPoints - j - 1]*psiFunctionRight[nPoints - j - 1]
    if xPosition < right_prob_limit:
        prob_final_exam += normFactor*psiFunctionRight[nPoints - j - 1]*psiFunctionRight[nPoints - j - 1]
    if(xPosition <= xMatch):
        break
    xPosition -= xDelta
probabilityIntegral *= xDelta
prob_final_exam *= xDelta * 100

print "  Probability of finding particle in between", left_prob_limit, "and", right_prob_limit, "is", str(round(prob_final_exam, 2)) + "%"

normFactorString = '  Wavefunction normalization factor %6.3e produces an integrated probability value %6.3e' % (normFactor, probabilityIntegral)
print normFactorString

yMax = 3.5
yMaxV = 2.0
yMin = -1.0
xMax = xRight
xMin = xLeft

yMinV = 1.0e5
deltaAxis = (xMax - xMin)/nPotentialPoints
xPosition = xMin
xPotential = []
yPotential = []
for j in range(nPotentialPoints):
    yValue = vRescale*vPotential(xPosition)
    if(yValue <= yMaxV):
        xPotential.append(xPosition)
        yPotential.append(yValue)

    xPosition += deltaAxis
    if(yValue < yMinV):
        yMinV = yValue

plt.plot(xPotential, yPotential, 'k', label=vString)

if(plotInitialGuess):
    xLeftPlotInitial = []
    yLeftPlotInitial = []
    leftPositionInitial = xLeft
    for j in range(leftPointsInitial):
        xLeftPlotInitial.append(leftPositionInitial)
        yLeftPlotInitial.append(psiFunctionLeftInitial[j])
        leftPositionInitial += xDelta

    xRightPlotInitial = []
    yRightPlotInitial = []
    rightPositionInitial = xRight
    for j in range(rightPointsInitial):
        xRightPlotInitial.append(rightPositionInitial)
        yRightPlotInitial.append(psiFunctionRightInitial[rightPointsInitial - j - 1])
        rightPositionInitial -= xDelta

    plt.plot(xLeftPlotInitial, yLeftPlotInitial, 'r', label='Not converged solution, wrong initial energy')
    plt.plot(xRightPlotInitial, yRightPlotInitial, 'r')

xLeftPlot = []
yLeftPlot = []
leftPosition = xLeft
for j in range(leftPoints):
    xLeftPlot.append(leftPosition)
    yLeftPlot.append(psiFunctionLeft[j])
    leftPosition += xDelta

xRightPlot = []
yRightPlot = []
rightPosition = xRight
for j in range(rightPoints):
    xRightPlot.append(rightPosition)
    yRightPlot.append(psiFunctionRight[nPoints - j - 1])
    rightPosition -= xDelta

plt.plot(xLeftPlot, yLeftPlot, 'b', label='Converged left side wavefunction')
plt.plot(xRightPlot, yRightPlot, 'g', label='Converged right side wavefunction')

xMatchLine = [xMatch, xMatch]
yMatchLine = [0.0, psiFunctionLeft[jMatch]]
matchString = 'Left-Right match position, x = %.1f' % xMatch
plt.plot(xMatchLine, yMatchLine, 'k:', label=matchString)

leftTurningPoint, rightTurningPoint = findTurningPoints(eGuess, nPoints)
xEnergy = [leftTurningPoint, rightTurningPoint]
yEnergy = [vRescale*eGuess, vRescale*eGuess]
energyString = 'Converged energy eigenvalue = %.3e' % eGuess
plt.plot(xEnergy, yEnergy, 'm', label=energyString)

yMinLeft = min(yLeftPlot)
yMinRight = min(yRightPlot)

yMinAll = min(yMinV, yMinLeft, yMinRight)
yMin = yMinAll - 0.10*(yMax - yMinAll)
plt.xlim(xMin, xMax)
plt.ylim(yMin, yMax)

if(vType == 1):
    plt.xlabel('Radial separation of Argon atoms (sigma units)')
else:
    plt.xlabel('x-position')                             # add axis labels
yLabelString = r'$\psi$' + ' wavefunction amplitude or rescaled potential V(x)'
plt.ylabel(yLabelString)
plt.plot(xPotential, yPotential, 'k')
plt.grid(True)

plt.title(titleString)

vFactor = 1.0/vRescale
vRescaleString = 'Energy rescale factor %.2e' % vFactor
xPositionText = xMin + 0.50*(xMax - xMin)
yPositionText = vRescale*eGuess - 0.03*(yMax - yMin)
plt.text(xPositionText, yPositionText, vRescaleString)

plt.legend(loc=1)
plt.show()
