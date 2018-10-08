# Program 5.2 Solution to LaPlace example in Figure 5.2 with the Gauss-Seidel algorithm modification, twoDGSLaPlace_Chapter5V1.py
#
#
#  Internal constraints and program operation
#  1) A rectangular grid of iRows and jColumns is used, i.e. a matrix of values
#  2) Only one copy of this grid (oldV) is needed
#  3) The initial distribution of values is printed first
#  4) The laplaceJacobiCalculate function is called.  This function causes the "new" matrix to
#     be updated from the "old" copy according to the first order iteration equation 5.14 (page 142)
#
#  Equation 5.14   V(i,j) = [V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1)]/4    Average around the point; boundary is not updated
#                  The left side is a new array while the right side is the old array replaced as the new array values are generated
#
#  5) The values in the new and the old matricies are compared point by point
#  6) The sum of absolute values of all the differences are returned by the updating function.
#  7) If that sum of differences is less than a convegence number, then the solution has been found.
#
#  Output from the program
#  1) Prints to the terminal of the initial guess and the converged solution
#  There is no plot output from this program
#

import  argparse                    # argument parser library

#
# Define the input parameter options and assign the default values and the variable types using the argument parser library
#
parser = argparse.ArgumentParser()

parser.add_argument('--iRows', default=15, type=int, help="Number of rows in the matrix of grid points; default 7")
parser.add_argument('--jColumns', default=15, type=int, help="Number of columns in the matrix of grid points; default 7")
parser.add_argument('--updateLimit', default=200, type=int, help="Maximum number of update calls allowed; default 100")
parser.add_argument('--leftBoundaryVoltage', default=+1.0, type=float, help="Voltage value at the left column boundary; default -1.0")
parser.add_argument('--rightBoundaryVoltage', default=+1.0, type=float, help="Voltage value at the right column boundary; default +1.0")
parser.add_argument('--interiorInitialVoltage', default=0.0, type=float, help="Initial value for all the interior grid points; default 0.0")
parser.add_argument('--epsilonConvergence', default=1.0e-04, type=float, help="Convergence limit between successive updates; default 1.0e-04")
parser.add_argument('--verbose', action='store_true', help="Give extra printout; default False")

args = parser.parse_args()
iRows = args.iRows
jColumns = args.jColumns
updateLimit = args.updateLimit
leftBoundaryVoltage = args.leftBoundaryVoltage
rightBoundaryVoltage = args.rightBoundaryVoltage
interiorInitialVoltage = args.interiorInitialVoltage
epsilonConvergence = args.epsilonConvergence
verbose = args.verbose

middleBoundaryVoltage = 0

def setBoundaryVoltage(array):       # Initializes the left and right column boundaries, and the top and lower insulating surfaces
    global iRows, jColumns, leftBoundaryVoltage, rightBoundaryVoltage
    jColumnsMinusOne = jColumns - 1
    iRowsMinusOne = iRows - 1
    #
    # Left and right side conductor initializations
    #
    print jColumnsMinusOne
    for iRow in range(iRows):
        array[iRow][0] = leftBoundaryVoltage
        array[iRow][jColumnsMinusOne] = rightBoundaryVoltage
        array[iRow][jColumnsMinusOne/2] = middleBoundaryVoltage

    #
    # Top and bottom row insulator initializations
    #
    stepV = (rightBoundaryVoltage - leftBoundaryVoltage)/float(jColumnsMinusOne)
    jColumn = 1
    while jColumn < jColumnsMinusOne:
        array[0][jColumn] = leftBoundaryVoltage + jColumn*stepV
        array[iRowsMinusOne][jColumn] = leftBoundaryVoltage + jColumn*stepV
        jColumn += 1
    return

def printVoltages(array):                          # Prints the voltage values at the grid points
    global iRows, jColumns
    for iRow in range(iRows):
        for jColumn in range(jColumns):
            voltageString = '%6.3f' % array[iRow][jColumn]
            print voltageString,
        print " "

    return

def laplaceJacobiCalculate(oldArray,  newArray):    # Does the LaPlace Jacobi method of solution with the Gauss-Seidel modification, i.e. oldArray and newArray are the same locations
    global iRows, jColumns, epsilonConvergence, verbose
    nCall = 0
    while nCall < updateLimit:
        diff1 = updateVoltages(oldArray, newArray)     # Gauss-Seidel modification uses the newest values as soon as they are available
        nCall += 1
        if(nCall==1):
            print "\n After first update"
            printVoltages(newArray)
        diff2 = updateVoltages(newArray, oldArray)         # Gauss-Seidel modification uses the newest values as soon as they are available
        if(abs(diff1) < epsilonConvergence or abs(diff2) < epsilonConvergence):
            return nCall
        if(verbose):
            print "\n At update iteration ", nCall, " the deltaV value is ", diff2
            printVoltages(newArray)

    print "\n  ***laplaceJacobiCalculate updating does not converge for ", nCall, " iterations***\n"
    return 0

def updateVoltages(oldArray, newArray):                   # Does the updating according to equation 5.14, newArray and oldArray are the same memory locations
    global iRows, jColumns

    #
    # Gauss-Jacobi update algorithm page 142
    # newV(i,j) = [oldV(i-1, j) + oldV(i+1, j), oldV(i, j-1) + oldV(i, j+1)]/4.
    # all (i,j) are updated except at the boundaries: jColumn=0, jColumn=jColumnsMinusOne, iRow = 0, and iRow = iRowsMinusOne
    #

    deltaV = 0.0
    jColumnsMinusOne = jColumns - 1
    iRowsMinusOne = iRows - 1
    jColumn = 1
    while jColumn < jColumnsMinusOne:
        iRow = 1
        if jColumn != 7:
            while iRow < iRowsMinusOne:
                oldValue = oldArray[iRow][jColumn]        # IMPORTANT CHANGE in Gauss-Seidel to save first the old value for deltaV comparison after the update
                newArray[iRow][jColumn] = 0.25*(oldArray[iRow+1][jColumn] + oldArray[iRow-1][jColumn] + oldArray[iRow][jColumn+1] + oldArray[iRow][jColumn-1])
                deltaV += abs(newArray[iRow][jColumn] - oldValue)
                iRow += 1
        jColumn += 1

    return deltaV

#
# Define the only voltage array in use
#
oldV = [[interiorInitialVoltage for jColumn in range(jColumns)] for iRow in range(iRows)]   # initialize all the old voltage grid points as the default interior voltage

#
# Set the boundary voltage for the voltage array (Gauss-Seidel has only one array in use)
#
setBoundaryVoltage(oldV)

print "\n Starting value of the voltages"
printVoltages(oldV)

successFail = laplaceJacobiCalculate(oldV, oldV)       # Gauss-Seidel algorithm uses the updated values as soon as they become available
if(successFail > 0):
    print "\n laplaceJacobiCalculate updating with Gauss-Seidel algorithm succeeded after ", successFail, " iterations"
else:
    print "\n laplaceJacobiCalculate updating did not converge"

printVoltages(oldV)

exit()