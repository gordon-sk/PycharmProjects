#  plotBaseballB2OverM.py, version January 26, 2017
#

import matplotlib                   # library of plotting functions
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt     # assign the main plotting function in the library to the name plt
import math as mp
import numpy as np                  # library of widely used numerical functions in python
import  argparse                    # argument parser library, to handle command line parameters

parser = argparse.ArgumentParser()   # obtain the object called parser which can interpret command line parameters
parser.add_argument('--minSpeed', type=float, default=0.0, help='Minimum speed in mph for plotting, default = 0.0')   # declare a parameter and its default value
parser.add_argument('--maxSpeed', type=float, default=200.0, help='Maximum speed in mph for plotting, default = 200.0') # declare a parameter and its default value
parser.add_argument('--deltaSpeed',type=float, default=1.0, help='Delta speed in mph for plot spacing, default = 1.0') # declare a parameter and its default value
args = parser.parse_args()  # retrieve the parameters into another object called arg

minSpeedMPH = args.minSpeed      # retrieve the minimum speed
maxSpeedMPH = args.maxSpeed      # retrieve the maximum speed
deltaSpeedMPH = args.deltaSpeed      # retrieve the increment in speed

#
# Indicate what were the two parameters
#
print "\n Plotting from minimum speed ", minSpeedMPH, " mph to maximum speed ", maxSpeedMPH, " mph in steps of ", deltaSpeedMPH, "mph"

xStart = minSpeedMPH

#
# Declare two empty lists to hold the (x,y) data points
#
xplot = []
yplot = []

speedConversionFactor = 0.4407

xValue = xStart
while xValue <= maxSpeedMPH:                            # note the : mark at the end
    speedMKS = xValue*speedConversionFactor
    yValue = 0.0039 + (0.0058/(1.0+mp.exp((speedMKS - 35.0)/5.0)))
    xplot.append(xValue)                                # an indentation must be given, here it is 4 spaces  yValue = np.sqrt(radiusSquared - xValue*xValue)     # all other statements in the loop must
    yplot.append(yValue)

    xValue = xValue + deltaSpeedMPH                     # the xValue is increased by one step, this is the last statement in the loop

#
# Plot this set of (x,y) coordinates as a continuous blue line
#
plt.plot(xplot,yplot, 'b', label='Normal baseball')      # put a label identification for this plot

#
# Add labels for the x and y coordinates, use a rectangular grid, use expanded ranges for the (x,y) axes
# Add the legend identifying the plots and put in a title for the figure
#
plt.xlabel('Velocity (miles/hour)')                             # add x-axis label
plt.ylabel('B2/m (in MKS units of inverse meters)')                                             # add y-axis label
plt.grid(True)

#
# Set the expanded plot limits according to the extremes of the data points
#
maximumPositionY = 1.2*max(yplot)     # getting 20% extra room the top for the text documentation
plt.xlim(minSpeedMPH, maxSpeedMPH)
plt.ylim(0.0,maximumPositionY)

plt.legend(loc=0)
plt.title('B2/m Value For A Baseball')

#
# The sequence order of the above plt commands is irrelevant
#

#
# Give the command to display the plot
#
plt.show()