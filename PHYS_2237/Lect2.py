#  plotCircleV1.py, version December 1, 2016
#  Type command line  python plotCircleV1.py
#  There are no input options on the command line

import matplotlib                   # library of plotting functions
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt     # assign the main plotting function in the library to the name plt
import numpy as np                  # library of widely used numerical functions in python


def his_1():

    # Python script to generate a set of (x,y) points which fall on a circle of a given radius
    # The points are first plotted discretely
    # A following plot command is given to plot the points as a continuous line
    # Other plot commands add labels to the x and y axes, add a figure title, put in a legend, and add grid lines
    #

    #
    # Generate the set of (x,y) points which fall on a circle
    # The radius of the circle and the number of points are hardcoded in this script
    #
    radius = 5.0   # radius of circle
    npltX = 100    # number of x points to generate; note each x value produces two y values

    xStart = -radius
    xStep = 2*radius/npltX  # the iteration over x starts at xStart and increases by xStep each time

    #
    # Declare two empty lists to hold the (x,y) data points
    #
    xplot = []
    yplot = []

    radiusSquared = radius*radius
    #
    # First iteration loop to produce the (x,y) points in the top half of the circle, starting at the left
    #
    xValue = xStart
    while xValue <= radius:                                 # note the : mark at the end
        xplot.append(xValue)                                # an indentation must be given, here it is 4 spaces
        yValue = np.sqrt(radiusSquared - xValue*xValue)     # all other statements in the loop must have the same iteration (except comments)
        yplot.append(yValue)    # positive square root for y

        xValue = xValue + xStep                             # the xValue is increased by one step, this is the last statement in the loop

    #
    # Second teration loop to produce the (x,y) points in the bottom half of the circle, starting at the right
    #
    xValue = radius - xStep
    while xValue >= -radius:                                 # note the : mark at the end
        xplot.append(xValue)                                 # an indentation must be given, here it is 4 spaces
        yValue = -np.sqrt(radiusSquared - xValue*xValue)     # all other statements in the loop must have the same iteration (except comments)
        yplot.append(yValue)    # negative square root for y

        xValue = xValue - xStep                             # the xValue is decreased by one step, this is the last statement in the loop

    #
    # Plot this set of (x,y) coordinates as a set of discrete, red data points
    #
    plt.plot(xplot,yplot, 'ro', label='Discrete data points')      # put a label identification for this plot

    #
    # Plot the same set of (x,y) coordinates as a continuous line
    #
    plt.plot(xplot,yplot, label='Continuous line')                 # pat a label identification for this plot

    #
    # Add labels for the x and y coordinates, use a rectangular grid, use expanded ranges for the (x,y) axes
    # Add the legend identifying the plots and put in a title for the figure
    #
    plt.xlabel('x coordinate')                             # add x-axis label
    plt.ylabel('y coordinate')                             # add y-axis label
    plt.grid(True)
    plt.xlim(-1.1*radius,1.1*radius)                       # make the x range 10% larger than the radius
    plt.ylim(-1.1*radius,1.1*radius)                       # make the y range 10% larger than the radius
    plt.legend(loc=0)
    plt.title('Plotting a circle as a set of discrete points and a continuous line')
    #
    # The sequence order of the above plt commands is irrelevant
    #

    #
    # Give the command to display the plot
    #
    plt.show()

his_1()