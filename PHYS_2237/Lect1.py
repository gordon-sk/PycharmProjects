import matplotlib                   # library of plotting functions
matplotlib.use('TkAgg')             # special code to make plots visible on Macintosh system
import matplotlib.pyplot as plt   # assign the main plotting function in the library to the name plt


def his_1():
#  plotPointsV1.py, version December 1, 2016
#  Type command line  python plotPointsV1.py
#  There are no input options on the command line
#  The script plots five pairs of (x,y) data points
#

# Python script to plot a set of 5 (x,y) values as discrete data points
#

#
# List (or array) of x values
#
    x = [1,2,3,4,5]

#
# List of y values
#
    y = [1,4,9,16,25]

#
# Plot this set of (x,y) coordinates as a set of red data points
#
    plt.plot(x,y, 'ro')      # the letter r means red and the letter o means a solid disk plot symbol

#
# Give the command to display the plot
#
    plt.show()

def mine_1():
    x=[1,2,3,4,5,6,7,8,9,10]
    y=[1,4,9,16,25,36,49,64,81,100]
    plt.plot(x,y,'o')
    plt.show()

def his_2():
#  plotPointsV2.py, version December 1, 2016
#  Type command line  python plotPointsV2.py
#  There are no input options on the command line
#

#
# Python script to plot a set of 5 (x,y) values as discrete data points
# A following plot command is given to plot the points as a continuous line
# Other plot commands add labels to the x and y axes, add a figure title, put in a legend, and add grid lines
#

#
# List (or array) of x values
#
    x = [1,2,3,4,5]

#
# List of y values
#
    y = [1,4,9,16,25]

#
# Plot this set of (x,y) coordinates as a set of discrete, red data points
#
    plt.plot(x,y, 'ro', label='Discrete data points')      # put a label identification for this plot

#
# Plot the same set of (x,y) coordinates as a continuous line
#
    plt.plot(x,y, label='Continuous line')                 # pat a label identification for this plot

#
# Add labels for the x and y coordinates, use a rectangular grid, use expanded ranges for the (x,y) axes
# Add the legend identifying the plots and put in a title for the figure
#
    plt.xlabel('x coordinate')                             # add x-axis label
    plt.ylabel('y coordinate')                             # add y-axis label
    plt.grid(True)
    plt.xlim(0,6)
    plt.ylim(0,30)
    plt.legend(loc=0)
    plt.title('Plotting a set of discrete points and a continuous line')
#
# The sequence order of the above plt commands is irrelevant
#

#
# Give the command to display the plot
#
    plt.show()


def mine_2():
    x=[1,2,3,4,5,6,7,8,9,10]
    y=[1,4,9,16,25,36,49,64,81,100]

    plt.plot(x,y,'ro',label='some schtuff')
    plt.plot(x,y,'r',label='same schtuff, more lines')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.xlim(0,6)
    plt.ylim(0,81)
    plt.legend(loc=0)
    plt.title("Goooo fuck urself")




    plt.show()

mine_2()