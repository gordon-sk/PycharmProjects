# Gordon Kiesling
# Physics 2237 -- Computational Physics
# Exam 1
# Dr. Charles Maguire
# Written in Python 2.7.10 to PEP 8 standard
#

import matplotlib.pyplot as mp
import math
import numpy as np
from scipy.integrate import odeint


# Main method calls helper methods for each part of the exam
# No inputs, outputs appear on screen or in terminal text box
#
def main():
    part_1()
    part_2()
    part_3()


# PART 1
# Performs part one of the exam -- a position and velocity dependent graph of the function dv/dx is produced
# Constants and functioned all defined inside, no output
#
def part_1():

    # Initial conditions
    timestep = .1  # seconds
    maxT = 10      # seconds
    v0 = 0         # m/s

    # defining our derivative function for odeint
    def fprime(variableList, t):
        v = variableList[0]
        dvdt = 5 - .5*v
        dxdt = v
        return [dvdt, dxdt]

    # Defining a timegrid for odeint
    nTimeSteps = int(maxT/timestep)
    timeGrid = np.linspace(0, maxT, nTimeSteps)
    initialValuesSpeedPosition = [v0, 0.0]

    # Solving the differential equation
    twoSolution = odeint(fprime, initialValuesSpeedPosition, timeGrid)

    # Pulling the data we want
    vel, pos = [], []
    for x in twoSolution:
        vel.append(x[0])
        pos.append(x[1])

    # Plotting
    mp.figure(1)
    mp.subplot(211)
    mp.plot(timeGrid, vel, 'ro', label="v(t)")
    mp.xlabel('Time (s)')
    mp.ylabel('Velocity (m/s)')

    mp.title('Motion Along A Hill With Quadratic Air Resistance')
    mp.grid(True)
    mp.legend(loc=2)

    mp.subplot(212)
    mp.plot(timeGrid, pos, 'bo', label="x(t)")
    mp.grid(True)
    mp.legend(loc=2)
    mp.xlabel('Time (s)')
    mp.ylabel('Position (m)')
    mp.show()

    # A sanity check done with Wolfram Alpha's analytic solution
    # their solution: v = -10e^-.5t + 10
    t = 0
    vel_check, time_check, diff = [], [], []
    while t < maxT:
        t += timestep
        vel_check.append(-10*math.exp(-.5*t)+10)
        time_check.append(t)

    # Plotting sanity check
    mp.plot(time_check, vel_check, 'ro', label="Sanity Check on Analytical Solution")
    mp.legend(loc=0)
    mp.xlabel("Time (s)")
    mp.ylabel("Position (m)")
    mp.xlim(0, 10)
    mp.ylim(0, 10)
    mp.grid(True)
    mp.show()


# PART 2
# Performs a numerical solution for the two coupled differential counting equations
# No input, outputs graph of counts over time
#
def part_2():
    # Initial conditions
    A0 = 5000     # count - N
    alpha = 1.2   # seconds - s
    B0 = 0        # count - N
    beta = .8     # seconds - s

    # Defining timegrid for odeint
    maxT = 10
    timestep = .1
    nTimeSteps = int(maxT/timestep)
    timeGrid = np.linspace(0, maxT, nTimeSteps)
    initial_values = [A0, B0]

    # define the derivatives for odeint
    def fprime(dependentVariables, t):
        NA = dependentVariables[0]
        NB = dependentVariables[1]
        dNAdt = -1*NA/alpha + NB/beta
        dNBdt = -1*NB/beta + NA/alpha
        return [dNAdt, dNBdt]

    # Solving the differential equations
    twoSolutions = odeint(fprime, initial_values, timeGrid)

    # Pulling the data we want
    series_A, series_B = [], []
    for x in twoSolutions:
        series_A.append(x[0])
        series_B.append(x[1])

    # Plotting
    mp.plot(timeGrid, series_A, 'ro', label="Series A")
    mp.plot(timeGrid, series_B, 'bo', label="Series B")
    mp.legend(loc=1)
    mp.xlabel('Time (s)')  # add axis labels
    mp.ylabel('Counts (N)')
    mp.title('Steady-state status of transmutable quantities')
    mp.grid(True)
    mp.legend(loc=1)
    mp.show()


# PART 3
# Projectile motion of a motorcyclist jumping off a ramp
# Differential equation is performed for all angles between 0 and 90 degrees
# Only theta=10 is graphed
#
def part_3():
    # position and dimensional intial values
    x0, y0 = 0, 0
    v0 = 60     # m/s

    # air resistance values
    A = 1       # m^2
    ro = 1.2    # kg/m^3
    mass = 300  # kg
    C = 1
    b2overM = .5*C*ro*A/mass

    # situational values and natural contants
    g = 9.81855392611

    # time
    maxT = 15
    timestep = .01
    nTimeSteps = int(maxT/timestep)
    timeGrid = np.linspace(0, maxT, nTimeSteps)

    # to satisfy the odeint method
    def fprime(variableList, t):
        vx = variableList[0]
        vy = variableList[1]
        vTotal = math.sqrt(pow(vx, 2) + pow(vy, 2))
        dvxdt = -1.0 * b2overM * vTotal * vx
        dvydt = -1.0 * (b2overM * vTotal * vy + g)
        return [dvxdt, dvydt, vx, vy]

    # Defining a dictionary object to save angles and associated ranges
    range_dict = {}
    for angle in range(0, 90, 2):               # Looping over theta
        theta_rads = math.radians(angle)        # Defining theta
        # defining theta-dependent initial values
        initial_values = [v0*math.cos(theta_rads), v0*math.sin(theta_rads), x0, y0]

        # Solving the ODE
        four_solutions = odeint(fprime, initial_values, timeGrid)   # returns [vx, vy, xpos, ypos]

        # pulling the data we're interested in
        xpos, ypos = [], []
        range_marker = False
        for x in range(four_solutions.__len__()):
            xpos.append(four_solutions[x][2])
            ypos.append(four_solutions[x][3])
            if four_solutions[x][3] < 0.0 and not range_marker:
                range_dict[angle] = four_solutions[x-1][2]       # Recording the maximum range before y becomes negative
                range_marker = True
        # Plotting for theta=10
        if angle == 10:
            mp.plot(xpos, ypos, 'bo', label="x(t)")
            mp.title("Projectile Motion with Quadratic Air Resistance")
            mp.xlabel("Horizontal Pos")
            mp.ylabel("Vertical Pos")
            mp.legend(loc=1)
            mp.ylim(0, 10)
            mp.xlim(0, 125)
            mp.grid(True)
            mp.show()

    # finding the maximum range recorded in the range_dict dictionary
    best_range = max(range_dict.values())
    best_angle = 0
    for x in range_dict.keys():
        if range_dict[x] == best_range:
            best_angle = x
    print "Part Three:"
    print "The best range is " + str(best_range) + " meters given by the angle " + str(best_angle) + " degrees"
    print "The angle 10 degrees gives a range of aproximately " + str(range_dict[10]) + " meters"


# Executing the file
main()
