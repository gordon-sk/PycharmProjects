# Gordon Kiesling
# PHYS 2237 Exam 2
# Prof. Charles Maguire
# 31 March 2017
# Python 2.7.10 in PEP 8

import matplotlib
matplotlib.use("TKAgg")

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint
import math as mp
import numpy as np
import sys



def problem_1(planet):

    sunPlot = False
    verbose = False
    nFrames = 500
    noAnimation = False
    doMovie = False
    nInterval = 200
    orbitSun = False

    FOURPISQUARE = 4*(mp.pi)*(mp.pi)
    AU = 1.496e+11                               # Earth's distance from the Sun in meters
    AUCUBE = AU*AU*AU                            # cube of AU since G*Msun is in units of AU**3/year**2
    YEAR = 3.156e+07                             # number of seconds in one year
    YEARSQUARE = YEAR*YEAR                       # square one year
    GSUNMASS = FOURPISQUARE*AUCUBE/YEARSQUARE    # used in potential energy formula
    SUNMASS = 1.991e+30                          # Sun's mass in kg

    maxT = 3
    deltaT = .002

    # 1
    sun_mass = SUNMASS  # kg

    # 2
    mercury_mass = .33011 * pow(10, 24) # kg
    mercury_x0 = .310       # Au
    mercury_vy0 = 12.433     # 2piAU/yr

    # 3
    if planet == 'Jupiter':
        planet_mass = 1898.19 * pow(10, 24) # kg
        planet_x0 = 4.95    # Au
        planet_vy0 = 2.8922       # Au / yr
    else:
        planet_mass = 4.8675 * pow(10, 24) # kg
        planet_x0 = .718    # Au
        planet_vy0 = 7.433  # Au / yr

    mass_sum = sun_mass + mercury_mass + planet_mass
    sun_x0 = (mercury_x0 * mercury_mass + planet_x0 * planet_mass) / mass_sum
    mercury_x0 = (mercury_x0 * (sun_mass + planet_mass) - planet_x0 * planet_mass) / mass_sum
    planet_x0 = (planet_x0 * (sun_mass * mercury_mass) - mercury_x0 * mercury_mass) / mass_sum

    sun_vy0 = -mercury_vy0 * mercury_mass / sun_mass - planet_vy0 * planet_mass / sun_mass

    sun_y0 = 0
    mercury_y0 = 0
    planet_y0 = 0

    sun_vx0 = 0
    mercury_vx0 = 0
    planet_vx0 = 0


    def fDerivative(variableList, t):                  # variableList dummy list array since there is more than one differential equation
        vx1 = variableList[0]                          # Sun speed in the x direction
        vy1 = variableList[1]                          # Sun speed in the y direction
        vx2 = variableList[2]                          # Mercury speed in the x direction
        vy2 = variableList[3]                          # Mercury speed in the y direction
        vx3 = variableList[4]                          # planet speed in the x direction
        vy3 = variableList[5]                          # planet speed in the y direction
        x1 = variableList[6]                           # Sun x coordinate
        y1 = variableList[7]                           # Sun y coordinate
        x2 = variableList[8]                           # Mercury x coordinate
        y2 = variableList[9]                           # Mercury y coordinate
        x3 = variableList[10]                          # planet x coordinate
        y3 = variableList[11]                          # planet y coordinate
        #
        # Compute distance parameters
        #
        d12Sq = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) # square of the distance between Sun and Mercury
        d13Sq = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) # square of the distance between the Sun and the planet
        d23Sq = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) # square of the distance between Mercury and the planet

        d12Cu = d12Sq * mp.sqrt(d12Sq)   # cube of the distance between the Sun and Mercury
        d13Cu = d13Sq * mp.sqrt(d13Sq)   # cube of the distance between the Sun and the planet
        d23Cu = d23Sq * mp.sqrt(d23Sq)   # cube of the distance between Mercury and the planet
        print "solving..."
        #
        # Compute the derivatives
        #
        dvx1dt = -FOURPISQUARE*(x1-x2)*mercury_mass/d12Cu - FOURPISQUARE*(x1-x3)*planet_mass/d13Cu    # the time derivative of velocity v1 in the x direction according to the universal gravity force component
        dvy1dt = -FOURPISQUARE*(y1-y2)*mercury_mass/d12Cu - FOURPISQUARE*(y1-y3)*planet_mass/d13Cu    # the time derivative of velocity v1 in the y direction according to the universal gravity force component
        dvx2dt = -FOURPISQUARE*(x2-x1)*sun_mass/d12Cu - FOURPISQUARE*(x2-x3)*planet_mass/d23Cu    # the time derivative of velocity v2 in the x direction according to the universal gravity force component
        dvy2dt = -FOURPISQUARE*(y2-y1)*sun_mass/d12Cu - FOURPISQUARE*(y2-y3)*planet_mass/d23Cu    # the time derivative of velocity v2 in the y direction according to the universal gravity force component
        dvx3dt = -FOURPISQUARE*(x3-x1)*sun_mass/d13Cu - FOURPISQUARE*(x3-x2)*mercury_mass/d23Cu    # the time derivative of velocity v3 in the x direction according to the universal gravity force component
        dvy3dt = -FOURPISQUARE*(y3-y1)*sun_mass/d13Cu - FOURPISQUARE*(y3-y2)*mercury_mass/d23Cu    # the time derivative of velocity v3 in the y direction according to the universal gravity force component
        netForceX = sun_mass*dvx1dt + mercury_mass*dvx2dt + planet_mass*dvx3dt
        netForceY = sun_mass*dvy1dt + mercury_mass*dvy2dt + planet_mass*dvy3dt
        return [dvx1dt, dvy1dt, dvx2dt, dvy2dt, dvx3dt, dvy3dt, vx1, vy1, vx2, vy2, vx3, vy3]    # return the twelve derivatives as a list object containing eight elements

    initialValuesSpeedsPositions = [sun_vx0, sun_vy0, mercury_vx0, mercury_vy0, planet_vx0, planet_vy0, sun_x0, sun_y0, mercury_x0, mercury_y0, planet_x0, planet_y0]   # starting values of vx1, vy1, vx2, vy2, x1, y1, x2, y2 for the iteration
    nTimeSteps = int(maxT / deltaT)
    timeGrid = np.linspace(0, maxT, nTimeSteps)           # time grid used by odeint method
    twelveSolution = odeint(fDerivative, initialValuesSpeedsPositions, timeGrid)   # odeint returns a list of values which is the RK4 solution
    vx1RK4 = twelveSolution[:,0]          # vx1 function of time obtained with RK4 solution
    vy1RK4 = twelveSolution[:,1]          # vy1 function of time obtained with RK4 solution
    vx2RK4 = twelveSolution[:,2]          # vx2 function of time obtained with RK4 solution
    vy2RK4 = twelveSolution[:,3]          # vy2 function of time obtained with RK4 solution
    vx3RK4 = twelveSolution[:,4]          # vx3 function of time obtained with RK4 solution
    vy3RK4 = twelveSolution[:,5]          # vy3 function of time obtained with RK4 solution

    x1RK4 = twelveSolution[:,6]           # x1 function of time obtained with RK4 solution
    y1RK4 = twelveSolution[:,7]           # y1 function of time obtained with RK4 solution
    x2RK4 = twelveSolution[:,8]           # x2 function of time obtained with RK4 solution
    y2RK4 = twelveSolution[:,9]           # y2 function of time obtained with RK4 solution
    x3RK4 = twelveSolution[:,10]          # x3 function of time obtained with RK4 solution
    y3RK4 = twelveSolution[:,11]          # y3 function of time obtained with RK4 solution


    ################################################################################


    #
    # Set up the radial position and the total velocity arrays for the first planet (Earth default)
    #
    radialPosition2 = []
    totalVelocity2= []
    radialPosition12 = []
    x12Position = []
    y12Position = []

    nTimeStep = 0
    while nTimeStep < nTimeSteps:
        radialPosition2.append(mp.sqrt(x2RK4[nTimeStep]*x2RK4[nTimeStep] + y2RK4[nTimeStep]*y2RK4[nTimeStep]))
        totalVelocity2.append(mp.sqrt(vx2RK4[nTimeStep]*vx2RK4[nTimeStep] + vy2RK4[nTimeStep]*vy2RK4[nTimeStep]))
        r12 = mp.sqrt((x2RK4[nTimeStep]-x1RK4[nTimeStep])*(x2RK4[nTimeStep]-x1RK4[nTimeStep]) + (y2RK4[nTimeStep]-y1RK4[nTimeStep])*(y2RK4[nTimeStep]-y1RK4[nTimeStep]))
        radialPosition12.append(r12)
        x12Position.append(x2RK4[nTimeStep] - x1RK4[nTimeStep])
        y12Position.append(y2RK4[nTimeStep] - y1RK4[nTimeStep])
        nTimeStep += 1

    meanRadialPosition2 = np.mean(radialPosition2)
    stdRadialPosition2 = np.std(radialPosition2)
    meanTotalVelocity2 = np.mean(totalVelocity2)
    stdTotalVelocity2 = np.std(totalVelocity2)

    meanRadialPosition12 = np.mean(radialPosition12)
    stdRadialPosition12 = np.std(radialPosition12)

    #
    # Set up the radial position and the total velocity arrays for the second planet (Jupiter default)
    #
    radialPosition3 = []
    totalVelocity3= []
    radialPosition13 = []
    x13Position = []
    y13Position = []

    nTimeStep = 0
    while nTimeStep < nTimeSteps:
        radialPosition3.append(mp.sqrt(x3RK4[nTimeStep]*x3RK4[nTimeStep] + y3RK4[nTimeStep]*y3RK4[nTimeStep]))
        totalVelocity3.append(mp.sqrt(vx3RK4[nTimeStep]*vx3RK4[nTimeStep] + vy3RK4[nTimeStep]*vy3RK4[nTimeStep]))
        r13 = mp.sqrt((x3RK4[nTimeStep]-x1RK4[nTimeStep])*(x3RK4[nTimeStep]-x1RK4[nTimeStep]) + (y3RK4[nTimeStep]-y1RK4[nTimeStep])*(y3RK4[nTimeStep]-y1RK4[nTimeStep]))
        radialPosition13.append(r13)
        x13Position.append(x3RK4[nTimeStep] - x1RK4[nTimeStep])
        y13Position.append(y3RK4[nTimeStep] - y1RK4[nTimeStep])
        nTimeStep += 1

    meanRadialPosition3 = np.mean(radialPosition3)
    stdRadialPosition3 = np.std(radialPosition3)
    meanTotalVelocity3 = np.mean(totalVelocity3)
    stdTotalVelocity3 = np.std(totalVelocity3)

    meanRadialPosition13 = np.mean(radialPosition13)
    stdRadialPosition13 = np.std(radialPosition13)
    radialPlot = False
    if(radialPlot):
        plt.figure(1)       # start a figure for a single plot of the orbit
        plt.subplot(211)    # this sets the upper half plot for the first planet
        radiusString2 = '%s %5.4f %s %5.4f %s' % ('Radius from c.m. = ', meanRadialPosition2, ' +/- ', stdRadialPosition2, ' AU')
        plt.plot(timeGrid, radialPosition2, label=radiusString2)    # plot as a continuous line
        radiusString12 = '%s %5.4f %s %5.4f %s' % ('Radius from Sun = ', meanRadialPosition12, ' +/- ', stdRadialPosition12, ' AU')
        plt.plot(timeGrid, radialPosition12, label=radiusString12)    # plot as a continuous line
        plt.ylabel('Radius (AU)')
        rMinimum = min(radialPosition2)
        rMaximum = max(radialPosition2)
        yMinimum = 0.995*rMinimum
        yMaximum = rMaximum + 0.6*(rMaximum - rMinimum)
        plt.xlim(0, timeGrid[nTimeSteps -1])
        plt.ylim(yMinimum, yMaximum)
        plt.title('Radius of the First Planet in a Sun + Two Planet System')
        plt.grid(True)
        plt.legend(loc=1)

        plt.subplot(212)    # this sets the lower half plot for second planet
        radiusString3 = '%s %5.4f %s %5.4f %s' % ('Radius from c.m. = ', meanRadialPosition3, ' +/- ', stdRadialPosition3, ' AU')
        plt.plot(timeGrid, radialPosition3, label=radiusString3)    # plot as a continuous line
        radiusString13 = '%s %5.4f %s %5.4f %s' % ('Radius from Sun = ', meanRadialPosition13, ' +/- ', stdRadialPosition13, ' AU')
        plt.plot(timeGrid, radialPosition13, label=radiusString13)    # plot as a continuous line
        plt.xlabel('Time (years)')
        plt.ylabel('Radius (AU)')
        rMinimum3 = min(radialPosition3)
        rMinimum13 = min(radialPosition13)
        rMaximum3 = max(radialPosition3)
        rMaximum13 = max(radialPosition13)
        yMinimum = 0.985*min(rMinimum3, rMinimum13)
        yMaximum = 1.04*max(rMaximum3, rMaximum13)
        plt.xlim(0, timeGrid[nTimeSteps -1])
        plt.ylim(yMinimum, yMaximum)
        plt.title('Radius of the Second Planet in a Sun + Two Planet System')
        plt.grid(True)
        plt.legend(loc=1)

        plt.show()          # show the complete figure
        exit()

    minimumPositionX1 = 1.5*min(x1RK4)
    minimumPositionX2 = 1.5*min(x2RK4)
    minimumPositionX3 = 1.5*min(x3RK4)

    maximumPositionX1 = 1.5*max(x1RK4)
    maximumPositionX2 = 1.5*max(x2RK4)
    maximumPositionX3 = 1.5*max(x3RK4)

    minimumPositionX = min(minimumPositionX1, minimumPositionX2, minimumPositionX3)
    maximumPositionX = max(maximumPositionX1, maximumPositionX2, maximumPositionX3)

    minimumPositionY1 = 1.5*min(y1RK4)
    minimumPositionY2 = 1.5*min(y2RK4)
    minimumPositionY3 = 1.5*min(y3RK4)

    maximumPositionY1 = 1.5*max(y1RK4)
    maximumPositionY2 = 1.5*max(y2RK4)
    maximumPositionY3 = 1.5*max(y3RK4)

    minimumPositionY = min(minimumPositionY1, minimumPositionY2, minimumPositionY3)
    maximumPositionY = max(maximumPositionY1, maximumPositionY2, maximumPositionY3)

    if(orbitSun):
        minimumPositionX13 = 1.5*min(x13Position)
        minimumPositionX12 = 1.5*min(x12Position)

        maximumPositionX13 = 1.5*max(x13Position)
        maximumPositionX12 = 1.5*max(x12Position)

        minimumPositionX = min(minimumPositionX13, minimumPositionX12)
        maximumPositionX = max(maximumPositionX13, maximumPositionX12)

        minimumPositionY13 = 1.5*min(y13Position)
        minimumPositionY12 = 1.5*min(y12Position)

        maximumPositionY13 = 1.5*max(y13Position)
        maximumPositionY12 = 1.5*max(y12Position)

        minimumPositionY = min(minimumPositionY13, minimumPositionY12)
        maximumPositionY = max(maximumPositionY13, maximumPositionY12)


    print "\n Minimum x position for the plot = ", minimumPositionX,
    print ",  maximum x position for the plot = ", maximumPositionX
    print " Minimum y position for the plot = ", minimumPositionY,
    print ",  maximum y position for the plot = ", maximumPositionY

    fig = plt.figure(1)       # start a figure for plots of the two orbits
    ax = fig.add_subplot(111)
    if(orbitSun):
        plt.title('Orbits Of Two Planets About the Sun')
        plt.plot(x12Position, y12Position, label='Orbit of the first planet')   # plot as a continuous line
        plt.plot(x13Position, y13Position, label='Orbit of the second planet')  # plot as a continuous line
    else:
        if(sunPlot):
            plt.title('Orbits of the Sun and Two Planets About the Center-of-Mass')
        else:
            plt.title('Orbits Of Two Planets About the Center-of-Mass')
        plt.plot(x2RK4, y2RK4, label='Orbit of the first planet')   # plot as a continuous line
        plt.plot(x3RK4, y3RK4, label='Orbit of the second planet')  # plot as a continuous line
        if(sunPlot):
            plt.plot(x1RK4, y1RK4, label='Orbit of the Sun')   # plot as a continuous line

    plt.xlabel('x Coordinate (AU)')                             # add axis labels
    plt.ylabel('y Coordinate (AU)')

    #
    # Center of mass point symbol
    #
    plt.scatter(0, 0, c='k', s=20)

    plt.grid(True)
    plt.xlim(minimumPositionX, maximumPositionX)
    plt.ylim(minimumPositionY, maximumPositionY)

    circularOrbit = False
    ellipticalOrbit = True
    if(stdRadialPosition2/meanRadialPosition2 < 1.e-05 and stdTotalVelocity2/meanTotalVelocity2):
        print "\n The orbit of the first planet is found to be circular (eccentricity = 0)"
        ellipticalOrbit = False
        circularOrbit = True
    else:
        print "\n The orbit of the first planet is found to be elliptical (eccentricity > 0)"

    eccentricity = 0
    orbitTime = 0
    numberOfOrbits = 0
    if(circularOrbit):
        orbitTime = 2*mp.pi*meanRadialPosition2/meanTotalVelocity2
        numberOfOrbits = int(maxT/orbitTime)

    lookingForNextPerihelion = True
    lookingForNextAphelion = False
    if(ellipticalOrbit):
        #
        # Algorithm is to accumulate the set of perihelion and aphelion radial values and their time values
        # The times to pass through successive perihelion points are stored from which a mean and standard deviation are computed
        # The first step is to confirm that the initial point is consistent with being a perihelion
        #
        if(radialPosition12[1] < radialPosition12[0]):
            print " Check that the initial position is consistent with being an aphelion is passed"
        else:
            lookingForNextAphelion = True
            lookingForNextPerihelion = False
            print " Looking for first aphelion position"


        radiusPerihelion = []
        radiusAphelion = []
        lastRadialPosition = radialPosition12[0]
        nTimeStep = 1
        timePerihelion = []
        timeAphelion = []
        nPerihelion = 0
        nAphelion = 0
        while nTimeStep < nTimeSteps:
            #
            # As the while loop iterates through the positions at each time step, it will be either looking for the next aphelion or the next perihelion
            #
            newRadialPosition = radialPosition12[nTimeStep]

            if(lookingForNextAphelion):
                #
                # If the new radial position is smaller than the previous radial position, then the next aphelion point has been crossed
                #
                if(newRadialPosition < lastRadialPosition):
                    nAphelion += 1
                    timeValue = 0.5*(timeGrid[nTimeStep] + timeGrid[nTimeStep-1])  # take an average of the current and previous times
                    timeAphelion.append(timeValue)     # store the time value of this aphelion
                    if(nAphelion == 1):
                         if(verbose): print "\n Found first aphelion at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                    else:
                        if(verbose): print "\n Found next aphelion ", nAphelion, " at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
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
                    nPerihelion += 1
                    timeValue = 0.5*(timeGrid[nTimeStep] + timeGrid[nTimeStep-1])  # take an average of the current and previous times
                    timePerihelion.append(timeValue)     # store the time difference last perihelion
                    if(nPerihelion == 1):
                        if(verbose): print "\n Found first Perihelion at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                    else:
                        if(verbose): print "\n Found next perihelion ", nPerihelion, " at position ", newRadialPosition, " at time ", timeGrid[nTimeStep]
                    radiusPerihelion.append(0.5*(lastRadialPosition+newRadialPosition))
                    lookingForNextPerihelion = False
                    lookingForNextAphelion = True

                lastRadialPosition = newRadialPosition

            nTimeStep += 1

        timeDifferencePerihelion = []
        timeDifferenceAphelion = []
        timeDifferenceBoth = []
        nTimesPerihelion = 1
        while nTimesPerihelion < len(timePerihelion):
            timeDifference = timePerihelion[nTimesPerihelion] - timePerihelion[nTimesPerihelion - 1]
            timeDifferencePerihelion.append(timeDifference)
            timeDifferenceBoth.append(timeDifference)
            nTimesPerihelion += 1

        nTimesAphelion = 1
        while nTimesAphelion < len(timeAphelion):
            timeDifference = timeAphelion[nTimesAphelion] - timeAphelion[nTimesAphelion - 1]
            timeDifferenceAphelion.append(timeDifference)
            timeDifferenceBoth.append(timeDifference)
            nTimesAphelion += 1

        if(verbose): print "\n Number of perihelion time differences ", len(timeDifferencePerihelion)
        if(verbose): print "\n Number of aphelion time differences ", len(timeDifferenceAphelion)

        orbitTimeAphelion = 0.0
        if(len(timeDifferenceAphelion) > 0):
            orbitTimeAphelion = np.mean(timeDifferenceAphelion)
            orbitTime = orbitTimeAphelion

        orbitTimePerihelion = 0.0
        if(len(timeDifferencePerihelion) > 0):
            orbitTimePerihelion = np.mean(timeDifferencePerihelion)
            orbitTime = orbitTimePerihelion

        orbitTimeBoth = 0.0
        if(len(timeDifferenceAphelion) > 0 and len(timeDifferencePerihelion) > 0):
            orbitTimeBoth = np.mean(timeDifferenceBoth)
            orbitTime = orbitTimeBoth
            print "\n Perihelion time ", orbitTimePerihelion, "  Aphelion time ", orbitTimeAphelion, " Both time ", orbitTimeBoth

        if(nPerihelion > 0 and nAphelion > 0):
            majorAxisLength = np.mean(radiusAphelion) + np.mean(radiusPerihelion)
            eccentricity = (np.mean(radiusAphelion) - np.mean(radiusPerihelion))/majorAxisLength
            print '%s %5.3f' % (' The orbit eccentricity = ', eccentricity)

        numberOfOrbits = max(nTimesPerihelion, nTimesAphelion)

    print '%s %d' % (' Number of completed orbits = ', numberOfOrbits)
    if(numberOfOrbits > 0):
        print '%s %5.3f %s' % (' Orbit time = ', orbitTime, ' years')

    xTextPosition1 = minimumPositionX + 0.53*(maximumPositionX - minimumPositionX)
    yTextPosition1 = minimumPositionY + 0.10*(maximumPositionY - minimumPositionY)
    mass2_text = ax.text(xTextPosition1, yTextPosition1, '',color='blue')
    mass2_text.set_text('Mass2 = %.3e solar masses' % mercury_mass)
    x20_text = ax.text(xTextPosition1, minimumPositionY + 0.06*(maximumPositionY - minimumPositionY), '',color='blue')
    x20_text.set_text('Initial x2 = %.3f AU' % mercury_x0)
    vy20_text = ax.text(xTextPosition1, minimumPositionY + 0.02*(maximumPositionY - minimumPositionY), '',color='blue')
    vy20_text.set_text('Initial vy2 = %.3f AU/Year' % mercury_vy0)

    xTextPosition2 = minimumPositionX + 0.03*(maximumPositionX - minimumPositionX)
    yTextPosition2 = yTextPosition1
    mass3_text = ax.text(xTextPosition2, yTextPosition2, '',color='green')
    mass3_text.set_text('Mass3 = %.3e solar masses' % planet_mass)
    planet_x0_text = ax.text(xTextPosition2, minimumPositionY + 0.06*(maximumPositionY - minimumPositionY), '',color='green')
    planet_x0_text.set_text('Initial x3 = %.3f AU' % planet_x0)
    vy30_text = ax.text(xTextPosition2, minimumPositionY + 0.02*(maximumPositionY - minimumPositionY), '',color='green')
    vy30_text.set_text('Initial vy3 = %.3f AU/Year' % planet_vy0)

    xTextPosition3 = minimumPositionX + 0.35*(maximumPositionX - minimumPositionX)
    yTextPosition3 = minimumPositionY + 0.30*(maximumPositionY - minimumPositionY)
    xTextPosition4 = xTextPosition3
    yTextPosition4 = minimumPositionY + 0.25*(maximumPositionY - minimumPositionY)

    orbitTime_text = ax.text(xTextPosition3, yTextPosition3, '', color='red')
    orbitTime_text.set_text('Orbit time = %.2f years' % orbitTime)
    eccentricity_text = ax.text(xTextPosition4, yTextPosition4, '', color='red')
    eccentricity_text.set_text('Eccentricity = %.3f' % eccentricity)

    xTextPosition5 = 0.5*(xTextPosition1 + xTextPosition2)
    yTextPosition5 = minimumPositionY + 0.65*(maximumPositionY - minimumPositionY)
    maximumTime_text = ax.text(xTextPosition5, yTextPosition5, '',color='red')
    maximumTime_text.set_text('Maximum time = %.1f years' % maxT)

    xTextPosition6 = 0.5*(xTextPosition1 + xTextPosition2)
    yTextPosition6 = minimumPositionY + 0.60*(maximumPositionY - minimumPositionY)
    timeStep_text = ax.text(xTextPosition6, yTextPosition6, '',color='red')
    timeStep_text.set_text('Time step = %.3e years' % deltaT)

    npts = nTimeStep - 1
    finalSpeed1 = np.sqrt(vx1RK4[npts]*vx1RK4[npts] + vy1RK4[npts]*vy1RK4[npts])
    finalSpeed2 = np.sqrt(vx2RK4[npts]*vx2RK4[npts] + vy2RK4[npts]*vy2RK4[npts])
    finalSpeed3 = np.sqrt(vx3RK4[npts]*vx3RK4[npts] + vy3RK4[npts]*vy3RK4[npts])

    finalSeparation12 = np.sqrt((x1RK4[npts] - x2RK4[npts])*(x1RK4[npts] - x2RK4[npts]) + (y1RK4[npts] - y2RK4[npts])*(y1RK4[npts] - y2RK4[npts]))
    initialSeparation12 = np.sqrt((sun_x0 - mercury_x0)*(sun_x0 - mercury_x0) + (sun_y0 - mercury_y0)*(sun_y0 - mercury_y0))
    finalSeparation13 = np.sqrt((x1RK4[npts] - x3RK4[npts])*(x1RK4[npts] - x3RK4[npts]) + (y1RK4[npts] - y3RK4[npts])*(y1RK4[npts] - y3RK4[npts]))
    initialSeparation13 = np.sqrt((sun_x0 - planet_x0)*(sun_x0 - planet_x0) + (sun_y0 - planet_y0)*(sun_y0 - planet_y0))
    finalSeparation23 = np.sqrt((x2RK4[npts] - x2RK4[npts])*(x2RK4[npts] - x3RK4[npts]) + (y2RK4[npts] - y3RK4[npts])*(y2RK4[npts] - y3RK4[npts]))
    initialSeparation23 = np.sqrt((mercury_x0 - planet_x0)*(mercury_x0 - planet_x0) + (mercury_y0 - planet_y0)*(mercury_y0 - planet_y0))

    initialPE123 = -GSUNMASS*SUNMASS*mercury_mass*sun_mass/(initialSeparation12*AU) - GSUNMASS*SUNMASS*sun_mass*planet_mass/(initialSeparation13*AU) - GSUNMASS*SUNMASS*mercury_mass*planet_mass/(initialSeparation23*AU)
    finalPE123 = -GSUNMASS*SUNMASS*mercury_mass*sun_mass/(finalSeparation12*AU) - GSUNMASS*SUNMASS*sun_mass*planet_mass/(finalSeparation13*AU) - GSUNMASS*SUNMASS*mercury_mass*planet_mass/(finalSeparation23*AU)
    initialEnergy123 = initialPE123 + 0.5*SUNMASS*sun_mass*pow(sun_vy0*AU/YEAR, 2) + 0.5*SUNMASS*mercury_mass*pow(mercury_vy0*AU/YEAR, 2) + 0.5*SUNMASS*planet_mass*pow(planet_vy0*AU/YEAR, 2)
    finalEnergy123 = finalPE123 + 0.5*SUNMASS*sun_mass*pow(finalSpeed1*AU/YEAR, 2) + 0.5*SUNMASS*mercury_mass*pow(finalSpeed2*AU/YEAR, 2) + 0.5*SUNMASS*planet_mass*pow(finalSpeed3*AU/YEAR, 2)

    print "\n Initial speed for sun_mass = ", sun_vy0, " AU/year,  final speed = ", finalSpeed1
    print " Initial speed for mercury_mass = ", mercury_vy0, " AU/year,  final speed = ", finalSpeed2
    print " Initial speed for planet_mass = ", planet_vy0, " AU/year,  final speed = ", finalSpeed3
    print " Initial system energy (1+2+3) = ", initialEnergy123, " (Joules)"
    print " Final system energy (1+2+3) = ", finalEnergy123, " (Joules)"
    print " Fractional energy change = ", (finalEnergy123 - initialEnergy123)/initialEnergy123

    #
    # Add the animation
    #
    line1, = ax.plot([], [], 'o-', lw=4)                                # this is an empty line object used in the initialization of the animation
    line2, = ax.plot([], [], 'o-', lw=4)                                # this is an empty line object used in the initialization of the animation
    line3, = ax.plot([], [], 'o-', lw=4)                                # this is an empty line object used in the initialization of the animation
    timeValue_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='red')   # this is a placeholder for where the time will be updated
    mass2Position_text = ax.text(0.02, 0.90, '', transform=ax.transAxes, color='blue')   # this is a placeholder for where the sun_mass (x,y) will be updated
    mass3Position_text = ax.text(0.02, 0.85, '', transform=ax.transAxes, color='green')   # this is a placeholder for where the mass3 (x,y) will be updated

    def init():                                # this initialization function tells what is being animiated, i.e. a single point "line", and the x and y coordinate displays
        """initialize animation"""
        line2.set_data([], [])
        line3.set_data([], [])
        timeValue_text.set_text('')
        mass2Position_text.set_text('')
        mass3Position_text.set_text('')
        if(sunPlot and orbitSun == False):
            line1.set_data([], [])
            return line1, line2, line3, timeValue_text, mass2Position_text, mass3Position_text
        else:
            return line2, line3, timeValue_text, mass2Position_text, mass3Position_text

    nSkipPoints = int(nTimeSteps/nFrames)
    if(nSkipPoints < 1):
        nSkipPoints = 1
    print "\n nSkipPoints = ", nSkipPoints

    def animate(i):                            # this is the function this being animated, the i index is increased according the range of the frames value below
        """perform animation step"""           # i ranges from 0 to (nFrames - 1) as given in the FuncAnimation function
        global x2RK4, y2RK4, x3RK4, y3RK4, nSkipPoints, orbitTime, orbitSun

        ii = nSkipPoints*i
        if(ii >= nTimeSteps - 1):
            ii = nTimeSteps - 1
        if(orbitSun):
            x2 = x2RK4[ii] - x1RK4[ii]
            y2 = y2RK4[ii] - y1RK4[ii]
        else:
            x2 = x2RK4[ii]
            y2 = y2RK4[ii]
        v2 = np.sqrt(vx2RK4[ii]*vx2RK4[ii] + vy2RK4[ii]*vy2RK4[ii])
        timeValue = ii*deltaT
        orbitNumber = int(timeValue/orbitTime) + 1
        line2.set_data(x2,y2)  # the "line" contains a single data point (x,y) for the first planet
        timeValue_text.set_text('Time = %.2f Years, First Planet Orbit Number %d' % (timeValue,orbitNumber))
        mass2Position_text.set_text('(x2,y2) = (%.2f %.2f) AU, v2 = %.2f AU/Year' % (x2,y2,v2))
        if(orbitSun):
            x3 = x3RK4[ii] - x1RK4[ii]
            y3 = y3RK4[ii] - y1RK4[ii]
        else:
            x3 = x3RK4[ii]
            y3 = y3RK4[ii]
        v3 = np.sqrt(vx3RK4[ii]*vx3RK4[ii] + vy3RK4[ii]*vy3RK4[ii])
        line3.set_data(x3,y3)                     # the "line" contains a single data point (x,y) for the second planet
        mass3Position_text.set_text('(x3,y3) = (%.2f %.2f) AU, v3 = %.2f AU/Year' % (x3,y3,v3))
        if(sunPlot and orbitSun == False):
            x1 = x1RK4[ii]
            y1 = y1RK4[ii]
            line1.set_data(x1,y1)  # the "line" contains a single data point (x,y) for the Sun orbiting about the c.m.
            return line1, line2, line3, timeValue_text, mass2Position_text, mass3Position_text
        else:
            return line2, line3, timeValue_text, mass2Position_text, mass3Position_text

    #
    # The blit=True means that the figure should be updated only for the parts which have changed, i.e. the data point position and the x and y value displays
    #
    if(noAnimation):
        print "\n Animation step is being skipped"
    else:
        print "\n Animation step is being done"
        ani = animation.FuncAnimation(fig, animate, frames=nFrames,
                                      interval=nInterval, blit=True, init_func=init)

    #
    # The ffmpeg binary must be put in the executable PATH, where the binary can be obtained from http://www.ffmpegmac.net/
    #
    if(doMovie):
        print "\n Movie production step is being done"
        FFwriter = animation.FFMpegWriter()
        ani.save('sunTwoPlanetsV1.mp4', writer = FFwriter, fps=30, extra_args=['-vcodec', 'libx264'])
        print " Movie production step is completed"
    else:
        print "\n Movie production step is being skipped"
    plt.show()

problem_1("Jupiter")