# Gordon Kiesling
# ASTR 3800 -- Structure Formation
# Assignment 3 -- xi(r) plotting
# Tuesday, 21 Mar 2017
# Written in Python 2.7.10 to PEP 8 standards
#

import math
import matplotlib.pyplot as mp
from scipy.spatial import ckdtree as kd
import sys

# a variable that won't change no matter what
c = 2.99792458*pow(10, 5)   # km/s


# MAIN METHOD
# A 'table of contents' method that organizes our data files and calls the necessary helper methods
# No input, output is detailed in helper method headers
#
def main():
    # Opening all of the data
    SDSS_Mr21_rspace = open_files("SDSS_Mr21_rspace.txt")   # [RA, dec, z]
    SDSS_Mr20_rspace = open_files("SDSS_Mr20_rspace.txt")
    SDSS_Mr20_zspace = open_files("SDSS_Mr20_zspace.txt")
    SDSS_random = open_files("SDSS_random.txt")
    DM = open_files("DM.txt")                               # [x, y, z] coordinates
    DM_random = open_files("DM_random.txt")

    # Converting SDSS spherical coordinate data to cartesian
    SDSS_random = cartesian_conversion(SDSS_random)
    SDSS_Mr21_rspace = cartesian_conversion(SDSS_Mr21_rspace)
    SDSS_Mr20_rspace = cartesian_conversion(SDSS_Mr20_rspace)
    SDSS_Mr20_zspace = cartesian_conversion(SDSS_Mr20_zspace)

    rest = False
    if not rest:
        # PART 1
        # Performing part 1, two graphs and three data sets along with the random data
        #
        xi21 = xi_function(SDSS_Mr21_rspace, SDSS_random)                   # r space
        xi20 = xi_function(SDSS_Mr20_rspace, SDSS_random)                  # r space
        # Plot a
        graphing(xi21, xi20)
        # Plot b
        xi20_zspace = xi_function(SDSS_Mr20_zspace, SDSS_random)                  # z space
        graphing(xi20, xi20_zspace)

    # PART 2
    # Performing part 2
    #
    xi_DM = xi_function(DM, DM_random)                                  # part two
    bias_20 = bias_function(xi20, xi_DM)
    bias_21 = bias_function(xi21, xi_DM)
    # plot c
    graphing(bias_20, bias_21)

    # PART 3
    # Extra credit: We'll see if we get to this
    # part_3()                                               # part three -- extra credit



# A helper method for reading the data and parsing it into a list of lists,
# where each sublist is the collected data in a specific file.
#
def open_files(filename):
    data_file = open(filename)
    readable = data_file.read().splitlines()  # Get the data and cut it into chunks by line
    data_file.close()
    galaxies = []
    for x in range(readable.__len__()):  # Here we separate galaxy data
        sub = readable[x].split(' ')  # splitting space-separated data
        while sub.__contains__(''):  # Getting rid of random empty spaces
            sub.remove('')
        for y in range(sub.__len__()):  # Converting string values to floats
            sub[y] = float(sub[y])
        galaxies.append(sub)  # Adding the sublist to our master list
    return galaxies
    # in the end, each sublist has format [Right Ascension, Declination, Redshift]


# XI FUNCTION
# A helper method that computes the xi function for survey data
# Input is a dataset and a randomly computed set of galaxy points
# Outputs a dictionary of r correlating to xi(r)
#
def xi_function(data, random):

    ND = data.data.__len__() * 1.0
    NR = random.data.__len__() * 1.0

    bins = []
    for x in range(16):
        r = -1 + x * (2.301 / 15)
        bins.append(pow(10.0, r))

    xi_graph = {}
    DD_former = [0]
    RR_former = [0]
    for x in range(bins.__len__()):
        DD = 0
        RR = 0
        DD = ((data.count_neighbors(data, (bins[x])) - data.data.__len__()) / 2) - DD_former[x]
        RR = ((random.count_neighbors(random, (bins[x])) - random.data.__len__()) / 2) - RR_former[x]
        DD_former.append(DD)
        RR_former.append(RR)
        try:
            xi_graph[bins[x]] = pow(ND / NR, 2) * (1.0 * DD / RR) - 1
        except ZeroDivisionError:
            xi_graph[bins[x]] = 0
    return xi_graph

    # KDtree.countpairs takes the form      number of pairs = a kdtree.count_pairs(another kdtree, radius)
    # a kdtree is created by listing tuples, ie tree = (0,0), (1,1), (2,2), etc
    # This is infuriating bc i just want it to count the pairs around one set of points


# CARTESIAN CONVERSION
# convert from spherical coordinates to cartesian coordinates
# Data in spherical coordinates appears as [right ascension, declination, redshift] or [RA, D, z]
# right ascension is phi, declination is theta, z is range
# x = r sin theta cos phi
# y = r sin theta sin phi
# z = r cos theta = zc/H cos theta
#
def cartesian_conversion(data):
    cartesian_data = []
    Ngal = data.__len__()
    for x in data:
        sub = []
        r = x[2] * c / 100
        sub.append(r * math.sin(math.radians(90-x[1])) * math.cos(math.radians(x[0])))
        sub.append(r * math.sin(math.radians(90-x[1])) * math.sin(math.radians(x[0])))
        sub.append(r * math.cos(math.radians(90-x[1])))
        cartesian_data.append(sub)
    sub = []
    for x in cartesian_data:
        sub.append([x[0], x[1], x[2]])
    return kd.cKDTree(sub)


# BIAS FUNCTION
# Computes and returns the bias function with xi(r) for galactic data and xi(r) for dark matter data
#
def bias_function(xi_gal, xi_DM):
    print xi_gal.__len__(), xi_DM.__len__()



# GRAPHING
# A helper method for graphing the data post-analysis
#
def graphing(graph1, graph2):
    mp.figure(1)
    mp.subplot(211)
    mp.plot(graph1.keys(), graph1.values(), 'go')
    mp.subplot(212)
    mp.plot(graph2.keys(), graph2.values(), 'bo')
    mp.xlabel("r (Mpc)")
    mp.ylabel("xi function")
    mp.grid(True)
    mp.show()


# Execution method
main()
print "Finished"