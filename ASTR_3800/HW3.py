# Gordon Kiesling
# ASTR 3800 -- Structure Formation -- Dr. Berlind
# Assignment 3 -- xi(r) plotting
# Tuesday, 21 Mar 2017
# Written in Python 2.7.10 to PEP 8 standards
#
from __future__ import division
import math
import matplotlib.pyplot as mp
from scipy.spatial import ckdtree as kd
import sys
# a variable that won't change no matter what
c = 2.99792458 * pow(10, 5)   # km/s


# MAIN METHOD
# A 'table of contents' method that organizes our data files and calls the necessary helper methods
# No input, output is detailed in helper method headers
#
def main():
    # Opening all of the data
    SDSS_Mr21_rspace = open_files("SDSS_Mr21_rspace.txt")   # [RA, dec, z]
    SDSS_Mr20_rspace = open_files("SDSS_Mr20_rspace.txt")   #
    SDSS_Mr20_zspace = open_files("SDSS_Mr20_zspace.txt")   #
    SDSS_random = open_files("SDSS_random.txt")             #
    DM = open_files("DM.txt")                               # [x, y, z] coordinates
    DM_random = open_files("DM_random.txt")                 #

    # Converting SDSS spherical coordinate data to cartesian
    SDSS_random = cartesian_conversion(SDSS_random)
    SDSS_Mr21_rspace = cartesian_conversion(SDSS_Mr21_rspace)   # Converting to xyz
    SDSS_Mr20_rspace = cartesian_conversion(SDSS_Mr20_rspace)   #
    SDSS_Mr20_zspace = cartesian_conversion(SDSS_Mr20_zspace)   #

    xi21, error21 = xi_function(SDSS_Mr21_rspace, SDSS_random)                  # r space
    xi20, error20 = xi_function(SDSS_Mr20_rspace, SDSS_random)                  # r space
    xi20_zspace, error20zspace = xi_function(SDSS_Mr20_zspace, SDSS_random)           # z space

    # PART 1
    # Performing part 1, two graphs and three data sets along with the random data
    #
    # Plot a
    graphing_a(xi21, xi20)
    # Plot b
    graphing_b(xi20, xi20_zspace)

    # PART 2
    # Performing part 2
    #
    xi_DM, DM_error = xi_function(kd.cKDTree(DM), kd.cKDTree(DM_random))
    bias_20 = bias_function(xi20, xi_DM)
    bias_21 = bias_function(xi21, xi_DM)
    # plot c
    graphing_c(bias_20, bias_21)

    # PART 3
    # Extra credit
    graphing_EC(xi20, xi21, error20, error21)


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
# Inputs are two cKDTrees of real data and a randomly computed set of galaxy points
# Outputs a dictionary of r correlating to xi(r), as well as a dictionary of r correlating to poisson data bar values
#
def xi_function(data, random):
    ND = data.data.__len__() * 1.0
    NR = random.data.__len__() * 1.0
    bins = []
    for x in range(16):
        r = -1 + x * (2.301 / 15)
        bins.append(pow(10.0, r))
    xi_graph = {}
    error_bars = {}
    DD_former, RR_former = [], []
    for x in range(bins.__len__()):
        DD = ((data.count_neighbors(data, bins[x])) - data.data.__len__()) / 2
        RR = ((random.count_neighbors(random, bins[x])) - random.data.__len__()) / 2
        for y in DD_former:
            DD -= y
        for z in RR_former:
            RR -= z
        DD_former.append(DD)
        RR_former.append(RR)
        try:
            xi_value = pow(NR / ND, 2) * (1.0 * DD / RR) - 1
        except ZeroDivisionError:
            xi_value = 0            # if RR(r) = 0
        xi_graph[bins[x]] = xi_value
        error_bars[bins[x]] = error(DD, RR, xi_value)
    return xi_graph, error_bars


# ERROR
# A function to compute Poisson error values
# Called internally via the xi_function method, takes DD, RR and xi values and outputs a float error value
#
def error(DD, RR, xi):
    try:
        sigma_over_xi = math.sqrt(pow((math.sqrt(DD) / DD), 2) + pow((math.sqrt(RR) / RR), 2))
    except ZeroDivisionError:
        sigma_over_xi = 0
    error = sigma_over_xi * xi
    return error


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
    if xi_gal.__len__() != xi_DM.__len__():
        print "error: xi function tables are unequal"
        sys.exit()
    else:
        xi_gal_k = sorted(xi_gal.keys())
        xi_gal_vals = [xi_gal[xx] for xx in xi_gal_k]
        DM_keys = sorted(xi_DM.keys())
        DM_vals = [xi_DM[xx] for xx in DM_keys]
        bias = {}
        for x in range(xi_gal_vals.__len__()):
            bias[x] = math.sqrt(xi_gal_vals[x] / DM_vals[x])
        return bias


# LOGGING
# Just logs both keys and values in each dictionary to be graphed
#
def logging(data):
    returned_data = {}
    for x in data.keys():
        if data[x] == 0:
            pass
        else:
            returned_data[math.log(x, 10)] = math.log(data[x], 10)
    return returned_data


# GRAPHING
# A set of helper methods for graphing the data post-analysis
#
def graphing_a(graph1, graph2):
    graph1 = logging(graph1)
    graph2 = logging(graph2)
    # Sort keys
    graph1_keys = sorted(graph1.keys())
    graph1_vals = [graph1[xx] for xx in graph1_keys]
    graph2_keys = sorted(graph2.keys())
    graph2_vals = [graph2[xx] for xx in graph2_keys]
    #
    mp.title("GRAPH A: logxi(r) vs logr for the two real-space galaxy samples, M<-20 and M<-21")
    mp.ylabel("log(xi(r))")
    mp.xlabel("log(r) (Mpc)")
    mp.plot(graph1_keys, graph1_vals, '-go', label="M<-21")
    mp.plot(graph2_keys, graph2_vals, '-ro', label="M<-20")
    mp.grid(True)
    mp.legend(loc=0)
    mp.show()


#
def graphing_b(graph1, graph2):
    graph1 = logging(graph1)
    graph2 = logging(graph2)
    graph1_keys = sorted(graph1.keys())
    graph1_vals = [graph1[xx] for xx in graph1_keys]
    graph2_keys = sorted(graph2.keys())
    graph2_vals = [graph2[xx] for xx in graph2_keys]
    mp.title("GRAPH B: logxi(r) vs logr for real space vs z space, M<-20")
    mp.ylabel("log(xi(r))")
    mp.xlabel("log(r) (Mpc)")
    mp.plot(graph1_keys, graph1_vals, '-go', label="xi(r) real space")
    mp.plot(graph2_keys, graph2_vals, '-ro', label="xi(r) z space")
    mp.legend(loc=0)
    mp.grid(True)
    mp.show()


#
def graphing_c(graph1, graph2):
    mp.title("GRAPH C: bias(r) vs log(r) for M<-20, M<-21")
    mp.ylabel("bias(r))")
    mp.xlabel("log(r) (Mpc)")
    graph1_keys = sorted(graph1.keys())
    graph1_vals = [graph1[xx] for xx in graph1_keys]
    graph2_keys = sorted(graph2.keys())
    graph2_vals = [graph2[xx] for xx in graph2_keys]
    for x in range(graph1_keys.__len__()):
        if graph1_keys[x] != 0:
            graph1_keys[x] = math.log(graph1_keys[x], 10)
        if graph2_keys[x] != 0:
            graph2_keys[x] = math.log(graph2_keys[x], 10)
    mp.plot(graph1_keys, graph1_vals, '-go', label="bias of M<-20")
    mp.plot(graph2_keys, graph2_vals, '-bo', label="bias of M<-21")
    mp.grid(True)
    mp.legend(loc=0)
    mp.show()


#
def graphing_EC(graph1, graph2, error1, error2):
    del error1[.1]
    del error2[.1]
    del graph1[.1]
    del graph2[.1]
    graph1_keys = sorted(graph1.keys())
    graph1_vals = [graph1[xx] for xx in graph1_keys]
    graph2_keys = sorted(graph2.keys())
    graph2_vals = [graph2[xx] for xx in graph2_keys]
    error1_keys = sorted(error1.keys())
    error1_vals = [error1[xx] for xx in error1_keys]
    error2_keys = sorted(error2.keys())
    error2_vals = [error2[xx] for xx in error2_keys]
    mp.title("EXTRA CREDIT GRAPH: logxi(r) vs logr for the two real-space galaxy samples, M<-20 and M<-21,"
             " with Poisson error bars")
    mp.grid(True)
    mp.ylabel("log(xi(r))")
    mp.xlabel("log(r) (Mpc)")
    mp.legend(loc=0)
    ax = mp.subplot()
    ax.set_xscale("log")
    ax.set_yscale("log")
    mp.errorbar(graph1_keys, graph1_vals, yerr=error1_vals, label="M<-20")
    mp.errorbar(graph2_keys, graph2_vals, yerr=error2_vals, label="M<-21")
    mp.show()


# Execution method
main()
print "Finished"
