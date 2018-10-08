# Gordon Kiesling
# ASTR 3800 -- Structure Formation -- Dr. Berlind
# Assignment 4 --
# Tuesday, 21 Mar 2017
# Written in Python 2.7.10 to PEP 8 standards
#
from __future__ import division
from scipy.integrate import quad
import matplotlib.pyplot as mp
import math

redshift_min = 0.0
redshift_max = 3.0
redshift_step_size = .01
redshift_array = []
for x in range(int(1/redshift_step_size * redshift_max)):   # to be iterated over for each function
    redshift_array.append(redshift_step_size*x)

c = 2.99792458 * pow(10, 5)     # km/s
h = 1.0                         # little h - Hubble
H0 = 100.0*h                    # Hubble constant
DH = c/H0                       # Hubble distance
tH = 9.78 * pow(10, 9)          # age of the universe * little h

# universes 1, 2, 3 and 4
# omega m, omega lambda, w
universe_1 = [1.0, 0.0, 0.0]
universe_2 = [.25, .75, -1.0]
universe_3 = [.25, .75, -.8]
universe_4 = [.25, .75, -1.2]
universe_list = [universe_1, universe_2, universe_3, universe_4]


# Main -- calls and graphs each function
# takes universe list and redshift array, no output
#
def main():
    DC_list = comoving_distance(universe_list)
    luminosity_distance(redshift_array, DC_list)
    angular_diameter(redshift_array, DC_list)
    comoving_volume(DC_list)
    lookback_time(redshift_array, universe_list)


# Comoving Distance
# performs the function on each redshift value, then graphs it
# returns the list of DC functions for use by other functions
#
def comoving_distance(universes):
    DC_list = []
    for universe in universes:
        DC_array = []
        for z in redshift_array:
            Dc = DH * quad(lambda z_prime: 1.0/math.sqrt((universe[0] * pow(1 + z_prime, 3) + universe[1] * pow(1 + z_prime, 3*(1 + universe[2])))), 0, z)[0]
            DC_array.append(Dc)
        DC_list.append(DC_array)

    mp.plot(redshift_array, DC_list[0], 'r', label="Universe 1 - No dark energy")
    mp.plot(redshift_array, DC_list[1], 'g', label="Universe 2 - Real dark energy")
    mp.plot(redshift_array, DC_list[2], 'b', label="Universe 3 - High w")
    mp.plot(redshift_array, DC_list[3], 'c', label="Universe 4 - Low w")
    mp.grid(True)
    mp.legend(loc=0)
    mp.xlabel("Redshift (z)")
    mp.ylabel("Comoviding Distance as function of z (h^-1Gpc)")
    mp.title("Comoving distances for four universe types")
    mp.show()
    return DC_list


# Luminosity Distance
# Takes redshifts and DC values for each universe
# Graphs the results, no output
#
def luminosity_distance(redshifts, DC_list):
    DL_list = []
    for universe in DC_list:
        Dl_array = []
        for z in range(redshifts.__len__()):
            Dl_array.append(universe[z] * (1 + redshifts[z]))
        DL_list.append(Dl_array)
    mp.plot(redshift_array, DL_list[0], 'r', label="Universe 1 - No dark energy")
    mp.plot(redshift_array, DL_list[1], 'g', label="Universe 2 - Real dark energy")
    mp.plot(redshift_array, DL_list[2], 'b', label="Universe 3 - High w")
    mp.plot(redshift_array, DL_list[3], 'c', label="Universe 4 - Low w")
    mp.grid(True)
    mp.legend(loc=0)
    mp.xlabel("Redshift (z)")
    mp.ylabel("Luminosity Distance as function of z (h^-1Gpc)")
    mp.title("Luminosity distances for four universe types")
    mp.show()


# Angular diameter
# Takes redshifts and DC values for each universe
# Graphs the function, no output
#
def angular_diameter(redshifts, DC_list):
    DA_list = []
    for universe in DC_list:
        DA_array = []
        for z in range(redshifts.__len__()):
            DA_array.append(universe[z] / (1 + redshifts[z]))
        DA_list.append(DA_array)
    mp.plot(redshift_array, DA_list[0], 'r', label="Universe 1 - No dark energy")
    mp.plot(redshift_array, DA_list[1], 'g', label="Universe 2 - Real dark energy")
    mp.plot(redshift_array, DA_list[2], 'b', label="Universe 3 - High w")
    mp.plot(redshift_array, DA_list[3], 'c', label="Universe 4 - Low w")
    mp.grid(True)
    mp.legend(loc=0)
    mp.xlabel("Redshift (z)")
    mp.ylabel("Angular Diameter as function of z (h^-1Gpc)")
    mp.title("Angular Diameter for four universe types")
    mp.show()


# Comoving Volume
# takes DC values
# Graphs result of function, no output
#
def comoving_volume(DC_list):
    Cv_list = []
    for universe in DC_list:
        CV_array = []
        for z in universe:
            CV_array.append(4.0 * math.pi / 3.0 * pow(z, 3))
        Cv_list.append(CV_array)
    mp.plot(redshift_array, Cv_list[0], 'r', label="Universe 1 - No dark energy")
    mp.plot(redshift_array, Cv_list[1], 'g', label="Universe 2 - Real dark energy")
    mp.plot(redshift_array, Cv_list[2], 'b', label="Universe 3 - High w")
    mp.plot(redshift_array, Cv_list[3], 'c', label="Universe 4 - Low w")
    mp.grid(True)
    mp.legend(loc=0)
    mp.xlabel("Redshift (z)")
    mp.ylabel("Comoving Volume of a sphere (h^-3Gpc^3")
    mp.title("Comoving Volumes for four universe types")
    mp.show()


# Lookback time
# Graphs function with only redshifts
# no output
#
def lookback_time(redshifts, universes):
    tL_list = []
    for universe in universes:
        tL_array = []
        for z in redshifts:
            tL = tH * quad(lambda z_prime: 1.0 / ((1.0 + z_prime) * math.sqrt((universe[0] * pow(1 + z_prime, 3) + universe[1]
                                                                        * pow(1 + z_prime, 3*(1 + universe[2]))))), 0, z)[0]
            tL_array.append(tL)
        tL_list.append(tL_array)
    mp.plot(redshift_array, tL_list[0], 'r', label="Universe 1 - No dark energy")
    mp.plot(redshift_array, tL_list[1], 'g', label="Universe 2 - Real dark energy")
    mp.plot(redshift_array, tL_list[2], 'b', label="Universe 3 - High w")
    mp.plot(redshift_array, tL_list[3], 'c', label="Universe 4 - Low w")
    mp.grid(True)
    mp.legend(loc=0)
    mp.xlabel("Redshift (z)")
    mp.ylabel("Lookback Time as function of z (Gyr)")
    mp.title("Lookback Time for four universe types")
    mp.show()


# Main executes
main()
