# Gordon Kiesling
# ASTR 3800
# Dr. Berlind
# Homework Assignment 2: Redshift Galaxy Survey Data Analysis
# Written in Python 2.7.10
# 23 February 2017
# Written to PEP 8 standard for python grammar & syntax... mostly

import matplotlib.pyplot as mp  # For plotting in Python
import math                     # For doing the calculations


# a variable that won't change no matter what
c = 2.99792458*pow(10, 5)   # km/s


# The table of contents
def main():
    data = open_file("SDSS_DR7.dat")
    part_A(data)
    part_B(data)
    Part_C(data)
    part_D_E(data, -20)       # Being code efficient in part D
    part_D_E(data, -19)
    part_D_E(data, -18)
    part_F(data, .1)


# A helper method for reading the data and parsing it into a list of lists,
# where each sublist is the collected data for a specific galaxy in the survey.
#
def open_file(filename):
    data_file = open(filename)
    readable = data_file.read().splitlines()  # Get the data and cut it into chunks by line
    data_file.close()  # One of my profs once told me if we didn't include this, it was sloppy code, so here it is
    galaxies = []  # Our master list of data
    for x in range(readable.__len__()):  # Here we separate galaxy data
        sub = readable[x].split(' ')  # splitting space-separated data
        while sub.__contains__(''):  # Getting rid of random empty spaces
            sub.remove('')
        for y in range(sub.__len__()):  # Converting string values to floats
            sub[y] = float(sub[y])
        galaxies.append(sub)  # Adding the sublist to our master list
    return galaxies
    # in the end, each sublist has format [Right Ascension, Declination, Redshift, Gband Abs Mag, Rband Abs Mag]


# PART A
# Plotting Declination vs Right Ascension and redshift vs brightness
# galaxies: The complete list of 550166 galaxy data lines
#
def part_A(galaxies):
    RA, dec, z, Mr = [], [], [], []
    for x in range(galaxies.__len__()):  # Splits up the data into nice separate lists
        temp = galaxies[x]
        RA.append(temp[0])
        dec.append(temp[1])
        z.append(temp[2])
        Mr.append(temp[4])

    # first scatter plot -- position of galaxies in sky
    mp.scatter(RA, dec, s=.001)
    mp.xlabel("Right Ascension")
    mp.ylabel("Declination")
    mp.title("Position of galaxies surveyed in sky")
    mp.grid(True)
    mp.show()

    # Second scatter plot -- redshift vs absolute magnitude
    mp.scatter(z, Mr, s=.01, c="red")
    mp.gca().invert_yaxis()
    mp.xlabel("redshift")
    mp.ylabel("Abs Mag, red band")
    mp.xlim(-.05, .5)
    mp.ylim(-10, -30)
    mp.title("Redshift vs absolute magnitude in the absolute magnitude of red band light")
    mp.grid(True)
    mp.show()


# PART B
# Plotting g-r color distribution, finding out how many galaxies are blue
# galaxies: The complete list of 550166 galaxy data lines
#
def part_B(galaxies):
    color_distro = dict()
    count = 0
    for x in range(galaxies.__len__()):
        color_val = round(galaxies[x][3] - galaxies[x][4], 4)
        if color_distro.has_key(color_val):
            color_distro[color_val] += 1
        else:
            color_distro[color_val] = 0
        if color_val > .75:
            count += 1
    blue_gals = str(round(100.0 * count / galaxies.__len__(), 3))
    print blue_gals + "% of galaxies surveyed are blue."
    # The value printed above is 58.625%

    mp.bar(color_distro.keys(), color_distro.values(), width=.01, color='g')
    mp.xlabel("g-r color")
    mp.figtext(.2, .8, blue_gals + "% of galaxies surveyed are blue.")
    mp.ylabel("Number of galaxies")
    mp.xlim(0, 1.5)
    mp.ylim(0, 21000)
    mp.show()


# PART C
# Plotting the r-band luminosity function... this method is used for multiple parts
# galaxies: The data list of galaxies to be analyzed
# dMr: The width of the bins that red-band Abs Mag will be placed in, val=.1
# z: the redshift, used to find the extent of space observed and thus the volume, val=.1
#
def Part_C(galaxies):
    r_band_lum_function(galaxies, .1, .1, "Part C")     # the helper method


# PARTs D AND E
# Creating volume limited subsamples and running the numbers on them
# Our magnitude limit is imported as lower_mag and used to determine which data we pull
# from the galaxies master list
#
def part_D_E(galaxies, lower_mag):

    max_z, sub_data = 0, []
    for x in galaxies:      # Finding our highest redshift value
        if x[4] == lower_mag:
            if x[2] > max_z:
                max_z = x[2]
    total_count, blue_count = 0, 0
    for x in galaxies:      # Counting with respect to our highest redshift value
        if x[2] < max_z:
            total_count += 1
            sub_data.append(x)
            if x[3] - x[4] > .75:
                blue_count += 1

    Vol = ((4.0/3) * math.pi * pow(max_z * c*pow(10, -2), 3)) * (2.295/(4 * math.pi))    # Same as part C

    print "For galaxies more luminous than " + str(lower_mag) + " in the r-band:"
    print "Redshift bound is z=" + str(max_z)
    print "Volume contained therein is " + str(Vol) + " Mpc^3/h^3"
    print "The number of galaxies in this volume limited sample is " + str(total_count)
    print str(100*round(blue_count*1.0/total_count, 3)) + "% of those galaxies are blue\n"

    r_band_lum_function(sub_data, .1, max_z, "Part E, lower mag limit = " + str(lower_mag))  # The helper method


# PART F
# Basically the helper method r_band_lum_function, but we calculate a different volume for each
# bin based on the average redshift
# Bin size remains imported at .1
# galaxies is the total sample of redshift survey data
# z is not imported, it is pulled from the total data for each bin and used to calculate an avg per bin
#
def part_F(galaxies, dMr):

    bins = {}       # declaring a dictionary object
    for x in range(galaxies.__len__()):
        val = galaxies[x][4]    # our r-band Magnitude
        val = round(val - round(val % dMr, 2), 1)
        # finding what bin it should be in, ex: 20.43 would be converted to 20.4 and count towards the 20.4-20.5 bin
        if bins.has_key(val):
            bins[val].append(galaxies[x][2])    # we keep the data for now, instead of just iterating over a count
        else:
            bins[val] = []  # unlike the helper method, we use a list instead of an int to accumulate z values for later

    for x in bins.keys():
        gal_count = bins[x].__len__()   # number of galaxies in the bin
        avg_counter = 0
        if gal_count == 0:
            bins.__delitem__(x)     # getting rid of empty data points
        else:
            for y in range(gal_count):
                avg_counter += bins[x][y]   # calculating the average redshift, pt 1
            Vol = ((4.0/3)*math.pi*pow((avg_counter/gal_count)*c*pow(10,-2), 3))*(2.295/(4*math.pi)) # integrating avg z
            bins[x] = math.log(((gal_count/Vol)/dMr), 10)   # do the math, make it plottable

    # Finally, we plot what we've got
    mp.plot(bins.keys(), bins.values(), 'ro')
    mp.title("Part F: Volume-weighted r-band Luminosity Function")
    mp.xlabel("r-band Abs Mag")
    mp.ylabel("Log(dn/dMr)")
    mp.grid(True)
    mp.show()


# A helper method to graph the r-band luminosity function
# takes a list of lists, with each sublist being individual galaxy data
# Also takes a specified binwidth (dMr) and redshift (z)
# Finally, for organizational purposes the string 'part' is used to title graphs
#
def r_band_lum_function(galaxies, dMr, z, part):
    # equations used:
    # Volume = (4.0/3.0)*pi*r^3
    # z = vr / c
    # vr = H0 * r
    Vol = ((4.0/3) * math.pi * pow(z * c * pow(10, -2), 3)) * (2.295/(4 * math.pi))
    # equation for volume with radius subbed -- c is converted to km/s, and volume is
    # multiplied by the ratio of sky observed to total sky at the end

    # First we find our Mr values, and place them into bins of width dMr, counting them by iterating
    # onto their associated values in the bins dictionary
    bins = {}
    for x in range(galaxies.__len__()):
        val = galaxies[x][4]
        val = round(val - round(val % dMr, 2), 1)
        if bins.has_key(val):
            bins[val] += 1
        else:
            bins[val] = 0
    # Now that bins[x] represents number of galaxies in a specific bin of red band magnitude,
    # we divide that by volume to get number density, divide that by binwidth, and log it
    for x in bins.keys():
        try:
            bins[x] = math.log((bins[x]/Vol)/dMr, 10)
        except ValueError:
            if bins[x] == 0:
                bins.__delitem__(x)

    # Finally it is graphed.
    mp.plot(bins.keys(), bins.values(), 'ro')
    mp.xlabel("Abs Mag red band")
    mp.ylabel("Log(dn/dMr)")
    mp.title(part + ": r-band Luminosity Function")
    mp.grid(True)
    mp.show()


# Finally, we call it all and watch it run
#
if __name__ == '__main__':
    main()
    print "Finished"
