# Gordon Kiesling
# ASTR 3800 Final Project
# Friends of Friends Dark Matter Halo Simulation
# Written in Python 2.7.10 to PEP 8 standards
#
#

from __future__ import division
import math
import periodic_kdtree as kd
import matplotlib.pyplot as plt

h = 1.0  # Hubble parameter
mass = 1.4 * pow(10, 10) / h  # solar masses
file_title = "DM_2p8M.txt"


# Main method
# no inputs, no outputs
# Defines global variables and executes helper methods
def main():
    # Step 1, read our data
    data = data_reader(file_title)

    # Some variables we'll be using all over the place
    global L, N, r_step
    L = 141.3 / h  # Mpc
    N = data.__len__()  # a count
    r_step = .2 * (L / pow(N, 1.0 / 3.0))  # .2 * L / N^1/3
    print "r_step is", r_step, "Mpc"

    # The rest of the program executes
    halos = halo_counting(data, r_step)
    output_writing(halos)
    plot1(halos)
    plot2(halos)


# method data_reader
# Input is the name of the file to be read
# output is a list of lists, where each sublist is the 3d coordinate data of
#
def data_reader(filename):
    data_file = open(filename)
    readable = data_file.read().splitlines()  # Get the data and cut it into chunks by line
    data_file.close()  # close memory leak
    particles = []  # Our master list of data
    for x in range(readable.__len__()):  # Here we separate galaxy data
        sub = readable[x].split(' ')  # splitting space-separated data
        while sub.__contains__(''):  # Getting rid of random empty spaces
            sub.remove('')
        for y in range(sub.__len__()):  # Converting string values to floats
            sub[y] = float(sub[y])
        particles.append(sub)  # Adding the sublist to our master list
    return particles
    # in the end, each sublist has format [x-coord, y-coord, z-coord]


# method halo_counting
# Uses PeriodicCKDTree to execute a friends-of-friends code in O(nlogn)
# Inputs are the raw data and the mean separation of particles, r_step
# returns a list of lists of lists where each sublist is the list of coordinates for particles in a specific halo
# and each sub-sublist is the x,y,z coordinate specs
def halo_counting(data, r_step):
    boundaries = [L, L, L]  # Setting boundaries for periodicCKDTree
    tree = kd.PeriodicCKDTree(boundaries, data)  #
    groups_list = []  # The list of lists for our final output
    for particle in data:  # Going over every particle in the data
        if particle is not None:  # Particles already slotted into Halos are replaced with a None object
            group = []  # Starting a fresh group for a new, untouched particle
            indexes = tree.query_ball_point(particle, r=r_step)  # returns a list of the indexes in data
            if indexes.__len__() == 1:  # If there are no particles within range to make a halo,
                group.append(data[indexes[0]])  # we simply call this particle the halo and quit,
                data[indexes[0]] = None  # then replace its spot in data with None
            else:
                del indexes[0]  # Don't include itself
                for x in indexes:
                    extra_indexes = tree.query_ball_point(data[x], r=r_step)  # finding more particles in range
                    if extra_indexes.__len__() > 1:  # if there is any particle in range,
                        for k in extra_indexes:
                            if not indexes.__contains__(k):  # and we haven't already counted it,
                                indexes.append(k)  # add it to the end of the list,
                for x in indexes:  # so we can call query on it, too
                    group.append(data[x])  # adding each particle to this group
                    data[x] = None  # removing particle from initial data list
            groups_list.append(group)

    # the following was just for my own reference
    count = 0
    for x in groups_list:
        if x.__len__() >= 10:
            count += 1
    print "There are", count, "halos with 10 or more particles"
    return groups_list


# output writing method
# takes data as input, outputs a text file with specs listed below
#
def output_writing(halos):
    print "starting writing"
    filler = file_title.split(".")[0]  # specifying which file was read
    writeup = open("Final_Project_Data_Writeup_file_" + filler + ".txt", 'w')  # writing a text file
    writeup.write("List of halos and their relevant data from raw dark matter data analysis\n")
    writeup.write("Format: Mass, center of mass coordinates, root mean square length\n")
    writeup.write("There are %.1e" % N + " particles in a box with length " + str(L) + " Mpc,\n")
    writeup.write("and there are " + str(halos.__len__()) + " halos in the data file " + file_title + "\n\n")
    for x in halos:
        writeup.write("%.3e" % (mass * x.__len__()) + " ")  # mass is just number of particles x particles mass
        center = center_of_mass(x)  # center of mass helper method
        writeup.write("[" + str(center[0]) + "  " + str(center[1]) + "  " + str(center[2]) + "] ")  # writing w/o ,
        writeup.write(str(round(root_mean_square(x, center), 3)) + "\n")  # root mean square helper method
    writeup.close()  # closing data leakage
    print "ending writing"


# First plotting method
# radius vs mass
#
def plot1(halos):
    masses, radii = [], []
    for x in halos:
        masses.append(x.__len__() * mass)
        radii.append(root_mean_square(x, center_of_mass(x)))
    plt.loglog(radii, masses, "go", label="Radius vs masses")
    plt.title("Radius vs Mass, dark matter halo data: " + file_title)
    plt.grid(True)
    plt.legend(loc=0)
    plt.xlabel("rms radius (Mpc)")
    plt.ylabel("mass (solar masses)")
    filler = file_title.split(".")[0]
    plt.savefig("figure-1 " + filler + ".png")


# second plotting method
# Mass vs number density of halos with mass greater than Mass
#
def plot2(halos):
    masses = {}
    for x in halos:
        x_mass = x.__len__() * mass
        num_dens = 0
        if not masses.keys().__contains__(x_mass):
            for halo in halos:
                if halo.__len__() * mass > x_mass:
                    num_dens += 1
        masses[x_mass] = num_dens / pow(L, 3)
    graph1_keys = sorted(masses.keys())
    graph1_vals = [masses[xx] for xx in graph1_keys]
    plt.loglog(graph1_keys, graph1_vals, 'go-')
    plt.grid(True)
    plt.title("Cumulative Halo Mass Function")
    plt.xlabel("Mass bins (Msol)")
    plt.ylabel("Number of Halos with mass greater than Mass")
    filler = file_title.split(".")[0]
    plt.savefig('figure-2' + filler + '.png')


# Root mean square helper method
# returns the rms radius value for a collection of points
def root_mean_square(halo, center):
    summed = 0
    for x in halo:
        distance = math.sqrt(pow(x[0] - center[0], 2) + pow(x[1] - center[1], 2) + pow(x[2] - center[2], 2))
        summed += pow(distance, 2)
    return math.sqrt(summed / halo.__len__())


# center of mass helper method
# returns the center of mass xyz coordinate for a halo
# simplified by all particles having the same mass
#
def center_of_mass(halo):
    x_coord, y_coord, z_coord = 0, 0, 0
    for particle in halo:
        x_coord += particle[0]
        y_coord += particle[1]
        z_coord += particle[2]
    x_coord /= halo.__len__()
    y_coord /= halo.__len__()
    z_coord /= halo.__len__()
    return [round(x_coord, 2), round(y_coord, 2), round(z_coord, 2)]


main()
