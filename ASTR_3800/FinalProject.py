# Gordon Kiesling
# ASTR 3800 Final Project
# Friends of Friends Dark Matter Halo Simulation
# Written in Python 2.7.10 to PEP 8 standards
#


from __future__ import division
import matplotlib.pyplot as plt
import math
from scipy.spatial import ckdtree as kd
import time

h = 1.0     # Hubble parameter
mass = 1.4 * pow(10, 10) / h        # solar masses
files = ["DM.txt", "DM_random.txt", "DM_2p8M.txt", "DM_5p5M.txt"]
print "Which file would you like to process? In order of size, choices are:"
for x in range(files.__len__()):
    print str(x) + " " + files[x]
file_title = files[input("Select file number: ")]


def main():
    global r_step, N, L
    data = data_reader(file_title)
    L = 141.3 / h                         # Mpc
    N = data.__len__()
    r_step = .2 * (L / pow(N, 1.0/3.0))  # .2 * L / N^1/3
    print "r_step is", r_step
    halos = brute_force(data, r_step)
    output_writing(halos)
    plots(halos)


def data_reader(filename):
    data_file = open(filename)
    readable = data_file.read().splitlines()    # Get the data and cut it into chunks by line
    data_file.close()                           # close memory leak
    particles = []                              # Our master list of data
    for x in range(readable.__len__()):         # Here we separate galaxy data
        sub = readable[x].split(' ')            # splitting space-separated data
        while sub.__contains__(''):             # Getting rid of random empty spaces
            sub.remove('')
        for y in range(sub.__len__()):          # Converting string values to floats
            sub[y] = float(sub[y])
        particles.append(sub)                   # Adding the sublist to our master list
    return particles
    # in the end, each sublist has format [x-coord, y-coord, z-coord]


def brute_force(data, r_step):
    tree = kd.cKDTree(data)
    groups_list = []
    for particle in data:
        if particle is not None:
            group = []
            indexes = tree.query_ball_point(particle, r_step)      # returns the indexes in data of particles within r_step from the initial particle
            if indexes.__len__() == 1:
                group.append(data[indexes[0]])
                data[indexes[0]] = None
            else:
                del indexes[0]
                for x in indexes:
                    extra_indexes = tree.query_ball_point(data[x], r_step)
                    if extra_indexes.__len__() > 1:
                        for k in extra_indexes:
                            if not indexes.__contains__(k):
                                indexes.append(k)
                for x in indexes:                                # excluding the particle counting itself
                    group.append(data[x])                        # adding each particle to this group
                    data[x] = None                               # removing particle from initial data list
            groups_list.append(group)
    count = 0
    for x in groups_list:
        if x.__len__() >= 10:
            count += 1
    print "There are", count, "halos with 10 or more particles"
    return groups_list


def output_writing(halos):
    writeup = open("Final_Project_Data_Writeup.txt", 'w')
    writeup.write("List of halos and their relevant data from raw dark matter data analysis\n")
    writeup.write("Format: Mass, center of mass coordinates, root mean square length\n")
    writeup.write("There are %.1e" % N + " particles in a box with length " + str(L) + " Mpc,\n")
    writeup.write("and there are " + str(halos.__len__()) + " halos in the data file " + file_title + "\n\n")
    for x in halos:
        writeup.write("%.3e" % (mass * x.__len__()) + " ")
        center = center_of_mass(x)
        writeup.write("[" + str(center[0]) + " " + str(center[1]) + " " + str(center[2]) + "] ")
        writeup.write(str(round(root_mean_square(x, center), 3)) + "\n")
    writeup.close()


def plots(halos):
    print "begin first plot"
    radii, masses = [], []
    for x in halos:
        masses.append(x.__len__() * mass)
        radii.append(root_mean_square(x, center_of_mass(x)))
        if halos[halos.__len__()-1] == x:
            print "finished analyizing plot one data"
    plt.loglog(radii, masses)
    plt.savefig("fuckme.png")


    print "begin second plot"
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
    plt.loglog(graph1_keys, graph1_vals, '-go')
    plt.grid(True)
    plt.title("Cumulative Halo Mass Function")
    plt.xlabel("Mass bins (Msol)")
    plt.ylabel("Number of Halos with mass greater than Mass")
    plt.show()

def root_mean_square(halo, center):
    summed = 0
    for x in halo:
        distance = math.sqrt(pow(x[0] - center[0], 2) + pow(x[1] - center[1], 2) + pow(x[2] - center[2], 2))
        summed += pow(distance, 2)
    return math.sqrt(summed / halo.__len__())


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
