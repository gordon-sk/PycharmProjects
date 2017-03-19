# Gordon Kiesling
# Lab 5 Analysis -- Gamma Spectroscopy
# Written in Python 2.7.10 to PEP8 standards... mostly
# Modern Physics Lab -- Dr. Velkomsky
#
#
import matplotlib.pyplot as mp
import math


# Main method calls the helper methods
#
def main():
    names = "Co60_GSK_MM.txt Na22_GSK_MM.txt Cs137_Bare_GSK_MM.txt Cs137_Al_TOP_GSK_MM.txt " \
            "Cs137_Al_Bottom_GSK_MM.txt Cs137_Lead_Top_GSK_MM.txt Cs137_GSK_MM_Lead_Bottom.txt".split(" ")
    count = 0
    for x in names:
        count += 1
        data, calib = open_file(x)
        print "File " + x + " opened and data read"
        graphable_data = condense(data, calib)
        print "Data condensed, graphing..."
        graphing(graphable_data, x, count)
        print "done\n"


# Generic data-reading method
#
def open_file(filename):
    data_file = open(filename)
    readable = data_file.read().splitlines()
    data_file.close()
    galaxies = []
    calibration_data = readable[readable.__len__()-4].split(" ")
    for x in range(12, readable.__len__()-17):
        sub = readable[x].split(' ')
        while sub.__contains__(''):
            sub.remove('')
        galaxies.append(int(sub[0]))  # Adding the sublist to our master list
    return galaxies, calibration_data


# Folding 16,0000 or so data points into 2048
#
def condense(data, calib):
    condensed_data = {}
    for x in range(0, data.__len__()-8, 8):
        sum = 0
        for y in range(8):
            sum += data[x+y]
        condensed_data[float(calib[0])+float(calib[1])*x] = sum
    # Logging the values so we can read the graphs a little better later on
    for x in condensed_data.keys():
        try:
            condensed_data[x] = math.log(condensed_data[x])
        except ValueError:
            condensed_data[x] = 0
    return condensed_data


# Graphing the data
#
def graphing(final_data, name, count):

    # We throw out the very low energy data, because really, who cares
    # It also screws up my graph and I'm not about that
    data = {}
    for x in final_data.keys():
        if x > 15:
            data[x] = final_data[x]

    # Just 6 lines of code for writing the title, nothing to see here
    name = name.split("_")
    try:
        name.remove("GSK")
        name.remove("MM.txt")
        name.remove("MM")
        name.remove(".txt")
    except ValueError:
        pass
    title = ""
    for x in name:
        title += (str(x) + " ")
    if name.__contains__("Bottom.txt"):
        title = "Cs137 Lead Bottom"

    # Graph that guy UP
    mp.bar(data.keys(), data.values(), color='g', width=.01, edgecolor='g')
    mp.xlim(0, 2050)
    mp.xlabel("Energy in KeV")

    xticks = []
    for x in range(0, 20, 2):
        xticks.append(x*100)

    mp.xticks(xticks)
    mp.ylim(0, max(data.values()) + 3)
    mp.ylabel("Counts, log scale")
    mp.title("Figure " + str(count) + ": " + title)
    mp.grid(True)
    mp.show()


main()
