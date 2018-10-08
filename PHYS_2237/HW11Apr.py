from __future__ import division
import math
import matplotlib.pyplot as plt
import time

def main():
    part_1()
    part_2()


def part_1():
    # given a function of the form y = a + b * ln(x),
    # b is given by n * sigma i to N of yi * ln(xi) minus sigma i to N yi * sigma i to N ln(xi)
    # over n * sigma i to N of ln(xi)**2 - (sigma i to N of ln(xi))**2
    #
    # a is given by sigma i to N of yi - b * sigma i to N of ln(xi) / n
    # C(L) = f(L)
    # L, max C value
    data = [[5, .948], [10, 1.258], [15, 1.597], [20, 1.667], [25, 1.828], [30, 2.057], [35, 1.894], [40, 2.045]]
    for x in data:
        x[1] = x[1] / (x[0]**2)
    n = 8
    bt1 = 0
    bt2 = 0
    bt3 = 0
    bb1 = 0
    bb2 = 0
    at1 = 0
    at2 = 0
    Cvals, Lvals = [], []
    for x in data:
        bt1 += x[1] * math.log(x[0])
        bt2 += x[1]
        bt3 += math.log(x[0])
        bb1 += math.log(x[0])**2
        bb2 += math.log(x[0])
        at1 += x[1]
        at2 += math.log(x[0])
        Cvals.append(x[1])
        Lvals.append(x[0])
    b = ((n * bt1) - (bt2 * bt3)) / ((n * bb1) - (bb2 ** 2))
    a = (at1 - (b * at2)) / n

    x = 1
    plots, Ls = [], []
    while x <= 40:
        plots.append(a + b * math.log(x))
        Ls.append(x)
        x += 1
    print a, b

    plt.plot(Ls, plots, '-g', label="fitted function")
    plt.plot(Lvals, Cvals, "r", label="data points")
    plt.text(12.5, .045, "Fitted equation is C/N: = " + str(round(a, 3)) + " " + str(round(b, 3)) + " * ln(L)", fontsize=15)
    plt.xlabel("Lattice size L")
    plt.ylabel("C/N")
    plt.title("C/N ~ L data and fitted function")
    plt.legend(loc=0)
    plt.grid(True)
    plt.show()


def part_2():
    filename = "HWpt2.txt"
    print "\n\n"
    data_file = open(filename)
    data = data_file.read()
    data_file.close()
    lines = data.splitlines()
    data = []
    for line in lines:
        data.append(line.split(" "))
    del data[0], data[0]
    # each list in data contains a list of the following values:
    # [Temp, Energy, Energy error, heat capacity C per spin particle]
    # lattice size was 40x40

    cMax, Tmax = 0, 0
    for set in data:
        del set[4]
        for x in range(set.__len__()):
            set[x] = float(set[x])
    for set in data:
        if set[3] > cMax:
            cMax = set[3]
            Tmax = set[0]
    print cMax, Tmax
    low_vals, high_vals = [], []
    for set in data:
        if set[0] < Tmax:
            low_vals.append([float(set[0]), float(set[3])])
        elif set[0] > Tmax:
            high_vals.append([float(set[0]), float(set[3])])

    # T is x is x[0], C is y is x[1]
    n = high_vals.__len__()
    bt1 = 0
    bt2 = 0
    bt3 = 0
    bb1 = 0
    bb2 = 0
    at1 = 0
    at2 = 0
    Cvals, Tvals = [], []
    for x in high_vals:
        print x
        bt1 += x[1] * math.log(x[0])
        bt2 += x[1]
        bt3 += math.log(x[0])
        bb1 += math.log(x[0])**2
        bb2 += math.log(x[0])
        at1 += x[1]
        at2 += math.log(x[0])
        Cvals.append(x[1])
        Tvals.append(x[0])
    b = ((n * bt1) - (bt2 * bt3)) / ((n * bb1) - (bb2 ** 2))
    a = (at1 - (b * at2)) / n

    x = Tmax
    plots, Ls = [], []
    while x <= 5.0:
        plots.append(a + b * math.log(x))
        Ls.append(x)
        x += .01
    print a, b

    plt.plot(Ls, plots, '-g', label="fitted function")
    plt.plot(Tvals, Cvals, "r", label="data points")
    plt.text(12.5, .045, "Fitted equation is C/N: = " + str(round(a, 3)) + " " + str(round(b, 3)) + " * ln(L)", fontsize=15)
    plt.xlabel("Lattice size L")
    plt.ylabel("C/N")
    plt.title("C/N ~ L data and fitted function")
    plt.legend(loc=0)
    plt.grid(True)
    plt.show()


main()
