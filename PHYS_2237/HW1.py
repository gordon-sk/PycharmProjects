import matplotlib
matplotlib.use('TkAgg')
import math
import matplotlib.pyplot as plt

def main():
    x = 1000
    t = 0
    xplot = []
    yplot = []
    tau = 1

    while t<5:
        yplot.append(x*math.exp(-t/tau))
        t+=.05
        xplot.append(t)

    plt.plot(xplot,yplot,'ro',label="Decay Chain")
    plt.legend(loc=0)
    plt.show()

main()