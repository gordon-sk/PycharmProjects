import os
import matplotlib.pyplot as plt
import novas.compat.eph_manager as eph_manager
import numpy as np
import time


def main():
    filename = "JPLEPH"
    file = eph_manager.ephem_open(filename)
    print(file)
    print((file[1]-file[0])/365.24)

    time.sleep(.1)
    input("What dates " + \
          "are you: ")

    print(eph_manager.planet_ephemeris((2444911.0, .1), center=2, target=9))

    eph_manager.ephem_close()

def year_shti():
    year = input("what year : ")
    while type(year) is not int and str(year).__len__() is not 4:
        print("Nice try lil shit")


if __name__ == '__main__':
    main()
    print("done")
