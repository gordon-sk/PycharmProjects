# Gordon Kiesling
# ASTR 3800
# Dr. Berlinds
# Homework Assignment 1: Passive stellar color / luminosity evolution
# Written in Python 2.7.10

import matplotlib                       # For plotting in Python
import matplotlib.pyplot as mp          #
import math                             # for doing math
import scipy.integrate as integrate     # for doing harder math
matplotlib.use("TkAgg")                 # Special case Mac OS requirement


# IMF for luminosity
# High mass stars will die first. The integral will be taken with a receding upper limit

# PART 1
# First we declares some lists, our lower mass limit, and our initial luminosity
#
time_pt1, mag = [], []
lower_mass = .08                                                      # lower mass limit defined
L0 = integrate.quad(lambda M: pow(M, -2.35+3.5), lower_mass, 50)[0]   # Luminosity at t=0

#
# Now we solve L(t) between t=0 and t=10Gyr, and save the values in the lists declared above
#
t = 0
while t < 10:
    if t == 0:
        upper_mass = 50
    else:
        upper_mass = math.pow(t/10.0, -.4)                      # setting our upper mass limit as a function of time
    x = (integrate.quad(lambda M: pow(M, -2.35+3.5), lower_mass, upper_mass)[0])/L0              # setting the ratio
    mag.append(-2.5*math.log(x, 10))                                                             # Convert to abs mag
    time_pt1.append(round(t, 3))                                                                 # both are appended
    t += .01

#
# Now, we plot it
#
mp.plot(time_pt1, mag, "b", label="L(t)")
mp.xlabel("Time in GYrs")
mp.ylabel("M(t) - M(0)")
mp.xlim(-1, 11)
mp.ylim(0, 10)
mp.grid(True)
mp.legend(loc=0)
mp.title("Luminosity (Abs Mag) of a star group as a function of time")
mp.show()


# PART 2
# g-r color as a function of time.... A three step process:
# 1. Calculate what fraction of stars exist in each mass bin for masses rising by .01Msol steps
# 2. Calculate luminosity of a star of that mass, and multiply step 1 by that value
# 3. Multiply that total luminosity by the color function
# Finally we redo this in time iterations up to 10 GYr
t = 0
total_color, time_pt2 = [], []
tp = lambda M: math.pow(M, 1.15)
bt = lambda R: math.pow(R, 1.15) * math.pow(10, -.4*(math.log(((R+2.0)/R)) - .65))

while t < 10:
    if t == 0:
        upper_mass = 50
    else:
        upper_mass = math.pow(t/10.0, -.4)   # setting our upper mass limit as a function of time
    bottom = integrate.quad(tp, lower_mass, upper_mass)[0]
    top = integrate.quad(bt, lower_mass, upper_mass)[0]
    total_color.append(-2.5*math.log(top/bottom, 10))
    time_pt2.append(round(t, 3))
    t += .01

#
# Plotting pt 2
#
mp.plot(time_pt2, total_color, "r", label="g-r(M) over time")
mp.xlabel("Time in GYrs")
mp.ylabel("g-r(M)")
mp.grid(True)
mp.legend(loc=0)
mp.title("Color (g-r) of a star group as a function of time")
mp.show()


# PART 3
# g-r plotted against absolute magnitude
# 10, 100 MYrs and 2, 5, and 10 GYrs are plotted

time = [.01, .1, 1.0, 2.0, 5.0, 10.0]
pt3_mag, pt3_color = [], []

for x in range(6):
    pt3_mag.append(mag[int(time[x]/.01)])               # isolate our special point values
    pt3_color.append(total_color[int(time[x]/.01)])

fig = mp.figure()
mp.plot(total_color, mag, "g")
mp.plot(pt3_color, pt3_mag, 'ro')
for x in range(time.__len__()):                         # This loop labels the special values
    mp.annotate(str(time[x]) + " GYrs", xy=([pt3_color[x], pt3_mag[x]+.5]))
mp.legend(loc=0)
mp.grid(True)
mp.xlabel("color")
mp.ylabel("magnitude (absolute)")
mp.title("Color vs Magnitude over time")
mp.legend()
mp.show()
