"""
 File name: project4.py
 Author: Gordon Kiesling
 VUnetid: Kiesligs
 Email: Gordon.s.kiesling@vanderbilt.edu
 Class: CS2204
 Date: 3/2/16
 Honor statement: I have neither given nor received unauthorized aid on this work
 Assignment Number: 4
 Description: Given a 3D MRI image, normalize intensities, 
 display the image, generate a histogram and CDF, use them to generate
 a histogram-equalized image, and display the new image.

This is another scripting problem, so you will not be using classes
for this project. Keep in mind then that you will not need to include
class-specific syntax like built-in functions or the self keyword.

You should not be editing your original image when calculating the
HEQ. The HEQ should be a stored in a new matrix.
"""


#
# Imports
import nibabel as nib
import matplotlib.pylab as pl
import numpy as np
import time

np.set_printoptions(threshold=10000, linewidth=10000, suppress=True)


##
# Timing info
localtime = time.asctime(time.localtime(time.time()))
print("Starting time is : ", localtime)


##
# Load the MR volume. It should be stored the same directory as project4.py
img = nib.load('KKI2009-42-MPRAGE.nii')

##
# Make sure it's reasonable
img.get_shape()
img.get_data_dtype()

##
# Convert the MR volume data to a 3D numpy array
data = img.get_data()
m, n, p = data.shape   # get the size of the 3 dimensions

#
# MR data is floating point with strange intensities. Convert to a 0...255
# grayscale image and get the shape
# YOUR CODE HERE (My code is 3 lines) #

for x in range(n):
    for y in range(p):
        i_n = data[85,x,y]
        data[85, x, y] = int(round((((data[85, x, y] - data.min()) / (data.max() - data.min())) - data.min()) * 255))

##
# Plot the middle slice as a reality check
pl.imshow(data[85, :, :], cmap=pl.cm.gray)
pl.show()

##
# Compute the histogram. Okay, this is really too easy: Python provides
# an iterator, but there you have it. I suppose I should check to make
# sure that the data is >= 0 and <= 255, but I am not going to do that.
#
# This will take a minute or two on a reasonable volume
localtime = time.asctime(time.localtime(time.time()))
print("Histogram start is : ", localtime)
#
# YOUR CODE HERE (My code is 3 lines) #
H = np.zeros((256))

for x in range(n-1):
    for y in range(p-1):
        H[data[85,x,y]] += 1

# Timing info
localtime = time.asctime(time.localtime(time.time()))
print("Histogram end is : ", localtime)

#
# Compute the cumulative distribution function. This is simply the integral of the histogram. If
# the histogram were a probability, this would integrate to 1. Instead it should sum to the total
# number of pixels (voxels) in the image.
# YOUR CODE HERE (My code is 5 lines)

CDF = np.zeros((256))
for x in range(1,n):
    for y in range(1,x):
        CDF[x] += H[y]

localtime = time.asctime(time.localtime(time.time()))

##
# Compute the histogram equalization. The function below produces a new image such that the CDF
# is linearized across the value range. This will be computed into a new array that is the same
# size/shape as the original data array. Call the new array: heq_image
# Timing info
localtime = time.asctime(time.localtime(time.time()))
print("HEQ start is : ", localtime)
##

heq_image = np.zeros((256,256))

min = 0
for x in range(n-1):
    if CDF[x] > 0:
       min = CDF[x]
       break
print(min)

for x in range(n-1):
    for y in range(p-1):
        heq_image[x, y] = int(round(((CDF[data[85,x,y]]-min)/(CDF[n-1] - min)) * 255))

##
# Timing info
localtime = time.asctime(time.localtime(time.time()))
print("Computation done at : ", localtime)

##
# Show the middle slice
# heq_image is a placeholder, replace with your variable name
pl.imshow(heq_image[:, :], cmap=pl.cm.gray)
pl.show()

# It's useful to see what we've done. That is to compute the histogram and CDF
# of the new image.
#
# Numpy has a histogram command and Matplotlib has a hist command.
# Only for this new calculation may you use them. It is up to you
# to look into the documentation to find the syntax.
##
# YOUR CODE HERE (My code is 8 lines) #

# numpy's histogram method wasn't working so I just used my own... can't trust the system
new_H = np.zeros((256))
for x in range(n-1):
    for y in range(p-1):
        new_H[heq_image[x,y]] += 1

new_CDF = np.zeros((256))
for x in range(1,n):
    for y in range(1,x):
        new_CDF[x] += new_H[y]

##

# Now plot them
# hist, newhist are placeholders for your original and new histogram
# replace them with whichever variable names you used.
x = np.arange(256)
print("check out these dope graphs")
pl.semilogy(x, H, 'r*', x, new_H, 'b+', linewidth=2)
pl.legend(['original histogram', 'equalized histogram'])
pl.show()

# cdf, newcdf are placeholders for your original and new cdfs
# replace them with whichever variable names you used.
pl.plot(x, CDF, 'r--', x, new_CDF, 'b-', linewidth=2)
pl.legend(['original cdf', 'equalized cdf'], loc=4)
pl.show()
