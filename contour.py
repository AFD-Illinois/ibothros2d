#!/usr/bin/python
#
# first argument is dump number.
# you must edit the file to change which variable is plotted!
#
# import
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
#
# open data file and get data 
data = np.loadtxt("ibothros2d.dat")
I = data[:,2]
#
# set font size
plt.rc('text', usetex=True) 
plt.rc('font', size = 16)
plt.rc('font', family='sans-serif')
#
# plot the variable
# reshape the array for input to contour
nx = 128
ny = 128
# make axes
# limits of field of view hardcoded; see main.c in ibothros2d
x = np.linspace(-25., 25., num=nx)
y = np.linspace(-25., 25., num=ny)
print len(x), len(y), len(I)
# set plotting array
I0 = 1.6e-8  # RJ limit for B_{\nu} at 10^6 K, 230 GHz
lI = np.log10(I + I0)
pvar = np.reshape(lI, (ny, nx)) 
pvar = np.transpose(pvar)
# prevent dashed negative contours
plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.contour(x, y, pvar, 20, colors = 'k')
# labels
plt.xlabel('$\Delta x \, c^2/(G M)$', weight='bold', fontsize=20)
plt.ylabel('$\Delta y \, c^2/(G M)$', weight='bold', fontsize=20)
# and, as in German, the verb comes last
plt.show()
#
