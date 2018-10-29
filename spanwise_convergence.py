#!/usr/bin/env/python
# -*- coding: utf-8 -*-

##################################################################
# Flageul and Tiselj
# https://doi.org/
##################################################################

# import modules
import os, sys
from numpy import *
from pylab import *

# define time step and empty arrays
dt = 0.027*10
t = []
Z_conv = []
Z_diff = []
err_EZ0 = []

# read t, Z-conv, Z-diff and EX0
f = open('./fort.3119.sym')
flines = f.readlines()
f.close()
for l in flines:
    d = l.split()
    t.append(dt * float(d[0]))
    Z_diff.append(float(d[1]))
    Z_conv.append(float(d[2]))
    err_EZ0.append(float(d[3]))

# Plot
plot(t,err_EZ0,'-.',color='g',label=r"${\parallel E_0^Z \parallel}^2$")
plot(t,Z_conv,'--',color='r',label=r"$Z-Conv$")
plot(t,Z_diff,'--',color='b',label=r"$Z-Diff$")
plot([1000,10000],[0.0002,0.0002/10],'-',color='k',label=r"$\frac{1}{T}$ and $\frac{1}{T^2}$")
plot([1000,10000],[0.00002,0.00002/10/10],'-',color='k')
text(4500,0.00007,r"$T^{-1}$")
text(4000,0.0000002,r"$T^{-2}$")

# Graph settings
matplotlib.rc('figure', figsize=(5.,4.13))
matplotlib.rc('text', usetex = True)
size=16
size_legend=14
size_label=20
linewidth=1.5
markersize=10
framealpha=1.
matplotlib.rc('lines', linewidth=linewidth,markersize=markersize)
matplotlib.rc('font', size=size)
matplotlib.rc('axes', labelsize=size_label, titlesize=size)
matplotlib.rc('legend', fontsize=size_legend, framealpha=framealpha)
#
axis([100.,dt*210000,0.000000001,0.02])
xscale('log')
yscale('log')
xlabel(r"$T^+$")
ylabel(r"Squared residuals")
legend(bbox_to_anchor=(0.43,0.45),numpoints=1)

# Save to file
savefig("spanwise_convergence.png",bbox_inches='tight')
savefig("spanwise_convergence.pdf",bbox_inches='tight')
