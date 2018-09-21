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
dt = 0.027
t = []
err_EX0 = []
err_ER = []

# read t and EX0
f = open('./fort.3115.sym')
flines = f.readlines()
f.close()
for l in flines:
    d = l.split()
    t.append(dt * float(d[0]))
    err_EX0.append(float(d[3]))

# read ER
f = open('./fort.3117.sym')
flines = f.readlines()
f.close()
for l in flines:
    d = l.split()
    err_ER.append(float(d[3]))

# Plot
plot(t,err_EX0,'-.',color='g',label=r"${\parallel E_0^X \parallel}^2$")
plot(t,err_ER,'--',color='r',label=r"${\parallel E_R \parallel}^2$")
plot([40,1000],[0.01,0.01*4/100],'-',color='k',label=r"$\frac{1}{T}$ and $\frac{1}{T^2}$")
plot([40,1000],[0.01,0.01*4*4/100/100],'-',color='k')
text(500,0.001,r"$T^{-1}$")
text(500,0.0001,r"$T^{-2}$")

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
axis([10.,dt*210000,0.000000001,0.02])
xscale('log')
yscale('log')
xlabel(r"$T$")
ylabel(r"Squared residuals")
legend(bbox_to_anchor=(0.75,0.95),numpoints=1)

# Save to file
savefig("streamwise_convergence.png",bbox_inches='tight')
savefig("streamwise_convergence.pdf",bbox_inches='tight')
