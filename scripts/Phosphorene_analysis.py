"""
Python analysis file for phosphorene results.

Author:       Denis Gagnon                  <denisg6@hotmail.com>
"""

import numpy as np
import argparse

import matplotlib as mpl
mpl.use('Agg') # No display
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Plot parameters
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{fourier}"]
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['figure.figsize'] = (4,4)

# Command line arguments
parser = argparse.ArgumentParser(description='Plot tight-binding probabilities')
parser.add_argument('--no-polygon', dest='polygon', action='store_false',
                    help='Disable hexagon in plot')
parser.set_defaults(polygon=True)

args = parser.parse_args()

# Load momentum map data
kx = np.loadtxt('xvec.dat')
ky = np.loadtxt('yvec.dat')
kx_mat,ky_mat=np.meshgrid(kx,ky)
prob_mat = np.loadtxt('probability.dat')

# Plot transition probability (v.1)
f = plt.figure()
ax1 = f.add_subplot(111)
ax1.pcolormesh(kx_mat,ky_mat,np.flipud(prob_mat),cmap='viridis',shading='gouraud')
ax1.set_xlabel(r'$k_x [\AA^{-1}]$')
ax1.set_ylabel(r'$k_y [\AA^{-1}]$')
ax1.set_aspect(1.0)

# Save figure
f.tight_layout()
plt.savefig('probability.png',dpi=300)

# Load vector potential data
gdata = np.loadtxt('potential.dat')

# Plot vector potential
mpl.rcParams['figure.figsize'] = (4,3)
f = plt.figure()
ax1 = f.add_subplot(111)
ax1.plot(gdata[:,0], gdata[:,1], label=r'$G_x$', color='k', ls='-')
ax1.plot(gdata[:,0], gdata[:,2], label=r'$G_y$', color='k', ls=':')
ax1.set_xlabel(r'$t$ [s]')
ax1.legend(loc='best')

f.tight_layout()
plt.savefig('potential.png',dpi=300)
