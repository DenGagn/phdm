"""
Python analysis file for GrapheneSine results.

The excitation
is a few-cycle sine wave with variable duration and amplitude, but zero
carrier-envelope phase

Author:       Denis Gagnon                  <denisg6@hotmail.com>
"""

import numpy as np
import argparse

import matplotlib as mpl
mpl.use('Agg') # No display
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Plot parameters
#mpl.rcParams['ps.useafm'] = True
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r"\usepackage{fourier}"]
#mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['figure.figsize'] = (4,4)

# Command line arguments
parser = argparse.ArgumentParser(description='Plot tight-binding probabilities')
parser.add_argument('--no-polygon', dest='polygon', action='store_false',
                    help='Disable hexagon in plot')
parser.set_defaults(polygon=True)

args = parser.parse_args()

# Load data
kx = np.loadtxt('xvec.dat')
ky = np.loadtxt('yvec.dat')
kx_mat,ky_mat=np.meshgrid(kx,ky)
prob_mat = np.loadtxt('probability.dat')

# Plot transition probability (v.1)
f = plt.figure()
ax1 = f.add_subplot(111)
ax1.pcolormesh(kx_mat,ky_mat,np.flipud(prob_mat),cmap='viridis',shading='gouraud')
ax1.set_xlabel(r'$k_x a$')
ax1.set_ylabel(r'$k_y a$')
ax1.set_aspect(1.0)

# Polygon patch for Brillouin zone
if (args.polygon == True):
    kappa = 4.0*np.pi/(3.0*np.sqrt(3.0)) # Coordinate of K point
    poly_list = np.array([[kappa, 0], [kappa*0.5, kappa*0.5*np.sqrt(3.0)],
                        [-kappa*0.5, kappa*0.5*np.sqrt(3.0)],
                        [-kappa, 0],
                        [-kappa*0.5, -kappa*0.5*np.sqrt(3.0)],
                        [kappa*0.5, -kappa*0.5*np.sqrt(3.0)]] )
    poly = Polygon(poly_list,ec='w',fc='none',lw=1.0)
    ax1.add_patch(poly)

# Save figure
f.tight_layout()
plt.savefig('GrapheneSine_probability.png',dpi=300)
