#! /usr/bin/env python


""" Create a 2x2 figure; for each one of them, scatter-plot the data
points (x, y, FWHM) of the detector and fit a regression plane. """

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from mpl_toolkits.mplot3d import Axes3D
import functools
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pylab
import scipy.optimize
import matplotlib


if __name__ == "__main__":


    def plot_plate_plane_fit(nplate):
        """ Subplot the plane fit of the 'plateN' data file """

        # Position of each plate
        positions = {1:3, 2:4, 3:2, 4:1}

        ax = pylab.subplot(2, 2, positions[nplate])
        points = np.loadtxt('./data/plate{0}'.format(nplate))
        x, y, fwhm = zip(*points)

        
        print("FWHM[{0}]={1}".format(nplate, fwhm))
        # Normalize FWHM
        fwhm /= np.max(np.abs(fwhm),axis=0)
        mymap = plt.cm.get_cmap('jet')
        my_colors = mymap(fwhm)
        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1, clip = False)
        ax.scatter(x, y, s=fwhm*100,  color=my_colors, cmap=mymap, norm=norm, alpha=1.0, edgecolors='None')
        #ax.plot(x, y, fwhm, 'o', **kwargs)
        ax.set_title("All plates")
        ax.set_title("Plate {0}".format(nplate))
        
        print("FHWM_n={0}".format(fwhm))


    for n in range(1, 5):
        plot_plate_plane_fit(n)

    fig = pylab.gcf()
    msg = "PANIC: Tilt determination during telescope commissioning"
    fig.suptitle(msg, fontsize=14)
    pylab.tight_layout()
    plt.show()

