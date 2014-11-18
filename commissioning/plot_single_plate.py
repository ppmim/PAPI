#! /usr/bin/env python


""" Scatter-plot the data points (x, y, FWHM) of the four detectors
and fit a regression plane."""

from __future__ import division
#from __future__ import print_function
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


    ax = pylab.subplot(1, 1, 1)

    all_points = []
    for nplate in range(1, 5):
        path = './data/plate{0}'.format(nplate)
        points = [list(x) for x in np.loadtxt(path)]

        all_points += points
        x, y, fwhm = zip(*points)
        mymap = plt.cm.get_cmap('jet')
        # Normalize FWHM
        print "FWHM[%d]=%s"%(nplate, fwhm)
        fwhm /= np.max(np.abs(fwhm),axis=0)
        my_colors = mymap(fwhm)
        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1, clip = False)
        ax.scatter(x, y, s=fwhm*100,  color=my_colors, cmap=mymap, alpha=1.0, edgecolors='None')
        #ax.plot(x, y, fwhm, 'o', **kwargs)
        ax.set_title("All plates")
        print "FHWM_n=%s"%fwhm


    
    pylab.tight_layout()
    pylab.legend(numpoints=1)
    pylab.show()

