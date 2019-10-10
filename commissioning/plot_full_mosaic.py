#! /usr/bin/env python


""" Scatter-plot the data points (x, y, FWHM) of the four detectors
and fit a regression plane."""

from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import matplotlib.pyplot as plt
import numpy as np
import pylab
import matplotlib
import sys


if __name__ == "__main__":

    data_file = sys.argv[1]
    
    ax = pylab.subplot(1, 1, 1)

    all_points = []
    path = data_file
    points = [list(x) for x in np.loadtxt(path)]

    all_points += points
    x, y, fwhm = zip(*points)
    mymap = plt.cm.get_cmap('jet')
    
    #print "FWHM[%d]=%s"%(nplate, fwhm)
    # Normalize FWHM
    fwhm /= np.max(np.abs(fwhm),axis=0)
    my_colors = mymap(fwhm)
    norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1, clip = False)
    ax.scatter(x, y, s=fwhm*100,  color=my_colors, cmap=mymap, alpha=1.0, edgecolors='None')
    ax.set_title("All plates")
    #print "FHWM_n=%s"%fwhm


    
    pylab.tight_layout()
    pylab.legend(numpoints=1)
    pylab.show()

