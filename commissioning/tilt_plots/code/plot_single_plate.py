#! /usr/bin/env python

# Author: Victor Terron (c) 2014
# Email: `echo vt2rron1iaa32s | tr 132 @.e`
# License: GNU GPLv3

""" Scatter-plot the data points (x, y, FWHM) of the four detectors
and fit a regression plane."""

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

# This code is taken from Stack Overflow:
# https://stackoverflow.com/a/20700063/184363
#
# The normal vector of a plane a*x + b*y +c*z = 0, equals (a,b,c)
def plane(x, y, params):
    a = params[0]
    b = params[1]
    c = params[2]
    z = a*x + b*y + c
    return z

def error(params, points):
    result = 0
    for (x,y,z) in points:
        plane_z = plane(x, y, params)
        diff = abs(plane_z - z)
        result += diff**2
    return result

def cross(a, b):
    return [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]

def get_labels():
    with open("./data/headers", "rt") as fd:
        return [x.strip() for x in fd.readline().split("|")]

def get_comparison_coords():
    with open("./data/comparison", "rt") as fd:
        return [float(x.strip()) for x in fd.readline().split()]

if __name__ == "__main__":

    colors = itertools.cycle(['blue', 'red', 'green', 'black'])

    ax = pylab.subplot(1, 1, 1, projection='3d')

    all_points = []
    for nplate in range(1, 5):
        path = './data/plate{0}'.format(nplate)
        points = [list(x) for x in np.loadtxt(path)]
        
        # To take into account the gap (167px)
        print("----Original POINTS {0}".format(points))
        for row in points:
            if row[0]<2048: row[0]-=84
            elif row[0]>=2048: row[0]+=84
            if row[1]<2048: row[1]-=84
            elif row[1]>=2048: row[1]+=84
        print("****New Points = {0}".format(points))
        # end_gap_correction
        
        all_points += points
        x, y, fwhm = zip(*points)
        kwargs = dict(color=next(colors),
                      label="Plate {0}".format(nplate))
        ax.plot(x, y, fwhm, 'o', **kwargs)
        ax.set_title("All plates")
        labels = get_labels()
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_zlabel(labels[2])

    # Now fit a plane to all the points
    x, y, fwhm = zip(*all_points)
    fun = functools.partial(error, points=all_points)
    params0 = [0, 0, 0]
    kwargs = dict(method='Powell')
    res = scipy.optimize.minimize(fun, params0, **kwargs)

    a = res.x[0]
    b = res.x[1]
    c = res.x[2]

    point  = np.array([0.0, 0.0, c])
    normal = np.array(cross([1,0,a], [0,1,b]))
    d = -point.dot(normal)

    print("Plane Equation:")
    print("A = {0}".format(a))
    print("B = {0}".format(b))
    print("C = {0}".format(-1))
    print("D = {0}".format(-d))

    xx, yy = np.meshgrid([min(x), max(x)], [min(y), max(y)])
    z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]
    ax.plot_surface(xx, yy, z, alpha=0.2, color=[0,1,0])

    # Now the comparison plane (Z = 0)
    a, b, c = get_comparison_coords()
    point  = np.array([0.0, 0.0, c])
    normal = np.array(cross([1,0,a], [0,1,b]))
    d = -point.dot(normal)

    print()
    print("Comparison plane equation:")
    print("A = {0}".format(a))
    print("B = {0}".format(b))
    print("C = {0}".format(-1))
    print("D = {0}".format(-d))

    x, y, fwhm = zip(*all_points)
    xx, yy = np.meshgrid([min(x), max(x)], [min(y), max(y)])
    z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]
    ax.plot_surface(xx, yy, z, alpha=0.15, color='blue')

    fig = pylab.gcf()
    msg = "PANIC: z plane determination-first detector focus cycle"
    fig.suptitle(msg, fontsize=14)
    pylab.tight_layout()
    pylab.legend(numpoints=1)
    pylab.show()

