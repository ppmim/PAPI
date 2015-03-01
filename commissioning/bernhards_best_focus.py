#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2015 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
#
# PAPI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt
import sys


def robust_mean(x):
    y = x.flatten()
    n = len(y)
    y.sort()
    ind_qt1 = round((n+1)/4.)
    ind_qt3 = round((n+1)*3/4.)
    IQR = y[ind_qt3]- y[ind_qt1]
    lowFense = y[ind_qt1] - 1.5*IQR
    highFense = y[ind_qt3] + 1.5*IQR
    ok = (y>lowFense)*(y<highFense)
    yy=y[ok]
    return yy.mean(dtype='double')

#-------Robust Standard Deviation---

def robust_std(x):
    y = x.flatten()
    n = len(y)
    y.sort()
    ind_qt1 = round((n+1)/4.)
    ind_qt3 = round((n+1)*3/4.)
    IQR = y[ind_qt3]- y[ind_qt1]
    lowFense = y[ind_qt1] - 1.5*IQR
    highFense = y[ind_qt3] + 1.5*IQR
    ok = (y>lowFense)*(y<highFense)
    yy=y[ok]
    return yy.std(dtype='double')

#-------Robust variance---

def robust_var(x):
    y = x.flatten()
    n = len(y)
    y.sort()
    ind_qt1 = round((n+1)/4.)
    ind_qt3 = round((n+1)*3/4.)
    IQR = y[ind_qt3]- y[ind_qt1]
    lowFense = y[ind_qt1] - 1.5*IQR
    highFense = y[ind_qt3] + 1.5*IQR
    ok = (y>lowFense)*(y<highFense)
    yy=y[ok]
    return yy.var(dtype='double')

def readStarfocus(log_file):
    """
    Read the results from the iraf.starfocus log file and compute
    the best focus for that execution. Only non-saturated data
    are read.
    """
    
    # Read the last lines 
    with open(log_file, "r") as f:
        f.seek (0, 2)           # Seek @ EOF
        fsize = f.tell()        # Get Size
        f.seek (max (fsize-2**15, 0), 0) # Set pos @ last chars
        lines = f.readlines()       # Read to end
    
    lines.reverse()
    my_lines = []
    
    # Look for the last execution (It must start with
    # 'NOAO/IRAF') 
    while not lines[0].strip().startswith('NOAO/IRAF'):
        my_lines.append(lines[0])
        lines.pop(0)

    my_lines.reverse()
    data = []
    # start to read the right columns
    for line in my_lines:
        if line.strip().startswith('Image'):
            # Heading line
            continue
        elif line.strip().startswith('Best'):
            # End reading
            break
        elif len(line.split()) == 0: 
            # Blank line
            continue
        else:
            # Read columns, only non-saturated data
            if line.split()[-1] != '*':
                if line.split()[0].startswith("/"):
                    data.append(line.split()[1:9])
                else:
                    data.append(line.split()[0:8])
     
    return data

def getBestFocus(data, output_file):
    """
    Fit the input data read from starfocus log file,
    to a parabola and find out the minimin (best focus estimation).
    
    Parameters
    ----------
    data: list
        A list with N rows x M (8) columns with the next correspondence:
        Column    Line     Mag   Focus   MFWHM Beta   Ellip      PA  
        
    output_file: str
        Filename of out plot with the fitting computed.
    """
    
    import socket
    hostname = socket.gethostname()
    if hostname == 'panic22':
        foclimits = [-1, 27]
    elif hostname == 'panic35':
        foclimits = [10, 60]
    else:
        foclimits = None 
        
    d = np.array(data, dtype=np.float32)
    good_focus_values = d[:, 3] # focus
    fwhm_values = d[:, 4] # PSF-value (MFWHM, GFWHM, FWHM, ...)
    
    
    print "\n---------"
    print "N_POINTS: ", len(fwhm_values)
    print "----------\n"
    
    
    m_foc = good_focus_values.mean()
    good_focus_values = good_focus_values - m_foc
    z = np.polyfit(good_focus_values, fwhm_values, 2)
    print "Fit = %s  \n"%str(z)
    # Note that poly1d returns polynomial’s coefficients, in increasing powers !
    pol = np.poly1d(z)
    
    xp = np.linspace(np.min(good_focus_values ) - 0.5, 
                     np.max(good_focus_values ) + 0.5, 500) # number or points to interpolate

    # best focus is derivative of parabola = 0
    # but check if it is correctly curved
    if pol[2] < 0:
        print "ERROR: Parabola fit unusable!"
    best_focus = - pol[1] / (2. * pol[2])
    min_fwhm = pol([best_focus])
    
    if foclimits and (best_focus + m_foc < foclimits[0] or best_focus + m_foc > foclimits[1]):
        print "ERROR: Best focus out of range!"
    print "BEST_FOCUS = ", best_focus + m_foc
    print "MIN_FWHM = ", min_fwhm
    
    # Plotting
    plt.plot(good_focus_values + m_foc, fwhm_values, '.')
    plt.plot(xp + m_foc, pol(xp), '-')
    plt.axvline(best_focus + m_foc, ls='--', c='r')
    plt.title("Focus serie - Fit: %f X^2 + %f X + %f\n Best Focus=%6.3f mm; FWHM=%6.3f (px)" 
        %(pol[2], pol[1], pol[0], best_focus + m_foc, min_fwhm))
    
    plt.xlabel("T-FOCUS (mm)")
    plt.ylabel("FWHM (pixels)")
    plt.xlim(np.min(good_focus_values + m_foc) - 0.1, np.max(good_focus_values + m_foc) + 0.1)
    plt.ylim(np.min(fwhm_values) - 1, np.max(fwhm_values) + 1 )
    if pol[2] < 0:
    	plt.figtext(0.5, 0.5, 'ERROR: Parabola fit unusable!', size='x-large', color='r', weight='bold', ha='center', va='bottom')
    if foclimits and (best_focus + m_foc < foclimits[0] or best_focus + m_foc > foclimits[1]):
    	plt.figtext(0.5, 0.5, 'ERROR: Best focus out of range!', size='x-large', color='r', weight='bold', ha='center', va='top')
    plt.grid()
    
    
    plt.savefig(output_file)
    
    show = True
    if show:
        plt.show(block=True)
            
    sys.stdout.write("\nGenerated plot file: %s\n"%output_file)
            
    # Print out Values
    #for idx, foc_value in enumerate(good_focus_values):
    #    sys.stdout.write("\nFoc. value: %s   -->  FWHM: %s"%(foc_value, fwhm_values[idx]))
        


# Testing
# -------
data = readStarfocus(sys.argv[1])
getBestFocus(data, "starfocus.pdf")
