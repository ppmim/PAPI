#! /usr/bin/env python


import numpy as np
import pyraf.iraf as iraf
import argparse
import sys

def read_psfmeasure_log(logfile):
        """
        Read read_psfmeasure log file and convert to x,y,fwhm simple files.
        """

        xout = np.array([])
        yout = np.array([])
        FWHM = np.array([])
        
        n_images = 0
        new_log = None
        for line in open(logfile, 'r'):
            # First line is a description of the file, second is a newline \n
            # and third the description of the columns.
            field_pattern = "/data1/PANIC"
            if line != "\n" and field_pattern in line.split()[0]:
                print "Found new Image"
                if new_log: new_log.close()
                new_log = open(logfile +".new.%0d"%n_images,"w")
                n_images += 1
                xout = np.append(xout, float(line.split()[1]))
                yout = np.append(yout, float(line.split()[2]))
                FWHM = np.append(FWHM, float(line.split()[4]))
                new_log.write("%f  %f  %f\n"%(xout[-1], yout[-1], FWHM[-1]))
            else:
                try:
                    xout = np.append(xout, float(line.split()[0])) 
                    yout = np.append(yout, float(line.split()[1]))
                    FWHM = np.append(FWHM, float(line.split()[3]))
                    new_log.write("%f  %f  %f\n"%(xout[-1], yout[-1], FWHM[-1]))
                except:
                    # First line
                    pass
        
        
	new_log.close()
        
        n_points = len(xout)/n_images
        print "N_Points",n_points
        print "N_images",n_images
        return
        for i in range(n_points):
            for im in range(n_images):
                print "FWHM_%d_%d = %f"%(im, i, FWHM[im*n_points+i])

if __name__ == "__main__":

    data_file = sys.argv[1]
    
    read_psfmeasure_log(data_file)



