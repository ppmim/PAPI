#! /usr/bin/env python

""" Get centers of stars regions. """

# Author: Jose M. Ibanez (c) 2014
# Email: jmiguel@iaa.es
# License: GNU GPLv3

import os
import os.path
import numpy as np
import shutil
import sys


def get_stars_regions_centers(path):
    """
    Computer the center of the set of stars of each of the 16 regions
    on the 4kx4k PANIC field.
    """
    
    coords_matrix = np.loadtxt(path)
    regions = np.zeros([16,3]) # region_i: xi, yi, counter_i
    
    # init struct to save the points
    p_regions = []
    for i in range(16):
        p_regions.append([])
        
    for x,y in coords_matrix:
        if x<1024 and y<1024:
            regions[0][0]+=x
            regions[0][1]+=y
            regions[0][2]+=1
            p_regions[0].append((x,y))
        elif x>=1024 and x<2048 and y<1024:
            regions[1][0]+=x
            regions[1][1]+=y
            regions[1][2]+=1
            p_regions[1].append((x,y))
        elif x>=2048 and x<3072 and y<1024:
            regions[2][0]+=x
            regions[2][1]+=y
            regions[2][2]+=1
            p_regions[2].append((x,y))
        elif x>=3072 and x<4096 and y<1024:
            regions[3][0]+=x
            regions[3][1]+=y
            regions[3][2]+=1
            p_regions[3].append((x,y))
        ##
        elif x<1024 and y>=1024 and y<2048:
            regions[4][0]+=x
            regions[4][1]+=y
            regions[4][2]+=1
            p_regions[4].append((x,y))
        elif x>=1024 and x<2048 and y>=1024 and y<2048:
            regions[5][0]+=x
            regions[5][1]+=y
            regions[5][2]+=1
            p_regions[5].append((x,y))
        elif x>=2048 and x<3072 and y>=1024 and y<2048:
            regions[6][0]+=x
            regions[6][1]+=y
            regions[6][2]+=1
            p_regions[6].append((x,y))
        elif x>=3072 and x<4096 and y>=1024 and y<2048:
            regions[7][0]+=x
            regions[7][1]+=y
            regions[7][2]+=1
            p_regions[7].append((x,y))
        ##
        elif x<1024 and y>=2048 and y<3072:
            regions[8][0]+=x
            regions[8][1]+=y
            regions[8][2]+=1
            p_regions[8].append((x,y))
        elif x>=1024 and x<2048 and y>=2048 and y<3072:
            regions[9][0]+=x
            regions[9][1]+=y
            regions[9][2]+=1
            p_regions[9].append((x,y))
        elif x>=2048 and x<3072 and y>=2048 and y<3072:
            regions[10][0]+=x
            regions[10][1]+=y
            regions[10][2]+=1
            p_regions[10].append((x,y))
        elif x>=3072 and x<4096 and y>=2048 and y<3072:
            regions[11][0]+=x
            regions[11][1]+=y
            regions[11][2]+=1
            p_regions[11].append((x,y))
        ##
        elif x<1024 and y>=3072 and y<4096:
            regions[12][0]+=x
            regions[12][1]+=y
            regions[12][2]+=1
            p_regions[12].append((x,y))
        elif x>=1024 and x<2048 and y>=3072 and y<4096:
            regions[13][0]+=x
            regions[13][1]+=y
            regions[13][2]+=1
            p_regions[13].append((x,y))
        elif x>=2048 and x<3072 and y>=3072 and y<4096:
            regions[14][0]+=x
            regions[14][1]+=y
            regions[14][2]+=1
            p_regions[14].append((x,y))
        elif x>=3072 and x<4096 and y>=3072 and y<4096:
            regions[15][0]+=x
            regions[15][1]+=y
            regions[15][2]+=1
            p_regions[15].append((x,y))

    #print p_regions
    reg2sg = {  'SG1-1': 4,'SG1-2':8,'SG1-3':7, 'SG1-4':3 ,
                'SG2-1': 12,'SG2-2':16,'SG2-3':15, 'SG2-4':11,
                'SG3-1': 10,'SG3-2':14,'SG3-3':13, 'SG3-4':9,
                'SG4-1': 2,'SG4-2':6,'SG4-3':5, 'SG4-4':1
            }
    
    # Compute the mean of each region
    for r in reg2sg:
        i = reg2sg[r] - 1
        if regions[i][2]>0:
            x = regions[i][0]/regions[i][2]
            y = regions[i][1]/regions[i][2]
            print "Center of region [%s] = (%f, %f)"%(r,x,y)
            #print "Points of region :",p_regions[i]
            # Write points into a file
            filename = "region_%s.reg"%r
            f = open(filename, "w")
            for p in p_regions[i]:
                f.write("%s %s\n"%(p[0],p[1]))
            f.close()
        else:
            print "Center of region [%s] = (%f, %f)"%(r, np.nan, np.nan)
    
        
    return
    
    """        
    for r in range(len(regions)):
        if regions[r][2]>0:
            x = regions[r][0]/regions[r][2]
            y = regions[r][1]/regions[r][2]
            print "Center of region [%d] = (%f, %f)"%(r,x,y)
        else:
            print "Center of region [%d] = (%f, %f)"%(r,np.nan,np.nan)
    """
def split_region_file():
    """
    Split a ds9 region file into N(=16) region files.
    """
    
    
if __name__ == "__main__":

    path = sys.argv[1]
    centers = get_stars_regions_centers(path)
    sys.exit()

