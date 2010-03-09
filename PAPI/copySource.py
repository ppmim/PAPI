#! /usr/bin/python

import os
import sys
import fileinput
import glob

__author__="panic"
__date__ ="$Apr 27, 2009 06:00:23 PM$"


def linkSourceFiles( source ):
"""Create a symbolic link to all the sources specified in 'source', which can be a file a dir"""
  
    if (os.path.isfile(source)):
        #We have a source-file with absolute path for the data files
        file=open(source, 'r')
        for line in file:
            print "FILE=", line
            if line.endswith('\n'):
                line=line.replace('\n','')
            os.symlink(line, os.path.basename(line))
    else:
        #We have a source-directory as input data
        for file in glob.glob(source+"*.fits"):
            print "FILE=", file
            os.symlink(file, os.path.basename(file))

#################################################################

if __name__ == "__main__":
    
    if len(sys.argv)>1:
        source = sys.argv[1]
        linkSourceFiles(source)
    else:
        print "Wrong number of arguments"
        sys.exit(0)      