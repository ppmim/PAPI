#! /usr/bin/python

import os
import sys
import fileinput
import glob

__author__="panic"
__date__ ="$Apr 27, 2009 06:00:23 PM$"


def linkSourceFiles( source, dest ):
    """Create a symbolic link to all the sources specified in 'source', which can be a file a dir"""
  
    if (type(source)==type(list())):
        #we suppose is a list
        for file in source:
            print "FILE=", file
            os.symlink(file, dest+"/"+os.path.basename(file))
    elif (os.path.isfile(source)):
        #We have a source-file with absolute path for the data files
        file=open(source, 'r')
        for line in file:
            #print "FILE=", line
            if line.endswith('\n'):
                line=line.replace('\n','')
            print 
            #print "DEST=", dest+"/"+os.path.basename(line)
            os.symlink(line, dest+"/"+os.path.basename(line))
    elif (os.path.isdir(source)):
        #We have a source-directory as input data
        for file in glob.glob(source+"*.fits"):
            print "FILE=", file
            os.symlink(file, dest+"/"+os.path.basename(file))
    else:
        print "Error reading source; cannot indentify the type of input"
#################################################################

if __name__ == "__main__":
    
    if len(sys.argv)>2:
        source = sys.argv[1]
        dest = sys.argv[2]
        linkSourceFiles(source, dest)
    else:
        print "Wrong number of arguments (=2)"
        sys.exit(0)      