#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# makeMaskedSkyList.py
#
# Created    : 25/03/2010    jmiguel@iaa.es
#
# TODO
################################################################################

################################################################################
# Import necessary modules

import getopt
import sys
import os
import logging
import fileinput
import time
from optparse import OptionParser

        
################################################################################
# Functions       
################################################################################
class MaskedSkyList:
    def __init__(self, files, offsets, type, output_file):
        
        self.source_file_list = files
        self.offsets_file = offsets
        self.type = type
        self.output_file = output_file
        
    def create (self):
        
        fileList=[line.replace( "\n", "") for line in fileinput.input(self..source_file_list)]
        offsetsList=[line.replace( "\n", "") for line in fileinput.input(self.offsets_file)]
        outF=open(self.output_file, "w+")
    
        j=0
        for i in range(0,len(fileList)):
            if self.type=='onoff':
                if not i%2: # even(par)
                    outF.write(fileList[i] + " " + offsetsList[j] + "\n")
                else: 
                    outF.write(fileList[i] + " " + offsetsList[j] + "\n")
                    j=j+1
            else: #offon
                if not i%2: # even(par)
                    outF.write(fileList[i] + " " + offsetsList[j] + "\n")
                else: 
                    outF.write(fileList[i] + " " + offsetsList[j] + "\n")
                    j=j+1
                
        outF.close()      
      
# main
if __name__ == "__main__":
    print 'Start makeMaskedSkyList....'
    # Get and check command-line options
    args = sys.argv[1:]
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
                  
    parser.add_option("-f", "--files",
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It has to be a file")
    
    parser.add_option("-o", "--offsets",
                  action="store", dest="offsets_file", help="Offsets file of object/target frames")
    
    parser.add_option("-t", "--type",
                  action="store", dest="type", default="onoff", help="type of dithering (onoff|offon)")
                  
    parser.add_option("-O", "--out",
                  action="store", dest="output_file",  help="Output file to be generated")
    
    (options, args) = parser.parse_args()
    
    if not options.source_file_list or not options.output_file or not options.offsets_file or not options.type or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )


    # 
    obj = MaskedSkyList(options.source_file, options.offsets_file, options.type, options.output_file)
    obj.create()
    #
    
    print 'End of makeMaskedSkyList....'
        
