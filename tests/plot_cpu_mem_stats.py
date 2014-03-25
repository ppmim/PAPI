#!/usr/bin/env python
"""
Module to do some plots of cpu and mem stats obteined with:

$> sar -P ALL 1 >> cpus_test1_all.txt 
$> sar -r 1 >> mem_test1_all.txt 
 
"""

import fileinput
import re
import numpy as np
import sys

import matplotlib
#matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import pylab


def do_plot(cpu_file, mem_file):

    #Read CPUs values
    lines_cpu = [line.replace("\n","") for line in fileinput.input(cpu_file)]
    prog = re.compile(".*([AP]M.*[all|\d])(\d+)\.(\d+).*")
    matrix = {'all':[],'0':[],'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],'7':[]}
    for line in lines_cpu:
        result = prog.match(line)
        if result:
            cpu_date = result.group(0).split()[0]
            cpu_number = result.group(0).split()[2]
            cpu_value = result.group(0).split()[3]
            matrix[cpu_number].append(cpu_value)
    
    # Read Memory values
    lines_mem = [line.replace("\n","") for line in fileinput.input(mem_file)]
    prog = re.compile(".*[AP]M(.*)(\d+)\.(\d+).*")
    mem = []
    for line in lines_mem:
        result = prog.match(line)
        if result:
            mem_date = result.group(0).split()[0]
            mem_value = result.group(0).split()[4]
            mem.append(mem_value)
    
    #Plot cpus
    #n = range(0,len(matrix['all']))
    my_len = 1000
    n = range(0,my_len)

    
    plt.subplot(10,1,1)
    plt.title("CPUs workload " )
    plt.xlabel("Time")
    plt.ylabel("CPU_all %")
    plt.ylim([0,100])

    plt.grid(True)
    #plt.plot(n, matrix['all'][:my_len],'bo')
    plt.fill_between(n, matrix['all'][:my_len],color='blue')
    
    for i in range(2,10):    
        plt.subplot(10,1,i)
        plt.grid(True)
        plt.ylim([0,110])
        plt.xlabel("Time (s)")
        plt.ylabel("CPU_%d %%"%(i-2))
        plt.plot(n, matrix["%d"%(i-2)][:my_len],'-')

    #Plot memory usage
    plt.subplot(10,1,10)
    plt.xlabel("Time")
    plt.ylabel("Memory %")
    plt.ylim([0,110])
    plt.xlim([0,my_len])
    plt.fill_between(n, mem[:my_len], color='green', )

    
    
    """
    #Full plot
    plt.plot(n, matrix['all'][:my_len], '*', n, matrix['1'][:my_len], '-',
             n, matrix['2'][:my_len], '-',n, matrix['3'][:my_len], '-',
             n, matrix['4'][:my_len], '-',n, matrix['5'][:my_len], '-',
             n, matrix['6'][:my_len], '-',n, matrix['7'][:my_len], '-')
    """
    plt.savefig("plot.eps")
    plt.show()        
    
    
################################################################################
# main
################################################################################
if __name__ == "__main__":
    
    do_plot("/home/panic/DEVELOP/PIPELINE/PANIC/trunk/tests/cpus_test1_all.txt",
            "/home/panic/DEVELOP/PIPELINE/PANIC/trunk/tests/mem_test1_all.txt")

    sys.exit(0)

    
