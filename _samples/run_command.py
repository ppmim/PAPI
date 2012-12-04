#/usr/bin/env python
"""This code is not currently used in PAPI, and popen2 is deprecated"""


import popen2  # deprecated module, we should use subprocess.Popen
import subprocess
import os
import select

def run_command(cmd, timeout=None, callback=None):
    """
    @summary: Run a command and return the text written to stdout and stderr, plus
    the return value.

    @return:  (int  value, string out, string err)
    """
    child = popen2.Popen3(cmd, True)
    
    #p = subprocess.Popen(cmd, shell=True, bufsize=0,
    #      stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
    #      stderr=subprocess.STDOUT, close_fds=True)
    
    #(fout, fin, ferr) = (p.stdin, p.stdout, p.stderr)
    
    (fout, fin, ferr) = (child.fromchild, child.tochild, child.childerr)
    fin.close()
    stdout = fout.fileno()
    stderr = ferr.fileno()
    outbl = []
    errbl = []
    outeof = erreof = False
    block = 1024
    while True:
        s=[]
        if not outeof:
            s.append(stdout)
        if not erreof:
            s.append(stderr)
        if not len(s):
            break
        (ready, nil1, nil2) = select.select(s, [], [], timeout)
        if stdout in ready:
            outchunk = os.read(stdout, block)
            if len(outchunk) == 0:
                outeof = True
            outbl.append(outchunk)
        if stderr in ready:
            errchunk = os.read(stderr, block)
            if len(errchunk) == 0:
                erreof = True
            errbl.append(errchunk)
        if callback:
            callback()
    fout.close()
    ferr.close()
    w = child.wait()
    if not os.WIFEXITED(w):
        return (-100, out, err)
    rtn = os.WEXITSTATUS(w)
    return (rtn, ''.join(outbl), ''.join(errbl))

