#! /usr/bin/env python

# Copyright (c) 2008-2015 IAA-CSIC  - All rights reserved. 
# Author: Jose M. Ibanez. 
# Instituto de Astrofisica de Andalucia, IAA-CSIC
#
# This file is part of PAPI (PANIC Pipeline)
#
# This program is free software: you can redistribute it and/or modify
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

from __future__ import print_function
import os
import time

major = 1
minor = 2

rlfile = 'version.py'
backup = 'version.py.bak'
path = os.path.dirname(os.path.abspath(__file__))

def make_release():

    try:
        release = time.strftime("%Y%m%d%H%M%S", time.gmtime(time.time()))
        release = int(release)

        cwd = os.getcwd()
        os.chdir(path)

        if os.path.exists(backup):
            os.remove(backup)

        if os.path.exists(rlfile):
            os.rename(rlfile, backup)

        with open(rlfile, 'w') as out_f:
            out_f.write("# this file was automatically generated\n")
            out_f.write("major = %d\n" % major)
            out_f.write("minor = %d\n" % minor)
            out_f.write("release = %d\n" % release)
            out_f.write("\n")
            out_f.write("__version__ = '%d.%d.%d' % (major, minor, release)\n")
            out_f.write("\n")
    finally:
        os.chdir(cwd)
        return "%d.%d.%d" % (major, minor, release)
    
if __name__ == "__main__":
    print(make_release())
    
