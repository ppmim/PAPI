#!/usr/bin/env python

# Copyright (c) 2009-2012 IAA-CSIC  - All rights reserved. 
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

pattern = "/data2/out/Q%02d/NorthA_H_%04d.Q%02d_D_F.skysub.ast.fits"
#pattern = "/data2/out/Q%02d/NorthA_J_%04d_coadd.Q%02d_D_F.skysub.ast.fits"
#pattern = "/data2/out/Q%02d/NorthA_Sky_J_%04d_coadd.Q%02d_F.skysub.ast.fits"

import misc.mef

def createMEFs(input_pattern, start, end):
    """
    Create a set of MEF from a given set of individual frames Q1,Q2,Q3,Q4
    It is done for a test, to get a SCAMP distorion plot with all detectors.
    """
    for i in range(start, end + 1):
        frame_list = []
        for det in range(1, 5):
            frame_list.append(input_pattern %(det, i, det))
        
        mef = misc.mef.MEF(frame_list)
        output_file = "mef_%04d.fits" % i
        mef.createMEF(output_file, out_dir="/data2/out/")
        
createMEFs(pattern, 1, 36)
