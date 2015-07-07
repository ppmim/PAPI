#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2015 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
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


import sys
import os
import shutil
import tempfile


# Logging
from misc.paLog import log
from misc.version import __version__
import misc.config
import montage_wrapper as montage

#
# Interface to Montage wrapper
#

def mosaic(files_to_mosaic, raw_dir, output_dir, tmp_dir,
                   background_match=True, out_mosaic=None):
    
    """
    Build the mosaic of the files provided using
    the mosaic function of Montage wrapper.
    
    Parameters
    ----------
    files_to_mosaic: list
        List of files to be mosaiced
    
    raw_dir: str
        Raw directory with input files to mosaic
    
    output_dir: str
        Out directory where mosaic file is created
        
    tmp_dir: str
        Temporal directroy to be used by montage
    
    out_mosaic: str
        Path to the out mosaic to be created
    
    
    Returns
    -------
    filename: str 
        Path to the mosaic file created.
    
    """
    
    log.debug("Running Montage.mosaic...")
    
    # Copy/link input files to the raw directory for montage.mosaic
    m_raw_dir = tempfile.mkdtemp(prefix='PAPI_mosaic_raw', dir=raw_dir)
    
    # Create output directory for montage.mosaic
    m_out_dir = tempfile.mkdtemp(prefix='PAPI_mosaic_out', dir=output_dir)
    # Because montage fails if output dir exits, we delete it in anycase
    os.rmdir(m_out_dir)
    
    # Create work directory for montage.mosaic
    m_work_dir = tempfile.mkdtemp(prefix='PAPI_mosaic_work', dir=tmp_dir)
    # Because montage fails if work dir exits, we delete it in anycase
    os.rmdir(m_work_dir)
    
    for f in files_to_mosaic:
        input_file = os.path.abspath(f)
        basen = os.path.basename(input_file)
        dest = os.path.join(m_raw_dir, basen)
        os.symlink(input_file, dest)
    
    # Call to montage wrapper
    montage.mosaic(input_dir=m_raw_dir,
                   output_dir=m_out_dir,
                   background_match=background_match, 
                   level_only=False,
                   work_dir=m_work_dir,
                   combine='mean')
    
    
    shutil.rmtree(m_raw_dir, ignore_errors=True)
    shutil.rmtree(m_work_dir, ignore_errors=True)
    
    
    # If out_mosaic path was provided, rename the mosaic
    # and remove the out directory.
    if out_mosaic != None:
        shutil.move(os.path.join(m_out_dir, "mosaic.fits"),
                    out_mosaic)
        shutil.rmtree(m_out_dir)
        return out_mosaic
    else:
        return os.path.join(m_out_dir, "mosaic.fits")