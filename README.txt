This is PAPI, the PANIC data reduction PIpeline. 

PAPI is distributed under GNU GPL, either version 3 of the License, 
or (at your option) any later version. See the file COPYING for details.

PAPI requires the following packages installed in order to
be able to be installed and work properly:

 - setuptools (http://peak.telecommunity.com/DevCenter/setuptools)
 - numpy (http://numpy.scipy.org/) 
 - scipy (http://www.scipy.org)
 - pyfits (http://www.stsci.edu/resources/software_hardware/pyfits)
 - pyraf (http://www.stsci.edu/institute/software_hardware/pyraf)
 - pywcs (http://stsdas.stsci.edu/astrolib/pywcs/)

PANIC is a general purpose Panoramic Near Infrared camera for Calar Alto. 
It is optimized for use at the 2.2m telescope, but can also be installed 
at the 3.5m telescope. It will work in the nIR bands Z, J, H and K. 

Installing
**********

1. Unzip and untar the papi files in a suitable location.

    * tar zxf papi_X.X.X.tgz
    * ln -s papi_X.X.X papi
    * cd papi


2. Build the source

        * cd papi/irdr/src
        * make all


3. Edit papi_setup.sh script

    * Modify the PAPI_HOME and PAPI_PROD in the papi_setup.[sh] file in the papi directory
    * Add the papi_setup.[sh] to your .bash_profile 

        * Bash: . $PAPI_HOME/papi_setup.sh
        * CSH: source $PAPI_HOME/papi_setup.csh


Documentation
*************

See the docs/ directory for full documentation.




Webpage: https://www.iaa.es/PANIC
Maintainer: jmiguel@iaa.es

     
