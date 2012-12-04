#!/usr/bin/env python

# script to calculate 'correct' Ra and Dec values taking
# into account the contents of a scamp header.

# module 're' is for regular expressions:
import re
import string
import math
import optparse 
import sys
import os

# check for Python 2.X with X >= 5; the 'optparse' module needs
# Python 2.5 for the used 'epilog' feature (see below):
version = string.split(string.split(sys.version) [0], ".")
# well, Python version 3 just gives us a syntax error at the
# first print statement :-)
if map(int, version) >= [3, 0, 0] or map(int, version) < [2, 5, 0]:
    sys.stderr.write("This script needs Python 2.Y.X (with Y >= 5)\n\n")
    sys.stderr.write("You have Python V%s.%s.%s\n" \
                      % (version[0], version[1], version[2]))
    sys.exit(1)
    
# import non-standard modules:
try:
    import pyfits, numpy
except ImportError:
    sys.stderr.write("This script needs the modules 'pyfits' and 'numpy'!\n")
    sys.stderr.write("see http://www.stsci.edu/resources/software_hardware/pyfits\n")
    sys.stderr.write("and/or http://numpy.scipy.org/\n")
    sys.exit(1)

# The following class allows us to define a usage message epilogue which does
# not skip newline characters:
class MyParser(optparse.OptionParser):
    def format_epilog(self, formatter):
        return self.epilog

usage = "%prog [Options] - apply scamp astrometry to objects"
parser = MyParser(usage=usage, epilog=
"""
Description:
Emmanuel Bertins 'scamp' program allows us to obtain astrometric
calibration of optical imaging data. It provides the resulting astrometric
information as ASCII header keywords which can be used by 'SExtractor'
or 'swarp'. Unfortunately, 'scamp' does not give us directly the
'correct' astrometric Ra and Dec positions of the sources it used
for its calibration. This script performs this task.
It calculates correct Ra, Dec values from a SExtractor catalogue
and an astrometric scamp header. Typically, the provided catalogue is
the same that 'scamp' used for its calibrations. The script takes
the image pixel positions from sources (by default the SExtractor
parameters X_IMAGE and Y_IMAGE) and outputs a new catalogue where
the objects binary LDAC table (default LDAC_OBJECTS) contains two new
entries for the 'correct' Ra and Dec object positions. If a file with
the provided output name already exists, the program overwrites the
old file.

Example:
- aplastrom_scamp.py -i 956233p_4C.cat -s 956233p_4.head -o 956233p_4C_corr.cat

- aplastrom_scamp.py -i 956233p_4C_ldac.cat -s 956233p_4.head
                     -t OBJECTS -p Xpos Ypos -o 956233p_4C_corr.cat

Known Bugs/Shortcomings:
The script only deals with astrometric solutions up to a polynomial degree
of five!

Authors:
  Patrick Kelly   (pkelly3@stanford.edu)
  Thomas Erben    (terben@astro.uni-bonn.de)
""")

parser.add_option("-i", "--input", dest="input",
                  help="name of input LDAC_FITS catalogue (default: %default)",
                  default="input.cat")
parser.add_option("-s", "--scamphead", dest="scamphead",
                  help="name of ASCII scamp header (default: %default)",
                  default="scamp.head")
parser.add_option("-o", "--output", dest="output",
                  help="name of output LDAC_FITS catalogue (default: %default)",
                  default="output.cat")
parser.add_option("-t", "--table", dest="table",
                  help="name of OBJECTS table (default: %default)",
                  default="LDAC_OBJECTS")
parser.add_option("-p", "--position-keys", dest="poskeys", nargs=2,
                  help="names of pixel position keys (default: %default)",
                  default=('X_IMAGE', 'Y_IMAGE'))
parser.add_option("-r", "--result-position-keys", dest="resultposkeys",
                  nargs=2,
                  help="names of result position keys (default: %default)",
                  default=('ALPHA_J2000_corr', 'DELTA_J2000_corr'))

(options, args) = parser.parse_args()

try:
    hdu = pyfits.open(options.input)
except:
    sys.stderr.write("Could not read input catalogue '%s'\n" % (options.input))
    parser.print_help()
    sys.exit(1)

# see how many relevant OBJECT tables we have in the input catalogue
# For later iteration we need the table indices within the catalogue:
tableindices = []
for i in range(len(hdu)):
    if hdu[i].name == options.table:
        tableindices.append(i)

print "catalogue %s has %d %s extensions" \
      % (options.input, len(tableindices), options.table)

try:
    hf = open(options.scamphead, 'r').readlines()
except IOError:
    sys.stderr.write("Could not read scamp header '%s'\n" % (options.scamphead))
    sys.exit(1)

# How many extensions do we have in the scamp header file? We get
# them by counting the 'END' lines. We need the indices of the END
# statements within the file so that we can easily iterate over individual
# extensions below.
# We precede the index list for the position of the END lines with
# a zero; see the 'for line ...' loop below
endindices = [ 0 ] + [ i for i, v in enumerate(hf) if v.strip() == 'END' ]

# sanity check for consistency of extensions in catalogue and scamp file:
if len(tableindices) != len(endindices) - 1:
    sys.stderr.write("The scamp catalogue has %d extensions\n" \
                    % (len(endindices) - 1))
    sys.stderr.write("This is inconsistent with the %d extensions of %s\n" \
                    %(len(tableindices), options.input))       
    sys.exit(1)
    
# now go through the OBJECT extensions:
for i in range(len(tableindices)):
    print "Calculating Ra and Dec for Object Extension %d" % (i + 1)
    # We save the headerkeywords within the scamp file in a dictionary;
    # the keys are the FITS keys and the dictionary values are the corresponding
    # header key values; we only need astrometry relevant quantities:
    headerkeys = {}
    for line in hf[endindices[i]:endindices[i + 1]]:
        if string.find(line, '=') != -1:
            res = re.split('=', line)
            name = res[0].strip()
            res = re.split('/', res[1])
            value = res[0].strip()
            if string.find(name, 'CD') != -1 or \
               string.find(name, 'PV') != -1 or \
               string.find(name, 'CR') != -1:
                headerkeys[name] = float(value)

    table = hdu[tableindices[i]].data
    try:
        x0 = table.field(options.poskeys[0]) - headerkeys['CRPIX1']
        y0 = table.field(options.poskeys[1]) - headerkeys['CRPIX2']
    except KeyError:
        sys.stderr.write("Keys %s and/or %s not present in object ext. %d\n" % \
                         (options.poskeys[0], options.poskeys[1], i + 1))
        sys.exit(1)

    x = x0 * headerkeys['CD1_1'] + y0 * headerkeys['CD1_2']
    y = x0 * headerkeys['CD2_1'] + y0 * headerkeys['CD2_2']
    r = (x**2. + y**2.)**0.5

    # the following lines evaluate polynomials of the form:
    # xi = PV1_0 + PV1_1 * x + PV1_2 * y + ...
    # VERY elegant solution with dictionaries, filter and reduce
    xi_terms = {'PV1_0' : numpy.ones(len(x)), \
                'PV1_1' : x, 'PV1_2' : y, 'PV1_3' : r, 'PV1_4' : x**2., \
                'PV1_5' : x * y, 'PV1_6' : y**2., 'PV1_7' : x**3., \
                'PV1_8' : x**2. * y, 'PV1_9' : x * y**2., 'PV1_10' : y**3., \
                'PV1_11' : r**3, 'PV1_12' : x**4, 'PV1_13' : x**3 * y, \
                'PV1_14' : x**2 * y**2, 'PV1_15' : x * y**3, 'PV1_16' : y**4, \
                'PV1_17' : x**5, 'PV1_18' : x**4 * y, 'PV1_19' : x**3 * y**2, \
                'PV1_20' : x**2 * y**3, 'PV1_21' : x * y**4, 'PV1_22' : y**5, \
                'PV1_23' : r**5}
    
    pv1_keys = filter(lambda x: string.find(x, 'PV1') != -1, headerkeys.keys())
    xi = reduce(lambda x, y : x + y, \
                [xi_terms[k] * headerkeys[k] for k in pv1_keys])

    # dito for eta:
    eta_terms = {'PV2_0' : numpy.ones(len(x)), \
                 'PV2_1' : y, 'PV2_2' : x, 'PV2_3' : r, 'PV2_4' : y**2., \
                 'PV2_5' : y * x, 'PV2_6' : x**2., 'PV2_7' : y**3., \
                 'PV2_8' : y**2. * x, 'PV2_9' : y * x**2., 'PV2_10' : x**3.,
                 'PV2_11' : r**3, 'PV2_12' : y**4, 'PV2_13' : y**3 * x, \
                 'PV2_14' : y**2 * x**2, 'PV2_15' : y * x**3, 'PV2_16' : x**4, \
                 'PV2_17' : y**5, 'PV2_18' : y**4 * x, 'PV2_19' : y**3 * x**2, \
                 'PV2_20' : y**2 * x**3, 'PV2_21' : y * x**4, 'PV2_22' : x**5, \
                 'PV2_23' : r**5}
    
    pv2_keys = filter(lambda x: string.find(x, 'PV2') != -1, headerkeys.keys())
    eta = reduce(lambda x, y : x + y, \
                 [eta_terms[k] * headerkeys[k] for k in pv2_keys])

    # create empty Numpy result arrays for calculated Ra and Dec values:
    RaCorr = numpy.empty(len(x0))
    DecCorr = numpy.empty(len(y0))
    
    CRVAL1 = headerkeys['CRVAL1'] / 180.0 * math.pi
    CRVAL2 = headerkeys['CRVAL2'] / 180.0 * math.pi

    # now calculate Ra and Dec according to the scamp
    # information:
    for k in range(len(xi)):
        XI = xi[k] / 180.0 * math.pi
        ETA = eta[k] / 180.0 * math.pi
        p = math.sqrt(XI**2. + ETA**2.) 
        c = math.atan(p)
        a = CRVAL1 + math.atan((XI * math.sin(c)) / \
                               (p * math.cos(CRVAL2) * math.cos(c) - \
                                ETA * math.sin(CRVAL2) * math.sin(c)))
        d = math.asin(math.cos(c) * math.sin(CRVAL2) + \
                      ETA * math.sin(c) * math.cos(CRVAL2) / p)
        
        RaCorr[k] = a * 180.0 / math.pi
        DecCorr[k] = d * 180.0 / math.pi

    # create new table columns and append them to the old
    # objects table:
    RaCorrColumn = pyfits.Column(name = options.resultposkeys[0],
                                 format = '1D', unit = 'deg',
                                 array = RaCorr)
    
    DecCorrColumn = pyfits.Column(name = options.resultposkeys[1],
                                  format = '1D', unit = 'deg',
                                  array = DecCorr)

    newObjectColumns = hdu[tableindices[i]].columns + \
                       RaCorrColumn + DecCorrColumn

    newObjectTableHDU =  pyfits.new_table(newObjectColumns,
                                          hdu[tableindices[i]].header)

    # just replace the old table by the new one:
    hdu[tableindices[i]] = newObjectTableHDU

# write result catalogue:
if os.path.isfile(options.output):
    print "WARNING: File %s already exists. I delete it!" % (options.output)
    os.remove(options.output)
    
print "Writing output catalogue %s" % (options.output)
try:    
    hdu.writeto(options.output)
except:
    sys.stderr.write("Could not write to '%s' (file exists?)\n" \
                     % (options.output))
    sys.exit(1)
    
# and bye:
sys.exit(0)

