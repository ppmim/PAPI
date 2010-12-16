#!/bin/csh 
#
# change ZPX WCS to ZPN WCS
#

if (! $?IRDR_BASEDIR) then
  echo please setenv IRDR_BASEDIR
  exit
endif

set path = ($path $IRDR_BASEDIR/extern/wcstools/bin)

foreach file ($*)
  sethead $file CTYPE1='RA---ZPN' CTYPE2='DEC--ZPN' PROJP1=1.0 PROJP3=220.0
end

exit
