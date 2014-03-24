! fittab
!
if "{p1}" .eq. "help" then
   write/out fittab to fit functions to table data
   use as @@ fittab xcol,ycol,weightcol fitfunction parameters mode,nitera eps0,eps new_column_with_fit
   goto exit
endif

DEFI/PAR P1 hist2 C "TABLE"
DEFI/PAR P2 :mR2,:lf,:weight c "xcol,ycol"
DEFI/PAR P3 pol C "FIT FUNCTION"
DEFI/PAR P4 -7.5,.4 N "PARAMETERS"
DEFI/PAR P5 -1,30 N " MODE,NITERA"
DEFI/PAR P6 1.E-4,1.E-4 N "EPS0,EPS"
DEFI/PAR P7 fita C "CREATE NEW COLUMN WITH FIT?"

WK INPUTC/C/21/30 'P3'
WK INPUTI/I/3/2 'P5'
WK INPUTR/R/1/2 'P6'
WK INPUTR/R/11/10 'P4'
WK INPUTC/C/1/20 'P7'


wk fitpar/r/30/1 0.
RUN /disk1/fried/midas/prog/fittab
!write/out fit done!
