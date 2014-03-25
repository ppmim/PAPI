!
! surf.prg
!

DEFI/PAR P1 x ima "file=?"
DEFI/PAR P2 one ima "weight file=?"
DEFI/PAR P3 CURSOR C "input [cursor/numbers]?"
DEFI/PAR P4 SAME TBL "output table name=?"
!DEFI/PAR P5 41.787,41.857 N "x,y center"
defi/par P5 1,0,10,10,3,3,.25 N " b/a, P.A., radmax, nbin, niter, kappa, arcsec/pix"

if "{p3}" .eq. "CURSOR" then
	write/out *** mark center ***
	cen/gauss
	cop/kk outputr/r/5/2 inputr/r/1/2
else
	write/key inputr/r/1/2 {p3} 
endif

WRITE/KEYW IN_A/C/1/20 {P1}
WRITE/KEYW IN_B/C/1/20 {P2}
WRITE/KEYW OUT_A/C/1/20 {P4}


WRITE/KEYW INPUTR/R/3/9/ {P5}

RUN FMP:surf

