!
! corfu.rpg
!
DEFI/PAR P1 f0003p00 C " table_in=?"
rd {p1}.tbl TSELTABL
DEFI/PAR P2 z C "table_out=?"
!            scale,rmin,rmax,dr,fac,seed
DEFI/PAR P3 1,1.,3.6,.4,1,63333 N  
DEFI/PAR P4 21,22  N " x_col,y_col=?"
DEFI/PAR P5 0.,0. N " x,y=?"
DEFI/PAR P6 AUTO C " MODE [AUTO/CROSS] = ? "

WRITE/KEYW INPUTR/R/1/6 {P3}
WRITE/KEYW INPUTI/I/1/10 0 ALL
WRITE/KEYW INPUTI/I/1/4 {P4}
WRITE/KEYW INPUTR/R/7/2 {P5}
WRITE/KEYW INPUTC/C/1/5 {P6}
RUN FMP:corfu
