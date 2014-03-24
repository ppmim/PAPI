DEFI/PAR P1 qso ima "input file=?"
DEFI/PAR P2 wei ima "weight file=?"
DEFI/PAR P3 psf ima "psf=?"
DEFI/PAR P4 ima ima "scale image=?"
DEFI/PAR P5 mod ima "output model=?"
DEFI/PAR P6 con ima "output convol=?"
DEFI/PAR P7 VV c "function=? [PL,GA,DG,LO,HS,VS,ES,PS,HF,VF,EF,VV,EV,GG,PV,SC]"
DEFI/PAR P8 2,.5,7,0,61.1,61.1,1 N " parameters=?"
!DEFI/PAR P8 0 N "IFAIL=?"
!DEFI/PAR P8 .8,.65,45.06,44.95,1.77,1.85,1.79,8.00,7.48 N "b/a,PA,xc,yc=?"
!DEFI/PAR P8 2.3,2.3,3.4,90,154
WRITE/KEYW IN_A/C/1/8 'P1'
WRITE/KEYW IN_B/C/1/8 'P2'
WRITE/KEYW IN_C/C/1/8 'P3'
WRITE/KEYW IN_D/C/1/8 'P4'
WRITE/KEYW OUT_A/C/1/8 'P5'
WRITE/KEYW OUT_B/C/1/8 'P6'
WRITE/KEYW INPUTC/C/1/40 'P7'
WRITE/KEYW INA/D/1/40 'P8'
! fixed parameters
!
WRITE/KEYW INPUTD/D/1/9 61.0,61.0,3,4,5,6,7,8,9
!
! logging , convolution on/of
write/key  action/c/1/2 on
write/key  action/c/3/2 of
!
! ifail
!
WRITE/KEYW INPUTI/I/1/1 0
!
RUN FMP:cnlf
