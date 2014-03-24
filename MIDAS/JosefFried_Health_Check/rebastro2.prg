!
!	Calls the programm rebastro.f to rebin images and perform an astrometric correction.
!	The center in delta is needed for computation purposes.
!	Be aware, that the step in RA has to be given a minus sign.
!
!	JWF 03.01
!



DEFI/PAR P1 ? C "input file a =?"


IF p1 .eq. "help" THEN
      write/out @@ rebastro2 input output start(ra),start(dec) npixx,npixy stepx,stepy
      return
ENDIF


DEFI/PAR P2 ? C "output file=?"
defi/par p3 ? n "start = start(ra),start(dec)? "
defi/par p4 2000,2000 n "npixx,npixy = ? "
DEFI/PAR P5 -6.2222e-5,6.2222e-5 N "stepx,y [degrees/pixel] = ?"


WRITE/KEYW IN_A/C/1/20 'P1'
WRITE/KEYW OUT_A/C/1/20 'P2'
WRITE/KEYW INPUTI/I/1/2 'P4'
WRITE/KEYW INPUTD/D/1/2 'p3'
write/keyw INPUTD/D/3/2 'p5'



RUN FMP:rebastro2.exe

cop/dd {p1} lhcuts/r/1/2 {p2} lhcuts/r/1/2
copy/dd {p1} *,3,history,lhcuts {p2} 
write/descr {p2} ident "{p1} rebinned"

