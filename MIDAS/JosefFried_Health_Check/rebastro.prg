!
!	Calls the programm rebastro.f to rebin images and perform an astrometric correction.
!	The center in delta is needed for computation purposes.
!	Be aware, that the step in RA has to be given a minus sign.
!
!	JWF 03.01
!



DEFI/PAR P1 ? C "input file a =?"


IF p1 .eq. "help" THEN
      write/out @@ rebastro input output step
      return
ENDIF


DEFI/PAR P2 ? C "output file=?"
!DEFI/PAR P3 -1.22e-4,1.22e-4 N "step_ra,step_dec in degrees/pixel = ?"
DEFI/PAR P3 -10e-5,10e-5 N "step_ra,step_dec in degrees/pixel = ?"






WRITE/KEYW IN_A/C/1/20 'P1'
WRITE/KEYW OUT_A/C/1/20 'P2'
write/key inputd/d/3/2 {p3}







RUN /disk2/fried/midas/prog/rebastro.exe

cop/dd {p1} lhcuts/r/1/2 {p2} lhcuts/r/1/2
copy/dd {p1} *,3,history,lhcuts {p2} 
write/descr {p2} ident "{p1} rebinned"

write/out
write/out output -> {p2}
write/out


exit:
