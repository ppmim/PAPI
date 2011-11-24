!
!	Calls the programm drizzle to rebin images and perform an astrometric correction.
!	The center in delta is needed for computation purposes.
!	Be aware, that the step in RA has to be given a minus sign.
!
!	JWF 03.01
!



DEFI/PAR P1 ? C "Input catalogue =?"


IF p1 .eq. "help" THEN
      write/out @@ drizzle catalogue output [start] [step] [npix]
      return
ENDIF


DEFI/PAR P2 ? C "output file=?"
defi/par p3 0,0 n "start_ra,start_dec = ? "
DEFI/PAR P4 0,0 N "step_ra,step_dec in degrees/pixel = ?"
defi/par p5 4100,4100 n "size of outframe in pixels [4100,4100] = ? "


defi/par p6 0,0 N "files from,to=?"


WRITE/KEYW IN_A/C/1/20 'P1'
WRITE/KEYW OUT_A/C/1/20 'P2'
write/key inputd/d/1/2 {p3}
write/key inputd/d/3/2 {p4}

write/key inputi/i/1/2 {p5}
write/key inputi/i/3/2 {p6}




! compute start if not given
!if {p3} .eq. 0,0 then
!   @@ astrocomp {p1} 1,1
!   WRITE/KEYW INPUTD/D/1/2 {Q1},{Q2}
!   write/out start values = {inputd(1)},{inputd(2)}
!endif



! compute step if not given
!if {p4} .eq. 0,0 then
!   step(2) = 6.225e-5
!   defi/local dec/r/1/1 {inputd(2)}
!   dec(1) = dec(1) + inputi(2)*step(2)/2
!   step(1) = step(2)/m$cos({dec(1)})
!   WRITE/KEYW INPUTD/D/3/2 -{step(1)},{step(2)}
!   write/out step values = {inputd(3)},{inputd(4)}
!endif


RUN /disk2/fried/midas/prog/drizzle.exe

!cop/dd {p1} lhcuts/r/1/2 {p2} lhcuts/r/1/2
!copy/dd {p1} *,3,history,lhcuts {p2} 
write/descr {p2} ident "{p1} rebinned"

write/out
write/out output -> {p2}
write/out


exit:
