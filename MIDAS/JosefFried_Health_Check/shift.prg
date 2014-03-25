!
! shift.prg
!
! to set star/step to align images
!

defi/par p1 raw1 ima "prefix" 
defi/par p2 ? n "first file # = ? "
defi/par p3 ? n "last file # = ? "

defi/local n/i/1/1 1
set/format i4
loa {p1}{p2} 

cen/gauss

defi/local start0/r/1/2 0.,0.
defi/local start/r/1/2 0.,0.

comp/key start0(1) = outputr(5)
comp/key start0(2) = outputr(6)

do n = {p2}+1 {p3}
	loa {p1}{n} 
	cen/gauss
	comp/key start(1) = -outputr(5)+start0(1)
	comp/key start(2) = -outputr(6)+start0(2)
	write/desc {p1}{n} start/d/1/2 {start(1)},{start(2)}
enddo 

