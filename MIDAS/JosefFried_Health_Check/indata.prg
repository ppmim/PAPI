!
! indata
! reads data from DAT, removes bias, creates mosaic
!
! jwf aug 99
!


defi/par p1 ? n " 1.file = ? "
defi/par p2 ? n " last file = ? "
defi/par p3 dat1 ima " prefix for output files = ? "
defi/par p4 yy ima " read from tape y/n  + reduce y/n [def= yy] = ? "

defi/local mfile/i/1/1 0
set/format i4

if p4(1:1) .eq. "y" then
	do mfile = {p1} {p2}
		inta/fits {mfile} inp /dev/rmt/0bn SOA
	enddo
	write/out data read ... 
endif



if p4(2:2) .eq. "y" then
	do mfile = {p1} {p2}
		write/out removing bias from files {mfile} ...
		@@ rmbias inp{mfile}
		@@ mosaic rbinp{mfile} {p3}{mfile} msk
		copy/dd inp{mfile} *,3 {p3}{mfile}
		cop/dd inp{mfile} ident {p3}{mfile} ident
		$rm inp{mfile}.bdf
		$rm inp{mfile}*.bdf
		$rm rbinp{mfile}*.bdf
	enddo
endif




exit: