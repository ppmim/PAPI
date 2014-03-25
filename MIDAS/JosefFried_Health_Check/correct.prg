!
! correct wfi data
!
! jwf aug 99
!

defi/par p1 ? n "1.file # = ? "
defi/par p2 ? n " last file # = ? "
defi/par p3 dat1 ima " prefix for infiles = ? "
defi/par p4 dat1c c " prefix for outfiles = ? "
defi/par p5 flat1 ima " flat field = ? "
defi/par p6 none ima " global sens = ? "
defi/par p7 y c " delete input file [def=y] = ? "


defi/local mfile/i/1/1 0
set/format i4

do mfile = {p1} {p2}
	write/out file {mfile} ...
	if p6 .eq. "none" then
	        comp/ima {p4}{mfile} = {p3}{mfile}*{p5}
	else
		comp/ima {p4}{mfile} = ({p3}{mfile}*{p5})/{p6}
	endif
	
	if p7 .eq. "y" then
		write/out delete file {p3}{mfile}.bdf ...
		$rm {p3}{mfile}.bdf
	endif

enddo