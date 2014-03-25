!
!  archive.prg
!

if p1(1:4) .eq. "help" then

	write/out use as 	@@ archive catalog outfile

	write/out to archive data:

	write/out 1.) create catalogue of frames using crea/icat 
	write/out "	example: crea/icat frames ../*/*.bdf"
	write/out "	this creates a catalogue of all bdfs in directories  ../*/"

	write/out 2.) use @@ archive frames frameslist

	write/out "	this will scan all files, store descriptors in 	"
	write/out "	a) an ascii file (frameslist.ascii in the example)"
	write/out "	b) an sylk file  (frameslist.sylk in the example)"
	write/out "	    which can be read with StarOffice"

	goto exit
endif


defi/par p1 frames char catalog

defi/par p2 list char output_file
write/keyw in_a/c/1/20 {p1}
write/keyw out_a/c/1/20 {p2}


run FMP:archive


exit: