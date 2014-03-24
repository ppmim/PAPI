!
!	astrometrie.prg
!	changed version of dwarf.prg by JW
!
! usage: astrometrie {frame} {standardcat} {chipnr} {matchingtol} {offset}
!
!	This program uses the USNO catalogue to extract an object catalogue from LAICA-chips 
!	complete with RA and DEC. "Objects" are what we get out of SExtractor, while "stars"
!	or "standard stars" will allways be taken from some catalogue. 
!	The catalogue of standard stars has to be provided beforehand. It is assumed, that it
!	is an ASCII file like the one you get, when you retrieve a part of the USNO out of SKYCAT. 
!	The header must be deleted. 
!	The WCS is set to an approximation using the asumption, that north is up and the
!	pointing determines the center of the field.
!	The keyword STARSEQ determines which stars are used for the interactive correction of 
!	the pointingoffset and the rotation of the field of view. "tol" is the tolerance in 
!	standard coordinates that you want to allow for when joining the standards and the 
!	objects table (0.00001 corresponds roughly to 2 arcseconds). 
!	Specifying big tolerances does not necessarily mean, that you get more stars, 
!	because each star that was matched two times is excluded.
!
!	JW 4.01
!
! major modifications JWF Sept/Oct 2002
!


set/format F20.10 I5
define/parameter p1 ? ? "Frame name:"
IF p1 .eq. "help" THEN
write/out @@ astrometrie frame standard_table chip_number tolerance offset backsub
return
ENDIF

define/parameter p2 ? ? "Name of the ASCII file:"
define/parameter p3 ? N "Number for the chip:"
define/parameter p4 0.00004,0.00004 ? "Tolerance in Xi,Eta:"
define/parameter p5 0,0 ? "Offset in x,y:"
define/parameter p6 n c "subtract background [y]/n ?"



defi/local sub/c/1/1 {p6}
define/local STARSEQ/i/1/1 5
define/local tol/d/1/2 {p4}
define/local N/i/1/1
define/local step/d/1/2 0 
define/local tols/d/1/2
define/local name/c/1/20
define/local OTIME/r/1/7
define/local ok/c/1/1
define/local exist/i/1/1
define/local descr/c/1/3 phi


write/key name {p1}


write/out
write/out
write/out "------------------------------------------------ "
write/out working on image: {name}  
write/out "------------------------------------------------ "
write/out


! Compute the standard coordinates in the standard catalogue.
!  {p2}.cat (=USNOcatalogue) -> {p1}std.tbl 
write/out
write/out Preparing the standard star catalogue ...
write/out
  
@@ laUSNOtbl {name} {p2}
!cco 
!comp/tab {p1}std :z = seq*1.0
!loa/tab {p1}std :xi :eta :z 1 3 1


! Compute the first approximation for the WCS.
! uses {p1}std to determine rotation of field and pointing offset  
! selects from {p1}std.tbl the stars which are on the chip -> {p1}temp.tbl
write/out
write/out Setting the WCS for frame {name} ...
write/out
@@ lawcs {name} {STARSEQ} {p3} {p5}
!  @@ lawcs2 {name} {STARSEQ} {p3} {p5}



write/out coarse positions ... blue

clea/chan 0
clear/chan over
loa {p1}
sel/tab {p1}temp all
loa/tab {p1}temp :xi :eta ? 1 4 5







! Apply distortion to xi and eta
! convert radius -> mm 
! note dist is in %

comp/tab {p1}temp :radius = 125/9.1e-3*sqrt(:xi**2+:eta**2)
comp/tab {p1}temp :cb = abs(:xi/(sqrt(:xi**2+:eta**2)))
comp/tab {p1}temp :sb = abs(:eta/(sqrt(:xi**2+:eta**2)))
comp/tab {p1}temp :dist = 9.7947e-3-8.2364e-4*:radius+1.0516e-4*:radius**2
comp/tab {p1}temp :xidist = :dist*:cb/100.
comp/tab {p1}temp :etadist = :dist*:sb/100.

if {p3} .eq. 1 then 
   comp/tab {p1}temp :xidis = :xi*(1.+:xidist)  
   comp/tab {p1}temp :etadis = :eta*(1.+:etadist)  
endif
if {p3} .eq. 2 then 
    comp/tab {p1}temp :xidis = :xi*(1.+:xidist)  
    comp/tab {p1}temp :etadis = :eta*(1.+:etadist)  
endif
if {p3} .eq. 3 then 
   comp/tab {p1}temp :xidis = :xi*(1.+:xidist)  
   comp/tab {p1}temp :etadis = :eta*(1.+:etadist)  
endif  
if {p3} .eq. 4 then 
   comp/tab {p1}temp :xidis = :xi*(1.+:xidist)  
   comp/tab {p1}temp :etadis = :eta*(1.+:etadist) 
endif

loa/tab {p1}temp :xidis :etadis ? 1 4 6
write/out distorted positions ...  yellow







! Create the object catalog
! output {p1}.tbl

inputi(1) = m$exist("{p1}.tbl")
if {inputi(1)} .eq. 0 then 
   write/out
   write/out {p1}.tbl does not exist - start SExtractor ...
   write/out
   @@ lasexzwerg {name} 
endif

if {inputi(1)} .eq. 1 then 
   write/out
   write/out {p1}.tbl exists - skip SExtractor ...
   write/out
endif



! Now we join the standards and the objects table
!  {p1}temp + {p1} -> {p1}cat 
write/out
write/out Creating a joined catalog for standards and the objects on frame {name} ...
write/out
@@ lajoin {name} {tol(1)} {tol(2)}

write/out
write/out Loading the standards which will be used for astrometry ... white squares
write/out
loa/tab {p1}cat :xidis_1 :etadis_1 ? 0 4 2

wk action/c/1/1 y
inquire/key action "Continue [y]/n"
if action(1:1) .eq. "n" then
   goto exit
endif

! Now we calculate the plate solutions for chip
! works on {p1}cat
! nstars <5 linear fit
! nstars <40 second order
! nstars >40 third order

write/out
write/out Calculating the plate solution 
write/out


set/cont astromet
@@ laastromet {name} y


! Now we produce the final object catalogue with coordinates in RA and DEC 
write/out
write/out Apply plate solution to the object catalogue ...
write/out
@@ laastrocomp {name}


! Clean up.
!$rm {name}cat.tbl {name}std.tbl 


! subtract background

a:


if "{p6}" .eq. "y" then 
   write/out   
   write/out  use background from Sextractor to subtract background ...
   write/out
   indisk/fits background.fits {p1}back.bdf
   compute/image {p1} = {p1} - {p1}back
endif





write/out
write/out All done ... 
write/out


exit:






