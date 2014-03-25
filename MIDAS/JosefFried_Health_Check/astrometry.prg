!
!	astrometrie.prg
!	changed version of dwarf.prg by JW
!
! usage: 
! @@ astrometrie {frame root} {chipnr} {standardcat} {standardepoch} {matchingtol} 
!
!	This program uses the USNO catalogue to extract an object catalogue 
!       from LAICA-chips complete with RA and DEC. 

!       "Objects" are what we get out of SExtractor, while 
!	"stars" or "standard stars" will allways be taken from some catalogue. 

!	The catalogue of standard stars has to be provided beforehand. 
!       It is assumed, that it is an ASCII file like the one you get, when 
!       you retrieve a part of the USNO out of SKYCAT. The header will be
!	treated properly. 

!	The WCS is set to an approximation using the asumption, that north is 
!       up and the pointing determines the center of the field.

!	The keyword STARSEQ determines which stars are used for the interactive 
!       correction of the pointingoffset and the rotation of the field of view. 
!	"tol" is the tolerance in standard coordinates that you want to allow 
!       for when joining the standards and the objects table (0.00001 corresponds 
!       roughly to 2 arcseconds). Specifying big tolerances does not necessarily 
!       mean, that you get more stars, because each star that was matched two times 
!       is excluded.
!
!	JW 4.01
!       modified JWF 11.04...03.05
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! use of tables:
!                in	       out
! lasUSNOtbl                 {root}std     standard stars

! laCOSMOStbl                {root}std     standard stars 

! lawcs        {root}std     {name}temp    linear wcs approx.      

! lasex                      {name}        object x/y measurements

! lajoin       {name}temp
!	       {name}        {name}cat     joined table with x/y for standards
 
! laastromet   {name}cat     keywords      TERMS,CX,CY,BX,BY,AL_DE0 -> frame 

! laastrocomp  {name}        {name}t       calibrated result
!					   :R_A :DEC -> {name}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



set/format F20.10 I1
define/parameter p1 ? ? "Frame root:"
IF p1 .eq. "help" THEN
   write/out
   write/out "astrometry                         JWF 22 march 2005" 
   write/out
   write/out Purpose: 
   write/out	      Use standard star coordinates from the output tables of laprep 
   write/out	      to compute an astrometric solution and apply this to an object
   write/out	      catalogue created by SExtractor. The astrometric
   write/out	      solution is stored as frame descriptors for further use.
   write/out
   write/out Syntax:
   write/out @@ astrometry frame_root chip_number [tolerance] [mode]
   write/out 
   write/out frame root = root of frame name  without chip#   no default
   write/out chip_number = 1,2,3, or 4                        no default
   write/out tolerance =  tolerance in generalized coords     def=0.00002,0.00002
   write/out mode = interactive/no background subtraction     def=yn
   write/out
   write/out Example:
   write/out @@ astrometry M33V0049 1 
   write/out @@ astrometry jun120033 2
   write/out 
   write/out Notes:
   write/out The frames have to be passed first through laprep!
   write/out
   write/out astrometry first sets an approximate WCS. Using this, coarse
   write/out coordinates of the standard stars are calculated. 
   write/out An object catalog using SEXtractor is created. If the coarse 
   write/out coordinates agree within a given tolerance with the
   write/out positions measured by SExtractor, they are used to calculate 
   write/out an  astrometric solution. This is done in 2 steps: a 2.order fit
   write/out is used to throw out stars deviating more than 1 arcsec and then
   write/out a 3.order fit is done. This fit is also applied to the 
   write/out object catalog and stored as descriptors of the frame for further
   write/out use.
   write/out astrometry creates some tables which may be deleted. The
   write/out astrometrically calibrated object catalog is called root_chipnumber.tbl	     
   return
ENDIF




define/parameter p2 ? N "Number for the chip:"
define/parameter p3 0.00002,0.00002 ? "Tolerance in Xi,Eta:"
define/parameter p4 "yn" C "mode = yn interactive/no background subtraction"

define/local root/c/1/20 {p1}
define/local chip/i/1/1 {p2}
define/local tol/d/1/2 {p3}
define/local mode/c/1/2 {p4}

define/local N/i/1/1
define/local step/d/1/2 0 
define/local name/c/1/10
define/local OTIME/r/1/7
define/local ok/c/1/1
define/local exist/i/1/1
define/local descr/c/1/3 phi


! 
! Compute the first approximation for the WCS.
! linear approximation+dist.correction 
! input  {root}std
! output {root}temp 
! relevant columns: :Xicor :Etacorr
!
write/out
write/out "-------------------------------------------------- "
write/out Setting up WCS in linear approximation
write/out "-------------------------------------------------- "
write/out

        set/for i1
	@@ lawcs {root} {chip} 





! 
! Create the object catalog
! this step must always be done anew even if it has been done before
! since the WCS might change after step2
! run Sextractor 
! outputtable {root}_{chip}
!

write/out
write/out "-------------------------------------------------- "
write/out SExtracting objects
write/out "-------------------------------------------------- "
write/out


	set/for i1
	@@ lasex {root} {chip}




!
! subtract sky using background of SExtractor [optional]

!
if mode(2:2) .eq. "y" then
   comp/key inputi = m$existd("{root}_{chip}.bdf","back")
 
  if inputi(1) .eq. 0 then 
     write/out
     write/out "-------------------------------------------------- "
     write/out Subtracting sky
     write/out "-------------------------------------------------- "

      set/for I1
      comp {root}_{chip} = {root}_{chip}-background.fits
      loa/ima {root}_{chip}
      @@ backdet
      wd {root}_{chip} back/i/1/1 1
  endif
endif


!
! Join the standards and the objects table
! input: {root}temp     the std table, with the linear approx. of WCS
!        {root}_{chip}         the output of SExtractor
! output {root}_{chip}cat      all common objects within the defined tolerance

write/out
write/out "-------------------------------------------------- "
write/out Join tables with standards and objects  
write/out "-------------------------------------------------- "

set/for I1
@@ lajoin {root} {chip} {tol(1)} {tol(2)}




!
! compute astrometric solution
!
! works on table {root}_{chip}cat
! the solution is stored as descriptors to this table


write/out
write/out "-------------------------------------------------- "
write/out Starting astrometry 
write/out "-------------------------------------------------- "

set/cont astromet
set/for I1


@@ laastromet {root} {chip} -3



! Clean up.
!$rm {root}_{chip}cat.tbl {root}_{chip}std.tbl 



exit:





