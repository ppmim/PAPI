!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! astroprep
!
! prepares the LAICA files for use in astrometry.prg
!
!
! gets coords of field center from descriptor O_POS 
!
! computes precession 
!
! checks pointing
!
! prepares standard star table
! 
!
!	JWF mar/2005
!
! input: frame root
!
! output: descriptor alpha2000,delta2000
!	  table {ascii}std
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


define/parameter p1 ? ? "root of names:"
define/parameter p2 ? ima "ascii file:"
define/parameter p3 2000 n "Equinox of standard stars [2000]: " 



defi/local root/c/1/20 {p1}
defi/local asci/c/1/20 {p2}
defi/local equinox/R/1/2 {p3}


set/for i1
defi/local file/c/1/20 "{p1}_1"
set/for I4




del/dis 
del/grap
crea/dis 0 840,840,430,200
crea/grap 0 600,480,0,500



write/out
write/out "-------------------------------------------------- "
write/out working on image: {file}
write/out "-------------------------------------------------- "
write/out

define/local delta2000/d/1/1
define/local alpha2000/d/1/1

defi/local ra/d/1/1 0.
defi/local de/d/1/1 0.



cop/dk {file} O_TIME/d/1/1 equi
cop/dk {file} O_POS/d/1/1 ra
cop/dk {file} O_POS/d/2/1 de



defi/local rah/i/1/1 0
defi/local ram/i/1/1 0
defi/local ras/r/1/1 0.
defi/local ddeg/i/1/1 0
defi/local dmin/i/1/1 0
defi/local dsec/r/1/1 0.

rah = ra/15
ram = (ra/15-rah)*60
ras = ra/15*3600.-rah*3600.-ram*60

if de .lt. 0 then
   de = -de
   ddeg = de
   dmin = (de-ddeg)*60.
   dsec = de*3600-ddeg*3600-dmin*60
   if ddeg .ge. 1 then
      ddeg = -ddeg
   else 
	if dmin .gt. 1 then
	   dmin=-dmin
	else 
	   dsec = -dsec
	endif
   endif
   de = -de 
else
   ddeg = de
   dmin = (de-ddeg)*60.
   dsec = de*3600-ddeg*3600-dmin*60
endif	


set/for I2 F5.2
write/out alpha = {rah},{ram},{ras} 
write/out delta = {ddeg},{dmin},{dsec}
write/out Epoch of observation = {equi}

comp/prec {rah},{ram},{ras} {ddeg},{dmin},{dsec} {equi} {equinox} 

defi/local alpha/d/1/1 0.
comp/key alpha = outputr(1)+outputr(2)/60.+outputr(3)/3600.
comp/key alpha = alpha*15.

defi/local delta/d/1/1 0.

if de .gt. 0 then 
   delta = outputr(4)+outputr(5)/60.+outputr(6)/3600.
else
   delta = m$abs(outputr(4))+m$abs(outputr(5)/60.)+m$abs(outputr(6)/3600.)
   delta = -delta
endif



comp/key inputi = m$existk("offpoint")
 
if inputi(1) .eq. 0 then 
      write/key offpoint/r/1/2 0.,0.
endif
write/out Offset {offpoint(1)},{offpoint(2)} used


comp/key alpha = alpha-{offpoint(1)}
comp/key delta = delta+{offpoint(2)}

set/for f15.6
write/out
write/out center of field  @ equinox {equinox(1)}: {alpha} {delta} [decimal degrees]




!
! note alpha2000/delta2000 are used in the following
! O_POS is not used anymore!
!

write/desc {root}_1 alpha2000/d/1/1 {alpha}
write/desc {root}_1 delta2000/d/1/1 {delta}
write/desc {root}_2 alpha2000/d/1/1 {alpha}
write/desc {root}_2 delta2000/d/1/1 {delta}
write/desc {root}_3 alpha2000/d/1/1 {alpha}
write/desc {root}_3 delta2000/d/1/1 {delta}
write/desc {root}_4 alpha2000/d/1/1 {alpha}
write/desc {root}_4 delta2000/d/1/1 {delta}


!
! set individual offsets of CCDs 
! positions ok f2004
!
write/key offccdx/r/1/4 -4.8,0.4,-3.7,-9.
write/key offccdy/r/1/4 -.5,-5.,-12.5,-16.5


!
! 2.step
! 
! Compute the standard coordinates in the standard catalogue.
! extracts relevant part of the total USNO table
! output {root}std.tbl  stored in key asc
! relevant columns: Rmag Bmag Ra DEC Xi Eta

write/out
write/out "------------------------------------------------ "
write/out Prepare standard star catalogue
write/out "------------------------------------------------ "
write/out
	set/for I1
	@@ laCOSMOS2tbl {root} {asci}

!
! 3.step
! correct pointing
!

write/out
write/out "------------------------------------------------ "
write/out Check + correct pointing
write/out "------------------------------------------------ "
write/out

	@@ lapointing {root} 1 
! we have to remake table
	@@ laCOSMOS2tbl {root} {asci}


exit: