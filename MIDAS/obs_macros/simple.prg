!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! .COPYRIGHT	(C) IAA-CSIC
! .IDENT	simple_expos.prg
!
! .CALL		@@ simple_expos icatalog itime=[tot_pointing,tot_single,single]
!				   object
! .PARAMETERS	
!		P1=itime=intergration time in seconds:
!			tot_pointing,tot_single,single
!			tot_pointing = total integration time for pointing
!			tot_single = total integration time for one image
!			single = integration time of a single readout
!		P2:object = object name for identifier
!               P3:filter = filter name for exposition
!               P4:pointing = pointing number to be stored in descriptor POINT_NO
!
! .PURPOSE	A procedure to take non dithered images for one pointing 
!		with PANIC
! .ENVIRONMENT	MIDAS
! .AUTHOR	jmiguel@iaa.es     
! .KEYWORDS	non dithering, pointing
! .COMMENTS	 
!
! .CONVENTIONS	!! --> command of line necessary
!		!* --> command not needed but could be useful for testing
!
! .VERSION	1.0	03.03.11
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!
!!!!!!!!!!!!!!!!!! SETUP AND HELP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! command line parameter setup
crossref  itime object 


	! help text
if p1(1:4) .eq. "help" then

write/out
write/out "simple_expos.prg"
write/out "call: @@ simple_expos  itime=[tot_pointing,tot_single,single]"
write/out "			      filter pointing"
write/out
write/out "This program takes dithered images for one pointing for OMEGA2000:"
write/out "The area with no overlap is maximally 26 arcsec on each side"
write/out 
write/out "The command line parameters are:"
write/out "P1:itime = [tot_p,tot_s,s]intergration time in seconds:"
write/out "	tot_pointing,tot_single,single"
write/out "	tot_pointing = total integration time for pointing"
write/out "	tot_single = total integration time for one image"
write/out "	single = integration time of a single readout"
write/out "P2:object = object name for identifier"
write/out "P3:filter = filter name for exposition"
write/out "P4:pointing = pointing number to be stored in descriptor POINT_NO"

write/out

return		! back to MIDAS session

endif	


	! define command line parameters
define/parameter P1 ? N "Enter the integration times in secs[tot,integ,single] :"
define/parameter P2 ? C "Enter an object name :"
define/parameter P3 ? C "Enter an filter name :"
define/parameter P4 1 N "Enter the identification number of the current pointing :"


	! extract integration times
define/local itime/r/1/3 {P1}
define/local tot_pointing/r/1/1 {itime(1)}
define/local tot_single/r/1/1 {itime(2)}
define/local single_time/r/1/1 {itime(3)}
define/local i_p/i/1/1 0
	
	! for return value of abort-file-check
define/local abort_check/i/1/1		

	! define keyword for path and image name
define/local pathname_ima/c/1/200

	! get pointing number
define/local point_no/i/1/1 {P4}

define/local loop/i/1/1 
define/local counter/i/1/1 

set/format I1 F4.2

	! display input parameters
write/out
write/out "** ANY ** dither pattern will be used to take images..."	
write/out "The total integration time for the poining will be:{tot_pointing} secs"
write/out "The total integration time for a single image is: {tot_single} secs"
write/out "The integration time for a single frame is: {single_time} secs"
write/out "The FILTER for expositions will be: {P3} "
write/out "The descriptor POINT_NO will contain the pointing number:{point_no}"
write/out

set/format

!!!!!!!!!!!!!!!!!!!!!! ENVIRONMENT VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

define/local geirslstabort/c/1/256  " "
geirslstabort = M$SYMBOL("GEIRSLSTABORT")

define/local tecs_script/c/1/256  " "
tecs_script = M$SYMBOL("TECS_SCRIPT")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

	! remove existing abort file
	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("{geirslstabort}")
if abort_check .eq. 1 then
  $rm {geirslstabort} 
endif

!!!!!!!!!!!!!!!!!! image / integration parameters  !!!!!!!!!!!!!!!!!!!!!!!!

	! calculate image repitions
define/local rep_image/i/1/1		! total number of images to be taken
     rep_image = m$nint(tot_pointing/tot_single)

define/local rep_integrate/i/1/1	!number of images for single-image-integration
rep_integrate = m$nint(tot_single/single_time)



!.............. camera commands now ............

	! write initial image descriptors
$cmd_panic_new counter DITH_NO set 1	! set dither counter to start position
$cmd_panic_new counter POINT_NO set {point_no}	! set pointing no to its value
$cmd_panic_new counter EXPO_NO clear		! clear exposure counter-->EXPO_NO=1

	! set single image parameters
set/format I1
$cmd_panic_new crep {rep_integrate}
$cmd_panic_new itime {single_time}
$cmd_panic_new filter {P3}
$cmd_panic_new sync

!
!
!-------------------------------------------------------------
	! do loop over number of images
do loop = 1 {rep_image} 1

set/format I1
write/out "Taking image {loop} of {rep_image}..."	

	! write object
$cmd_panic_new object {P2}:{loop}/{rep_image}


$cmd_panic_new read
$cmd_panic_new sync


	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("{geirslstabort}")
if abort_check .eq. 1 then
  write/out "Program is aborted..."
  $rm {geirslstabortp}
  $play -q $GEIRS_DIR/SOUNDS/crash.au
  goto exit
endif


$cmd_panic_new save 


	! handle image descriptors
$cmd_panic_new counter DITH_NO inc			
$cmd_panic_new counter EXPO_NO clear		! reset exposure counter 


	! increment counter
counter = counter+1


enddo
!-------------------------------------------------------------

$cmd_panic_new sync    ! wait for last save



!!!!!!!!!!!!!!!!! CLOSE DOWN THINGS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write/out
write/out "All images for pointing are finished..." 
write/out
$play -q $GEIRS_DIR/SOUNDS/gong.au

exit:


write/out
write/out "Observing of current pointing done..."

return

