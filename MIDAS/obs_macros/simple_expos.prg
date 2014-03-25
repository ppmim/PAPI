!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! .COPYRIGHT	(C) IAA-CSIC
! .IDENT	single_expos.prg
!
! .CALL		@@ simple_expos icatalog itime=[tot_pointing,tot_single,single]
!				   object
! .PARAMETERS	P1=icatatlog = name of image catalog for images
!		P2=itime=intergration time in seconds:
!			tot_pointing,tot_single,single
!			tot_pointing = total integration time for pointing
!			tot_single = total integration time for one image
!			single = integration time of a single readout
!		P3:object = object name for identifier
!
! .PURPOSE	A procedure to take non dithered images for one pointing 
!		with PANIC
! .ENVIRONMENT	MIDAS
! .AUTHOR	jmiguel@iaa.es     
! .KEYWORDS	non dithering, pointing
! .COMMENTS	An image catalog is opened, which will contain all images taken.
!		The procedure will use the image names, filter etc specified in the
!		GUI. At the current observing position, a non dithered image 
!		for the specified total integration time will be taken. 
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
crossref icatalog itime object 


	! help text
if p1(1:4) .eq. "help" then

write/out
write/out "dither_pointing.prg"
write/out "call: @@ dither_pointing icatalog itime=[tot_pointing,tot_single,single]"
write/out "			     start_pos pointing	"
write/out
write/out "This program takes dithered images for one pointing for OMEGA2000:"
write/out "The area with no overlap is maximally 26 arcsec on each side"
write/out 
write/out "The command line parameters are:"
write/out "P1:icatatlog = name of image catalog for images"
write/out "P2:itime = [tot_p,tot_s,s]intergration time in seconds:"
write/out "	tot_pointing,tot_single,single"
write/out "	tot_pointing = total integration time for pointing"
write/out "	tot_single = total integration time for one image"
write/out "	single = integration time of a single readout"
write/out "P3:object = object name for identifier"
write/out "P4:pointing = pointing number to be stored in descriptor POINT_NO"
write/out

return		! back to MIDAS session

endif	


	! define command line parameters
define/parameter P1 ? C "Enter a name for the image catalog :"
define/parameter P2 ? N "Enter the integration times in secs[tot,integ,single] :"
define/parameter P3 ? C "Enter an object name :"
define/parameter P4 1 N "Enter the identification number of the current pointing :"

	! define keyword which contains the path where pipe_files are stored
define/local file_path/c/1/80 "/home/panic/tmp"
define/local isodate/c/1/32 " "
define/local fctrl/i/1/2 0,0
isodate = m$isodate()
open/file tel_pos_{isodate}.log WRITE fctrl
write/file {fctrl(1)} {isodate} Log of telescope movements for {P5} (pointing {P7})

	! extract catalog name
define/local icatalog/c/1/80 {P1}

	! extract integration times
define/local itime/r/1/3 {P2}
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
write/out "A dither pattern will now be used to take images..."	
write/out "All images will be listed in the image catalog: {icatalog}"
write/out "The total integration time for the poining will be:{tot_pointing} secs"
write/out "The total integration time for a single image is: {tot_single} secs"
write/out "The integration time for a single frame is: {single_time} secs" 
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

!!!!!!!!!!!!!!!!! PASS INFORMATION TO PIPELINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! activate image catalog: all new *.fits images will be added 
create/icat {icatalog} null	! no entries in catalog
set/icat {icatalog}
	

define/local icat_path/c/1/180			! holds path and name of icat
define/local string_length/i/1/1		! length of path
define/local file_id_1/i/1/2			! file id for open/file
define/local file_id_2/i/1/2			! file id for open/file
define/local file_id_3/i/1/2			! file id for open/file


$pwd | write/keyword icat_path

string_length = m$len(icat_path)
string_length = string_length+1			! start at this position

icat_path({string_length}:) = "/{icatalog}"

	! file contains path and name of active image catalog
open/file {file_path}/MacroLstIcat write file_id_1
write/file {file_id_1(1)},key icat_path
close/file {file_id_1(1)}

set/format f4.1

	! file which contains the integration time per image
open/file {file_path}/MacroIntTime write file_id_2
write/file {file_id_2(1)} {tot_single}
close/file {file_id_2(1)}


	! file contains only name of active image catalog
open/file {file_path}/MacroLstIcatName write file_id_3
write/file {file_id_3(1)},key icatalog
close/file {file_id_3(1)}


 

!!!!!!!!!!!!!!!!!! image / integration parameters  !!!!!!!!!!!!!!!!!!!!!!!!

	! calculate image repitions
define/local rep_image/i/1/1		! total number of images to be taken
     rep_image = m$nint(tot_pointing/tot_single)

define/local rep_integrate/i/1/1	!number of images for single-image-integration
rep_integrate = m$nint(tot_single/single_time)



!.............. camera commands now ............

	! write initial image descriptors
$cmd_panic_new counter DITH_NO set 0	! set dither counter to start position
$cmd_panic_new counter POINT_NO set {point_no}	! set pointing no to its value
$cmd_panic_new counter EXPO_NO clear		! clear exposure counter-->EXPO_NO=1


	! set single image parameters
set/format I1
$cmd_panic_new crep {rep_integrate}
$cmd_panic_new itime {single_time}
$cmd_panic_new sync

!
!
!-------------------------------------------------------------
	! do loop over number of images
do loop = 1 {rep_image} 1

set/format I1
write/out "Taking image {loop} of {rep_image}..."	

	! write object
$cmd_panic_new object {P3}:{loop}/{rep_image}


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

write/out "VOY POR AQUI (a)"

	! add file to image catalog
if loop .ge. 2 then
  wait/sec 5
  !set/midas output=logonly
  write/out "VOY POR AQUI (B1)"
  $cmd_panic_new last	| awk '{print $2}' | write/keyword pathname_ima 
	! add file to icat
  write/out "VOY POR AQUI (B2) IMAGEN= {pathname_ima}"      
  add/icat {icatalog} {pathname_ima}
  write/out "VOY POR AQUI (B3)"
  
  set/midas output=yes
endif

write/out "VOY POR AQUI (c)"

$cmd_panic_new save -i


	! handle image descriptors
$cmd_panic_new counter DITH_NO incr			! increment dither counter by 1
$cmd_panic_new counter EXPO_NO clear		! reset exposure counter 


	! increment counter
counter = counter+1


enddo
!-------------------------------------------------------------

$cmd_panic_new sync    ! wait for last save

set/midas output=logonly

	! add last file to image catalog
$cmd_panic_new last | awk '{print $2}' | write/keyword pathname_ima 
	! add file to icat
add/icat {icatalog} {pathname_ima}

set/midas output=yes

        ! close tele-log-file
isodate = m$isodate()
write/file {fctrl(1)} {isodate} done


!!!!!!!!!!!!!!!!! CLOSE DOWN THINGS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write/out
write/out "All images for pointing are finished..." 
write/out
$play -q $GEIRS_DIR/SOUNDS/gong.au

exit:

close/file {fctrl}

	! deactivate image catalog
clear/icat

	! give out catalog
write/out
write/out "The following images were taken and are listed in the image catalog {icatalog}:"
write/out

read/icat {icatalog}

write/out
write/out "Observing of current pointing done..."

return

