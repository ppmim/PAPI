!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! .COPYRIGHT	(C) MPIA
! .IDENT	dither_pointing.prg
!
! .CALL		@@ dither_pointing icatalog itime=[tot_pointing,tot_single,single]
!				   object start_pos pointing offsets
! .PARAMETERS	P1=icatatlog = name of image catalog for images
!		P2=itime=intergration time in seconds:
!			tot_pointing,tot_single,single
!			tot_pointing = total integration time for pointing
!			tot_single = total integration time for one image
!			single = integration time of a single readout
!		P3:object = object name for identifier
!		P4:start_pos,tel_pos
!                       start_pos = continue dither pattern at this position (N)
!			tel_pos   = AQ   telescope is at position 1 (aequisition)
!                                 = PREV telescope is at position N-1  (default)
!		P5:pointing = pointing number to be stored in descriptor POINT_NO
!		P6:offsets = flag for integer or non-integer pixels offsets
!			     0 = integer pixel offsets ; 1 = non-integer offsets
!
! .PURPOSE	A procedure to take dithered images for one pointing 
!		with OMEGA2000
! .ENVIRONMENT	MIDAS
! .AUTHOR	Rene Fassbender     
! .KEYWORDS	dithering, pointing
! .COMMENTS	An image catalog is opened, which will contain all images taken.
!		The procedure will use the image names, filter etc specified in the
!		GUI. At the current observing position, a dithered image pattern
!		for the specified total integration time will be taken. A single
!		dither pattern contains 20 observing positions, after that the
!		pattern will be repeated with a slight offset for the first image.
!		The area with no overlap is at most 26 arcsec on each side.
!
! .CONVENTIONS	!! --> command of line necessary
!		!* --> command not needed but could be useful for testing
!
! .VERSION	1.0	16.12.02
!		2.0 	16.01.03	adjusted to CA environment on o2k
!		3.0 	17.01.03	implementation of more features P3, P4
!		4.0	12.03.03	implementation of: abort check, telescope
!					return check, write next file in icat,
!					write pointing keywords, integer pixel 
!					offsets, writing icat in file
!		4.1	16.03.03	use awk for return check
!		4.2	10.04.03	write catalog name in seperate file,
!					changed integer repetition offsets,
!					added key file_path, change awk to 1.line
!               5.0     03.04.04        back to old repitions pattern and
!                                       drizzle-compatible offsets  (HJR)
!                                       add 1 / N-1  option for starting position
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!
!!!!!!!!!!!!!!!!!! SETUP AND HELP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! command line parameter setup
crossref icatalog itime object start_pos pointing offsets


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
write/out "P4:start/telescope position = start_pos,tel_pos"
write/out "             start_pos = continue dither pattern at this position (#n, def=1)"
write/out "             tel_pos   = flag to indicate telescope position"
write/out "                       = AQ   telescope is at position 1 (aequisition)"
write/out "                       = PREV telescope is at position n-1  (default)"
write/out "P5:pointing = pointing number to be stored in descriptor POINT_NO"
write/out "P6:offsets = flag for integer or non-integer pixels offsets"
write/out "	        0 = integer pixel offsets ; 1 = non-integer offsets"
write/out

return		! back to MIDAS session

endif	


	! define command line parameters
define/parameter P1 ? C "Enter a name for the image catalog :"
define/parameter P2 ? N "Enter the integration times in secs[tot,integ,single] :"
define/parameter P3 ? C "Enter an object name :"
define/parameter P4 1,PREV C "Enter the sequence number n of start position :"
define/parameter P5 1 N "Enter the identification number of the current pointing :"
define/parameter P6 0 N "Enter the offset flag :"

	! define keyword which contains the path where pipe_files are stored
define/local file_path/c/1/80 "/disk-a/o2k/tmp"
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
define/local point_no/i/1/1 {P5}

	! get start positions
i_p = m$index(P4,",")
if i_p .eq. 0 then
   define/local start_position/i/1/1 {P4}
   define/local tel_flag/c/1/8 "PREV"
else
   define/local start_position/i/1/1 {P4(1:{i_p})}
   i_p = i_p + 1
   define/local tel_flag/c/1/8 {P4({i_p}:)}
endif
tel_flag = m$upper(tel_flag)
if tel_flag(1:2) .ne. "AQ"  then
   if tel_flag(1:4) .ne. "PREV" then
      write/out
      write/out "Input error:"
      write/out "Flag for current telescope position (parameter 4b) has to be AQ or PREV"
      write/out "   ... abort    "
      write/out
      $auplay /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
      goto exit
   endif
endif	

set/format I1 F4.2

	! display input parameters
write/out
write/out "A dither pattern will now be used to take images..."	
write/out "All images will be listed in the image catalog: {icatalog}"
write/out "The total integration time for the poining will be:{tot_pointing} secs"
write/out "The total integration time for a single image is: {tot_single} secs"
write/out "The integration time for a single frame is: {single_time} secs" 
write/out "The dither pattern will be started at position: {start_position}"
write/out "The identifier will contain the object name: {P3}"
write/out "The dither pattern will be started at position: {start_position}"
write/out "The descriptor POINT_NO will contain the pointing number:{point_no}"
write/out "The dither offsets will be (0=integer pixels ; 1=non-integer):{P6}"
write/out "Telescope is assumed to be at {tel_flag} position"
write/out "Telescope movements are logged into file tel_pos_{isodate}.log"
write/out

set/format

	! remove existing abort file
	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
if abort_check .eq. 1 then
  $rm /disk-a/o2k/tmp/geirsLstAbort 
endif
!
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

!
!!!!!!!!!!!!!!!!! DITHER SEQUENCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! define dither offsets in units of 1/10 arcsecs
	! all positions are defined relative to last position
set/format I1
define/local X_offset/i/1/20 0 all
define/local Y_offset/i/1/20 0 all

	! for in units of 1/3 pixel
if {P6} .eq. 1 then

!  define/local x_offset/i/1/20 0,101,-150,180,-227,-48,198,128,-375,225
!  write/keyword x_offset/i/11/10 -209,150,-90,317,-42,-216,143,-237,36,254
  
   write/keyword x_offset/i/1/10 0,100,-150,180,-226,-48,198,127,-348,198
   write/keyword x_offset/i/11/10 -208,150,-90,289,-15,-216,142,-237,36,253
 

!  define/local y_offset/i/1/20 0,101,57,-216,-29,150,-198,182,-87,252
!  write/keyword y_offset/i/11/10 -74,-309,162,-98,291,-99,-56,-192,360,-395

   write/keyword y_offset/i/1/10 0,100,57,-216,-29,150,-198,181,-87,216
   write/keyword y_offset/i/11/10 -38,-318,171,-97,291,-99,-55,-192,324,-322
 
  
  define/local X_back/I/1/1 -136
  define/local Y_back/i/1/1 160

  write/out "Fractional pixel offsets for drizzle were defined..."
!..............................................................................

else	! integer pixel offsets	: multiples of 9 = 2 pixels

!  define/local x_offset/i/1/20 0,99,-153,180,-225,-45,198,126,-378,225
!  write/keyword x_offset/i/11/10 -207,153,-90,315,-45,-216,144,-234,36,252

   write/keyword x_offset/i/1/10 0,99,-153,180,-225,-45,198,126,-351,198
   write/keyword x_offset/i/11/10 -207,153,-90,288,-18,-216,144,-234,36,252
      
!  define/local y_offset/i/1/20 0,99,54,-216,-27,153,-198,180,-90,252
!  write/keyword y_offset/i/11/10 -72,-306,162,-99,288,-99,-54,-189,360,-396

   write/keyword y_offset/i/1/10 0,99,54,-216,-27,153,-198,180,-90,216
   write/keyword y_offset/i/11/10 -36,-315,171,-99,288,-99,-54,-189,324,-324

  define/local X_back/I/1/1 -135
  define/local Y_back/i/1/1 162

  write/out "Integer pixel offsets were defined..."

endif

  ! define repetion offsets (always integer pixels) --> modified to yield sum=0
define/local x_repetition/i/1/20 0,-54,72,27,-90,27,54,-72,72,-54
write/keyword x_repetition/i/11/10 27,54,-27,-54,45,-54,72,-99,99,-45

define/local y_repetition/i/1/20 54,45,-99,99,-72,72,-54,27,27,-72
write/keyword y_repetition/i/11/10 45,-72,72,-72,27,54,-54,27,27,-81



!!!!!!!!!!!!!!!!!! image / integration parameters  !!!!!!!!!!!!!!!!!!!!!!!!

	! calculate image repitions
define/local rep_image/i/1/1		! total number of images to be taken
rep_image = m$nint(tot_pointing/tot_single)

define/local rep_integrate/i/1/1	!number of images for single-image-integration
rep_integrate = m$nint(tot_single/single_time)
!
!!!!!!!!!!!!!!!!! TELESCOPE MACRO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! keywords
define/local loop/i/1/1 
define/local counter/i/1/1 1		! for counting to 20
define/local pattern_reps/i/1/1 2	! for counting how many times 20
define/local tel_return/i/1/1 0		! character keyword for telescop returns:0,-1


 	! find the proper start position
if start_position .gt. 400  then
   start_position = start_position-(start_position/400)*400
   if start_position .eq. 0  then
      start_position = 400
   endif
endif
if start_position .le. 20 then
	counter = start_position
else
	pattern_reps = start_position/20
	counter = start_position - (20*pattern_reps)
	pattern_reps = pattern_reps+2			! pattern for next 20 block
endif

if counter .eq. 0 then
	counter = 20
	pattern_reps = pattern_reps-1	! set to next block
endif


!!!!!!!!!!!!!!!!!!  move telescope if not at PREV position  !!!!!!!!!!!!!!!!

if tel_flag(1:2) .eq. "AQ"  then
   set/format I3
   define/local x_move/i/1/1 0
   define/local y_move/i/1/1 0
   define/local ic/i/1/1 0
   define/local ir/i/1/1 0
   define/local xy_test/i/1/1 0

   ir = pattern_reps - 1
   if ir .ge. 2 then
      do ic = 2 {ir}
         x_move = x_move + x_repetition({ic})
         y_move = y_move + y_repetition({ic})
      enddo
   endif
   ir = counter - 1
   if ir .ge. 1 then
      do ic = 1 {ir}
         x_move = x_move + x_offset({ic})
         y_move = y_move + y_offset({ic})
      enddo
   endif
   xy_test = m$abs(x_move)+m$abs(y_move)
   if xy_test .ne. 0 then
      write/out "To position telescope correctly we need to offset by"
      write/out "      dx = {x_move}    dy = {y_move}  [1/10 arcsec]"

      !*goto exit		! test without hardware commands
      set/format I5 	! for telescope command

	! offset telescope
      $ $TECS_SCRIPT/t_coord_system xy 
      $ $TECS_SCRIPT/t_offset {x_move} {y_move} -
        | awk '{if(NR==1){print $1}}' | write/keyword tel_return	! pipe return

      if tel_return .ne. 0 then
         write/out "ERROR: Telescope return value for t_offset signals an error..."
         write/out "...the program is aborted"
         $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
         goto exit
      else
         isodate = m$isodate()
         write/file {fctrl(1)} {isodate} {x_move} {y_move} preset
      endif	
      write/out "                  done !"
   endif
endif
   
if pattern_reps .gt. 20 then
  pattern_reps = 1		! start at beginning if start_pos > 400
endif

!*goto exit		! test without hardware commands


!.............. telescope and camera commands now ............

	! write initial image descriptors
$cmd_o2000 counter DITH_NO set {start_position}	! set dither counter to start position
$cmd_o2000 counter POINT_NO set {point_no}	! set pointing no to its value
$cmd_o2000 counter EXPO_NO clear		! clear exposure counter-->EXPO_NO=1


	! set telescope in XY-mode
$ $TECS_SCRIPT/t_coord_system xy

	! set single image parameters
set/format I1
$cmd_o2000 crep {rep_integrate}
$cmd_o2000 itime {single_time}
$cmd_o2000 sync

!
!
!-------------------------------------------------------------
	! do loop over number of images
do loop = 1 {rep_image} 1

set/format I1
write/out "Taking image {loop} of {rep_image}..."	

	! write object
$cmd_o2000 object {P3}:{loop}/{rep_image}


set/format I5 	! for telescope command

	! offset telescope
$ $TECS_SCRIPT/t_offset {x_offset({counter})} {y_offset({counter})} -
| awk '{if(NR==1){print $1}}' | write/keyword tel_return	! pipe return

if tel_return .ne. 0 then
  write/out "ERROR: Telescope return value for t_offset signals an error..."
  write/out "...the program is aborted"
  $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
  goto exit
else
  isodate = m$isodate()
  write/file {fctrl(1)} {isodate} {x_offset({counter})} {y_offset({counter})} {loop} obj
endif	


$cmd_o2000 read
$cmd_o2000 sync


	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
if abort_check .eq. 1 then
  write/out "Program is aborted..."
  $rm /disk-a/o2k/tmp/geirsLstAbort
  $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
  goto exit
endif


	! add file to image catalog
if loop .ge. 2 then
  set/midas output=logonly
  $cmd_o2000 last	! writes last filename in file geirsLstFile
	! writes last filename in keyword pathname_ima
  write/keyword pathname_ima </disk-a/o2k/tmp/geirsLstFile 
	! add file to icat
  add/icat {icatalog} {pathname_ima}
  set/midas output=yes
endif


$cmd_o2000 save -i


	! handle image descriptors
$cmd_o2000 counter DITH_NO incr			! increment dither counter by 1
$cmd_o2000 counter EXPO_NO clear		! reset exposure counter 


	! check wether 20 images were taken
if counter .eq. 20 then

	set/format I5
		! set telescope back to first position of last dither pattern
		! calculated for integer pixel offsets
	$ $TECS_SCRIPT/t_offset {X_back} {Y_back} -
	 | awk '{if(NR==1){print $1}}' | write/keyword tel_return

	if tel_return .ne. 0 then
  	  write/out "ERROR: Telescope return value for t_offset signals an error..."
  	  write/out "...the program is aborted"
        $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
  	  goto exit
        else
          isodate = m$isodate()
          write/file {fctrl(1)} {isodate} {X_back} {Y_back} {loop} back_to_start
	endif	
	
		! set telescope to next starting position
	$ $TECS_SCRIPT/t_offset {x_repetition({pattern_reps})} {y_repetition({pattern_reps})} -
	| awk '{if(NR==1){print $1}}' | write/keyword tel_return

	if tel_return .ne. 0 then
  	  write/out "ERROR: Telescope return value for t_offset signals an error..."
  	  write/out "...the program is aborted"
        $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
  	  goto exit
        else
          isodate = m$isodate()
          write/file {fctrl(1)} {isodate} {x_repetition({pattern_reps})} {y_repetition({pattern_reps})} {loop} {pattern_reps} next_rep
	endif		


		! after 400 images, start with first one
	if pattern_reps .eq. 20 then
	  pattern_reps = 0
	endif


	pattern_reps = pattern_reps+1
	counter = 0

endif


	! increment counter
counter = counter+1


enddo
!-------------------------------------------------------------

$cmd_o2000 sync    ! wait for last save

set/midas output=logonly

	! add last file to image catalog
$cmd_o2000 last	! writes last filename in file geirsLstFile
	! writes last filename in keyword pathname_ima
write/keyword pathname_ima </disk-a/o2k/tmp/geirsLstFile 
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
$auplay /disk-a/staff/GEIRS/SOUNDS/gong.au

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

