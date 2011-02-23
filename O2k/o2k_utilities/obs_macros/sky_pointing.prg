!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! .COPYRIGHT	(C) MPIA
! .IDENT	sky_pointing.prg
!
! .CALL		@@ sky_pointing icatalog itime=[tot_pointing,tot_single,single,tot_single_sky]
!				move  direction  object  start_pos,tel_pos  pointing offsets
! 
!.PARAMETERS	P1=icatatlog : name of image catalog for images
!		P2=itime : intergration time in seconds:
!			tot_pointing,tot_single,single,tot_single_sky
!			tot_pointing = total integration time for pointing
!			tot_single = total integration time for one image
!			single = integration time of a single readout
!			tot_single_sky = total integration time for sky exposure [default=tot_single]
!		P3=move : arcmin from object to move telescope to 
!		P4=direction : direction for sky --> N,S,W,E or ALL
!			     all will use 8 different direction at "move distance"
!		P5:object = object name for identifier
!		P6:start_pos,tel_pos
!                       start_pos = continue dither pattern at this position (N)
!			tel_pos   = AQ   telescope is at position 1 (aequisition)
!                                 = PREV telescope is at position N-1  (default)
!		P7:pointing = pointing identication number for descriptor POINT_NO
!		P8:offsets = flag for integer or non integer pixel offsets
!			     0 = integer pixel offsets ; 1 = non-integer offsets
!
! .PURPOSE	A procedure to take alternately dithered images and sky 
!		with OMEGA2000
! .ENVIRONMENT	MIDAS
! .AUTHOR	Rene Fassbender     
! .KEYWORDS	dithering, pointing, sky
! .COMMENTS	An image catalog is opened, which will contain all images taken.
!		The procedure will use the image names, filter etc specified in the
!		GUI. At the current observing position, a dithered image pattern
!		for the specified total integration time will be taken.
!		After each saved image (tot_single secs) the telescope is moved
!		to the designated sky position. A sky image is taken with the
!		designated integration time and the same dither pattern. A single
!		dither pattern contains 20 observing positions, after that the
!		pattern will be repeated with a slight offset for the first image.
!		The area with no overlap is maximally 26 arcsec on each side.
!               A log-file will be created with all telescope offsets issued.
!
! .CONVENTIONS	!! --> command of line necessary
!		!* --> command not needed but could be useful for testing
!
! .VERSION	1.0	16.12.02
!		2.0	16.01.03	adjusted to CA o2k system
!		3.0	18.01.03	implementation of features P5, P6
!		4.0	12.03.03	implemented telescope return value check,
!					abort check, write descriptors, integer
!					pixel offsets
!		4.01	18.03.03	corrected some typos and missing ()
!		4.1	10.04.03	write file for catalog name only,
!					changed integer repetition offsets,
!					added key file_path, changed awk to 1.line
!		4.2	27.05.03	add separate integration time for sky
!               4.3      9.07.04        adjust repetition pattern and drizzle-
!                                       compatible offset to dither_pointing.
!                                       add 1 / N-1  option for starting position
!               4.4     12.07.04        telescope log file added.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!!!!!!!!!!!!!!!!!! SETUP AND HELP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! command line parameter setup
crossref icatalog itime move direction object start_pos pointing offsets


	! help text
if p1(1:4) .eq. "help" then
   write/out
   write/out "sky_pointing.prg"
   write/out "call:@@ sky_pointing icatalog itime=[tot_pointing,tot_single,single,tot_single_sky]"
   write/out "		        [move] [direction] object start_pos pointing"
   write/out 
   write/out "The command line parameters are:"
   write/out
   write/out "P1=icatatlog = name of image catalog"
   write/out "P2=itime=intergration time in seconds:"
   write/out "	tot_pointing,tot_single,single,tot_single_sky"
   write/out "	tot_pointing = total integration time for pointing"
   write/out "	tot_single = total integration time for one image"
   write/out "	single = integration time of a single readout"
   write/out "	tot_single_sky = total integration time for a single sky image [default=tot_single]"
   write/out "P3=move : arcmin from object to move telescope to "
   write/out "P4=direction : direction for sky --> N,S,W,E or ALL"
   write/out "		  all will use 8 different direction at move_distance"
   write/out "P5:object = object name for identifier"
   write/out "P6:start/telescope position = start_pos,tel_pos"
   write/out "            start_pos = continue dither pattern at this position (#n, def=1)"
   write/out "            tel_pos   = flag to indicate telescope position:"
   write/out "                          = AQ   telescope is at dither position 1 (aequisition)"
   write/out "                          = PREV telescope is at position n-1 (default)"
   write/out "P7:pointing = pointing identication number for descriptor POINT_NO"
   write/out "P8:offsets  = flag for integer or non integer pixel offsets"
   write/out "		     0 = integer pixel offsets ; 1 = non-integer offsets"
   write/out
   write/out "This program takes alternately dithered images and sky frames for one "
   write/out "pointing for OMEGA2000:"
   write/out "The area with no overlap is maximally 26 arcsec on each side"
   write/out
   return		! back to MIDAS session
endif	

define/parameter P1 ? C "Enter a name for the image catalog:"
define/parameter P2 ? N "Enter the integration times in secs[tot,integ,single,tot_single_sky]:"
define/parameter P3 30 N "Enter move parameter in arcmin within limits" 0,150
define/parameter P4 A C "direction"
define/parameter P5 ? C "Enter an object name:"
define/parameter P6 1,PREV C "Sequence number of start position :"
define/parameter P7 1 N "Enter an identification number for the descriptor POINT_NO:"
define/parameter P8 0 N "offsets"


!
	! define keyword which contains the path where pipe_files are stored
define/local file_path/c/1/80 "/disk-a/o2k/tmp"
define/local fctrl/i/1/2 0,0
define/local i_p/i/1/1 0
define/local isodate/c/1/32 " "
isodate = m$isodate()
open/file tel_pos_{isodate}.log WRITE fctrl
write/file {fctrl(1)} {isodate} Log of telescope movements for {P5} (pointing {P7})

	! extract catalog name
define/local icatalog/c/1/80 {P1}

	! extract integration times
define/local itime/r/1/4 400,20,2,0		! defaults; sky_int_time = object int time
write/key itime/r/1/4 {P2}

define/local tot_pointing/r/1/1 {itime(1)}
define/local tot_single/r/1/1 {itime(2)}
define/local single_time/r/1/1 {itime(3)}
define/local tot_sky/r/1/1 {itime(4)}		! new sky parameter

define/local sky_check/i/1/1 0			! if set to 0-->sky and image have same int_time,
						! if set to 1--> different int times


	! if tot_single_sky was not specified, use same sky integration time as for object
if itime(4) .gt. 0 then				! time was specified

  if itime(4) .ne. itime(2) then		! if sky and object integration time differ 
    sky_check = 1		
  endif

else						! set tot_single_sky to tot_single_image

  tot_sky = tot_single

endif



	! initialize pointing number
define/local point_no/i/1/1 {P7}

	! get start position

i_p = m$index(P6,",")
if i_p .eq. 0 then
   define/local start_position/i/1/1 {P6}
   define/local tel_flag/c/1/8 "PREV"
else
   define/local start_position/i/1/1 {P6(1:{i_p})}
   i_p = i_p + 1
   define/local tel_flag/c/1/8 {P6({i_p}:)}
endif
tel_flag = m$upper(tel_flag)
if tel_flag(1:2) .ne. "AQ"  then
   if tel_flag(1:4) .ne. "PREV" then
      write/out
      write/out "Input error:"
      write/out "Flag for current telescope position (parameter 6b) has to be AQ or PREV"
      write/out "   ... abort    "
      write/out
      $auplay /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
      goto exit
   endif
endif	

	! for return value of abort-file-check
define/local abort_check/i/1/1 0	

	! define keyword for path and image name
define/local pathname_ima/c/1/200

	! extract sky direction and distance
define/local sky_move/r/1/1 {P3}
define/local dir_flag/i/1/1 0
	! 0=all direction; 1=N ; 2=E ; 3=S ; 4=W


P4 = m$lower(P4)	! change keyword to lower case

if P4(1:1) .eq. "n" then
	dir_flag = 1
elseif P4(1:1) .eq. "e" then
	dir_flag = 2
elseif P4(1:1) .eq. "s" then
	dir_flag = 3	
elseif P4(1:1) .eq. "w" then
	dir_flag = 4
endif
	
	! display input parameters
set/format I1 F4.2

write/out
write/out "A dither pattern will now be used to take images..."	
write/out "All images will be listed in the image catalog: {icatalog}"
write/out "The total integration time for the poining will be:{tot_pointing} secs"
write/out "The total integration time for a single image is: {tot_single} secs"
write/out "The total integration time for a single sky exposure is: {tot_sky} secs"
write/out "The integration time for a single frame is: {single_time} secs" 
write/out "Sky frames will be taken {sky_move} arcminutes from object,"
write/out "at direction {dir_flag}"
write/out "(0=all directions ; 1=N ; 2=E ; 3=S ; 4=W)"
write/out "The identifier will contain the object name: {P5}"
write/out "The dither pattern will be started at position: {start_position}"
write/out "Telescope is assumed to be at {tel_flag} position"
write/out "The descriptor POINT_NO will contain the pointing number:{point_no}"
write/out "The dither offsets will be (0=integer pixels ; 1=non-integer):{P8}"
write/out "Telescope is assumed to be at {tel_flag} position"
write/out "Telescope movements are logged into file tel_pos_{isodate}.log"
write/out

!write/out "itime={itime(1)},{itime(2)},{itime(3)},{itime(4)} ; P2 = {P2}"
!*goto exit	! for test purposes


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

	! for offsets of 1/3 pixel
if {P8} .eq. 1 then


  define/local x_offset/i/1/20 0,101,-150,180,-227,-48,198,128,-375,225
  write/keyword x_offset/i/11/10 -209,150,-90,317,-42,-216,143,-237,36,254

  define/local y_offset/i/1/20 0,101,57,-216,-29,150,-198,182,-87,252
  write/keyword y_offset/i/11/10 -74,-309,162,-98,291,-99,-56,-192,360,-395
  
  define/local X_back/I/1/1 -137
  define/local Y_back/i/1/1 197

  write/out "Fractional pixel offsets for drizzle were defined..."
!..............................................................................

else	! integer pixel offsets	: multiples of 9 = 2 pixels

  define/local x_offset/i/1/20 0,99,-153,180,-225,-45,198,126,-378,225
  write/keyword x_offset/i/11/10 -207,153,-90,315,-45,-216,144,-234,36,252

  define/local y_offset/i/1/20 0,99,54,-216,-27,153,-198,180,-90,252
  write/keyword y_offset/i/11/10 -72,-306,162,-99,288,-99,-54,-189,360,-396

  define/local X_back/I/1/1 -135
  define/local Y_back/i/1/1 198

  write/out "Integer pixel offsets were defined..."

endif

  ! define repetion offsets (always integer pixels) --> modified to yield sum=0
define/local x_repetition/i/1/20 0,-54,72,27,-90,27,54,-72,72,-54
write/keyword x_repetition/i/11/10 27,54,-27,-54,45,-54,72,-99,99,-45

define/local y_repetition/i/1/20 54,45,-99,99,-72,72,-54,27,27,-72
write/keyword y_repetition/i/11/10 45,-72,72,-72,27,54,-54,27,27,-81


	! calculate repitions
define/local rep_image/i/1/1		! total number of images to be taken
rep_image = m$nint(tot_pointing/tot_single)

define/local rep_integrate/i/1/1	!number of images for single-image-integra
rep_integrate = m$nint(tot_single/single_time)


	! for sky integration
define/local rep_int_sky/i/1/1		! number of images for a total single sky exposure
rep_int_sky = m$nint(tot_sky/single_time)

write/out "For the object a total number of {rep_image} dither positions with"
write/out "                 {rep_integrate} exposures per position will be taken."
write/out "For the sky      {rep_int_sky} exposures per position will be taken."
write/out

	! calculate sky offsets in units of 1/10 arsecs
define/local alpha_sky/I/1/1 0
define/local delta_sky/I/1/1 0

if dir_flag .eq. 1 then
	delta_sky = 600*({sky_move})
elseif dir_flag .eq. 2 then
	alpha_sky = 600*({sky_move})
elseif dir_flag .eq. 3 then
	delta_sky = -600*({sky_move})	
elseif dir_flag .eq. 4 then
	alpha_sky = -600*({sky_move})
else 
	alpha_sky = 600*({sky_move})
	delta_sky = 600*({sky_move})
endif


	! dummy variables for possible --alpha_sky entries
define/local m_alpha_sky/I/1/1					! contains the minus sign
m_alpha_sky = -({alpha_sky})

define/local m_delta_sky/I/1/1
m_delta_sky = -({delta_sky})

	! define vectors for sky positions
define/local alpha_all/i/1/8 {alpha_sky},{alpha_sky},{alpha_sky},0,{m_alpha_sky},{m_alpha_sky},{m_alpha_sky},0

define/local delta_all/i/1/8 {delta_sky},0,{m_delta_sky},{m_delta_sky},{m_delta_sky},0,{delta_sky},{delta_sky}
!
!!!!!!!!!!!!!!!!! TELESCOPE MACRO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! keywords
define/local loop/i/1/1 
define/local counter/i/1/1 1		! for counting to 20
define/local pattern_reps/i/1/1 2	! for counting how many times 20
define/local move_counter/i/1/1 1	! for sky positions
define/local tel_return/i/1/1 0		! for telescope return values:0,-1

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
  pattern_reps = 1		! start at beginning if start_pos>400
endif

!*goto exit		! test without hardware commands


!''''''''''''''' telescope and camera commands now '''''''''''''

	! write initial image descriptors
$cmd_o2000 counter DITH_NO set {start_position}	! set dither counter to start position
$cmd_o2000 counter POINT_NO set {point_no}	! set pointing no to its value
$cmd_o2000 counter EXPO_NO clear		! --> EXPO_NO = 1


	! set telescope in XY-mode
$ $TECS_SCRIPT/t_coord_system xy 


	! set single image parameters
set/format I1
$cmd_o2000 crep {rep_integrate}
$cmd_o2000 itime {single_time}
$cmd_o2000 sync

!
!-------------------------------------------------------------
	! do loop over number of images
do loop = 1 {rep_image} 1


  	! if sky has different integration time ---> set back
if sky_check .eq. 1 then
	! set single image parameters
  set/format I1
  $cmd_o2000 crep {rep_integrate}
  $cmd_o2000 sync
endif

set/format I1
write/out "Taking image {loop} of {rep_image}..."	

	! write object
$cmd_o2000 object {P5}: {loop}/{rep_image}


set/format I5 	! for telescope command

	! offset telescope
$ $TECS_SCRIPT/t_offset {x_offset({counter})} {y_offset({counter})} -
| awk '{if(NR==1){print $1}}' | write/keyword tel_return	

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

!*$cmd_o2000 sync	! test, not needed


!..................................................
	! take sky image

  	! if sky has different integration time ---> set to sky paremeters
if sky_check .eq. 1 then
	! set single image parameters
  set/format I1
  $cmd_o2000 crep {rep_int_sky}
endif



set/format I1
write/out "       sky"

	! write object name
$cmd_o2000 object sky for {P5}:{loop}/{rep_image}


set/format I5 	! for telescope command

	! first move telescope to sky position
$ $TECS_SCRIPT/t_offset {alpha_sky} {delta_sky} -
| awk '{if(NR==1){print $1}}' | write/keyword tel_return

if tel_return .ne. 0 then
  write/out "ERROR: Telescope return value for t_offset signals an error..."
  write/out "...the program is aborted"
  $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
  goto exit
else
  isodate = m$isodate()
  write/file {fctrl(1)} {isodate} {alpha_sky} {delta_sky} {loop} move_to_sky
endif


	! take sky
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
set/midas output=logonly
$cmd_o2000 last	! writes last filename in file geirsLstFile
	! writes last filename in keyword pathname_ima
write/keyword pathname_ima </disk-a/o2k/tmp/geirsLstFile 
	! add file to icat
add/icat {icatalog} {pathname_ima}
set/midas output=yes


$cmd_o2000 save -i


	! move telescope back to object position
alpha_sky = -1*alpha_sky
delta_sky = -1*delta_sky

$ $TECS_SCRIPT/t_offset {alpha_sky} {delta_sky} -
| awk '{if(NR==1){print $1}}' | write/keyword tel_return

if tel_return .ne. 0 then
  write/out "ERROR: Telescope return value for t_offset signals an error..."
  write/out "...the program is aborted"
  $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
  goto exit
else
  isodate = m$isodate()
  write/file {fctrl(1)} {isodate} {alpha_sky} {delta_sky} {loop} back_to_obj
endif


	! set offsets back for next sky
alpha_sky = -1*alpha_sky
delta_sky = -1*delta_sky


	! if all directions for sky are used, get next sky direction
if dir_flag .eq. 0 then

	move_counter = move_counter+1

	if move_counter .eq. 9 then 	! set offsets back to first
	move_counter = 1 
	endif

	alpha_sky = alpha_all({move_counter})
	delta_sky = delta_all({move_counter})
	
endif

!..................................................

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


	! handle image descriptors
	! object EXPO_NO = 1; sky EXPO_NO = 2 ; DITH_NO is the same
$cmd_o2000 counter DITH_NO incr			! increment dither counter by 1
$cmd_o2000 counter EXPO_NO clear		! reset exposure counter 


	! increment counter
counter = counter+1


enddo
!-------------------------------------------------------------


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
write/out "All images for pointing finished..." 
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

return
