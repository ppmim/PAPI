!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! .COPYRIGHT	(C) MPIA
! .IDENT	grid_pointing.prg
!
! .CALL		@@ grid_pointing icatalog itime=[tot_pointing,tot_single,single]
!				   object start_pos pointing grid step
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
!               P6:grid = number of grid cells (only 4, 9 or 16)"
!               P7:step = step in arcsec between the centers of the cells (def=60)"
!
! .PURPOSE	A procedure to take dithered images for one pointing using a 
!		grid pattern.
! .ENVIRONMENT	MIDAS
! .AUTHOR	Rene Fassbender  (modified for grid pattern, 21.2.08 U. Thiele)   
! .KEYWORDS	dithering, pointing
! .COMMENTS	An image catalog is opened, which will contain all images taken.
!		The procedure will use the image names, filter etc specified in the
!		GUI. At the current observing position, a dithered image pattern
!		for the specified total integration time will be taken. A single
!		dither pattern contains 4, 9 or 16 observing positions, after that the
!		pattern will be repeated. 
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
!               6.0     01.03.11        Adapted to PANIC by jmiguel@iaa.es
!
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!
!!!!!!!!!!!!!!!!!! SETUP AND HELP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! command line parameter setup
crossref icatalog itime object start_pos pointing offsets


	! help text
if p1(1:4) .eq. "help" then

write/out
write/out "grid_pointing.prg"
write/out "call: @@ grid_pointing icatalog itime=[tot_pointing,tot_single,single]"
write/out "			     start_pos pointing	"
write/out
write/out "This program takes dithered images following a grid pattern"
write/out "The number of grid cells and the step can be choosen"
write/out "Default is a grid of 2x2 with a step of 60arcsec"
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
write/out "P6:grid = number of grid cells (only 4, 9 or 16)"
write/out "P7:step = step in arcsec between the centers of the cells (def=60)"

return		! back to MIDAS session

endif	


	! define command line parameters
define/parameter P1 ? C "Enter a name for the image catalog :"
define/parameter P2 ? N "Enter the integration times in secs[tot,integ,single] :"
define/parameter P3 ? C "Enter an object name :"
define/parameter P4 1,PREV C "Enter the sequence number n of start position :"
define/parameter P5 1 N "Enter the identification number of the current pointing :"
define/parameter P6 4 N "Enter the number of grid cells (4,9,16):"
define/parameter P7 60 N "Enter the step between the grid cells (in arcsec):"

        ! common definitions 
define/local x_off/i/1/1
define/local y_off/i/1/1
define/local nn/i/1/1
define/local notallowed/i/1/1

	! define keyword which contains the path where pipe_files are stored
define/local file_path/c/1/80 "/home/panic/tmp"  ! PENDIENTE
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
      $play -q $GEIRS_DIR/SOUNDS/sorrydave.au
      goto exit
   endif
endif	

        ! get number of grid cells
define/local grid/i/1/1 {P6}

if grid .ne. 4 .and. grid .ne. 9 then
   notallowed = 1
endif
if grid .ne. 16 .and. notallowed .eq. 1 then
      write/out " "
      write/out "only 4, 9 or 16 allowed!  ...abort "
      write/out " "
      $play -q $GEIRS_DIR/SOUNDS/sorrydave.au
      goto exit
endif

        ! get step between cells
define/local step/i/1/1 {P7}


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
write/out "Telescope is assumed to be at {tel_flag} position"
write/out "Telescope movements are logged into file tel_pos_{isodate}.log"
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
	
	! define dither offsets in units of 1, later multiplied by step factor
	! all positions are defined relative to last position
set/format I1

if {grid} .eq. 4 then

  define/local x_off1/r/1/1/ 0.5
  define/local y_off1/r/1/1/ 0.5

  define/local x_offset/r/1/4 0,-1,0,1
  define/local y_offset/r/1/4 0,0,-1,0
 
  define/local X_back/r/1/1 0
  define/local Y_back/r/1/1 1

endif

if {grid} .eq. 9 then

  define/local x_off1/r/1/1/ 0
  define/local y_off1/r/1/1/ 0

  define/local x_offset/r/1/9 0,1,0,-1,-1,0,0,1,1
  define/local y_offset/r/1/9 0,0,1,0,0,-1,-1,0,0
 
  define/local X_back/r/1/1 -1
  define/local Y_back/r/1/1 1

endif

if {grid} .eq. 16 then

  define/local x_off1/r/1/1/ 0.5
  define/local y_off1/r/1/1/ 0.5

  define/local x_offset/r/1/16 0,-1,0,1,1,0,0,-1,-1,-1,0,0,0,1,1,1
  define/local y_offset/r/1/16 0,0,-1,0,0,1,1,0,0,0,-1,-1,-1,0,0,0
 
  define/local X_back/r/1/1 -1
  define/local Y_back/r/1/1 2

endif

!..............................................................................


  ! define repetition offsets (always integer pixels) --> modified to yield sum=0
define/local x_repetition/r/1/20 0 all

define/local y_repetition/i/1/20 0 all



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
define/local counter/i/1/1 1		! for counting to {grid}
define/local pattern_reps/i/1/1 2	! for counting how many times {grid}
define/local tel_return/i/1/1 0		! character keyword for telescop returns:0,-1


	! offset telescope to start position
      x_off = x_off1 * step * 10
      y_off = y_off1 * step * 10
      $ {tecs_script}/t_coord_system xy 
      $ {tecs_script}/t_offset {x_off} {y_off} -
        | awk '{if(NR==1){print $1}}' | write/keyword tel_return	! pipe return

      if tel_return .ne. 0 then
         write/out "ERROR: Telescope return value for t_offset signals an error..."
         write/out "...the program is aborted"
         $play -q $GEIRS_DIR/SOUNDS/crash.au
         goto exit
      else
         isodate = m$isodate()
         write/file {fctrl(1)} {isodate} {x_off} {y_off} preset
      endif	
      write/out " "
      write/out "Telescope at start position (first grid cell)"
      write/out " "


 	! find the proper start position
do nn = 1 100
  if start_position .gt. grid  then
     start_position = start_position - grid
     if start_position .eq. 0  then
        start_position = 1
     endif
  endif
enddo

	counter = start_position



!!!!!!!!!!!!!!!!!!  move telescope if not at PREV position  !!!!!!!!!!!!!!!!

if tel_flag(1:2) .eq. "AQ"  then
   set/format I3
   define/local x_move/i/1/1 0
   define/local y_move/i/1/1 0
   define/local x_mov/i/1/1 0
   define/local y_mov/i/1/1 0
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
      x_mov = x_move * step * 10
      y_mov = y_move * step * 10
      write/out "To position telescope correctly we need to offset by"
      write/out "      dx = {x_mov}    dy = {y_mov}  [1/10 arcsec]"

      !*goto exit		! test without hardware commands
      set/format I5 	! for telescope command

	! offset telescope
      $ {tecs_script}/t_coord_system xy 
      $ {tecs_script}/t_offset {x_mov} {y_mov} -
        | awk '{if(NR==1){print $1}}' | write/keyword tel_return	! pipe return

      if tel_return .ne. 0 then
         write/out "ERROR: Telescope return value for t_offset signals an error..."
         write/out "...the program is aborted"
         $play -q $GEIRS_DIR/SOUNDS/crash.au
         goto exit
      else
         isodate = m$isodate()
         write/file {fctrl(1)} {isodate} {x_mov} {y_mov} preset
      endif	
      write/out "                  done !"
   endif
endif
   
if pattern_reps .gt. grid then
  pattern_reps = 1		! start at beginning if start_pos > 400
endif

!*goto exit		! test without hardware commands


!.............. telescope and camera commands now ............

	! write initial image descriptors
$cmd_panic_new counter DITH_NO set {start_position}	! set dither counter to start position
$cmd_panic_new counter POINT_NO set {point_no}	! set pointing no to its value
$cmd_panic_new counter EXPO_NO clear		! clear exposure counter-->EXPO_NO=1


	! set telescope in XY-mode
$ {tecs_script}/t_coord_system xy

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


set/format I5 	! for telescope command

	! offset telescope
x_off = x_offset({counter}) * step * 10
y_off = y_offset({counter}) * step * 10
$ {tecs_script}/t_offset {x_off} {y_off} -
| awk '{if(NR==1){print $1}}' | write/keyword tel_return	! pipe return

if tel_return .ne. 0 then
  write/out "ERROR: Telescope return value for t_offset signals an error..."
  write/out "...the program is aborted"
  $play -q $GEIRS_DIR/SOUNDS/crash.au
  goto exit
else
  isodate = m$isodate()
  write/file {fctrl(1)} {isodate} {x_off} {y_off} {loop} obj
endif	


$cmd_panic_new read
$cmd_panic_new sync


	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("{geirslstabort}")
if abort_check .eq. 1 then
  write/out "Program is aborted..."
  $rm {geirslstabort}
  $play -q $GEIRS_DIR/SOUNDS/crash.au
  goto exit
endif


	! add file to image catalog
if loop .ge. 2 then
  set/midas output=logonly
  $cmd_panic_new last | awk '{print $2}' | write/keyword pathname_ima
	! add file to icat
  add/icat {icatalog} {pathname_ima}
  set/midas output=yes
endif


$cmd_panic_new save -i


	! handle image descriptors
$cmd_panic_new counter DITH_NO incr			! increment dither counter by 1
$cmd_panic_new counter EXPO_NO clear		! reset exposure counter 


	! check wether {grid} images were taken
if counter .eq. grid then

	set/format I5
		! set telescope back to first position of last dither pattern
		! calculated for integer pixel offsets
        x_off = X_back * step * 10
        y_off = Y_back * step * 10
 	$ {tecs_script}/t_offset {x_off} {y_off} -
 	 | awk '{if(NR==1){print $1}}' | write/keyword tel_return

 	if tel_return .ne. 0 then
   	  write/out "ERROR: Telescope return value for t_offset signals an error..."
   	  write/out "...the program is aborted"
         $play -q $GEIRS_DIR/SOUNDS/crash.au
   	  goto exit
         else
           isodate = m$isodate()
           write/file {fctrl(1)} {isodate} {X_back} {Y_back} {loop} back_to_start
 	endif	
	
!		! set telescope to next starting position
!	$ {tecs_script}/t_offset {x_repetition({pattern_reps})} {y_repetition({pattern_reps})} -
!	| awk '{if(NR==1){print $1}}' | write/keyword tel_return

!	if tel_return .ne. 0 then
!  	  write/out "ERROR: Telescope return value for t_offset signals an error..."
!  	  write/out "...the program is aborted"
!        $play -q $GEIRS_DIR/SOUNDS/crash.au
!  	  goto exit
!        else
!          isodate = m$isodate()
!          write/file {fctrl(1)} {isodate} {x_repetition({pattern_reps})} {y_repetition({pattern_reps})} {loop} {pattern_reps} next_rep
!	endif		


!		! after 400 images, start with first one
!	if pattern_reps .eq. grid then
!	  pattern_reps = 0
!	endif


!	pattern_reps = pattern_reps+1
	counter = 0

endif


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

