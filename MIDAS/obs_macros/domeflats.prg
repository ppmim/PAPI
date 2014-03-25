!     d o m e f l a t s . p r g                29-Jul-04    HJR
!
!   Take domeflats with flatfield lamps
!   Lamp is to be specified in call
!

crossref ident eq lamp time coadds number save saturation

define/par P1 ? ? "Identifier : "
if P1(1:4) .eq. "help"  then
   write/out
   write/out "Purpose: Take two series of dome flats,"
   write/out "         first with lamp on then with lamp off"
   write/out
   write/out "P1 identifier [no default]"
   write/out "P2 ="
   write/out "P3 Lamp(,level) to select the lamp [no default]"
   write/out "         and in case of lamp 5 also the intesitiy level in Watt"
   write/out "P4 Exposure time of a single read in sec [no default]"
   write/out "P5 Number images to be of added in memory [1]"
   write/out "P6 Number of read cycles [1]"
   write/out "P7 Save option, i.e. save individual images (# images = # of coadds"
   write/out "                or integrated images (one file per read-cycle) [single]"
   write/out "P8 Saturation level in counts [30000]"
   write/out
   goto exit
endif

define/par P2 =
define/par P3 ? N/A "Lamp [1...5], level [1...10]: "
define/par P4 ? N/A "Integration time per exposure : "
define/par P5 1 N/A "Number of coadds : "
define/par P6 1 N/A "Number of read cycles : "
define/par P7 i C/A "Save option [single/integrated] : "
define/par P8 30000 N/A "Saturation level : "

define/local lamp/i/1/2 {P3}		! lamp_ID,level
define/local exp_time/r/1/1 {P4}
define/local coadds/i/1/1 {P5}
define/local n_ima/i/1/1 {P6}
P7 = m$lower(P7)
define/local save/c/1/8 {P7}
define/local i/i/1/1 0
define/local act_level/r/1/1 0
define/local curr_level/r/1/2 0,0
define/local diff_level/r/1/1 0
define/local saturation/r/1/1 {P8}

define/local geirslstabort/c/1/256  " "
geirslstabort = M$SYMBOL("GEIRSLSTABORT")

lamp(2) = lamp(2) - 1

write/out "GEIRSLSTABORT = {geirslstabort}"

define/local abort_check/i/1/1 0	! for return value of abort-file-check
	! remove existing abort file
	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("{geirslstabort}")


if abort_check .eq. 1 then
  $rm {geirslstabort}
endif

$ $CA_FFLTEST | write/key i

if i .eq. 1 then
   write/out "GUI for flatfield lamps is running. Please shut it down to run this macro."
   $play -q $GEIRS_DIR/SOUNDS/sorrydave.au
   goto exit
endif
if i .eq. 2 then
   write/out "Error message from flatfield lamps. Please contact staff."
   $play -q $GEIRS_DIR/SOUNDS/sorrydave.au
   goto exit
endif


if coadds .eq. 1  then
   save = "single"
endif

set/format i1

if lamp(1) .eq. 5  then
   $ $FLAT_COMMAND L{lamp(1)} on {lamp(2)}
   lamp(2) = lamp(2) + 1
   write/out "         Using lamp 5 with {lamp(2)}W"
else
   $ $FLAT_COMMAND L{lamp(1)} on
   write/out "         Using lamp {lamp(1)}"
endif

wait/sec 2

$cmd_panic_new itime {exp_time}
$cmd_panic_new crep {coadds}
$cmd_panic_new object {P1} lamp on
$cmd_panic_new sync

do i = 1 {n_ima}
   write/out "         Taking image {i} with lamp on ..."
   $cmd_panic_new read
   $cmd_panic_new sync
		! abort check: 0=does not exist ; 1=abort file exists
	abort_check = m$exist("{geirslstabort}")
	if abort_check .eq. 1 then
  	  write/out "         Program is aborted..."	
	  $rm {geirslstabort} 	! remove file again
        $play -q $GEIRS_DIR/SOUNDS/crash.au
  	  goto exit
	endif

   if i .eq. 1  then
      $cmd_panic_new median -stdout -raw | awk '{if(NR==1){print $1}}'| write/key curr_level/r/1/2
      act_level = curr_level(2)-curr_level(1)
      diff_level = act_level
      if curr_level(1) .gt. 5000  then
         write/out
         write/out "         Image saturated. Sequence aborted !"
         write/out "         Level = {curr_level(2)}!"
         write/out
         $play -q $GEIRS_DIR/SOUNDS/crash.au
         goto exit
      endif
      if act_level .gt. saturation  then
         write/out
         write/out "         Image saturated. Sequence aborted !"
         write/out "         Level = {diff_level}!"
         write/out
         $play -q $GEIRS_DIR/SOUNDS/crash.au
         goto exit
      else
         write/out "         Level of images with lamp on = {act_level} counts"
      endif
   endif
   if save(1:1) .eq. "s" then
      $cmd_panic_new save
   else
      $cmd_panic_new save -i
   endif
enddo

$ $FLAT_COMMAND ALLOFF
wait/sec 10

$cmd_panic_new object {P1} lamp off
$cmd_panic_new sync

do i = 1 {n_ima}
   write/out "         Taking image {i} with lamp off ..."
   $cmd_panic_new read
   $cmd_panic_new sync

		! abort check: 0=does not exist ; 1=abort file exists
	abort_check = m$exist("{geirslstabort}")
	if abort_check .eq. 1 then
  	  write/out "         Program is aborted..."	
	  $rm {geirslstabort}	! remove file again
        $play -q $GEIRS_DIR/SOUNDS/crash.au
  	  goto exit
	endif


   if save(1:1) .eq. "s" then
      $cmd_panic_new save
   else
      $cmd_panic_new save -i
   endif
   if i .eq. 1  then
      $cmd_panic_new median -stdout -raw | awk '{if(NR==1){print $1}}'| write/key curr_level/r/1/2
      act_level = curr_level(2)-curr_level(1)
      diff_level = diff_level - act_level
      write/out "         Level of images with lamp off = {act_level} counts"
      write/out "         Level difference (on - off) = {diff_level} counts"
   endif
enddo
$cmd_panic_new sync
write/out "         done ..."
$play -q $GEIRS_DIR/SOUNDS/gong.au

exit:
return
