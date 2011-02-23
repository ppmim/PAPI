!	c a l i b r a t i o _ s e r i e s . p r g		28-Jul-04	HJR
!
!  Series of exposures for calibration (flat, dark)
!
!  Filter and read-out mode have to be set up manually before call!
!
!  25.9.04:   reset introduced to deal (hopefully) with fifo overflow problem
!  2.10.04:   min_time introduced

crossref ident eq time coadds number spacing save reset

define/par P1 ? ? "Identifier : "
if P1(1:4) .eq. "help" then
   write/out
   write/out "darks.prg"
   write/out "call: o2k/calser ident = [max_<>t] [coadds] [n_exp] [lin/log] [save_mode] [reset]" 
   write/out 
   write/out "The command line parameters are:"
   write/out "  P1 = Identifier (enclose in quotes if blank is to be used!)"
   write/out "  P2 = "
   write/out "  P3 = maximum  or  minumum,maximum exposure time"
   write/out "       If only maximum is given, minimum is taken as zero"
   write/out "       (If time < shortest possible time, shortest time is taken"
   write/out "  P4 = number of coadds, i.e. number of exposures to be added up in memory"
   write/out "  P5 = number of images to be taken"
   write/out "  P6 = linear, logarithmic or zero increase of exposure times [lin]"
   write/out "  P7 = save mode, single images or integrated images"
   write/out "  P8 = reset (def.) / no_reset : set exptiem to 1 sec and read detector twice"
   write/out "                                before increasing exposure time"
   write/out
   write/out "To abort a running sequence press ABORT in the GEIRS gui! Never use ^C"
   write/out
   goto exit
endif
define/par P2 =
define/par P3 ? N/A "Minimun, maximum exposure time [sec] : "
define/par P4 1 N/A "Number of coadds in memory : "
define/par P5 5 N/A "Number of images : "
define/par P6 lin ? "Linear, logarithmic or zero increase of intervalls : "
define/par P7 s C/A "Save mode [single/integrated] : "
define/par P8 no_reset C/A "Reset integration time : "

P6 = m$upper(P6)
P7 = m$lower(P7)
P8 = m$lower(P8)


define/local i/i/1/1 0
define/local n_exp/i/1/1 {P5}
define/local times/r/1/2 {P3}
define/local max_time/r/1/1 0.
define/local min_time/r/1/1 0.
define/local coadds/i/1/1 {P4}
define/local dt/r/1/1 0.
define/local exp_time/r/1/1 0.
define/local true_time/r/1/1 0.
define/local lexp_time/r/1/1 0.

i = m$index(P3,",")
if i .eq. 0 then
   min_time = 0.
   max_time = times(1)
else
  min_time = times(1)
  max_time = times(2)
endif

define/local abort_check/i/1/1 0	! for return value of abort-file-check
	! remove existing abort file
	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
if abort_check .eq. 1 then
  $rm /disk-a/o2k/tmp/geirsLstAbort 
endif

if P6(1:3) .eq. "LOG" then
   dt = m$log10(max_time-min_time)/(n_exp-1)
   exp_time = min_time
else
   if P6(1:4) .eq. "ZERO"  then
      dt = 0
      exp_time = max_time
   else
      dt = (max_time-min_time)/(n_exp-1)
      exp_time = min_time
   endif
endif

set/format I1 F5.1

! setup camera

$cmd_o2000 object {P1}

$cmd_o2000 counter DITH_NO clear	! clear dither counter 
$cmd_o2000 counter POINT_NO clear	! clear pointing no 
$cmd_o2000 counter EXPO_NO clear	! clear exposure counter-->EXPO_NO=1

do i = 1 {n_exp}
   if P8(1:3) .ne. "no_"  then
	write/out "            now resetting exposure time to avoid fifo-overflow ..."
	write/out "            next 2 images are dummies and will not be saved!"
      $cmd_o2000 itime 0
	$cmd_o2000 crep 1
      $cmd_o2000 sync
      $cmd_o2000 read
      $cmd_o2000 sync
      $cmd_o2000 read
      $cmd_o2000 sync
   endif

   true_time = m$nint(exp_time*10.)/10.

   $cmd_o2000 itime {true_time}
   $cmd_o2000 crep {coadds}
   $cmd_o2000 sync
   write/out "         Exposing step {i} with {exp_time}sec ..."
   $cmd_o2000 read
   $cmd_o2000 sync
		! abort check: 0=does not exist ; 1=abort file exists
	abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
	if abort_check .eq. 1 then
  	  write/out "         Program is aborted..."	
	  $rm /disk-a/o2k/tmp/geirsLstAbort 	! remove file again
        $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
  	  goto exit
	endif

   if coadds .eq. 1 then
      $cmd_o2000 save
   else
      if P7(1:1) .eq. "i" then
         $cmd_o2000 save -i
      else
         $cmd_o2000 save
      endif
   endif

   if P6(1:3) .eq. "LOG"  then
      lexp_time = lexp_time + dt
      exp_time = min_time + 10**lexp_time
   else
      if P6(1:4) .eq. "ZERO"  then
         exp_time = max_time
	else
         exp_time = exp_time + dt
      endif
   endif
enddo

$cmd_o2000 sync
write/out "         calibration series finished ..."
$auplay /disk-a/staff/GEIRS/SOUNDS/gong.au

exit:
return
