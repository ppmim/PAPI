! rel_cal				HJR      21. September 2003
!
! modified for PANIC:   jmiguel@iaa.es  01-March-11
!
! relative calibration of an PANIC mosaic of 1(PENDIENTE) square degree (4x4 pointings)
! 3x3 pointings in intersection of deep exposures
! Telescope assumed to be positioned in center of 4x4 mosaic, spacing 14.3arcmin

crossref ident eq expose

define/par P1 ? ? "Identifier for images : "
define/par P2 =
define/par P3 ? N/A "Exposure time of single frame, number of exposures : "

if p1(1:4) .eq. "help" then
   write/out
   write/out "rel_cal.prg"
   write/out "call: panic/relcal [object] = [exposure]" 
   write/out 
   write/out "The command line parameters are:"
   write/out "  P1 = identifier for the images taken (keyword OBJECT or IDENT)"
   write/out "  P2 = "
   write/out "  P3 = single frame exposure times, number of exposures"
   write/out
   write/out "To abort a sequence press ABORT in the GEIRS gui! Never use ^C"
   write/out
   write/out "Filter and read-out mode have to be set correctly before call !"
   write/out "Telescope assumed to be in center of 
   write/out "                     area (4x4 pointings) to be calibrated!"
   write/out
   goto ende
endif

write/key outputr/r/1/2 {P3}
define/local dit/r/1/1 0
define/local ndit/i/1/1 0
define/local ans/c/1/32 ""
define/local tel_return/i/1/1 0	! keyword for telescop returns: 0,-1
define/local abort_check/i/1/1 0	! for return value of abort-file-check

dit = outputr(1)
ndit = outputr(2)


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

! telescope on first field
$ {tecs_script}/t_offset -8580 +8580 -
| awk '{if(NR==1){print $1}}' | write/keyword tel_return	! pipe return

if tel_return .lt. 0 then
  write/out "ERROR: Telescope signals error setting to first field ..."
  write/out "...the program is aborted"
  $play -q $GEIRS_DIR/SOUNDS/crash.au
  goto ende
endif	

! set image parameters
$cmd_panic_new crep {ndit}
$cmd_panic_new itime {dit}
$cmd_panic_new sync

define/local ix/i/1/1 0
define/local iy/i/1/1 0

set/format i1

do iy = 1 3
   do ix = 1 3
	write/out "Taking image ({ix},{iy})"
	! write object
	$cmd_panic_new object {P1} relcal {ix},{iy}
	$cmd_panic_new read
	$cmd_panic_new sync
	if ndit .gt. 1 then
		$cmd_panic_new save -i
	else
		$cmd_panic_new save
	endif

	! abort check: 0=does not exist ; 1=abort file exists
	abort_check = m$exist("{geirslstabort}")
	if abort_check .eq. 1 then
	   write/out "Program was aborted from GUI ..."	
	   $rm {geirslstabort} 	! remove file again
         $play -q $GEIRS_DIR/SOUNDS/crash.au
 	   goto ende
	endif
      if ix .lt. 3 then
	   $ {tecs_script}/t_offset 8580 0 -
           | awk '{if(NR==1){print $1}}' | write/keyword tel_return

         if tel_return .lt. 0 then
            write/out "ERROR: Telescope signals error offsetting in X to {ix}+1 ..."
            write/out "...the program is aborted"
            $play -q $GEIRS_DIR/SOUNDS/crash.au
            goto ende
         endif	
	endif
   enddo
   if iy .lt. 3 then
      $ {tecs_script}/t_offset -17160 -8580 -
            | awk '{if(NR==1){print $1}}' | write/keyword tel_return  ! PENDIENTE

      if tel_return .lt. 0 then
         write/out "ERROR: Telescope signals error offsetting in Y to {iy}+1 ..."
         write/out "...the program is aborted"
         $play -q $GEIRS_DIR/SOUNDS/crash.au
         goto ende
      endif	
   else
	$ {tecs_script}/t_offset -8580 +8580 -
            | awk '{if(NR==1){print $1}}' | write/keyword tel_return    ! PENDIENTE

      if tel_return .lt. 0 then
         write/out "ERROR: Telescope signals error going back to origin ..."
         write/out "...the program is aborted"
         $play -q $GEIRS_DIR/SOUNDS/crash.au
         goto ende
      endif	
   endif
enddo

$cmd_panic_new sync

write/out
write/out "         calibration sequence done ..."
write/out
$play -q $GEIRS_DIR/SOUNDS/gong.au

ende:
return
