!	d e t _ t i l t . p r g		22-Jul-05	HJR
!
!  Series of exposures with telescope not tracking to determine detector tilt
!
!  Filter and read-out mode have to be set up manually before call!
!

crossref ident eq exp_time coadds number intervall

define/par P1 ? ? "Identifier : "
if P1(1:4) .eq. "help" then
   write/out
   write/out "darks.prg"
   write/out "call: o2k/det_tilt ident = [exp_time] [coadds] [n_exp]  time_intervall" 
   write/out 
   write/out "The command line parameters are:"
   write/out "  P1 = Identifier (enclose in quotes if blank is to be used!)"
   write/out "  P2 = "
   write/out "  P3 = exposure time"
   write/out "  P4 = number of coadds, i.e. number of exposures to be added up in memory"
   write/out "  P5 = number of images to be taken"
   write/out "  P6 = time intervall between exposures"
   write/out
   write/out "To abort a running sequence press ABORT in the GEIRS gui! Never use ^C"
   write/out
   goto exit
endif
define/par P2 =
define/par P3 ? N/A "Exposure time [sec] : "
define/par P4 1 N/A "Number of coadds in memory : "
define/par P5 5 N/A "Number of images : "
define/par P6 10 ? "Time intervall : "
define/par P7 i ? "Save mode [single or integrated] : "

define/local i/i/1/1 0
define/local n_exp/i/1/1 {P5}
define/local coadds/i/1/1 {P4}
define/local exp_time/r/1/1 {P3}
define/local wait_time/r/1/1 0
define/local time_int/r/1/1 {P6}

P7 = m$lower(P7)

wait_time = time_int - coadds * exp_time

if wait_time .le. 0 then
  write/out "Time intervall too short for combination of exp_time and coadds!"
  $play -q /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
  goto exit
endif

define/local abort_check/i/1/1 0	! for return value of abort-file-check
	! remove existing abort file
	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("{geirslstabort}")
if abort_check .eq. 1 then
  $rm {geirslstabort} 
endif

set/format I1 F5.1

! setup camera

$cmd_panic_new object {P1}

$cmd_panic_new counter DITH_NO clear	! clear dither counter 
$cmd_panic_new counter POINT_NO clear	! clear pointing no 
$cmd_panic_new counter EXPO_NO clear	! clear exposure counter-->EXPO_NO=1

$cmd_panic_new itime {exp_time}
$cmd_panic_new crep {coadds}
$cmd_panic_new sync


do i = 1 {n_exp}
   write/out "         Exposing step {i} with {exp_time}sec ..."
   $cmd_panic_new read
   $cmd_panic_new sync
		! abort check: 0=does not exist ; 1=abort file exists
	abort_check = m$exist("{geirslstabort}")
	if abort_check .eq. 1 then
  	  write/out "         Program is aborted..."	
	  $rm {geirslstabort} 	! remove file again
        $play -q /disk-a/staff/GEIRS/SOUNDS/crash.au
  	  goto exit
	endif

   if coadds .eq. 1 then
      $cmd_panic_new save
   else
      if P7(1:1) .eq. "i" then
         $cmd_panic_new save -i
      else
         $cmd_panic_new save
      endif
   endif

   wait/sec {wait_time}
enddo

$cmd_panic_new sync
write/out "         det_tilt series finished ..."
$play -q /disk-a/staff/GEIRS/SOUNDS/gong.au

exit:
return
