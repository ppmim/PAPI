!	o f f s e t _ t e s t . p r g		23-Jul-05	HJR
!
!  Series of exposures with telescope being offsetted
!
!  Filter and read-out mode have to be set up manually before call!
!

crossref ident eq exp_time coadds number intervall save

define/par P1 ? ? "Identifier : "
if P1(1:4) .eq. "help" then
   write/out
   write/out "darks.prg"
   write/out "call: o2k/det_tilt ident = [exp_time] [coadds] [n_exp]  time_intervall  save_mode" 
   write/out 
   write/out "The command line parameters are:"
   write/out "  P1 = Identifier (enclose in quotes if blank is to be used!)"
   write/out "  P2 = "
   write/out "  P3 = exposure time"
   write/out "  P4 = number of coadds, i.e. number of exposures to be added up in memory"
   write/out "  P5 = number of images to be taken"
   write/out "  P6 = time intervall between exposures"
   write/out "  P7 = save mode (integrated / single)"
   write/out "  P8 = offset flag"
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
define/par P8 offset ? "Offset flag : "

P7 = m$lower(P7)
P8 = m$lower(P8)

define/local i/i/1/1 0
define/local n_exp/i/1/1 {P5}
define/local coadds/i/1/1 {P4}
define/local exp_time/r/1/1 {P3}
define/local wait_time/r/1/1 0
define/local time_int/r/1/1 {P6}
define/local isodate/c/1/32 " "
define/local tel_return/i/1/1 0
define/local test/i/1/1 0

define/local fctrl/i/1/2 0,0
isodate = m$isodate()
open/file tel_pos_{isodate}.log WRITE fctrl
write/file {fctrl(1)} {isodate} Log of telescope movements for {P1}

define/local x_offset/i/1/20 0 all
define/local y_offset/i/1/20 0 all

if P8(1:3) .eq. "off"  then
   write/keyword x_offset/i/1/10 0,99,-153,180,-225,-45,198,126,-351,198
   write/keyword x_offset/i/11/10 -207,153,-90,288,-18,-216,144,-234,36,252
   write/keyword y_offset/i/1/10 0,99,54,-216,-27,153,-198,180,-90,216
   write/keyword y_offset/i/11/10 -36,-315,171,-99,288,-99,-54,-189,324,-324
endif

$ {tecs_script}/t_coord_system xy 

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
   set/format I1
   test = m$abs(x_offset({i})) + m$abs(y_offset({i})) 
   if test .ne. 0 then
      $ {tecs_script}/t_offset {x_offset({i})} {y_offset({i})} -
         | awk '{if(NR==1){print $1}}' | write/keyword tel_return	! pipe return
      if tel_return .ne. 0 then
        write/out "ERROR: Telescope return value for t_offset signals an error..."
        write/out "...the program is aborted"
        $play -q /disk-a/staff/GEIRS/SOUNDS/crash.au
        goto exit
      else
         isodate = m$isodate()
         write/file {fctrl(1)} {isodate} {x_offset({i})} {y_offset({i})} {i}
         write/out "         Telescope offsetted by {x_offset({i})} {y_offset({i})} 1/10 arcsec"
      endif
   endif	

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

   if i .ne. n_exp  then
      write/out "         Waiting for {wait_time} seconds"
      wait/sec {wait_time}
   endif
enddo

$cmd_panic_new sync
write/out "         offset series finished ..."
$play -q /disk-a/staff/GEIRS/SOUNDS/gong.au

exit:
close/file {fctrl(1)}
return
