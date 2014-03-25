!   a c q u i s i t i o n . p r g           HJR    14-Jul-04
!
!   modified by jmiguel@iaa.es  for PANIC   2-March-2011
!
!   acquire a filed for PANIC
!
!
define/par P1 ? ? "Telescope position and equinox : "
define/par P2 =
define/par P3 0,0 N/A "Position of alignment star [Xpix,Ypix] : "
define/par P4 2 N/A "Single frame integration time : "
define/par P5 10 N/A "Number of frames to be added in memory : "
define/par P6 NO_FLAT ? "Flat field frame : "

define/local last_ima/c/1/64 " "
define/local answer/c/1/16 " "                         !
define/local back/c/1/8 " "
define/local i/i/1/1 0

!!!!!!!!!!!!!!!!!!!!!! ENVIRONMENT VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

define/local geirslstabort/c/1/256  " "
geirslstabort = M$SYMBOL("GEIRSLSTABORT")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

$cmd_panic_new tele absolute {P1}         ! move telescope to first field
$cmd_panic_new sync tele

$cmd_panic_new crep  {P5}
$cmd_panic_new itime {P4}
$cmd_panic_new sync
take_AQ:
$cmd_panic_new read
$cmd_panic_new sync
$cmd_panic_new save -i

! load flatfield corrected image into display for pixel-accurate telescope alignment

set/midas output=logonly        ! suppress MIDAS output to screen
$cmd_panic_new last | awk '{print $2}' | write/keyword last_ima
set/midas output=yes            ! re-activate MIDAS output

if P6(1:7) .ne. "NO_FLAT" then
   i = m$exist("{P6}.fits")
   if i .le. 0  then
      write/out "Flatfield not found.  Exiting procedure...!"
      goto exit
   endif
   i = m$len(P6)
endif

if P3(1:3) .eq. "0,0"  then
   write/out
   write/out "Skipping pixel-alignment!"
   write/out "Displaying image for field verification only"
   if P6(1:7) .ne. "NO_FLAT"  then
      compute AQ = {last_ima}.fits / {P6(1:{i})}.fits  ! flatfield acquisition frame
      del/descr AQ lhcuts              ! force new calculation of cut levels
      load AQ ce=c sc=-3,a             ! display image, centered on alignment star
   else
      load {last_ima}.fits sc=-3,a
   endif

   inquire/key answer "Was acquistion successful and shall we proceed to take data ?"
   answer = m$lower(answer)
   if answer(1:1) .eq. "n"  then
      back(1:4) = "stop"
   else
      back(1:7) = "proceed"
   endif
else
   compute AQ = {last_ima}.fits / {P6(1:{i})}.fits  ! flatfield acquisition frame
   del/descr AQ lhcuts              ! force new calculation of cut levels
   load AQ ce={P3} sc=1             ! display image, centered on alignment star
   back/det
   panic/offset {P3}                  ! position alignment star to this position
   inquire/key answer "Do you want to check the alignment [yes/no/stop] : "
   answer = m$lower(answer)
   if answer(1:1) .eq. "y" then
      goto take_AQ
   else
      if answer(1:1) .eq. "s" then
         back(1:4) = "stop"
      else
         back(1:7) = "proceed"
      endif
      goto exit
   endif
endif

exit:
return {back}
