!  t r a n s m o n . p r g                                 HJR   16-Jul-04
!
! monitoring atmospheric transmission with OMEGA2000
!

define/par P1 ? N/A "Magnitude of star : "
define/par P2 ? ? "Filter used : "

P1 = m$lower(P1)
if p1(1:4) .eq. "help" then
   write/out
   write/out "transmon.prg"
   write/out "call: o2k/transmon" 
   write/out 
   write/out "Command line parameter:"
   write/out "   P1 = magnitude of the objekt to be used for monitoring"
   write/out "   P2 = filter used for the exposure"
   write/out "        supported are J, H, Kprime and K"
   write/out
   write/out "Result is appended to log file Transmission.log"
   write/out
   got to exit
endif

define/local fctrl/i/1/2 0,0
define/local sat_value/r/1/1 30000.          !    <--------  SATURATION LEVEL
open/file Transmission.log append fctrl

clear/chan over
cen/gauss ? ? ? 2,1,1

define/local fwhm_x/r/1/1 0
define/local fwhm_y/r/1/1 0
define/local peak/r/1/1 0
define/local signal/r/1/1 0

fwhm_x = outputr(12)
fwhm_y = outputr(13)
peak = outputr(14)

define/local expected/r/1/1 0
define/local ratio/r/1/1 0

signal = 1.13*fwhm_x*fwhm_y*peak

! write log
@@ PM:tvdisply
define/local l/i/1/1 0
define/local disp_frame/c/1/24 {tvdisply}
define/local isodate/c/1/17 " "
define/local filter/c/1/32 " "
define/local exp_time/d/1/1 0.
define/local coadds/i/1/1 0
define/local max_signal/r/1/1 0
copy/dk {disp_frame} o_time/d/7/1 exp_time/d/1/1
copy/dk {disp_frame} filters/c/1/32 filter/c/1/32
copy/dk {disp_frame} ncoadds/i/1/1 coadds/i/1/1
max_signal = peak / coadds
if max_signal .gt. {ccd_para(4)} then
   write/out
   write/out "The star seems to be saturated!"
   write/out "No log entry possible!"
   write/out
   $play -q /disk-a/staff/GEIRS/SOUNDS/crash.au
   goto exit
endif

isodate = m$isodate()

l = m$len(P2)
if P2(1:{l}) .eq. "J"  then
   define/local count_rate/r/1/1 4.78e9
endif
if P2(1:{l}) .eq. "H"  then
   define/local count_rate/r/1/1 2.22e9
endif
if P2(1:{l}) .eq. "Kprime"  then
   define/local count_rate/r/1/1 1.67e9
endif
if P2(1:{l}) .eq. "K"  then
   define/local count_rate/r/1/1 1.88e9
endif

expected = count_rate * 10 ** (-{P1}/2.5)
ratio = signal / exp_time / expected * 100

set/format F5.1
write/file {fctrl(1)} {isodate} {filter} {exp_time}sec {ratio}%
write/out
write/out "Atmospheric transmission in {filter} is {ratio}% of expected value."

exit:
close/file {fctrl(1)}
return
