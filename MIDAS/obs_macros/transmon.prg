!  t r a n s m o n . p r g                                 HJR   16-Jul-04
!
! monitoring atmospheric transmission with OMEGA2000
!
! modified for PANIC:   jmiguel@iaa.es  01-March-11

define/par P1 ? N/A "Magnitude of star : "
define/par P2 ? ? "Filter used : "

P1 = m$lower(P1)
if p1(1:4) .eq. "help" then
   write/out
   write/out "transmon.prg"
   write/out "call: panic/transmon" 
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

if mid_session .ne. 31  then
   write/out "Please use QUICKLOOK (green) MIDAS window to start transmon !"
   $ play -q $GEIRS_DIR/SOUNDS/sorrydave.au
   goto exit
endif

define/local fctrl/i/1/2 0,0
define/local sat_value/r/1/1 30000.          !    <--------  SATURATION LEVEL
open/file Transmission.log append fctrl

write/out "Mark object first, then sky background ..."

clear/chan over
get/curs CR_mon N ? 2,2

      set/midas output=logonly
stat ? CR_mon ? ? FNNN CR_mon,A
      set/midas output=yes

define/local xstart/i/1/1 0
define/local ystart/i/1/1 0
define/local xend/i/1/1 0
define/local yend/i/1/1 0
define/local n_star/i/1/1 0
define/local n_back/i/1/1 0

define/local star/r/1/1 0
define/local back/r/1/1 0
define/local signal/r/1/1 0
define/local expected/r/1/1 0
define/local ratio/r/1/1 0

xstart = m$value(CR_mon,:Xstartpix,@1)
ystart = m$value(CR_mon,:Ystartpix,@1)
xend = m$value(CR_mon,:Xendpix,@1)
yend = m$value(CR_mon,:Yendpix,@1)
n_star = (xstart-xend+1)*(ystart-yend+1)

xstart = m$value(CR_mon,:Xstartpix,@2)
ystart = m$value(CR_mon,:Ystartpix,@2)
xend = m$value(CR_mon,:Xendpix,@2)
yend = m$value(CR_mon,:Yendpix,@2)

n_back = (xstart-xend+1)*(ystart-yend+1)

star = m$value(CR_mon,:tot_intens,@1)
back = m$value(CR_mon,:tot_intens,@2)

signal = star - back*n_star/n_back

! write log
@@ PM:tvdisply   ! PENDIENTE 
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
max_signal = m$value(CR_mon,:MAX,@1)
max_signal = max_signal / coadds
if max_signal .gt. {ccd_para(4)} then
   write/out
   write/out "The star seems to be saturated!"
   write/out "No log entry possible!"
   write/out
   $play -q $GEIRS_DIR/SOUNDS/sorrydave.au
   goto exit
endif

isodate = m$isodate()

l = m$len(P2)
if P2(1:{l}) .eq. "J"  then
   define/local count_rate/r/1/1 4.78e9   ! PENDIENTE
endif
if P2(1:{l}) .eq. "H"  then
   define/local count_rate/r/1/1 2.22e9   ! PENDIENTE 
endif
if P2(1:{l}) .eq. "Kprime"  then
   define/local count_rate/r/1/1 1.67e9   ! PENDIENTE
endif
if P2(1:{l}) .eq. "K"  then
   define/local count_rate/r/1/1 1.88e9   ! PENDIENTE
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
