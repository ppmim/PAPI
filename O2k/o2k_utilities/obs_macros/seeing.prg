!  s e e i n g . p r g                                 HJR   9-Sep-03
!
! measure seeing with OMEGA2000
!

define/par P1 zoom ? "Zoom or not : "

P1 = m$lower(P1)
if p1(1:4) .eq. "help" then
   write/out
   write/out "seeing.prg"
   write/out "call: o2k/seeing" 
   write/out 
   write/out "Command line parameter:"
   write/out  "   P1 = zoom or no-zoom"
   write/out
   write/out "If zoom-option is selection, measurement is taken in zoom=4 window"
   write/out "Image to be analysed has to be loaded into MIDAS display before call!"
   write/out "Measure with left, terminate with right mouse-button"
   write/out "Result is stored in MIDAS-table seeing.tbl"
   write/out
   return
endif

if mid_session .ne. 31  then
   write/out "Please use QUICKLOOK (green) MIDAS window to start seeing !"
   $ auplay /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
   goto exit
endif

define/local seeing_x/r/1/2 0,0
define/local seeing_y/r/1/2 0,0
define/local seeing/r/1/2 0,0
define/local help_x/r/1/1 0
define/local help_y/r/1/1 0
define/local fctrl/i/1/2 0,0
define/local ncoadds/i/1/1 0
define/local saturat/r/1/1 {ccd_para(4)}

open/file Seeing.log append fctrl

if P1(1:4) .eq. "zoom"  then
   center/gauss cursor seeing ? ? ? W,4
else
   center/gauss cursor seeing
endif
@@ PM:tvdisply
define/local disp_frame/c/1/24 {tvdisply}
define/local isodate/c/1/17 " "
define/local filter/c/1/32 " "
define/local exp_time/d/1/1 0.
copy/dk {disp_frame} o_time/d/7/1 exp_time/d/1/1
copy/dk {disp_frame} filters/c/1/32 filter/c/1/32
copy/dk {disp_frame} ncoadds/i/1/1 ncoadds/i/1/1
saturat = saturat*ncoadds

compute/tab seeing :xfwhm = :xfwhm*0.44942
compute/tab seeing :yfwhm = :yfwhm*0.44942

name/col seeing :xfwhm F5.2 "arcsec"
name/col seeing :yfwhm F5.2 "arcsec"

set/midas output=logonly
sel/tab seeing :max.le.{saturat}
sel/tab seeing sel.and.:status.eq.0
sel/tab seeing sel.and.:xfwhm.gt.0.5
sel/tab seeing sel.and.:yfwhm.gt.0.5
stat/tab seeing :xfwhm
seeing_x(1) = outputr(3)
seeing_x(2) = outputr(4)
stat/tab seeing :yfwhm
seeing_y(1) = outputr(3)
seeing_y(2) = outputr(4)
seeing(1) = m$sqrt(seeing_x(1)*seeing_y(1))
help_x = seeing_x(1)**2 / (seeing_x(1)**2 + seeing_y(1)**2)*seeing_x(2)**2
help_y = seeing_y(1)**2 / (seeing_x(1)**2 + seeing_y(1)**2)*seeing_y(2)**2
seeing(2) = m$sqrt(help_x + help_y)
set/midas output=yes

set/format F5.2
write/out
write/out "Seeing in X = ({seeing_x(1)} +/- {seeing_x(2)}) arcsec"
write/out "Seeing in Y = ({seeing_y(1)} +/- {seeing_y(2)}) arcsec"
write/out "------------------------------------------------------------------------------"
write/out "Seeing =  ({seeing(1)} +/- {seeing(2)}) arcsec"

! write log
isodate = m$isodate()

set/format F5.1
write/file {fctrl(1)} {isodate} {filter} {exp_time}sec {seeing(1)} +/- {seeing(2)}

exit:
clear/chan overlay
close/file {fctrl(1)}
return
