!		LOGIN_quick.PRG        for OMEGA2000 quicklook

define/local ans/c/1/80
inquire/key ans "Specify data directory : "
ch/dir {ans}

monitpar(20) = 2000

WRITE/KEY ERR_CTRL/I/1/3 0,2,1
write/key mid$sys/c/21/10 "$debu     "
write/key user/c/1/20 "panic"
write/key INSTR_ID/C/1/8 PANIC
write/key DETECTOR/C/1/8 PANIC-HAWAI

!@@ PM:init OMEGA2k OMEGA2k 1,1 20 ? ? ? silent
!@@ PM:maincmd
set/cont panic

!create/gui help

!create/gra 0 687,476,45,544

!init/display p5=RGBQ
!create/disp 0 1024,1024,750,468 4
!load/lut rainbow3
!display/lut

!load/itt neg 0
!load/itt neg 1
!load/itt neg 2
!load/itt neg 3

set/midas prompt=QUICK
write/key mid_session/i/1/1 31

return
