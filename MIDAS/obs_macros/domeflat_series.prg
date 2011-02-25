defi/loc z/i/1/1/ 0
defi/loc Y/i/1/1/ 0
defi/loc J/i/1/1/ 0
defi/loc H/i/1/1/ 0
defi/loc K/i/1/1/ 0
defi/loc Ks/i/1/1/ 0
defi/loc Kp/i/1/1/ 0
defi/loc nb2122/i/1/1/ 0
defi/loc nb2144/i/1/1/ 0
defi/loc mon/i/1/1/ 0
defi/loc moff/i/1/1/ 0
defi/loc abort_check/i/1/1/ 0
defi/loc answ/c/1/1/ 

write/out
write/out "Programme to calibrate a series of O-2000 filters"
write/out

write/out " "
write/out "     Programme can be aborted with the command - touch abortfile -"  
write/out "     To this purpose it is recommended to open an xterminal from" 
write/out "     the same directory where the programme is started"
write/out " "
write/out " "
write/out " "
inq/key answ " - Set correct output filename - done >"
write/out " "
write/out " - Set the telescope to the domewall flatfield position and"
inq/key answ " make sure that no FFL is running - OK >"
write/out " "

$flats ALLOFF

write/out
write/out "Select filters to calibrate [y/n(def)]"
write/out

answ = "n"
inq/key answ " z: "
if answ .eq. "Y" .or. answ .eq. "y" then
  z = 1
endif

answ = "n"
inq/key answ " Y: "
if answ .eq. "Y" .or. answ .eq. "y" then
  Y = 1
endif

answ = "n"
inq/key answ " J: "
if answ .eq. "Y" .or. answ .eq. "y" then
  J = 1
endif

answ = "n"
inq/key answ " H: "
if answ .eq. "Y" .or. answ .eq. "y" then
  H = 1
endif

answ = "n"
inq/key answ " K: "
if answ .eq. "Y" .or. answ .eq. "y" then
  K = 1
endif

answ = "n"
inq/key answ " Ks: "
if answ .eq. "Y" .or. answ .eq. "y" then
  Ks = 1
endif

answ = "n"
inq/key answ " Kp: "
if answ .eq. "Y" .or. answ .eq. "y" then
  Kp = 1
endif

answ = "n"
inq/key answ " NB2122: "
if answ .eq. "Y" .or. answ .eq. "y" then
  nb2122 = 1
endif

answ = "n"
inq/key answ " NB2144: "
if answ .eq. "Y" .or. answ .eq. "y" then
   nb2144 = 1
endif

answ = "n"
inq/key answ " meth_on: "
if answ .eq. "Y" .or. answ .eq. "y" then
  mon = 1
endif

answ = "n"
inq/key answ " meth_off: "
if answ .eq. "Y" .or. answ .eq. "y" then
  moff = 1
endif

! --- Taking flatfields ---

if z .eq. 1 then
   write/out "calibrating z ..."
   $cmd_panic_new filter Z
   $cmd_panic_new sync
   o2k/domeflats "dome_z" = 5,10 6 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if Y .eq. 1 then
   write/out "calibrating Y ..."
   $cmd_panic_new filter Y
   $cmd_panic_new sync
   o2k/domeflats "dome_Y" = 5,9 10 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if J .eq. 1 then
   write/out "calibrating J ..."
   $cmd_panic_new filter J
   $cmd_panic_new sync
   o2k/domeflats "dome_J" = 5,9 2.0 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if H .eq. 1 then
   write/out "calibrating H ..."
   $cmd_panic_new filter H
   $cmd_panic_new sync
   o2k/domeflats "dome_H" = 5,1 2.0 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if K .eq. 1 then
   write/out "calibrating K ..."
   $cmd_panic_new filter K
   $cmd_panic_new sync
   o2k/domeflats "dome_K" = 5,1 1.6 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if Ks .eq. 1 then
   write/out "calibrating Ks ..."
   $cmd_panic_new filter KS
   $cmd_panic_new sync
   o2k/domeflats "dome_Ks" = 5,1 1.6 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if Kp .eq. 1 then
   write/out "calibrating K-prime ..."
   $cmd_panic_new filter K-PRIME
   $cmd_panic_new sync
   o2k/domeflats "dome_Kp" = 5,1 1.6 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if nb2122 .eq. 1 then
   write/out "calibrating NB2122 ..."
   $cmd_panic_new filter NB2122 
   $cmd_panic_new sync
   o2k/domeflats "dome_NB2122" = 3 1.6 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if nb2144 .eq. 1 then
   write/out "calibrating NB2144 ..."
   $cmd_panic_new filter NB2144 
   $cmd_panic_new sync
   o2k/domeflats "dome_NB2144" = 3 1.6 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if mon .eq. 1 then
   write/out "calibrating methane_on ..."
   $cmd_panic_new filter METHANE_ON
   $cmd_panic_new sync
   o2k/domeflats "dome_meth_on" = 5,9 4.0 10 10 i
endif   

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

if moff .eq. 1 then
   write/out "calibrating methane_off ..."
   $cmd_panic_new filter METHANE_OFF
   $cmd_panic_new sync
   o2k/domeflats "dome_meth_off" = 5,10 4.0 10 10 i
endif   


exit:
return
