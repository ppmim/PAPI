defi/loc abort_check/i/1/1/ 0
defi/loc answ/c/1/1/ 

write/out " "
write/out "Programme to take a series of dark exposures with O-2000"
write/out "(1.6  2.0  3.0  4.0  5.0  6.0 seconds "
write/out " "

write/out " "
write/out "     Programme can be aborted with the command - touch abortfile -"  
write/out "     To this purpose it is recommended to open an xterminal from" 
write/out "     the same directory where the programme is started"
write/out " "
write/out " "
write/out " "
inq/key answ " - Set correct output filename - done >"
write/out " "

! --- Taking darks ---

$cmd_panic_new filter BLANK
$cmd_panic_new sync

write/out " Dark 1.6 s "
o2k/calser dark1.6 = 1.6,1.6 10 5 0 i 

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

write/out " Dark 2.0 s "
o2k/calser dark2.0 = 2.0,2.0 10 5 0 i
 
abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

write/out " Dark 3.0 s "
o2k/calser dark3.0 = 3.0,3.0 10 5 0 i

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

write/out " Dark 4.0 s "
o2k/calser dark4.0 = 4.0,4.0 10 5 0 i

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

write/out " Dark 5.0 s "
o2k/calser dark5.0 = 5.0,5.0 10 5 0 i

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif

write/out " Dark 6.0 s "
o2k/calser dark6.0 = 6.0,6.0 10 5 0 i

abort_check = m$exist("abortfile")
if abort_check .eq. 1 then
   write/out "         Program is aborted..."
   $rm -f abortfile
goto exit
endif


exit:
return
