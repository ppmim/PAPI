! This prg should speed up observation !

define/parameter P1 ? C "Enter a name for the image catalog :"
define/parameter P2 ? N "Enter the integration times in secs[tot,integ,single] :"
define/parameter P3 ? C "Enter an object name :"
define/parameter P4 ? N "Enter position of offset star :"
define/parameter P5 ? ? "Enter name of flatfield :"

define/local int_time/r/1/3 {P2}
define/local pointings/r/1/1 {int_time(1)}
define/local integrated/r/1/1 {int_time(2)}
define/local single/r/1/1 {int_time(3)}
define/local framename/c/1/200
define/local abort_check/i/1/1  
define/local offset/r/1/2 {P4} 
define/local abort/c/1/80  
    
    ! remove existing abort file
    ! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
if abort_check .eq. 1 then
$rm /disk-a/o2k/tmp/geirsLstAbort 
endif

$cmd_o2000 object "{P3} (AQ)"
$cmd_o2000 crep 10
$cmd_o2000 itime 2

$cmd_o2000 read
$cmd_o2000 sync

    ! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
if abort_check .eq. 1 then
  write/out "Program is aborted..."
  $rm /disk-a/o2k/tmp/geirsLstAbort
  $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
  goto exit
endif

$cmd_o2000 last   ! writes last filename in file geirsLstFile
    ! writes last filename in keyword framename
write/keyword framename </disk-a/o2k/tmp/geirsLstFile 


compute/ima test = {framename} / {P5}

load/ima test sc=2 ce={offset(1)},(offset(2)}
back/det ? 2,5
o2k/offset {offset(1)},(offset(2)}

inquire/key abort "Everything allright? y/n"

if abort .eq. n then
 goto exit
endif 

o2k/dither {P1}_I {pointings},{integrated},{single} {P3}

$cmd_o2000 last   ! writes last filename in file geirsLstFile
    ! writes last filename in keyword framename
write/keyword framename </disk-a/o2k/tmp/geirsLstFile 


compute/ima test = {framename} / {P5}
$auplay /disk-a/staff/GEIRS/SOUNDS/whistle.au

load/ima test sc=2 ce={offset(1)},(offset(2)}
back/det ? 2,5
o2k/offset {offset(1)},(offset(2)}

inquire/key abort "Everything allright? y/n"

if abort .eq. n then
 goto exit
endif 

o2k/dither {P1}_II {pointings},{integrated},{single} {P3} 26,AQ



exit:
return
