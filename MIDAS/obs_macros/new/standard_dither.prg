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
abort_check = m$exist("{geirslstabort}")
if abort_check .eq. 1 then
$rm {geirslstabort} 
endif

$cmd_panic_new object "{P3} (AQ)"
$cmd_panic_new crep 10
$cmd_panic_new itime 2

$cmd_panic_new read
$cmd_panic_new sync

    ! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("{geirslstabort}")
if abort_check .eq. 1 then
  write/out "Program is aborted..."
  $rm {geirslstabort}
  $play -q $GEIRS_DIR/SOUNDS/crash.au
  goto exit
endif

$cmd_panic_new last   ! writes last filename in file geirsLstFile
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

$cmd_panic_new last   ! writes last filename in file geirsLstFile
    ! writes last filename in keyword framename
write/keyword framename </disk-a/o2k/tmp/geirsLstFile 


compute/ima test = {framename} / {P5}
$play -q $GEIRS_DIR/SOUNDS/whistle.au

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
