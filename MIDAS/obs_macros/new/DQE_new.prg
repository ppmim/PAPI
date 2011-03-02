! Measure relative DQE of OMEGA2000  (September 2003   HJR)
!
!
!     --------------------------------
!     |   <---- -RA              (R) |      star (arrows) and 
!     |             +RA        ^     |      telescope  (RA/DEC signs)
!     |             -DEC     /    |  |      movements 
!     |                    /      |  |
!     |                  /        v  |
!     |                /        +DEC |
!     |             (S)              |
!     |                              |
!     |                              |
!     |                              |
!     |                              |
!     |                              |
!     |                              |
!     --------------------------------
!
!  Measure relative DQE by positioning the telescope from the 
!  start position (S) first into the reference position (R) by 
!  offsetting with ref_off(1,2). Then the object is shifted first 
!  in -X direction (telescope in -RA). 
!  After each row the telescope is moved back to the start position 
!  for flux check with reference to the first image. Then the 
!  next row is started by moving the star in -Y (telescope in +DEC)
!  direction. 
!  In the end a final image is taken again at the start position (S).
!
! ===================================================================!
define/par P1 " " C/A "RA (J2000) of reference position [Format hh:mm:ss.s]: " 
define/par P2 " " C/A "DEC (J2000) of reference position [Format +dd:mm:ss]: " 
define/par P3 3,2 N/A "Number of pointings in X,Y : "
define/par P4 " " ? "Filter used : "
define/par P5 1,2.0 N/A "Integration parameters crep,itime : "
define/par P6 SAVE C/A "Save images [if not give NO_SAVE] : "
define/par P7 1,1 N/A "Start position in X,Y [1,1] : "

! define local variables
define/local ans/c/1/32 " "
define/local save/c/1/32 "{P6}"
save = m$upper(save)
define/local X_step/i/1/1 0
define/local Y_step/i/1/1 0
define/local x_offset/r/1/1 0
define/local y_offset/i/1/1 0
define/local ix/i/1/1 0
define/local iy/i/1/1 0
define/local crep/i/1/1 0
define/local itime/r/1/1 0.
define/local first_x/i/1/1 1
define/local first_y/i/1/1 1
define/local jump_x/i/1/1 0
define/local jump_y/i/1/1 0
define/local test_jump/i/1/1 0
DEFINE/LOCAL imaname/c/1/80 " "
DEFINE/LOCAL starposx/r/1/1 0
DEFINE/LOCAL starposy/r/1/1 0
DEFINE/LOCAL corr_off/i/1/2 0,0
DEFINE/LOCAL cg_flag/i/1/1 0

define/local abort_check/i/1/1 0    ! for return value of abort-file-check
  ! abort check: 0=does not exist ; 1=abort file exists
  abort_check = m$exist("{geirslstabort}")
  if abort_check .eq. 1 then
    $rm {geirslstabort}  ! remove file to reset
  endif

write/key outputr/r/1/2 {P5}
crep = m$nint(outputr(1))
itime = outputr(2)

write/key outputi/i/1/2 {P7}
first_x = outputi(1)
first_y = outputi(2)

! Raster festlegen
define/local abs_RA/C/1/16 "{P1}"
define/local abs_DEC/C/1/16 "{P2}"
write/key outputi/i/1/2 {P3}
define/local n_X/I/1/1 {outputi(1)}
define/local n_Y/I/1/1 {outputi(2)}
define/local ref_off/i/1/2 4400,4400     ! in 1/10arcsec = 7.5arcmin

define/local tel_wait/i/1/1 0.  ! wait [sec] after telescope movement

X_step = 14.9*600/(n_X-1) ! telescope offsets [in arsec if given with decimal point!]
Y_step = 14.9*600/(n_Y-1)

Set/format I1 F6.1
! set image parameters
$cmd_panic_new crep {crep}
$cmd_panic_new itime {itime}
$cmd_panic_new sync

$cmd_panic_new filter {P4}
$cmd_panic_new sync

!-------------------------------------------------------------

inquire/key ans "Are you ready to preset telescope to -
RA(J2000) {abs_RA}, DEC(J2000) {abs_DEC} [no] : "

ans = m$upper(ans)

if ans(1:1) .eq. "Y"  then
   ! telescope to reference position (R) via start (S)
   $ {tecs_script}/t_posit {abs_RA} {abs_DEC} 2000.0
      $play -q $GEIRS_DIR/SOUNDS/doorbell.au
   inquire/key ans "All set with telescope centered on start position [no] ?"
   ans = m$upper(ans)
   if ans(1:1) .ne. "Y"  then
      write/out
      write/out "Program stopped ..."
      Write/out
      goto ende
   endif
   $ {tecs_script}/t_offset +{ref_off(1)} -{ref_off(2)}     ! position of first image
   wait/secs {tel_wait}
else
   write/out
   write/out "Program stopped ..."
   write/out
   goto ende
endif

write/out "Taking DQE image {P4} [{ix},{iy}]"

  ! write object
$cmd_panic_new object DQE reference {P4} start

$cmd_panic_new read
$cmd_panic_new sync
if save(1:7) .ne. "NO_SAVE" then
   $cmd_panic_new save -i
endif

$cmd_panic_new sync
! get position of star
$cmd_panic_new last 
write/keyword imaname </disk-a/o2k/tmp/geirsLstFile

@@ O2K_UTIL:/obs_macros/new/cuts {imaname}
starposx = {Q1}
starposy = {Q2} 
WRITE/OUT {starposx} {starposy}

! Set telescope to start position if not [1,1]

test_jump = first_x + first_y
if test_jump .gt. 2 then
   jump_x = (first_x - 1) * x_step
   jump_y = (first_y - 1) * y_step
   $ {tecs_script}/t_offset -{jump_x} {jump_y}
endif

do iy = {first_y} {n_Y}
 do ix = {first_x} {n_X}
  write/out "Taking DQE image {P4} [{ix},{iy}]"

  ! write object
  $cmd_panic_new object DQE measurement {P4} {ix},{iy}

  $cmd_panic_new read
  $cmd_panic_new sync
  if save(1:7) .ne. "NO_SAVE" then
     $cmd_panic_new save -i
  endif

  if ix .lt. n_X  then
    ! offset telescope in RA  (move star in +alpha direction)
    $ {tecs_script}/t_offset -{X_step} 0
    wait/secs {tel_wait} ! wait for telescope  before read
  endif

  ! abort check: 0=does not exist ; 1=abort file exists
  abort_check = m$exist("{geirslstabort}")
  if abort_check .eq. 1 then
      write/out "Program was aborted from GUI ..." 
    $rm {geirslstabort}  ! remove file again
              $play -q $GEIRS_DIR/SOUNDS/crash.au
      goto ende
  endif
 enddo
 first_x = 1

 if iy .lt. n_Y  then
     ! telescope to start, again via reference position (R)
  $ {tecs_script}/t_posit {abs_RA} {abs_DEC} 2000.0
            $play -q $GEIRS_DIR/SOUNDS/doorbell.au
  inquire/key ans "Telescope at start position [no] ? "
  ans = m$upper(ans)
  if ans(1:1) .ne. "Y" then
     write/out
     write/out " Program stopped"
     write/out
     goto ende
  endif
   
  $ {tecs_script}/t_offset +{ref_off(1)} -{ref_off(2)}
  wait/secs {tel_wait}

inquire/key ans "Move telescope by 10arcsec in RA and 20arcsec in DEC "

  write/out "Taking image at reference position {iy} ..."
  $cmd_panic_new object DQE reference {P4} {iy}
  $cmd_panic_new read
  $cmd_panic_new sync
  if save(1:7) .ne. "NO_SAVE" then 
    $cmd_panic_new save -i
  endif
  
  $cmd_panic_new sync      !get filename of last saved image
  $cmd_panic_new last
  write/keyword imaname </disk-a/o2k/tmp/geirsLstFile
  
  @@ O2K_UTIL:/obs_macros/new/autoguide {imaname} {starposx} {starposy}
  
  corr_off(1) = {Q1}
  corr_off(2) = {Q2}
  write/key cg_flag/i/1/1 {Q3}
  if cg_flag .eq. 0 then
     WRITE/OUT Corrections from autoguide: {corr_off(1)} {corr_off(2)}
  endif

  !starpos not changed, should be in correct position if everything ok

  ! offset telescope in Y (move star in -delta direction)
  y_offset = iy*Y_step + corr_off(2)
  $ {tecs_script}/t_offset {corr_off(1)} {y_offset}
  wait/secs {tel_wait}
 endif

 ! abort check: 0=does not exist ; 1=abort file exists
 abort_check = m$exist("{geirslstabort}")
 if abort_check .eq. 1 then
     write/out "Program was aborted from GUI ..." 
     $rm {geirslstabort}  ! remove file again
        $play -q $GEIRS_DIR/SOUNDS/crash.au
     goto ende
 endif
enddo

write/out "Taking last image at starting position..."

! telescope to start for final reference
$ {tecs_script}/t_posit {abs_RA} {abs_DEC} 2000.0
 ans = m$upper(ans)
 if ans(1:1) .ne. "Y" then
    write/out
    write/out " Program stopped"
    write/out
    goto ende
 endif
$ {tecs_script}/t_offset +{ref_off(1)} -{ref_off(2)}
wait/secs {tel_wait}
$cmd_panic_new object DQE reference {P4} final
$cmd_panic_new read
$cmd_panic_new sync
if save(1:7) .ne. "NO_SAVE" then 
 $cmd_panic_new save -i
endif

write/out
write/out "         all done ..."
write/out
$play -q $GEIRS_DIR/SOUNDS/gong.au

ende:
return
