!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! !
! offset.prg
! purpose: computes offset between two positions
! use as @@ offset x,y     or panic/offset  x,y
!       x,y = the position [pixels] you want to go to
! this procedure works with rotated cass flange too!
! original: jwf mar 97
! modified for OMEGA2000:   HJR  15-Sep-03
! modified for PANIC:   jmiguel@iaa.es  01-March-11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

defi/par p1 ?  n "New  position [X_pix,Y_pix] : "

P1 = m$lower(P1)
if p1(1:4) .eq. "help" then
   write/out
   write/out "offset.prg"
   write/out "call: panic/offset [new_X,new_Y]" 
   write/out 
   write/out "The command line parameters are:"
   write/out "  P1=new_X,new_Y: destination position in pixels for object to be marked"
   write/out
   write/out "Telescope will be moved to place object (current position measured"
   write/out "     via CENTER/GAUSS on display) to new position."
   write/out
   goto exit
endif

if mid_session .ne. 31  then
   write/out "Please use QUICKLOOK (green) MIDAS window to start offset !"
   $ play -q $GEIRS_DIR/SOUNDS/sorrydave.au
   goto exit
endif


defi/local owanted/r/1/2 {p1}
defi/local pixsize/r/1/1 18.              ! pixel size in mue
defi/local tel_scal/r/1/1/ 40.0516221     ! telescope scale mue/arcsec
defi/local pix/r/1/1 0.
defi/local step/d/1/2 1,1
defi/local start/d/1/2 1,1


!!!!!!!!!!!!!!!!!!!!!! ENVIRONMENT VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

define/local tecs_script/c/1/256  " "
tecs_script = M$SYMBOL("TECS_SCRIPT")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

comp/key pix = {pixsize}/{tel_scal}       ! arcsec / pixel    PENDIENTE

set/format f6.4

defi/local opos/r/1/2 0.,0.      ! position of object [pixel]
defi/local off/r/1/2 0.,0.

write/out " "

write/out "             mark object with cursor box"
write/out "             modify box size with arrow keys"
write/out "             left mouse button to input position"
write/out "             right mouse button to terminate input"
write/out " "

set/midas output=logonly
cen/gauss ? ? ? 2,1,1
set/midas output=yes
write/out " "

set/format f10.2
copy/kk outputr/r/5/2 opos/r/1/2

write/out " the object is at pixel {opos(1)},{opos(2)}"


! compute offset

comp/key off(1) = opos(1)-owanted(1)
comp/key off(2) = opos(2)-owanted(2)

set/format f10.2
write/out  the offset to pixel {owanted(1)},{owanted(2)}  is :
write/out "      x = {off(1)}  y = {off(2)} pixel"
comp/key off(1) = -off(1)*pix
comp/key off(2) = off(2)*pix
write/out " "
write/out " or  RA = {off(1)}  Del = {off(2)} arcsec"
write/out " "

$ play -q $GEIRS_DIR/SOUNDS/doorbell.au
defi/local move/c/1/1 a
inquire/key move "move by hand or automatically [def]? [h/a]"


if move .eq. "a" then
    write/out selecting xy coord system, offsetting the telescope ...
    $ {tecs_script}/t_coord_system  xy
    $ {tecs_script}/t_offset {off(1)} {off(2)}
    define/local offset/i/1/2 
    owanted(1) = owanted(1) *10
    owanted(2) = owanted(2) *10
    offset(1) = M$NINT(owanted(1))
    offset(2) = M$NINT(owanted(2))
    define/local file_id_1/i/1/2            ! file id for open/file
    open/file /home/panic/tmp/offset_ok.dat write file_id_1  ! PENDIENTE
    set/format I5
    write/file {file_id_1(1)} {owanted(1)},{owanted(2)}
    close/file {file_id_1(1)}
    set/format f10.2
    
    
  
endif

if move .eq. "h" then
    write/out Select the x/y coordinate system , then move the telescope manually !
endif
$ play -q $GEIRS_DIR/SOUNDS/gong.au

exit:
return
