!********************************************************
!
!     maskalign.prg   03-July-97   E.T.
!
!     new version 27-Sept-97   HJR+JWF
!
!*******************************************************

define/par p1 ? ? "Name of field : "
define/par p2 ? ? "Aperture position for mask : "
define/par p3 ? ? "Exposure time for final full field verification of alignment : "
define/par p4 1 ? "Exposure time for holes with internal lamp : "
define/par p5 10 ? "Exposure time for stars in holes : "

if p1(1:4) .eq. "help" then
	write/out maskalign takes 
	write/out 1.) an image of the mask with the build in lamp
	write/out 2.) a short exposure of the sky to identify the 2 bright stars
	write/out "   calculate tel.offset, change in PA of Cassflange"
	write/out 3.) take 2 short exposures around the stars to improve upon
	write/out "    tel.offset and change in PA of Cassflange"
	write/out step 3 may be iterated until optimum is reached
	write/out
endif


defi/loc pe/r/1/2 1247.5,1029.5

!--------------------------------------
defi/loc scale/r/1/2 0.32358,0.32358
!-------------------------------------

defi/loc ans/c/1/2 xx

defi/loc dec/r/1/1 0
defi/loc PA_C_old/r/1/1 0
defi/loc PA_C_new/r/1/1 0

defi/loc cut/r/1/2 500,1000
defi/loc winx/r/1/2 800,1200
defi/loc winy/r/1/2 1400,1800

defi/loc delta_a/r/1/1 0
defi/loc delta_d/r/1/1 0

defi/loc x/r/1/2 0 all
defi/loc y/r/1/2 0 all
defi/loc loch1/i/1/4 0 all
defi/loc loch2/i/1/4 0 all

defi/loc xs/r/1/2 0 all
defi/loc ys/r/1/2 0 all

defi/loc xss/r/1/2 0 all
defi/loc yss/r/1/2 0 all

defi/loc b/r/1/2 0 all
defi/loc bb/r/1/1 0
defi/loc bs/r/1/2 0 all
defi/loc bsb/r/1/1 0

defi/loc a1/d/1/1 0
defi/loc a2/d/1/1 0
defi/loc alpha/d/1/1 0
defi/loc k/r/1/1 0

defi/loc c/r/1/2 0 all
defi/loc cb/r/1/1 0
defi/loc cs/r/1/2 0 all
defi/loc csb/r/1/1 0

defi/loc as1/d/1/1 0
defi/loc as2/d/1/1 0
defi/loc alphas/d/1/1 0
defi/loc absalpha/d/1/1 0


defi/loc hx/d/1/2 0 all
defi/loc hy/d/1/2 0 all
defi/loc px/d/1/2 0 all
defi/loc py/d/1/2 0 all
defi/loc cofr/d/1/2 0 all
defi/loc dcofr/d/1/2 0 all


defi/loc dx/r/1/2 0 all
defi/loc dy/r/1/2 0 all
defi/loc d/r/1/2 0 all
defi/loc dd/r/1/2 0 all

defi/loc n/i/1/1 0
defi/loc text/c/1/12

defi/loc loop/i/1/1 0
defi/loc newfil/i/1/1 0

defi/loc freepos/i/1/1 5
!========================================================================

! counter to keep aquisition frames despite test mode
n = m$existk("NR_AQ") 
if n .eq. 0 then
	write/key NR_AQ/i/1/1 0
endif 

NR_AQ = NR_AQ + 1

! setup acquisition frame to locate alignment holes

write/out "Configuring MOSCA and CCD for mask exposure with internal lamp ..."
write/out
$ camera_chip1_x1 2
$ camera_chip1_y1 1
$ camera_chip1_x2 1149
$ camera_chip1_y2 1024
$ camera_chip1_hbin 2
$ camera_chip1_vbin 2
$ camera_dewar1_exptime {p4}
$ camera_dewar1_exptype science
$ camera_objectname "Mask {P1}"
$ camera_expmode test 

! set filter to free
$ /disk-a/staff/MOSCA/scripts/mosca_filter_1 0 wait
! set desired mask
newfil = {p2}
$ /disk-a/staff/MOSCA/scripts/mosca_mask {newfil} wait
! mirror in, continuum lamp on
$ /disk-a/staff/MOSCA/scripts/mosca_calib_mirror in wait
$ /disk-a/staff/MOSCA/scripts/mosca_lamp_cont on 40% wait


write/out
write/out ****************************************
write/out Exposing mask with internal lamp to
write/out measure positions of alignment holes ...
write/out ****************************************
write/out

$ camera_dewar1_start wait
set/format I4
$ cp test0001.fits AQ{NR_AQ}a.fits

write/out "Moving MOSCA mirror out of beam ..."
write/out
$ /disk-a/staff/MOSCA/scripts/mosca_calib_mirror out wait

do loop = 1  100
	newfil = m$exist("newdata")
	if {newfil} .eq. 0 then
		goto endloop1
	endif
	wait/secs 1.
enddo
endloop1:

set/format F10.5
displ/chan 0
load test0001 chan=0 scale=-2 center=c,c cuts=0,1000

step1:
clear/chan over
set/curs 2 rectangle 250,250,260,260

writ/out " Mark small hole 1 ..."

center/gauss cursor
x(1) = outputr(5)+15.5
x(2) = outputr(6)+61.8

writ/out
writ/out " Mark small hole 2 ..."

center/gauss cursor
y(1) = outputr(5)+15.5
y(2) = outputr(6)+61.8

write/key text/c/1/1 y
inquire/key text "Are you happy with this [y[def]/n] ?"
text = m$lower(text)
if text(1:1) .eq. "n" goto step1



write/out
write/out "Configuring MOSCA and CCD for exposure of field ..."
write/out
$ camera_dewar1_exptime {p5}
$ camera_dewar1_exptype science
$ camera_objectname "1. aquisition for {P1}"

! set mask to free
newfil = freepos
$ /disk-a/staff/MOSCA/scripts/mosca_mask {newfil} wait


write/out
write/out ******************************************
write/out Taking aquisition exposure without mask to
write/out measure positions of alignment stars ....
write/out ******************************************
write/out

$ camera_dewar1_start wait

set/format I4
$ cp test0001.fits AQ{NR_AQ}b.fits
do loop = 1  100
	newfil = m$exist("newdata")
	if {newfil} .eq. 0 then
		goto endloop2
	endif
	wait/secs 1.
enddo
endloop2:
 

stat/ima test0001 [400,400:600,600] ? ? FN
cut(1) = outputr(8) - 2.*outputr(4)
cut(2) = outputr(8) + 5.*outputr(4)
load test0001 chan=0 scale=-2 cuts={cut(1)},{cut(2)} center=c,c

write/key text/c/1/1 n
inquire/key text "Do you want to change the low/high cuts [y/n[def]] ?"
text = m$lower(text)
if text(1:1) .eq. "y" then 
	inquire/key text "New low,high cut ?"
	cuts test0001 {text}
	load test0001 chan=0 scale=-2 center=c,c
endif


copy/dk test0001 CAHA.ADA.POSANG PA_C_old


step2:
clear/chan over
draw/circle {x(1)},{x(2)},10 F intens=red
draw/circle {y(1)},{y(2)},10 F intens=red
writ/out " Mark position of star 1 ... "

center/gauss cursor
xs(1) = outputr(5)
xs(2) = outputr(6)

writ/out " Mark position of star 2 ..."

center/gauss cursor
ys(1) = outputr(5)
ys(2) = outputr(6)

write/key text/c/1/1 y
write/out  "Are you happy with this [y[def]/n/exit] ?"
text = m$lower(text)
if text(1:1) .eq. "n" then
	goto step2
endif

if text(1:4) .eq. "exit" then 
	goto exit
endif




b(1) = x(1) - y(1)
b(2) = x(2) - y(2)

bs(1) = xs(1) - ys(1)
bs(2) = xs(2) - ys(2)

bb = m$sqrt( b(1)*b(1) + b(2)*b(2) )
bsb = m$sqrt( bs(1)*bs(1) + bs(2)*bs(2) )

a1 = ( b(1)*bs(1) + b(2)*bs(2) ) / ( bb*bsb )

a2 = ( b(1)*bs(2) - b(2)*bs(1) ) / ( bb*bsb )



alpha = m$atan( a2/a1 )
k = a1 / m$cos(alpha)


if alpha .lt. 0.01 then
          dx(1) = x(1) - xs(1)
          dx(2) = x(2) - xs(2)
          dy(1) = y(1) - ys(1)
          dy(2) = y(2) - ys(2)
else


dx(1) = x(1) - a1*( xs(1) - pe(1) ) - a2*( xs(2) - pe(2) ) - pe(1)
dx(2) = x(2) + a2*( xs(1) - pe(1) ) - a1*( xs(2) - pe(2) ) - pe(2)

dy(1) = y(1) - a1*( ys(1) - pe(1) ) - a2*( ys(2) - pe(2) ) - pe(1)
dy(2) = y(2) + a2*( ys(1) - pe(1) ) - a1*( ys(2) - pe(2) ) - pe(2)

end if

d(1) = 0.5 * ( dx(1) + dy(1) )
d(2) = 0.5 * ( dx(2) + dy(2) )

dd(1) = dx(1) - dy(1)
dd(2) = dx(2) - dy(2)

dd(1) = dd(1) * scale(1)
dd(2) = dd(2) * scale(2)

d(1) = d(1) * scale(1)
d(2) = -d(2) * scale(2)

set/format F10.2

PA_C_new = PA_C_old + alpha

write/out
write/out "********************************************************"
write/out "* First rotate instrument by {alpha} degrees ..."
write/out *
write/out *   PA_CASS = {PA_C_old}  ----> {PA_C_new}
write/out *
write/out *   k     = {k}
write/out "*-------------------------------------------------------"
write/out "
write/out *Then apply telescope offsets
write/out *
write/out *  du = {d(1)} " arcsec      " (@{dd(1)})
write/out *  dv = {d(2)} " arcsec      " (@{dd(2)})
write/out *
write/out *                (based on estimated center of rotation!)
write/out "********************************************************"
write/out

inquire/key text "Hit <CR> when rotation has been applied ..."
write/out 
write/key text/c/1/1 y
inquire/key text "Move telescope du = {d(1)} dv = {d(2)} y[def]/n :"
text = m$lower(text)
if text(1:1) .eq. "y" then
	$/disk-a/staff/TECS35/scripts/t3_offset {d(1)} {d(2)}
else 
	write/out Move telescope by hand!
endif

write/out
inquire/key text "Turn on TV-Autoguider and hit <CR> when done ..."

clear/chan over
write/out

loch1(1) = 2300-{x(1)}-40.
loch1(2) = {x(2)}-100.
loch1(3) = 2300-{x(1)}+40.
loch1(4) = {x(2)}+40.
set/format I4

write/out
write/out "Configuring MOSCA and CCD for aquisition exposures of alignment holes ...."
write/out

$ camera_chip1_x1 {loch1(1)}
$ camera_chip1_y1 {loch1(2)}
$ camera_chip1_x2 {loch1(3)}
$ camera_chip1_y2 {loch1(4)}
$ camera_chip1_hbin 1
$ camera_chip1_vbin 1
$ camera_dewar1_exptime {p5}
$ camera_dewar1_exptype science
$ camera_objectname "Hole 1 for {P1}"

! set desired mask
newfil = {p2}
$ /disk-a/staff/MOSCA/scripts/mosca_mask {newfil} wait


write/out "------------------------------------------"
write/out " Taking aquisition exposure for hole 1 ..."
write/out "------------------------------------------"

$ camera_dewar1_start wait

set/format I4
$ cp test0001.fits AQ{NR_AQ}c.fits
do loop = 1  100
	newfil = m$exist("newdata")
	if {newfil} .eq. 0 then
		goto endloop3
	endif
	wait/secs 1.
enddo
endloop3:
!stat/ima test0001 [<,<:>,>] ? 125,10000 FN
!cut(1) = outputr(8)-2*outputr(4)
!cut(2) = outputr(8)+5*outputr(4)
!set/format f10.2
!load test0001 chan=0 scale=4 cuts={cut(1)},{cut(2)} center=c,c
load test0001 chan=0 scale=4 cuts=0,1000 center=c,c

copy/dk test0001 CAHA.ADA.POSANG PA_C_old
set/curs 2 rectangle 250,250,280,280

writ/out " Mark position of small hole 1 ... "
center/gauss cursor
x(1) = outputr(5)+15.5
x(2) = outputr(6)+61.8
write/out
write/out "Set cursor on sky background for background determination ..."
set/curs 2 rectangle 250,250,270,270
@@ backdet

writ/out " Mark position of star 1 ..."
set/curs 2 rectangle 250,250,300,300
center/gauss cursor
xss(1) = outputr(5)
xss(2) = outputr(6)


write/out "------------------------------------------"
write/out " Taking aquisition exposure for hole 2 ..."
write/out "------------------------------------------"

loch2(1) = 2300-{y(1)}-40.
loch2(2) = {y(2)}-100.
loch2(3) = 2300-{y(1)}+40.
loch2(4) = {y(2)}+40.
set/format I4
$ camera_chip1_x1 {loch2(1)}
$ camera_chip1_y1 {loch2(2)}
$ camera_chip1_x2 {loch2(3)}
$ camera_chip1_y2 {loch2(4)}
$ camera_chip1_hbin 1
$ camera_chip1_vbin 1
$ camera_dewar1_exptime {p5}
$ camera_dewar1_exptype science
$ camera_objectname "Hole 2 for {P1}"

$ camera_dewar1_start wait


set/format I4
$ cp test0001.fits AQ{NR_AQ}d.fits
do loop = 1  100
	newfil = m$exist("newdata")
	if {newfil} .eq. 0 then
		goto endloop4
	endif
	wait/secs 1.
enddo
endloop4:
!stat/ima test0001 [<,<:>,>] ? 125,10000 FN
!cut(1) = outputr(8)-2*outputr(4)
!cut(2) = outputr(8)+5*outputr(4)
!set/format f10.2
!load test0001 chan=0 scale=4 cuts={cut(1)},{cut(2)} center=c,c
load test0001 chan=0 scale=4 cuts=0,1000 center=c,c

writ/out " Mark position of small hole 2 ..."
set/curs 2 rectangle 250,250,280,280
center/gauss cursor
y(1) = outputr(5)+15.5
y(2) = outputr(6)+61.8
write/out
write/out "Set cursor on sky background for background determination ..."
set/curs 2 rectangle 250,250,270,270
@@ backdet

writ/out " Mark position of star 2 ..."
set/curs 2 rectangle 250,250,300,300
center/gauss cursor
yss(1) = outputr(5)
yss(2) = outputr(6)

!=========================================================================

! try to improve the center of rotation

c(1) = xs(1) - ys(1)
c(2) = xs(2) - ys(2)

cs(1) = xss(1) - yss(1)
cs(2) = xss(2) - yss(2)

cb = m$sqrt( c(1)*c(1) + c(2)*c(2) )
csb = m$sqrt( cs(1)*cs(1) + cs(2) * cs(2) )

as1 = ( c(1)*cs(1) + c(2)*cs(2) ) / ( cb * csb )
as2 = ( c(1)*cs(2) - c(2)*cs(1) ) / ( cb * csb )

alphas = m$atan( as2/as1)

absalpha = m$abs(alphas)

if  {absalpha} .ge. 0.1 then

	hx(1) = xss(1) - ( as1*xs(1) - as2*xs(2) )
	hx(2) = xss(2) - ( as2*xs(1) + as1*xs(2) )

	px(1) = 0.5 * ( hx(1) - hx(2)*( as1 + 1.0 )/as2 )
	px(2) = 0.5 * ( hx(1)*( as1 + 1.0 )/as2 + hx(2) )


	hy(1) = yss(1) - ( as1*ys(1) - as2*ys(2) )
	hy(2) = yss(2) - ( as2*ys(1) + as1*ys(2) )

	py(1) = 0.5 * ( hy(1) - hy(2)*( as1 + 1.0 )/as2 )
	py(2) = 0.5 * ( hy(1)*( as1 + 1.0 )/as2 + hy(2) )


	cofr(1) = 0.5 * ( px(1) + py(1) )
	cofr(2) = 0.5 * ( px(2) + py(2) )

	dcofr(1) = px(1) - py(1)
	dcofr(2) = px(2) - py(2)

	write/out ****************************************************
	write/out *  Total rotation angle  {alphas} degrees
	write/out *
	write/out *       New center of rotation
	write/out *
	write/out *   x-coordinate: {cofr(1)} "(Pixel) " ({dcofr(1)})
	write/out *   y-coordinate: {cofr(2)} "(Pixel) " ({dcofr(2)})
	write/out *
	write/out * Attention: This center of rotation is
	write/out *            reliable only, if no significant
	write/out "*            linear offset had been applied !"
	write/key ans/c/1/1 y
	set/format f10.2
	write/out *
	write/out * Do you want to use these values [y=def]"
	inq/key ans "*   or keep old center [{pe(1)},{pe(2)}] [n] ?"
	write/out ****************************************************
   	if ans(1:1) .eq. "n" then
 	       cofr(1) = pe(1)
	       cofr(2) = pe(2)
     	else
   	end if

else

	write/out ****************************************************
	write/out *  Rotation angle was not large enough to improve
	write/out *  rotation center
	Write/out *  ----> estimated center of rotation will be used !
	write/out ****************************************************

	cofr(1) = pe(1)
	cofr(2) = pe(2)

end if


!======================================================================
loop1:
!======================================================================

bs(1) = xss(1) - yss(1)
bs(2) = xss(2) - yss(2)

bsb = m$sqrt( bs(1)*bs(1) + bs(2)*bs(2) )

a1 = ( b(1)*bs(1) + b(2)*bs(2) ) / ( bb*bsb )

a2 = ( b(1)*bs(2) - b(2)*bs(1) ) / ( bb*bsb )

alpha = m$atan( a2/a1 )

if alpha .lt. 0.01 then
	dx(1) = x(1) - xss(1)
	dx(2) = x(2) - xss(2)
	dy(1) = y(1) - yss(1)
	dy(2) = y(2) - yss(2)
else

	dx(1) = x(1) - a1*( xss(1) - cofr(1) ) - a2*( xss(2) - cofr(2) ) -cofr(1)
	dx(2) = x(2) + a2*( xss(1) - cofr(1) ) - a1*( xss(2) - cofr(2) ) -cofr(2)

	dy(1) = y(1) - a1*( yss(1) - cofr(1) ) - a2*( yss(2) - cofr(2) ) -cofr(1)
	dy(2) = y(2) + a2*( yss(1) - cofr(1) ) - a1*( yss(2) - cofr(2) ) -cofr(2)

end if

d(1) = 0.5 * ( dx(1) + dy(1) )
d(2) = 0.5 * ( dx(2) + dy(2) )

dd(1) = dx(1) - dy(1)
dd(2) = dx(2) - dy(2)


dd(1) = dd(1) * scale(1)
dd(2) = dd(2) * scale(2)

d(1) = d(1) * scale(1)
d(2) = -d(2) * scale(2)

PA_C_new = PA_C_old + alpha

set/format f10.3
write/out *************************************************************
write/out * First rotate by {alpha} degrees
write/out *
write/out *   PA_CASS = {PA_C_old}  ----> {PA_C_new}
write/out *
write/out * then apply offsets
write/out *
write/out *  du = {d(1)} " arcsec    " (ddu = @{dd(1)})
write/out *  dv = {d(2)} " arcsec    " (ddv = @{dd(2)})
write/out *
write/out * to offset use  mg {d(1),{d(2)}  on TV-Guider!
write/out *
write/out *************************************************************


write/key ans/c/1/1 I

write/out   "       Further iteration             > I(terate)   [def] "
write/out   "       Verification with full field  > V(erification)    "
inq/key ans "       End                           > E(nd)            :"
ans = m$upper(ans)

if ans(1:1) .eq. "I" then
	write/out "----------------------------"
	write/out " Next aquisition exposures :"
	write/out "----------------------------"
	write/out

	loch1(1) = 2300-{x(1)}-40.
	loch1(2) = {x(2)}-100.
	loch1(3) = 2300-{x(1)}+40.
	loch1(4) = {x(2)}+40.
	set/format I4
	$ camera_chip1_x1 {loch1(1)}
	$ camera_chip1_y1 {loch1(2)}
	$ camera_chip1_x2 {loch1(3)}
	$ camera_chip1_y2 {loch1(4)}
	$ camera_chip1_hbin 1
	$ camera_chip1_vbin 1
	$ camera_dewar1_exptime {p5}
	$ camera_dewar1_exptype science
	$ camera_objectname "Hole 1 for {P1}"

	$ camera_dewar1_start wait

	set/format I4
	$ cp test0001.fits AQ{NR_AQ}e.fits
	do loop = 1  100
		newfil = m$exist("newdata")
		if {newfil} .eq. 0 then
			goto endloop5
		endif
		wait/secs 1.
	enddo
	endloop5:
!	stat/ima test0001 [<,<:>,>] ? 125,10000 FN
!	cut(1) = outputr(8)-2*outputr(4)
!	cut(2) = outputr(8)+5*outputr(4)
!	set/format f10.2
!	load test0001 chan=0 scale=4 cuts={cut(1)},{cut(2)} center=c,c
	load test0001 chan=0 scale=4 cuts=0,1000 center=c,c
	write/out "Set cursor on sky background for background determination ..."
	@@ backdet


        copy/dk test0001 CAHA.ADA.POSANG PA_C_old

	writ/out " Mark position of star 1 with cursor box ! "
        center/gauss cursor
	write/out OK ... continuing
        xss(1) = outputr(5)
        xss(2) = outputr(6)

	loch2(1) = 2300-{y(1)}-40.
	loch2(2) = {y(2)}-100.
	loch2(3) = 2300-{y(1)}+40.
	loch2(4) = {y(2)}+40.
	set/format I4
	$ camera_chip1_x1 {loch2(1)}
	$ camera_chip1_y1 {loch2(2)}
	$ camera_chip1_x2 {loch2(3)}
	$ camera_chip1_y2 {loch2(4)}
	$ camera_chip1_hbin 1
	$ camera_chip1_vbin 1
	$ camera_dewar1_exptime {p5}
	$ camera_dewar1_exptype science
	$ camera_objectname "Hole 2 for {P1}"

	$ camera_dewar1_start wait
	
	set/format I4
	$ cp test0001.fits AQ{NR_AQ}f.fits
	do loop = 1  100
		newfil = m$exist("newdata")
			if {newfil} .eq. 0 then
			goto endloop6
		endif
		wait/secs 1.
	enddo
	endloop6:
!	stat/ima test0001 [<,<:>,>] ? 125,10000 FN
!	cut(1) = outputr(8)-2*outputr(4)
!	cut(2) = outputr(8)+5*outputr(4)
!	set/format f10.2
!	load test0001 chan=0 scale=4 cuts={cut(1)},{cut(2)} center=c,c
	load test0001 chan=0 scale=4 cuts=0,1000 center=c,c
	write/out "Set cursor on sky background for background determination ..."
	@@ backdet


        write/out " Mark position of star 2 with cursor box ! "

        center/gauss cursor
        yss(1) = outputr(5)
        yss(2) = outputr(6)

	goto loop1
endif

if ans(1:1) .eq. "V" then

	write/out
	write/out "-----------------------------------------"
	write/out " Exposing full field for verification ..."
	write/out "-----------------------------------------"

	set/format I4
	$ camera_chip1_x1 2
	$ camera_chip1_y1 1
	$ camera_chip1_x2 2298
	$ camera_chip1_y2 2048
	$ camera_chip1_hbin 1
	$ camera_chip1_vbin 1
	$ camera_dewar1_exptime {P3}
	$ camera_dewar1_exptype science
	$ camera_objectname "Final AQ {P1}"

	$ camera_dewar1_start wait
	
	set/format I4
	$ cp test0001.fits AQ{NR_AQ}g.fits
	do loop = 1  100
		newfil = m$exist("newdata")
		if {newfil} .eq. 0 then
			goto endloop7
		endif
		wait/secs 1.
	enddo
	endloop7:

	loch1(1) = x(1)-15.
	loch1(2) = x(2)-15.
	loch1(3) = x(1)+15.
	loch1(4) = x(2)+15.
	set/format I4
	stat/ima test0001 [{loch1(1)},{loch1(2)}:{loch1(3)},{loch1(4)}] ? ? FN
	cut(1) = outputr(8)-2*outputr(4)
	cut(2) = outputr(8)+5*outputr(4)

	set/format f10.2
        load test0001 chan=0 scale=-4 cuts={cut(1)},{cut(2)} center=c,c

endif

if ans(1:1) .eq. "E" then

	goto end

end if

end:

! reset CCD and MOSCA


write/out
write/out "Get out of CCD TEST mode to save next frames !"
write/out "Choose your grism and start exposure of spectrum !"
write/out




$camera_expmode normal
return





exit:




