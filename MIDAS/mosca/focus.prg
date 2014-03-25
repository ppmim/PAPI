!
! focus.prg
!
! evaluate focus series
!
! HJR+JWF on a cloudy night oct 1996
!
! modifed jwf oct 97 oct 98 nov2000
!

if p1(1:4) .eq. "help" then
	write/out @@ focus [label of plot]
	write/out "         to evaluate a focus series"
	write/out You have to mark the images of an object with the cursor
	write/out center/gauss is used to determine the sigmas of the images
	write/out polynomes are fitted to x and y sigmas  to determine 
	write/out the best focus. 
	write/out
	return
endif


defi/local z/d/1/1 0.
set/format f10.2

defi/local fokus/r/1/2 0,0
defi/local focus/r/1/2 0.,0.


!
! firstcheck if key focser exists
!

defi/local n/i/1/1 0
n = m$existk("focser")
if n .eq. 0 then 
	write/out ERROR !
	write/out You have to take a focus series first (@@ autofocus)
	write/out exiting ...
	goto exit
endif


! 
! define device to be focussed
!
defi/local device/c/1/1 m

if {focser(1)} .gt. 40 then 
	device = "m"
	write/out focussing MOSCA ...
else
	device = "t"
	write/out focussing telescope ...
	write/out NOTE telecope focus has to be entered in microns!
endif







inquire/key z "First focus value  [{focser(1)}] : "
if z .ne. 0 then
		fokus(1) = z
	else
		fokus(1) = focser(1)
endif

z = 0.
inquire/key z "Focus step  [{focser(2)}] : "
if z .ne. 0 then
		fokus(2) = z
	else
		fokus(2) = focser(2)
endif







crea/grap
write/out
write/out "Mark with cursor box"
write/out "NOTE:the last focus step is  doubly separated"
write/out " "
 
clear/chan over
set/graph
set/format
center/gauss ? focus

write/out "okay ... determining best focus now!"





defi/local pixsize/r/1/1 15.             !pixsize in mue
defi/local inst_scal/r/1/1 0.273         !instrument magnification factor 
defi/local tel_scal/r/1/1/ 169.7         !telescope scale mue/arcsec
defi/local pix/r/1/2 0.,0. 


cop/dk {idimemc} step/d/1/2 step/d/1/2
cop/dk {idimemc} pixsize/r/1/1 pixsize/r/1/1 


comp/key pix(1) = {pixsize}/({tel_scal}*{inst_scal})
comp/key pix(2) = {pixsize}/({tel_scal}*{inst_scal})






write/key log/i/4/1 2

comp/tab focus :xfwhm = :xfwhm*{pix(1)}
comp/tab focus :yfwhm = :yfwhm*{pix(2)}
compute/tab focus :focus = (seq-1)*{fokus(2)}+{fokus(1)}
regression/poly focus :xfwhm :focus 2
create/tab fitfocus 3 100
save/regre focus coeffx
compute/regr focus :fitx = coeffx
compute/tab fitfocus :fitf = (seq-1)*({fokus(2)}/10)+({fokus(1)})
compute/tab fitfocus :fitx = {outputd(1)}+({outputd(2)})*:fitf+({outputd(3)})*:fitf**2

define/local fwhmx/d/1/1 0.
define/local focusx/d/1/1 0.
focusx = -{outputd(2)}/(2*{outputd(3)})
fwhmx = ({outputd(1)}-({outputd(2)})**2/(4*({outputd(3)})))


regression/poly focus :yfwhm :focus 2
save/regre focus coeffy
compute/regr focus :fity = coeffy
compute/tab fitfocus :fity = {outputd(1)}+({outputd(2)})*:fitf+({outputd(3)})*:fitf**2
define/local fwhmy/d/1/1 0.
define/local focusy/d/1/1 0.
focusy = -{outputd(2)}/(2*{outputd(3)})
fwhmy = ({outputd(1)}-({outputd(2)})**2/(4*({outputd(3)})))



stat/tab focus :xfwhm
defi/loc fwhmmin/r/1/1 0.
comp/key fwhmmin = outputr(1)
defi/loc fwhmmax/r/1/1 0.
comp/key fwhmmax = outputr(2)
stat/tab focus :yfwhm

if fwhmmin .gt. outputr(1) then
	fwhmmin = outputr(1)
endif
if fwhmmax .lt. outputr(2) then
	fwhmmax = outputr(2)
endif

fwhmmin = fwhmmin-0.1
fwhmmax = fwhmmax+0.1

stat/tab focus :focus
defi/loc focmin/r/1/1 0.
comp/key focmin = outputr(1)-{focus(2)}
defi/loc focmax/r/1/1 0.
comp/key focmax = outputr(2)+{focus(2)}

crea/tab dummy 2 2 NULL
crea/col dummy :focus " "
crea/col dummy :FWHM "arc sec"
dummy,1,1 = {focmin}
dummy,1,2 = {focmax}
dummy,2,1 = {fwhmmin}
dummy,2,2 = {fwhmmax}

set/grap stype=0 ltype=0 pmode=1 colour=1
plot/tab dummy #1 #2



set/grap stype=6 colour=2 ltype=0

over/tab focus :focus :xfwhm
set/grap ltype=1 stype=0
over/tab fitfocus :fitf :fitx
set/grap ltype=0 stype=11  colour=4
over/tab focus :focus :yfwhm
set/graph ltype=2 stype=0
over/tab fitfocus :fitf :fity

set/grap colour=1

defi/local fwhm/r/1/1 0.
comp/key focus = (focusx+focusy)*0.5
comp/key fwhm = M$sqrt(fwhmx*fwhmy)

write/key log/i/4/1 0
set/format i2 f8.3

write/out
write/out "Best focus = {focus} (in X = {focusx}, in Y = {focusy})"
write/out
write/out "Seeing = {fwhm} arcsec FWHM"
write/out


!label/graph {p1} 80,90,mm 0 1 0
label/graph "Focus = {focus}" 80,95,mm 0 1 0
label/graph "Seeing = {fwhm} arcsec" 80,90,mm 0 1 0
set/grap color=2
over/symbol 6 65,85,mm
label/grap "x-focus" 80,85,mm 0 1 0
set/grap color=4
over/symbol 11 65,82,mm
label/grap "y-focus" 80,82,mm 0 1 0



! focus telescope

if device .eq. "t" then
	defi/local move/c/1/1 n
	inquire/key move "focus telescope  [def]? [y/n] (def=n)"

	if move .eq. "y" then
		write/out --- focussing the telescope 
		set/format f6.3	
		$/disk-a/staff/TECS35/scripts/t3_afocus {focus}
		goto exit
	endif
endif


! focus MOSCA

if device .eq. "m" then
	wk move/c/1/1 n
	inquire/key move "Set focus zero point in MOSCA ? [y/n] (def=n)"
	if move .eq. "y" then
		inquire/key z "Focus value  [{focus}] : "
		if z .ne. 0 then
			focus(1) = z
		endif
		write/out --- focussing MOSCA 
		set/format f10.0	
		$/disk-a/staff/MOSCA/scripts/mosca_afocus {focus}
		$/disk-a/staff/MOSCA/scripts/mosca_set_focus_zero {focus}
		write/out	
		write/out MOSCA focus zero point set
		write/out	
	endif
endif


exit:

$rm dummy.tbl
return

