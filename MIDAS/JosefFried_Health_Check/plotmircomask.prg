!
! plotmask.prg
!
! plots MOS mask
!
! jwf mar/apr 1998
!

defi/par p1 ? ima "table=?"

!crea table with mask dimension if it does not exist
defi/local ex/i/1/1 0
ex = m$existc("mask_dim"," ")

if {ex} .eq. -2 then
	write/out create table with mech. dimensions of mask

	crea/tab mask_dim 2 4
	crea/col mask_dim :xm f10.3 "mm"
	crea/col mask_dim :ym f10.3 "mm"

	wri/tab mask_dim :xm @1 0.
	wri/tab mask_dim :ym @1 0.

	wri/tab mask_dim :xm @2 0.
	wri/tab mask_dim :ym @2 125.

	wri/tab mask_dim :xm @3 -105.
	wri/tab mask_dim :ym @3 135.

	wri/tab mask_dim :xm @4 -125.
	wri/tab mask_dim :ym @4 115.

	wri/tab mask_dim :xm @5 -120.
	wri/tab mask_dim :ym @5 0.

	wri/tab mask_dim :xm @6 0.
	wri/tab mask_dim :ym @6 0.

endif

set/grap stype=0 lwidth=3
plot/axes -140,10,50,10 -10,140,50,10 -100,-100 "x \  [mm]" "y \  [mm]"
set/grap stype=0 lwidth=3
over/tab mask_dim :xm :ym
over/line 2 -2.5,2.5 -2.5,117.5
over/line 2 -2.5,117.5 -117.5,117.5
over/line 2 -117.5,117.5 -117.5,2.5
over/line 2 -117.5,2.5 -2.5,2.5

set/grap stype=0 lwidth=3

defi/local yps/r/1/1 0.
defi/local ype/r/1/1 0.
defi/local xp/r/1/1 0.
defi/local loop/i/1/1 1
defi/local nslits/i/1/4 0,0



! draw slits

sel/tab {p1}_mask :type.eq.0
cop/tab {p1}_mask z
nslits = {outputi(1)}
set/format i2
write/out  {nslits}  slits

do loop = 1 {nslits}
	set/grap lwidth=3
	xp = m$value(z,:x_m,@{loop}) 
	yps = m$value(z,:y_m,@{loop}) 
!	ype = m$value(z,:y_m_e,@{loop}) 
!	over/line 1  {xp},{yps} {xp},{ype}
	ype = m$value(z,:y_m,@{loop}) 
	set/grap lwidth=1
!	over/symbol 12 {xp},{ype} 0.5
	over/symbol 2 {xp},{yps} 0.2
enddo
set/grap lwidth=2

! large holes for stars
write/out holes for stars

sel/tab {p1}_mask :type.eq.1
defi/local nstars/i/1/1 {outputi(1)}
cop/tab {p1}_mask z
loop = {nstars}-1
xp = m$value(z,:x_m,@{loop}) 
yps = m$value(z,:y_m,@{loop}) 
set/grap color=2
over/symbol 2 {xp},{yps} .5
over/symbol 2 {xp},{yps} .7

loop = {nstars}
xp = m$value(z,:x_m,@{loop}) 
yps = m$value(z,:y_m,@{loop}) 
over/symbol 2 {xp},{yps} .5
over/symbol 2 {xp},{yps} .7





! ref. holes
sel/tab {p1}_mask :type.eq.2
cop/tab {p1}_mask z
defi/local xo/r/1/1 0.
defi/local yo/r/1/1 0.
loop = {nstars}-1
xp = m$value(z,:x_m,@{loop})
yps = m$value(z,:y_m,@{loop})
over/symbol 2 {xp},{yps} .2

loop = {nstars}
xp = m$value(z,:x_m,@{loop}) 
yps = m$value(z,:y_m,@{loop}) 
over/symbol 2 {xp},{yps} .2
set/grap color=1

! draw windrose

defi/local x/r/1/4 0.,0.,-10.,10.
defi/local y/r/1/4 -10.,10.,0.,0.
defi/local xx/r/1/4 0,0,0,0
defi/local yy/r/1/4 0,0,0,0
defi/local n/i/1/1 0
defi/local xoff/r/1/1 45.
defi/local yoff/r/1/1 20.defi/local alpcen/r/1/3 0.,0.,0.
defi/local delcen/r/1/3 0.,0.,0.
cop/dk {p1}_mask.tbl alpcen alpcen/r/1/3
cop/dk {p1}_mask.tbl delcen delcen/r/1/3
cop/dk {p1}_mask.tbl posang pa
cop/dk {p1}_mask.tbl offhole offhole


do n = 1 4
	xx({n}) = xoff(1)+x({n})*m$cos({pa})+y({n})*m$sin({pa})
	yy({n}) = yoff(1)-x({n})*m$sin({pa})+y({n})*m$cos({pa})
enddo
set/format i1 f3.0
over/line 1 {xx(1)},{yy(1)} {xx(2)},{yy(2)}
over/line 1 {xx(3)},{yy(3)} {xx(4)},{yy(4)}
label/grap E {xx(2)},{yy(2)} 0 1 1 
label/grap N {xx(4)},{yy(4)} 0 1 1

label/grap " {p1}" 20.,120. 0. 1.5 1
label/grap " (as seen from below)" 20.,110. 0. .6 1
label/grap " alp cen : {alpcen(1)} {alpcen(2)}" 20.,80 0. .7 1
label/grap " del cen : {delcen(1)} {delcen(2)}" 20.,70. 0. .7 1
set/format f5.3
label/grap "{alpcen(3)}" 63.,80 0. .7 1
label/grap "{delcen(3)}" 63.,70. 0. .7 1
set/format f8.3
label/grap " Cass.fl.ang : {pa}" 20.,60. 0. .7 1
label/grap " hole offsets : {offhole(1)} {offhole(2)}" 20.,50. 0. .7 1


label/grap "{user}" -160.,-35. 0. .7 1
comp/key inputc = m$time()
label/grap "{inputc(1:30)}" -110.,-35. 0. .7 1
