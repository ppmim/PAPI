!
! compmask.prg
!
! JWF apr/may 98, jan 99
! 
! computes mask slits from table with alphas and deltas
!
! the ref.stars must have type=1 , the objects type=0
!
!


if "{p1}" .eq. "help"  then
	write/out
	write/out @@ compmask table_in AHR,AMIN,ASEC DDEG,DMIN,DSEC PA masktype 
	write/out "	[grism] [hole_off] [hole_diam]"
	write/out
	write/out where
	write/out table_in = input table whith the columns 
	write/out "           :AHR,:AMIN,:ASEC,:DDEG,:DMIN,:DSEC,:type"
	write/out "           where type=0 for objects and type=1 for the 2 ref.stars"
	write/out AHR,AMIN,ASEC = RA of field center
	write/out DDEG,DMIN,DSEC = Dec of field center
	write/out PA = cass.flange angle (note: the slits are at PA=90 deg for 
	write/out "       cass.flange angle=0 deg!)
	write/out masktype = slit or hole
	write/out grism = grism name [def = green_500] used for calculating approximate
	write/out "        spectral coverage"
	write/out hole_off = offset of holes [def=5,-20] arcsec from holes for ref.stars
	write/out
	write/out hole_diam = diameter of holes for ref.stars [def = 10 (arcsecs)]
	write/out 
	write/out hints for mask-makers
	write/out =====================
	write/out 
	write/out 1.) most important is good astrometry (error < 0.3 ")
	write/out "    all coordinates must be in the same equinox!" 
	write/out 2.) a position angle not =0 may help to get more objects per mask
	 

	goto exit
endif

defi/par p1 ? ima "table=?"
defi/par p2 ? n "alpha cen=? [hh,mm,ss.s]"
defi/par p3 ? n "delta cen=?  [dd,mm,ss.s]"
defi/par p4 ? n "pos.angle=? [deg]"
defi/par p5 ? ima "masktype (slit/hole) = ? " 
defi/par p6 green_500 ima "grism [def.=green_500] = ?"
defi/par p7 5.,-20. n "offsets of ref. holes [def = 5.,-20. (arcsec)]=?"
defi/par p8 10. n "diameter of holes for stars [def = 10. (arcsec)]=?"

defi/local masktype/c/1/4 {p5}

if "{masktype(1:4)}" .ne. "slit" .and.  "{masktype(1:4)}" .ne. "hole" then
	write/out ERROR: wrong masktype entered -> exiting
	goto exit
endif


defi/local mrows/i/1/1 0
defi/local typem/i/1/1 0
defi/local typel/i/1/1 0

defi/local zi/i/1/1 0
defi/local zr/r/1/2 0.
defi/local mslit/i/1/1 0
defi/local m/i/1/1 0
defi/local n/i/1/4 0,0,0,0
defi/local l/i/1/1 0
defi/local spa/d/1/1 0.
defi/local cpa/d/1/1 0.
defi/local sdel0/d/1/1 0.
defi/local cdel0/d/1/1 0.
defi/local al0/r/1/1 0.
defi/local del0/r/1/1 0.
defi/local z/r/1/3 0,0,0
defi/local text/c/1/20 " "
defi/local xm/r/1/1 0.
defi/local xl/r/1/1 0.

defi/local len/r/1/1 0.
defi/local dis/r/1/1 0.
defi/local lamcen/r/1/1 0.

!tel focal length in mm
defi/local f/r/1/1 35003.1376

! offsets for holes in arcsec
defi/local offhole/r/1/2 {p7}

! mechanical center of mask in mm
defi/local maskcen/r/1/2 60.,-60.

! diameter of holes for stars 
defi/local dia/r/1/1 {p8}

set/format i2
cop/dk {p1}.tbl tblcontr/i/1/4 n
zi = {n(4)}
mrows = zi+2

write/out "  "
write/out " input table = {p1} contains {zi} rows"
write/out " alpha center = {p2} "
write/out " delta center = {p3} "
write/out " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!  check input
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! check column :type

defi/local ex/i/1/1 0
comp/key ex(1) = M$EXISTC("{p1}",":type")
if {ex} .eq. -1 then
	write/out column :type does not exist -> exiting
	goto exit
endif

! check column :type
comp/key ex(1) = M$EXISTC("{p1}",":type")
if {ex} .eq. -1 then
	write/out column :type does not exist -> exiting
	goto exit
endif

!check number of ref.stars
sel/tab {p1} :type.eq.1
if outputi(1) .ne. 2 then
	write/out # of ref.stars unequal 2 -> exiting
	goto exit
endif
sel/tab {p1} all



! check if columns R_A and DEC exist
defi/local exa/i/1/1 0
defi/local exb/i/1/1 0
exa(1) = M$EXISTC("{p1}",":R_A")
exb(1) = M$EXISTC("{p1}",":DEC")


! columns do not exist -> ahr,amin,asec,ddeg,dmin,dsec required

if {exa} .eq. -1 .or. {exb} .eq. -1 then
	write/out INFO: input columns R_A,DEC do not exist-use ahr,amin,asec,ddeg,dmin,dsec
	comp/key ex(1) = M$EXISTC("{p1}",":AHR")
	if {ex} .eq. -1 then
     	   write/out column :AHR does not exist -> exiting
     	   goto exit
	endif
	comp/key ex(1) = M$EXISTC("{p1}",":AMIN")
	if {ex} .eq. -1 then
 	       write/out column :AMIN does not exist -> exiting
 	       goto exit
	endif
	comp/key ex(1) = M$EXISTC("{p1}",":ASEC")
	if {ex} .eq. -1 then
 	       write/out column :ASEC does not exist -> exiting
 	       goto exit
	endif
	
	comp/key ex(1) = M$EXISTC("{p1}",":DDEG")
	if {ex} .eq. -1 then
	        write/out column :AHR does not exist -> exiting
 		goto exit
	endif
	comp/key ex(1) = M$EXISTC("{p1}",":DMIN")
	if {ex} .eq. -1 then
	        write/out column :DMIN does not exist -> exiting
	        goto exit
	endif
	comp/key ex(1) = M$EXISTC("{p1}",":DSEC")
	if {ex} .eq. -1 then
        	write/out column :DSEC does not exist -> exiting
 	       goto exit
	endif
endif


! check grism name
if "{p6}" .eq. "blue_1000" then 
	dis = 1.229
	lamcen = 3923.
	goto start
endif

if "{p6}" .eq. "green_1000" then 
	dis = 1.312
	lamcen = 5314.
	goto start
endif


if "{p6}" .eq. "red_1000" then
	dis = 1.259
	lamcen = 6982.
	goto start
endif

if "{p6}" .eq. "blue_500" then
	dis = 2.92
	lamcen = 3807.
	goto start
endif

if "{p6}" .eq. "green_500" then
	dis = 2.883
	lamcen = 6982.
	goto start
endif

if "{p6}" .eq. "green_250" then
	dis = 6.00
	lamcen = 5086.
	goto start
endif

write/out grism {p6} unknown - check name given ---> exiting
goto exit





start:

! everything ok - continue
! field center -> decimal degrees
set/format f15.6
write/key z/r/1/3 {p2}
al0 = (z(1)+z(2)/60.+z(3)/3600.)*15.

write/key z/r/1/3 {p3}
if z(1) .lt. 0 then
		del0 = z(1)-z(2)/60.-z(3)/3600.
	else
		del0 = z(1)+z(2)/60.+z(3)/3600.
endif
sdel0 = m$sin(del0)
cdel0 = m$cos(del0)



write/out create output table and compute generalized coordinates...
sel/tab {p1} all
cop/tab  {p1} {p1}_mask
crea/col {p1}_mask :delalp R*8
crea/col {p1}_mask :denom R*8


! create columns ad and dd if not given
if {exa} .eq. -1 .or. {exb} .eq. -1 then
	crea/col {p1}_mask :ad R*8
	crea/col {p1}_mask :dd R*8
 	! convert alpha,delta -> degrees
	comp/tab {p1}_mask :ad = (:AHR+:AMIN/60.+:ASEC/3600.)*15.
	comp/tab {p1}_mask :dd = :DDEG+:DMIN/60.+:DSEC/3600.
	! convert to -:dd if necessary
	
endif 


write/desc {p1}_mask.tbl alpcen/r/1/3 {p2}
write/desc {p1}_mask.tbl delcen/r/1/3 {p3}
write/desc {p1}_mask.tbl posang/r/1/1 {p4}
write/desc {p1}_mask.tbl offhole/r/1/2 {p7}
write/desc {p1}_mask.tbl dia/r/1/1 {p8}



! comp XY coordinates (see Montenbruck+Pfleger, Astronomie mit dem pers.Comp.)
!  X1||alpha  Y1||delta
comp/tab {p1}_mask :delalp = :ad-{al0}
comp/tab {p1}_mask :denom = {cdel0}*cos(:dd)*cos(:delalp)+{sdel0}*sin(:dd)
comp/tab {p1}_mask :X1 = -cos(:dd)*sin(:delalp)/:denom
comp/tab {p1}_mask :Y1 = -({sdel0}*cos(:dd)*cos(:delalp)-{cdel0}*sin(:dd))/:denom

		
! rotate coordinate system
write/out rotate coordinate system by {p4} degrees...
spa = m$sin({p4})
cpa = m$cos({p4})
comp/tab {p1}_mask :X = :X1*{cpa}+:Y1*{spa}
comp/tab {p1}_mask :Y = -:X1*{spa}+:Y1*{cpa}



! create rows for reference holes, write data

set/format i2
do m = 1 {zi}
	typem = m$value({p1}_mask,:type,@{m})
	if {typem} .eq. 1 then
		zi = zi+1
!		write/out object at row {m} = star -> create ref. hole
		zr(1) = m$value({p1}_mask,:X,@{m})
		zr(2) = m$value({p1}_mask,:Y,@{m})
		zr(1) = zr(1)-offhole(1)*4.848e-6
		zr(2) = zr(2)+offhole(2)*4.848e-6
		l = zi-1
		crea/row {p1}_mask @{l} 1
		l = m+1
		{p1}_mask,:X,@{zi} = zr(1)
		{p1}_mask,:Y,@{zi} = zr(2)
		write/tab {p1}_mask :type @{zi} 2
!		text = m$value({p1}_mask,:number,@{m})
!		write/tab {p1}_mask :number @{zi} -{text(1:5)}
	endif
enddo
		
write/out "Info: reference holes are offset by {offhole(1)},{offhole(2)} [arcsec] from stars"



! sort table. compute lengths of slits
write/out compute lengths of slits...
sort/tab {p1}_mask :X
crea/col {p1}_mask :x_e
crea/col {p1}_mask :x_s

len = dia*0.5*4.848e-6

set/format i2
zi = mrows-1


! :x_s(1)  and :x_e(last)  = 20 arcsec from object
do m = 1 zi
	typem = m$value({p1}_mask,:type,@{m})
	xm = m$value({p1}_mask,:X,@{m})
	l = m+1
	typel = m$value({p1}_mask,:type,@{l})
	xl = m$value({p1}_mask,:X,@{l})
	n = m-1

	if {typem} .eq. 0 .and. {typel} .eq. 0 then
		if {m} .eq. 1 then
			zr(1) = xm-20.*4.848e-6
		else
			zr(1) =  m$value({p1}_mask,:x_e,@{n})   
		endif
		{p1}_mask,:x_s,@{m} = zr(1)
		zr(1) = (xm+xl)*0.5
		{p1}_mask,:x_e,@{m} = zr(1)
	endif

	if {typem} .eq. 0 .and. {typel} .eq. 1 then
		if {m} .eq. 1 then
			zr(1) = xm-20.*4.848e-6
		else
			zr(1) =  m$value({p1}_mask,:x_e,@{n})   
		endif
		{p1}_mask,:x_s,@{m} = zr(1)
		zr(1) = xl-len
		{p1}_mask,:x_e,@{m} = zr(1)
	endif

	if {typem} .eq. 0 .and. {typel} .eq. 2 then
		if {m} .eq. 1 then
			zr(1) = xm-20.*4.848e-6
		else
			zr(1) =  m$value({p1}_mask,:x_e,@{n})   
		endif
		{p1}_mask,:x_s,@{m} = zr(1)
		zr(1) = xl-4.848e-6
		{p1}_mask,:x_e,@{m} = zr(1)
	endif

	if {typem} .eq. 1  then
		{p1}_mask,:x_s,@{m} = xm-len
		{p1}_mask,:x_e,@{m} = xm+len
	endif
	
	if {typem} .eq. 2  then
		{p1}_mask,:x_s,@{m} = xm-4.848e-6
		{p1}_mask,:x_e,@{m} = xm+4.848e-6
	endif
enddo


! write last data point	
typem = m$value({p1}_mask,:type,@{mrows})
if {typem} .eq. 1  then
	{p1}_mask,:x_s,@{mrows} = xm-len
	{p1}_mask,:x_e,@{mrows} = xm+len
endif
if {typem} .eq. 0  then
	n = mrows-1
	zr(1) =  m$value({p1}_mask,:x_e,@{n})
	{p1}_mask,:x_s,@{mrows} = zr(1)
        zr(1) =  m$value({p1}_mask,:x,@{mrows})+20.*4.848e-6
	{p1}_mask,:x_e,@{mrows} = zr(1)
endif




! compute astronomical x/y positions N -> +y   E -> +x 
! [x_a ...] = arcsec
write/out compute astronomical x/y positions...
comp/tab {p1}_mask :x_a = :X*206264.8
comp/tab {p1}_mask :y_a = :Y*206264.8
comp/tab {p1}_mask :x_a_s = :x_s*206264.8
comp/tab {p1}_mask :x_a_e = :x_e*206264.8
name/col {p1}_mask :x_a  "arcsec" f10.3
name/col {p1}_mask :y_a  "arcsec" f10.3
name/col {p1}_mask :x_a_s  "arcsec" f10.3
name/col {p1}_mask :x_a_e  "arcsec" f10.3


! compute mechanical x/y positions N -> +x_m   E -> +y_m 
write/out  compute  mechanical x/y positions... 
comp/tab {p1}_mask :x_m = :Y*{f}+{maskcen(2)}
comp/tab {p1}_mask :y_m = -:X*{f}+{maskcen(1)}
comp/tab {p1}_mask :y_m_s = -:x_s*{f}+{maskcen(1)}
comp/tab {p1}_mask :y_m_e = -:x_e*{f}+{maskcen(1)}

! set slitlength = 0 if masktype = hole
if "{masktype(1:4)}" .eq. "hole" then
	comp/tab {p1}_mask :x_a_s = :x_a
	comp/tab {p1}_mask :x_a_e = :x_a
	comp/tab {p1}_mask :y_m_s = :y_m
	comp/tab {p1}_mask :y_m_e = :y_m
endif


name/col {p1}_mask :x_m  "mm" f10.3
name/col {p1}_mask :y_m  "mm" f10.3
name/col {p1}_mask :y_m_s  "mm" f10.3
name/col {p1}_mask :y_m_e  "mm" f10.3


! compute length of slits
comp/tab {p1}_mask :slen = :x_a_e-:x_a_s
name/col {p1}_mask :slen f8.1 "arcsec"







! compute spectral range
dis = -{dis}/0.015*0.27
comp/tab {p1}_mask :lam1 =  {lamcen}+57.*{dis}+(60.+:x_m)*{dis}
comp/tab {p1}_mask :lam2 =  {lamcen}-57.*{dis}+(60.+:x_m)*{dis}
	

write/out 
sel/tab {p1}_mask :type.eq.0
name/col {p1}_mask :lam1 f6.0 "Angstroem"
name/col {p1}_mask :lam2 f6.0 "Angstroem"

write/out
write/out the result of the computation is:
write/out
sel/tab {p1}_mask all
read/tab {p1}_mask :ident,:R_A,:DEC,:slen,lam1,lam2,type
write/out slen = length of slits [arcsec]
write/out lam1,lam2 = approx. wavelength range covered [Angstroem] by grism {p6}




assi/pri file {p1}.log
write/out
pri/tab {p1}_mask :ident,:R_A,:DEC,:mag,:slen,lam1,lam2,type
assi/pri
write/out log-file -> {p1}.log
write/out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! here comes the plotting stuff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



write/out start plotting...
write/out
!crea table with mask dimension if it does not exist

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
if "{masktype(1:4)}" .eq. "slit" then 
	write/out  {nslits}  slits
endif
if "{masktype(1:4)}" .eq. "hole" then 
	write/out  {nslits}  holes
endif



if "{masktype(1:4)}" .eq. "slit" then
	do loop = 1 {nslits}
		set/grap lwidth=3
		xp = m$value(z,:x_m,@{loop}) 
		yps = m$value(z,:y_m_s,@{loop}) 
		ype = m$value(z,:y_m_e,@{loop}) 
		over/line 1  {xp},{yps} {xp},{ype}
		ype = m$value(z,:y_m,@{loop}) 
		set/grap lwidth=1
		over/symbol 12 {xp},{ype} 0.5
	enddo
endif


if "{masktype(1:4)}" .eq. "hole" then
	do loop = 1 {nslits}
		set/grap lwidth=3
		xp = m$value(z,:x_m,@{loop}) 
		ype = m$value(z,:y_m,@{loop}) 
		over/symbol 2 {xp},{ype} 0.5
	enddo
endif






set/grap lwidth=2 color=2

! large holes for stars
write/out holes for stars

sel/tab {p1}_mask :type.eq.1
defi/local nstars/i/1/1 {outputi(1)}
cop/tab {p1}_mask z
loop = {nstars}-1
xp = m$value(z,:x_m,@{loop}) 
yps = m$value(z,:y_m,@{loop}) 
over/symbol 2 {xp},{yps} .5

loop = {nstars}
xp = m$value(z,:x_m,@{loop}) 
yps = m$value(z,:y_m,@{loop}) 
over/symbol 2 {xp},{yps} .5

set/grap color=1



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


! draw windrose

defi/local x/r/1/4 0.,0.,-10.,10.
defi/local y/r/1/4 -10.,10.,0.,0.
defi/local xx/r/1/4 0,0,0,0
defi/local yy/r/1/4 0,0,0,0
defi/local xoff/r/1/1 45.
defi/local yoff/r/1/1 20.defi/local alpcen/r/1/3 0.,0.,0.
defi/local delcen/r/1/3 0.,0.,0.
cop/dk {p1}_mask.tbl alpcen alpcen/r/1/3
cop/dk {p1}_mask.tbl delcen delcen/r/1/3
cop/dk {p1}_mask.tbl posang pa
cop/dk {p1}_mask.tbl offhole offhole

set/format i1 f3.0
do n = 1 4
	xx({n}) = xoff(1)+x({n})*m$cos({pa})+y({n})*m$sin({pa})
	yy({n}) = yoff(1)-x({n})*m$sin({pa})+y({n})*m$cos({pa})
enddo

over/line 1 {xx(1)},{yy(1)} {xx(2)},{yy(2)}
over/line 1 {xx(3)},{yy(3)} {xx(4)},{yy(4)}
label/grap E {xx(2)},{yy(2)} 0 1 1 
label/grap N {xx(4)},{yy(4)} 0 1 1

over/symbol 16 -80,120 1.5
label/grap "\lambda" -83.,123.5 0. .7

label/grap " mask {p1}" 20.,120. 0. 1.5 1
label/grap " (as seen from below)" 20.,110. 0. .6 1
label/grap " alp cen : {alpcen(1)} {alpcen(2)}" 20.,80 0. .7 1
label/grap " del cen : {delcen(1)} {delcen(2)}" 20.,70. 0. .7 1
set/format f5.3
label/grap "{alpcen(3)}" 63.,80 0. .7 1
label/grap "{delcen(3)}" 63.,70. 0. .7 1
set/format f8.3
label/grap " Cassfl.angle : {pa}" 20.,60. 0. .7 1
label/grap " hole offsets : {offhole(1)} {offhole(2)}" 20.,50. 0. .7 1


label/grap "{user}" -160.,-35. 0. .7 1
comp/key inputc = m$time()
label/grap "{inputc(1:30)}" -110.,-35. 0. .7 1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! here comes the CNC stuff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
abc:
write/out
n = 0
inquire/key n "enter 4 digit number to create CNC machine code or CR to exit :"

if {n} .le. 1000  then
	 goto exit 
endif
if {n} .gt. 9999  then
	write/out ERROR: file number must be between 1000 - 9999
	goto abc 
endif



defi/local numb/i/1/1 10

defi/local file/i/1/2 0,0

! used to correct length of slits for min. overlap
defi/local fras/r/1/1 0.25
if "{masktype(1:4)}" .eq. "slit" then
	fras = fras/4.
endif
if "{masktype(1:4)}" .eq. "hole" then
	fras = 0.
endif



open/file {n} write file

write/file {file} %{n} ( -D11-000000- -06.05.98- Maske{n})
write/file {file} ?
write/file {file} 0000
write/file {file} %{p2}*T

!define tools
write/file {file} T1 I/M:M R0.5 L0 FE300 FZ100 S2000 (Stichel 0.93)
write/file {file} T2 I/M:M R0.8 L11 FE300 FZ100 S1500 (Bohrer 1.7)
write/file {file} T3 I/M:M R0.2 L26.356 FE300 FZ100 S2000 (Bohrer 0.35)

write/file {file} ?
write/file {file} 0000
write/file {file} " "
write/file {file} %{n}*% ( -D11-000000- -06.05.98- Maske{n})


write/file {file} N10 G0 Z100 G17 T1 M3 M62

!  slits

sel/tab {p1}_mask :type.eq.0
cop/tab {p1}_mask z
nslits = {outputi(1)}
set/format i2 f6.3
write/out  {nslits}  slits

do loop = 1 {nslits}
	numb = {numb}+10
	xp = m$value(z,:x_m,@{loop}) 
	yps = m$value(z,:y_m_s,@{loop})-fras
	ype = m$value(z,:y_m_e,@{loop})+fras
	write/file {file} N{numb} G0 X{xp} Y{ype}
        numb = {numb}+10
	write/file {file} N{numb} G1 Z0
	numb = numb+10
	write/file {file} N{numb} G1 X{xp} Y{yps}
	numb = {numb}+10
	write/file {file} N{numb} G1 X{xp} Y{ype}
	numb = numb+10
	write/file {file} N{numb} G0 Z5
enddo

! Bohre Loecher an

sel/tab {p1}_mask :type.ne.0
cop/tab {p1}_mask z
nslits = {outputi(1)}
set/format i2 f6.3
write/out  {nslits}  holes

do loop = 1 {nslits}
	numb = {numb}+10
	xp = m$value(z,:x_m,@{loop}) 
	yps = m$value(z,:y_m,@{loop}) 
	write/file {file} N{numb} G0 X{xp} Y{yps}
        numb = {numb}+10
	write/file {file} N{numb} G1 Z0.3
	numb = numb+10
	write/file {file} N{numb} G0 Z5
enddo







! large holes for stars
! change tools
numb = numb+10
write/file {file} N{numb} G0 Z100
numb = numb+10
write/file {file} N{numb} M6
numb = numb+10
write/file {file} N{numb} G17 T2 M3 M62 [1.7]

sel/tab {p1}_mask :type.eq.1
nstars = {outputi(1)}
write/out large holes for {nstars} stars
cop/tab {p1}_mask z
do loop = 1 {nstars}
	xp = m$value(z,:x_m,@{loop}) 
	yps = m$value(z,:y_m,@{loop}) 
	numb = numb+10
	write/file {file} N{numb} G0 X{xp} Y{yps}
	numb = numb+10
	write/file {file} N{numb} G0 Z2
	numb = numb+10
	write/file {file} N{numb} G1 Z-2
	numb = numb+10
	write/file {file} N{numb} G0 Z10
enddo

! ref. holes
! change tools

numb = numb+10
write/file {file} N{numb} G0 Z100
numb = numb+10
write/file {file} N{numb} M6
numb = numb+10
write/file {file} N{numb} G17 T3 M3 M62 [0.35]


sel/tab {p1}_mask :type.eq.2
nstars = {outputi(1)}
write/out {nstars}  ref.holes...
cop/tab {p1}_mask z
xo = 0.
yo = 0.

do loop = 1 {nstars}
	xp = m$value(z,:x_m,@{loop}) 
	yps = m$value(z,:y_m,@{loop}) 
	numb = numb+10
	write/file {file} N{numb} G0 X{xp} Y{yps}
	numb = numb+10
	write/file {file} N{numb} G0 Z2
	numb = numb+10
	write/file {file} N{numb} G1 Z-0.3
	numb = numb+10
	write/file {file} N{numb} G0 Z10
enddo

numb = numb+10
write/file {file} N{numb} G0 Z100 M30
write/file {file} ?
write/file {file} A44C
close/file {file}
write/out
write/out all done !
write/out CNC code -> file {n}
write/out 


exit:
