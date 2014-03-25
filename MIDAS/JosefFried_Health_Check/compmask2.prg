!
! compmask.prg
!
! JWF apr/may 98
!
! computes mask slits from table with alphas and deltas
!
! the ref.stars must have type=1 , the objects type=0
!
!



if "{p1}" .eq. "help"  then
	write/out
	write/out @@ compmask table_in AHR,AMIN,ASEC DDEG,DMIN,DSEC PA grism 
	write/out where
	write/out table_in = input table whith the columns 
	write/out "           :AHR,:AMIN,:ASEC,:DDEG,:DMIN,:DSEC,:type"
	write/out "           where type=0 for objects and type=1 for the 2 ref.stars"
	write/out AHR,AMIN,ASEC = RA of field center
	write/out DDEG,DMIN,DSEC = Dec of field center
	write/out PA = position angle
	write/out grism = grism name (e.g. green_500) used for calculating approximate
	write/out "        spectral coverage"
	write/out 
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
defi/par p5 green_500 ima "grism [def.=green_500] = ?"
defi/par p6 5.,-20. n "offsets of ref. holes [arcsec]=?"
defi/par p7 10. n "diameter of holes for stars [arcsec]=?"



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
defi/local offhole/r/1/2 {p6}

! mechanical center of mask in mm
defi/local maskcen/r/1/2 60.,-60.

! diameter of holes for stars 
defi/local dia/r/1/1 {p7}

set/format i2
cop/dk {p1}.tbl tblcontr/i/1/4 n
zi = {n(4)}
mrows = zi+2

write/out "  "
write/out " input table = {p1} contains {zi} rows"
write/out " alpha center = {p2} "
write/out " delta center = {p3} "
write/out " "
!
! check if columns exist
!
write/out special version for input from astromet
 
write/out check if input columns exist, # of ref stars ...

defi/local ex/i/1/1 0
comp/key ex(1) = M$EXISTC("{p1}",":dd")
if {ex} .eq. -1 then
	write/out column :ad does not exist -> exiting
	goto exit
endif
comp/key ex(1) = M$EXISTC("{p1}",":dd")
if {ex} .eq. -1 then
	write/out column :dd does not exist -> exiting
	goto exit
endif

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


! check grism name
if "{p5}" .eq. "blue_1000" then 
	dis = 1.229
	lamcen = 3923.
	goto start
endif

if "{p5}" .eq. "green_1000" then 
	dis = 1.312
	lamcen = 5314.
	goto start
endif


if "{p5}" .eq. "red_1000" then
	dis = 1.259
	lamcen = 6982.
	goto start
endif

if "{p5}" .eq. "blue_500" then
	dis = 2.92
	lamcen = 3807.
	goto start
endif

if "{p5}" .eq. "green_500" then
	dis = 2.883
	lamcen = 6982.
	goto start
endif

if "{p5}" .eq. "green_250" then
	dis = 6.00
	lamcen = 5086.
	goto start
endif

write/out grism {p5} unknown - check name given ---> exiting
goto exit


start:

! everything ok - continue


! field center in degrees
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
write/desc {p1}_mask.tbl alpcen/r/1/3 {p2}
write/desc {p1}_mask.tbl delcen/r/1/3 {p3}
write/desc {p1}_mask.tbl posang/r/1/1 {p4}
write/desc {p1}_mask.tbl offhole/r/1/2 {p6}
write/desc {p1}_mask.tbl dia/r/1/1 {p7}



! comp XY coordinates (see Montenbruck+Pfleger, Astronmie mit dem pers.Comp.)
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

name/col {p1}_mask :x_m  "mm" f10.3
name/col {p1}_mask :y_m  "mm" f10.3
name/col {p1}_mask :y_m_s  "mm" f10.3
name/col {p1}_mask :y_m_e  "mm" f10.3


! compute length of slits

comp/tab {p1}_mask :slen = :x_a_e-:x_a_s
name/col {p1}_mask :slen f8.1 "arcsec"

! compute spectral range

dis = {dis}/0.015

comp/tab {p1}_mask :lam1 = -(57.+(60.+:x_m))*0.27*{dis}+{lamcen}
comp/tab {p1}_mask :lam2 =  (57.+(60.+:x_m))*0.27*{dis}+{lamcen}
	

write/out 
sel/tab {p1}_mask :type.eq.0
name/col {p1}_mask :lam1 f6.0 "Angstroem"
name/col {p1}_mask :lam2 f6.0 "Angstroem"

write/out
write/out the result of the computation is:
write/out
read/tab {p1}_mask :ident,:R_A,:DEC,slen,lam1,lam2
write/out slen = length of slits [arcsec]
write/out lam1,lam2 = approx. wavelength range covered [Angstroem] by grism {p5}

write/out
write/out " use @@ plotmask {p1} to plot the mask"
write/out " use @@ mechmask {p1} to create CNC program"
write/out " output is stored in table  {p1}_mask "
write/out

assi/pri file {p1}.log
write/out
pri/tab {p1}_mask :ident,:R_A,:DEC,:mag,slen,lam1,lam2
assi/pri
write/out log-file -> {p1}.log

$more {p1}.log


exit:
