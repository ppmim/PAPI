!
! shuffle.prg
!
! shifts spectra in maskfiles
!
!
! jwf may/jul 98
!
!
defi/par p1 t0153 ima "inputfile=?"
defi/par p2 CNC5023_mask ima "table file=?"
defi/par p3 z ima outputfile=?

defi/local iwork/i/1/3 0,0,0
defi/local rwork/r/1/3 0,0,0
defi/local nrows/i/1/1 0
defi/local loop/i/1/1 0
defi/local xc/r/1/1 0.
defi/local yc/r/1/1 0.
defi/local slen/r/1/1 0.
defi/local scale/r/1/1 0.33
defi/local xs/r/1/1 0.
defi/local xe/r/1/1 0.
defi/local ys/r/1/1 0.
defi/local ystart/r/1/1 0.

defi/local sl/r/1/1 0.
defi/local xoff/r/1/1 1250.
defi/local yoff/r/1/1 1024.
defi/local nshift/i/1/1 0. 

! create outframe 

sel/tab {p2} :type.eq.0
nrows = {outputi(1)}
stat/tab {p2} :y_a
nshift = (outputi(2)-outputi(1))/2.
npix(2) = npix(2)+nshift

crea/ima {p3} 2,{npix(1)},{npix(2)} {start(1)},{start(2)},{step(1)},{step(2)}
cop/dk {p1} start/d/1/2 start



! extract subframes and insert them in out frame

set/format i2
do loop = 1 {nrows}
	xc = m$value({p2},:x_a,@{loop})
	yc = m$value({p2},:y_a,@{loop})
	slen = m$value({p2},:slen,@{loop})
	! convert arcsec to pixels
	xs = {yc}*{scale}+{xoff}
	ys = {xc}*{scale}+{yoff}
	sl = {slen}*{scale}
	xs = {ys} 
	xe = {xs}+{sl}
	ystart = start(2)+0.
!rk xc,yc,slen,xs,ys,sl,xs,xe
	write/out extracting subframe {p3}{loop} = {p1}[{xs},<:{xe},>]
	extra/ima z{loop} = {p1}[{xs},<:{xe},>]
	wd z{loop} start/d/2/1 {ystart}	
	insert/ima z{loop} z 
enddo

write/out done - now rotating image...
rota/counter_clock z {p3} 1

