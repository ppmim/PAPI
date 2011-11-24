!
! check the astrometry
!


defi/par p1 ? ima "frame1 name=?"
defi/par p2 ? ima "frame2 name=?"

defi/par p3 1.5e-3 n "box width = ? "


defi/local box/r/1/1 {p3}

loa {p1} 

cen/gaus cursor  check 


loa/tab check :xcen :ycen ? ? 10
comp/tab check :xc1 = :xcen
comp/tab check :yc1 = :ycen

cen/gauss {p2},check check

comp/tab check :diff_r_a = (:xc1-:xcen)*3600.
comp/tab check :diff_dec = (:yc1-:ycen)*3600.
name/col check :diff_r_a "arcsec"
name/col check :diff_dec "arcsec"
sel/tab check :R_MAG1.gt.15.and.:R_mag1.lt.18.
write/out
write/out Errors in RA
stat/tab check :diff_r_a
write/out
write/out Errors in DEC
write/out
stat/tab check :diff_dec
write/out
set/grap pmode=1

set/grap xaxis=-2,2,0.5,.1 
plot/hist check :diff_r_a ? 0.05
label/grap  "R_A" 20,20
set/grap color=2
label/grap "DEC" 20,30
overplot/hist check :diff_dec
set/grap

name/col check :diff_r_a f10.6
name/col check :diff_dec f10.6
rt check :diff_r_a :diff_dec

set/grap xaxis=-1,1,0.5,.1
plot/hist check :diff_r_a ? 0.05
label/grap  "R_A" 20,20
set/grap color=2
label/grap "DEC" 20,30
overplot/hist check :diff_dec
set/grap

stat/tab check :diff_r_a
stat/tab check :diff_dec

