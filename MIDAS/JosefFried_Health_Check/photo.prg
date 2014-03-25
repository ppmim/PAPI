!
! photo.prg  
!
defi/par p2 0.4068 n pixelscale=?
defi/local inr/r/1/1 0.
defi/local outr/r/1/1 0.
defi/local outs/r/1/1 0.
defi/local pix/r/1/1 0.
set/format f10.3

copy/dk {idimemc} CAL/d/1/1 inr
copy/dk {idimemc} pixscale/d/1/1 pix

!load {idimemc}

stat/ima ? cursor option=fn

comp/key outr = {inr}-2.5*M$log10({outputr(7)})
write/out "in window     {outputr(7)} counts   <=>  {outr} mag "


comp/key outr = {inr}+5.*M$log10({pix})-2.5*M$log10({outputr(8)})
write/out "median intensity {outputr(8)} counts  <=>  {outr} mag/sqas "
write/out "(image scale = {pix} arcsec/pix)"


