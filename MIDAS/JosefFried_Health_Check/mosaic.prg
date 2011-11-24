! km_mosaic.prg
!
! NEW VERSION:  13/Jan/99
!

defi/par p1 ? ?    " input images >"
defi/par p2 mask ? " mask/nomask >"

!-------------------------------------------------------------------------------
if P1 .eq. "ZYX" then
	writ/out " call:"
	writ/out "       MOSAIC/wfi  frame(of which _,a,b,c,... must exist)  mask/nomask"
	writ/out " e.g.  MOSAIC/wfi  prep7542  ... will produce prep7542_all"
	writ/out 
	return
endif
!-------------------------------------------------------------------------------

defi/loc I/I/1/1 0

writ/out 
writ/out " create mosaic frame ..."

i = M$exist("mosaic0.bdf")

if I .eq. 1 then
	$cp mosaic0.bdf {p1}_all.bdf
else 
	create/ima {p1}_all 2,8574,8256
	writ/out " save empty mosaic frame ..."
	$cp {p1}_all.bdf mosaic0.bdf
endif
	
writ/out " insert images ..."
writ/out " a -> all"
insert/ima {p1}a {p1}_all 3,4128

writ/out " b -> all"
insert/ima {p1}b {p1}_all 2145,4125

writ/out " c -> all"
insert/ima {p1}c {p1}_all 4289,4124

writ/out " d -> all"
insert/ima {p1}d {p1}_all 6432,4126

writ/out " e -> all"
insert/ima {p1}e {p1}_all 6430,2

writ/out " f -> all"
insert/ima {p1}f {p1}_all 4285,0

writ/out " g -> all"
insert/ima {p1}g {p1}_all 2145,0

writ/out " h -> all"
insert/ima {p1}h {p1}_all 0,4

if p2(1:4) .eq. "mask" then
	comp/ima {p1}_all = {p1}_all * {actmask}
endif

copy/dd {p1}a *,3 {p1}_all

return





