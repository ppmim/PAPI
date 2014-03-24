!
! analyses series of dark frames 
!
! JWF sep2011
!

! input files as indis/fits *rdoff*.fits - this stores files as toto0001...
! or like indis/fits *xyz*fits root=nov23 ...

defi/par p1 toto ima "root= ?"
defi/par p2 1,10 n "files start,end in sequence=?"
defi/par p3 sp c "window in frame [Q1,Q2,...]=?" 
defi/par p4 n c "skip stacking="

defi/loc root/c/1/20 {p1}
defi/loc mf/i/1/2 {p2}
defi/loc wind/c/1/20 {p3}

defi/loc area/c/1/30 " " all

if "{wind(1:2)}" .eq. "Q1" then
   area = "[@10,@10:@2030,@2030]"
endif


if "{wind(1:2)}" .eq. "Q2" then
   area = "[@2100,@10:@4080,@2100]"
endif


if "{wind(1:2)}" .eq. "Q3" then
   area = "[@2100,@2100:@4080,@4080]"
endif


if "{wind(1:2)}" .eq. "Q4" then
   area = "[@10,@2100:@2100,@4080]"
endif


if "{wind(1:2)}" .eq. "al" then
   area = "[@10,@10:@4080,@4080]"
endif

if "{wind(1:2)}" .eq. "sp" then
   area = "[@970,@3086:@1222,@3127]"
  area = "[@768,@3350:@1268,@3850]"
  area = "[@1912,@2200:@1950,@2230]"
  area = "[@1770,@2280:@1805,@2300]"
endif



defi/loc n/i/1/1 0
defi/loc m/i/1/1 0
defi/loc m1/i/1/1 0
defi/loc m2/i/1/1 0
defi/loc nrow/i/1/1 0


crea/icat ana.cat {p1}*.bdf
sho/icat ana.cat >NULL
defi/loc nf/i/1/1 {outputi(1)}

write/out 
wri/out files from {mf(1)} to {mf(2)} out of sequence with {nf} files 
write/out 


   
if "{p4}" .eq. "y" then
   goto analyse
endif

@@ stack ana {mf(1)},{mf(2)},1 mean sig

analyse:


write/out area used = {area}
stat/ima mean {area} >NULL
write/out mean level in stack =  {outputr(3)} +- {outputr(4)}     
stat/ima sig {area} >NULL
write/out mean sigma in stack =  {outputr(3)} +- {outputr(4)}  

! remove : in area to draw rect in display!

inputi(1) = m$index(area,":")
write/key area/c/{inputi(1)}/1 ","

ad0
clea/chan over
loa mean sc=f cuts=-20,20
draw/rect {area} F
ad1
loa sig sc=1 cuts=0,10
draw/rect {area} F