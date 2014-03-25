!
! adjust levels in overlap of 2 frames
!

defi/par p1 M33V0049_1 ima "frame 1?"
defi/par p2 M33V0050_1 ima "frame 2?"


cop/dk {p1} start start1
cop/dk {p1} step step1
cop/dk {p1} npix npix1

cop/dk {p2} start start2
cop/dk {p2} step step2
cop/dk {p2} npix npix2

rk start1,start2
rk step1,step2



! find overlap in x 
defi/local x11/r/1/1 0
defi/local x12/r/1/1 0
defi/local x21/r/1/1 0
defi/local x22/r/1/1 0

defi/local xs/r/1/1 0
defi/local xe/r/1/1 0

x11 = start1(1)
x12 = start1(1)+(npix1(1)-1)*step1(1)
x21 = start2(1)
x22 = start2(1)+(npix2(1)-1)*step2(1)

xs = x11
if x21 .gt. x11 then
   xs = x21
endif

xe = x12
if x22 .lt. x12 then 
   xe = x22
endif
   
set/for f10.3
rk x11,x12,x21,x22,xs,xe

! find overlap in y
defi/local y11/r/1/1 0
defi/local y12/r/1/1 0
defi/local y21/r/1/1 0
defi/local y22/r/1/1 0

defi/local ys/r/1/1 0
defi/local ye/r/1/1 0

y11 = start2(2)
y12 = start2(2)+(npix1(2)-1)*step1(2)
y21 = start2(2)
y22 = start2(2)+(npix2(2)-1)*step2(2)

ys = y11
if y21 .gt. y11 then
   ys = y21
endif

ye = y12
if y22 .lt. y12 then 
   ye = y22
endif
   

set/for f10.3
rk y11,y12,y21,y22,ys,ye

rk xs,ys,xe,ye





defi/local area1/c/1/20 [{xs},{ys}:{xe},{ye}]
ef
stat/ima {p1} [{xs},{ys}:{xe},{ye}]
eo
defi/local level/r/1/4 0,0,0,0

!stat/ima {p1} {area} >NULL
level(1) = outputr(3)
write/out Level 1 = {level(1)}

!stat/ima {p2} {area} >NULL
level(2) = outputr(3)

write/out Level 2 = {level(2)}