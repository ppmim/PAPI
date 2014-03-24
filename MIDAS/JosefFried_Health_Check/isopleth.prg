!
! isopleth.prg
!
!
! converts table to bdf
!
! jwf nov 2000
!


defi/par p1 car0093aobj ima "input table = ? "
defi/par p2 z ima "output frame = ? "

defi/par p3 1,2000,1,4000 n "x1,x2,y1,y2 = ? "
defi/par p4 1000,2000 n "npix_out_x,npix_out_y = ? "
defi/par p5 #4,#5,#15,#16


write/key in_a/c/1/20 {p1}
write/key out_a/c/1/20 {p2}

write/key inputd/d/1/10 {p3}
write/key inputi/i/1/10 {p4}
write/key inputc/c/1/20 {p5}

run isopleth 
