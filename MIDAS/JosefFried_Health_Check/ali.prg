!
! ali.prg
!
! align im1 to im2    result is stored in {p3}
!
! jwf aug 96
!


defi/par p1 ? ima "image = ? "
defi/par p2 ? ima "ref.image = ? "
defi/par p3 ? ima "result = ? "
defi/par p4 unit c "type of alignment=?"
load {p1}
write/out measure positions of >3 stars 
cen/gauss ? ima1  
   
load {p2}
write/out  measure positions of same stars 
cen/gauss ? ima2

align/ima ima1 ima2 {p4}
rebin/rotate {p1} {p3} KEYWORD {p2}

load {p3}
write/out image {p3} is aligned to image {p2}
