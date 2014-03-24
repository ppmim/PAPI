!
!	Creates a LAICA Mosaic.
!	We assume that flatfielding etc. have been done. Also that the files have been through 
!	astrometrie.prg allready. The filename root has to have eight characters.
!	We assume you know wich part of the sky you want.
!
!	JW 5.01
!
set/format i1
define/parameter p1 ? ? "Four input frames:"
IF p1 .eq. "help" THEN
   write/out @@ mosaiclaica input output step start,end
   write/out Each frame name has to have eight characters, no _chip#!
   write/out
   write/out example: @@ moscailaica 01day01,01day02,01day03,01day04 ngc0815 1.322e-4 241.2
   write/out step = degree/pixel default=6.25e-5
   return
ENDIF

define/parameter p2 ? ? "Output frame:"
define/parameter p3 6.25e-5  N "Step size in the output frame:"
define/parameter p4 ? ? "start(ra),start(dec),end(ra),end(dec)"

define/local step/d/1/2 -{p3},{p3}
define/local opos/d/1/2 
define/local i/i/1/1 
define/local names/c/1/35 {p1}
define/local npix/i/1/2
define/local alpha1/d/1/4
define/local delta1/d/1/4
define/local alpha2/d/1/4
define/local delta2/d/1/4
define/local alpha3/d/1/4
define/local delta3/d/1/4
define/local alpha4/d/1/4
define/local delta4/d/1/4
define/local start/d/1/4 'p4'
define/local end/d/1/2 {start(3)},{start(4)}


! Compute the rebinned images for alpha,delta.

cop/dk {names(1:8)}_1 O_POS/d/2/1 opos(2)
do i = 1 4
 @@ rebastro {names(1:8)}_{i} temp1_{i} {step(1)},{step(2)} {opos(2)}
 comp z = {names(1:8)}_{i}*0.+1.
 @@ rebastro z z1 {step(1)},{step(2)} {opos(2)}
 comp temp1_{i} = temp1_{i}/z1
enddo








cop/dk {names(10:17)}_1 O_POS/d/2/1 opos(2)
do i = 1 4
 @@ rebastro {names(10:17)}_{i} temp2_{i} {step(1)},{step(2)} {opos(2)}
 comp z = {names(10:17)}_{i}*0.+1.
 @@ rebastro z z1 {step(1)},{step(2)} {opos(2)}
 comp temp2_{i} = temp2_{i}/z1

enddo


cop/dk {names(19:26)}_1 O_POS/d/2/1 opos(2)
do i = 1 4
 @@ rebastro {names(19:26)}_{i} temp3_{i} {step(1)},{step(2)} {opos(2)}
 comp z = {names(19:26)}_{i}*0.+1.
 @@ rebastro z z1 {step(1)},{step(2)} {opos(2)}
 comp temp3_{i} = temp3_{i}/z1
enddo


cop/dk {names(28:35)}_1 O_POS/d/2/1 opos(2)
do i = 1 4
 @@ rebastro {names(28:35)}_{i} temp4_{i} {step(1)},{step(2)} {opos(2)}
 comp z = {names(28:35)}_{i}*0.+1.
 @@ rebastro z z1 {step(1)},{step(2)} {opos(2)}
 comp temp4_{i} = temp4_{i}/z1
enddo






! Create the appropriate output frame.
!$ rm ngc6240.bdf

opos(2) = {temp1_1,O_POS(2:2)}
step(1) = step(1) / m$cos({opos(2)})
npix(1) = M$ABS((end(1)-start(1)) / step(1))
npix(2) = M$ABS((end(2)-start(2)) / step(2))
$rm {p2}.bdf
write/out
write/out create output image ...
create/image {p2} 2,{npix(1)},{npix(2)} {start(1)},{start(2)},{step(1)},{step(2)}

write/descr {p2} CUNIT "            DEC         RA         "
write/descr {p2} IDENT "                                "


! Insert the subimages
write/out insert the subimages ...
do i = 1 4 
   insert/image temp1_{i} {p2} 
   insert/image temp2_{i} {p2} 
   insert/image temp3_{i} {p2} 
   insert/image temp4_{i} {p2} 
enddo


write/out all done - set the Identifier of your field correctly!
$date
exit:
