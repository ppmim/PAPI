!"auto-guide" computes the real position of a star 
! that is expected at P2,P3 in frame P1
! it returns the difference in 1/10 arcsec in x(Q1) und y(Q2)
! and return status of center/gauss (Q3)

DEFINE/PAR P1 ? I/A "Name of image file: "
DEFINE/PAR P2 ? N/A "x-value of expected star: "
DEFINE/PAR P3 ? N/A "y-value of expected star: "

DEFINE/LOCAL offset/i/1/2 0,0
DEFINE/LOCAL size/i/1/2 20,20

SET/MIDAS output=logonly
CENTER/GAUSS {P1},{P2},{P3} ? ? ? {size(1)},{size(2)}
SET/MIDAS output=yes	
offset(1) = M$NINT(({P2}-outputr(5))* 0.44942 * 10)
offset(2) = M$NINT(-({P3}-outputr(6))* 0.44942 * 10 )

IF outputr(11) .NE. 0 then
    WRITE/OUT star tracking not possible
ENDIF

return {offset(1)} {offset(2)} {outputr(11)}