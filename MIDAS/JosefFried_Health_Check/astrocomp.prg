! 
!	Computes RA and DEC with the plate solution form astro/trans.
!
!	JW 03.01
set/format F15.10,F15.13
define/parameter p1 ? ? "Image:"
define/parameter p2 ? ? "Pixel:"

define/local aldel0/d/1/2
define/local i/i/1/1
define/local cx/D/1/1
define/local cy/D/1/1
define/local bx/D/1/9
define/local by/d/1/9
define/local cen/d/1/2 'p2'
define/local Xi/d/1/1
define/local Eta/d/1/1
define/local alpha/d/1/1
define/local delta/d/1/1

copy/dk {p1} cx cx
copy/dk {p1} cy cy
copy/dk {p1} bx bx
copy/dk {p1} by by
copy/dk {p1} al_de0 aldel0

cx = cx/57.2957795
cy = cy/57.2957795

DO i = 1 9
bx({i}) = bx({i})/57.2957795
by({i}) = by({i})/57.2957795
ENDDO

compute/key cen(1) = cen(1)/1000
compute/key cen(2) = cen(2)/1000

compute/key Xi = cx + bx(1)*cen(1) 
compute/key Xi = Xi + bx(2)*{cen(2)} + bx(3)*cen(1)*{cen(2)} + bx(4)*cen(1)**2
compute/key Xi = Xi + bx(5)*{cen(2)}**2 + bx(6)*cen(1)**3 + bx(7)*{cen(2)}**3
compute/key Xi = Xi + bx(8)*cen(1)*{cen(2)}**2 + bx(9)*cen(1)**2*{cen(2)}
compute/key Xi = - Xi

compute/key Eta = cy + by(1)*cen(1) 
compute/key Eta = Eta + by(2)*{cen(2)} + by(3)*cen(1)*{cen(2)} + by(4)*cen(1)**2
compute/key Eta = Eta + by(5)*{cen(2)}**2 + by(6)*cen(1)**3 + by(7)*{cen(2)}**3
compute/key Eta = Eta + by(8)*cen(1)*{cen(2)}**2 + by(9)*cen(1)**2*{cen(2)}

aldel0(1) = aldel0(1) * 57.2957795
aldel0(2) = aldel0(2) * 57.2957795
compute/key alpha = aldel0(1)+m$atan(-Xi/(m$cos(aldel0(2))-Eta*m$sin(aldel0(2))))
compute/key delta = m$asin((m$sin(aldel0(2))+Eta*m$cos(aldel0(2)))/(1+Xi**2+Eta**2)**0.5) 

return {alpha} {delta}




