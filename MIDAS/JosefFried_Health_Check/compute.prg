! km_compute.prg
!
!  08/Jan/99

def/par p1 ZYX C " Enter OUTPUT filename >"
def/par p2 =
!-------------------------------------------------------------------------------
if P1 .eq. "ZYX" then
	writ/out " call:"
	writ/out "       COMPUTE/wfi OUT =    IN_1  operator  IN_2"
	writ/out " e.g.  COMPUTE/wfi result = wfi07481 - wfi07480"
	writ/out 
	return
endif
!-------------------------------------------------------------------------------

def/par p3 ? C " Enter first  operand >"
def/par p4 ? ? " operator: "
def/par p5 ? C " Enter second operand >"

comp/ima {p1}  = {p3} {p4}{p5}
comp/ima {p1}a = {p3}a{p4}{p5}a
comp/ima {p1}b = {p3}b{p4}{p5}b
comp/ima {p1}c = {p3}c{p4}{p5}c
comp/ima {p1}d = {p3}d{p4}{p5}d
comp/ima {p1}e = {p3}e{p4}{p5}e
comp/ima {p1}f = {p3}f{p4}{p5}f
comp/ima {p1}g = {p3}g{p4}{p5}g

return







