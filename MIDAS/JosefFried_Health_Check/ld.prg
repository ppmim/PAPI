defi/par p1 z ima "input file=?"
defi/par p2 -1,-1  n "lhcuts=?"
defi/par p3 99 N "scale=?"

if "{p2}" .ne. "-1,-1" then
	wd {p1} lhcuts/r/1/2 {p2}
endif


if "{p3}" .eq. "99" then
	load/ima {p1} 
else
	load/ima {p1} scale={p3}
endif


