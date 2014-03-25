! gets star coordinates with user interaction

DEFINE/PAR P1 ? I/A "Name of image file: "
DEFINE/LOCAL mean/I/1/1 0
DEFINE/LOCAL sdev/I/1/1 0
DEFINE/LOCAL cute/R/1/2 0,0
DEFINE/LOCAL answ/C/1/32 "Y"
DEFINE/LOCAL cent/I/1/2 1950,1950
DEFINE/LOCAL scal/I/1/1 2


STAT/IMA {P1} option=FNYN
mean = outputr(3)
sdev = outputr(4)
cute(1) = mean - (6 * sdev)
cute(2) = mean + (6 * sdev)

LOOP:

LOAD/IMA {P1} cuts={cute(1)},{cute(2)} center={cent(1)},{cent(2)} sc={scal}
 $play -q $GEIRS_DIR/SOUNDS/doorbell.au
INQUIRE/KEY answ "Ready to choose star(y) or change cuts(n)? " 
answ = m$upper(answ)
IF answ(1:1) .EQ. "N" THEN
 $play -q $GEIRS_DIR/SOUNDS/doorbell.au
	INQUIRE/KEY cute "lower,upper cut? "
	INQUIRE/KEY cent "center (x,y)?"
	INQUIRE/KEY scal "scale? "
	answ = "Y"
	GOTO LOOP
ENDIF


WRITE/OUT "Please choose star...(use arrow key to change cursor size)"
SET/MIDAS output=logonly
CENTER/GAUSS CURSOR ? ? 2,1,1 
SET/MIDAS output=yes

return {outputr(5)} {outputr(6)}



