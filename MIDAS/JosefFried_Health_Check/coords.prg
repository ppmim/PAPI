!
! reads coordinates, converts to RA,DEC
!

cen/gauss ? cursor

defi/local r_a/r/1/1
defi/local dec/r/1/1

r_a = outputr(5)
dec = outputr(6)


define/local R_AH/i/1/1
define/local R_AM/i/1/1
define/local R_AS/r/1/1
define/local DECD/i/1/1
define/local DECM/i/1/1
define/local DECS/r/1/1

R_A = R_A/15.
R_AH = R_A
R_AM = (R_A-R_AH)*60.
R_AS = (R_A-R_AH)*3600.-R_AM*60.

DECD = DEC
DECM = (DEC-DECD)*60.
DECS = (DEC-DECD)*3600.-DECM*60.


set/for I3 f5.3
write/out RA = {R_AH} {R_AM} {R_AS}
write/out DEC = {DECD} {DECM} {DECS} 
