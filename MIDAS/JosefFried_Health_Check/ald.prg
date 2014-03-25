defi/par p1 z ima "input file=?"
defi/par p2 0 N "channel=?"
defi/par p3 -2 N "scale=?"
stat/ima 'p1' [@10,@10:@1000,@1000] option=fs 
COPY/DK 'P1' statistic/R/8/2 INPUTR/R/1/2
COMP/KEY INPUTR(3) = INPUTR(1)-2.*M$ABS(INPUTR(2))
COMP/KEY INPUTR(4) = INPUTR(1)+10.*M$ABS(INPUTR(2))
COPY/KD INPUTR/R/3/2 'P1' LHCUTS/R/1/2
RD 'P1' LHCUTS
LOAD 'P1' chanl='p2' SCALE='p3' cent=c,c
disp/chan 'p2'
!
