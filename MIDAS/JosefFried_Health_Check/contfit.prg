DEFI/PAR P1 ? C "input file=?"
DEFI/PAR P2 ? C "output file=?"
DEFI/PAR P3 3,3 n "degree =[3,3] ?" 
defi/par p4 5. n "kappa = ?"
WRITE/KEYW IN_A/C/1/20 'P1'
WRITE/KEYW OUT_A/C/1/20 'P2'
WRITE/KEYW INPUTI/I/1/2 'P3'
write/key inputr/r/1/1 {p4}
plot {p1}
RUN FMP:contfit
set/grap color=2
over {p2}
set/grap
