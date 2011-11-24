DEFI/PAR P1 ? C "input file=?"
DEFI/PAR P2 ? C "weight file=?"
DEFI/PAR P3 ? C "output file=?"
DEFI/PAR P4 3 N "order of fit=?"
DEFI/PAR P5 10,10,1024,1024 N "lower left, upper right=?"
WRITE/KEYW IN_A/C/1/8 'P1'
WRITE/KEYW IN_B/C/1/8 'P2'
WRITE/KEYW OUT_A/C/1/8 'P3'
WRITE/KEYW INPUTI/I/1/1 'P4'
WRITE/KEYW INPUTI/I/2/5 'P5'
RUN FMP:backfit
