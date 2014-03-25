DEFI/PAR P1 ? C "input file=?"
DEFI/PAR P2 ? C "output file=?"
DEFI/PAR P3 L C "mode [Lin/Mag]= ?" 
WRITE/KEYW IN_A/C/1/20 'P1'
WRITE/KEYW OUT_A/C/1/20 'P2'
WRITE/KEYW INPUTC/C/1/1 'P3'
RUN FMP:extinct
WRITE/DESC 'P2' EXTINCT/C/1/100 "extinct 'P1' 'P2' 'P3'"
