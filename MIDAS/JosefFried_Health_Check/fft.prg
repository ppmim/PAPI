DEFI/PAR P1 1 N "direction=?"file=?"
DEFI/PAR P2 IN C "input real=?"
DEFI/PAR P3 zero C "input imag=?"
DEFI/PAR P4 real C "output real=?"
DEFI/PAR P5 imag C "output imag=?"
WRITE/KEYW IN_A/C/1/20 'P2'
WRITE/KEYW IN_B/C/1/20 'P3'
WRITE/KEYW OUT_A/C/1/20 'P4'
WRITE/KEYW OUT_B/C/1/20 'P5'
WRITE/KEYW INPUTI/I/1/1 'P1'
RUN FMP:fft
