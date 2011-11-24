defi/par p1 qso ima "outputframe = ?"
defi/par p2 1.,2 n "z, Av = ?"
defi/par p3 3 n " mode[1=lam*flam,2=flam,3=fnu,4=mJ] = ?"
defi/par p4 2300,1 n "npx,npy = ?"
defi/par p5 3000.,5. n "start, step = ?"
defi/par p6 qsospect.dat ima "inputfile = ?"
WRITE/KEYW IN_A/C/1/60 'P6'
WRITE/KEYW OUT_A/C/1/60 'P1'
WRITE/KEYW INPUTR/R/3/2 'P2'
WRITE/KEYW INPUTI/I/3/1 'P3'
WRITE/KEYW INPUTI/I/1/2 'P4'
WRITE/KEYW INPUTR/R/1/2 'P5'
RUN FMP:qsospect

!plot 'p1'
cop/it 'p1' 'p1' :lambda
comp/tab 'p1' :lambda = :lambda*0.1
name/col 'p1' #1 "nm"
name/col 'p1' #2 :qso "fnu"
