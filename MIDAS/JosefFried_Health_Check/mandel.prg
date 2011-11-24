!
!++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  man out a,b,side        npix   nloop 
!
!  def: z   -2.,-1.25,2.5   500    100
!--------------------------------------------------
!
defi/par p1 z ima "outfile=?"
defi/par p2 -2.,-1.25,2.5 n "a,b,side=?"
defi/par p3 1000 n npix=?
defi/par p4 10000 n "nloop=?"
WRITE/KEYW OUT_A/C/1/60 'p1'
WRITE/KEYW INPUTD/D/1/3 'p2'
WRITE/KEYW INPUTI/I/1/1 'p3'
WRITE/KEYW INPUTI/I/2/2 'p4'
RUN FMP:mandel
COPY/KD INPUTI/I/2/2 'p1' LOOP/I/1/1 
