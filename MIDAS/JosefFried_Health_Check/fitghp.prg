!
! fit gauss-hermite polynomes
!
defi/par p1 power_fcq ima "input file=?"

defi/par p2 0.1,50,3,0.,0. n "ampl,v,sig,h3,h4=?"


cop/it {p1} zkl1
comp/tab zkl :pos = SEQ*1.
pt zkl1 :pos #1
