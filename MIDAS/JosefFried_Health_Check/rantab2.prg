DEFI/PAR P1 z C " table_out=?"
DEFI/PAR P2 4,1000 N "nser,nobj=?"
defi/par P3 0.5,1 n "effective radius, excentricity =?"
defi/par P4 100000 n "seed=?"

WRITE/KEYW IN_A/C/1/8 {P1}
WRITE/KEYW INPUTI/I/1/6 {P2}
write/keyw inputr/r/1/2 {P3}

if {P4} .eq. -1 then
   com/key inputi(3) = m$secs()
   write/out seed for random number generator = {inputi(3)}
endif


RUN FMP:rantab2

!pt 'p1' #1 #2


comp/tab 'p1' :r = sqrt(:x**2+:y**2)
!plot/hist 'p1'.tbl :r


defi/local alfa/r/1/1 12.0
defi/local delta/r/1/1 30.
defi/local size/r/1/1 0.5

comp/tab 'p1' :al = :x*{size}+{alfa}
comp/tab 'p1' :cd = {delta}
comp/tab 'p1' :al = :al/cos(:cd)

comp/tab 'p1' :del = :y*{size}+{delta}
name/col 'p1' :al R12.6
name/col 'p1' :del S12.6
plot/tab 'p1' :al :del
