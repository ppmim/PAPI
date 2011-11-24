DEFI/PAR P1 red C "red file=?"
DEFI/PAR P2 green  C "green file=?"
DEFI/PAR P3 blue C "red file=?"
DEFI/PAR P4 test C "out file=?"

!DEFI/PAR P5 0,0 n "cuts r = ?"
!DEFI/PAR P6 0,0 n "cuts g = ?"
!DEFI/PAR P7 0,0 n "cuts b = ?"

WRITE/KEYW IN_A/C/1/20 {p1}
WRITE/KEYW IN_B/C/1/20 {p2}
WRITE/KEYW IN_C/C/1/20 {p3}
WRITE/KEYW OUT_A/C/1/20 {p4}


defi/par p5 8,0.05,1000000 n "Q,alpha,m=?"

!WRITE/KEYW inputr/r/1/2 {p5}
!WRITE/KEYW inputr/r/3/2 {p6}
!WRITE/KEYW inputr/r/5/2 {p7}

WRITE/KEYW inputr/r/7/3 {p5}


$rm {p4}.pmm
!RUN /disk1/fried/midas/prog/RGB.exe

run RGB.exe

write/out the file {p4}.pmm is very large --- convert  it 
! $xv {p4}.pmm