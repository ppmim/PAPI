!
!  cosmic.prg
!
defi/par p1 z char catalog
defi/par p2 2,5 n frames
defi/par p3 med ima med.frame
defi/par p4 sig ima sig.frame
defi/par p5 3.5 n kappa
! mode=1 => comp med,sig
! mode=2 => stop after med,sig computation
defi/par p6 0 n mode
write/keyw in_a/c/1/20 {p1}
write/keyw inputi/i/2/2 {p2}  
write/keyw out_a/c/1/20 {p3}  
write/keyw out_b/c/1/20 {p4}  
write/keyw inputr/r/1/1 {p5}
write/keyw inputi/i/1/1 {p6}
! replace 1. character in file for result files by
write/keyw inputc/c/1/1 h
run FMP:cosmic 
