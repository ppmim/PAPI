!
! apmcorfu.rpg
!
defi/par p1 acorr inputcatalog=?
defi/par p2 -1,-1 n "# of selected tables"
defi/par p3 16,3,4,5 n #rows,icol(1),icol(2),icol(3)
defi/par p4 out char "out table"
write/keyw in_a/c/1/20 {p1}
write/keyw inputi/i/1/2 {p2}  
write/keyw inputi/i/3/4 {p3}  
write/keyw out_a/c/1/20 {p4}
run FMP:apmcorfu
