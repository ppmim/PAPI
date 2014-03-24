!
! fluxcalibration
!
defi/par p1 spz ima "input image=?"
defi/par p2 z ima "output image=?"
defi/par p3 feige34 c "standard star data file (output from @@ stand) =?"
defi/par p4 stand  c "standard star measurement file=?"


defi/local extimo/r/1/7 0
defi/local airmo/r/1/1/ 1.

defi/local extims/r/1/7 0
defi/local airms/r/1/1/ 1.

defi/local fac/r/1/1 1.
 
cop/dk {p1}  O_TIME extimo
cop/dk {p1} O_AIRM airmo

cop/dk {p4}  O_TIME extims
cop/dk {p4} O_AIRM airms


comp/key fac = extims(7) / extimo(7)

comp/ima {p2} = {fac} * {p1} * {p3} / {p4}

copy/dd {p4} *,3 {p2}
