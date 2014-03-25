
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! normalize frames

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



defi/par p1 in ima  "input root = ? "
defi/par p2 out ima " output root = ? "
defi/par p3 1,1 n "from, to = ? "


defi/local frame/i/1/2 {p3}


defi/local lo/r/1/1 0.
defi/local hi/r/1/1 0.

set/format i4



stat/ima {p1} [@400,@400:@2000,@4000] 
hi = {outputr(8)}+5*{outputr(4)}
lo = {outputr(8)}-5*{outputr(4)}

stat/ima {p1} [@400,@400:@2000,@4000] #1000 {lo},{hi}
hi = {outputr(8)}+3*{outputr(4)}
lo = {outputr(8)}-3*{outputr(4)}

stat/ima {p1} [@400,@400:@2000,@4000] #1000 {lo},{hi}


comp/ima {p2} = {p1}/{outputr(8)}
   

write/out frame {p1} normalized by factor {outputr(8)}




