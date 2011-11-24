!
! fit plane to data points
!

crea/tab test 3 20
crea/col test :x
crea/col test :y
crea/col test :z


! setze A=2 B=3.5 C=2


comp/tab test :x = (SEQ-10)*35
comp/tab test :y = (SEQ-10)*35
comp/tab test :z = 0.*:x+3.5*:y+2.5

crea/ran noise 1,20
cop/it noise noise
cop/tt noise #1 test :noise

comp/tab test :z = :z+10*:noise

comp/tab test :x = :x+:noise
comp/tab test :y = :y-2*:noise

regr/lin test :z :x,:y
save/regr test fit
comp/regr test :fit = fit
rt test


goto exit
exit:
