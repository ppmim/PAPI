! iso.prg
!
! isocontour program
!




defi/local level/r/1/3 0,0,0
 
inquire/key level " enter levels: min,max,inc :"
write/out mark lower left, upper right corner


set/curs 1,2 rect 100,100,200,200
set/grap frame=square

plot/cont {idimemc(1:20)} C ? {level(1)}:{level(2)}:{level(3)}


set/grap
set/curs
