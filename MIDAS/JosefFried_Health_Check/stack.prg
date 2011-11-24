!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  stack.prg
!
! computes average,median od mode for stack of frames
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!


if p1 .eq. "help" then 
   
   write/out stack files with different modes
   write/out @@ stack incat frame_from,to outframe sigmaframe data_range_low,high,kappa
   write/out mode=0 => output = average
   write/out mode=1 => output = median
   write/out mode=2 => output = mode = 3*mean-2*median
   write/out  mode=3 => output = biweight
   goto exit
endif



defi/par p1 test char "input catalog = ? "
defi/par p2 1,9 n "frames from,to = ? "
defi/par p3 out ima "output frame = ? "
defi/par p4 sig ima "sigma frame = ? "
defi/par p5 -100000,100000,3 n "data range low,high,kappa = ? "

! mode=0 => output = average
! mode=1 => output = median
! mode=2 => output = mode = 3*mean-2*median
! mode=3 => output = biweight

defi/par p6 0 n mode
defi/par p7 1000000000 n "max.memory = ? "

write/keyw in_a/c/1/20 {p1}
write/keyw inputi/i/2/3 {p2}  
write/keyw out_a/c/1/20 {p3}  
write/keyw out_b/c/1/20 {p4}  
write/keyw inputr/r/1/3 {p5}
write/keyw inputi/i/1/1 {p6}
!write/keyw inputi/i/4/1 {p7}

run /home/staff/prog/stack 

write/desc {p3} stack_par/r/1/3 {p5}
write/desc {p3} stack_mode/i/1/1 {p6}
 

exit: