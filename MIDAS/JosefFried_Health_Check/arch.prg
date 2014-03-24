!
! archiving programm
!
defi/par p1 t0016 ima "frame=?"
defi/par p2 Test.sylk ima " output file [ number only]= ? "
 
defi/local file/i/1/2 0,0
open/file {p2} write file

cop/dk {p1} ident ident
write/file {file} {ident}

cop/dk {p1} O_AIRM O_AIRM
write/file {file} C;X2;Y1;K{O_AIRM}
close/file {file}

