!
! login.prg 
!
! jwf may 97
!



echo/on
crea/comm ald @@ :ald
crea/comm ld @@ :ld
crea/comm ad0 assi/disp d,0
crea/comm ad1 assi/disp d,1
!crea/comm cdc copy/disp color P6=LC8B
!crea/comm cdb copy/disp color_b  P6=LC8B
crea/comm cdis copy/dis laser_b ? N
crea/comm cg center/gauss
crea/comm dp0 disp/chan 0
crea/comm dp1 disp/chan 1
crea/comm hlp @@ hlp
crea/comm rk read/key
crea/comm wk write/key
crea/comm rd read/desc
crea/comm wd write/desc
crea/comm rt read/tab
crea/comm pt plot/tab

crea/comm agl assi/graph laps_e
crea/comm agp assi/grap postscript
crea/comm ag0 assi/graph g,0
crea/comm ag1 assi/graph g,1
crea/comm ef echo/full 
crea/comm eo echo/off

crea/comm cco clea/chan over
crea/comm stc stat ? cursor

!
!
!   
!
!
wk user "JWF"
set/grap stype=6
set/grap pmode=2
set/grap ssize=1.0

set/midas_system editor=$edt
set/midas_system prompt=Midas
set/midas_system user=novice
set/midas_system environment=MidHost
set/midas_system path=+/disk1/fried/079.A-0644A.1/prg/
set/midas_system path=+/disk1/fried/midas/prog/
!set/midas_system path=+/disk1/fried/midas/prog/laica/


!
!
! create displays
!
crea/disp 0 
crea/grap 0 
load/lut heat
disp/lut on
!
eo
