!
! analyses series of frames 
!
! JWF sep2011
!

! do before:
! input files as indis/fits *rdoff*.fits - this stores files as toto0001...
! or like indis/fits *xyz*fits root=nov23 ...

defi/par p1 ds ima "root= ?"
defi/par p2 1,10 n "1., last file in sequence=?" 
defi/par p3 10,2,10 n "files pack size, start,end in each pack=?"
defi/par p4 sp c "window in frame [Q1,Q2,...]=?" 
defi/par p5 dsub06 ima "output root=?"

defi/loc root/c/1/20 {p1}
defi/loc nf/i/1/2 {p2}
defi/loc np/i/1/3 {p3}
defi/loc out/c/1/20 {p5}
!defi/loc nfil/i/1/2 {p6}


defi/loc wind/c/1/20 {p4}
defi/loc area/c/1/30 " " all

if "{wind(1:2)}" .eq. "Q1" then
   area = "[@10,@10:@2030,@2030]"
endif


if "{wind(1:2)}" .eq. "Q2" then
   area = "[@2100,@10:@4080,@2100]"
endif


if "{wind(1:2)}" .eq. "Q3" then
   area = "[@2100,@2100:@4080,@4080]"
endif


if "{wind(1:2)}" .eq. "Q4" then
   area = "[@10,@2100:@2100,@4080]"
endif


if "{wind(1:2)}" .eq. "al" then
   area = "[@10,@10:@4080,@4080]"
endif

if "{wind(1:2)}" .eq. "sp" then
   area = "[@970,@3086:@1222,@3127]"
  area = "[@768,@3350:@1268,@3850]"
  area = "[@1900,@2190:@1950,@2240]"
  area = "[@1940,@2185:@1980,@2250]"
endif



defi/loc n/i/1/1 0
defi/loc m/i/1/1 0
defi/loc m1/i/1/1 0
defi/loc m2/i/1/1 0
 

crea/icat ana.cat {p1}*.bdf
sho/icat ana.cat >NULL
defi/loc nftot/i/1/1 {outputi(1)}

analyse:

write/out use area {area} for analysis ...

crea/tab {out} 6 5
crea/col {out} :itime "sec" 
crea/col {out} :signal "ADU"
crea/col {out} :stdev "ADU"




defi/loc nrow/i/1/1 0
do n = nf(1) nf(2)
	nrow = nrow+1
      stat/ima {root}{n} {area} plot=p >NULL
      {out},:signal,@{nrow} = m$value(outputr(3))
      {out},:stdev,@{nrow} = m$value(outputr(4))
      cop/dk {root}{n} itime itime
      {out},:itime,@{nrow} = m$value(itime) 

enddo

rt {out}


set/gra
assi/gra g,1
plot/tab {out} :itime :signal 
regr/lin {out} :signal :itime
save/regr {out} fit
comp/regr {out} :fitsig = fit

set/gra color=2 stype=7
over/tab {out} :itime :fitsig
set/gra stype=0 ltype=1
over/tab {out} :itime :fitsig
set/gra

set/for f10.4
label/grap slope={outputd(2)} 50,100,mm ? ? 1           

   
exit:
set/gra
set/for

rea/tab {out}



