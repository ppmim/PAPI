!
! analyses series of frames 
!
! JWF sep2011
!

! do before:
! input files as indis/fits *rdoff*.fits - this stores files as toto0001...
! or like indis/fits *xyz*fits root=nov23 ...

defi/par p1 toto ima "root= ?"
defi/par p2 1,10 n "1., last file in sequence=?" 
defi/par p3 10,2,10 n "files pack size, start,end in each pack=?"
defi/par p4 sp c "window in frame [Q1,Q2,...]=?" 
defi/par p5 dsub06 ima "output root=?"
defi/par p6 -1,-1 c "start, end in sequence  -1=do stacking" 

defi/loc root/c/1/20 {p1}
defi/loc nf/i/1/2 {p2}
defi/loc np/i/1/3 {p3}
defi/loc out/c/1/20 {p5}
defi/loc nfil/i/1/2 {p6}


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
  area = "[@1756,@2317:@1944,@2351]"
endif



defi/loc n/i/1/1 0
defi/loc m/i/1/1 0
defi/loc m1/i/1/1 0
defi/loc m2/i/1/1 0
 

crea/icat ana.cat {p1}*.bdf
sho/icat ana.cat >NULL
defi/loc nftot/i/1/1 {outputi(1)}

defi/loc npack/i/1/1 
npack = (nf(2)-nf(1)+1)/np(1)

write/out 
write/out files from {nf(1)} to {nf(2)} out of sequence with {nftot} files 
write/out {npack} packs with {np(1)} files each, use files {np(2)}...{np(3)} in each pack
write/out


   
if {nfil(1)} .ne. -1 then
   goto analyse
else
    nfil(1) = 1	
    nfil(2) = npack
endif



do n = 1 npack
      m1 = nf(1)+(n-1)*np(1)+np(2)-1
      m2 = m1+np(3)-np(2)
      write/out stack pack # {n} files from {m1}...{m2}
      @@ stack ana {m1},{m2},1 {out}_m{n} {out}_s{n}
      cop/dd {root}{m1} itime {out}_m{n} 
      cop/dd {root}{m1} itime {out}_s{n}
      assi/dis 0
      loa/ima {out}_m{n} sc=-2 cuts=F
      assi/dis 1
  !    loa/ima {out}_s{n} sc=-2 cuts=F
enddo



analyse:

write/out use area {area} for analysis ...

crea/tab {out} 6 5
crea/col {out} :itime "sec" 
crea/col {out} :signal "ADU"
crea/col {out} :stdev "ADU"
crea/col {out} :var "ADU"





do n = nfil(1) nfil(2)
      stat/ima {out}_m{n} {area} plot=p >NULL
      {out},:signal,@{n} = m$value(outputr(3))
      stat/ima {out}_s{n} {area} plot=p >NULL
      {out},:stdev,@{n} = m$value(outputr(3))
      cop/dk {out}_m{n} itime itime
      {out},:itime,@{n} = m$value(itime) 

enddo

rt {out}

comp/tab {out} :var = :stdev**2

set/gra
assi/gra g,0
plot/tab {out} :signal :var
regr/lin {out} :var :signal
save/regr {out} fit
comp/regr {out} :fitvar = fit

set/gra color=2 stype=7
over/tab {out} :signal :fitvar 
set/gra stype=0 ltype=1
over/tab {out} :signal :fitvar
set/gra


set/for f6.2
label/grap "slope = {outputd(2)}" 40,100,mm ? ? 1           
label/grap "intercept = {outputd(1)}" 40,95,mm ? ? 1
defi/loc w/r/1/1 0.
w = 1./outputd(2)
label/grap "gain = {w}  [e-/ADU]" 80,30,mm ? ? 1
w = m$sqrt(outputd(1)/(outputd(2)**2))  
label/grap "read noise = {w} [e-]" 80,25,mm ? ? 1

! gain = 1/slope
! intercept = (R/gain)**2 => R = gain*sqrt(intercept)


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



