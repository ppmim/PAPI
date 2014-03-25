defi/par p1 test2edi ima

defi/loc tab/c/1/20 {p1}

comp/tab {tab} :var = :stdev**2

set/gra
assi/gra g,0
plot/tab {tab} :signal :var
regr/lin {tab} :var :signal
save/regr {tab} fit
comp/regr {tab} :fit = fit

set/gra color=2 
over/tab {tab} :signal :fit
set/gra stype=0 ltype=1
over/tab {tab} :signal :fit
set/gra
set/for f6.2
label/grap slope={outputd(2)} 50,100,mm ? ? 1           

   

set/gra
assi/gra g,1
plot/tab {tab} :itime :signal 
regr/lin {tab} :signal :itime
save/regr {tab} fit
comp/regr {tab} :fit = fit

set/gra color=2 
over/tab {tab} :itime :fit
set/gra stype=0 ltype=1
over/tab {tab} :itime :fit
set/gra
set/for f10.4
label/grap slope={outputd(2)} 50,100,mm ? ? 1           

   
exit:
set/gra
set/for
