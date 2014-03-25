defi/par p1 planck ima "output table = ? "

defi/par p2 5000 n "T = ? "

defi/local T/r/1/1 {p2}
defi/local k/r/1/1 1.3807e-16 
defi/local c/r/1/1 3e10
defi/local h/r/1/1 6.62e-27

crea/tab {p1} 2 100
crea/col {p1} :Wellenlaenge R "cm" F10.4
crea/col {p1} :Intensitaet R F10.4


comp/tab {p1} :Log_Wellenlaenge = -6+SEQ*0.1
comp/tab {p1} :lam = 10**:Log_Wellenlaenge

comp/tab {p1} :Intensitaet = (2*{h}*{c}**2)/(:lam**5*(exp({h}*{c}/(:lam*{k}*{T}))-1)) 

set/grap stype=0
plot/tab {p1} :Log_Wellenlaenge log(#2)
set/grap