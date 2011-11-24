! cosmology prg

defi/par p1 0.5 n "z  " 
defi/par p2 0.72,0.27,0.73 n " h0/100,OMmat,OMlam=? [0.72,0.27,0.73] " 
defi/par p3 on ima "log on/off [on] =?"

defi/local z/r/1/1 {p1}
defi/local h/r/1/3 {p2}



 
 
defi/local h0/r/1/1 {h(1)}   !hubble constant
defi/local om/r/1/1 {h(2)}    !Omega matter
defi/local ol/r/1/1 {h(3)}   !Omega lambda



defi/local w/r/1/1 0.
defi/local e/r/1/1 0.
defi/local i/i/1/1 0
defi/local imax/i/1/1 0
defi/local zs/r/1/1 0.

defi/local dc/r/1/1 0.
defi/local dh/r/1/1 0.
defi/local da/r/1/1 0.
defi/local dl/r/1/1 0.
defi/local dm/r/1/1 0.
defi/local dvc/r/1/1 0.


set/for

! comoving distance dc

imax = z(1)*1000
do i = 1 {imax}
	zs = i*0.001
	e = m$sqrt(om*(1+zs)**3+ol)
	dc = dc + 0.001/e
enddo

dh = 3000/h0
dc = dc*dh     !now in Mpc

e = m$sqrt(om*(1+z)**3+ol)

da = dc/(1+z)
dl = (1+z)**2*da
dm = 5*m$log10(dl*1e6/10)
dvc = dh*(1+z)**2*da**2/e


if "{p3}" .eq. "on" then
write/out
write/out H0 = {h0}  Omatter = {om}  Olambda = {ol}
write/out E = {e}
write/out D_H = {dh} [Mpc]  -> outputr(1) 
write/out D_C = {dc} [Mpc]  -> outputr(2) 
write/out D_A = {da}  [Mpc] -> outputr(3) 
write/out D_L = {dl}  [Mpc] -> outputr(4) 
write/out D_M = {dm}  [mag] -> outputr(5) 
write/out dVc = {dvc}  [Mpc**3/dOmega[sr]/dz] -> outputr(6) 
endif

outputr(1) = {dh}
outputr(2) = {dc}
outputr(3) = {da}
outputr(4) = {dl}
outputr(5) = {dm}
outputr(6) = {dvc}

exit:
