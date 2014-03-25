!
! cosmology.prg
!
! computes various comological quantities
!
! JWF oct 99
!




defi/par p1 cosmology  c "output table = ? "
defi/par p2 75,0.3,0.7 N "H0, Omega_m, Omega_lambda = ? "
defi/par p3 0.1,8.,0.1 n "z1,z2,delz = ? "

write/key inputd/d/1/3 {p2}
write/key inputd/d/4/3 {p3}
run /disk2/fried/midas/prog/cosmology


crea/tab {p1} 6 100 cos.dat
name/col {p1} #1 :z " " f6.2
name/col {p1} #2 :d_ang "Mpc" 
name/col {p1} #3 :d_lum "Mpc"
name/col {p1} #4 :dvdz "Mpc**3"
name/col {p1} #5 :Vol "Mpc**3"
name/col {p1} #6 :t_lookback "years" 
comp/tab {p1} :d_mod = 5*log10(:d_lum)+25.
name/col {p1} :d_mod f6.2
write/out output -> {p1}.tbl
