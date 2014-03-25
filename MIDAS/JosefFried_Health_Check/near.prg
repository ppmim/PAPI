defi/par p1 ot566a C " table_in=?"
defi/par p2 z c "table out=?"
defi/par p3 2,3 N "x_col, y_col=?"
defi/par p4 75.,1,0.5,1.,150. n "H0, Omega, z_lens,z_source, vel_disp = "
defi/par p5 .40 N "scale=? [arcsec/pix]"
write/keyw in_a/c/1/12 {p1}
write/keyw out_a/c/1/12 {p2}
write/keyw inputi/i/1/8 {p3}
write/keyw inputr/r/1/5 {p4}
write/keyw inputr/r/6/1 {p5}


RUN FMP:near
