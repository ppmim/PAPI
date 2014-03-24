! 
! sex.prg
!

defi/par p1 sumaBR2_F090 ima "frame name=?"
defi/par p2 sumaBR2_sex ima "out table name=?"

defi/par p3 32. n "zero_magnitude=?"
defi/par p4 19. n "lim_mag_low=?"
defi/par p5 22.5 n "lim_mag_hi=?"

outdisk/fits {p1}.bdf z.fits

$/disk-b/fried/sex/sextractor2.0.21/source/sex z.fits




crea/tab {p2} 13 9999  test.cat

write/out table created - renaming columns...

name/col {p2} #1 :ident I i4
name/col {p2} #2 :fluxerr_iso R f10.2
name/col {p2} #3 :flux_auto R f10.2
name/col {p2} #4 :fluxerr_auto R f10.2
name/col {p2} #5 :xcen R f10.3
name/col {p2} #6 :ycen R f10.3
name/col {p2} #7 :cxx R f10.4
name/col {p2} #8 :cyy R f10.4
name/col {p2} #9 :cxy R f10.4
name/col {p2} #10 :elong R f10.4
name/col {p2} #11 :ellipt R f10.4
name/col {p2} #12 :flags I i4
name/col {p2} #13 :classification F f6.2
 

write/out computing magnitudes...
comp/tab {p2} :mag_auto = {p3}-2.5*log10(:flux_auto)
name/col {p2} :mag_auto f8.2
write/out selecting objects within {p4} and {p5} mags ...
sel/tab {p2} :mag_auto.gt.{p4}.and.:mag_auto.lt.{p5}

write/out selecting stars ...
sel/tab {p2} selection.and.:classification.gt.0.5
loa/tab {p2} :xcen :ycen ? 1 6 6
write/out

write/out selecting galaxies ...
sel/tab {p2} :mag_auto.gt.{p4}.and.:mag_auto.lt.{p5}
sel/tab {p2} selection.and.:classification.lt.0.5
loa/tab {p2} :xcen :ycen ? 0 6 3
write/out 



exit:






