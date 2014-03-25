! Subtracts the background
define/parameter p1 ? ima file=?

!loa {p1} sc=-6 cuts=600,1200
fit/flat_sky {p1} cursor 2,2 {p1}back
comp {p1}sub = {p1}-{p1}back
loa {p1}sub sc=-6 cuts=-50,400
cco
goto exit
 




define/local fwhm/d/1/1
define/local aperture/d/1/1
define/local exptim/r/1/1

outdisk/fits {p1} {p1}.fits

write/out starting SExtractor 
$/disk-b/fried/sex/sextractor2.0.21/source/sex {p1}.fits -c rebin.sex

indisk/fits background.fits {p1}back.bdf

compute/image {p1}sub = {p1} - {p1}back
write/out created {p1}sub ...

!opy/dk {p1} O_TIME/d/7/1 exptim

!$ rm {p1}back.bdf {p1}.fits {p1}.bdf background.fits

!compute/image {p1} = {p1}sub / {exptim}
!$ rm {p1}sub.bdf



exit:
