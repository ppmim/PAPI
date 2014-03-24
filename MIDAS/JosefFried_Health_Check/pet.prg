
goto a





!
wd 1831 start/d/1/3 1,1
wd 1831 step/d/1/3 1,1

@@ sex 1831 1831all
sel/tab 1831all :mag_auto.lt.20
cop/tab 1831all 1831

loa 1831
loa/tab 1831 :xcen :ycen ? ? 10


a:



crea/tab usno 12 1000 usno.cat USNO2Midas.fmt



del/col usno :PPM
comp/tab usno :PPM = 100000+SEQ
comp/tab usno :MAG = :B_mag
name/col usno :RAH :R_A F15.6
name/col usno :GRAD :DEC f15.6
comp/tab usno :PMA = 0.
comp/tab usno :PMD = 0.
comp/tab usno :ident = 1000+SEQ
comp/tab usno :std = 1
!comp/tab usno :xerr = 0.
!comp/tab usno :yerr = 0.

sel/tab usno :R_A.gt.240.5496.and.:R_A.lt.240.62583
sel/tab usno selection.and.:DEC.gt.17.8299.and.:DEC.lt.17.90277
sort/tab usno :DEC
sel/tab usno all

goto exit

sort/tab 1831 :ycen


set/cont astromet
astro/transfo usno 1831 mean 2000,2000 n,n  11000000,11000000

!? 111000000,111000000 A



