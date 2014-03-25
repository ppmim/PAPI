!
! prepare flats for WFI run july 1999
!
!

write/key mfile/i/1/1 0
set/format i4



goto a


! create mask

! 1. night eve flat - read files on DAT1

do mfile = 31 31
	inta/fits {mfile} inp /dev/rmt/0bn SOA
enddo

! create mask

@@ mosaic inp0031 inp0031_all nomask
replace/ima inp0031_all msk 0,300=0.
comp/ima msk = msk/msk







! flat1       files 31-40  dat1 


@@ indata 31 40 dat1 yy


do mfile = 31 40
        @@ norm dat1{mfile} nor1{mfile}
enddo

crea/icat nor1 nor1*.bdf
aver/ima flat1 = nor1.cat ? ? median 0.9,3.








! 1. night mor flat   files 88+89  on DAT2

!@@ indata 89 90 dat2 yy

do mfile = 89 90
        @@ norm dat2{mfile} nor2{mfile}
enddo

crea/icat nor2 nor2*.bdf
aver/ima flat2 = nor2.cat ? ? median 0.9,3.











! 2. night eve flat - read files from DAT3


@@ indata 15 22 dat3 yy



do mfile = 17 22
        @@ norm dat3{mfile} nor3a{mfile}
enddo

crea/icat nor3a nor3a*.bdf
aver/ima flat3a = nor3a.cat ? ? median 0.9,3.







! 2. night mor flat - files from DAT3

a:
@@ indata 113 119 dat3 yy 

do mfile = 113 119
        @@ norm dat3{mfile} nor3b{mfile}
enddo

crea/icat nor3b nor3b*.bdf
aver/ima flat3b = nor3b.cat ? ? median 0.9,3.





exit:
