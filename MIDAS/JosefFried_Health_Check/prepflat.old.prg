!
! prepare flats for WFI run july 1999
!
!

defi/local mfile/i/1/1 0
set/format i4

goto a

! 1. night eve flat - read files on DAT1

do mfile = 31 40
	inta/fits {mfile} inp /dev/rmt/0bn SOA
enddo

! create mask

@@ mosaic inp0031 nomask
replace/ima inp0031_all msk 0,300=0.
comp/ima msk = msk/msk



! remove bias, mosaic images + normalize

a:

do mfile = 31 40
        @@ norm rbinp{mfile}_mos mfile}
enddo

crea/icat nor nor*.bdf
aver/ima flat1 = nor.cat ? ? median 0.9,1.6




! 1. night mor flat - read files on DAT2

do mfile = 89 90
	inta/fits {mfile} inp /dev/rmt/0 
enddo


do mfile = 89 90
	@@ rmbias inp{mfile}
	@@ mosaic rbinp{mfile} msk
        @@ norm rbinp{mfile}_all nor2{mfile}
enddo

crea/icat nor2 nor2*.bdf
aver/ima flat2 = nor2.cat ? ? median 0.9,1.6








! 2. night eve flat - read files from DAT3

a:

do mfile = 10 17
	inta/fits {mfile} inp /dev/rmt/0
	stat/ima inp{mfile} 
enddo

goto exit






do mfile = 10 17
	@@ rmbias inp{mfile}
	@@ mosaic rbinp{mfile} msk
        @@ norm rbinp{mfile}_all nor3{mfile}

enddo

crea/icat nor3 nor3*.bdf
aver/ima flat3 = nor3.cat ? ? median 0.9,1.6






! 2. night mor flat - files from DAT3



do mfile = 102 112
	@@ rmbias inp{mfile}
	@@ mosaic rbinp{mfile} msk
        @@ norm rbinp{mfile}_all nor4{mfile}
!	$rm *inp{mfile}.bdf
!	$rm *inp{mfile}*.bdf
enddo

crea/icat nor4 nor4*.bdf
! del/icat nor4 ? 1,20  delete all except 111,112

aver/ima flat4 = nor4.cat ? ? median 0.9,1.6



goto exit








! 3. night eve flat - read files from DAT4

do mfile = 8 11
	inta/fits {mfile} inp /dev/rmt/0
enddo
do mfile = 8 11
	@@ rmbias inp{mfile}
	@@ mosaic rbinp{mfile} msk
        @@ norm rbinp{mfile}_all nor5{mfile}
	$rm *inp{mfile}.bdf
	$rm *inp{mfile}*.bdf
enddo

crea/icat nor5 nor5*.bdf
aver/ima flat5 = nor5.cat ? ? median 0.9,1.6




exit:
