! Program for correcting the cross-talking of NICS,
! based on the original work of  Leonardo Testi
! translated into a fortran program by E. Oliva 
!
! Version 1.3 March 2006 by F. Ghinassi
! (after Fasti installation)
! Version for SUN, HP and similar platforms (without byte-swapping)
!
! The source file (crt_nics6.f) is compiled by
!
!   f77 -o crt_nics crt_nics6.f
! or
!   g77 -o crt_nics crt_nics6.f
!
! To run it use
!
!   ./crt_nics in_file [out_file]
!
! where in_file and out_file  are the names of the input and output
! files, respectively. If out_file is omitted it will take the name
! of the input file with a "c" in front of it.
!
! Example:
!
!   ./crt_nics EDXZ0039.fts
!
! will produce a corrected image called cEDXZ0039.fts
! 
      integer*2 ima2(1024,1024)
      character*30 filein
      character*40 fileout

      common/images/araw(512,512),braw(512,512),
     ,              craw(512,512),draw(512,512),
     ,              acor(512,512),bcor(512,512),
     ,              ccor(512,512),dcor(512,512)
      common/cross_talk/rfAqB,rfAqC,rfAqD,dxAqB,dxAqC,dxAqD,
     ,                  rfBqA,rfBqC,rfBqD,dxBqA,dxBqC,dxBqD,
     ,                  rfCqA,rfCqB,rfCqD,dxCqA,dxCqB,dxCqD,
     ,                  rfDqA,rfDqB,rfDqC,dxDqA,dxDqB,dxDqC,
     ,                  factor

!
! valori determinati da Leo Testi il 12/11/2000
!
      factor=0.8e-5  ! fattore di riga, effetto della somma di cio' che c'e' 
                     ! sulla riga a Y=coste su ogni pixel dello stesso Y degli
                     ! altri 3 quadranti. NB questo e' positivo, cioe' ha 
                     ! effetto opposto a quello del cross-talking
!
! I quadranti sono:
!
!                  C  D
!                  A  B
!

      rfAqB=0.040    ! effetto di A su B
      rfAqC=0.040
      rfAqD=0.040
      dxAqB=0.40    ! parte di A(x+1)*rfAqB che va su B(x)
      dxAqC=0.40
      dxAqD=0.40
   
      rfBqA=0.041
      rfBqC=0.039
      rfBqD=0.039
      dxBqA=0.35
      dxBqC=0.35
      dxBqD=0.35
   
      rfCqA=0.039
      rfCqB=0.039
      rfCqD=0.039
      dxCqA=0.42
      dxCqB=0.42
      dxCqD=0.42

      rfDqA=0.038
      rfDqB=0.038
      rfDqC=0.037
      dxDqA=0.38
      dxDqB=0.38
      dxDqC=0.38

!
! Per i nomi files cf. 
!       crt_nis
!
      call getarg (1,filein)
      if(filein.eq." ") then
	write(6,'(a)') " "
	write(6,'(a)') 
     +     "Purpose:  correct a NICS frame for cross-talking"
	write(6,'(a)') " "
	write(6,'(a)') "Usage:  crt_nics  in_file [out_file]"
	write(6,'(a)') 
     +   'if out_file is omitted it will be "c"in_file'
	write(6,'(a)') " "
	write(6,'(a)') "Example: crt_nics EXAZ0001.fts"
	write(6,'(a)') 
     +  "           writes corrected frame in cEXAZ0001.fts"
	write(6,'(a)') " "
	stop
      endif
      call getarg (2,fileout)
      open(unit=1,file=filein, 
     ,     form='unformatted',recl=2880,access='direct')
      if(fileout.eq." ") write(fileout,'(a,a)') "c",filein
      open(unit=2,file=fileout,
     ,     form='unformatted',recl=2880,access='direct')
      call get_ima(filein,fileout,nrec_header,ima2)
      close(unit=1)
!
! Initialize a,b,c,d subimages
!
      do j=1,512
	do i=1,512
	  araw(i,j)=ima2(i,j)
	  braw(i,j)=ima2(i+512,j)
	  craw(i,j)=ima2(i,j+512)
	  draw(i,j)=ima2(i+512,j+512)
	  acor(i,j)=araw(i,j)
	  bcor(i,j)=braw(i,j)
	  ccor(i,j)=craw(i,j)
	  dcor(i,j)=draw(i,j)
	enddo
      enddo
!
! Apply cross-talking corrections
!
      call apply_cross_talking
!
! Copy results in outside image and write it in output file
!
      do j=1,512
        do i=1,512
	  ima2(i,j)=acor(i,j)
	  ima2(i+512,j)=bcor(i,j)
	  ima2(i,j+512)=ccor(i,j)
	  ima2(i+512,j+512)=dcor(i,j)
        enddo
      enddo
      call out_ima(ima2,nrec_header)
      close(unit=2)
      stop ''
      end

      subroutine apply_cross_talking

      common/images/araw(512,512),braw(512,512),
     ,              craw(512,512),draw(512,512),
     ,              acor(512,512),bcor(512,512),
     ,              ccor(512,512),dcor(512,512)
      common/cross_talk/rfAqB,rfAqC,rfAqD,dxAqB,dxAqC,dxAqD,
     ,                  rfBqA,rfBqC,rfBqD,dxBqA,dxBqC,dxBqD,
     ,                  rfCqA,rfCqB,rfCqD,dxCqA,dxCqB,dxCqD,
     ,                  rfDqA,rfDqB,rfDqC,dxDqA,dxDqB,dxDqC,
     ,                  factor

      do iter=1,2
        do j=1,512
	  asum=0.
	  bsum=0.
	  csum=0.
	  dsum=0.
!
! Correzione cross-talking negativo
!
  	  do i=1,512
            if(i.lt.512) then
	      ash=acor(i+1,j)
	      bsh=bcor(i+1,j)
	      csh=ccor(i+1,j)
	      dsh=dcor(i+1,j)
	    else
	      ash=0.
	      bsh=0.
	      csh=0.
	      dsh=0.
	    endif
	    acor(i,j)=araw(i,j)+
     +                rfBqA*((1-dxBqA)*bcor(i,j)+dxBqA*bsh)+
     +                rfCqA*((1-dxCqA)*ccor(i,j)+dxCqA*csh)+
     +                rfDqA*((1-dxDqA)*dcor(i,j)+dxDqA*dsh)
	    bcor(i,j)=braw(i,j)+
     +                rfAqB*((1-dxAqB)*acor(i,j)+dxAqB*ash)+
     +                rfCqB*((1-dxCqB)*ccor(i,j)+dxCqB*csh)+
     +                rfDqB*((1-dxDqB)*dcor(i,j)+dxDqB*dsh)
	    ccor(i,j)=craw(i,j)+
     +                rfAqC*((1-dxAqC)*acor(i,j)+dxAqC*ash)+
     +                rfBqC*((1-dxBqC)*bcor(i,j)+dxBqC*bsh)+
     +                rfDqC*((1-dxDqC)*dcor(i,j)+dxDqC*dsh)
	    dcor(i,j)=draw(i,j)+
     +                rfAqD*((1-dxAqD)*acor(i,j)+dxAqD*ash)+
     +                rfBqD*((1-dxBqD)*bcor(i,j)+dxBqD*bsh)+
     +                rfCqD*((1-dxCqD)*ccor(i,j)+dxCqD*csh)

	    asum=asum+acor(i,j)
	    bsum=bsum+bcor(i,j)
	    csum=csum+ccor(i,j)
	    dsum=dsum+dcor(i,j)
	  enddo
!
! Correzione effetto di riga (positivo)
!
	  do i=1,512
	    acor(i,j)=acor(i,j)-factor*(bsum+csum+dsum)
	    bcor(i,j)=bcor(i,j)-factor*(asum+csum+dsum)
	    ccor(i,j)=ccor(i,j)-factor*(asum+bsum+dsum)
	    dcor(i,j)=dcor(i,j)-factor*(asum+bsum+csum)
	  enddo
	enddo
      enddo

      return
      end

      subroutine get_ima(filein,fileout,nrec_header,ima_in)
!
! 1) read header records (up to END) and copies them into output
!    file adding an history string with the name of the files
!    Store # of last written record in "nrec_header"
! 2) Read the 1024^2 I*2 image matrix and store it in ima_in
!
      character*30 filein,fileout
      character*2880 crecord
      character*80 string
      integer*2 ima_in(1)

      do i=1,100
        irec=i
	read(1,rec=irec,err=99) crecord
        do j=1,36
	  string(1:80)=crecord(((j-1)*80+1):((j-1)*80+80))
          if(string(1:3).EQ.'END') then
	    crecord(((j-1)*80+1):((j-1)*80+80))=' '
	    write(crecord(((j-1)*80+1):((j-1)*80+80)),'(4a)')
     ,      'HISTORY crt_nics ',filein,' ',fileout
            if(j.le.35) then
	      write(crecord(((j-0)*80+1):((j-0)*80+80)),'(a)')
     ,        'END'
              write(2,rec=irec) crecord
	      nrec_header=irec
	      goto 10
            else
	      write(2,rec=irec) crecord
	      nrec_header=irec+1
	      crecord='END'
	      write(2,rec=nrec_header) crecord
	      goto 10
            endif
	  endif
	enddo
	write(2,rec=irec) crecord
      enddo
99    stop 'END not found in header'

10    continue
      n=1
      do i=irec+1,irec+729
        n2=n+1439
        n2=min(n2,1024*1024) ! do not go beyond 1024**2
        read(1,rec=i,err=98) (ima_in(k),k=n,n2)
        n=n2+1
      enddo
      return
98    stop 'Premature end-of-file in input fits frame'
      end

      subroutine out_ima(ima_out,nrec_header)
!
! Write I*2 matrix starting from record=nrec_header+1
!
      integer*2 ima_out(1)

      n=1
      do i=nrec_header+1,nrec_header+729
        irec=i
        write(2,rec=i) (ima_out(k),k=n,n+1439)
        n=n+1440
      enddo
      return
      end
