C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  PROGRAM STACK
C
C
C  J.FRIED  feb 2003
C
C  modifed Jun 2005
C
C----------------------------------------------------------
C

c vers 1.3.2004
c vers 6/2005



	CHARACTER*20	CATA,FRAME(500),SIGMA,LOCATION
	CHARACTER	CUNIT*64
	CHARACTER	IDENT*72
	CHARACTER	HISTORY*72
	BYTE            CR

	INTEGER  	STAT,NPIX(2),IPAR(10),npixcube(3),npixo(2)
	INTEGER	        MADRID(1)
	INTEGER         KNULL,KUNIT(1)
	
	integer*8       ip,ipc,iploc,ipsig

	REAL*4		PAR(6)
	REAL*8		START(2),STEP(2),dval(10)
	real*8          startcube(3),stepcube(3)




	INCLUDE		'st_def.inc'
	INCLUDE		'st_dat.inc'

	COMMON /VMR/ MADRID


	DATA  CR/13/

C
C-----  INITIALIZE MIDASENVIRONMENT
C
	CALL STSPRO('stack')


C
C-----  GET PARAMETERS
C
	CALL STKRDR('INPUTR',1,6,IAV,PAR,KUNIT,KNULL,STAT)
	rlow = par(1)
	rhigh = par(2)
	rkappa = par(3)


C
C------ GET MODE, max.memory
C
	CALL STKRDI('INPUTI',1,6,IAV,IPAR,KUNIT,KNULL,STAT)

	MODE=IPAR(1)
	NFR1=IPAR(2)
	NFR2=IPAR(3)
	NFR3=IPAR(4)




C
C----- OPEN CATALOG, LOOP OVER OBJECTS
C
	CALL STKRDC('IN_A',1,1,20,IAV,CATA,KUNIT,KNULL,STAT)
	CALL STCSHO(CATA,NFRAMES,LAST,STAT)

	PRINT *,' >>> program stack vers. 9/2011<<< '

	IF(MODE.EQ.0) THEN
		PRINT *,' compute average + sigma  '
	END IF
	IF(MODE.EQ.1) THEN
		PRINT *,' compute median + sigma '
	END IF
	IF(MODE.EQ.2) THEN
		PRINT *,' compute mode + sigma '
	END IF
	IF(MODE.EQ.3) THEN
		PRINT *,' compute biweight + sigma '
	END IF
	IF(MODE.EQ.4) THEN
		PRINT *,' compute aver  + ',rkappa,' * sigma clipping'
	END IF



	PRINT *,' catalogue ',CATA,' contains',NFRAMES,'  entries '
	PRINT *,' selected are frames',nfr1,' to',nfr2,' step = ',nfr3
	
	NFRAMES=(NFR2-NFR1)/NFR3+1
	
	PRINT *,' to process ' , nframes,' frames'

	PRINT *,' data range from ',rlow,' to ',rhigh


	MFR=0 
	
	   
	DO NFR=NFR1,NFR2,NFR3
		MFR=MFR+1
		CALL STCFND(CATA,NFR,FRAME(NFR),STAT)


                CALL STIGET(FRAME(NFR),D_R4_FORMAT,F_I_MODE,F_IMA_TYPE
     1             ,2,NAXIS,NPIX,START,STEP,IDENT,CUNIT,IP,IM,STAT)




		IF (MFR.EQ.1) THEN
		   do i = 1,2
		      startcube(i)=1
		      stepcube(i) = 1
		      npixcube(i)=npix(i)
                   end do

		   startcube(3) =1
      	           stepcube(3)=1
		   npixcube(3)=nframes


	
C
C----- OPEN OUTPUT DATA CUBE
C



		   CALL STIPUT('CUBE',D_R4_FORMAT,F_O_MODE,F_IMA_TYPE
     1		   ,3,NPIXCUBE,STARTCUBE,STEPCUBE,'cube',CUNIT,IPC,IMC,STAT)


C
C-----  MAP SIGMA FILE
C
		   
c		   print*, ' map sigma file',sigma
		   CALL STKRDC('OUT_B',1,1,20,IAV,SIGMA,KUNIT,KNULL,STAT)
		   CALL STIPUT(SIGMA,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE
     1                 ,2,NPIX,start,step,IDENT,CUNIT,IPSIG 
     1                  ,IMSIG,STAT)


C
C-----  MAP OUTPUT FILE
C

		   
c		   print *,' map output file',location
		   CALL STKRDC('OUT_A',1,1,20,IAV,LOCATION,KUNIT,KNULL,STAT)	
		   CALL STIPUT(LOCATION,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE
	1		,2,NPIX,START,STEP,ident,CUNIT,IPLOC,IMLOC,STAT)
		   

	   end if

C
C----- FILL DATA CUBE
C
	 

		CALL FILL(MADRID(IP),NPIX(1),NPIX(2)
     1          ,madrid(ipc),npixcube(1),npixcube(2),npixcube(3),MFR)
	
c
c---- close input frame
c	      
		PRINT *, ' frame  # ',nfr,'  ',FRAME(NFR)


		call stfclo(IM,stat)	  

	END DO


	PRINT *,' cube filled - start computation '



	CALL MEAN(MADRID(IPC),MADRID(IPLOC),MADRID(IPSIG)
	1		,NPIX(1),NPIX(2),NFRAMES,MFR,rlow,rhigh,rkappa,mode)




C
C-----  END
C

C---- CLOSE CUBE

	CALL STFCLO(IMC,STAT)

 1000	print *,' output -> ',location
	print *,' sigma -> ', sigma

	CALL STSEPI
	STOP  


	END



	SUBROUTINE FILL(A,NPX,NPY,CUBE,NPXCUBE,NPYCUBE,npzcube,NFR)

	INTEGER    I,J,NPX,NPY,npxcube,npycube,npzcube,nfr


	REAL*4	A(NPX,NPY),CUBE(NPXCUBE,NPYCUBE,npzcube)

	

	DO J=1,NPY
	   DO I=1,NPX
	      CUBE(I,J,NFR) = a(i,j)
	   END DO
	END DO


	RETURN
	END




	SUBROUTINE MEAN(CUBE,LOC,SIG,NPX,NPY,NFRAMES,NFR,rlow,rhigh
     1	,rkappa,mode)

	BYTE CR

	REAL*4	CUBE(NPX,NPY,NFRAMES),LOC(NPX,NPY),SIG(NPX,NPY)
	REAL*4	WORK(2000),work2(2000)
	integer mode,npts
	real*4  rlow,rhigh,XMED,AVER

	DATA    CR/13/


c
c------ average 
c
	if (mode .eq. 0 ) then
	   DO J=1,NPY
	      IDONE = NINT((100.*J)/NPY)
c	      PRINT 11, IDONE, CR
 11	      FORMAT(5x,'average + sigma  ',i3,' % done ',a1,$)
	      DO I=1,NPX

		 npts = 0
		 DO K=1,NFRAMES
		    if (cube(i,j,k) .ge. rlow .and. cube(i,j,k).le. rhigh) then
		       npts = npts+1
		       WORK(npts)=CUBE(I,J,K)
		    endif
		 END DO

		 LOC(I,J) = 0.
		 if(npts .eq. 0) then
		    loc(i,j) = 0.
		    sig(i,j) = 0.
		 endif
		 if(npts.eq.1) then
		    loc(i,j) = work(1)
		    sig(i,j) = 0.
		 endif
		 if(npts.ge.2) then
		    call ms_meansigma (work,npts,aver,sigma)
		    loc(i,j) = aver
		    sig(i,j) = sigma
		 endif
	      END DO
	   END DO
	endif




c
c---- median
c


	if (mode .eq. 1 ) then
	   DO J=1,NPY
	      IDONE = NINT((100.*J)/NPY)

	      PRINT 12, IDONE, CR
 12	      FORMAT(5x,'median + sigma',i3,' % done ',a1,$)

	      DO I=1,NPX

		 npts = 0
		 DO K=1,NFRAMES
		    if (cube(i,j,k) .ge. rlow .and. cube(i,j,k).le. rhigh) then
		       npts = npts+1
		       WORK(npts) = CUBE(I,J,K)
		    endif
		 END DO

		 LOC(i,j) = 0.
		 sig(i,j) = 0.
		 if (npts .ge. 2 ) then 
		    call ms_median1 (WORK,NPTS,XMED)
		    LOC(I,J) = XMED
		    sum1 = 0.0
		    sum2 = 0.0
		    do  iw = 1,npts 
		       s = xmed - work (iw)
		       sum1 = sum1 + (s * s)
		       sum2 = sum2 + s
		    end do
		    vari = (sum1 - (sum2 / real (npts))) / real (npts-1)
		    sig(i,j) = sqrt (vari)
		 endif
		 if (npts .eq. 1) then
		    LOC(i,j) = work(1)
		    sig(i,j) = 0.
		 endif

	      END DO
	   END DO
	endif



c
c---- mode = 3*mean-2*median
c

	if (mode .eq. 2) then
	   DO J=1,NPY
	      IDONE = NINT((100.*J)/NPY)
	      PRINT 13, IDONE, CR
 13	      FORMAT(5x,'mode = 3*mean-2*median ',i3,' % done ',a1,$)

	      DO I=1,NPX

		 npts = 0
		 DO K=1,NFRAMES
		    if (cube(i,j,k) .ge. rlow .and. cube(i,j,k).le. rhigh) then
		       npts = npts+1
		       WORK(npts) = CUBE(I,J,K)
		    endif
		 END DO

		 LOC(i,j) = 0.
		 sig(i,j) = 0.

		 if (npts .eq. 1) then
		    LOC(i,j) = work(1)
		    sig(i,j) = 0.
		 endif
		 if (npts .ge. 2 ) then 
		    call ms_median1 (WORK,NPTS,XMED)
		    call ms_mean (WORK,NPTS,xmean)
		    LOC(I,J) = 3.*xmean-2.*xmed
		    
		    sum1 = 0.0
		    sum2 = 0.0
		    do  iw = 1,npts 
		       s = loc(i,j) - work (iw)
		       sum1 = sum1 + (s * s)
		       sum2 = sum2 + s
		    end do
		    vari = (sum1 - (sum2 / real (npts))) / real (npts-1)
		    sig(i,j) = sqrt (vari)
		 endif

	      END DO
	   END DO

	endif	


c
c------ biweight
c
	if (mode .eq. 3 ) then
	   DO J=1,NPY
	      IDONE = NINT((100.*J)/NPY)
	      PRINT 14, IDONE, CR
 14	      FORMAT(5x,i3,' % done ',a1,$)
	      DO I=1,NPX

		 npts = 0
		 DO K=1,NFRAMES
		    if (cube(i,j,k) .ge. rlow .and. cube(i,j,k).le. rhigh) then
		       npts = npts+1
		       WORK(npts)=CUBE(I,J,K)
		    endif
		 END DO

		 LOC(I,J) = 0.
		 if(npts .eq. 0) then
		    loc(i,j) = 0.
		    sig(i,j) = 0.
		 endif
		 if(npts.eq.1) then
		    loc(i,j) = work(1)
		    sig(i,j) = 0.
		 endif
		 if(npts.ge.2) then
		    call biwgt (work,npts,aver,sigma)
		    loc(i,j) = aver
		    sig(i,j) = sigma
		 endif
	      END DO
	   END DO
	endif


c
c------ kappa*sigma
c
	if (mode .eq. 4 ) then

	   nclip = 0
	   DO J=1,NPY
	      IDONE = NINT((100.*J)/NPY)
	      PRINT 16, IDONE, CR
 16	      FORMAT(5x,'kappa*sigma',i3,' % done ',a1,$)
	      DO I=1,NPX

		 npts = 0
		 DO K=1,NFRAMES
		    if (cube(i,j,k) .ge. rlow .and. cube(i,j,k).le. rhigh) then
		       npts = npts+1
		       WORK(npts) = CUBE(I,J,K)
		    endif
		 END DO



		 if(npts .eq. 0) then
		    loc(i,j) = 0.	
		    sig(i,j) = 0.
		 endif
		 if(npts .eq. 1) then
		    loc(i,j) = work(1)	
		    sig(i,j) = 0.
		 endif
		 if(npts.ge.2) then
		    call ms_meansigma (work,npts,aver,sigma)
		    loc(i,j) = aver
		    sig(i,j) = sigma

c---    kappa+sig clipping
		    ik =0
		    do kl = 1,npts
		       if (abs(work(kl)-aver) .le. rkappa*sigma) then
			  ik = ik+1
			  work2(ik) = work(kl)
		       endif
		    enddo

		    if(ik.lt.2) then
		       loc(i,j) = work2(1)
		    endif
		    if(ik.ge.2) then
		       call ms_meansigma(work2,ik,aver,sigma)
		       loc(i,j) = aver
		    endif
		 nclip = nclip + (npts-ik)
		    
		 endif
	      END DO
	   END DO
	   mpx = npx*npy
	   clip = nclip*100/mpx
	   print *,' '
	   print *,nclip,' points of ',mpx ,' = ',clip,' % clipped'
	endif











	PRINT *,' '

	RETURN
	END







	SUBROUTINE REMCOS(IN,MED,SIG,OUT,NPX,NPY,MFR,PAR,FSTAT,FAC)


	INTEGER*4	NPIX(2),IPAR(4)

	REAL*4		IN(NPX,NPY),MED(NPX,NPY)
	REAL*4		OUT(NPX,NPY),SIG(NPX,NPY)
	REAL*4		PAR(6),FSTAT(3,100),FAC(100)
	REAL*4		KAPPA,NOISE




	KAPPA = PAR(1)
	NOISE = FSTAT(3,MFR)**2


C
C
C


	IREP=0
	DO J=1,NPY
		DO I=1,NPX

			sigma = sqrt(abs( sig(i,j)**2 + noise ))
			diff = in(i,j)*FAC(MFR) - med(i,j)

			if(diff.gt.kappa*sigma) then
				IREP = IREP+1
c                       out(i,j) = med(i,j)    !! muss das nicht hier herein?
			ELSE
				out(I,J) = in(I,J)
			END IF
		END DO
	END DO
	PRINT *,IREP,'  pixels replaced'

	RETURN
	END
		





	SUBROUTINE MDIAN1(X,N,XMED)
c
c   aus numrecip
c    
	DIMENSION X(N)
 	CALL SORT(N,X)
     	N2=N/2
	IF(2*N2.EQ.N)THEN
        	XMED=0.5*(X(N2)+X(N2+1))
      	ELSE
	        XMED=X(N2+1)
      	ENDIF
      	RETURN
     	 END






	SUBROUTINE SORT(N,RA)
      	DIMENSION RA(N)
    	L=N/2+1
    	IR=N
10    	CONTINUE
        IF(L.GT.1)THEN
        	L=L-1
          	RRA=RA(L)
        ELSE
	        RRA=RA(IR)
        	RA(IR)=RA(1)
          	IR=IR-1
          	IF(IR.EQ.1)THEN
            		RA(1)=RRA
            		RETURN
          	ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
        	IF(J.LT.IR)THEN
            		IF(RA(J).LT.RA(J+1))J=J+1
	          ENDIF
        	  IF(RRA.LT.RA(J))THEN
            		RA(I)=RA(J)
          		I=J
 		        J=J+J
	          ELSE
      		      J=IR+1
        	  ENDIF
      	  GO TO 20
        ENDIF
        RA(I)=RRA
      	GO TO 10
      	END
	


 
c-----------------------------------------------------------------------------
      subroutine biwgt (x, n, xbl, xbs)
c-----------------------------------------------------------------------------

c  This subroutine calculates the Biweight estimator of the location and
c  scale of a data set.  This is a robust estimator, and is preferred to
c  the median by all rational folks.  The Biweight function is given by:
c
c                                  u((1-u*u)**2)     abs(u) <= 1
c                         f(u) =
c                                  0                 abs(u) >  1
c
c  where u is defined by
c
c                         u = (X(I) - Med) / (c * Mad)  .
c
c  Med, Mad, and c are the median, the median absolute deviation from
c  the median, and the tuning constant respectively. The tuning constant
c  is a parameter which is chosen depending on the sample size and the
c  specific function being used for the scale estimate.  Here the tuning
c  constant is set to 6.0 for calculation of the location, and to 9.0 for
c  calculation of the scale, as recommended by Tukey.
c
c  The input data are provided in the single precision real array X, of
c  length N.  The location and scale are returned in the real variables
c  XBL and XBS, respectively.  On input errors, the scale XBS is returned
c  with negative values:
c
c        for N < 1     XBS = -2.0  and  XBL = 0.0
c        for N = 1     XBS = -1.0  and  XBL = X(1)
c
c  WARNING: This version of the Biweight routine has the side effect of
c  sorting the input array X, and returning the sorted data in place.


      real*8  sum1, sum2, t0, t1
      real    x(n), xbl, xbs, xmed, xmad, cmad, cmadsq, delta
      integer n, i



      if (n .lt. 1) then
         xbl =  0.0
         xbs = -2.0
         return
      else if (n .eq. 1) then
         xbl =  x(1)
         xbs = -1.0
         return
      else if (n .eq. 2) then
         xbl = (x(1) + x(2)) / 2.0
         xbs = abs(x(1) - xbl)
         return
      end if

      call medmad (x, n, xmed, xmad)

      if (xmad .eq. 0.0) then
         xbl = xmed
         xbs = xmad
         return
      end if

      xbl = xmed
      delta = xmad
      cmad = 6.0 * xmad
      cmadsq = cmad * cmad

      icnt=0
      do while ((abs(delta) .ge. 1.0e-4 * xmad).and.(icnt.lt.20))

         sum1 = 0.0
         sum2 = 0.0

         do i=1,n
            t0 = x(i) - xbl
            if (abs(t0) .lt. cmad) then
               t1 = cmadsq - t0 * t0
               t1 = t1 * t1
               sum1 = sum1  + (t0 * t1)
               sum2 = sum2 + t1
            endif
         end do

         delta = sum1 / sum2
         xbl = xbl + delta
         icnt=icnt+1
      end do

      sum1 = 0.0
      sum2 = 0.0
      cmad = 9.0 * xmad
      cmadsq = cmad * cmad

      do i=1,n
         t0 = x(i) - xbl
         if (abs(t0) .lt. cmad) then
            t0 = t0 * t0
            t1 = cmadsq - t0
            sum1 = sum1 + (t0 * t1 * t1 * t1 * t1)
            sum2 = sum2 + t1 * (cmadsq - 5.0 * t0)
         endif
      end do

      xbs = n * sqrt(sum1/(n-1.0)) / abs(sum2)

      return
      end




c------------------------------------------------------------------------------
      subroutine medmad (x,n,xmed,xmad)
c------------------------------------------------------------------------------

c  This routine calculates the median and the median absolute deviation of
c  a data set.  The data are supplied in the array X, of length N. The
c  median is returned in XMED, and the median absolute deviation in XMAD.
c
c  If N is odd, the median is the value from the data set that has equal
c  numbers of values greater than it and less than it.  If N is even, the
c  median is the mean of the two central values of the data set.
c
c  The median absolute deviation, as the name implies, is the median of
c  the distribution of absolute values of the deviation of the data about
c  the median.
c
c  WARNING: This routine has the side effect of sorting the input array X.


      integer i, j, k, n, n2
      real x(n), xmed, xmad, xi, xj
      logical even

      even = mod(n,2) .eq. 0
      n2 = n / 2

c  sort the data

      call sort (n,x)

c  calculate the median

      if (even) then
         xmed = 0.5 * (x(n2) + x(n2+1))
      else
         xmed = x(n2+1)
      endif

c  calculate the mad

      i = n2
      j = n2 + 2
      if (even) j = n2 + 1
      xi = xmed - x(i)
      xj = x(j) - xmed

      do k=1,n2
         if (xi .lt. xj) then
            xmad = xi
            if (i .gt. 1) then
               i = i - 1
               xi = xmed - x(i)
            else
               j = j + 1
               xj = x(j) - xmed
            end if
         else
            xmad = xj
            if (j .lt. n) then
               j = j + 1
               xj = x(j) - xmed
            else
               i = i - 1
               xi = xmed - x(i)
            end if
         end if
      end do

      if (even) then
         if (xi .lt. xj) then
            xmad = 0.5 * (xmad + xi)
         else
            xmad = 0.5 * (xmad + xj)
         end if
      end if

      return
      end





C
C    REAL FUNCTION HLQEST
C
C    PURPOSE       COMPUTES THE HODGES-LEHMANN LOCATION ESTIMATOR:
C                  MEDIAN OF ( X(I) + X(J) ) / 2   FOR 1 LE I LE J LE N
C
C    USAGE         RESULT = HLQEST(X,N,LB,RB,Q)
C
C    ARGUMENTS  X   REAL ARRAY OF OBSERVATIONS  (INPUT)
C                 * VALUES OF X MUST BE IN NONDECREASING ORDER *
C
C               N   INTEGER NUMBER OF OBSERVATIONS  (INPUT)
C                 * N MUST NOT BE LESS THAN 1 *
C
C               LB  INTEGER ARRAY OF LENGTH N FOR WORKSPACE
C
C               RB  INTEGER ARRAY OF LENGTH N FOR WORKSPACE
C
C               Q   INTEGER ARRAY OF LENGTH N FOR WORKSPACE
C
C         NOTE ---  ONLY LB,RB, AND Q ARE CHANGED IN COMPUTATION
C
C   EXTERNAL ROUTINE
C              RAN  FUNCTION PROVIDING UNIFORM RANDOM VARIABLES
C                   IN THE INTERVAL (0,1)
C                   RAN REQUIRES A DUMMY INTEGER ARGUMENT
C
C   NOTES           HLQEST HAS AN EXPECTED TIME COMPLEXITY ON
C                   THE ORDER OF N * LG( N )
C
C  J F MONAHAN, APRIL 1982, DEPT OF STAT, N C S U, RALEIGH, N C 27650
C  FINAL VERSION  JUNE 1983
C

      real function fast_hodges_lehmann (x, n, lb, rb, q)

C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
C.. In/Out Status: Maybe Read, Not Written ..
      real x(*)
C.. In/Out Status: Read, Not Written ..
      integer n
C.. In/Out Status: Maybe Read, Maybe Written ..
      integer lb(*)
C.. In/Out Status: Maybe Read, Maybe Written ..
      integer rb(*)
C.. In/Out Status: Maybe Read, Maybe Written ..
      integer q(*)
C
C.. Local Scalars ..
      integer i,ipiq,iq,j,k,k1,k2,l,lbi,mdll,mdlrow,mdlu,nn,rbi,sm,sq
      real am,amn,amx
C
C.. Intrinsic Functions ..
      intrinsic int, max, min, real

      real      ranran
      external  ranran
C
C
C  TAKE CARE OF SPECIAL CASES: N=1 AND N=2
C
      if (n .gt. 2) then
C
C  FIND THE TOTAL NUMBER OF PAIRS (NN) AND THE MEDIAN(S) (K1,K2) NEEDED
C
         nn = (n*(n+1)) / 2
         k1 = (nn+1) / 2
         k2 = (nn+2) / 2
C
C  INITIALIZE LEFT AND RIGHT BOUNDS
C
         do i = 1,n
            lb(i) = i
            rb(i) = n
         end do
C  SM = NUMBER IN SET S AT STEP M
         sm = nn
C  L = NUMBER OF PAIRS LESS THAN THOSE IN SET S AT STEP M
         l = 0
C
C
C  USE THE MEDIAN OF X(I)'S TO PARTITION ON THE FIRST STEP
C
         mdll = (n+1) / 2
         mdlu = (n+2) / 2
         am = x(mdll) + x(mdlu)
         do 20 while (.true.)
C
C       *****   PARTITION STEP   *****
C
C  USE AM TO PARTITION S0 INTO 2 GROUPS: THOSE .LT. AM, THOSE .GE. AM
C  Q(I)= HOW MANY PAIRS (X(I)+X(J)) IN ROW I LESS THAN AM
            j = n
C                              START IN UPPER RIGHT CORNER
            sq = 0
C                              I COUNTS ROWS
            do 10 i = 1,n
               q(i) = 0
               do while (.true.)
C                              HAVE WE HIT THE DIAGONAL ?
                  if (j .lt. i) goto 10
C                              SHALL WE MOVE LEFT ?
                  if (x(i)+x(j) .lt. am) goto 1000
                  j = j - 1
               end do
C                              WE'RE DONE IN THIS ROW
 1000          q(i) = j - i + 1
C  SQ = TOTAL NUMBER OF PAIRS LESS THAN AM
               sq = sq + q(i)
 10         continue
C
C  ***  FINISHED PARTITION --- START BRANCHING  ***
C
C  IF CONSECUTIVE PARTITIONS ARE THE SAME WE PROBABLY HAVE TIES
            if (sq .ne. l) then
C
C  ARE WE NEARLY DONE, WITH THE VALUES WE WANT ON THE BORDER?
C  IF(WE NEED  MAX OF THOSE .LT. AM -OR- MIN OF THOSE .GE. AM) GO TO 90
C
               if (sq .eq. k2-1) goto 1200
C
C  THE SET S IS SPLIT, WHICH PIECE DO WE KEEP?
C  70  =  CUT OFF BOTTOM,   90  =  NEARLY DONE,   60  =  CUT OFF TOP
C
               if (sq .lt. k1) then
                  do i = 1,n
C                            RESET LEFT BOUNDS FOR EACH ROW
                     lb(i) = i + q(i)
                  end do
               else
                  if (sq .eq. k1) goto 1200
                  do i = 1,n
C                            RESET RIGHT BOUNDS FOR EACH ROW
                     rb(i) = i + q(i) - 1
                  end do
               end if
C
C  COUNT   SM = NUMBER OF PAIRS STILL IN NEW SET S
C           L = NUMBER OF PAIRS LESS THAN THOSE IN NEW SET S
               l = 0
               sm = 0
               do i = 1,n
                  l = l + lb(i) - i
                  sm = sm + rb(i) - lb(i) + 1
               end do
C
C        *****   NORMAL RESTART JUMP   *****
C
C  CAN ONLY GET TO 2 LEFT IF K1.NE.K2  -- GO GET THEIR AVERAGE
               if (sm .gt. 2) then
C                        USE RANDOM ROW MEDIAN AS PARTITION ELEMENT
C
C   *****   RESTART HERE UNLESS WORRIED ABOUT TIES   *****
C
                  k = int(real(sm)*ranran(sm))
C                        K IS A RANDOM INTEGER FROM O TO SM-1
                  do i = 1,n
                     j = i
                     if (k .le. rb(i)-lb(i)) goto 1100
                     k = k - rb(i) + lb(i) - 1
                  end do
C                        J IS A RANDOM ROW --- NOW GET ITS MEDIAN
 1100             mdlrow = (lb(j)+rb(j)) / 2
                  am = x(j) + x(mdlrow)
                  goto 20
               end if
            end if
C
C  USE THE MIDRANGE OF SET S AS PARTITION ELEMENT WHEN TIES ARE LIKELY
C   -- OR GET THE AVERAGE OF THE LAST 2 ELEMENTS
C
            amx = x(1) + x(1)
            amn = x(n) + x(n)
            do i = 1,n
C   SKIP THIS ROW IF NO ELEMENT IN IT IS IN SET S ON THIS STEP
               if (lb(i) .le. rb(i)) then
                  lbi = lb(i)
C                             GET THE SMALLEST IN THIS ROW
                  amn = min(amn,x(lbi)+x(i))
                  rbi = rb(i)
C                             GET THE LARGEST IN THIS ROW
                  amx = max(amx,x(rbi)+x(i))
               end if
            end do
            am = (amx+amn) / 2.
C  BE CAREFUL TO CUT OFF SOMETHING -- ROUNDOFF CAN DO WIERD THINGS
            if (am.le.amn .or. am.gt.amx) then
               am = amx
            end if
C  UNLESS FINISHED, JUMP TO PARTITION STEP
            if (amn.eq.amx .or. sm.eq.2) then
C  ALL DONE IF ALL OF S IS THE SAME -OR- IF ONLY 2 ELEMENTS ARE LEFT
               fast_hodges_lehmann = am / 2.
               return
            end if
 20      continue
C
C  FIND:   MAX OF THOSE .LT. AM
C          MIN OF THOSE .GE. AM
 1200    amn = x(n) + x(n)
         amx = x(1) + x(1)
         do i = 1,n
            iq = q(i)
            ipiq = i + iq
            if (iq .gt. 0) then
               amx = max(amx,x(i)+x(ipiq-1))
            end if
            ipiq = i + iq
            if (iq .lt. n-i+1) then
               amn = min(amn,x(i)+x(ipiq))
            end if
         end do
         fast_hodges_lehmann = (amn+amx) / 4.
C  WE ARE DONE, BUT WHICH SITUATION ARE WE IN?
         if (k1 .ge. k2) then
            if (sq .eq. k1) then
               fast_hodges_lehmann = amx / 2.
            end if
            if (sq .eq. k1-1) then
               fast_hodges_lehmann = amn / 2.
            end if
         end if
      else
         fast_hodges_lehmann = x(1)
         if (n .ne. 1) then
            fast_hodges_lehmann = (x(1)+x(2)) / 2.
         end if
      end if
      end


C
C
      real function ranran (ixx)
C  UNIFORM PSEUDORANDOM NUMBER GENERATOR
C  FORTRAN VERSION OF LEWIS, GOODMAN, MILLER
C  SCHRAGE,  ACM TOMS V.5 (1979) P132
C  FIRST CALL SETS SEED TO IXX, LATER IXX IGNORED
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
C.. In/Out Status: Maybe Read, Not Written ..
      integer ixx
C
C.. Local Scalars ..
      integer a,b15,b16,fhi,ix,k,leftlo,p,xalo,xhi
C
C.. Intrinsic Functions ..
      intrinsic real
C
C.. Data Declarations ..
      data a/16807/ b15/32768/ b16/65536/ p/2147483647/
      save a, b15, b16, p

      data ix/0/
      save ix
C
C ... Executable Statements ...
C
      if (ix .eq. 0) then
         ix = ixx
      end if
      xhi = ix / b16
      xalo = (ix-xhi*b16) * a
      leftlo = xalo / b16
      fhi = xhi*a + leftlo
      k = fhi / b15
      ix = (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16) + k
      if (ix .lt. 0) then
         ix = ix + p
      end if
      ranran = real(ix) * 4.656612875e-10
      end
C
C






      subroutine  ms_simple_myriad  (datavec, ndata, scaleconstant,
     $                               myriadvalue, myriaderror)
c     ------------------------------------------------------------

      implicit          none

      integer        ndata
      real           datavec (ndata)

      real           scaleconstant , scaleconst
      real           myriadvalue , myriaderror , datamedian
      real           myriadcost, ave , sdev
      integer        myriadindex


c     /* set scaleconstant */

      scaleconst = scaleconstant
      if (scaleconstant .lt. 0.0)   scaleconst = 0.01

c     /* use robust sigma as myriad error */

      datamedian = 0.0
      call   ms_robust_sigma (datavec, ndata, datamedian,
     $                        .false., .true., myriaderror)
      myriaderror = myriaderror / sqrt (float(ndata))

      if (myriaderror .lt. 0.0)
     $   then
            call   ms_meansigma (datavec, ndata, ave, sdev)
            myriaderror = sdev
         endif

c     /* compute sample myriad */

      call  ms_simplemyriadcost  (datavec, ndata, scaleconst,
     $                            myriadvalue, myriadindex, myriadcost)

      return
      end







      subroutine  ms_simplemyriadcost   (data, ndata, scaleconstant,
     $                            myriadvalue, myriadindex, myriadcost)
c     -------------------------------------------------------------------

      implicit       none

      integer        ndatamax
      parameter      (ndatamax = 10000)

      real           costfunction (ndatamax)

      integer        ndata
      real           data (ndata)

      real           scaleconstant , scaleconstant2
      real           maxcost, myriadvalue, myriadcost , cost
      integer        maxindx , myriadindex
      integer        i , j

      character * 1     carriage return


      carriage return = char (13)

      if (ndata .gt. ndatamax)
     $   then
           print * , ' <ms_simplemyriadcost> ... workspace exceeded !'
           stop
         end if

      scaleconstant2 = scaleconstant ** 2

      if (scaleconstant .gt. 0.0)
     $   then
            do i = 1 , ndata , 1
               costfunction (i) = 0.0
            end do
            do i = 1 , ndata , 1

               if (mod (i,1000) .eq. 0)
     $              then
                       print 1122 , i , carriage return
 1122                  format  ('                          ' //
     $                          '                . . . running ' //
     $                          'data point ' , i5 , ' ' , a1 , $)
                    endif

               do j = 1 , ndata , 1
                  cost = scaleconstant2 + ((data(j)-data(i))**2)
                  if (cost .gt. 0.0) costfunction(i) = costfunction(i) +
     $                                                 alog10(cost)
               end do
               if (i .eq. 1)
     $            then
                     maxcost = costfunction (i)
                     maxindx = i
                  else
                     if (costfunction (i) .lt. maxcost)
     $                  then
                           maxindx = i
                           maxcost = costfunction (i)
                        endif
                  endif
            end do
         else
            do i = 1 , ndata , 1
               costfunction (i) = 0.0
            end do
            do i = 1 , ndata , 1

               if (mod (i,1000) .eq. 0)
     $              then
                       print 1122 , i , carriage return
                    endif

               do j = 1 , ndata , 1
                  cost = (data(j)-data(i))**2
                  if (cost .ne. 0.0) costfunction(i) = costfunction(i) +
     $                                                 cost
               end do
               if (i .eq. 1)
     $            then
                     maxcost = costfunction (i)
                     maxindx = i
                  else
                     if (costfunction (i) .lt. maxcost)
     $                  then
                           maxindx = i
                           maxcost = costfunction (i)
                        endif
                  endif
            end do
         endif

      Myriadindex = maxindx
      Myriadvalue = Data (Myriadindex)
      Myriadcost  = Costfunction (Myriadindex)

      return
      end






      subroutine    ms_meansigma   (vec , n , aver , sdev)
c     ----------------------------------------------------

c     /*  find mean and standard deviation of vector elements  */

      implicit          none

      real        vec (*) , aver , sdev
      real        s , vari , sum1 , sum2
      integer     n , i

      aver = 0.0
      sdev = 0.0
      if (n .le. 0)  return

      do 10000  i = 1 , n , 1
         aver = aver + vec (i)
10000 continue

      aver = aver / float (n)

      sum1 = 0.0
      sum2 = 0.0
      do  i = 1 , n , 1
         s = aver - vec (i)
         sum1 = sum1 + (s * s)
         sum2 = sum2 + s
      end do

      vari = (sum1 - (sum2 / real (n))) / real (n-1)
      sdev = sqrt (vari)

      return
      end







      subroutine   ms_robust_sigma  (xdata, ndata, medianvalue,
     $                               basezero, domedian, robustsigma)
c     ---------------------------------------------------------------

      implicit          none

      integer     ndatamax
      parameter      (ndatamax = 10000)

      real        tmp1 (ndatamax)

      integer     i , ndata , nlow
      real        xdata (ndata) , medianvalue , robustsigma
      logical     basezero , domedian

      real        epsilon , xmed , x0 , medmad , meanmad
      real        numerator , denominator , sigmavalue


      if (ndata .gt. ndatamax)
     1  then
            print * , ' <ms_robust_sigma> ... workspace exceeded !'
            stop
         end if

      epsilon = 1.0E-20
      if (basezero)
     $   then
            x0 = 0.0
         else
            if (domedian)
     $         then
                  call  ms_median1 (xdata, ndata, xmed)
                  x0 = xmed
               else
                  x0 = medianvalue
               end if
         end if

      do i = 1 , ndata , 1
         tmp1 (i) = abs (xdata(i) - x0)
      end do
      call  ms_median1 (tmp1, ndata, medmad)
      medmad = medmad / 0.6745

      if (medmad .lt. epsilon)
     $   then
            call  ms_mean (tmp1, ndata, meanmad)
            medmad = meanmad / 0.8
         end if

      if (medmad .lt. epsilon)
     $   then
            robustsigma = 0.0
            return
         end if

      nlow        = 0
      numerator   = 0.0
      denominator = 0.0
      do i = 1 , ndata , 1
         tmp1 (i) = ((xdata(i) - x0) / (6.0 * medmad))**2
         if (tmp1 (i) .le. 1.0)
     $      then
               nlow        = nlow + 1
               numerator   = numerator   +
     $                       ((1.0-tmp1(i))**4) * ((xdata(i) - x0)**2)
               denominator = denominator +
     $                       ((1.0-tmp1(i)) * (1.0 - (5.0*tmp1(i))))
            end if
      end do

      if (nlow .le. 3)
     $   then
            print * , ' <ms_robust_sigma> ... weird distribution !'
            print * , ' <ms_robust_sigma>     returning -1 !'
            robustsigma = -1.0
            return
         end if

      sigmavalue = real (nlow) * numerator /
     $             (denominator*(denominator-1.0))

      if (sigmavalue .gt. 0.0)   sigmavalue = sqrt (sigmavalue)

      robustsigma = sigmavalue

      return
      end









      subroutine  ms_median1 (x,n,xmed)
c     ---------------------------------

      integer     n, n2
      real        x (n) , xmed

      call sort(n,x)
      n2=n/2
      if(2*n2.eq.n)then
        xmed=0.5*(x(n2)+x(n2+1))
      else
        xmed=x(n2+1)
      endif
      return
      end






      subroutine   ms_mean   (vec, n, aver)
c     -------------------------------------

c     /*  find mean of vector elements  */

      implicit          none

      real        vec (*) , aver
      integer     n , i

      aver = 0.0
      if (n .le. 0)  return

      do 10000  i = 1 , n , 1
         aver = aver + vec (i)
10000 continue

      aver = aver / real (n)

      return
      end




