!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!... Identification 	:  BACKDET.PRG
!
!... Version		:  4.4	2-Nov-87	HJR
!    Revisions		:  8-Sep-88 HJR   processing without TV
!
!... Purpose		:  set CUTS for frame using STATISTIC-information
!
!... Call	 BACK/DET  frame  minus,plus  window
!
!
!... Action		: image will be displayed (if not yet done) and cursor
!			  is used to define clean background area for STATISTIC
!
!			Set cursor board to :	TRACK off
!						RATE  on
!						CURSOR 1 + 2 on
!
! Parameter definition and explanation :
!
DEFINE/PAR P1 XXX
!			Input frame (default : displayed image)
DEFINE/PAR P2 3,15
!			Cuts in units of sigma. Value will be subtracted from
!			median for lower and added for upper cut value.
!			Default : 3,15 
DEFINE/PAR P3 CURS ?  "Window for statistics : "
!
!-----------------------------------------------------------------------------!
!



if p1(1:4) .eq. "help" then

	write/out @@ backdet [ima] [lfac,hfac]
	write/out determine low/high cuts from statistics within cursor box:
	write/out low  = median-lfac*stddev   high = median+hfac*stddev
	write/out ima = image defaulted to loaded image
	write/out lfac,hfac defaulted to 3,15
	return
endif


define/local DIF/R/1/2 {P2}
define/local CT1/R/1/1 0.
define/local CT2/R/1/1 0.
define/local MED/R/1/1 0.
define/local SIG/R/1/1 0.
define/local FRAME/c/1/24 " "
!
IF P1(1:3) .EQ. "XXX" THEN
	@@ tvdisply
	WRITE/KEY FRAME/C/1/24 {TVDISPLY}
ELSE
	WRITE/KEY FRAME/C/1/24 {P1}
ENDIF
!
IF P3(1:4) .EQ. "CURS"  THEN
	LOG/OFF
	IF P1 .NE. "XXX" THEN
		LOAD {FRAME}
	ENDIF
	write/out "Click left   mouse button for ENTER"
	write/out "Click center mouse button for EXIT"
	STAT/IMA ? CURSOR ? ? FS
	LOG/ON
ELSE
	STAT/IMA {tvdisply} {P3} ? ? FS
ENDIF
!
med = outputr(8)
sig = outputr(4)
COMPUTE/KEY CT1 = MED-DIF(1)*SIG
COMPUTE/KEY CT2 = MED+DIF(2)*SIG
CUTS {FRAME} {CT1},{CT2}
!
IF P3(1:4) .EQ. "CURS"  THEN
	CLEAR/CHAN OVERLAY
	LOAD {FRAME}
ELSE
	CLEAR/CHAN OVERLAY
	LOAD {FRAME}
ENDIF
!
RETURN
