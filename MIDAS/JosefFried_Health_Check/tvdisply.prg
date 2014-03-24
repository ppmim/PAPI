!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!									   !
!... Identification:  TVDISPLY.PRG     					   !
!									   !
!... Version:	4.4	date: 09-Mar-88   author: HJR	  		   !
!									   !
!... Purpose: Fill keyword TVDISPLY with the name of the image displayed   !
!		on the TV						   !
!									   !
!... Call: @@ RM:TVDISPLY						   !
!									   !
!... Action to be taken by user during execution : none			   !
!									   !
!... Requirements for execution :					   !
!									   !
!           PRG        EXE        COM        DATA       keywords	   !
!									   !
!									   !
!... Parameter definition and explanation :				   !
!									   !
!	no parameters							   !
!									   !
!--------------------------------------------------------------------------!
!									   !
!
define/local IHELP1/I/1/1 0
define/local IHELP2/I/1/1 0
WRITE/KEY TVDISPLY/C/1/24 " "
IHELP1 = IDIDEV(15)*20+1
IHELP2 = IHELP1 + 19
set/format i2
!
! WRITE/KEY TVDISPLY/C/1/24 {IDIMEMC({IHELP1}:{IHELP2})}
WRITE/KEY TVDISPLY/C/1/24 {IDIMEMC(1:20)}
!
RETURN
