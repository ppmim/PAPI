C @(#)tabled.inc	6.1.1.1 (ESO-IPG) 7/16/93 20:24:02
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C.IDENTIFICATION       TABLED.INC
C.AUTHOR              J.D.Ponz       ESO - Garching
C.PURPOSE
C       Defines table related parameters
C
C.VERSION              3.0 6-NOV-1987
C
C.COMMENTS
C
C       ERRILL       ! illegal expression
C       ERRSEM       ! stack empty
C       ERRSOV       ! stack overflow
C       ERRSOP       ! error in operation code
C       ERRPAR
C       ERRCOD
C       ERRFMT
C       ERRRED
C       ERRSCL
C       ERRCOL
C       ERRREF
C       ERRALG
C       ERRGNP ! Regression, No of points
C       ERRGAL ! Regression, algorithm
C
C------------------------------------------------------
C
C ... File access mode
C
C       DATA              FIMODE/0/
C       DATA              FOMODE/1/
C       DATA              FUMODE/2/
C       DATA              FXMODE/9/
C       DATA              FDMODE/2/
C       DATA              FDMODE/10/
C
C ... Table store mode
C
C       DATA              FTRANS/0/
C       DATA              FRECOR/1/
C
C ... Data types 
C
C       DATA DI1FOR/1/
C       DATA DI2FOR/2/
C       DATA DI4FOR/4/
C       DATA DR4FOR/10/
C       DATA DR8FOR/18/
C       DATA DL1FOR/21/
C       DATA DL2FOR/22/
C       DATA DL4FOR/24/
C       DATA DCFORM/30/
C       DATA DXFORM/40/
C       DATA DPFORM/50/
       INCLUDE     'MID_INCLUDE:ST_DAT.INC'
C       
C ... Error codes at level 2
C
       DATA              ERRILL/-100/       
       DATA              ERRSEM/-104/       
       DATA              ERRSOV/-107/       
       DATA              ERRSOP/-108/       
       DATA              ERRPAR/-500/
       DATA              ERRCOD/-501/
       DATA              ERRFMT/-502/
       DATA              ERRRED/-503/
       DATA              ERRSCL/-504/
       DATA              ERRCOL/-505/
       DATA              ERRREF/-506/
       DATA              ERRALG/-507/
       DATA              ERRGNP/-600/ 
       DATA              ERRGAL/-601/ 
C
C -----------------------------------------------------------------------
C

