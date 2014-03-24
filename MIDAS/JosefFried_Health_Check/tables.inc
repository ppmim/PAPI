C @(#)tables.inc	6.1.1.1 (ESO-IPG) 7/16/93 20:24:22
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C.IDENTIFICATION       TABLES.INC
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
C ... File mode IO
C
C       INTEGER              FIMODE
C       INTEGER              FOMODE
C       INTEGER              FUMODE
C       INTEGER              FXMODE
C       INTEGER              FDMODE
C
C ... Data types 
C
C       INTEGER DI1FOR
C       INTEGER DI2FOR
C       INTEGER DI4FOR
C       INTEGER DR4FOR
C       INTEGER DR8FOR
C       INTEGER DL1FOR
C       INTEGER DL2FOR
C       INTEGER DL4FOR
C       INTEGER DCFORM
C       INTEGER DXFORM
C       INTEGER DPFORM
C
C ... Table store format
C
C       INTEGER              FTRANS
C       INTEGER              FRECOR
       INCLUDE             'MID_INCLUDE:ST_DEF.INC'
C       
C ... Error codes at level 2
C
       CHARACTER*64             TPARBF(8)
       INTEGER              ERRILL       
       INTEGER              ERRSEM       
       INTEGER              ERRSOV       
       INTEGER              ERRSOP       
       INTEGER              ERRPAR
       INTEGER              ERRCOD
       INTEGER              ERRFMT
       INTEGER              ERRRED
       INTEGER              ERRSCL
       INTEGER              ERRCOL
       INTEGER              ERRREF
       INTEGER              ERRALG
       INTEGER              ERRGNP 
       INTEGER              ERRGAL 
C
       COMMON/TBLCM1/           TPARBF
C
C -----------------------------------------------------------------------
C
