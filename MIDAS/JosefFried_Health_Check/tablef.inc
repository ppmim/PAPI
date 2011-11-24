C @(#)tablef.inc	6.1.1.1 (ESO-IPG) 7/16/93 20:24:21
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C.IDENTIFICATION           TABLEF.INC
C.AUTHOR                   J.D.Ponz       ESO - Garching
C
C.PURPOSE
C
C  Defines Format structures for FORTRAN programs
C
C.VERSION                  3.0 6-Nov-1987
C
C------------------------------------------------------
       INTEGER              FMNBUF
       PARAMETER           (FMNBUF=256)
       CHARACTER*16         FMLABE (FMNBUF)
       CHARACTER*16         FMUNIT (FMNBUF)
       CHARACTER*5          FMTYPE (FMNBUF)
       CHARACTER*8          FMFORM (FMNBUF)
       INTEGER              FMNCOL
       INTEGER              FMNROW
       INTEGER              FMIDEN
       INTEGER              FMLLEN
       INTEGER              FMFIRS (FMNBUF)
       INTEGER              FMLAST (FMNBUF)
       COMMON/TBLCM2/              FMNCOL,
     -                            FMNROW,
     -                            FMIDEN,
     -                            FMLLEN,
     -                            FMFIRS,
     -                            FMLAST
       COMMON/TBLCM3/               FMFORM,
     -                            FMLABE,
     -                            FMUNIT,
     -                            FMTYPE
C
C -----------------------------------------------------------------------
C
