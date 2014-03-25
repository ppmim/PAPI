C @(#)tablet.inc	6.1.1.1 (ESO-IPG) 7/16/93 20:24:23
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C         Control Data Structure Compile Expressions
C
C       TABLET.INC
C
C       Stack Handling
C
C       EXPPTR              Top of stack
C       EXPOVF                     Overflow control
C       EXPCAA                     Character stack
C       EXPPOO                     Integer stack
C
C       Polish Conversion
C
C       EXPNFU              Number of allowed functions (16)
C        EXPIFU              Function names
C       EXPNAR              Number of arguments of the functions
C       EXPNSY              Number of allowed symbols (16)
C       EXPSYM              List of allowed symbols
C       EXPPRE              Table with precedence rules
C
C       Operation Control
C
C       EXPPR1              Buffer with polish expression
C       EXPPR2              Auxiliar buffer
C       EXPATO              List of atoms
C       EXPNAT              Number of atoms in list
C       EXPLAT              Length of atoms
C       EXPOPE              Operation to execute
C       EXPOPC              Operation control
C       EXPRES              Operation result
C       EXPIPT              Pointers to atoms
C       EXPOPT              Pointers to atoms
C       
C       J.D.Ponz       ESO - Garching  25may82
C
C---------------------------------------------------------
       INTEGER              EXPPTR
       INTEGER              EXPOVF
       INTEGER              EXPPOO(40)
       INTEGER              EXPNFU
       INTEGER              EXPNAR(40)
       INTEGER              EXPNSY
       INTEGER              EXPPRE(21,3)
       INTEGER              EXPOPC(6)
       INTEGER              EXPNAT
       INTEGER              EXPLAT(40)
       INTEGER              EXPIPT(40)
       INTEGER              EXPOPT(40)
       INTEGER              EXPCNU
       INTEGER              EXPVNU
       CHARACTER*1              EXPCAA(40)
       CHARACTER*6              EXPIFU(40)
       CHARACTER*1              EXPSYM(21)
       CHARACTER*4              EXPREL(6)
       CHARACTER*17              EXPATO(40)
       CHARACTER*17              EXPOPE
       CHARACTER*128              EXPPR1
       CHARACTER*128              EXPPR2
       CHARACTER*17              EXPRES
C
       COMMON/EXPCO1/              EXPPTR,
     .                            EXPOVF,
     .                            EXPPOO,
     .                            EXPNFU,
     .                            EXPNAR,
     .                            EXPNSY,
     .                            EXPPRE,
     .                            EXPIPT,
     .                            EXPOPT,
     .                            EXPOPC,
     .                            EXPNAT,
     .                            EXPLAT,
     .                            EXPVNU,
     .                            EXPCNU
       COMMON/EXPCO2/              EXPCAA,
     .                            EXPSYM,
     .                            EXPIFU,
     .                            EXPREL,
     .                            EXPATO,
     .                            EXPOPE,
     .                            EXPRES,
     .                            EXPPR1,
     .                            EXPPR2
C
C ... data definition in (routine TEPOLS old)
C ... data definition in TABLEU.INC       
C
