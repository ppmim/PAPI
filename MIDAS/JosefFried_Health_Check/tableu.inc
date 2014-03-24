C @(#)tableu.inc	6.1.1.1 (ESO-IPG) 7/16/93 20:24:23
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C         Control Data Structure Compile Expressions
C
C       TABLEU.INC
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
C
C ... INITIALIZE DATA IN COMMON BLOCK
C
       DATA              EXPNFU/16/
       DATA              EXPNSY/21/
       DATA              EXPNAR / 16*1,24*0/
       DATA              EXPIFU
     . / 'SQRT  ','LN    ','LOG10 ','EXP   ','SIN   ','COS   ',
     .   'TAN   ','ASIN  ','ACOS  ','ATAN  ','SINH  ','COSH  ',
     .   'TANH  ','ABS   ','INT   ','ANINT ',      24*'      '/
       DATA              EXPREL
     . /'.LE.','.LT.','.GE.','.GT.','.EQ.','.NE.'/
       DATA              EXPSYM
     . /'@','`',';','|','<','{','>','}','~','+',
     .  '-','*','/','^','P','F','C','L','D','(',')'/
C     . /'o','a','n','l','<','g','>','e','~','+',
C     .  '-','*','/','p','P','F','C','L','D','(',')'/
       DATA              EXPPRE
     . / 1,  3,  5,  7,  7,  7,  7,  7,  7,  9,
     .   9, 11, 11, 14, 15, 15, 15, 15, 15, 17,  0,
     .   2,  4,  6,  8,  8,  8,  8,  8,  8, 10,
     .  10, 12, 12, 13, 16, 16, 16, 16, 16,  0,  0,
     .  -1, -1,  0, -1, -1, -1, -1, -1, -1, -1,
     .  -1, -1, -1, -1,  0,  1,  1,  1,  1,  0,  0/
C
C -----------------------------------------------------------------------
