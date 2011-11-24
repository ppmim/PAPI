C @(#)st_dat.inc	6.1.1.1 (ESO-IPG) 7/16/93 20:23:48
C 
C 
C --------------------- MIDAS definition Table ---------------------------
C 	
C                       on file ST_DAT.INC
C 
C  to set the variables to specific values
C  version [2.00]  880414:  creation  (KB)
C  version [2.20]  880908:  add filetype no's (KB)
C  version [2.30]  900607:  update according to midas_def.h  (KB)
C  version [2.35]  910125:  add F_OLD_TYPE, F_ASC_TYPE (KB)
C 
C 
      DATA F_O_MODE   /1/
      DATA F_IO_MODE  /2/
      DATA F_I_MODE   /0/
      DATA F_U_MODE   /2/
      DATA F_X_MODE   /9/
      DATA F_D_MODE   /2/
C
C 
      DATA F_XD_PARM  /0/
C
C 
      DATA F_OLD_TYPE  /0/
      DATA F_IMA_TYPE  /1/
      DATA F_TBL_TYPE  /3/
      DATA F_FIT_TYPE  /4/
      DATA F_ASC_TYPE  /9/
C
C 
      DATA D_OLD_FORMAT  /0/
      DATA D_I1_FORMAT   /1/
      DATA D_I2_FORMAT   /2/
      DATA D_UI2_FORMAT  /102/
      DATA D_I4_FORMAT   /4/
      DATA D_R4_FORMAT   /10/
      DATA D_R8_FORMAT   /18/
      DATA D_L1_FORMAT   /21/
      DATA D_L2_FORMAT   /22/
      DATA D_L4_FORMAT   /24/
      DATA D_C_FORMAT    /30/
      DATA D_X_FORMAT    /40/
      DATA D_P_FORMAT    /50/
C
C 
      DATA F_TRANS       /0/
      DATA F_RECORD      /1/
C 
C  ------------------------------------------------------------------  
