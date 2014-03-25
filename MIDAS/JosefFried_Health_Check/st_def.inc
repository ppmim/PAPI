C @(#)st_def.inc	6.1.1.1 (ESO-IPG) 7/16/93 20:23:41
C 
C 
C --------------------- MIDAS definition Table ---------------------------
C 	
C                       on file ST_DEF.INC
C 
C   to declare the variables
C   version [2.00]  880414:  creation  (KB)
C   version [2.10]  880908:  add definitions for filetype no's (KB)
C   version [2.20]  890203:  add definitions for table formats (JDP)
C   version [2.30]  900607:  update according to midas_def.h (KB)
C   version [2.35]  910125:  add F_OLD_TYPE, F_ASC_TYPE  (KB)
C 
C 
C   Mode definitions:   
C 
C
      INTEGER F_I_MODE                     ! map file for input only    
      INTEGER F_O_MODE                     ! map file for output        
      INTEGER F_IO_MODE                    ! map file for updating      
      INTEGER F_U_MODE                     ! map file for updating      
      INTEGER F_X_MODE                     ! create/map virtual memory  
      INTEGER F_D_MODE                     ! table descriptor access
C
C 
C   File parameter definitions:   
C 
C
      INTEGER F_XD_PARM                    ! create file with Xlarge dsc area
C
C 
C   File type definitions:   
C 
C
      INTEGER F_OLD_TYPE                   ! old file type
      INTEGER F_IMA_TYPE                   ! image file type
      INTEGER F_TBL_TYPE                   ! table file type
      INTEGER F_FIT_TYPE                   ! fit file type
      INTEGER F_ASC_TYPE                   ! ASCII file type
C
C 
C   Data Format definitions:   
C 
C
C
      INTEGER D_OLD_FORMAT     
      INTEGER D_I1_FORMAT     
      INTEGER D_I2_FORMAT     
      INTEGER D_UI2_FORMAT     
      INTEGER D_I4_FORMAT     
      INTEGER D_R4_FORMAT    
      INTEGER D_R8_FORMAT
      INTEGER D_L1_FORMAT
      INTEGER D_L2_FORMAT
      INTEGER D_L4_FORMAT
      INTEGER D_C_FORMAT
      INTEGER D_X_FORMAT
      INTEGER D_P_FORMAT
C
C     table file format
C
      INTEGER F_TRANS                 ! transposed format (by columns)
      INTEGER F_RECORD                ! record format  
C  ------------------------------------------------------------------  
