C @(#)st_def.inc	18.1.1.1 (ESO-IPG) 09/03/02 18:37:44
C 
C 
C --------------------- MIDAS definition Table ---------------------------
C 	
C   to declare the variables
C 
C   which are then really assigned to values in st_dat.inc
C 
C   880414:  creation  (KB)
C   010410		last modif
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
      INTEGER F_FITS_PARM                  ! set flags for FITS header proc.
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
C 
      INTEGER X_Y_PLANE  
      INTEGER X_Z_PLANE 
      INTEGER Z_Y_PLANE
C 
C  ------------------------------------------------------------------  
