/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		flat_correction.h
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.PURPOSE	      Header File for flat_correction.c 
.COMMENTS
.VERSION		13.08.02
--------------------------------------------------------------------------------------*/


/* Function Prototype */
int Flat(float *p_inframe, float *p_flat, float *p_outframe, int total_pix, int image_no_out,\
	 char *name_inframe, char *flat_frame, int dark_flag, float *p_dark, char *dark_name,\
	 double dark_itime, int inframe_imno);


