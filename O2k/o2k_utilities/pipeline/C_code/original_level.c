/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		original_level.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		MODULE 10 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		un_normalization
.PURPOSE	      	set the count level back to original values  
.IN/OUTPUT		IN: pointer to first element of input frame
			IN: image number of input frame
			IN: total number of pixels of frame
.RETURNS		Status: 0 = ok
.COMMENTS		This function sets back the effect of 'normalize_frame.c'.
			The values of the normalization process were:
			NORMALIZATION(1)= normalization factor that the image was multiplied by
		     	NORMALIZATION(2)= NORM_LEVEL (count level the median was scaled to)
			The function will divide the input frame by the normalization factor
			to restore the original count level.
.CALL			int Original(float *p_inframe, int image_no, int tot_pix)
.VERSION		1.00		06.09.02
			1.1		06.03.03	changed all output to SCTPUT
-----------------------------------------------------------------------------------------*/


#include<stdio.h>
#include<midas_def.h>		/* MIDAS definitions */
		
/* Header Files of Secundary Modules */
#include "original_level.h"



int Original(float *p_inframe, int image_no, int tot_pix)
{
  /* local variables */
  int fstat;		/* function status */
  int count;		/* loop variable */
  int actvals;		/* actual number of values returned from function */
  int unit = 0;		/* unused dummy variables for interface */
  int null = 0;

  float norm_factor;  	/* the factor the image was multiplied by */

  char aux_string[200];	/* auxilary string buffer for SCTPUT */


  /* info on screen */
  fstat = SCTPUT("\nSetting count level back...");


  /* get the normalization factor from descriptor NORMALIZATION */
  fstat = SCDRDR(image_no,"NORMALIZATION",1,1,&actvals,&norm_factor,&unit,&null);


  if(actvals == 0)
    {
  fstat = SCTPUT("ERROR in function Original(): Descriptor NORMALIZATION does not exist");
  return 1;		/* Error status for: Descriptor not present */
    }


  sprintf(aux_string,"\nNorm_factor = %f\n", norm_factor);
  fstat = SCTPUT(aux_string);


  /* set count level of image back */
  for(count=1 ; count<=tot_pix ; count++)
    {
      *p_inframe = (*p_inframe) / norm_factor;

      p_inframe++;
    }


  /* update descriptor HISTORY */
  fstat = SCDWRC(image_no,"HISTORY",1,\
		 "Set back to original count level: divided by factor NORMALIZATION(1) ",-1,80,&unit);

 
  /* back to calling function */
  return 0;
}
