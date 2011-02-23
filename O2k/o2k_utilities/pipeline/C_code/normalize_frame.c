/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		normalize_frame.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		MODULE 3 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		normalization
.PURPOSE	      	normalize the median value of a frame to given level  
.IN/OUTPUT		IN: pointer to first element of input frame
			IN: image number of input frame
			IN: total number of pixels of frame
.RETURNS		Status: 0 = ok
.COMMENTS		The normalization level is set to TOT_IMAGES * median_of_first_image,
			which is roughly the count level of the sum of the images on the stack.
 			The input frame is expected to have a descriptor SUB_STATISTICS.
			The median (SUB_STATISTICS(3)) is than taken and scaled to the above 
			level. Every pixel of the image is than multiplied by this factor.
			A (real) descriptor NORMALIZATION is created, which stores:
			NORMALIZATION(1)= normalization factor that the image was multiplied by
			NORMALIZATION(2)= NORM_LEVEL (count level the median was scaled to)
.CALL			int Normalize(float *p_inframe, int image_no, int tot_pix)
.VERSION		1.00		22.08.02
			1.1		05.09.02	The normalization level is determined
			     		interactivlely during the first call of the function.
					LEVEL is set to: TOT_IMAGES*median_of_first_image,
					which is about the sum_level of all images on the stack 
			1.2		23.10.02	adjustment to external TOT_IMAGES
			1.3		14.11.02	update HISTORY and cancel LHCUTS setting
			1.4		06.03.03	changed all output to SCTPUT
--------------------------------------------------------------------------------------*/


#include<stdio.h>
#include<stdlib.h>
#include<midas_def.h>		/* MIDAS definitions */
		
/* Header Files of Secundary Modules */
#include "normalize_frame.h"

/* Pipeline Parameters, including NORM_LEVEL */
#include "pipe_para.h"




int Normalize(float *p_inframe, int image_no, int tot_pix)
{
  /* make global variables visible */
  extern int TOT_IMAGES;				


  /* local variables */
  int fstat;		/* function status */
  int count;		/* loop variable */

  static float level = 0;	/* !!!!! different from version 1.00 !!!!! */

  float median;		/* median level taken from descriptor SUB_STATISTICS */
  float norm_factor;  	/* the factor the image is multiplied by */

  int actvals;		/* actual number of values returned from function */
  int unit = 0;		/* unused dummy variables for interface */
  int null = 0;

  char aux_string[200];	/* auxilary buffer for MIDAS output */



  /* write info on screen */
  fstat =  SCTPUT("Normalization is done ...");


  /* get the median level from SUB_STATISTICS */
  fstat = SCDRDR(image_no,"SUB_STATISTICS",3,1,&actvals,&median,&unit,&null);

  if(actvals == 0)
    {
  fstat = SCTPUT("ERROR in function Normalize(): Descriptor SUB_STATISTICS does not exist");
  return 1;		/* Error status for: Descriptor not present */
    }


  /* if this is first function call, initialize normalization level */
  if(level == 0)
    {
      level = median * TOT_IMAGES;

      sprintf(aux_string,"\nThe normalization level was set to: %f counts\n", level);
      fstat = SCTPUT(aux_string);
    }
  


  /* get normalization factor */
  norm_factor = level / median;


  /* do normalization of image */
  for(count=1 ; count<=tot_pix ; count++)
    {
      *p_inframe = (*p_inframe) * norm_factor;

      p_inframe++;
    }



  /* write descriptor NORMALIZATION */
  fstat = SCDWRR(image_no,"NORMALIZATION",&norm_factor,1,1,&unit);
  fstat = SCDWRR(image_no,"NORMALIZATION",&level,2,1,&unit);

  /* help text for descriptor*/
  fstat = SCDWRH(image_no,"NORMALIZATION","NormalizationFactor(1), MedianCountLevel(2)",1,50 );

  /* update descriptor HISTORY */
  fstat = SCDWRC(image_no,"HISTORY",1,"Frame normalized as specified in descriptor NORMALIZATION ",-1,80,&unit);


  /* back to calling function */
  return 0;
}
