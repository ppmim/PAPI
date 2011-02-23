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
.COMMENTS		Input frame must be a .BDF frame that was opened in read-write-mode
			(F_IO_MODE). The normalization level is NORM_LEVEL (in pipe_para.h).
			The input frame is expected to have a descriptor SUB_STATISTICS.
			The median (SUB_STATISTICS(3)) is than taken and scaled to NORM_LEVEL.
			Every pixel of the image is than scaled with the same factor.
			A (real) descriptor NORMALIZATION is created, which stores:
			NORMALIZATION(1)= normalization factor that the image was multiplied by
			NORMALIZATION(2)= NORM_LEVEL (count level the median was scaled to)
.CALL			int Normalize(float *p_inframe, int image_no, int tot_pix)
.VERSION		1.00		22.08.02
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
  /* local variables */
  int fstat;		/* function status */
  int count;		/* loop variable */

  float median;		/* median level taken from descriptor SUB_STATISTICS */
  float norm_factor;  	/* the factor the image is multiplied by */
  float level = NORM_LEVEL;
  float sigma;		/* stddev of image */
  float low, high; 	/* for LHCUTS */


  int actvals;		/* actual number of values returned from function */
  int unit = 0;		/* unused dummy variables for interface */
  int null = 0;



  /* get the median level from SUB_STATISTICS */
  fstat = SCDRDR(image_no,"SUB_STATISTICS",3,1,&actvals,&median,&unit,&null);

  if(actvals == 0)
    {
  fstat = SCTPUT("ERROR in function Normalize(): Descriptor SUB_STATISTICS does not exist");
  return 1;		/* Error status for: Descriptor not present */
    }


  printf("\nMedian = %f\n", median);	// Test



  /* get normalization factor */
  norm_factor = NORM_LEVEL / median;


  /* do normalization of image */
  for(count=1 ; count<=tot_pix ; count++)
    {
      *p_inframe = (*p_inframe) * norm_factor;

      p_inframe++;
    }



  /* write descriptor NORMALIZATION */
  fstat = SCDWRR(image_no,"NORMALIZATION",&norm_factor,1,1,&unit);
  fstat = SCDWRR(image_no,"NORMALIZATION",&level,2,1,&unit);

  fstat = SCDWRH(image_no,"NORMALIZATION","NormalizationFactor(1), MedianCountLevel(2)",1,50 );		/* help text */

  /* update descriptor HISTORY */
  fstat = SCDWRC(image_no,"HISTORY",1,"Normalized ",-1,14,&unit);


  /* set LHCUTS to (-3*sigma, +10*sigma) */
  fstat = SCDRDR(image_no,"SUB_STATISTICS",5,1,&actvals,&sigma,&unit,&null);  		/* get stdev */

  sigma = sigma * norm_factor;		/* scale the stddev */
  low = NORM_LEVEL - 3*sigma;
  high = NORM_LEVEL + 10*sigma;

  fstat = SCDWRR(image_no,"LHCUTS",&low,1,1,&unit);
  fstat = SCDWRR(image_no,"LHCUTS",&high,2,1,&unit);


  /* back to calling function */
  return 0;
}
