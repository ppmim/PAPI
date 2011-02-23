/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		sum_preparation.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 14 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		rejection, normalization, copying 
.PURPOSE	     	do preparation steps for cosmics filtering and summation
.IN/OUTPUT		IN: p_frame = pointer to input frame
			IN: n_average = number of frames for cosmics median cleaning
.RETURNS		1 --> frame was rejected from summation
			Status: 0 = ok
.COMMENTS		When the function Preparation() is called, it will take the frame
			information of the top-level in sum.info and do the following:
			1) check wether the frame belongs to the same ditherpattern, if not
			--> put frame on rejection stack
			2) determine the approximate pixel offsets from the observing position
			   and store them in sum.info
			3) copy frame to the clean_stack and do median normalization to the
			   count_level specified in NORM_LEVEL (pipe_para.h)
.CALL			int Preparation(float *p_frame, int n_average)
.VERSION		1.00   		26.11.02
			1.1			02.12.02	adjusted RA to be in hours
			1.2			04.12.02	take cosmics_stack out of clean.
			1.3			06.03.03	changed all output to SCTPUT
			1.4			07.04.03	save original offsets
			2.0			17.04.04	do telescope drift correction
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<string.h>
#include<math.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "sum_preparation.h"

/* Pipeline related parameters/constants */
#include "pipe_para.h"		

/* Header file with structure definitions */
#include "struct_def.h"

/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif

/* make gloabel variables accessible */
extern DEEP_STRUCT sum;
extern COSMICS_STACK clean;					
extern SUM_PAGE rejection[];				
extern int reject_count;	
extern float cosmics_stack[N_AVERAGE_MAX][PIX_AXIS][PIX_AXIS];



int Preparation(float *p_frame, int n_average)
{
  /* local variables */
  int fstat;
  int last_pos;						/* position of relevant info on info stack */
  int loop_x , loop_y;

  double del_alpha;		       			/* difference in observing position ...*/
  double del_dec;					/* ...for rejection */
  double PI = acos(-1);

  float norm_factor;					/* normalization factor */

  char aux_string[200];					/* string buffer */


  /* initialize last_pos */
  last_pos = sum.total_opened - 1;

  fstat = SCTPUT("\nPreperation for summation process...");


  /* check wether last frames belongs to the same ditherpattern */
  if(sum.total_opened > 1)	/* at least second image */
    {
      /* get difference in observing position in arc_minutes */
      del_alpha = sum.info[last_pos].alpha - sum.reference.alpha;
      del_dec = sum.info[last_pos].declination - sum.reference.declination;

      /* conversion degrees --> arc_minutes */
      del_alpha = 15 * 60 * del_alpha * cos(PI * sum.info[last_pos].declination / 180);	
      /* for RA in hours */
      del_dec = 60 * del_dec;

      /* take absolute value */
      del_alpha = fabs(del_alpha);
      del_dec = fabs(del_dec);


      /* do rejection of observing position is too for away */
      if(del_alpha > DELTA_REJECTION  ||  del_dec > DELTA_REJECTION)
	{
	  fstat = SCTPUT("\n\nWARNING (in summation process)");
	  fstat = SCTPUT("The following image belongs to a different dither pattern and is rejected: ");
	  fstat = SCTPUT(sum.info[last_pos].ima_name);

	  /* copy info to rejection stack */
	  rejection[reject_count] = sum.info[last_pos];

	  /* get pixel offsets */
	  rejection[reject_count].delta_x = 60 * del_alpha / PIX_SCALE;
	  rejection[reject_count].delta_y = 60 * del_dec / PIX_SCALE;

	  /* display rejection */
	  sprintf(aux_string,"x,y-pixel-offsets were to big: %f, %f\n", \
		  rejection[reject_count].delta_x, rejection[reject_count].delta_y );
	  fstat = SCTPUT(aux_string);
	  sprintf(aux_string,"del_alpha and del_dec in arc_min are: %f, %f \n", del_alpha, del_dec);
	  fstat = SCTPUT(aux_string);

	  /* set back top stack position in sum-structure */
	  sum.total_opened--;
	  reject_count++;

	  /* return to calling function */
	  return 1;	/* signal to calling function that frame was rejected */
	}
    }



  /* get pixel offsets from observing position */
  /* get offsets in arcsecs (with sign) */
  del_alpha = (sum.info[last_pos].alpha - sum.reference.alpha) * 3600 * 15;	/* for RA in hours */
  del_dec = (sum.info[last_pos].declination - sum.reference.declination) * 3600;

  /* store offsets in sum.info.delta_x , _y */
  /* pixel_offset = arcsec_offset / plate_scale (arcsec/pix) */
  /* correct del_RA with cos(declination) */
  sum.info[last_pos].delta_x = del_alpha * cos(PI * sum.info[last_pos].declination / 180) / PIX_SCALE;
  sum.info[last_pos].delta_y = del_dec / PIX_SCALE;

  /* save these offsets as original offsets */
  sum.info[last_pos].org_delta_x = sum.info[last_pos].delta_x;
  sum.info[last_pos].org_delta_y = sum.info[last_pos].delta_y;


  /* do telescope drift correction by using the measured drift from last frame */
  if(last_pos >= 3)
    {
      sum.info[last_pos].delta_x +=  sum.info[last_pos-1].delta_x -  sum.info[last_pos-1].org_delta_x;
      sum.info[last_pos].delta_y +=  sum.info[last_pos-1].delta_y -  sum.info[last_pos-1].org_delta_y;
    }


  /* get normalization factor */
  norm_factor = NORM_LEVEL / sum.info[last_pos].countlevel;


  sprintf(aux_string,"\nloading image to clean stack at position: %d",clean.top_stack);
  fstat = SCTPUT(aux_string);


  /* load image to stack and do normalization */ 
  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
    {
      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	{
	  cosmics_stack[clean.top_stack][loop_y][loop_x] = (*p_frame) * norm_factor;
	  p_frame++;

	}
    }

  /* bookkeeping */
  clean.info[clean.top_stack] = sum.info[last_pos];


  /* incrementation modulo n_average */
  clean.top_stack = (clean.top_stack+1)%n_average;


  return 0;
}
