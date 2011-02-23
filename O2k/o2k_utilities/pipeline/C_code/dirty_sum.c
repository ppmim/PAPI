/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		dirty_sum.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 16 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		real sum 
.PURPOSE	     	sum up pictures in world coordinates without cosmics removal
.IN/OUTPUT		IN: p_inframe = pointer to input image
.RETURNS		Status: 0 = ok
.COMMENTS		Take input frame and get relevant information on frame from top_level
			of sum.info stack: x_offset, y_offset ,name etc.
			Get integer pixel offsets.
			Add input frame to sum.dirty_sum output frame. Where the input frame and
			the reference frame do not overlap, set value to 0.	Thus cosmics are 
			not removed.
.CALL			int Dirty(float *p_inframe)
.VERSION		1.00   		27.11.02
			1.1			02.12.02 	changed setting to get 0 for no overlap
			1.2			03.12.02	corrected for inverse RA orientation
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<string.h>
#include<math.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "dirty_sum.h"

/* Header file with structure definitions */
#include "struct_def.h"

/* Parameter and Symbolic Constant Files */
#include "pipe_para.h"		

/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif


 /* make gloabel variables accessible */
extern DEEP_STRUCT sum;



int Dirty(float *p_inframe)
{
  /* local variables */
  double next_int;							/* buffer for offsets */

  int fstat;								/* function status */
  int last_pos = sum.total_opened -1;					/* position of frame info */
  int del_x;								/* pixel offsets */
  int del_y;
  int loop_x, loop_y;							/* loop parameters */

 

  fstat = SCTPUT("\nMaking dirty sum...");


  /* get integer pixels offsets for frame */
  /* find closest integer */
  if(sum.total_opened > 1)	/* if not reference frame */
    {
      /* x_offset */ 
      next_int = ceil(sum.info[last_pos].delta_x);		/* next larger integer */

      if(fabs(sum.info[last_pos].delta_x - next_int) < 0.5)	/* if next_int is closest */
	del_x = next_int;
      else							/* else take the smaller one */
	del_x = next_int-1;

      /* y_offset */
      next_int = ceil(sum.info[last_pos].delta_y);		/* next larger integer */

      if(fabs(sum.info[last_pos].delta_y - next_int) < 0.5)	/* if next_int is closest */
	del_y = next_int;
      else							/* else take the smaller one */
	del_y = next_int-1;

    }


  /* Test */
  /*  printf("\n...using offsets in x and y: %d , %d", del_x, del_y); */


  /* if this is the reference frame, do straight copying */
  if(sum.total_opened <= 1)		/* reference frame */
    {
 
      for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	{
	  for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	    {
	      sum.dirty_sum[loop_y][loop_x] = *p_inframe;

	      /* pointer incrementation */
	      p_inframe++;     	      
	    }
	}

      return 0;		/* back to calling function */
    }



  /* if this is a new frame, do summation in world coordinates */
  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
    {
      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	{
	  /* check wether world coordinate pixel is available */
	  /* in x-direction */
	  if((loop_x + del_x) < 0  || (loop_x + del_x) >= PIX_AXIS )
	    {
	      sum.dirty_sum[loop_y][loop_x] = -1;	/* no value available */
	      continue;		/* go to next pixel */
	    }

	  /* in y-direction */
	  if((loop_y - del_y) < 0  || (loop_y - del_y) >= PIX_AXIS )
	    {
	      sum.dirty_sum[loop_y][loop_x] = -1;	/* no value available */
	      continue;		/* go to next pixel */
	    }

	  
	  /* add value to existing one, if it is within overlap of all frames */
	  if(sum.dirty_sum[loop_y][loop_x] > -0.5)
	    sum.dirty_sum[loop_y][loop_x] += *(p_inframe + (PIX_AXIS * (loop_y - del_y)) + (loop_x + del_x));

	}
    }


  return 0;
}
