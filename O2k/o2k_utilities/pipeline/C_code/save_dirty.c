/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		save_dirty.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 18 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		save dirty_sum 
.PURPOSE	     	save the frames dirty_sum and difference and update HISTORY
.IN/OUTPUT		none
.RETURNS		Status: 0 = ok
.COMMENTS		The function takes the frame sum.dirty_sum and writes it to the output
			frame specified in sum.dirty. The difference of sum.master_sum and
			sum.dirty_sum is written in the frame specified in sum.difference.
			HISTORY is updated for the output frames and frames are closed.
.CALL			int Savedirty()
.VERSION		1.00   		28.11.02
			1.01			29.11.02	corrected diff_ima_number for HISTORY
			1.1			02.12.02	changed to get 0 for no overlap
			1.2			04.12.02	switched difference around
			1.3			04.12.02	update LHCUTS
			1.31			11.02.03	only conditional include of pipe_para.h
			1.4			26.03.03	write pixel offsets in HISTORY
			1.5			07.04.03	write descriptores X_ORG_OFFSETS and
								X_REAL_OFFSETS
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "save_dirty.h"
#include "subframe_statistics.h"

/* Header file with structure definitions */
#include "struct_def.h"

/* Parameter and Symbolic Constant Files */
#if !defined( DETECTOR )				/* include only if necessary */
#include "pipe_para.h"		
#endif

/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif


/* make gloabel variables accessible */
extern DEEP_STRUCT sum;



int Savedirty()
{
  /* local variables */
  int fstat;								/* function status */
  int last_pos = sum.total_opened -1;					/* position of frame info */
  int loop_x, loop_y;							/* loop parameters */
  int count;
  int unit = 0;								/* dummy */

  float *p_dirty;							/* pointer to output image */
  float *p_difference;

  char aux_string[200];							/* buffer for SCTPUT */



  fstat = SCTPUT("Savin dirty_sum now...");


  /* initialize output pointers */
  p_dirty = sum.dirty.p_frame;					/* pointer to mapped output for dirty */
  p_difference = sum.difference.p_frame;			/* pointer to difference output */



  /* write output image */
  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
    {
      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	{
	  if(sum.dirty_sum[loop_y][loop_x] < -0.5)
	    *p_dirty = 0;

	  else
	    *p_dirty = sum.dirty_sum[loop_y][loop_x];

	  /* write difference_frame = dirty-master */
	  *p_difference = sum.dirty_sum[loop_y][loop_x] - sum.master_sum[loop_y][loop_x];

	  p_dirty++;
	  p_difference++;
	}
    }


  /* update LHCUTS for frames */
  fstat = Statistics(sum.dirty.p_frame, sum.dirty.ima_number,3,10);  
  fstat = Statistics(sum.difference.p_frame, sum.difference.ima_number,3,5); 


  /* update HISTORY for dirty_sum*/
  fstat = SCDWRC(sum.dirty.ima_number,"HISTORY",1,"This image is the real sum of:",-1,80,&unit);

  for(count=0 ; count<= last_pos ; count++)
    {
      fstat = SCDWRC(sum.dirty.ima_number,"HISTORY",1,sum.info[count].ima_name,-1,80,&unit);

      sprintf(aux_string,"pixel offsets for the last frame were: %4.1f in x ; %4.1f in y",\
	      sum.info[count].delta_x, sum.info[count].delta_y);
      fstat = SCDWRC(sum.dirty.ima_number,"HISTORY",1,aux_string,-1,80,&unit);


      /* write descriptors for pixel offsets */
      /* origianl offsets */
      fstat = SCDWRR(sum.dirty.ima_number,"X_ORG_OFFSETS",&sum.info[count].org_delta_x,-1,1,&unit);
      fstat = SCDWRR(sum.dirty.ima_number,"Y_ORG_OFFSETS",&sum.info[count].org_delta_y,-1,1,&unit);

      /* determined offsets */
      fstat = SCDWRR(sum.dirty.ima_number,"X_REAL_OFFSETS",&sum.info[count].delta_x,-1,1,&unit);
      fstat = SCDWRR(sum.dirty.ima_number,"Y_REAL_OFFSETS",&sum.info[count].delta_y,-1,1,&unit);
    }


  /* write help text for new descriptors */
  fstat = SCDWRH(sum.dirty.ima_number,"X_ORG_OFFSETS","x_offset determined from telescope position",-1,80);
  fstat = SCDWRH(sum.dirty.ima_number,"Y_ORG_OFFSETS","y_offset determined from telescope position",-1,80);
  fstat = SCDWRH(sum.dirty.ima_number,"X_REAL_OFFSETS","x_offset determined from found objects",-1,80);
  fstat = SCDWRH(sum.dirty.ima_number,"Y_REAL_OFFSETS","y_offset determined from found objects",-1,80);




  /* close this frame */
  fstat = SCTPUT("\nWriting out dirty_sum frame now:");
  fstat = SCTPUT(sum.dirty.name);

  fstat = SCFCLO(sum.dirty.ima_number);
    
 
  /* update HISTORY for difference_sum*/
  fstat = SCDWRC(sum.difference.ima_number,"HISTORY",1,\
		 "This image is difference of the two sum images:",-1,80,&unit);
  fstat = SCDWRC(sum.difference.ima_number,"HISTORY",1,sum.dirty.name,-1,80,&unit);
  fstat = SCDWRC(sum.difference.ima_number,"HISTORY",1,sum.master.name,-1,80,&unit);
 
  /* close this frame */
  fstat = SCTPUT("\nWriting out difference_sum frame now:");
  fstat = SCTPUT(sum.difference.name);

  fstat = SCFCLO(sum.difference.ima_number);


  return 0;
}
