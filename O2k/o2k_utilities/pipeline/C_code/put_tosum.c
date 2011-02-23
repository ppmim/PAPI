/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		put_tosum.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 17 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		cleaned sum 
.PURPOSE	     	sum up pictures in world coordinates without with removal
.IN/OUTPUT		IN: sum_save_flag = flag for saving summed frames
			IN: n_average = number of images on clean.cosmics_stack
			IN: n_real = N_AVERAGE ; number of images for median
.RETURNS		Status: 0 = ok
.COMMENTS		The functions takes the current frame in sum.scratch_sum and normalizes
			the frame values to the actual count level, which is the sum of the 
			countlevels of all frames on clean.cosmics_stack.
			This renormalized frame is then added to the frame sum.master_sum.
			If intermediate results should be saved (SUM_SAVE_FLAG = 0), the current
		    	sum frame is copied to an output frame, HISTORY updated. The offsets
			determined from the telescope position are written in X/Y_ORG_OFFSETS,
			the real offsets from the object finding in descriptor X_Y_REAL_OFFSETS
.CALL			int Tosum(int sum_save_flag, int n_average, int n_real)
.VERSION		1.00   		27.11.02
			1.01			28.11.02	implement use of sum.total_sum
			1.1			02.12.02	change to make value 0 if no overlap
								and take out frame closing
			1.2			03.12.02	corrected summation
			1.3			06.03.03	changed for case if n_average=0
			1.4			06.03.03	changed all output to SCTPUT
			1.5			26.03.03	write offsets in HISTORY in descriptors:
								X_ORG_OFFSETS, X_REAL_OFFSETS
			2.0			05.06.03	add n_real to get HISTORY right
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<string.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "put_tosum.h"

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
extern COSMICS_STACK clean;



int Tosum(int sum_save_flag, int n_average, int n_real)
{				      
  /* local variables */
  int fstat;								/* function status */
  int last_pos = sum.total_opened -1;					/* position of frame info */
  int loop_x, loop_y;							/* loop parameters */
  int count;
  int test = 0;
  int unit = 0;								/* dummy */

  float real_factor;							/* factor to get real level */
  float real_level = 0;							/* real countlevel */
  float *p_outframe;							/* pointer to output image */

  char aux_string[200];							/* string buffer */


  fstat = SCTPUT("Put_to_sum now...");


  /* initialize output pointer */
  p_outframe = sum.master.p_frame;					/* pointer to mapped output */


  /* do only if new frames are on stack --> n_average>0 */
  if(n_average >= 1)
    {
      /* get real count physical level */
      for(count=0 ; count<n_average ; count++ )
	real_level += clean.info[count].countlevel;
    
      /* save this in sum.sum_level */
      sum.sum_level += real_level;

      /* get factor for renormalization */
      /* median level was NORM_LEVEL --> real_factor = real_level / NORM_LEVEL */
      real_factor = real_level / NORM_LEVEL;

      sprintf(aux_string,"The factor for setting countlevel right is: %f", real_factor);
      fstat = SCTPUT(aux_string);


      /* do renormalization and add to sum.master_sum*/
      /* if this is first summation, do copying : sum.total_sum == 0*/
      /* if not, check for overlap */ 
      if(sum.total_sum == 0)
	{
	  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	    {
	      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
		{
		  /* do sum for the first time */
		  if(sum.scratch_sum[loop_y][loop_x] > -0.5)
		    sum.scratch_sum[loop_y][loop_x] *= real_factor;

		  /* add it */
		  sum.master_sum[loop_y][loop_x] = sum.scratch_sum[loop_y][loop_x];
		}
	    }
	}


      /* new frames */
      else
	{ 
	  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	    {
	      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
		{
		  /* do sum only if there is overlap for all frames */
		  if(sum.master_sum[loop_y][loop_x] > -0.5)
		    {
		  
		      if( sum.scratch_sum[loop_y][loop_x] > -0.5)
			{
			  sum.scratch_sum[loop_y][loop_x] *= real_factor;
			  
			  /* add it */
			  sum.master_sum[loop_y][loop_x] += sum.scratch_sum[loop_y][loop_x];
			}

		      
		      /* if no overlap with new frame, set to -1 */
		      else
			sum.master_sum[loop_y][loop_x] = -1;

		    }

		}
	    }
	}


      /* update sum.total_sum */
      sum.total_sum = sum.total_opened;

    }




  /* if sum_save_flag = 0, write sum.master_sum in output frame */    
  if(sum_save_flag == 0)
    {

      for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	{
	  for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	    {

	      if(sum.master_sum[loop_y][loop_x] >= 0)
		{	
		  *p_outframe = sum.master_sum[loop_y][loop_x];
		  test++;
		}

	      else
		*p_outframe = 0;


	      p_outframe++;
	    }
	}


      /* update HISTORY */
      fstat = SCDWRC(sum.master.ima_number,"HISTORY",1,"This image is the cosmics cleaned sum of:",-1,80,&unit);
  
      sprintf(aux_string,"The reference median-image was determined from %d images.", n_real);
      fstat = SCDWRC(sum.master.ima_number,"HISTORY",1,aux_string,-1,80,&unit);


      for(count=0 ; count<= last_pos ; count++)
	{
	  fstat = SCDWRC(sum.master.ima_number,"HISTORY",1,sum.info[count].ima_name,-1,80,&unit);

	  sprintf(aux_string,"pixel offsets for the last frame were: %4.1f in x ; %4.1f in y",\
		  sum.info[count].delta_x, sum.info[count].delta_y);
	  fstat = SCDWRC(sum.master.ima_number,"HISTORY",1,aux_string,-1,80,&unit);

	  /* write descriptors for pixel offsets */
	  /* origianl offsets */
	  fstat = SCDWRR(sum.master.ima_number,"X_ORG_OFFSETS",&sum.info[count].org_delta_x,-1,1,&unit);
	  fstat = SCDWRR(sum.master.ima_number,"Y_ORG_OFFSETS",&sum.info[count].org_delta_y,-1,1,&unit);

	  /* determined offsets */
	  fstat = SCDWRR(sum.master.ima_number,"X_REAL_OFFSETS",&sum.info[count].delta_x,-1,1,&unit);
	  fstat = SCDWRR(sum.master.ima_number,"Y_REAL_OFFSETS",&sum.info[count].delta_y,-1,1,&unit);
	}


      /* write help text for new descriptors */
      fstat = SCDWRH(sum.master.ima_number,"X_ORG_OFFSETS","x_offset determined from telescope position",-1,80);
      fstat = SCDWRH(sum.master.ima_number,"Y_ORG_OFFSETS","y_offset determined from telescope position",-1,80);
      fstat = SCDWRH(sum.master.ima_number,"X_REAL_OFFSETS","x_offset determined from found objects",-1,80);
      fstat = SCDWRH(sum.master.ima_number,"Y_REAL_OFFSETS","y_offset determined from found objects",-1,80);




      sprintf(aux_string,"\n%d values were written in sum_output_frame...\n", test);
      fstat = SCTPUT(aux_string);
    }



  return 0;
}
