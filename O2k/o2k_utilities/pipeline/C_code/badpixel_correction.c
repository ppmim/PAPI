/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		badpixel_correction.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 21 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		bad pixel correction
.PURPOSE	     	replace bad pixel by average of surrounding good pixels
.IN/OUTPUT		IN: pointer to frame to be BP-corrected
.RETURNS		Status: 0 = ok
.COMMENTS		Bad pixels are searched by looking for 1's in the badpixel mask.
			Once a bad pixel is found, the neighboring good pixels are searched.
			The bad pixel is then replaced by the average of the linearly
			scaled (good) pixel values in the horizontal and vertical direction. 
			The procedure is done in the image_stack at position new_pos.
.CALL			int Badpixel()
.VERSION		1.00   	04.03.03
			2.0		28.03.03	added debugging mode
			3.0		17.04.03 	corrected whole BPM process; now use bpm_aux
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "badpixel_correction.h"

/* Parameter and Symbolic Constant Files */
#include "pipe_para.h"		

/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif

/* check wether TOT_IMAGES_MAX is defined */
#if !defined(TOT_IMAGES_MAX)
#define TOT_IMAGES_MAX (2*SKY_FRAMES_MAX+1)
#endif


/* Make global variables visible to module */
extern float badpixel_mask[PIX_AXIS][PIX_AXIS];
extern float bpm_aux[PIX_AXIS][PIX_AXIS];

			
 

int Badpixel(float *p_inframe)
{
  /* local variables */
  int fstat;						   	/* function status */
  int loop_x, loop_y;					     	/* loop parameters */
  int test = 0;						     	/* count corrected pixels */
  int go_left;					      		/* steps to next good pixel*/
  int go_right;
  int go_up;
  int go_down;

  float left_value;					    	/* value of good pixel */
  float right_value;
  float up_value;
  float down_value;

  float vert_value;						/* good value in vertical direction */
  float horiz_value;						/* in horizontal direction */
  float good_value;						/* value that is replaced */
  float *p_aux = p_inframe;					/* auxiliary pointer to inframe*/

  char aux_string[200];



  fstat = SCTPUT("\nDoing bad pixel correction ...");




  /* copy input frame to auxiliary array */
  for(loop_y = 0 ; loop_y < PIX_AXIS ; loop_y++)
    {
      for(loop_x = 0 ; loop_x < PIX_AXIS ; loop_x++)
	{
	  bpm_aux[loop_y][loop_x] = *p_inframe;		/* copy value */

	  p_inframe++;					/* increment pointer */
	}
    }


 

  /* do bad pixel correction in bpm_aux[][] */
  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
    {
      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	{

	  /* test wether pixel is bad --> value == 1 */
	  if(badpixel_mask[loop_y][loop_x] > 0.5)
	    {
	      /* pixel bad */
	      test++;

	      /* reset search variables */
	      go_left = 1;
	      go_right = 1;
	      go_up = 1;
	      go_down = 1;


	      /* look for next good pixel to left*/
	      while( (loop_x - go_left) >= 0  )
		{		
		  if( badpixel_mask[loop_y][loop_x - go_left] < 0.5)		/* good pixel */
		    break;

		  go_left++;
		}

	      if( (loop_x - go_left) < 0  )		/* check wether inside frame */
		go_left = 0;

	      /* save good pixel if available */
	      if(go_left > 0)
		left_value =  bpm_aux[loop_y][loop_x-go_left];



	      /* look for next good pixel to right*/
	      while( (loop_x + go_right) < PIX_AXIS)
		{
		  if( badpixel_mask[loop_y][loop_x + go_right] < 0.5)		/* good pixel */
		    break;

		  go_right++;
		}

	      if( (loop_x + go_right) >= PIX_AXIS  )	/* check wether inside frame */
		go_right = 0;

	      /* save good pixel if available */
	      if(go_right > 0)
		right_value =  bpm_aux[loop_y][loop_x+go_right];



	      /* look for next good pixel up*/
	      while((loop_y + go_up) < PIX_AXIS)
		{
		  if( badpixel_mask[loop_y+go_up][loop_x] < 0.5)		/* good pixel */
		    break;

		  go_up++;
		}

	      if( (loop_y + go_up) >= PIX_AXIS  )	/* check wether inside frame */
		go_up = 0;

	      /* save good pixel if available */
	      if(go_up > 0)
		up_value =  bpm_aux[loop_y + go_up][loop_x];



	      /* look for next good pixel down*/
	      while((loop_y - go_down) >= 0)
		{
		  if( badpixel_mask[loop_y-go_down][loop_x] < 0.5)		/* good pixel */
		    break;

		  go_down++;
		}

	      if( (loop_y - go_down) < 0 )		/* check wether inside frame */
		go_down = 0;

	      /* save good pixel if available */
	      if(go_down > 0)
		down_value =  bpm_aux[loop_y - go_down][loop_x];




	      /* get scaled horzontal value */
	      if(go_left == 0)			/* no left value available */
		horiz_value = right_value;

	      else if(go_right == 0)		/* no right value available */
		horiz_value = left_value;
				
	      else				/* both values available */
		{
		  horiz_value = left_value + (right_value - left_value)*go_left/(go_left+go_right);
		  /* linear interpolation between the two values */
		}


	      /* get scaled vertical value */
	      if(go_up == 0)			/* no up value available */
		vert_value = down_value;

	      else if(go_down == 0)		/* no down value available */
		vert_value = up_value;
				
	      else				/* both values available */
		{
		  vert_value = down_value + (up_value - down_value)*go_down/(go_down+go_up);
		  /* linear interpolation between the two values */
		}


	      /* get good value as average*/
	      good_value = (horiz_value+vert_value)/2;

	      /* replace bad pixel by good value */
	       bpm_aux[loop_y][loop_x] = good_value;


	      /* debugging mode */
	      if(DEBUG == 2)
		{
		  sprintf(aux_string,"..det. good values: left %6.1f, right: %6.1f up: %6.1f, down: %6.1f", \
		     left_value, right_value, up_value, down_value);
		  fstat = SCTPUT(aux_string);


		  sprintf(aux_string,"Bad pixel replaced by %f, calcutlated from hor: %f and vert: %f", \
		     good_value,horiz_value,vert_value);
		  fstat = SCTPUT(aux_string);
		}


	    }	/* end of if bad pixel */


	}
    }		/* end for loop */





 /* write corrected frame back to input frame */
  for(loop_y = 0 ; loop_y < PIX_AXIS ; loop_y++)
    {
      for(loop_x = 0 ; loop_x < PIX_AXIS ; loop_x++)
	{
	  *p_aux = bpm_aux[loop_y][loop_x];		/* copy value */

	  p_aux++;					/* increment pointer */
	}
    }






  if(DEBUG == 1)
    {
      sprintf(aux_string,"...%d bad pixels corrected", test);
      fstat = SCTPUT(aux_string);
    }


  return 0;
}
