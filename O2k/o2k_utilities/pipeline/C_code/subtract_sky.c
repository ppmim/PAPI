/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		subtract_sky.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 8 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		subtract sky
.PURPOSE	     	sutract modelled sky from the master_image and add a constant that
			corresponds to the average sky level  
.IN/OUTPUT		IN: SKY_MODE (everthing for descriptor HISTORY)
			IN: N_SMALLEST
			IN: KAPPA
.RETURNS		Status: 0 = ok
.COMMENTS		The mean sky value is determined in the subarray specified in the 
			detector parameter file. The following operations are then done on 
			the frame at position 'master_pos' in 'stack_book':
			1) the sky value is subtracted pixelwise, taken from 'sky_frame'
			2) the mean sky value is added (to maintain the correct statistics) 
			The descriptor HISTORY is updated.
.CALL			int Subtract(int sky_mode, int n_smallest, float kappa)
.VERSION		1.00		04.09.02
			2.0		14.11.02 	get sky_mode parameters and update HISTORY
			2.1		06.03.03	changed all output to SCTPUT
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "subtract_sky.h"

/* Header File for Structure Definitions */
#include "struct_def.h"

/* Pipeline Parameters, including DETECTOR */
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




int Subtract(int sky_mode, int n_smallest, float kappa)
{
  /* Make global variables visible to module */
  extern BOOKPAGE stack_book[TOT_IMAGES_MAX];		      	
  extern float sky_frame[PIX_AXIS][PIX_AXIS];
  extern int master_pos;			      
  extern int TOT_IMAGES;
	

  /* Local Variables */
  int fstat;								/* function status */
  int count = 0;
  int unit = 0;								/* dummy */
  int loop_x = 0;							/* loop variable */
  int loop_y = 0;
  int i,j;

  float *p_frame;							/* pointer variable for frame */
  float sum = 0;							/* variables for statistics */
  float mean;								/* mean sky level */
  float min;
  float max = 0;

  int xpix = X_END - X_START;				    	/* length of subarray in x */
  int ypix = Y_END -Y_START;  				     	/* length of subarray in y */
  int tot_pix = (xpix+1)*(ypix+1);				/* total number of pixels in subarray */

  char string_buffer[80];					/* buffer for HISTORY update */
  char aux_string[200];						/* auxilary string buffer for SCTPUT */


  /* info on screen */
  fstat = SCTPUT("Subtracting sky now...");


  /* get minimum, maximum and average = mean of sky values */
  min = sky_frame[Y_START][X_START];		    	/* initialize min with first array value */

  for(i=Y_START ; i <= Y_END ; i++)
    {
    for(j=X_START ; j <= X_END ; j++)  
      {
	sum += sky_frame[i][j];

	if( sky_frame[i][j] < min)
	  min =  sky_frame[i][j];

	if( sky_frame[i][j] > max)
	  max =  sky_frame[i][j];
      }
    }

  mean = sum / tot_pix;

  sprintf(aux_string,"mean sky value = %f, max = %f, min = %f\n", mean, max, min);
  fstat = SCTPUT(aux_string);


  /* intitialize pointer variable with adress of frame to be manipulated */
  p_frame = stack_book[master_pos].p_frame;

  /* do the sky subtraction */
  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
    {
      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	{
	  *p_frame = *p_frame - sky_frame[loop_y][loop_x] + mean;	 

	  p_frame++;
	}
    }

 

  /* update descriptor HISTORY */
  fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,"Sky was subtracted from image.",-1,80,&unit);

  if(sky_mode == 1)		/* median mode */
    {     
      fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,\
		     "Sky was determined and subtracted in pixel-median-mode",-1,80,&unit);
    }

  else if(sky_mode == 0) 	/* minimum mode */
    {
      fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,\
		     "Sky was determined and subtracted in minimum-average-mode",-1,80,&unit);

      fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,\
		     "Pixel-sky-value = average of n_smallest values",-1,80,&unit);

      /* convert input parameters in string */
      sprintf(string_buffer,"The value for n_smallest was set to: %d", n_smallest);

      fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,string_buffer,-1,80,&unit);
    }

  else if(sky_mode == 2) 	/* clipping mode */
    {
      fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,\
		     "Sky was determined and subtracted in clip-outlier-mode",-1,80,&unit);

      fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,\
		     "Pixel-sky-value = median of values without outliers",-1,80,&unit);

      /* convert input parameters in string */
      sprintf(string_buffer,"The value for the clipping parameter kappa was set to: %f", kappa);

      fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,string_buffer,-1,80,&unit);
    }


  /* write added constant in HISTORY */
  sprintf(string_buffer,"A constant sky value was added to each pixel:: %f", mean);
  fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,string_buffer,-1,80,&unit);



  /* write all images for sky determination in HISTORY */
  fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,\
		 "The sky was determined from the following images:",-1,80,&unit);


  for(count=0 ; count < TOT_IMAGES ; count++)
    {
      fstat = SCDWRC(stack_book[master_pos].frame_imno,"HISTORY",1,\
		     stack_book[count].frame_name,-1,80,&unit);
    }



  return 0;
}
