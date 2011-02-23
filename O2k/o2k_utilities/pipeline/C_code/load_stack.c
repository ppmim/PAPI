/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		load_stack.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 5 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		sky cube, image stack
.PURPOSE	     	load an input frame onto the image_stack (from top)
.IN/OUTPUT		IN: pointer to input image
			IN: image number of input image
			IN: pointer to image name
.RETURNS		Status: 0 = ok
.COMMENTS		The global variable 'top_stack' is checked and accordingly...
			If the stack is not full yet, the input image is put on top.
			If stack is full, the oldest image is copied to the array 'old_image'
			and its image number to 'old_number'. Then the new input image is
			copied into the stack in place of the old image.
			The global variables 'stack_book[]','master_pos','new_pos' and
			'old_pos' are updated.
.CALL			int Load(float *p_input, int image_no, char frame_name[], int sky_mode)
.VERSION		1.00		29.08.02
			1.1		29.10.02	corrected error concerning 'old_number' 
			2.0		14.11.02	add sky_mode parameter and check wether old 
							image is needed (not needed for clipping_mode)
		     	2.1		04.03.03	comment out printf()'s
-----------------------------------------------------------------------------------------*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header File for Structure Definitions */
#include "struct_def.h"
		
/* Header Files of Secundary Modules */
#include "load_stack.h"

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



int Load(float *p_input, int image_no, char frame_name[], int sky_mode)
{
  /* Make global variables visible to module */
  extern float image_stack[PIX_AXIS][PIX_AXIS][TOT_IMAGES_MAX];		
  extern float old_image[PIX_AXIS][PIX_AXIS];
  extern BOOKPAGE stack_book[TOT_IMAGES_MAX];		      	
  extern int old_number;	
  extern int top_stack;				
  extern int master_pos;			      
  extern int new_pos;
  extern int old_pos;				

  extern int TOT_IMAGES;				   	/* actual # of images for stack */


  /* Local Variables */
  int fstat;								/* function status */
  int test = 0;
  int loop_x = 0;							/* loop variable */
  int loop_y = 0;


  /* info on screen */
  fstat = SCTPUT("Loading image stack now...");


  /* Case 1: stack not filled up yet */
  if(top_stack < TOT_IMAGES)					/* test wether unfilled levels exist */
    {
      /* do stack_bookkeeping */
      stack_book[top_stack].frame_imno = image_no;	 	/* write image # into stack book */
      stack_book[top_stack].p_frame = p_input;			/* pointer to image */
      strcpy(stack_book[top_stack].frame_name, frame_name);	/* copy image name to stack_book */

      /* load image to stack */
      for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	{
	  for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	    {
	      image_stack[loop_y][loop_x][top_stack] = *p_input;
	      p_input++;
	      test++;
	    }
	}

/* printf("\n%d pixel values were written into unfilled stack level %d\n",\
 test, top_stack); */


      if(test == TOT_PIX)
	fstat = SCTPUT("Image was properly loaded onto stack");
      else
	{
	fstat = SCTPUT("Error: Image was not copied correctly to stack");
	return 14;						/* 14 = file handling error */
	}

      /* do final bookkeeping */
      new_pos = top_stack;					/* last image on stack */
      top_stack++;						/* next open stack level */

      return 0;							/* done, back to main() */
    }



  /* Case 2: stack is full */
  else if(top_stack == TOT_IMAGES)
    {
      /* do stack_bookkeeping */
      old_number = stack_book[old_pos].frame_imno;	   	/* copy old image number */
      stack_book[old_pos].frame_imno = image_no;	 	/* write new image # into stack book */
      stack_book[old_pos].p_frame = p_input;			/* pointer to image */
      strcpy(stack_book[old_pos].frame_name, frame_name);	/* copy image name to stack_book */


      /* write old frame if needed (for sky_mode = 0 or 1) */
      if(sky_mode==0 || sky_mode==1)	
	{
	  /* oldest image is copied into 'old_image' with image number 'old_number' */
	  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	    {
	      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
		{
		  old_image[loop_y][loop_x] = image_stack[loop_y][loop_x][old_pos];
		}
	    }

	  /* printf("\nOld_image loaded with image number %d\n", old_number); */
	}



      /* new frame is copied to position of old frame */
      for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	{
	  for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	    {
	      image_stack[loop_y][loop_x][old_pos] = *p_input;
	      p_input++;

	      test++;
	    }
	}
      /* printf("\nNew frame was written into full stack at level %d\n", old_pos); */

      if(test == TOT_PIX)
	fstat = SCTPUT("Image was properly copied into stack");
      else
	{
	fstat = SCTPUT("Error: Image was not copied correctly into stack");
	return 14;						/* 14 = file handling error */
	}


      /* final bookkeeping */
      new_pos = old_pos;					/* position of newest image on stack */
      old_pos = (old_pos+1)%TOT_IMAGES;			   	/* incr. old_pos modulo 'tot_images' */
      master_pos = (master_pos+1)%TOT_IMAGES;		   	/* increment master_pos */

      return 0;							/* done, back to main() */
    }



  /* case 3: top_stack has undefined level */
  else
  {
    fstat = SCTPUT("Error: Stack level overflow");
    /* printf("\nThe stack level is: %d ; should be between 0 and %d\n",\
       top_stack, TOT_IMAGES); */

    return 7;							/* 7 = input invalid */
  }

}
