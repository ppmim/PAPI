/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		copy_frame.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 9 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS	      	copy frames
.PURPOSE		copy an input frame to an output frame
.IN/OUTPUT		IN: pointer to input frame 
		     	IN: pointer to the 	output frame
			IN: total number of pixels of frames	
			IN: image number of output frame (for updating descriptors)
			IN: pointer to name of input frame (for updating descriptor HISTORY)
.RETURNS		Status: 0 = ok
.COMMENTS		This function is similar to the function Flat() of Module 1 and is
			used to setup the output frame in case no flatfield correction is
			done in the beginning of the pipeline. The contents of the input frame
			is copied pixelwise to the output frame.
.CALL			int Copy(float *p_inframe, float *p_outframe, int total_pix, 
				   int image_no_out, char *name_inframe)
.VERSION		1.00		05.09.02
			1.01		15.11.02	minor changes on HISTORY
--------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<midas_def.h>		/* MIDAS definitions */


/* Header Files of Secundary Modules */
#include "copy_frame.h"


int Copy(float *p_inframe, float *p_outframe, int total_pix, int image_no_out, char *name_inframe)
{

  /* local variables */
  int count;
  int fstat; 		/* test Function STATus */
  int unit = 0;		/* dummy variable for interface functions */

  time_t time_secs;	/* time variable */

  char *p_timestring;	/* pointer for local time */


  /* display message */
  fstat =  SCTPUT("\n\nCopying frame now...");


  /* loop that does the copying */
  for (count=1 ; count <= total_pix ; count++)
    {
     
      *p_outframe = *p_inframe;

      /* pointer incrementation */
      p_inframe++;     
      p_outframe++;
    }



  /* update descriptors HISTORY of outframe ; create it first if necessary */
  /* get time and date of reduction */
  time(&time_secs);
  p_timestring = ctime(&time_secs);

  /* write time into descriptor */
  fstat = SCDWRC(image_no_out,"HISTORY",1,p_timestring,-1,26,&unit);

  /* original frame */
  fstat = SCDWRC(image_no_out,"HISTORY",1,"Frame was copied from original image:",-1,80,&unit);
  fstat = SCDWRC(image_no_out,"HISTORY",1,name_inframe,-1,80,&unit);


  fstat = SCTPUT("...copying done\n\n");

  return 0;

}

