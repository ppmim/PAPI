/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		save_sky.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 7 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		sky frame
.PURPOSE	     	save sky frame
.IN/OUTPUT		IN: image number of sky frame
			IN: pointer to sky frame
.RETURNS		Status: 0 = ok
.COMMENTS		The functions copies the current contents of the array sky_frame to
			the destinated sky frame and unmaps and closes this frame afterwards.
.CALL			int Savesky(int sky_imno, float *p_skyframe)
.VERSION		1.00		03.09.02
			1.1		15.11.02	update HISTORY
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<midas_def.h>							/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "save_sky.h"

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




int Savesky(int sky_imno, float *p_skyframe)
{
  /* Make global variables visible to module */
  extern BOOKPAGE stack_book[TOT_IMAGES_MAX];		      	
  extern float sky_frame[PIX_AXIS][PIX_AXIS];
  extern int TOT_IMAGES;

 	      	
  /* Local Variables */			
  int fstat;								/* function status */
  int count = 0;
  int loop = 0;								/* loop variable */
  int unit = 0;								/* dummy */

  float *p_sky_array = &sky_frame[0][0];				/* pointer variable to array */


  /* info on screen */
  fstat = SCTPUT("Saving sky now...\n");


  /* Copy array sky_frame to destinated image frame */
  for(loop=1 ; loop<=TOT_PIX ; loop++)
    {
      *p_skyframe = *p_sky_array;

      p_skyframe++;
      p_sky_array++;
    }


  /* write descriptor HISTORY */
  /* write all images for sky determination in HISTORY */
  fstat = SCDWRC(sky_imno,"HISTORY",1,\
		 "The sky was determined from the following images:",-1,80,&unit);

  for(count=0 ; count < TOT_IMAGES ; count++)
    {
      fstat = SCDWRC(sky_imno,"HISTORY",1,\
		     stack_book[count].frame_name,-1,80,&unit);
    }



  /* unmap and close saved sky frame */
  fstat = SCFCLO(sky_imno);


  return 0;
}

