/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		frame_info.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 13 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		image information 
.PURPOSE	     	get relevant information on input image for summation process
.IN/OUTPUT		IN: name of input image
			IN: image number of input image
.RETURNS		Status: 0 = ok
.COMMENTS		The input frame must have a current Statistics() update.
			Frame name, image number, observing positons, and median countlevel
			are then stored at the next open position of the info_stack inside the 
			sum structure. If this is the first input image, the information is also
			stored in the 'reference' info page
.CALL			int Info(char *p_name, int image_no)
.VERSION		1.00   		22.11.02
			1.1			02.12.02	take out icrementation of total_opened
			1.2			03.12.02	update O_POS for output frame
			1.3			06.03.03	changed all output to SCTPUT
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<string.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "frame_info.h"

/* Header file with structure definitions */
#include "struct_def.h"

/* make gloabel variables accessible */
extern DEEP_STRUCT sum;	



int Info(char *p_name, int image_no)
{

  /* local variables */
  int fstat;								/* function status */
  int actvals;								/* return from read descriptor */
  int unit = 0;								/* dummies */
  int null = 0;

  char aux_string[200];							/* auxilary string buffer */


  fstat = SCTPUT("\nGet frame information for summation...");


  /* copy name and ima_number to open info-stack-level */
  strncpy(sum.info[sum.total_opened].ima_name, p_name, 84);
  sum.info[sum.total_opened].ima_number = image_no;


  /* get median countlevel from SUB_STATISTICS */
  fstat =  SCDRDR(image_no,"SUB_STATISTICS",3,1,&actvals,&sum.info[sum.total_opened].countlevel,&unit,&null );


  /* get observing positions */
  /* RA */
  fstat = SCDRDD(image_no,"O_POS",1,1,&actvals,&sum.info[sum.total_opened].alpha,&unit,&null);	  
  /* DEC */
  fstat = SCDRDD(image_no,"O_POS",2,1,&actvals,&sum.info[sum.total_opened].declination,&unit,&null); 

  sprintf(aux_string,"new observing positions are: alpha = %f , dec = %f",sum.info[sum.total_opened].alpha,\
	  sum.info[sum.total_opened].declination);
  fstat = SCTPUT(aux_string);


  /* check wether this is first image */
  if(sum.total_opened == 0)		/* first opened image */
    {
      /* write info into 'reference' */
      sum.reference = sum.info[0];

      /* copy descriptor O_POS from reference to outout frame */
      fstat = SCDCOP(sum.reference.ima_number,sum.master.ima_number,4,"O_POS");

      /* give out reference frame */
      fstat = SCTPUT("\nThe reference frame for the summation process is:");
      fstat = SCTPUT(sum.reference.ima_name);
      sprintf(aux_string,"image number: %d", sum.reference.ima_number);
      fstat = SCTPUT(aux_string);
      sprintf(aux_string,"median countlevel: %f", sum.reference.countlevel);      
      fstat = SCTPUT(aux_string);     
      sprintf(aux_string,"observing position: RA = %f", sum.reference.alpha);
      fstat = SCTPUT(aux_string);      
      sprintf(aux_string,"observing position: DEC = %f\n",\
	      sum.reference.declination);
      fstat = SCTPUT(aux_string);
    }



  return 0;
}
