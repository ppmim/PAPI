/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		flat_correction.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 1 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		flatfield correction
.PURPOSE		flatfield correction by dividing input frame by normalized flatfield 
.IN/OUTPUT		IN: pointer to input frame to be flatfield corrected
			IN: pointer to flatfield frame
			IN: pointer to the flatfield-corrected	output frame
			IN: total number of pixels of frames	
			IN: image number of output frame (for updating descriptors)
			IN: pointer to name of input frame (for updatin descriptor HISTORY)
			IN: pointer to name of flatfield frame
			IN: dark_flag=0-->no dark current subtraction; =1-->do dark subtraction
			IN: pointer to dark current frame
			IN: pointer to name of dark frame
			IN: integration time of dark frame for proper scaling
			IN: image number of inframe
.RETURNS		Status: 0 = ok
.COMMENTS		The function assumes that all relevant frames have been opened
			inside the main module. The pointers to the first pixel of each frame
			is then passed on to the function, where pixel by pixel division is done
			over all pixels.
			A decent flatfield frame has to be used, e.g. there should be no values
 			close to 0, otherwise the "corrected values" -> infinity! 
			The quality of the flatfield is not checked (only that values are 
			> 0.01). If dark_flag is set to 1, the scaled dark current is 
			subtracted.
.CALL			int Flat(float *p_inframe, float *p_flat, float *p_outframe, 
			  int total_pix, int image_no_out, char *name_inframe, char *flat_frame,
		        int dark_flag, float *p_dark, char *dark_name,
			  double dark_itime, int inframe_imno)
.VERSION		1.00		13.08.02
			1.1		06.09.02 	Update function, so it can handle the same
			      				input and output frame
			1.2		14.11.02	Update descriptor HISTORY, with time and 
							flatfield name
			1.3		06.03.03	changed output to SCTPUT
			1.4		26.03.03	add check that flatfield value > 0.01
			2.0		28.03.03	added dark current correction at same time
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Pipeline Parameters, including DETECTOR */
#include "pipe_para.h"

/* Header Files of Secundary Modules */
#include "flat_correction.h"


int Flat(float *p_inframe, float *p_flat, float *p_outframe, int total_pix, int image_no_out,\
	 char *name_inframe, char *flat_frame, int dark_flag, float *p_dark, char *dark_name,\
	 double dark_itime, int inframe_imno)
{

  /* local variables */
  int count;
  int fstat; 		/* test Function STATus */
  int unit = 0;		/* dummy variable for interface functions */
  int null = 0;
  int actvals;
  int n_rep;		/* number of integrated images with ITIME each */

  double int_time;	/* integration time of the frame to be dark corrected */

  float scale_factor;	/* factor for integration time scaling */

  time_t time_secs;	/* time variable */

  char *p_timestring;	/* pointer for local time */
  char aux_string[200];	/* string buffer */
  


  /* Do only flatfield correction if dark_flag=0*/
  if(dark_flag == 0)
    {

      /* Start */
      fstat =  SCTPUT("\nDoing flatfield correction ...");


      /* loop that does the pixel by pixel division */
      /* test wether input and output is the same */
      if(p_inframe == p_outframe)
	{
	  /* input = output image */

	  for (count=1 ; count <= total_pix ; count++)
	    {
	      /* do test wether vflatfield value is ok, otherwise do not correct */
	      if((*p_flat) > 0.01)
		*p_outframe = (*p_outframe) / (*p_flat);

	      else
		*p_outframe = (*p_outframe);


	      /* pointer incrementation */
	      p_flat++;
	      p_outframe++;
	    }
	}


      else
	{
	  /* input and output are different */

	  for (count=1 ; count <= total_pix ; count++)
	    {
	      /* do test wether flatfield value is ok, otherwise do not correct */
	      if((*p_flat) > 0.01)
		*p_outframe = (*p_inframe) / (*p_flat);

	      else
		*p_outframe = (*p_inframe);


	      /* pointer incrementation */
	      p_inframe++;
	      p_flat++;
	      p_outframe++;
	    }
	}


      /* update descriptors HISTORY of outframe ; create it first if necessary */
      /* get time and date of reduction */
      time(&time_secs);
      p_timestring = ctime(&time_secs);

      /* write time into descriptor */
      fstat = SCDWRC(image_no_out,"HISTORY",1,p_timestring,-1,26,&unit);

      /* original frame */
      fstat = SCDWRC(image_no_out,"HISTORY",1,"Original Image:",-1,80,&unit);
      fstat = SCDWRC(image_no_out,"HISTORY",1,name_inframe,-1,80,&unit);

      /* flatfield frame used*/
      fstat = SCDWRC(image_no_out,"HISTORY",1,"Flatfield correction done with frame: ",-1,80,&unit);
      fstat = SCDWRC(image_no_out,"HISTORY",1,flat_frame,-1,80,&unit);


      return 0;

    }	/* end of do flatfield only */





  /* ------------------------------------------------------------------------------------ */
  /* if dark_flag > 0, do also dark subtraction */



  /* Start */
  fstat =  SCTPUT("\nDoing flatfield correction and dark current subtraction ...");


  /* get integration time of input image */
  fstat = SCDRDD(inframe_imno,"ITIME",1,1,&actvals,&int_time,&unit,&null);

  /* get coadds = number of integrated images with ITIME*/
  fstat = SCDRDI(inframe_imno,"NCOADDS",1,1,&actvals,&n_rep,&unit,&null);

  /* calculate scale factor */
  scale_factor =  (int_time * n_rep) / dark_itime;


  if(DEBUG == 1)
    {
      sprintf(aux_string,"Scale factor for dark current subtraction in frame is set to %5.2f", scale_factor);
      fstat = SCTPUT(aux_string);
    }



  /* do flatfield and dark current correction */
  for (count=1 ; count <= total_pix ; count++)
    {
      /* do test wether flatfield value is ok, otherwise do not correct */
      if((*p_flat) > 0.01)
	*p_outframe = ((*p_inframe) - (scale_factor * (*p_dark))) / (*p_flat);
      
      else
	*p_outframe = (*p_inframe) - (scale_factor * (*p_dark));


      /* pointer incrementation */
      p_inframe++;
      p_flat++;
      p_outframe++;
      p_dark++;
    }



  /* update descriptors HISTORY of outframe ; create it first if necessary */
  /* get time and date of reduction */
  time(&time_secs);
  p_timestring = ctime(&time_secs);

  /* write time into descriptor */
  fstat = SCDWRC(image_no_out,"HISTORY",1,p_timestring,-1,26,&unit);

  /* original frame */
  fstat = SCDWRC(image_no_out,"HISTORY",1,"Original Image:",-1,80,&unit);
  fstat = SCDWRC(image_no_out,"HISTORY",1,name_inframe,-1,80,&unit);

  /* flatfield frame used*/
  fstat = SCDWRC(image_no_out,"HISTORY",1,"Flatfield correction done with frame: ",-1,80,&unit);
  fstat = SCDWRC(image_no_out,"HISTORY",1,flat_frame,-1,80,&unit);

  /* dark frame used*/
  fstat = SCDWRC(image_no_out,"HISTORY",1,"Dark current subtraction done with frame: ",-1,80,&unit);
  fstat = SCDWRC(image_no_out,"HISTORY",1,dark_name,-1,80,&unit);



  return 0;

}

