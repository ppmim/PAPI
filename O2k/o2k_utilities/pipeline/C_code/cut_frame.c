/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		cut_frame.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 20 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		cut frames 
.PURPOSE	     	cuts and copies the sum frames into an actual size frame
.IN/OUTPUT		IN: flag = action flag: 0=get cuts ; 1=copy frames
			IN: p_outframe = pointer to output frame
			IN: p_xstart, p_xend, p_ystart, p_yend = pointers to the variables
				which contain the start and end points of the subframes
.RETURNS		Status: 0 = ok
.COMMENTS		for flag=0:
			The frame sum.master_sum is used to look for the first and the last
			data value. The pixel coordinates for these values are then written
			into the corresponding variables.
			for flag=1
			The determined subframe of sum.master_sum is copied into the real 
			size output frame. Done only for the masterframe.
.CALL			int Cut(int flag, float *p_outframe,
					int *p_xstart, int *p_xend, int *p_ystart, int *p_yend)
.VERSION		1.00		05.03.03
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "cut_frame.h"

/* Pipeline related parameters/constants */
#include "pipe_para.h"		

/* Header file with structure definitions */
#include "struct_def.h"

/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif

/* make gloabel variables accessible */
extern DEEP_STRUCT sum;
					
		  
int Cut(int flag, float *p_outframe, int *p_xstart, int *p_xend, int *p_ystart, int *p_yend)
{
  /* local variables */
  int fstat;
  int loop_x , loop_y;

  int break_guard = 0;



  /* do action according to the flag */ 

  /* ---------------- FIND CUTS -------------------------------------------------------- */
  if(flag == 0)
    {
      fstat = SCTPUT("\nDetermining cuts for sum frame...");

      /* where no value is stored, sum.master_frame contains a -1 */

      /* find lower left corner of data subframe */
      for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	{
	  for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	    {
	      if(sum.master_sum[loop_y][loop_x] > -0.5)
		{
		  /* store start points */
		  *p_xstart = loop_x;	
		  *p_ystart = loop_y;

		  /* set break guard */
		  break_guard = 1;
		  break;
		}
	    }

	  if (break_guard == 1)
	    break;		/* leave outer loop */
	}	/* outer for loop */




      break_guard = 0;	/* reset */

      /* find upper right corner of data subframe */
      /* start at coordinate 2048,2048 and go backwards */
      for(loop_y=(PIX_AXIS-1) ; loop_y>=0 ; loop_y--)
	{
	  for(loop_x=(PIX_AXIS-1) ; loop_x>=0 ; loop_x--)
	    {
	      if(sum.master_sum[loop_y][loop_x] > -0.5)
		{
		  /* store start points */
		  *p_xend = loop_x;	
		  *p_yend = loop_y;

		  /* set break guard */
		  break_guard = 1;
		  break;
		}
	    }

	  if (break_guard == 1)
	    break;		/* leave outer loop */
	}	/* outer for loop */

      return 0;
    }



  /* ------------------ COPY CUT FRAME ------------------------------------------------- */

  else if(flag == 1)
    {
      fstat = SCTPUT("\nCutting out sum frames...");


      /* do copying */
      for(loop_y = *p_ystart ; loop_y <= *p_yend ; loop_y++)
	{
	  for(loop_x = *p_xstart ; loop_x <= *p_xend ; loop_x++)
	    {
	      *p_outframe = sum.master_sum[loop_y][loop_x];

	      p_outframe++;	/* increment output pointer */
	    }
	}	/* outer for loop */


      return 0;
    }



  else 
    {
      fstat = SCTPUT("\nERROR: wrong action call in module cut_frame...");

      return 1;
    }


}
