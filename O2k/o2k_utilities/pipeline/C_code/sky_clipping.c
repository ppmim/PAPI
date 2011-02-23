/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		sky_clipping.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 12 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		sky determination, kappa-sigma clipping process
.PURPOSE	     	determine the sky for an image from a stack of images taken before and
			after the masterimage
.IN/OUTPUT		IN: kappa --> cut-off = median + kappa*sigma
			IN: end_flag = flag for end_reduction (called loose_end in MAIN MODULE)
.RETURNS		Status: 0 = ok
.COMMENTS		Whenever the function Clipping() is called, the current images on the 
			image stack are taken and a sky frame created via a kappa_sigma 
			clipping process. The pixel column is sorted and then the real median
			taken as the local sky level. The squareroot of this level (calculated
			on the electron level) is taken as the sigma of the sky. The 
			cut-off is then set to kappa*sigma above the median level. All values 
			higher then the cut-off are clipped off. The real sky is determined as 
			the median of the values lower than the cut-off.
.CALL			int Clipping(float kappa, int end_flag)
.VERSION		1.00		05.11.02
			1.01		12.11.02	corrected the sigma determination
			2.0		12.11.02	implement SHELL's method for sorting instead
							of qsort
			2.1		06.03.03	changed all output to SCTPUT
			3.0		21.03.03	implemented local stddev
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<midas_def.h>		/* MIDAS definitions */
#include<stdlib.h>		/* for qsort() */
#include<math.h>		/* for sqrt() */
#include<time.h>		/* for CPU time determination */

/* Header File for Structure Definitions */
#include "struct_def.h"

/* Header Files of Secundary Modules */
#include "sky_clipping.h"

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


/* prototype for sorting function */
void shell(unsigned long n, float a[]);


int Clipping(float kappa, int end_flag)
{
  /* Make global variables visible to module */
  extern float image_stack[PIX_AXIS][PIX_AXIS][TOT_IMAGES_MAX];		
  extern float sky_frame[PIX_AXIS][PIX_AXIS];
  extern BOOKPAGE stack_book[TOT_IMAGES_MAX];
  extern int master_pos;
  extern int new_pos;

  extern float sqrt_qe_map[PIX_AXIS][PIX_AXIS];

  extern int SKY_FRAMES;					/* number of images on both sides */
  extern int TOT_IMAGES;					/* total number of images on stack */

		      
  /* Local Variables */
  float column[TOT_IMAGES];			/* array for temporary storage of a pixel column */
  float norm_level = 0;				/* variable for normalization level */
  static float elec_gain = 0;		      	/* electron to counts conversion factor */
  float sigma = 0;				/* standard deviation of sky level */
  float cut_off = 0;				/* cut-off value */

  long time_1, time_2;				/* for CPU-time determination */

  int fstat;								/* function status */
  int test = 0;
  int loop_x = 0;							/* loop variable */
  int loop_y = 0;
  int loop_z = 0;
  int clip_pos = TOT_IMAGES -1;					     	/* position of cut */
  int count = 0;
  int actvals = 0;
  int unit = 0;								/* dummy */
  int null = 0;

  char aux_string[200];							/* string buffer */


  /* Start CPU-time recording */
  time_1 = clock();


  fstat = SCTPUT("\nDetermining sky using kappa-sigma clipping now...");


  /* if this is the first function call, get the electron to count conversion factor */
  /* from the descriptor ELECGAIN */
  if(elec_gain == 0)
    {
      fstat = SCDRDR(stack_book[master_pos].frame_imno,"ELECGAIN",1,1,&actvals,&elec_gain,&unit,&null);
      sprintf(aux_string,"\nELECGAIN is: %f", elec_gain);
      fstat = SCTPUT(aux_string);
    }


  /* get the current normalization factor of image in master_pos */
  /* if masterframe is still open */
  if(end_flag != 2)
    fstat = SCDRDR(stack_book[master_pos].frame_imno,"NORMALIZATION",1,1,&actvals,&norm_level,&unit,&null);

  /* otherwise take norm_level from last image */
  else
    fstat = SCDRDR(stack_book[new_pos].frame_imno,"NORMALIZATION",1,1,&actvals,&norm_level,&unit,&null);



  /* start the kappa-sigma clipping sky determination process */

  /* load the pixel column and sort it */
  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
    {
      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	{
	  for(loop_z=0 ; loop_z<TOT_IMAGES ; loop_z++)
	    {
	      /* store a pixel column in the temporary array */		     
	      column[loop_z] = image_stack[loop_y][loop_x][loop_z];
	    }
	      
	  /* sort pixel column with qsort */
	  shell(TOT_IMAGES, column);


	  /* calculate sigma from median */
	  /* the prefactors give the sigma as the statistical standard deviation */
	  sigma = (float) sqrt((double)(column[SKY_FRAMES] * norm_level / elec_gain));

	  /* calculate local cut-off as median + local stddev */
	  cut_off = column[SKY_FRAMES] + (kappa * sigma / sqrt_qe_map[loop_y][loop_x]);


	  /* clip off values above cut-off */
	  /* start checking from top of sorted column */
	  for(clip_pos = TOT_IMAGES -1 ; clip_pos >= SKY_FRAMES ; clip_pos--)
	    {
	      /* test wether value is above cut-off */
	      if(column[clip_pos] <= cut_off)
		break;

	      test++;							/* test clips */
	    }


	  /* clip_pos now has the position of the last good value */

	  /* if clip_pos is even --> # of good values is odd  */
	  /* --> median of good values is at position: clip_pos / 2 */
	  if((clip_pos%2) == 0)						/* true if clip_pos is even */
	    {		
	      /* store cleaned median in sky_frame */
	      sky_frame[loop_y][loop_x] = column[clip_pos/2];
	      continue;
	    }


	  /* if clip_pos is odd --> # of good values is even */
	  /* --> clean median is average of two middle values */
	  else
	    {
	      count = (clip_pos -1) / 2;			       	/* index of lower middle value */

	      sky_frame[loop_y][loop_x] = (column[count] + column[count+1]) / 2;
	    }


	}	/* end of loop_x */
    }		/* end of loop_y  */



  /* final comments */
  sprintf(aux_string,"The last cut_off was: %f", cut_off);
  fstat = SCTPUT(aux_string); 
  sprintf(aux_string," %d values were clipped, that is %f clips per pixel column\n", test, ((float)test)/TOT_PIX);
  fstat = SCTPUT(aux_string); 
  fstat = SCTPUT("...kappa-sigma clipping sky done");


  /* get used CPU time*/
  time_2 = clock();

  sprintf(aux_string,\
	  "The CPU-time for the sky determination (kappa-sigma clipping mode mode) was: %6.2f seconds\n", \
	  ((float)time_2-(float)time_1)/1000000);
  fstat = SCTPUT(aux_string);


  return 0;


}	     	/* end of function Minimum() */






/* This is an implementation of the SHELL'S METHOD for sorting an array of n elements
   in ascending numerical order
   n is the number of array elements
  'a' is an array that is replaced on output by its sorted rearrangement
   Source: Numerical Recipes in C; Press, Teukolsky, Vetterling, Flannery, 1992 
   For N < 50 this algorithm is very fast

Modification (RF 12.11.02): use regular array indices a[0],...,a[n-1]

CALL: void shell(unsigned long n, float a[])						   	*/

void shell(unsigned long n, float a[])
{
  unsigned long i,j,inc;
  float v;
  inc = 1;	

do
  {
    inc *= 3;
    inc++;
  }while(inc <= n);	

do
  {
    inc /= 3;

    for(i=inc;i<n;i++)	/* changed from: i=inc+1 ; inc <= n */
      {
	v=a[i];
	j=i;

	while(a[j-inc] > v)
	  {
	    a[j] = a[j-inc];
	    j -= inc;
	    
	    if(j < inc)	/* changed from <= */
	      break;
	  }

	a[j] = v;

      }

  }while(inc > 1);

}

