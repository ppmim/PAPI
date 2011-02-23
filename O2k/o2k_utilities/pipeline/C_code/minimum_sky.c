/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		minimum_sky.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 11 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		sky determination, minimum process
.PURPOSE	     	determine the sky for an image from a stack of images taken before and
			after the masterimage
.IN/OUTPUT		IN: n_smallest = number of minimum values for average
.RETURNS		Status: 0 = ok
.COMMENTS		Whenever the function Minimum() is called, the current images on the 
			image stack are taken and a sky frame created via a minimum process.
		 	The sky value is determined as the average of the n smallest values 
			in a column.	The result is stored in the global array 'sky_frame'. 
			2 different algorithms are implemented:
			a) for n_smallest=1: a very fast algotithm is used, the sky is re-
			   determined only if necessary
			b) for n_smallest>1: the column is sorted and then the average taken 
.CALL			int Minimum(int n_smallest)
.VERSION		1.00		31.10.02
			2.00		12.11.02	implementation of SHELL's method as sorting
							algorithm instead of qsort()
			2.1		06.03.03	changed all output to SCTPUT
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<midas_def.h>		/* MIDAS definitions */
#include<stdlib.h>		/* for qsort() */
#include<time.h>		/* for CPU time determination */

/* Header Files of Secundary Modules */
#include "minimum_sky.h"

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


/* prototype of the sorting routine */
/* index _1 for distinction from shell() in clipping() */
void shell_1(unsigned long n, float a[]);	


int Minimum(int n_smallest)
{
  /* Make global variables visible to module */
  extern float image_stack[PIX_AXIS][PIX_AXIS][TOT_IMAGES_MAX];		
  extern float old_image[PIX_AXIS][PIX_AXIS];
  extern float sky_frame[PIX_AXIS][PIX_AXIS];
  extern int old_number;	
  extern int new_pos;

  extern int TOT_IMAGES;				/* total number of images on stack */

		      
  /* Local Variables */
  float column[TOT_IMAGES];			/* array for temporary storage of a pixel column */
  float min = 0;				/* variable for minimum */
  float sum = 0;			      	/* variable to store sum of n_smallest elements */

  long time_1, time_2;				/* for CPU-time determination */

  int fstat;								/* function status */
  int test = 0;
  int loop_x = 0;							/* loop variable */
  int loop_y = 0;
  int loop_z = 0;
  int count = 1;

  char aux_string[200];				/* auxilary string buffer for SCTPUT */

  /* Start CPU-time recording */
  time_1 = clock();



  /*------------------ a) Minimum determination for n_smallest=1 ------------------------*/

  if(n_smallest == 1)
    {

      /* Minimum determination for each pixel, ...done if */
      /* no old image stored yet (imno<0 do not exist) */
      /* or total image number <= 3 */
      if((old_number < 0) || (TOT_IMAGES < 4))    	
	{
	  fstat = SCTPUT("\nDetermining minimum sky now...");

	  /* Get Minimum for each pixel column and store it in sky_frame */
	  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	    {
	      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
		{
		  for(loop_z=0 ; loop_z<TOT_IMAGES ; loop_z++)
		    {
		      /* for first element of column, set min=element */
		      if(loop_z == 0)
			{		
			  min = image_stack[loop_y][loop_x][loop_z];	    	     
			  continue;			/* continue with loop_z=1 */
			}

		      /* if current column element is smaller, adjust min */
		      if(image_stack[loop_y][loop_x][loop_z] < min)
			min = image_stack[loop_y][loop_x][loop_z];

		    }
	      
		  /* store minimum in sky_frame */
		  sky_frame[loop_y][loop_x] = min;

		}
	    }

	  fstat = SCTPUT("...sky done");

	  /* get used CPU time*/
	  time_2 = clock();

	  sprintf(aux_string,"\nCPU-time for the sky determination (in minimum mode) was: %6.2f seconds\n", \
		  ((float)time_2-(float)time_1)/1000000);
	  fstat = SCTPUT(aux_string);
	

	  return 0;

	}



      /* compare new pixel with old pixel value and then decide wether a new */
      /* minimum determination is necessary */
      fstat = SCTPUT("\nDetermining sky now using the fast minimum method...");

      for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	{
	  for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	    {

	      /* if the sky value is still valid, the marked block is omitted */
	      /* and the next pixel column is analysed */

	      /* Case 1: old value is larger than sky */
	      /* --> minimum is still in column */
	      if(old_image[loop_y][loop_x] > sky_frame[loop_y][loop_x] )
		{
		  /* if new value larger than sky, then sky is still good */
		  if(image_stack[loop_y][loop_x][new_pos] > sky_frame[loop_y][loop_x])
		    continue;	    

		  /* if new value smaller, then this is the new sky value */
		  else
		    {
		      sky_frame[loop_y][loop_x] = image_stack[loop_y][loop_x][new_pos];
		      continue;	   
		    }
		}


	      /* Case 2: old value was minimum */
	      /* --> new minimum has to be determined as in above */
	      else
		{	  
		  /*......................................................*/
		  for(loop_z=0 ; loop_z<TOT_IMAGES ; loop_z++)
		    {
		      /* for first element of column, set min=element */
		      if(loop_z == 0)
			{		
			  min = image_stack[loop_y][loop_x][loop_z];	    	     
			  continue;			/* continue with loop_z=1 */
			}

		      /* if current column element is smaller, adjust min */
		      if(image_stack[loop_y][loop_x][loop_z] < min)
			min = image_stack[loop_y][loop_x][loop_z];

		    }

		  /*......................................................*/	      

		  /* store minimum it in sky_frame */
		  sky_frame[loop_y][loop_x] = min;

		  test++;	/* check how often minimum is taken */
		}

	    }		/* end of loop_x for-loop */
	}		/* end of loop_y for-loop */


      /* final comments */
      fstat = SCTPUT("...minimum sky done");

      /* get used CPU time*/
      time_2 = clock();

      sprintf(aux_string,"The CPU-time for the sky determination (in minimum mode) was: %6.2f seconds", \
	      ((float)time_2-(float)time_1)/1000000);
      fstat = SCTPUT(aux_string);

      sprintf(aux_string,"The minimum was %d times determined, this is a percentage of: %d\n",\
	      test, test*100/TOT_PIX);
      fstat = SCTPUT(aux_string);
 

      return 0;

    }		/* end of if(n_smallest == 1) */




  /*------------------- b) Determine average of n_smallest values as sky ----------------*/
  fstat = SCTPUT("\nDetermining sky as average of the smallest n values now...");


  /* load the pixel column, sort it and take average of n_smallest values */
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
	  shell_1(TOT_IMAGES, column);
	 

	  /* sum the n_smallest elements */
	  for(count=0 ; count < n_smallest ; count++)
	    sum = sum + column[count];

	  /* store average in sky_frame */
	  sky_frame[loop_y][loop_x] = sum / n_smallest;

	  /* set sum back to 0 */
	  sum = 0;

	}	/* end of loop_x */
    }		/* end of loop_y  */



  /* final comments */
  fstat = SCTPUT("...minimum sky done");

  /* get used CPU time*/
  time_2 = clock();

  sprintf(aux_string,"The CPU-time for the sky determination (in minimum average mode) was: %6.2f seconds", \
	  ((float)time_2-(float)time_1)/1000000);
  fstat = SCTPUT(aux_string);


  return 0;


}	     	/* end of function Minimum() */




/* This is an implementation of the SHELL'S METHOD for sorting an array of n elements
   in ascending numerical order.
   n is the number of array elements
  'a' is an array that is replaced on output by its sorted rearrangement
   Source: Numerical Recipes in C; Press, Teukolsky, Vetterling, Flannery, 1992 
   For N < 50 this algorithm is very fast

Modification (RF 12.11.02): use regular array indices a[0],...,a[n-1]

CALL: void shell(unsigned long n, float a[])						   	*/

void shell_1(unsigned long n, float a[])
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
