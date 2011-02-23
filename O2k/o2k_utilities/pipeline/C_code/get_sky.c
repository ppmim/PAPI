/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		get_sky.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 6 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		sky determination, median process
.PURPOSE	     	determine the sky for an image from a stack of images taken before and
			after the masterimage
.IN/OUTPUT		none
.RETURNS		Status: 0 = ok
.COMMENTS		Whenever the function Sky() is called, the current images on the 
			image stack are taken and a sky frame created via a median process.
			The function 'select()' is used to determine the median.
			The result is stored in the global array 'sky_frame'. 
			Whenever a comparison with the old_images is possible (e.g, if 
			old_number > 0), a faster process is used: The new and the old pixel
			value are both compared to the median==sky_pixel. If both are on the 
			same side of the median, the sky_pixel is still valid, a new median
			is not needed. If not, the median of the column has to be determined. 
.CALL			int Sky()
.VERSION		1.00		02.09.02
			1.1		23.10.02	modifications for command line parameters
			1.2		06.03.03	changed all output to SCTPUT
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<midas_def.h>		/* MIDAS definitions */
#include<time.h>		/* for CPU time determination */

/* Header Files of Secundary Modules */
#include "get_sky.h"

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

/* Setup for function SELECT() */
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;		/* Macro definition for function SELECT()*/

float select_1(unsigned long k, unsigned long n, float arr[]);		/* SELECT() prototype */
/* index _1 to distinguish it from select() in other modules */




int Sky()
{
  /* Make global variables visible to module */
  extern float image_stack[PIX_AXIS][PIX_AXIS][TOT_IMAGES_MAX];		
  extern float old_image[PIX_AXIS][PIX_AXIS];
  extern float sky_frame[PIX_AXIS][PIX_AXIS];

  extern int old_number;	
  extern int new_pos;

  extern int TOT_IMAGES;				/* total number of images on stack */
  extern int SKY_FRAMES;				/* # of images used before/after masterframe*/
			      
		       

  /* Local Variables */
  const unsigned long median = SKY_FRAMES+1;  	/* define which element is to be selected */		
  const unsigned long total = TOT_IMAGES;	/* input variable for select() */

  float column[TOT_IMAGES];			/* array for temporary storage of a pixel column */

  long time_1, time_2;				/* for CPU-time determination */

  char aux_string[200];				/* auxilary string buffer for SCTPUT */

  int fstat;								/* function status */
  int test = 0;
  int loop_x = 0;							/* loop variable */
  int loop_y = 0;
  int loop_z = 0;


  /* Start CPU-time recording */
  time_1 = clock();


  /* Median determination for each pixel, ...done if */
  /* no old image stored yet (imno<0 do not exist) */
  /* or total image number <= 5 */
  if((old_number < 0) || (TOT_IMAGES < 6))    	
    {
      fstat = SCTPUT("\nDetermining sky now...");

      /* Get Median for each pixel column and store it in sky_frame */
      for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
	{
	  for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	    {
	      for(loop_z=0 ; loop_z<TOT_IMAGES ; loop_z++)
		{
		  /* store a pixel column in temporary array */
		  column[loop_z] = image_stack[loop_y][loop_x][loop_z];	    	     
		}
	      
	      /* determine median and store it in sky_frame */
	      sky_frame[loop_y][loop_x] = select_1(median, total, column);

	    }
	}

      fstat = SCTPUT("...sky done");

      /* get used CPU time*/
      time_2 = clock();

      sprintf(aux_string,"CPU-time for the sky determination was: %6.2f seconds\n\n",\
	      ((float)time_2-(float)time_1)/1000000);
      fstat = SCTPUT(aux_string);

      return 0;
      

    }



  /* compare new pixel with old pixel value and then decide wether a new */
  /* median determination is necessary */
  fstat = SCTPUT("\nDetermining sky now using the fast method...");

  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
    {
      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	{

	  /* if the sky value is still valid, the marked block is omitted */
	  /* and the next pixel column is analysed */

	  /* Case 1: both values smaller than median */
	  if(image_stack[loop_y][loop_x][new_pos] <= sky_frame[loop_y][loop_x] )
	    {
	      if(old_image[loop_y][loop_x] < sky_frame[loop_y][loop_x])
		continue;	    
	    }

	  /* Case 2: both values larger than median */
	  else
	    {	  
	      if(old_image[loop_y][loop_x] > sky_frame[loop_y][loop_x])
		continue;	    
	    }


	  /* ------------------------------------------- */
	  for(loop_z=0 ; loop_z<TOT_IMAGES ; loop_z++)
	    {
	      /* store a pixel column in temporary array */
	      column[loop_z] = image_stack[loop_y][loop_x][loop_z];	    	     
	    }
	  test++;	/* check how often median is taken */
	  
	  /* determine median and store it in sky_frame */
	  sky_frame[loop_y][loop_x] = select_1(median, total, column);
	  /* ------------------------------------------- */

	}
    }

  fstat = SCTPUT("...sky done");

  /* get used CPU time*/
  time_2 = clock();

  sprintf(aux_string,"The CPU-time for the sky determination was: %6.2f seconds", \
	  ((float)time_2-(float)time_1)/1000000);
  fstat = SCTPUT(aux_string); 
  sprintf(aux_string,"\nThe median was %d times determined, this is a percentage of: %d\n",\
	  test, test*100/TOT_PIX);
  fstat = SCTPUT(aux_string);


  return 0;

}






/*-----------------------------------------------------------------------------------------
Function SELECT() 
from book: Numerical Recipes in C, p. 341, (by W.H. Press); Second Edition

SELECT() finds the k_th smallest element of an array[1...n] by rearranging the array
in a way that arr[k] holds the k_th smallest element, with all smaller elements moved to
arr[1...k-1] (in arbitrary order) and all larger elements moved to arr[k+1...n] (in 
arbitrary order).
This is the fastest general method for selection, O(n) (rather than O(n*logn) as 
for quicksort). The method rearranges the array and uses partitioning.

float select(unsigned long k, unsigned long n, float arr[])
returns the k_th smallest element out of an array of n elements  

MODIFIED VERSION (RF:21.8.02):
The code was changed to comply with C-array-indices: arr[0], arr[1],...,arr[n-1]
Now, k is really the k_th smallest array element: k=1 being the smallest,..., k=n the 
largest.
n is now the number of array elements, as it is supposed to. 
-----------------------------------------------------------------------------------------*/

float select_1(unsigned long k, unsigned long n, float arr[])
{
  unsigned long i, ir, j, l, mid;
  float a, temp;

  l=0;			/* modified from: l=1 */
  ir = n-1;		/* modified from: ir=n */

  for(;;)
    {
      if(ir <= l+1)
	{
	  if(ir == l+1 && arr[ir] < arr[l])
	    {
	      SWAP(arr[l],arr[ir])
	    }
	  return arr[k-1];		/* modified from: arr[k] */
	}else{
	  mid=(l+ir) >> 1;
	  SWAP(arr[mid],arr[l+1])
	    if(arr[l] > arr[ir])
	      {
		SWAP(arr[l],arr[ir])
	      }
	  if(arr[l+1] > arr[ir])
	    {
	      SWAP(arr[l+1],arr[ir])
	    }
	  if(arr[l] > arr[l+1])
	    {
	      SWAP(arr[l],arr[l+1])
	    }
	  i=l+1;
	  j=ir;
	  a=arr[l+1];
	  for(;;)
	    {
	      do i++; while(arr[i] < a);
	      do j--; while(arr[j] > a);
	      if(j < i) break;
	      SWAP(arr[i],arr[j])
	    }
	  arr[l+1]=arr[j];
	  arr[j]=a;
	  if(j >= k-1) ir=j-1;		/* modified from: if(j>0k) */
	  if(j <= k-1) l=i;		/* modified from: if(j<=k) */
	}
    } 
}
