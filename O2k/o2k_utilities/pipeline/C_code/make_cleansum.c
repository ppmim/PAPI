/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		make_cleansum.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 15 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		cosmics removal
.PURPOSE	     	summation with cosmics removal via median process
.IN/OUTPUT		IN: n_average = number of images for median
.RETURNS		Status: 0 = ok
.COMMENTS		When Cleansum() is called, a cosmics cleaned sum image of all frames
			on clean.image_stack will be determined. 
			Pixels with the same world coordinates are stored in a pixel-column. 
			If all images contribute to the world coordinate pixel, 
			the median is taken and stored in sum.scratch_sum.
		  	If not all images contribute, the pixel is set to 0. The pixel shifting
 			is done with respect to the reference frame.
.CALL			int Cleansum(int n_average)
.VERSION		1.00		26.11.02
			1.1		02.12.02	changed to get 0 for no overlap
			1.2		03.12.02	corrected for inverse RA orintation on frames
			1.3	 	04.12.02	take cosmics stack out of clean.
			1.4		06.03.03	changed all output to SCTPUT
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include<midas_def.h>		/* MIDAS definitions */
#include<time.h>		/* for CPU time determination */

/* Header Files of Secundary Modules */
#include "make_cleansum.h"

/* Structure Definitions */
#include "struct_def.h"

/* Pipeline Parameters, including DETECTOR */
#include "pipe_para.h"

/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif


/* Setup for function SELECT() */
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;		/* Macro definition for function SELECT()*/

float select_2(unsigned long k, unsigned long n, float arr[]);		/* SELECT() prototype */
/* index _2 to distinguish it from select() in other modules */


 /* Make global variables visible to module */
extern DEEP_STRUCT sum;
extern COSMICS_STACK clean;
extern float cosmics_stack[N_AVERAGE_MAX][PIX_AXIS][PIX_AXIS];



int Cleansum(int n_average)
{

  /* Local Variables */
  const unsigned long median = (n_average-1) / 2;  	/* define which element is to be selected */
  const unsigned long total = n_average;		/* input variable for select() */

  float column[n_average];			/* array for temporary storage of a pixel column */

  long time_1, time_2;				/* for CPU-time determination */

  int fstat;				      	/* function status */
  int loop_x = 0;			      	/* loop variables */
  int loop_y = 0;
  int loop_z = 0;
  int count = 0;
  int x_offset[n_average];		     	/* arrays for storage of integer pixel offsets */
  int y_offset[n_average];
  int next_int;					/* variable for closest integer determination */
  int guard = 0;				/* guards medium-launch */
  int test = 0;					/* test number of selection processes */

  char aux_string[200];				/* string buffer */


  /* Start CPU-time recording */
  time_1 = clock();


  fstat = SCTPUT("\nMaking clean_sum now...");

 
  /* get integer pixel offsets for all frames in clean.image_stack */
  for(count=0 ; count<n_average ; count++)
    {

      /* find closest integer */
      /* x_offset */
      next_int = ceil(clean.info[count].delta_x);		/* next larger integer */

      if(fabs(clean.info[count].delta_x - next_int) < 0.5)	/* if next_int is closest */
	x_offset[count] = next_int;
      else							/* else take the smaller one */
	x_offset[count] = next_int-1;

      /* y_offset */
      next_int = ceil(clean.info[count].delta_y);		/* next larger integer */

      if(fabs(clean.info[count].delta_y - next_int) < 0.5)	/* if next_int is closest */
	y_offset[count] = next_int;
      else							/* else take the smaller one */
	y_offset[count] = next_int-1;

    }


  fstat = SCTPUT("\nworking on summation...");


  /* determine world coordinate pixel column for each pixel in reference frame */
  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
    {
      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	{

	  for(loop_z=0 ; loop_z<n_average ; loop_z++)
	    {

	      /* check wether world coordinate pixel is available */
	      /* in x-direction */
	      if((loop_x + x_offset[loop_z]) < 0  || (loop_x + x_offset[loop_z]) >= PIX_AXIS )
		{
		  sum.scratch_sum[loop_y][loop_x] = -1;	/* no value available */
		  guard = 1;
		  break;	/* go to next pixel */
		}

	      /* in y-direction */
	      if((loop_y - y_offset[loop_z]) < 0  || (loop_y - y_offset[loop_z]) >= PIX_AXIS )
		{
		  sum.scratch_sum[loop_y][loop_x] = -1;	/* no value available */
		  guard = 1;
		  break;	/* go to next pixel */
		}


	      /* store pixel_value in column[] */
	      column[loop_z] = cosmics_stack[loop_z][loop_y-y_offset[loop_z]][loop_x+x_offset[loop_z]];
	    	     
	    }
	  

	  /* determine median and store it in sum.scratch_sum */
	  /* do only if column was fully loaded */
	  if(guard == 0)
	    {
	    sum.scratch_sum[loop_y][loop_x] = select_2(median, total, column);
	    test++;
	    }

	  guard = 0;	/* set back guard for next pixel */
	}
    }



  sprintf(aux_string,"The median was taken %d times\n", test);
  fstat = SCTPUT(aux_string);

  fstat = SCTPUT("...cleansum done");

  /* get used CPU time*/
  time_2 = clock();

  sprintf(aux_string,"CPU-time for the clean_sum determination was: %6.2f seconds\n\n", \
	  ((float)time_2-(float)time_1)/1000000);
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

float select_2(unsigned long k, unsigned long n, float arr[])
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
