/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		subframe_statistics.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 2 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		subframe statistics
.PURPOSE	      	get standard statistics of a subframe
.IN/OUTPUT		IN: pointer to input frame
			IN: image number of input frame
			IN: low_cut in units of sigma
			IN: high_cut in units of sigma
.RETURNS		Status: 0 = ok
.COMMENTS		The function will use the subframe specified in the detector-parameter
			file to compute fundamental statistical properties of the frame:
			MINIMUM, MAXIMUM, MEDIAN, MEAN, STD DEVIATION
			The (real)descriptor SUB_STATISTICS is created and the results stored:
			(The name FR_STAT is not used to avoid confusion with MPIAPHOT)
			SUB_STATISTICS(1): MIN
			SUB_STATISTICS(2): MAX
			SUB_STATISTICS(3): MEDIAN
			SUB_STATISTICS(4): MEAN
			SUB_STATISTICS(5): STD DEVIATION
			The descriptor LHCUTS is set to:
					   	[med - low_cut*stddev , med + high_cut*stddev].
			Furthermore, the results are stored in the keyword OUTPUTR.
			Indices as above. Works for .bdf and .fits frames.
.CALL			int Statistics(float *p_inframe, int image_no,...
 				      		float low_cut, float high_cut)
.VERSION		1.00		15.08.02
			2.0		13.11.02	implementation of flexible LHCUTS
							and updating of descriptor HISTORY
			2.1		26.03.03	initialize npix with PIX_AXIS
			2.2		22.04.04	change name of select to select_0
							to avoid clashes with compilers
							of other systems
--------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<midas_def.h>		/* MIDAS definitions */
		
/* Header Files of Secundary Modules */
#include "subframe_statistics.h"

/* Pipeline Parameters, including DETECTOR */
#include "pipe_para.h"

/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif


/* Macro definition for function SELECT()*/
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;	

/* SELECT() prototype */
float select_0(unsigned long k, unsigned long n, float arr[]);		




int Statistics(float *p_inframe, int image_no, float low_cut, float high_cut)
{
  /* local variables */
  int fstat;			/* test Function STATus */

  int naxis;			/* value of descriptor NAXIS */
  int npix[2];			/* value of descriptor NPIX */

  int actvals;			/* actual number of values returned from function */
  int unit;			/* dummy variable for function */
  int null;			/* number of null values in data (returned from function) */

  int xpix = X_END - X_START;		/* length of subarray in x */
  int ypix = Y_END -Y_START;  		/* length of subarray in y */
  int tot_pix = (xpix+1)*(ypix+1);	/* total number of pixels in subarray */

  float subarray[xpix+1][ypix+1];	/* array the holds the values of the subarray */
  float *p_subframe;	    		/* pointer to subframe elements */
  int i,j;				/* loop variables */

  float med_array[tot_pix];		/* 1d array, that is used for the median selection */
  unsigned long tot_unsign = tot_pix;	/* argument for SELECT(), total number of pixels */
  unsigned long med = tot_unsign/2;	/* med = element that is used as median */
  int index = 0;			/* index used to initialize med_array */

  float low_cuts;			/* variable to update LHCUTS */
  float high_cuts;
  float sum = 0;			/* variables for statistics */
  float min;
  float max = 0;
  float mean;
  float median;
  float stddev;


  /* give out info on screen */
  fstat =  SCTPUT("\nFrame Statistics is determined ...");


  /* check some descriptors of input frame */
  /* check wether NAXIS = 2 */
  fstat = SCDRDI(image_no,"NAXIS",1,1,&actvals,&naxis,&unit,&null);

  if(naxis != 2)
    fstat =  SCTPUT("Warning: Descriptor NAXIS of image is not 2! Statistics is done in first plane.");

  /* get NPIX = number of pixels of each axis */
  fstat = SCDRDI(image_no,"NPIX",1,2,&actvals,npix,&unit,&null);

  if(npix[0]<X_END || npix[0]<Y_END || npix[1]<X_END || npix[1]<Y_END)
    {
      fstat = SCTPUT("Error in function Statistics: Subframe exceeds input frame dimension");
      return 7;		/* error code for: INPUT INVALID */
    }


  /* set npix to known PIX_AXIS value */
  npix[0] = PIX_AXIS;
  npix[1] = PIX_AXIS;


  /* store subframe values in array */
  /* initialize subframe_pointer to address of first element of subframe */
  /* indices of subarray: columns(j) 0 to xpix ; rows(i) 0 to ypix */
  p_subframe = p_inframe + (Y_START-1)*npix[0] + X_START - 1;

  /* rows i ; columns j */
  for(i=0 ; i <= ypix ; i++)
    {

      for(j=0 ; j <= xpix ; j++)  		/* one row */
	{
	subarray[i][j] = *p_subframe;

	med_array[index] = *p_subframe;		/* write all values in 1d array */
	index++;

	p_subframe++;				/* next element in row */
	}

      /* set pointer to first element of next row */
      p_subframe = p_inframe + (Y_START+i)*npix[0] + X_START - 1;
    }



  /* now get the statistics in subarray */

  /* first minimum, maximum and average = mean  */
  min = subarray[0][0];		/* initialize min with first array value */

  for(i=0 ; i <= ypix ; i++)
    {
    for(j=0 ; j <= xpix ; j++)  
      {
	sum += subarray[i][j];

	if( subarray[i][j] < min)
	  min =  subarray[i][j];

	if( subarray[i][j] > max)
	  max =  subarray[i][j];
      }
    }

  mean = sum / tot_pix;


  /* now the standard deviation */
  sum = 0;

  for(i=0 ; i <= ypix ; i++)
    {

    for(j=0 ; j <= xpix ; j++)  
      sum += (subarray[i][j] - mean)*(subarray[i][j] - mean);
      
    }

  sum = sum / (tot_pix-1);
  stddev = (float) sqrt( (double) sum);		/* do appropriate typecasting for sqrt() function */


  /* finally the median: use function SELECT() */
  median = select_0(med, tot_unsign, med_array);


  /* write results in keyword OUTPUTR */
  fstat = SCKWRR("OUTPUTR",&min,1,1,&unit);
  fstat = SCKWRR("OUTPUTR",&max,2,1,&unit);
  fstat = SCKWRR("OUTPUTR",&median,3,1,&unit);
  fstat = SCKWRR("OUTPUTR",&mean,4,1,&unit);
  fstat = SCKWRR("OUTPUTR",&stddev,5,1,&unit);


  /* update the descriptor LHCUTS */
  /* low_cut = function argument ; low_cuts = function variable */
  low_cuts = median - low_cut*stddev;
  high_cuts = median + high_cut*stddev;

  fstat = SCDWRR(image_no,"LHCUTS",&low_cuts,1,1,&unit);
  fstat = SCDWRR(image_no,"LHCUTS",&high_cuts,2,1,&unit);
  fstat = SCDWRR(image_no,"LHCUTS",&min,3,1,&unit);
  fstat = SCDWRR(image_no,"LHCUTS",&max,4,1,&unit);


  /* write the results in the (real) descriptor SUB_STATISTICS */
  fstat = SCDWRR(image_no,"SUB_STATISTICS",&min,1,1,&unit);
  fstat = SCDWRR(image_no,"SUB_STATISTICS",&max,2,1,&unit);
  fstat = SCDWRR(image_no,"SUB_STATISTICS",&median,3,1,&unit);
  fstat = SCDWRR(image_no,"SUB_STATISTICS",&mean,4,1,&unit);
  fstat = SCDWRR(image_no,"SUB_STATISTICS",&stddev,5,1,&unit);

  /* write descriptor help text */
  fstat = SCDWRH(image_no,"SUB_STATISTICS","VALUES:MIN(1),MAX(2),MEDIAN(3),MEAN(4),STDDEV(5)",1,80);


  /* back to calling function */
  return 0;
}




/*-------------------------------------------------------------------------------------
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
Now, k is really the k_th smallest array element: 
							k=1 being the smallest,..., k=n the largest.
n is now the number of array elements, as it is supposed to. 
---------------------------------------------------------------------------------------*/






float select_0(unsigned long k, unsigned long n, float arr[])
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
