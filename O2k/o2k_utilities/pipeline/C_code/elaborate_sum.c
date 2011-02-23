/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER	      	elaborate_sum.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 22 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		elaborate cosmics removal
.PURPOSE	     	summation with cosmics removal via comparison to median image
.IN/OUTPUT		IN: n_average = number of images for median
			IN: kappa_sum = clipping parameter for cosmics removal
.RETURNS		Status: 0 = ok
.COMMENTS		This is the more elaborate version of make_cleansum. The median image
			produced by cleansum in sum.scratch_sum is taken as reference frame.
			From the real statistical level of a pixel and all its neighbors the
			sigma is calculated. The cut_off is then set to max_level+kappa*sigma.
			The median value is then replaced by (sum_of real values)/n_average,
			where the outlier values above the cut_off are replaced by the median
			value. This method is using, that cosmics are always outliers at the
			high value end. The centers of bright stars need special treatment:
			If a cosmics candidate is discovered, it is checked wether the outlier
			might belong to a star. Pixels in the neighborhood with
			radius SEARCH_STAR  around candidates are checked and counted if they
			are above the background level. If >COSMICS_CUT, star pixels are found,
			the value is taken as star and not rejected,--> otherwise cosmic and
			rejected.
			WARNING: Pipeline will crash, if flatfield contains zeros, since values 
			are divided by the sensitivity map, which is sqrt(flatfield).
.CALL			int Elaborate(int n_average, float kappa_sum)
.VERSION		1.00		19.03.03
			2.0		20.03.03	debugged and added local sensitivity sigma,
							made median_max finding faster
			2.1		08.04.03	corrected hanging during ELECGAIN det.
			2.11		02.06.04	added comment on flatfield zeros
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include<midas_def.h>		/* MIDAS definitions */
#include<time.h>		/* for CPU time determination */

/* Header Files of Secundary Modules */
#include "elaborate_sum.h"

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


 /* Make global variables visible to module */
extern DEEP_STRUCT sum;
extern COSMICS_STACK clean;
extern float cosmics_stack[N_AVERAGE_MAX][PIX_AXIS][PIX_AXIS];
extern float sqrt_qe_map[PIX_AXIS][PIX_AXIS];


int Elaborate(int n_average, float kappa_sum)
{

  /* Local Variables */
  float column[n_average];			/* array for temporary storage of a pixel column */
  float median;					/* median stored in sum.scratch_sum */
  float median_max;				/* maximum median in neighborhood */
  float sigma;					/* global stddev */
  float back_sigma;				/* local stddev of background for star detection  */
  float pixel_sum;				/* sum for the replacement of the median */
  float scale_factor;				/* scale factor for real pixel statistics */
  float cut_off;				/* cut off value for cosmics rejection */
  float back_cut;				/* background cut off value */
  static float elec_gain = 0;		      	/* electron to counts conversion factor */

  long time_1, time_2;				/* for CPU-time determination */

  int fstat;				      	/* function status */
  int loop_x = 0;			      	/* loop variables */
  int loop_y = 0;
  int loop_z = 0;
  int neighbor_x = 0;				/* loop variables for neighboring pixels */
  int neighbor_y = 0;
  int repl_z;					/* loop variable for sum replacement */
  int star_x;					/* loop variables for star check */
  int star_y;
  int star_check = 0;				/* count star pixels */
  int mid_level = (n_average-1)/2;  	   	/* middle stack level for the unscaling paras */
  int max_guard = 0;				/* guard for median_max finding */
  int max_x = 0;				/* column were last max was found */
  int count = 0;
  int x_offset[n_average];		     	/* arrays for storage of integer pixel offsets */
  int y_offset[n_average];
  int next_int;					/* variable for closest integer determination */
  int guard = 0;				/* guards medium-launch */
  int test = 0;					/* test number of selection processes */
  int test_2 = 0;
  int test_3 = 0;
  int actvals = 0;				/* dummies */
  int unit = 0;							
  int null = 0;


  char aux_string[200];				/* string buffer */


  /* Start CPU-time recording */
  time_1 = clock();


  fstat = SCTPUT("\nMaking elaborate_sum now...");


  /* if this is the first function call, get the electron to count conversion factor */
  /* from the descriptor ELECGAIN --> only if stack is loaded*/
  if(elec_gain == 0)
    {
      fstat = SCDRDR(clean.info[n_average-1].ima_number,"ELECGAIN",1,1,&actvals,&elec_gain,&unit,&null);
      sprintf(aux_string,"\nELECGAIN is: %f", elec_gain);
      fstat = SCTPUT(aux_string);
    }

      
  /* find scaling factor for proper statistics[counts/electron] = norm_fact/elec_gain */ 	  
  /* use the real countlevel of the middle image */
  scale_factor = (NORM_LEVEL/clean.info[mid_level].countlevel)/elec_gain;

  sprintf(aux_string,"\nThe scale_factor = norm_factor/elec_gain [counts per electron] is: %f", scale_factor);
  fstat = SCTPUT(aux_string);


  /* find standard deviation of background */
  back_sigma =  (float) sqrt((double)(NORM_LEVEL * scale_factor));

  sprintf(aux_string,"The stddev of the scaled background of %d is: %f", NORM_LEVEL, back_sigma);
  fstat = SCTPUT(aux_string);




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




  fstat = SCTPUT("\nworking on elaborate summation and cosmics removal...");


  /* determine world coordinate pixel column for each pixel in reference frame */
  for(loop_y=0 ; loop_y<PIX_AXIS ; loop_y++)
    {
      for(loop_x=0 ; loop_x<PIX_AXIS ; loop_x++)
	{
	  /* reset max_guard at beginning of row */
	  if(loop_x == 0)
	    max_guard = 0;


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
	    	     
	    }	/* end of z-loop */

	  /* ++++++++++++++++++++ pixel column complete ++++++++++++++++++++++++++++++++ */	  


	  /* determine cut_off and the real sum that replaces the median */
	  /* do only if column was fully loaded */
	  if(guard == 0)
	    {
	      median = sum.scratch_sum[loop_y][loop_x];	/* median value as determined by cleansum */

		
	      /* find maximum count level in immediate 3x3 neighborhood */
	      /* if this is the first time for this row or median_max was in */
	      /* discarded pixels --> check full neighborhood */
	      if(max_guard == 0  ||  max_x < loop_x-1)
		{
		  median_max = median;			/* initialize the maximum countlevel */
		  max_x = loop_x;

		  /* set max_guard */
		  max_guard = 1;			/* next time, take faster method */
		  test_3++;


		  for(neighbor_y = loop_y-1 ; neighbor_y <= loop_y+1 ; neighbor_y++)
		    {
		      for(neighbor_x = loop_x-1 ; neighbor_x <= loop_x+1 ; neighbor_x++)
			{
			  /* check wether neighbor pixels are available */
			  if(neighbor_x>=0 && neighbor_y>=0 && neighbor_x<PIX_AXIS && neighbor_y<PIX_AXIS)
			    {
			      /* all neighbor pixels are compared */
			      if( sum.scratch_sum[neighbor_y][neighbor_x] > median_max)
				{			  
				  median_max = sum.scratch_sum[neighbor_y][neighbor_x];
				  max_x = neighbor_x;
				}	
			    }
			  
			}
		    }	/* end of for */

		} 	/* end of if max_guard */


	      /* otherwise only 3 new pixels to left have to be checked for max_median */
	      else
		{
		  /* only new pixels on right end have to be checked */
		  neighbor_x = loop_x+1;
		  
		  for(neighbor_y = loop_y-1 ; neighbor_y <= loop_y+1 ; neighbor_y++)
		    {
		     
		      /* check wether neighbor pixels are available */
		      if(neighbor_y>=0 && neighbor_x<PIX_AXIS && neighbor_y<PIX_AXIS)
			{
			  /* update median_max if higher value is found */
			  if( sum.scratch_sum[neighbor_y][neighbor_x] > median_max)
			    {			  
			      median_max = sum.scratch_sum[neighbor_y][neighbor_x];
			      max_x = neighbor_x;
			    }	
			  
			}
		    }	/* end of for */

		}	/* end of else */

	      /* end of maximum finding */



	      /* set cut off parameters */
	      /* find real statistics */
	      sigma = (float) sqrt((double)(median_max * scale_factor));

	      /* set local cut off */
	      cut_off = median_max + (kappa_sum * sigma / sqrt_qe_map[loop_y][loop_x]);

	      /* set local background cutoff as 3 sigma above background*/
	      back_cut = NORM_LEVEL + (BACK_CUT * back_sigma / sqrt_qe_map[loop_y][loop_x]);


	      
	      /* do the real pixel sum over the whole pixel column*/
	      pixel_sum = 0;			/* reset */



	      for(repl_z=0 ; repl_z<n_average ; repl_z++)
		{
		  /* use real value if under cut_off, else use median if cosmic */
		  if(column[repl_z] < cut_off)
		    pixel_sum += column[repl_z];


		  else	/* there is a cosmics candidate */
		    {
		      star_check = 0;	/* reset */	
	
		      /* check neighborhood for star and count pixels above background */
		      for(star_y = loop_y-SEARCH_STAR ; star_y <= loop_y+SEARCH_STAR ; star_y++)
			{
			  for(star_x = loop_x-SEARCH_STAR ; star_x <= loop_x+SEARCH_STAR ; star_x++)
			    {
			      /* count pixel if above background */
			      if(sum.scratch_sum[star_y][star_x] > back_cut)
				star_check++;
			    }
			}	/* end star check */

		      
		      if(star_check < COSMICS_CUT)	/* this is a cosmic */
			{
			  pixel_sum += median;
			  test++;
			}

		      else				/* this is star flux */
			{
			  pixel_sum += column[repl_z];
			  test_2++;
			}


		    }	/* end of else */

		}		/* end of for() */
		
	      /*replace median value with new and better value,scaled down to norm_level*/
	      sum.scratch_sum[loop_y][loop_x] = pixel_sum / n_average;


	    }		/* end of if guard */


	  guard = 0;	/* set back guard for next pixel */
	  /* ++++++++++++++++++++ next pixel +++++++++++++++++++++++++++++++++++++++++++++ */
	}	/* x-loop */
    }		/* y-loop */



  sprintf(aux_string,"The full 3x3 local neighborhood was checked %d times", test_3);
  fstat = SCTPUT(aux_string);
  sprintf(aux_string,"%d values were clipped for cosmics removal, with kappa_sum = %f ", test, kappa_sum);
  fstat = SCTPUT(aux_string);
  sprintf(aux_string,"%d cosmics-candidates were classified as star flux and not clipped ", test_2);
  fstat = SCTPUT(aux_string);


  fstat = SCTPUT("\n...elaborate sum done");

  /* get used CPU time*/
  time_2 = clock();

  sprintf(aux_string,"CPU-time for the elaborate sum with cosmics removal was: %6.2f seconds\n\n", \
	  ((float)time_2-(float)time_1)/1000000);
  fstat = SCTPUT(aux_string);


  return 0;

}

