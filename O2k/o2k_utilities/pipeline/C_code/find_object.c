/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		find_object.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		Module 19 ; IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS
.KEYWORDS		find objects, find x-y-move 
.PURPOSE	     	find object positions and determine best x-y-move for summation
.IN/OUTPUT		IN: n_average = number of frames to take median 			
.RETURNS		Status: 0 = ok ; 1 = no object was found
.COMMENTS		FOR MASTERFRAME:
			If this is the first frame to be analysed, an appropriate threshold 
			will be determined by adjusting a test threshold until a number of
			objects in the interval [N_OBJECT,2*N_OBJECT] is found.
			Then the location and FWHM of about n_object stellar objects will be 
			determined be the function Analyze() and stored in structure candidates.
			The average FWHM is determined and with appropriate cuts stellar objects
			are selected and stored in the global structure objects. The seeing of
			the masterframe is determined as the average FWHM of the selected 
			objects and stored in the global variable master_seeing.
			FOR ALL OTHER FRAMES: 
			The predetermined pixel offsets are taken and the object in the object 
			list are searched for around the expected position with an increasing
			search radius up to the maximum SEARCH_RADIUS pixels. 
		   	A found object is then analyzed and checked wether it is the
			same (total intensity). The new pixel offsets are stored in an array.
			The array with the offsets is then analyzed an outliers excluded.
			The real offsets are the averages of the selected offsets, which are 
			then stored in the corresponding frame structure.
.CALL			int Object(int n_average)
.VERSION		1.00   		22.12.02
			2.0			12.02.03 	masterframe part finished		
			3.0			13.02.03	finished working version
			3.1			17.02.03	change intensity to total intensity
								put search parameters in pipe_para.h
								corrected wrong intitial. of stack_pos
			3.2			03.03.03	returns 1 if no objects left after
								object selection
			3.3			05.03.03	made dirty median more stable,with check
			3.4			06.03.03	changed all output to SCTPUT
			3.41			20.03.03	added BACK_CUT para for back_cut
			3.5			20.03.03	corrected no-return-bug
			4.0			21.03.03	updated for local background for Analize
			4.1			26.03.03 	corrected tot_int by subtracting 
								NORM_LEVEL, added parameter ANA_RAD,
								and THRESH_MAX	
			4.2			28.03.03	store real intensity = tot_int/norm_fac
								to make less sensitive to background
								variations --> in Analyze 
-----------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<string.h>
#include<math.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Header Files of Secundary Modules */
#include "find_object.h"

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
extern COSMICS_STACK clean;					
extern OBJECT_INFO objects[];				
extern int top_object;	
extern float cosmics_stack[N_AVERAGE_MAX][PIX_AXIS][PIX_AXIS];
extern float object_threshold;
extern float master_seeing;
extern float sqrt_qe_map[PIX_AXIS][PIX_AXIS];


/* prototypes */
int Analyze(int x_pix, int y_pix, int stack_pos, float back_cut, float norm_fac);


/* Module Variables */
/* define structure for object candidates visible to whole module */
OBJECT_INFO candidates[2*N_OBJECT];
int top_cand = 0;			/* filled top level of candidate structure */





int Object(int n_average)
{
  /*------------------ Local Variables --------------------------------------------------*/
  int fstat;
  int last_pos;						/* position of relevant info on info stack */
  int loop_x , loop_y;					/* loop variables */
  int count;
  int aux;  					
  int x_start;						/* cornerpoints for relevant search area */
  int x_end;
  int y_start;
  int y_end;
  int unit = 0;						/* dummies */
  int null = 0;
  int actvals;
  int n_candidates = 0;					/* number of object candidates */
  int stack_pos;					/* holds stack position */
  int over_flag = 0;					/* for threshold convergence */
  int under_flag = 0;
  int delta_x;						/* stores predetermined pixel offsets */
  int delta_y;
  int x_exp;						/* expected pbject positions */
  int y_exp;
  int search_radius;					/* for object search */
  int median_pos;					/* position of median value */

  float norm_factor;					/* normalization factor for stack */
  float stddev;						/* standard deviation for frame */
  float threshold;					/* threshold for object selection */
  float old_thresh;					/* auxilary variable for threshold det */
  float background_cut;					/* backgound level for object analysis */
  float sum_fwhm = 0;					/* for object analysis */
  float ave_fwhm = 0;					/* average of fwhm */
  float intensity_ratio = 0;				/* ratio of central intensities for object */
 							/* identification */
  float x_off_median = 0;				/* for offset selection */
  float y_off_median = 0;
  float x_off_best = 0;					/* final offsets to be stored */
  float y_off_best = 0;
  float x_diff = 0;					/* difference in offsets */
  float y_diff = 0;

  char aux_string[200];					/* string buffer */


  /* arrays to store newly determined pixel offsets and fwhm of found objects*/
  float x_offset[2*N_OBJECT];
  float y_offset[2*N_OBJECT];
  float fwhm_object[2*N_OBJECT];




  /*------------------- Setup -----------------------------------------------------------*/
  /* initialize last_pos */
  last_pos = sum.total_opened - 1;

  /* initialize stack position */
  stack_pos = (clean.top_stack-1);

    if(stack_pos < 0)
      stack_pos = n_average - 1;	/* top level */


  /* display routine */
  fstat = SCTPUT("\n\nFinding objects now...");
  sprintf(aux_string,"\nnumber of objects to be found: %d", N_OBJECT);
  fstat = SCTPUT(aux_string); 
  sprintf(aux_string,"last_pos is: %d", last_pos);
  fstat = SCTPUT(aux_string);
  sprintf(aux_string,"stack_pos = %d\n", stack_pos);
  fstat = SCTPUT(aux_string);


  /* get cornerpoints for relevant search area */
  /* cut_off 100 pixel=45arcsecs from frame, so that objects are within overlap */
  x_start = 100;
  y_start = 100;
  x_end = PIX_AXIS - 100;
  y_end = PIX_AXIS -100;


  /* get frame statistics and scale it to NORM_LEVEL on stack */
  /* Remark: statistics was done in subframe --> ok, because images are reduced,thus flat*/
  /* median countlevel is set to NORM_LEVEl */
  /* get normalization factor */
  norm_factor = NORM_LEVEL / sum.info[last_pos].countlevel;

  /* read stddev from SUBFRAME_STATISTICS */
  fstat =  SCDRDR(sum.info[last_pos].ima_number,"SUB_STATISTICS",5,1,&actvals,&stddev,&unit,&null );
  stddev = stddev * norm_factor;	/* sttdev for normalized countlevel */




  /* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
  /* &&&&&&&&&&&&&&&&&&&& FIND OBJECT SECTION for MASTERFRAME &&&&&&&&&&&&&&&&&&&&&&&&&& */

  /*--------------------- Select Reference Objects --------------------------------------*/
  /* check wether this is masterframe */
  if(sum.total_opened <= 1)	
    {
      fstat = SCTPUT("\nLooking for appropriate threshold and reference objects...");


      /*.................. Find Approriate Threshold .....................................*/
      /* start with 10 sigma above median */
      threshold = NORM_LEVEL + 10*stddev;

      /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
      do		/* until appropriate number of objects is found */
	{
	  n_candidates = 0;	/* reset */

	  /* find appropriate threshold */    
	  for(loop_y = y_start ; loop_y <= y_end ; loop_y++)
	    {
	      for(loop_x = x_start ; loop_x <= x_end ; loop_x++)
		{

		  if(cosmics_stack[0][loop_y][loop_x] > threshold)
		    {
		      n_candidates++;
		    
		      /* look for next object in different line */
		      loop_x += 20;	/* offset of 20 pixels = 8-9 arcsecs */
		      loop_y += 20;

		      if(n_candidates > 2*N_OBJECT)
			break;
		    }

		}

	      if(n_candidates > 2*N_OBJECT)
		break;
	    }


	  sprintf(aux_string,"\nNumber of candidates: %d at threshold %f", n_candidates, threshold);
	  fstat = SCTPUT(aux_string);


	  /* set new threshold if necessary */
	  if(n_candidates > 2*N_OBJECT)
	    {
	      /* note that threshold was to high */
	      over_flag = 1;

	      /* do only if real threshold was not embraced yet */
	      if(under_flag == 0)
		{
		  /* save old threshold */
		  old_thresh = threshold;

		  /* need higher cutoff */
		  threshold += 2.5*stddev;	/* increase by 2.5 sigma */ 
		}

	    }

	  else if(n_candidates < N_OBJECT)
	    {
      	      /* note that threshold was to low */
	      under_flag = 1;

	      /* do only if real threshold was not embraced yet */
	      if(over_flag == 0)
		{
		  /* save old threshold */
		  old_thresh = threshold;

		  /* need lower cutoff */
		  threshold -= stddev;	/* decrease by 1 sigma */
		}
	    }


	  /* if real threshold was embraced, set threshold in middle and break*/
	  if(under_flag == 1 && over_flag == 1)
	    {
	      /* a good enough threshold to be used */
	      threshold = (threshold + old_thresh)/2;

	      break;	/* leave threshold determination */
	    }


	  /* break if threshold gets too low */
	  if(threshold < NORM_LEVEL)
	    break;

	  /* break if threshold gets too high */
	  /* if too many bright objects are in field, the brightest objects could  */
	  /* all be saturated --> break at 3*background level */
	  if(threshold >  (THRESH_MAX * NORM_LEVEL))
	    break;


	} while(n_candidates < N_OBJECT  || n_candidates > 2*N_OBJECT);
      /* search interval = [n_object , 2*n_object] */
      /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	
 
     /* set global threshold */
      object_threshold = threshold;
      
      sprintf(aux_string,"\n\nObject_threshold in masterframe was set to: %f", object_threshold);
      fstat = SCTPUT(aux_string);



      /*...................... Find Reference Objects in Masterframe .....................*/
      n_candidates = 0;

      for(loop_y = y_start ; loop_y <= y_end ; loop_y++)
	{
	  for(loop_x = x_start ; loop_x <= x_end ; loop_x++)
	    {

	      if(cosmics_stack[0][loop_y][loop_x] > object_threshold)
		{
		  n_candidates++;
		      
		  /* set local background cuts to median+3sigma */
		  background_cut = NORM_LEVEL + (BACK_CUT * stddev / sqrt_qe_map[loop_y][loop_x]);

		  /* call analize function */
		  /* return values: 0 --> object_cand ; 1 --> cosmic ; 2 --> extended object*/
		  /* background cuts used too identify cosmics by counting pixels above cut */
		  fstat = Analyze(loop_x , loop_y, 0, background_cut, norm_factor);

		  /* increment top_cand if object was found */
		  if(fstat == 0)
		    top_cand++;


		  /* look for next object in different line */
		  loop_x += 20;	/* offset of 20 pixels to go beyond object*/
		  loop_y += 20;

		  if(top_cand >= 2*N_OBJECT)
		    break;
		}

	    }

	  if(top_cand >= 2*N_OBJECT)
	    break;			/* break if enough objects were found */
	}



      /*...................... Save Stellar Objects in Object_List .......................*/
      /* top_cand has the real number of objects in candidate structure */
      /* the highest structure index is therefore: (top_cand-1) */

      /* display candidate objects */
      fstat = SCTPUT("\n\nThe object candidates are (cosmics already excluded): \n");
      fstat = SCTPUT("object#    x_pos      y_pos      fwhm");

      for( count=0 ; count <= (top_cand-1) ; count++ )
	{
	  sprintf(aux_string,"%4d %12.3f %10.3f %8.3f",\
		  count, candidates[count].x_pos, candidates[count].y_pos, candidates[count].fwhm );
	  fstat = SCTPUT(aux_string);
	}


      /* do FWHM statistics with candidates */
      /* find average FWHM */
      for( count=0 ; count <= (top_cand-1) ; count++ )
	{
	  sum_fwhm += candidates[count].fwhm;
	}

      ave_fwhm = sum_fwhm / top_cand;

      top_object = 0;		/* reset bookkeeping variable for object list */


      /* write stellar objects in object list */
      /* exclude the extended objects that have a FWHM > 1.5*average */
      /* and objects with FWHM < 0.5*average */
      for( count=0 ; count <= (top_cand-1) ; count++ )
	{
	  if(candidates[count].fwhm > (FWHM_MIN*ave_fwhm) && candidates[count].fwhm < (FWHM_MAX*ave_fwhm))
	    {
	      objects[top_object] = candidates[count];	/* copy structure entry */
	      top_object++;
	    }
	}
      /* top_object = number of objects ; highest index = top_object-1 */


      /* display object results */
      fstat = SCTPUT("\n\n\nThe selected objects are:");
      sprintf(aux_string,"Average FWHM was: %f arcsecs ; object cuts were [%f , %f]",\
	      ave_fwhm, 0.5*ave_fwhm, 1.5*ave_fwhm);
      fstat = SCTPUT(aux_string); 
      sprintf(aux_string,"%d object from %d candidates were selected...", top_object, top_cand);
      fstat = SCTPUT(aux_string);
      fstat = SCTPUT("\n\nobject#    x_pos      y_pos      fwhm    tot_intensity(normalized)");

      for( count=0 ; count <= (top_object-1) ; count++ )
	{
	  sprintf(aux_string,"%4d %12.3f %10.3f %8.3f %12.1f", count, objects[count].x_pos, \
		  objects[count].y_pos, objects[count].fwhm, objects[count].intensity );
	  fstat = SCTPUT(aux_string);
	}


      /* determine seing in masterframe as average of the fwhm of selected objects */
      sum_fwhm = 0; 	/* reset */

      for( count=0 ; count <= (top_object-1) ; count++ )
	{
	  sum_fwhm += objects[count].fwhm;
	}

      /* save seeing in global variable */
      master_seeing = sum_fwhm / top_object;

      sprintf(aux_string,"\nThe seeing in the master frame for the summation process is: %f\n", master_seeing);
      fstat = SCTPUT(aux_string);

      return 0; 		/* done with masterframe */
    }

 




  /* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
  /* &&&&&&&&&&&&&&&&& FIND MOVE SECTION FOR ALL OTHER FRAMES */


  /*------------------ Search for Reference Objects -------------------------------------*/
  /* if this is not the masterframe, look for same objects as in list */
  
  /* get predetermined pixel offsets of current frame as integers offsets */
  delta_x = sum.info[last_pos].delta_x;
  delta_y = sum.info[last_pos].delta_y;


  /* reset the global variable top_cand */
  top_cand = 0;



  /* """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" */
  /* try to find all objects in the object list at the expected position */
  /* the search will be within SEARCH_RADIUS pixels around expected position */
  for(count=0 ; count <= (top_object-1) ; count++)
    {
      /* expected object position */
      /* RA offset */
      x_exp = objects[count].x_pos + delta_x;
      /* declination offset */
      y_exp = objects[count].y_pos - delta_y;



      /* if object is too close to edge of image --> continue with next one */
      /* rejected if closer than 20pix = 8-9 arcsecs from edge */
      if( x_exp < 20  ||  x_exp > (PIX_AXIS-20)  ||  y_exp < 20  ||  y_exp > (PIX_AXIS-20) )
	continue;


      /* set starting search radius */
      search_radius = 5;


  
      /* ====================================================================== */
      /* adjust search radius, in steps of 5 pixels */
      do
	{
	  /* set fstat to 3 to activate the break condition */
	  fstat = 3;


	  /* do search */	
	  for(loop_y = (y_exp-search_radius) ; loop_y <= (y_exp+search_radius) ; loop_y++)
	    {
	      for(loop_x = (x_exp-search_radius) ; loop_x <= (x_exp+search_radius) ; loop_x++)
		{

		  if(cosmics_stack[stack_pos][loop_y][loop_x] > object_threshold)
		    {
		      /* set local background cuts to median+3sigma */
		      background_cut = NORM_LEVEL + (BACK_CUT * stddev / sqrt_qe_map[loop_y][loop_x]);

		      /* call analize function */
		      /* return values: 0-->object_cand ; 1-->cosmic ; 2-->extended object*/
		      /* background cuts used too identify cosmics */
		      /* by counting pixels above cut */
		      fstat = Analyze(loop_x , loop_y, stack_pos, background_cut, norm_factor);
		     

		      /* ############################################################## */
		      /* test wether found object is the same by comparing tot_intensity */
		      if(fstat == 0)	/* object found */
			{		
			  /* calculate intensity ratio */
			  intensity_ratio = objects[count].intensity / candidates[top_cand].intensity;


			  /* accept object if intensities are within a factor of 2 */
			  if(intensity_ratio > INTENSITY_MIN  &&  intensity_ratio < INTENSITY_MAX)
			    {
			      /* save new offsets and seeing in appropriate array */
			      /* X-OFFSET_new = x_found - objects[count].x_pos */
			      x_offset[top_cand] = candidates[top_cand].x_pos - objects[count].x_pos;
			      /* Y_OFFSET_new = -y_found + objects[count].y_pos */
			      y_offset[top_cand] = objects[count].y_pos - candidates[top_cand].y_pos;
			      /* save object fwhm */
			      fwhm_object[top_cand] =  candidates[top_cand].fwhm;

			      /* increment top_cand if object was found and break */
			      top_cand++;
			      break;
			    }

			  /* otherwise continue search and inactivate break condition */
			  else
			    {
			      if(DEBUG == 3)
				{			     
				  sprintf(aux_string,\
					  "Intensity rejection of object %d...intensity ratio was: %4.2f",\
					  count, intensity_ratio);
				  fstat = SCTPUT(aux_string);
				}

			      loop_x += 5;	/* continue search */
			      loop_y += 5;

			      fstat = 4;	/* inactivate break */
			    }

			}	


		      /* if no object was found...*/
		      /* ...continue search at advanced position*/
		      else
			{
			  loop_x += 5;	/* continue search */
			  loop_y += 5;
			}
		      /* ############################################################## */

		    }	/* end of above threshold loop */


		}	/* inner for loop */

	      if(fstat == 0)	/* object found */
		break;

	    }	     	/* outer for loop */


	  if(fstat == 0)	/* object found --> break */
	    break;


	  /* increase radius if object was not found */
	  search_radius += 5;
	  /* printf("\nsearch radius set to %d", search_radius); */

	}while(search_radius <= SEARCH_RADIUS);
      /* ==================================================================== */


    }	/* end of loop over all object in object list */
  /* """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""" */




  /* -------------- Reject frame of no objects were found ------------------------------ */
  /* indicate no objects by returning 1 */
  if(top_cand == 0)
    {
      fstat = SCTPUT("...returning without any objects found");
      return 1;
    }


  /* -------------- Analyze the newly found offsets ------------------------------------ */
  /* top_cand = # of found objects ; top_can-1 = highest index */

  /* display results */
  sprintf(aux_string,"\n\nThe predetermined telescope offsets in pixels in x and y were: %f , %f\n",\
	  sum.info[last_pos].delta_x,  sum.info[last_pos].delta_y);
  fstat = SCTPUT(aux_string);
  sprintf(aux_string,"%d out of %d objects were found\n", top_cand, top_object);
  fstat = SCTPUT(aux_string); 
  sprintf(aux_string,"\nThe found offset are:");
  fstat = SCTPUT(aux_string); 
  sprintf(aux_string,"number  x_offset  y_offset");
  fstat = SCTPUT(aux_string);

  for(count=0 ; count <= (top_cand-1) ; count++)
    {  
      sprintf(aux_string,"%4d %8.2f %8.2f", count, x_offset[count], y_offset[count]);
      fstat = SCTPUT(aux_string);
    }


  /* determine frame seeing from average fwhm of decent stellar objects */
  sum_fwhm = 0; 	/* reset */

  for( count=0 ; count <= (top_cand-1) ; count++ )
      sum_fwhm += fwhm_object[count];
    
  ave_fwhm = sum_fwhm / top_cand;	/* average fwhm of all objects */


  n_candidates = 0;			/* reset count variable */
  sum_fwhm = 0; 			/* reset */


  /* exclude possible outliers */
  for( count=0 ; count <= (top_cand-1) ; count++ )
    {
      /* take fwhm if within 0.5 to 1.5 * average_fwhm */
      if( fwhm_object[count] > (FWHM_MIN*ave_fwhm) &&  fwhm_object[count] < (FWHM_MAX*ave_fwhm))
	{
	  n_candidates++;	/* count selected objects */
	  sum_fwhm += fwhm_object[count];
	}
    }

  /* seeing is the cleaned ave_fwhm */
  ave_fwhm = sum_fwhm / n_candidates;

  /* save seeing in corresponding frame info structure */
  sum.info[last_pos].seeing = ave_fwhm;

  sprintf(aux_string,"\nThe seeing in the current frame number %d for the summation process is: %f\n arcsecs",\
	  last_pos, sum.info[last_pos].seeing);
  fstat = SCTPUT(aux_string);


  /* determine cuts for offset determination for the top_cand found objects */
  /* get median of offsets */
  median_pos = top_cand/2;  	/* initialize position of median value */

  loop_x = 0;	/* reset */
  loop_y = 0;



  if(top_cand >= 3)	/* at least 3 candidates */
    {
      /* x-median */
      do	/* dirty way to find median */
	{
	  n_candidates = 0;		/* reset */

	  for( count=0 ; count <= (top_cand-1) ; count++ )
	    {
	      if(x_offset[count] <= x_offset[loop_x])	/* <= because outliers tend to be larger */
		n_candidates++;      
	    }

	  loop_x++;

	  /* stability check to prevent index to overshoot */
	  if(loop_x >= top_cand  &&  n_candidates != median_pos) /* index too high and no median found*/
	    {
	      loop_x = 0;	/* reset */
	      median_pos--;	/* try lower median position */
	    }

	}while(n_candidates != median_pos);


      x_off_median = x_offset[loop_x-1];		/* approximate median in x (position +- 1) */


      /* reset median_pos */
      median_pos = top_cand/2;

      /* y-median */
      do	/* dirty way to find median */
	{
	  n_candidates = 0;		/* reset */

	  for( count=0 ; count <= (top_cand-1) ; count++ )
	    {
	      if(y_offset[count] <= y_offset[loop_y])
		n_candidates++;      
	    }

	  loop_y++;

	  /* stability check to prevent index to overshoot */
	  /* index too high and no median found*/
	  if(loop_y >= top_cand  &&  n_candidates != median_pos)      
	    {
	      loop_y = 0;	/* reset */
	      median_pos--;	/* try lower median position */
	    }

	}while(n_candidates != median_pos);

      y_off_median = y_offset[loop_y-1];		/* approximate median in y (position +- 1) */


      /* exclude all offsets that deviate more than 3 pixels from median-->wrong objects */
      aux = 0;	/* reset count variable */


      for( count=0 ; count <= (top_cand-1) ; count++ )
	{
	  /* determine difference in offsets */
	  x_diff = x_offset[count] - x_off_median;
	  y_diff = y_offset[count] - y_off_median;

	  /* select offsets */
	  if( x_diff > (-OFFSET_CUT) && x_diff < OFFSET_CUT && y_diff > (-OFFSET_CUT)  && y_diff < OFFSET_CUT )
	    {
	      aux++;	/* count */

	      x_off_best += x_offset[count];
	      y_off_best += y_offset[count];
	    } 

	}

    } 	/* end of if top_can >= 3 */

  /* if less than 3 candidates were found, reject image */
  else
    aux = 0;

  /* if no objects are left, return 0 */
  if(aux == 0)
    {
      fstat = SCTPUT("...returning without any objects found");
      return 1;
    }


  /* determine best offsets as average of selected offsets */
  x_off_best = x_off_best / aux;
  y_off_best = y_off_best / aux;


  sprintf(aux_string,"\nThe old pixel offsets in x and y: %f and %f ",\
	  sum.info[last_pos].delta_x, sum.info[last_pos].delta_y );
  fstat = SCTPUT(aux_string); 
  sprintf(aux_string,"...are now replaced with the newly determined offsets taken from %d objects", aux);
  fstat = SCTPUT(aux_string); 
  sprintf(aux_string,"The best offsets in x and y are: %f and %f\n\n", x_off_best, y_off_best);
  fstat = SCTPUT(aux_string);


  /* replace offsets */
  sum.info[last_pos].delta_x = x_off_best;
  sum.info[last_pos].delta_y = y_off_best;


  /* and copy new offsets also to clean.info */
  clean.info[stack_pos] = sum.info[last_pos];



  return 0;	/* end of main() */
}






/* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
/* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */


/*----------------------- Functions -----------------------------------------------------*/

int Analyze(int x_pix, int y_pix, int stack_pos, float back_cut, float norm_fac)
{
  /* PARAMETERS */
  /* x_pix = x position of detected object pixel ; y_pix = y position */
  /* stack_pos = image stack position of image */
  /* back_cut = background cut for cosmics filtering */
  /* norm_fac = normalization factor to get real total intensity */

  /* function to determine the center position of the object and the approximate FWHM */
  /* FWHM and central intensity are stored in current position of structure candidates */
  /* objects which are too extended (FWHM > 4arcsecs) and cosmics are marked by a return */
  /* value > 0 */


  /* function variables */
  int func_x;				/* function loop variable */
  int func_y;				
  int count_pix = 0;			/* count pixels above threshold */
  int x_high = 0;		       	/* position of highest intensity pixel */
  int y_high = 0;
  int x_cm_coord;			/* center of mass coordinates */
  int y_cm_coord;

  double PI = acos(-1);

  float high_level = 0;			/* holds highest intensity level */
  float cm_x = 0;			/* center of mass in x */
  float cm_y = 0;			/* in y */
  float tot_intensity  = 0;	      	/* total intensity in object */
  float intensity = 0;			/* pixel intensity */
  float fwhm_cut;			/* cut_level for FWHM determination */
  float fwhm_obj;			/* determined fwhm for object */



  /* look for highest intensity pixel in 15x10 pixel environment to upper right */
  /* and count pixels above threshold */
  for(func_y = y_pix ; func_y <= (y_pix+10) ; func_y++)		/* 10 pix up */
    {
      for(func_x = (x_pix-5) ; func_x <= (x_pix+10) ; func_x++)	/* 5 to left and 10 pix to right */
	{
	  /* count object pixels if above background cut*/
	  if(cosmics_stack[stack_pos][func_y][func_x] > back_cut)
	    count_pix++;
	  
	  /* save position of highest intensity pixel */
	  if(cosmics_stack[stack_pos][func_y][func_x] > high_level)
	    {
	      high_level = cosmics_stack[stack_pos][func_y][func_x];	/* save new high_level */
	      x_high = func_x;
	      y_high = func_y;
	    }		
	}
    }


  /* if object was a cosmic, return 1 */
  /* cosmic == less than 10 pixels are above backgound cut */
  if(count_pix < COSMICS_CUT)
    return 1;


  /* if a real object was found, take the intensity maximum position as origin */
  /* do a momentum analysis in a 10 pixel radius region around origin */
  /* center of mass == 1. momentum = sum_over(intensity*coordinate)/ total_intensity */
  /* size of object == 2. momentum */
  /* background = Norm_level , to be subtracted for intensity analysis */
  for(func_y = (y_high-ANA_RAD) ; func_y <= (y_high+ANA_RAD) ; func_y++)
    {
      for(func_x = (x_high-ANA_RAD) ; func_x <= (x_high+ANA_RAD) ; func_x++)
	{
	  /* pixel intensity above background */
	  intensity = cosmics_stack[stack_pos][func_y][func_x] - NORM_LEVEL;

	  /* if intensity is > 0, use for center of mass determination */
	  if(intensity > 0)
	    {
	      /* total intensity */	 
	      tot_intensity += intensity;

	      /* preliminary center of mass position in x and y */
	      cm_x += (intensity * (x_high - func_x));
	      cm_y += (intensity * (y_high - func_y));
	    }

	}
    }


  /* real center of mass position relative to origin */
  cm_x = cm_x / tot_intensity;
  cm_y = cm_y / tot_intensity;


  /* write results (absolute positions) in top level of candidate structure */
  candidates[top_cand].x_pos = cm_x + (float)x_high;		/* absolute position in x */
  candidates[top_cand].y_pos = cm_y + (float)y_high;		/* absolute position in y */

  x_cm_coord = cm_x + x_high;		/* integer coordinate position in x */
  y_cm_coord = cm_y + y_high;		/* integer coordinate position in y */


  /* do a rough FWHM analysis */
  /* set central intensity to center_of_mass_intensity */
  high_level =  cosmics_stack[stack_pos][ y_cm_coord][ x_cm_coord] - NORM_LEVEL;

  /* set half maximum level for FWHM */
  fwhm_cut = high_level / 2;


  /* count pixels above this level in area around center of mass */
  count_pix = 0;		/* reset */

  for(func_y = (y_cm_coord-10) ; func_y <= (y_cm_coord+10) ; func_y++)
    {
      for(func_x = (x_cm_coord-10) ; func_x <= (x_cm_coord+10) ; func_x++)
	{
	  /* pixel intensity above background */
	  intensity = cosmics_stack[stack_pos][func_y][func_x] - NORM_LEVEL;

	  /* if intensity is > 0, use for center of mass determination*/
	  if(intensity > fwhm_cut)
	    count_pix++;
	}
    }


  /* the approximate FWHM is the radius of the area (in square-arcsecs) */
  /* where the intensity was above the fwhm_cut */
  fwhm_obj = 2 * (PIX_SCALE) * sqrt( (double)count_pix / PI );

  /* preselect extended objects with fwhm > 4arcsecs */
  if(fwhm_obj > FWHM_CUT)
    return 2;    

  /* save fwhm and total intensity */
  candidates[top_cand].fwhm = fwhm_obj;
  candidates[top_cand].intensity = tot_intensity / norm_fac;	
  /* get real flux to make less sensitive to background variations*/


  /* 0 is returned if real object candidate was found */
  return 0;
}
