/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		OMEGA_pipeline.c
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		MAIN MODULE: IR-Pipeline for OMEGA2000
.LANGUAGE   		C ; MIDAS-C-Application
.ENVIRONMENT		MIDAS: #include<midas_def.h>
.KEYWORDS		pipeline, main funktion
.PURPOSE		control the processing of the icoming images 
.IN/OUTPUT		IN: parameters from PRG
.RETURNS		Program Status
.COMMENTS
.VERSION		1.00		26.08.02
			2.00		22.10.02	use of command line parameters in PRG and
							create FITS files
		 	2.1		29.10.02	handle end images and different sky modes
			2.2		13.11.02	update descriptors
			3.0		20.11.02	implementation of summation of reduced images
			3.1		04.12.02	define cosmics_stack seperately
			3.11		11.02.03	place MIDAS connection after variable decs
			3.2		03.03.03	impelemented summation in reduction process
			3.3		04.03.03	implemented bad-pixel-correction
							& and start using SCTPUT instead of printf()
			3.4		05.03.03	implemented cutting of master_sum
			3.5		06.03.03	corrected sum saving for case where master-sum
							is saved at end and no more images on stack
			4.0		20.03.03	added elaborate sum and parameter KAPPA_SUM
			4.1		21.03.03	added local detector sensitivity map
			4.11		26.03.03	check badpixel mask return values after map.
			4.2		28.03.03	add dark_current correction
			4.3		08.04.03	implemented online mode, do filter checking
-----------------------------------------------------------------------------------------*/



/*++++++++++++++++ 1 INCLUDE HEADER FILES +++++++++++++++++++++++++++++++++++++++++++++++*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#include<time.h>
#include<midas_def.h>		/* MIDAS definitions */

/* Structure Definitions */
#include "struct_def.h"

/* Parameter and Symbolic Constant Files */
#include "pipe_para.h"		/* pipeline related parameters/constants */

/* detector related constants: OPRIME_para.h for OMEGA PRIME */
                                   /* O2000_para.h  for OMEGA 2000  */
/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif


/* Header Files of Secundary Modules */
#include "flat_correction.h"		    		/* Module 1 */
#include "subframe_statistics.h"			/* Module 2 */
#include "normalize_frame.h"				/* Module 3 */
#include "extract_name.h"				/* Module 4 */
#include "load_stack.h"				/* Module 5 */
#include "get_sky.h"					/* Module 6 */
#include "save_sky.h"					/* Module 7 */
#include "subtract_sky.h"				/* Module 8 */
#include "copy_frame.h"			      	/* Module 9 */
#include "original_level.h"				/* Module 10 */
#include "minimum_sky.h"				/* Module 11 */
#include "sky_clipping.h"				/* Module 12 */
#include "frame_info.h"				/* Module 13 */
#include "sum_preparation.h"				/* Module 14 */
#include "make_cleansum.h"				/* Module 15 */
#include "dirty_sum.h"					/* Module 16 */
#include "put_tosum.h"					/* Module 17 */
#include "save_dirty.h"				/* Module 18 */
#include "find_object.h"				/* Module 19 */
#include "cut_frame.h"					/* Module 20 */
#include "badpixel_correction.h"			/* Module 21 */
#include "elaborate_sum.h"				/* Module 22 */



/* Prototypes of Local Functions */
int Delay(int secs);
int Copy_badpix(float *p_inframe);
int QE_map(float *p_flatfield);
int QE_one();



/*++++++++++++++++ 2 GLOBAL VARIABLES +++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* All arrays and variables that are needed for the storage and bookkeeping of the image */
/* stack are defined as global, e.g. all functions can access them */

/*---------------- 2.1 Structure Definitions --------------------------------------------*/
/* Structure definitions are stored in header file: struct_def.h */
/* Secundary modules that use these structures need this header file */


/*---------------- 2.2 Symblolic Constants and Makros -----------------------------------*/

#define TOT_IMAGES_MAX (2*SKY_FRAMES_MAX+1)		 	/* total number of images for stack */




/*---------------- 2.3 Variables and Arrays----------------------------------------------*/
/* bad-pixel-mask */
float badpixel_mask[PIX_AXIS][PIX_AXIS];		/* 2D array for bad pixel mask */

/* sensitivity mask */
/* contains squareroot of global flatfield and is normalize to 1 in statistic window */
float sqrt_qe_map[PIX_AXIS][PIX_AXIS];			/* contains local sigma corrections */


/* single image reduction */
float image_stack[PIX_AXIS][PIX_AXIS][TOT_IMAGES_MAX];	/* stack that can hold an image and the  */
							/* SKY_FRAMES_MAX images before and after  */
float sky_frame[PIX_AXIS][PIX_AXIS];			/* image array that holds the sky values */
						       	/* for the current masterframe */
float old_image[PIX_AXIS][PIX_AXIS];		      	/* image frame to temporarely store */
							/* image to be deleted */
BOOKPAGE stack_book[TOT_IMAGES_MAX];		      	/* vector of structures that hold all */
							/* relevant information on the images */
							/* in the corresponding stack position */	  
int old_number = -1;				       	/* image number of old image */
int top_stack = 0;					/* level to which stack is filled */
int master_pos = 0;		      			/* position of masterframe on stack (index) */
int new_pos = 0, old_pos = 0;			       	/* position of newest/oldest image on stack */
int SKY_FRAMES = 0;					/* #of frames for sky determination */
							/* after init. constant throughout pipeline */
int TOT_IMAGES = 0;					/* total number of stack positions used */
							/* constant after real initialization */
/* summation process */
DEEP_STRUCT sum;					/* declare 'sum' as a sructure of form */
							/* DEEP_STRUCT */
COSMICS_STACK clean;					/* declare clean as the structure for  */ 
							/* cosmics removal */
OBJECT_INFO objects[2*N_OBJECT];			 /* structure for storage of object info */
							/* up to 2*N_OBJECT objects can be stored */
float cosmics_stack[N_AVERAGE_MAX][PIX_AXIS][PIX_AXIS]; /* image stack: images stored in a series*/
							/* instead of pixels in series */
SUM_PAGE rejection[100];				/* declare 'rejection' as a structure for */
						     	/* saving info on rejected images */
int reject_count = 0;				      	/* number of rejected frames for summation */
int top_object = 0;					/* holds highest filled object level */

float object_threshold;					/* threshold for objects selection */
float master_seeing;					/* seeing in master frame for summation */



/*+++++++++++++++ 3 MAIN FUNCTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int main()

{

  /*------------- 3.1 Local Variables ---------------------------------------------------*/
  /*............ 3.1.1 General Variables ................................................*/
  int stat;  	 					/* holds current return status */
  int count;						/* loop variable */


  /*............ 3.1.2 Error handling variables (3.2)....................................*/
  int continue_flag = 0;				/* continuation after error: */
  							/* 1=always,0=stop on error only,-1=on warn.*/
  int log_flag = 2;					/* 2 = write warnings and errors in logfile */
  						      	/* 1 = log only errors */
  int display_flag = 1;					/* 2 = display warnings and errors */
  							/* 1 = display errors */


  /*............ 3.1.3 Variables first used in (3.3).....................................*/
  /* numerical command line parameters are capitalized */
  int no_of_char = 0;					/* number of characters of read in names */
  int actvals = 0;					/* number of returned keyword values */
  int unit = 0, null = 0;				/* dummies for funktions */

  char image_cat[81];					/* path and name of icat: last_icat */
  char image_cat_name[81];				/* only name of icat: P1 */
  char flat_name[81];				    	/* flatfield name: P8*/
  char badpixel_name[81];				/* name of bad-pixel-mask */
  char darkframe_name[81];				/* name for dark current name */
  char aux_string[200];					/* auxilary buffer for SCIPUT output */

  float KAPPA = 0;					/* parameter for outlier mode */
  float CUT_MIN = 1;					/* parameter for LHCUTS:... */
  float CUT_MAX = 3;					/* ...default: -1*sigma , +3*sigma */
  float KAPPA_SUM = 0;					/* cut parameter for elaborate sum */
  float WAIT_TIME = 0;					/* integration time for online reduction */

  int SKY_MODE = 1;					/* mode for sky determination; default=1 */
  int N_SMALLEST = 1;					/* parameter for minimum mode */
  int N_AVERAGE = 7;					/* parameter for summation process: P4 */
  int ACTION_FLAG = 0;					/* flag for summation action: P4 */
  int SUM_SAVE_FLAG = 0;				/* flag for summation saving: P4 */
  int FLAT_FLAG = 0;					/* flag for flatfield correction: P7 */
  int SAVE_SKY = 0;					/* flag for sky saving: P7 */
  int OUTPUT_FLAG = 0;					/* flag for screen output in P7 */
  int ONLINE_FLAG = 0;					/* flag for online reduction: online_flag */   
 

  /*............ 3.1.6 Variables first used in (3.6).....................................*/

  char flat_ident[73];					/* IDENT of flatfield */
  char badpixel_ident[73];				/* IDENT for bad pixel mask */
  char darkframe_ident[73];				/* IDENT for dark current frame */

  float *p_flatfield;					/* pointer to mapped flatfield data */
  float *p_badpixel;					/* pointer for bad pixel mask */
  float *p_darkframe;					/* pointer for dark current frame */

  double dark_itime;					/* variable for dark integration time */

  int flat_imno;					/* image number of flatfield */
  int badpixel_imno;					/* image number of bad pixel mask */
  int darkframe_imno;					/* image number of dark current frame */


  /*............ 3.1.7 Variables first used in (3.7).....................................*/
  /* dummy buffers for name prefixes */
  char sum_buf[5] = {'s','u','m','_'};
  char ave_buf[5] = {'a','v','e','_'};
  char dif_buf[5] = {'d','i','f','_'};
  char fits_buf[6] = {'.','f','i','t','s'};


  /*............ 3.1.10 Variables first used in (3.10)...................................*/
  char inframe_1[81];      				/* names of next image in catalog */
  char ident_1[73];					/* identifier of that image */
  char ident_old[73];					/* identifier of last image */
  char outframe_1[85] = {'r','e','d','_'}; 	     	/* name of output frame */
  char ident_outframe_1[81] = "Reduced Image";	/* ascii identifier of image */
  char skyframe_name[85] = {'s','k','y','_'};		/* name of current sky frame */
  char ident_sky[73] = "Modelled Sky";		/* identifier of sky frame */
  char cunit[40];				       	/* unit of each axis */
  char filter_new[10];					/* buffer for new filter name */
  char filter_old[10];					/* buffer for old filter name */


  float *p_inframe_1;					/* pointer to mapped input image */
  float *p_outframe_1;					/* pointer to created output image */
  float *p_skyframe;					/* pointer to sky frame */

  int cat_no = 0;					/* current number of catalog image */
  int imno_inframe_1;					/* image number of inframe 1 */
  int imno_outframe_1;					/* image number of outframe 1 */
  int imno_skyframe;					/* image number of sky frame */
  int naxis;						/* # dimens of image */
  int npix[2] = {PIX_AXIS,PIX_AXIS};		       	/* size of each dim */
  int wait_time = 0;					/* delay time for online reduction */
  int online_check = 0;					/* count variable for online reduction */
  int current_images = 0;			      	/* current number of images in catalog */
  int last_number = 0;					/* last entry number of catalog */

  double start[2] = {1,1};			    	/* start coordinate, step : ARRAYS! */ 
  double step[2] = {1,1};


  /*........... 3.1.15 Variables used in summation process (3.10.15).....................*/
  int reject_check = 0;					/* check wether frame was rejected */
  int name_index = 2;					/* index for new output image */
  int cut_x_start;					/* start point of cut sum frames */
  int cut_x_end;					/* end point of cut sum frames */
  int cut_y_start;					/* start point of cut sum frames */
  int cut_y_end;					/* end point of cut sum frames */
  int cutframe_imno;					/* image number of cut sum frame */
  int cut_npix[2];					/* pixels per dimension for cut sum frame */ 

  char new_name[92];					/* new name for output image */
  char index_char[2];					/* help buffer for indexing */
  char cutframe_ident[73];				/* IDENT of cut sum frame */
  char cutframe_name[95]= {'c','u','t','_'};	       	/* name for cut sum frame */ 

  float *p_cutframe;					/* pointer to mapped cut_sum  data */
  

  /*........... 3.1.20 Variables first used in (3.20)....................................*/
  int last_image = 0;					/* counts down reduction of last image */
  int loop_counter = 0;				      	/* loop_count variable */
  int loose_end = 0;					/* handles reduction of start/end images */
  							/* 0=first images still to be reduced */
  							/* 1=normal mode ; 2=do end ; 3=end pipe */ 



  /*------------- 3.2 Setup MIDAS connection & Error Handling ---------------------------*/
  SCSPRO("OMEGA_pipeline.c");				/* connect to MIDAS environment */

  /* error handling */
  stat = SCECNT("PUT",&continue_flag,&log_flag,&display_flag);			/* default: 0, 2, 1 */



  /*------------- 3.3 Get Input Parameters from PRG -------------------------------------*/

  /* character strings */
  stat = SCKGETC("P1",1,80,&no_of_char,image_cat_name); 	      	/* get name of image catalog */
  stat = SCKGETC("last_icat",1,80,&no_of_char,image_cat); 	      	/* get name and path of icat */
  stat = SCKGETC("flatfield_key",1,80,&no_of_char,flat_name); 		/* get name of input flatfield */
  stat = SCKGETC("badpixel_key",1,80,&no_of_char,badpixel_name);	/* name of bad-pixel-mask */
  stat = SCKGETC("darkframe_key",1,80,&no_of_char,darkframe_name);	/* name of dark current frame */


  /* integers */
  /* get sky_frames from P3 */
  stat = SCKRDI("frames_key",1,1,&actvals,&SKY_FRAMES,&unit,&null);

  /* get flags from P7 */
  stat = SCKRDI("flag_key",1,1,&actvals,&SAVE_SKY,&unit,&null);
  stat = SCKRDI("flag_key",2,1,&actvals,&FLAT_FLAG,&unit,&null);
  stat = SCKRDI("flag_key",3,1,&actvals,&OUTPUT_FLAG,&unit,&null);

  /* get summation process parameters */
  stat = SCKRDI("n_average_key",1,1,&actvals,&N_AVERAGE,&unit,&null);
  stat = SCKRDI("action_flag_key",1,1,&actvals,&ACTION_FLAG,&unit,&null);
  stat = SCKRDI("sum_save_key",1,1,&actvals,&SUM_SAVE_FLAG,&unit,&null);

  /* flags for online reduction */
  stat = SCKRDI("online_flag",1,1,&actvals,&ONLINE_FLAG,&unit,&null);


  /* real parameters */
  /* get KAPPA_SUM */
  stat = SCKRDR("kappa_sum_key",1,1,&actvals,&KAPPA_SUM,&unit,&null);

  /* get cuts for LHCUTS */
  stat = SCKRDR("cuts_key",1,1,&actvals,&CUT_MIN,&unit,&null); 
  stat = SCKRDR("cuts_key",2,1,&actvals,&CUT_MAX,&unit,&null);

  /* integration time for online reduction */
  stat = SCKRDR("int_time",1,1,&actvals,&WAIT_TIME,&unit,&null);

  /* mixed parameters */
  /* get sky mode */
  stat = SCKRDI("mode_key_sky",1,1,&actvals,&SKY_MODE,&unit,&null);

  /* get sky_mode specifications */
  if(SKY_MODE == 0)							/* test if sky_mode == minimum */
    stat = SCKRDI("mode_key_n",1,1,&actvals,&N_SMALLEST,&unit,&null);	/* set N_SMALLEST */
	
  else if(SKY_MODE == 2)						/* test if sky_mode == outlier */
    stat = SCKRDR("mode_key_kappa",1,1,&actvals,&KAPPA,&unit,&null);   	/* set KAPPA */
       
  else if(SKY_MODE != 1)
	stat = SCTPUT("\nWARNING: Invalid sky_mode parameter! Sky_mode is set to median");
    




  /*------------- 3.4 Check for Invalid Input Parameters --------------------------------*/
  /* check image catalog name */

  /* check: SKY_FRAMES */
  if(SKY_FRAMES > SKY_FRAMES_MAX)
    {
      stat = SCTPUT("\nWARNING: Invalid parameter input");
      sprintf(aux_string,"Maximum allowed number of frames for sky determination is: %d\n", SKY_FRAMES_MAX);
      stat = SCTPUT(aux_string);
      sprintf(aux_string,"Input number of frames for sky determination is: %d\n", SKY_FRAMES);
      stat = SCTPUT(aux_string);
      stat = SCTPUT("SKY_FRAME parameter is set to maximum value ...");
    
      SKY_FRAMES = SKY_FRAMES_MAX;
    }

  if(SKY_FRAMES < 1)
    {
      stat = SCTPUT("\nWARNING: Invalid parameter input");
      sprintf(aux_string,"Minimum allowed number of frames for sky determination is: 1\n");
      stat = SCTPUT(aux_string);
      sprintf(aux_string,"Input number of frames for sky determination is: %d\n", SKY_FRAMES);
      stat = SCTPUT(aux_string);
      stat = SCTPUT("SKY_FRAME parameter is set to 1 ...");
    
      SKY_FRAMES = 1;
    }

  /* check FLAT_FLAG */
  if(FLAT_FLAG < 0 || FLAT_FLAG >= 3)
    {
      stat = SCTPUT("\nWARNING: Invalid parameter input");
      sprintf(aux_string,"Allowed parameters for FLAT_FLAG are: 0,1,2\n");
      stat = SCTPUT(aux_string);
      sprintf(aux_string,"Input parameter is: %d\n", FLAT_FLAG);
      stat = SCTPUT(aux_string);
      stat = SCTPUT("FLAT_FLAG parameter is set to 1 ...");
    
      FLAT_FLAG = 1;
    }

  /* check SAVE_SKY */
  if(SAVE_SKY != 0 && SAVE_SKY != 1)
    {
      stat = SCTPUT("\nWARNING: Invalid parameter input");
      sprintf(aux_string,"Allowed parameters for SAVE_SKY are : 0,1\n");
      stat = SCTPUT(aux_string);
      sprintf(aux_string,"Input parameter is: %d\n", SAVE_SKY);
      stat = SCTPUT(aux_string);
      stat = SCTPUT("SKY_FRAME parameter is set to 0 ...");
    
      SAVE_SKY = 0;
    }

  /* check sky mode */
  if(N_SMALLEST >= 2*SKY_FRAMES || N_SMALLEST < 1)
    {
      stat = SCTPUT("\nWARNING: Invalid parameter input for N_SMALLEST of sky_mode MINIMUM");
      sprintf(aux_string,"Input parameter is outside allowed range: %d\n", N_SMALLEST);
      stat = SCTPUT(aux_string);
      sprintf(aux_string,"Input parameter is set to %d\n", SKY_FRAMES);
      stat = SCTPUT(aux_string);

      N_SMALLEST = SKY_FRAMES;
    }


  /* check CUTS */
  /* take absolute values of cuts */
  CUT_MIN = (float) fabs((double) CUT_MIN);
  CUT_MAX = (float) fabs((double) CUT_MAX);


  /* check N_AVERAGE */
  if(N_AVERAGE > N_AVERAGE_MAX || N_AVERAGE < 3)
    {
      stat = SCTPUT("\nWARNING: Invalid parameter input for summation parameter N_AVERAGE: ");
      sprintf(aux_string,"Input parameter is outside allowed range: 3- %d\n", N_AVERAGE_MAX);
      stat = SCTPUT(aux_string);
      sprintf(aux_string,"Input parameter is set to %d\n", N_AVERAGE_DEFAULT);
      stat = SCTPUT(aux_string);

      N_AVERAGE = N_AVERAGE_DEFAULT;
    }
  /* make sure N_AVERAGE is odd */
  if(N_AVERAGE%2 == 0)	/* if even */
    {
      N_AVERAGE++;
      stat = SCTPUT("\nN_AVERAGE was increased by 1 to get an odd number");
    }


  /* check ACTION_FLAG */
  if(ACTION_FLAG > 2 || ACTION_FLAG < 0)
    {
      stat = SCTPUT("\nWARNING: Invalid parameter input for summation parameter ACTION_FLAG: ");
      sprintf(aux_string,"Input parameter is set to 0\n");
      stat = SCTPUT(aux_string);

      ACTION_FLAG = 0;
    }

  /* check SUM_SAVE_FLAG */
  if(SUM_SAVE_FLAG > 2 || SUM_SAVE_FLAG < 0)
    {
      stat = SCTPUT("\nWARNING: Invalid parameter input for summation parameter SUM_SAVE_FLAG: ");
      sprintf(aux_string,"Input parameter is set to 1\n");
      stat = SCTPUT(aux_string);
      
      SUM_SAVE_FLAG = 1;
    }


  /* check flatfield */



  /* initialize total number of images used for stack */
  TOT_IMAGES = 2*SKY_FRAMES + 1;			/* initialize total image number */




  /*------------- 3.5 Give out the Parameters used for Pipeline -------------------------*/
  stat = SCTPUT("\nWelcome to OMEGA_pipeline");
  stat = SCTPUT("\nThe parameters used for the pipeline are:\n");

  sprintf(aux_string,"input image catalog: %s", image_cat);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"number of frames for sky determination before and after masterframe: %d", SKY_FRAMES);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"total number of images on sky stack: %d", TOT_IMAGES);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"sky mode is: %d", SKY_MODE);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"n = N_SMALLEST is: %d", N_SMALLEST);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"k = KAPPA is: %4.2f", KAPPA);
  stat = SCTPUT(aux_string);
  stat = SCTPUT("[0,n=minimum mode (takes average of n smallest); 1=median mode; 2,k=outlier mode (k-clipping)]");
  sprintf(aux_string,"flag for flatfield correction: %d", FLAT_FLAG);
  stat = SCTPUT(aux_string);
  stat = SCTPUT("[0=no correction ; 1=correction at start; 2=correction at end]");
  sprintf(aux_string,"flag for sky saving: %d", SAVE_SKY);
  stat = SCTPUT(aux_string);
  stat = SCTPUT("[0=no sky saving; 1=sky for each frame will be saved]");
  sprintf(aux_string,"Take median of #n_average for summation process: %d", N_AVERAGE);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"action_flag for summation process: %d", ACTION_FLAG);
  stat = SCTPUT(aux_string);
  stat = SCTPUT("[0=do reduction+summation; 1=summation only; 2=reduction only]");
  sprintf(aux_string,"sum_save_flag for summation process: %d", SUM_SAVE_FLAG);  
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"kappa_sum for laborate cosmics removal is: %4.2f", KAPPA_SUM);  
  stat = SCTPUT(aux_string);
  stat = SCTPUT("[0=write out every n_average_th frame; 1=save only final clean_sum; 2=save master+real+differ]");
  sprintf(aux_string,"cuts for LHCUTS are: -%4.2f , +%4.2f", CUT_MIN ,CUT_MAX);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"online reduction mode is (0=not online ; 1=online): %d ", ONLINE_FLAG);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"flatfield frame used: %s ", flat_name);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"bad-pixel-mask used: %s ", badpixel_name);
  stat = SCTPUT(aux_string);
  sprintf(aux_string,"dark current frame used: %s\n\n", darkframe_name);
  stat = SCTPUT(aux_string);



  /* wait for 8 seconds to check parameters on display */
  stat = Delay(8);



  /*------------- 3.6 Open the Flatfield Image & Bad Pixel Mask & Dark Frame -------------*/
  /* only needed if FLAT_FLAG is not 0 and ACTION_FLAG is not 1 */
  if(FLAT_FLAG != 0 && ACTION_FLAG != 1)
    {
      /* open dark current frame if DARK_SUBTRACTION = 1*/
      if(DARK_SUBTRACTION == 1)
	{
	  stat = SCIGET(darkframe_name,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,2,&naxis,npix,start,step,darkframe_ident,\
			cunit,(char *)&p_darkframe,&darkframe_imno);

	  /* check return status */
	  if(stat != 0)
	    {
	      stat = SCTPUT("\nERROR: return status of opened darkframe image is not 0!");
	      sprintf(aux_string,"Return status of SCIGET is: %d\n", stat);
	      stat = SCTPUT(aux_string);
	      stat = SCTPUT("The pipeline is aborted ...");
	      exit(EXIT_FAILURE);
	    }
	  /* give out warning of dimension of image is not 2 */
	  if(naxis != 2)
	    stat = SCTPUT("\nWARNING: Dimension of darkframe image is not 2!");
	  /* abort pipeline if pixels per dimension differ from detector parameters */
	  if((npix[0] != PIX_AXIS) || (npix[1] != PIX_AXIS))
	    {
	      stat = SCTPUT("\nERROR: Opened darkframe image has not the number of pixels as specified in the");
	      stat = SCTPUT("detector parameter file! The pipeline is aborted ...");
	      exit(EXIT_FAILURE);
	    }
	  /* give out warning if 'start' is not (1,1)*/
	  if((start[0] != 1) || (start[1] != 1))
	    stat = SCTPUT("\nWARNING: Start coordinates for darkframe image are not 1,1!");


	  /* get integration time of darkframe */
	  stat = SCDRDD(darkframe_imno,"ITIME",1,1,&actvals,&dark_itime,&unit,&null);


	  sprintf(aux_string,"\nDark frame %s opened with total integration time %5.2f sec...",\
		  darkframe_name, dark_itime);
	  stat = SCTPUT(aux_string);

	}	/* end of if DARK_SUBTRACTION = 1 */


      /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
      /* open bad-pixel-mask */
      stat = SCIGET(badpixel_name,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,2,&naxis,npix,start,step,badpixel_ident,cunit,\
		    (char *)&p_badpixel,&badpixel_imno);

     /* check return status */
      if(stat != 0)
	{
	  stat = SCTPUT("\nERROR: return status of opened badpixel image is not 0!");
	  sprintf(aux_string,"Return status of SCIGET is: %d\n", stat);
	  stat = SCTPUT(aux_string);
	  stat = SCTPUT("The pipeline is aborted ...");
	  exit(EXIT_FAILURE);
	}
      /* give out warning of dimension of image is not 2 */
      if(naxis != 2)
	stat = SCTPUT("\nWARNING: Dimension of badpixel image is not 2!");
      /* abort pipeline if pixels per dimension differ from detector parameters */
      if((npix[0] != PIX_AXIS) || (npix[1] != PIX_AXIS))
	{
	  stat = SCTPUT("\nERROR: Opened badpixel image has not the number of pixels as specified in the");
	  stat = SCTPUT("detector parameter file! The pipeline is aborted ...");
	  exit(EXIT_FAILURE);
	}
      /* give out warning if 'start' is not (1,1)*/
      if((start[0] != 1) || (start[1] != 1))
	stat = SCTPUT("\nWARNING: Start coordinates for badpixelmask are not 1,1!");


      /* copy to 2D array badpixel_mask */
      stat = Copy_badpix(p_badpixel);


      /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
      /* open .FITS-flatfield in read mode */
      stat = SCIGET(flat_name,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,2,&naxis,npix,start,step,flat_ident,cunit,\
		    (char *)&p_flatfield,&flat_imno);


      /* Check the return values for flatfield:*/
      /* check return status */
      if(stat != 0)
	{
	  stat = SCTPUT("\nERROR: return status of opened flatfield image is not 0!");
	  sprintf(aux_string,"Return status of SCIGET is: %d\n", stat);
	  stat = SCTPUT(aux_string);
	  stat = SCTPUT("The pipeline is aborted ...");
	  exit(EXIT_FAILURE);
	}
      /* give out warning of dimension of image is not 2 */
      if(naxis != 2)
	stat = SCTPUT("\nWARNING: Dimension of input image is not 2!");
     /* abort pipeline if pixels per dimension differ from detector parameters */
      if((npix[0] != PIX_AXIS) || (npix[1] != PIX_AXIS))
	{
	  stat = SCTPUT("\nERROR: Opened flatfield image has not the number of pixels as specified in the");
	  stat = SCTPUT("detector parameter file! The pipeline is aborted ...");
	  exit(EXIT_FAILURE);
	}
      /* give out warning if 'start' is not (1,1)*/
      if((start[0] != 1) || (start[1] != 1))
	stat = SCTPUT("\nWARNING: Start coordinates for flatfield are not 1,1!");



      /* create the local detector sensitivity map */
      stat = QE_map(p_flatfield);

    }	/* end of if flatfield */


  else	/* no flatfield opened */
    {
      stat = QE_one();			/* copy 1's in sensitivity map */
    }




  /*------------- 3.7 Open Output Frames for Summation images ---------------------------*/
  /* define names for summation output images */
  /* the catalog_name is used as name, with corresponding prefixes */
  /* name: sum_image-cat.fits or ave_ or dif_ */

  /* combine name parts */
  strncat(sum.master.name, sum_buf, 4);
  strncat(sum.master.name, image_cat_name, 80);
  strncat(sum.master.name, fits_buf, 5);

  strncat(sum.dirty.name, ave_buf, 4);
  strncat(sum.dirty.name, image_cat_name, 80);
  strncat(sum.dirty.name, fits_buf, 5);

  strncat(sum.difference.name, dif_buf, 4);
  strncat(sum.difference.name, image_cat_name, 80);
  strncat(sum.difference.name, fits_buf, 5);


  /* open frames if summation is done */
  if (ACTION_FLAG != 2)		/* 2 = reduction only */
    {

      stat = SCIPUT(sum.master.name,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,npix,start,step,\
		    "combined sum image of input catalog","",\
		    (char *)&sum.master.p_frame,&sum.master.ima_number); 

      stat = SCTPUT("\nNew frame for the combined and cosmics cleaned sum image created:");
      stat = SCTPUT(sum.master.name);


      /* open frames for average and difference if needed: SUM_SAVE_FLAG = 2 */
      if(SUM_SAVE_FLAG == 2)
	{
	  stat = SCIPUT(sum.dirty.name,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,npix,start,step,\
			"sum image, without cosmics removal","",\
			(char *)&sum.dirty.p_frame,&sum.dirty.ima_number); 

	  stat = SCTPUT("\nNew frame for the combined sum image created (without cosmics removal):");
	  stat = SCTPUT(sum.dirty.name);


	  stat = SCIPUT(sum.difference.name,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,npix,start,step,\
			"difference dirty_sum-cleaned_sum","",\
			(char *)&sum.difference.p_frame,&sum.difference.ima_number); 

	  stat = SCTPUT("\nNew frame for difference of dirty_sum-cleaned image created:");
	  stat = SCTPUT(sum.difference.name);
	}
    }



  /*------------- 3.8 Initialize Relevant Variables -------------------------------------*/
  /* set bookkeeping variables to 0 */
   sum.total_sum = 0;
   sum.total_opened = 0;
   sum.overflow = 99;

   clean.top_stack = 0;
   reject_count = 0;



 





  /*------------- 3.10 Reduce Single Images and do Summation ----------------------------*/

  /*################## loop over images #################################################*/
  for( ; ; )								/* endless loop */
    {
      stat = SCTPUT("\nNew start of the image loop...");


      /*.......... 3.10.1 Online reduction handling ......................................*/
      /* check wether in online mode (=1) and not end of catalog handling */
      if(ONLINE_FLAG == 1  &&  loose_end < 2 )
	{
	  /* get current number of images in icat */
	  stat = SCCSHO(image_cat,&current_images,&last_number);

	  /* display current images in icat */
	  sprintf(aux_string,"\n%d images in active catalog, %d processed...", last_number, cat_no);
	  stat = SCTPUT(aux_string);

	  /* set integer wait time to integration_time/4 */
	  wait_time = WAIT_TIME/4;

	  /* reset online_check */
	  online_check = 0;


	  /* if last_number in icat is not > cat_no --> wait and check again */
	  while(last_number <= cat_no)
	    {
	      stat = SCTPUT("Waiting for next image...");

	      /* wait for a unit */
	      stat = Delay(wait_time);

	      /* increments wait counter */
	      online_check++;

	      /* break if no new image is coming in --> do end reduction */
	      if(online_check > WAIT_MAX)
		break;

	      /* test for new images again */
	      stat = SCCSHO(image_cat,&current_images,&last_number);

	    }	/* end of while() */

	}	/* end of if ONLINE_FLAG */




      /*.......... 3.10.1 Get name and identifier of next image in catalog ...............*/
      stat = SCCGET(image_cat,1,inframe_1,ident_1,&cat_no);	

      /* check wether EOF */
      if(isspace((int) inframe_1[0]))					/* use string-test-function */
	{
	  stat = SCTPUT("\nLast image of catalog loaded...");
	  loose_end = 2;						/* set up final image handling */
	}


      /* if EOF is reached, skip setup steps and go to 3.10.16/20 */
      if(loose_end != 2)						/* true if new image availabe */
	{


	  /*....... 3.10.2 Open this image in read_mode ...................................*/
	  stat = SCIGET(inframe_1,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,2,&naxis,npix,start,step,ident_1,cunit,\
			(char *)&p_inframe_1,&imno_inframe_1);



	  /*....... 3.10.3 Check the return values ........................................*/
	  /* check return status of SCIGET */
	  if(stat != 0)
	    {
	      stat = SCTPUT("\nERROR: return status of opened image is not 0!");
	      sprintf(aux_string,"Return status of SCIGET is: %d\n", stat);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"Error occurred while opening image: %s \n", inframe_1);
	      stat = SCTPUT(aux_string);
	      stat = SCTPUT("The pipeline is aborted ...");
	      exit(EXIT_FAILURE);
	    }

	  /* give out warning of dimension of image is not 2 */
	  if(naxis != 2)
	    stat = SCTPUT("\nWARNING: Dimension of input image is not 2!");

	  /* abort pipeline if pixels per dimension differ from detector parameters */
	  if((npix[0] != PIX_AXIS) || (npix[1] != PIX_AXIS))
	    {
	      stat = SCTPUT("\nERROR: Opened images have not the number of pixels as specified in the");
	      stat = SCTPUT("detector parameter file! The pipeline is aborted ...");
	      exit(EXIT_FAILURE);
	    }

	  /* give out warning if 'start' is not (1,1)*/
	  if((start[0] != 1) || (start[1] != 1))
	    stat = SCTPUT("\nWARNING: Start coordinates are not 1,1!");

	  /* check if identifier has changed */
	  if(top_stack >= 1)
	    {
	      if(strcmp(ident_1, ident_old) != 0)			     	/* compare identifiers */
		{
		  stat = SCTPUT("\nWARNING: Descriptor IDENT has changed!");
		  stat = SCTPUT("IDENT of last and new images are: ");
		  stat = SCTPUT(ident_old);
		  stat = SCTPUT(ident_1);
		}
	    }
	  strcpy(ident_old, ident_1);						/* store new identifier */


	  /* check filter */
	  /* fet filter of new frame */
	  SCDGETC(imno_inframe_1,"FILTER",1,9,&actvals,filter_new);

	  if(top_stack >= 1)
	    {
	      if(strcmp(filter_new, filter_old) != 0)			     	/* compare filters */
		{
		  stat = SCTPUT("\nWARNING: Filter in new image has changed!");
		  stat = SCTPUT("FILTER of last and new images are: ");
		  stat = SCTPUT(filter_old);
		  stat = SCTPUT(filter_new);
		}
	    }
	  strcpy(filter_old, filter_new);				      	/* store new filter */




	  /*....... 3.10.4 Set up output frame ............................................*/
	  /* not needed for summation only: ACTION_FLAG = 1 */
	  if(ACTION_FLAG != 1)
	    {

	      /* initialize name for output frame */
	      stat = Name(inframe_1,outframe_1);

	      /* create an output frame:  FITS_image in write_mode */
	      stat = SCTPUT("\nNew frame to be created:");
	      stat = SCTPUT(outframe_1);

	      /* SCIPUT works only with F_IMA_TYPE */
	      stat = SCIPUT(outframe_1,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,npix,start,step,ident_outframe_1,cunit,\
			    (char *)&p_outframe_1,&imno_outframe_1); 
	  
	      /* check return status of SCIPUT */
	      if(stat != 0)
		{
		  stat = SCTPUT("\nERROR: return status of created image is not 0!");
		  sprintf(aux_string,"Return status of SCIPUT is: %d\n", stat);
		  stat = SCTPUT(aux_string);
		  sprintf(aux_string,"\nError occurred while creating frame: %s \n", outframe_1);
		  stat = SCTPUT(aux_string);
		  stat = SCTPUT("The pipeline is aborted ...");
		  exit(EXIT_FAILURE);
		}

	      /* copy all other descriptors (non standard descr.) to output frame */
	      stat = SCDCOP(imno_inframe_1,imno_outframe_1,3,inframe_1);



	      /*....... 3.10.5 Initialize output frame .....................................*/
	      /* if flatfield correction is supposed to be done in the beginning,initialize */
	      /* output frame with flatfield corrected input frame */
	      /* do scale dark frame subtraction at the same time */
	      if(FLAT_FLAG == 1)
		{
		  stat = Flat(p_inframe_1, p_flatfield, p_outframe_1, TOT_PIX, imno_outframe_1,\
			      inframe_1, flat_name, DARK_SUBTRACTION, p_darkframe, darkframe_name,\
			      dark_itime, imno_inframe_1);
		}
	      /* otherwise the input frame is just copied to output frame */
	      else
		{
		  stat = Copy(p_inframe_1, p_outframe_1, TOT_PIX, imno_outframe_1, inframe_1);
		}


	      /*....... 3.10.6 Unmap and close input frame .................................*/
	      stat = SCFCLO(imno_inframe_1);


	      /*....... 3.10.7 Get frame statistics ........................................*/
	      stat = Statistics(p_outframe_1, imno_outframe_1, CUT_MIN, CUT_MAX);


	      /*....... 3.10.8 Normalize image .............................................*/
	      stat = Normalize(p_outframe_1, imno_outframe_1, TOT_PIX);


	      /*....... 3.10.9 Load image onto stack .......................................*/
	      stat = Load(p_outframe_1,imno_outframe_1, outframe_1, SKY_MODE);


	      /*....... 3.10.10 Do Bad Pixel Correction ....................................*/
	      /* bad pixel correction always done here; better to do flat_correction first */
	      /* if flat_flag = 0 --> no correction is done */
	      /* last loaded image is corrected at level new_pos*/
	      if(FLAT_FLAG != 0)
		stat =Badpixel();
		


	      /*.......... 3.10.11 Setup for saving the sky frames .........................*/
	      if(SAVE_SKY == 1)					/* check wether sky should be saved */
		{
		  stat = Name(inframe_1,skyframe_name);		     	/* initialize name */

		  strcpy(stack_book[new_pos].sky_name, skyframe_name);	/* write name in stack_book */

		  stat = SCTPUT("\nNew sky frame to be created:");
		  stat = SCTPUT(skyframe_name);
      
		  stat = SCIPUT(skyframe_name,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,npix,start,step,ident_sky,cunit,\
				(char *)&p_skyframe,&imno_skyframe);	/* FITS_image in write_mode */

		  /* update stack_book */
		  stack_book[new_pos].sky_imno = imno_skyframe;
		  stack_book[new_pos].p_sky = p_skyframe;
		}

	    }	/* end of if(ACTION_FLAG != 1) */



	  /* (((((((((( BEGIN SUMMATION ONLY ((((((((((((((((((((((((((((((((((((((((((((( */
	  /*.......... 3.10.15 Do summation process only ..................................*/
	  /* if only summation is done (ACTION_FLAG = 1), the whole summation process is */
	  /* done in part 3.10.15, as long as new images are coming in */
	  /* the end-of-catalog-handling is done in 3.10.16 */
	 
	  else if(ACTION_FLAG == 1)
	    {
	      /* get statistics */
	      stat = Statistics(p_inframe_1, imno_inframe_1, CUT_MIN, CUT_MAX);

	      /* get frame info */
	      stat = Info(inframe_1,imno_inframe_1);	
	      sum.total_opened++;


	      /* preparation for cosmics filtering */
	      /* rejection, pixel offset, normalization, and copying to clean stack */
	      reject_check = Preparation(p_inframe_1, N_AVERAGE);
	      /* if = 1 --> frame was rejected */


	      /* find objects and xy-move, if not rejected */
	      /* if 1 is returned, no objects were found --> reject */
	      if(reject_check == 0)
		{
		  reject_check = Object(N_AVERAGE);


		  /* check return value */
		  if(reject_check == 1)
		    {
		      stat = SCTPUT("\nWARNING: No objects to match were found.The following frame is rejected:");
		      stat = SCTPUT(inframe_1);

		      /* add to rejection list */
		      /* copy info to rejection stack */
		      rejection[reject_count] = sum.info[sum.total_opened-1];

		      /* set back top stack position in sum-structure */
		      sum.total_opened--;
		      reject_count++;

		      /* set back clean.top_stack to right level */
		      clean.top_stack--;

		      if(clean.top_stack < 0)
			clean.top_stack = N_AVERAGE - 1;	/* top level of clean stack */
		    }

		}


	      /* do the dirty sum */
	      if(SUM_SAVE_FLAG == 2  && reject_check == 0)		/* do dirty sum for new frame */
		stat = Dirty(p_inframe_1);


	      /* close input frame */
	      stat = SCTPUT("\nClosing frame:");
	      stat = SCTPUT(inframe_1);
	      stat = SCFCLO(imno_inframe_1);



	      /* continue with next frame if reject_check = 1 */
	      if(reject_check == 1)
		continue;



	      /* do clean median if clean_stack is loaded */
	      if(clean.top_stack == 0)		/* true when stack is loaded */
		{	
		  stat = Cleansum(N_AVERAGE);

		  /* do the elaborate sum */
		  stat = Elaborate(N_AVERAGE, KAPPA_SUM);

		  /* set countlevel to real sum of counts */
		  /* do summation if new median image is available */
		  stat = Tosum(SUM_SAVE_FLAG, N_AVERAGE);



		  /* close frame and open new output frame if necessary */
		  if(SUM_SAVE_FLAG == 0)
		    {
		      /* update LHCUTS */
		      stat = Statistics(sum.master.p_frame, sum.master.ima_number, CUT_MIN, CUT_MAX);

		      /* close this frame */
		      stat = SCTPUT("\nWriting out sum frame now:");
		      stat = SCTPUT(sum.master.name);

		      stat = SCFCLO(sum.master.ima_number);


		      /* name: with new index */
		      if(name_index == 2)
			sprintf(new_name,"%d_%s", name_index,sum.master.name);

		      else		
			{		
			  sprintf(index_char,"%d", name_index);
			  new_name[0] = index_char[0];
			}

		      name_index++;
		      
		      /* write new name in sum.master.name */
		      strncpy(sum.master.name,new_name,89);

		      /* open new output frame */
		      stat = SCIPUT(sum.master.name,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,npix,start,step,\
				    "combined sum image of input catalog","",\
				    (char *)&sum.master.p_frame,&sum.master.ima_number); 

		      stat = SCTPUT("\nNew frame for the combined and cosmics cleaned sum image created:");
		      stat = SCTPUT(sum.master.name);
		    }

		}	/* end of clean median */


	      /* test variables */
	      sprintf(aux_string,"\nsum.total_opened = %d ", sum.total_opened);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"sum.total_sum = %d ", sum.total_sum);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"reject_count = %d ", reject_count);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"x-pixel-offset of last image = %f ", sum.info[sum.total_opened-1].delta_x);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"y-pixel-offset of last image = %f ", sum.info[sum.total_opened-1].delta_y);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"ncount_level of last image = %f ", sum.info[sum.total_opened-1].countlevel);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"sum.sum_level = %f ", sum.sum_level);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"clean.top_stack = %d", clean.top_stack);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"test overflow = %d \n", sum.overflow);
	      stat = SCTPUT(aux_string);


	      /* go to beginnig of image loop */
	      continue;

	    }	/* end of summation_only */


	  else
	    {
	      SCESIG("OMEGA_pipeline","ERROR: undefined ACTION (reduction or summation)",stat);
	      SCETER(ERR_INPINV,"The pipeline is aborted...");
	    }

	}	/* end of if(loose_end != 2) */



      /*.......... 3.10.16 Summation Process: End-of-Catalog .............................*/
      /* only for summation only: ACTION_FLAG = 1 */
      /* loose_end = 2 --> no more new input images */
      if(ACTION_FLAG == 1)
	{
	  stat = SCTPUT("\nDoing summation of last images now...");

	  /* take clean_stack median */
	  if(clean.top_stack > 0)	/* unused images on cosmics_stack */	
	    {
	      stat = Cleansum(N_AVERAGE);

	      /* do the elaborate sum */
	      stat = Elaborate(N_AVERAGE, KAPPA_SUM);

	      /* set countlevel to real sum of counts and do summation */
	      /* arg1 = 0 --> frame is written to output */
	      /* new frames in position: 0,1,...,top_stack-1 */
	      stat = Tosum(0,clean.top_stack);	

	      /* update LHCUTS */
	      stat = Statistics(sum.master.p_frame, sum.master.ima_number, CUT_MIN, CUT_MAX);

	      /* close this frame */
	      stat = SCTPUT("\nWriting out sum frame now:");
	      stat = SCTPUT(sum.master.name);

	      stat = SCFCLO(sum.master.ima_number);

	    }

	  /* if clean.top_stack == 0 and master_sum is only saved at end */
	  /*(sum_save_flag = 1 or 2) */
	  else if(clean.top_stack == 0  && SUM_SAVE_FLAG >= 1)
	    {
	      /* arg1 = 0 --> frame is written to output */
	      /* new frames in position: 0,1,...,top_stack-1 */
	      stat = Tosum(0,0);	

	      /* update LHCUTS */
	      stat = Statistics(sum.master.p_frame, sum.master.ima_number, CUT_MIN, CUT_MAX);

	      /* close this frame */
	      stat = SCTPUT("\nWriting out sum frame now:");
	      stat = SCTPUT(sum.master.name);

	      stat = SCFCLO(sum.master.ima_number);
	    }


	  /* give out warning if total number of images is smaller than n_average */
	  if(sum.total_opened < N_AVERAGE)
	    {	   
	      stat = SCTPUT("\nWARNING: Number of input images for summation is smaller than");
	      stat = SCTPUT("number of images used for median process. Cosmics removal did not work...");
	    }


	  /* write out last dirty_sum and difference */
	  if(SUM_SAVE_FLAG == 2)
	    stat = Savedirty();


	  /* ------------- CUT SUM IMAGE ------------------------------------------------- */
	  /* cut sum images to actual size */
	  /* find appropriate size; flag = 0 */
	  stat = Cut(0 ,p_cutframe, &cut_x_start, &cut_x_end, &cut_y_start, &cut_y_end);

	  /* write dimensions */
	  cut_npix[0] = cut_x_end - cut_x_start + 1;	  
	  cut_npix[1] = cut_y_end - cut_y_start + 1;

	  /* create new name */
	  /* --> cut_lastsumname.fits */
	  strncpy(cutframe_name+4,sum.master.name,90);

	  /* create new cut frame */
	  stat = SCTPUT("\nNew frame for the cut sum image is created with name:");
	  stat = SCTPUT(cutframe_name);
	  sprintf(aux_string,"The pixel dimension of the new frame are: %d , %d", cut_npix[0], cut_npix[1]);
	  stat = SCTPUT(aux_string);

	  stat = SCIPUT(cutframe_name,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,cut_npix,start,step,\
			"cut and combined sum image of input catalog",cunit,\
			(char *)&p_cutframe,&cutframe_imno); 

	  /* open last master_sum image again for the descriptors */
	  stat = SCIGET(sum.master.name,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,2,&naxis,npix,start,step,\
			cutframe_ident,cunit,\
			(char *)&sum.master.p_frame,&sum.master.ima_number); 

	  /* copy descriptors from master_sum frame */
	  /* copy all other descriptors (non standard descr.) to output frame */
	  stat = SCDCOP(sum.master.ima_number,cutframe_imno,3,inframe_1);
	  /* copy descriptor LHCUTS */
	  stat = SCDCOP(sum.master.ima_number,cutframe_imno,4,"LHCUTS");

	  /* close master_sum frame */
	  stat = SCFCLO(sum.master.ima_number);

	  /* add descriptor CUT_SUBFRAME with subframe coordinates */
	  sprintf(aux_string,"The frame is a subframe of %s, with low-left and up-right coord: %d,%d ; %d.%d "\
		  , sum.master.name, cut_x_start, cut_y_start, cut_x_end, cut_y_end);
	  SCDWRC(cutframe_imno,"CUT_SUBFRAME",1,aux_string,1,200,&unit);

	  /* do copying of subframe; flag = 1 */
	  stat =  Cut(1 ,p_cutframe, &cut_x_start, &cut_x_end, &cut_y_start, &cut_y_end);

	  /* close this frame */
	  stat = SCTPUT("\nWriting out cut sum frame now:");
	  stat = SCTPUT(cutframe_name);

	  stat = SCFCLO(cutframe_imno);

	  /* ----------------------------------------------------------------------------- */



	  /* display rejected images */
	  stat = SCTPUT("\n\nThe following images were rejected from summation: ");

	  for(count=0 ; count < reject_count ; count++)
	    {
	    stat = SCTPUT(rejection[count].ima_name);
	    sprintf(aux_string,"x,y-pixel-offsets were to big or no objects were found: %f, %f\n", \
		   rejection[count].delta_x, rejection[count].delta_y );
	    stat = SCTPUT(aux_string);
	    }


	  /* back to PRG */
	  stat = SCTPUT("\nSummation process finished, returning to PRG...");
	  
	  SCSEPI();
	  exit(0);
	}
      /* )))))))))) END OF SUMMATION ONLY ))))))))))))))))))))))))))))))))))))))))))))))) */





      /*%%%%%%%%%%% MAIN REDUCTION 3.10.20 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
      /* not needed for summation only: ACTION_FLAG = 1 */
      /* start the main reduction of images when stack is fully loaded */
      if(top_stack == TOT_IMAGES  &&  ACTION_FLAG != 1)  	/* fully loaded stack for reduction*/
	{

	  /*.......... 3.10.20.1 Determine sky for masterframe .........................*/
	  /* switch is selecting the sky_mode */
	  switch(SKY_MODE)
	    {
	    case 0:		/* minimum mode: sky = average of n-smallest elements */
	      {
		stat = Minimum(N_SMALLEST);
		break;
	      }

		  
	    case 1:	 	/* median mode: sky = median */
	      {
		stat = Sky();
		break;
	      }
	      

	    case 2:    		/* outlier mode: sky = median of elements without outliers */
	      {
		stat = Clipping(KAPPA, loose_end);
		break;
	      }
		  
	    }		/* end of switch */



	  /* rest of 3.10.20 is embedded in a do-while loop to allow for the reduction */
	  /* of images at beginning and end of catalog */
	  do	/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/	
	    {
	      /* get new master_pos for end-image-reduction */
	      if(loose_end == 2)
		master_pos = (master_pos +1)%TOT_IMAGES;


	      /*.......... 3.10.20.2 Save skyframe and then unmap and close frame ..........*/
	      if(SAVE_SKY == 1)
		{
		  /* update LHCUTS and statistics */
		  /* use the sky_frame for statistics, because the frame to be saved is */
		  /* still empty; results are saved in target frame */
		  stat = Statistics(&sky_frame[0][0], stack_book[master_pos].sky_imno, CUT_MIN, CUT_MAX);
		  
		  stat = Savesky(stack_book[master_pos].sky_imno, stack_book[master_pos].p_sky);
		}

		  
	      /*.......... 3.10.20.3 Subtract sky in master frame ..........................*/
	      stat = Subtract(SKY_MODE, N_SMALLEST, KAPPA);



	      /*.......... 3.10.20.10 Do a late flatfield correction .......................*/
	      /* flat correction of sky subtracted image */
	      if(FLAT_FLAG == 2)	
		{
		  stat = Flat(stack_book[master_pos].p_frame, p_flatfield, stack_book[master_pos].p_frame,\
			      TOT_PIX, stack_book[master_pos].frame_imno, inframe_1,flat_name,0,
			      p_flatfield, flat_name, 0, 0);
		  /* last 4 arguments only dummies, since no dark subtraction is done */
		}


	      /*.......... 3.10.20.20 Set count level back .................................*/
	      stat = Original(stack_book[master_pos].p_frame, stack_book[master_pos].frame_imno , TOT_PIX);
	


	      /*.......... 3.10.20.30 Update statistic and LHCUTS of reduced frame .........*/
	      stat = Statistics(stack_book[master_pos].p_frame, stack_book[master_pos].frame_imno,\
				CUT_MIN, CUT_MAX);





	      /* §§§§§§§§§§ START SUMMATION IN REDUCTION §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ */
	      /*.......... 3.10.20.40  Summation Process ...................................*/
	      /* summation in combination with reduction */
	      /* as long as new images are coming in: --> always */
	      /* the end-of-catalog-handling is done in 3.10.20.50 */

	      if(ACTION_FLAG == 0)	
		{
		  /* do all procedures with objects in stack_book[master_pos]. */
		  /* p_frame (* float) ; frame_imno (int) ; frame_name (* char) */

		  /* statistics is already up to date */

		  /* get frame info */
		  stat = Info(stack_book[master_pos].frame_name, stack_book[master_pos].frame_imno);	
		  sum.total_opened++;


		  /* preparation for cosmics filtering */
		  /* rejection, pixel offset, normalization, and copying to clean stack */
		  reject_check = Preparation( stack_book[master_pos].p_frame, N_AVERAGE);
		  /* if = 1 --> frame was rejected */


		  /* find objects and xy-move, if not rejected */
		  /* if 1 is returned, no objects were found --> reject */
		  if(reject_check == 0)
		    {
		      reject_check = Object(N_AVERAGE);


		      /* check return value */
		      if(reject_check == 1)
			{
			  stat = SCTPUT("\nWARNING: No objects to match were found.");
			  stat = SCTPUT("The following frame is rejected:");
			  stat = SCTPUT(stack_book[master_pos].frame_name);

			  /* add to rejection list */
			  /* copy info to rejection stack */
			  rejection[reject_count] = sum.info[sum.total_opened-1];

			  /* set back top stack position in sum-structure */
			  sum.total_opened--;
			  reject_count++;

			  /* set back clean.top_stack to right level */
			  clean.top_stack--;

			  if(clean.top_stack < 0)
			    clean.top_stack = N_AVERAGE - 1;	/* top level of clean stack */
			}
		    }		/* end of reject_check == 0 */




		  /* continue with next frame if reject_check = 1 */
		  /* no continue commands allowed because of variable handling at end */
		  if(reject_check == 0)
		    {

		      /* do the dirty sum */
		      if(SUM_SAVE_FLAG == 2)		/* do dirty sum for new frame */
			stat = Dirty(stack_book[master_pos].p_frame);


		      /* do clean median if clean_stack is loaded */
		      if(clean.top_stack == 0)		/* true when stack is loaded */
			{	
			  stat = Cleansum(N_AVERAGE);
		  
			  /* do the elaborate sum */
			  stat = Elaborate(N_AVERAGE, KAPPA_SUM);

			  /* set countlevel to real sum of counts */
			  /* do summation if new median image is available */
			  stat = Tosum(SUM_SAVE_FLAG, N_AVERAGE);



			  /* close frame and open new output frame if necessary */
			  if(SUM_SAVE_FLAG == 0)
			    {
			      /* update LHCUTS */
			      stat = Statistics(sum.master.p_frame, sum.master.ima_number, CUT_MIN, CUT_MAX);

			      /* close this frame */
			      stat = SCTPUT("\nWriting out sum frame now:");
			      stat = SCTPUT(sum.master.name);

			      stat = SCFCLO(sum.master.ima_number);


			      /* name: with new index */
			      if(name_index == 2)
				sprintf(new_name,"%d_%s", name_index,sum.master.name);

			      else		
				{		
				  sprintf(index_char,"%d", name_index);
				  new_name[0] = index_char[0];
				}

			      name_index++;
		      
			      /* write new name in sum.master.name */
			      strncpy(sum.master.name,new_name,89);

			      /* open new output frame */
			      stat = SCIPUT(sum.master.name,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,npix,start,step,\
					    "combined sum image of input catalog","",\
					    (char *)&sum.master.p_frame,&sum.master.ima_number); 

			      stat = SCTPUT("\nNew frame for the cosmics cleaned sum image created:");
			      stat = SCTPUT(sum.master.name);
			    }

			}	/* end of clean median */


		      /* test variables */
		      sprintf(aux_string,"\nsum.total_opened = %d ", sum.total_opened);
		      stat = SCTPUT(aux_string);
		      sprintf(aux_string,"sum.total_sum = %d ", sum.total_sum);
		      stat = SCTPUT(aux_string);
		      sprintf(aux_string,"reject_count = %d ", reject_count);
		      stat = SCTPUT(aux_string);
		      sprintf(aux_string,"x-pixel-offset of last image = %f ",\
			      sum.info[sum.total_opened-1].delta_x);
		      stat = SCTPUT(aux_string);
		      sprintf(aux_string,"y-pixel-offset of last image = %f ", \
			      sum.info[sum.total_opened-1].delta_y);
		      stat = SCTPUT(aux_string);
		      sprintf(aux_string,"count_level of last image = %f ",\
			      sum.info[sum.total_opened-1].countlevel);
		      stat = SCTPUT(aux_string);
		      sprintf(aux_string,"sum.sum_level = %f ", sum.sum_level);
		      stat = SCTPUT(aux_string);
		      sprintf(aux_string,"clean.top_stack = %d", clean.top_stack);
		      stat = SCTPUT(aux_string);
		      sprintf(aux_string,"test overflow = %d \n", sum.overflow);
		      stat = SCTPUT(aux_string);
		      
		    } 	/* end of if reject_check == 0 */
		}	/* end of if ACTION_FLAG == 0 */
	      



	      /*.......... 3.10.20.50  Summation End of Catalog ............................*/
	      /* summation in combination with reduction */
	      /* test for rejection (0 = not rejected) and loose_end (2 = end of catalog)*/
	      if((ACTION_FLAG == 0)  && (loose_end == 2) && (last_image == SKY_FRAMES-1))       	
		{
		  /* loose_end = 2 --> no more new input images */
		  /* Careful: there are still reduced images coming in */
		  /* no more images <--> last_image == SKY_FRAMES-1 */

		  stat = SCTPUT("\nDoing summation of last images now...");

		  /* take clean_stack median */
		  if(clean.top_stack > 0)	/* unused images on cosmics_stack */	
		    {
		      stat = Cleansum(N_AVERAGE);

		      /* do the elaborate sum */
		      stat = Elaborate(N_AVERAGE, KAPPA_SUM);

		      /* set countlevel to real sum of counts and do summation */
		      /* arg1 = 0 --> write values in output frame */
		      /* new frames in position: 0,1,...,top_stack-1 */
		      stat = Tosum(0,clean.top_stack);	

		      /* update LHCUTS */
		      stat = Statistics(sum.master.p_frame, sum.master.ima_number, CUT_MIN, CUT_MAX);

		      /* close this frame */
		      stat = SCTPUT("\nWriting out sum frame now:");
		      stat = SCTPUT(sum.master.name);

		      stat = SCFCLO(sum.master.ima_number);
	    
		    }


		  /* if clean.top_stack == 0 and master_sum is only saved at end */
		  /*(sum_save_flag = 1 or 2) */
		  else if(clean.top_stack == 0  && SUM_SAVE_FLAG >= 1)
		    {
		      /* arg1 = 0 --> frame is written to output */
		      /* new frames in position: 0,1,...,top_stack-1 */
		      stat = Tosum(0,0);	

		      /* update LHCUTS */
		      stat = Statistics(sum.master.p_frame, sum.master.ima_number, CUT_MIN, CUT_MAX);
	      
		      /* close this frame */
		      stat = SCTPUT("\nWriting out sum frame now:");
		      stat = SCTPUT(sum.master.name);

		      stat = SCFCLO(sum.master.ima_number);
		    }



		  /* give out warning if total number of images is smaller than n_average */
		  if(sum.total_opened < N_AVERAGE)
		    {	   
		      stat = SCTPUT("\nWARNING: Number of input images for summation is smaller than");
		      stat = SCTPUT("number of images used for median process. Cosmics removal did not work...");
		    }


		  /* write out last dirty_sum and difference */
		  if(SUM_SAVE_FLAG == 2)
		    stat = Savedirty();



		  /* --------- CUT SUM IMAGE ---------------------------------------------- */
		  /* cut sum images to actual size */
		  /* find appropriate size; flag = 0 */
		  stat = Cut(0 ,p_cutframe, &cut_x_start, &cut_x_end, &cut_y_start, &cut_y_end);

		  /* write dimensions */
		  cut_npix[0] = cut_x_end - cut_x_start + 1;	  
		  cut_npix[1] = cut_y_end - cut_y_start + 1;

		  /* create new name */
		  /* --> cut_lastsumname.fits */
		  strncpy(cutframe_name+4,sum.master.name,90);

		  /* create new cut frame */
		  stat = SCTPUT("\nNew frame for the cut sum image is created with name:");
		  stat = SCTPUT(cutframe_name);
		  sprintf(aux_string,"The pixel dimension of the new frame are: %d , %d",\
			  cut_npix[0], cut_npix[1]);
		  stat = SCTPUT(aux_string);

		  stat = SCIPUT(cutframe_name,D_R4_FORMAT,F_O_MODE,F_IMA_TYPE,2,cut_npix,start,step,\
				"cut and combined sum image of input catalog",cunit,\
				(char *)&p_cutframe,&cutframe_imno); 

		  /* open last master_sum image again for the descriptors */
		  stat = SCIGET(sum.master.name,D_R4_FORMAT,F_I_MODE,F_IMA_TYPE,2,&naxis,npix,start,step,\
				cutframe_ident,cunit,\
				(char *)&sum.master.p_frame,&sum.master.ima_number); 

		  /* copy descriptors from master_sum frame */
		  /* copy all other descriptors (non standard descr.) to output frame */
		  stat = SCDCOP(sum.master.ima_number,cutframe_imno,3,inframe_1);
		  /* copy descriptor LHCUTS */
		  stat = SCDCOP(sum.master.ima_number,cutframe_imno,4,"LHCUTS");

		  /* close master_sum frame */
		  stat = SCFCLO(sum.master.ima_number);

		  /* add descriptor CUT_SUBFRAME with subframe coordinates */
		  sprintf(aux_string,\
			  "The frame is a subframe of %s, with low-left and up-right coord: %d,%d ; %d.%d "\
			  , sum.master.name, cut_x_start, cut_y_start, cut_x_end, cut_y_end);
		  SCDWRC(cutframe_imno,"CUT_SUBFRAME",1,aux_string,1,200,&unit);

		  /* do copying of subframe; flag = 1 */
		  stat =  Cut(1 ,p_cutframe, &cut_x_start, &cut_x_end, &cut_y_start, &cut_y_end);

		  /* close this frame */
		  stat = SCTPUT("\nWriting out cut sum frame now:");
		  stat = SCTPUT(cutframe_name);

		  stat = SCFCLO(cutframe_imno);

		  /* ---------------------------------------------------------------------- */





		  /* display rejected images */
		  stat = SCTPUT("\n\nThe following images were rejected from summation: ");
		  
		  for(count=0 ; count < reject_count ; count++)
		    {
		      stat = SCTPUT(rejection[count].ima_name);
		      sprintf(aux_string,"x,y-pixel-offsets were to big or no objects were found: %f, %f\n", \
			     rejection[count].delta_x, rejection[count].delta_y );
		      stat = SCTPUT(aux_string);
		    }


		  stat = SCTPUT("\nSummation process finished...");
	 		  
		}
	      /* §§§§§§§§§§§ END SUMMATION IN REDUCTION §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ */



	      
	      /*.......... 3.10.20.60 Close down reduced frame .............................*/
	      sprintf(aux_string,"Image %s with image number %d is now closed...\n" \
		     ,stack_book[master_pos].frame_name, stack_book[master_pos].frame_imno);
	      stat = SCTPUT(aux_string);

	      stat = SCFCLO(stack_book[master_pos].frame_imno);


	      
	      /*.......... 3.10.20.70 Print out current loop parameters ....................*/
	      /* increment loop counter */
	      loop_counter++;				   	/* count number of reduced images */

	      /* test variables */
	      sprintf(aux_string,"\nThe variable values after the image loop are: \n");
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"Total number of images: %d", TOT_IMAGES);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"Top_Stack = %d", top_stack);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"Master_pos = %d", master_pos);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"New_pos = %d", new_pos);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"Old_pos = %d", old_pos);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"Old_number = %d", old_number);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"Last opened output frame has image number = %d", imno_outframe_1);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"Number of reduced images = %d\n\n", loop_counter);
	      stat = SCTPUT(aux_string);
	      sprintf(aux_string,"------------------- %d --------------------------------\n", loop_counter);
	      stat = SCTPUT(aux_string);


	      /*.......... 3.10.20.80 Organize handling of start/end images ................*/
	      /* set break for reduction of first SKY_FRAMES images at beginning */
	      if(loose_end == 0 && master_pos == SKY_FRAMES)
		{
		  stat = SCTPUT("\nFirst image block is reduced...");
		  loose_end = 1;
		  break;
		}

	      /* new master_pos for start_image reduction*/
	      if(loose_end == 0)
		master_pos++; 


	      /* handle last image block */
	      if(loose_end == 2)
		{
		  last_image++;

		  /* if last SKY_FRAMES images are reduced: break */
		  if(last_image == SKY_FRAMES)
		    {
		      stat = SCTPUT("\nLast image block is reduced...");
		      loose_end = 3;						/* close down pipeline */
		      break;
		    }
		}




	    }while(loose_end == 0 || loose_end == 2);		    		/* end of do-while loop */
	  /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/



	}	/* end of if(top_stack == TOT_IMAGES  &&  ACTION_FLAG != 1) */
      /*%%%%%%%%%%% END OF MAIN REDUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/





      /*.......... 3.10.60 Test break condition ..........................................*/
      /* break condition: loose_end == 3 */
      if(loose_end == 3)
	{
	  stat = SCTPUT("\nReduction of image catalog is done...");

	  /* break the endless image loop */
	  break;
	}


    }
  /*################# end of endless image loop #########################################*/






  /*------------- 3.99 Close Down Everything --------------------------------------------*/
  SCSEPI();
  exit(0);
}								/* END of main() */





/*+++++++++++++++ 4 FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*--------------- 4.1 Delay Function ----------------------------------------------------*/
/* delays programm for 'secs' seconds */

int Delay(int secs)
{
  time_t starttime;

  starttime = time(0);

do
  {
    /* do nothing */
  }while(time(0) <= starttime+secs);

 return 0;
}




/*--------------- 4.2 Copy Bad_Pixel_Mask -----------------------------------------------*/
/* copies the bad pixel mask to 2D frame "badpixel_mask" */
/* IN: pointer to mapped bad pixel mask */

int Copy_badpix(float *p_inframe)
{
  /* local count variables */
  int loop_x;
  int loop_y;
  int test = 0;
  int fstat;

  char aux_string[200];


  fstat = SCTPUT("\nCopying bad pixel mask...");

  /* do copying */
  for(loop_y = 0 ; loop_y < PIX_AXIS ; loop_y++)
    {
      for(loop_x = 0 ; loop_x < PIX_AXIS ; loop_x++)
	{
	  badpixel_mask[loop_y][loop_x] = *p_inframe;	/* copy value */

	  if((*p_inframe) > 0.5)
	    test++;	

	  p_inframe++;					/* increment pointer */
	}
    }

  sprintf(aux_string,"... %d bad pixels in bad pixel mask", test);
  fstat = SCTPUT(aux_string);


  return 0;
}





/*--------------- 4.3 Create Local Sensitivity Map --------------------------------------*/
/* copies the flatfield frame to 2D array sqrt_qe_map, normalizes it to one in  */
/* statistics window and takes squareroot for all pixels */
/* IN: pointer to flatfield frame */

int QE_map(float *p_flatfield)
{
  /* local count variables */
  int loop_x;
  int loop_y;
  int pix_no = 0;				/* number of pixels in statistics window */
  int fstat = 0;

  float mean_sum = 0;				/* sum up values */
  float mean = 0;					/* average in statistics window */


  fstat = SCTPUT("\nCreating local sensitivity map...");


  /* do copying */
  for(loop_y = 0 ; loop_y < PIX_AXIS ; loop_y++)
    {
      for(loop_x = 0 ; loop_x < PIX_AXIS ; loop_x++)
	{
	  sqrt_qe_map[loop_y][loop_x] = *p_flatfield;			 	/* copy value */
	
	  p_flatfield++;						  	/* increment pointer */
	}
    }


  /* find average in statistics window */
  for(loop_y = Y_START ; loop_y <= Y_END ; loop_y++)
    {
      for(loop_x = X_START ; loop_x <= X_END ; loop_x++)
	{
	  mean_sum += sqrt_qe_map[loop_y][loop_x];			 	/* sum value */
	  pix_no++;						  		/* increment counter */
	}
    }

  mean = mean_sum / pix_no;							/* average */


  /* do normalization to statistics window and take squareroot */
  for(loop_y = 0 ; loop_y < PIX_AXIS ; loop_y++)
    {
      for(loop_x = 0 ; loop_x < PIX_AXIS ; loop_x++)
	{
	  sqrt_qe_map[loop_y][loop_x] = (float) sqrt((double)(sqrt_qe_map[loop_y][loop_x] / mean));	   
	}
    }



  return 0;
}





/*--------------- 4.4 Create dummie Sensitivity Map --------_----------------------------*/
/* if no flatfield is available 2D array sqrt_qe_map with only 1's is created  */
/* IN: nothing */

int QE_one()
{
  /* local count variables */
  int loop_x;
  int loop_y;

  /* do copying */
  for(loop_y = 0 ; loop_y < PIX_AXIS ; loop_y++)
    {
      for(loop_x = 0 ; loop_x < PIX_AXIS ; loop_x++)
	{
	  sqrt_qe_map[loop_y][loop_x] = 1;			/* initialize with 1 */	
	}
    }

  return 0;
}
