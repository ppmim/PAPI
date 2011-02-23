/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		pipe_para.h
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		IR-Pipeline for OMEGA2000
.LANGUAGE   		C 
.ENVIRONMENT		MIDAS
.KEYWORDS		Parameters, Symbolic Constants
.PURPOSE		Define Pipeline related constants/parameters
.COMMENTS
.VERSION		1.00	12.08.02
			2.00 	22.10.02	modifications to use command line parameters
			3.00 	20.11.02	add summation process parameters
			3.1	12.02.03	add find object parameter
			3.2	26.03.03	added THRESH_MAX and ANA_RAD
			3.3	28.03.03	added DEBUG and DARK_SUBTRACTION
			3.4	08.04.03	added WAIT_MAX
			4.0	22.04.04	reduce memory requirements for other systems
-----------------------------------------------------------------------------------------*/


/* Define Detector: 0 = OMEGA2000 ; 1 = OMEGA_PRIME */
#define DETECTOR 0


/* Parameters for SINGLE IMAGE REDUCTION */
#define SKY_FRAMES_MAX 3  	/* maximun allowed number of frames before and after */
				/* main frame used for sky */
				/* # of stack levels = 2*SKY_FRAMES_MAX + 1 */


/* Parameters for SUMMATION PROCESS */
#define N_AVERAGE_MAX 7     	/* maximum allowed number of total frames for median process*/
				/* during image summation */

#define N_AVERAGE_DEFAULT 5 	/* default value for N_AVERAGE */


#define NORM_LEVEL 100000	/* count level that all images are normalized to */
				/* before summation  */	



/* Debugging parameter for additional function output */
#define DEBUG 0			/* 0=only normal output ; 1=additional output from OMEGA_pipe */
				/* 2=badpixel_correction ; 3=find_object*/


/* Turn dark current subtraction on and off */
#define DARK_SUBTRACTION 1	/* 0=off ; 1=do dark current subtraction */



/* TUNABLE PARAMETERS */
#define WAIT_MAX 10	      	/* maximal wait time for next image in online reduction mode */
				/* in units of 1/4*integration_time */


#define N_OBJECT 60		/* number of objects to be found for findobj-routine */
				/* has to be within [1,100] */

#define SEARCH_RADIUS 20	/* maximum search radius in pixels for the matching of objects */
				/* the object will be searched within this radius around */
				/* the expected position in steps of 5*/
				/* corresponds to the maximum tolerated pointing error of */
				/* the telescope */

#define SEARCH_STAR 2		/* pixel radius from intensity center of a cosmics candidate */
				/* were there is searched for star flux. If >COSMICS_CUT */
				/* pixels are found, value is taken as star value */


/* Parameters used in find_object.c to identify objects */
#define THRESH_MAX 3		/* maximal set threshold in units of the background */
#define ANA_RAD 8		/* pixel radius around intensity maximum in which star */
				/* is analyzed by subfunction Analyze*/
#define FWHM_MIN 0.5		/* factor for minimal fwhm selection:fwhm > fwhm_master * fact */
#define FWHM_MAX 1.5		/* factor for maximal fwhm selection */

#define INTENSITY_MIN 0.4		/* factor for object selection: object is accepted, if */
#define INTENSITY_MAX 2.5		/* total insity > MIN*masterframe_inten < MAX*masterint */
					/* if background changes, the scale factors change */
					/* --> should be set to allow for variations */

#define OFFSET_CUT 3			/* maximum pixel difference from median to be accepted */
					/* as the right object */
#define COSMICS_CUT 4			/* minimal number of pixels above background to be ac- */
					/* cepted as object. Used in find_object & elaborate_sum */
#define FWHM_CUT 4			/* maximal fwhm in arcsecs to pass the first object */
					/* selection ; objects > are classified as extended */
#define BACK_CUT 3			/* # of stddev above background to be accepted as object */
