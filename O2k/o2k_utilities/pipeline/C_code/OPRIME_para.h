/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		OPRIME_para.h
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE	   	IR-Pipeline for OMEGA2000
.LANGUAGE   		C 
.ENVIRONMENT		MIDAS
.KEYWORDS		Parameters, Symbolic Constants
.PURPOSE		Define OMEGA PRIME detector related constants
.COMMENTS
.VERSION		1.00	12.08.02
--------------------------------------------------------------------------------------*/


/* OMEGA PRIME Detector Constants */

#define TOT_PIX 1048576		 	/* =1024*1024 TOTal number of PIXels on detector*/
#define PIX_AXIS 1024			/* PIXels per AXIS = 1024 */


#define PIX_SIZE 18.5			/* PIXel SIZE of detector is 18.5 microns  */
#define PIX_SCALE 0.3961		/* PIXel SCALE is 0.4 arcsec/pix */
#define FIELD_OF_VIEW 6.8		/* FIELD OF VIEW is 6.8*6.8 arcmin */
#define FOCAL_RATIO 2.6			/* ? FOCAL RATIO at prime focus is 2.6  */


#define DELTA_REJECTION	 2		/* arc_minutes beyond which incoming frames are rejected */
					/* from summation process: in RA and DEC */



/* specified frame area to do statistcs in*/

#define X_START 350			/* x-start-coordinate in pixels */
#define X_END 450			/* x-end-coordinate in pixels */
#define Y_START 350			/* y-start-coordinate in pixels */
#define Y_END 450			/* y-end-coordinate in pixels */
