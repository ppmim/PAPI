/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		O2000_para.h
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE		IR-Pipeline for OMEGA2000
.LANGUAGE   		C 
.ENVIRONMENT		MIDAS
.KEYWORDS		Parameters, Symbolic Constants
.PURPOSE		Define OMEGA2000 detector related constants
.COMMENTS
.VERSION		1.00	12.08.02
--------------------------------------------------------------------------------------*/


/* OMEGA2000 Detector Constants */

#define TOT_PIX 4194304			/* =2048*2048 TOTal number of PIXels on detector*/
#define PIX_AXIS 2048			/* PIXels per AXIS = 2048 */


#define PIX_SIZE 18			/* PIXel SIZE of detector is 18 microns  */
#define PIX_SCALE 0.45			/* PIXel SCALE is 0.45 arcsec/pix */
#define FIELD_OF_VIEW 15.4		/* FIELD OF VIEW is 15.4*15.4 arcmin */
#define FOCAL_RATIO 2.35		/* FOCAL RATIO at prime focus is 2.35 */


#define DELTA_REJECTION	 3		/* arc_minutes beyond which incoming frames are rejected */
					/* from summation process: in RA and DEC */



/* specified frame area to do statistcs in*/

#define X_START 700			/* x-start-coordinate in pixels */
#define X_END 900			/* x-end-coordinate in pixels */
#define Y_START 700			/* y-start-coordinate in pixels */
#define Y_END 900			/* y-end-coordinate in pixels */
