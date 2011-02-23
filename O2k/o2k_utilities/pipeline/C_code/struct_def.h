/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.COPYRIGHT 		(C) MPIA 
.IDENTIFIER		struct_def.h
.AUTHOR    		Rene Fassbender , MPIA - Heidelberg
.MODULE	   	IR-Pipeline for OMEGA2000
.LANGUAGE   		C 
.ENVIRONMENT		MIDAS
.KEYWORDS		Structure Definitions
.PURPOSE		Define Structures to make them accessible to all modules
.COMMENTS
.VERSION		1.00	04.09.02
			2.0   	22.11.02	update for summation process
			2.1	29.11.02	changed order of arrays to avoid errors
			2.2	03.12.02	increased info stack size
			2.3	04.12.02	take clean.cosmic_stack out of structure
			2.4	22.12.02	add find_object structure
			2.5	11.02.02	only conditional include for pipe_para
			2.6	07.04.03	added original offsets in SUM_PAGE
			2.7	05.06.03	changed max sum images from 200 to 1000
-----------------------------------------------------------------------------------------*/

/* include headers */
/* Parameter and Symbolic Constant Files */
#if !defined( DETECTOR )
#include "pipe_para.h"		/* pipeline related parameters/constants */
#endif


/* Parameter File of specified Detector */
#if DETECTOR==1
#include "OPRIME_para.h"	
#elif DETECTOR==0
#include "O2000_para.h"
#endif





/*-------------------- Structure Definitions --------------------------------------------*/
/* Secundary modules that use these structures need the this header file */

/* SINGLE IMAGE REDUCTION */
/* stucture with synonym BOOKPAGE for the bookkeeping of the image stack */
typedef struct {
  int frame_imno;					/* image number of frame */
  float *p_frame;					/* pointer to frame */
  char frame_name[85];					/* name of frame */

  /* the next part is used only if sky frames are saved */
  int sky_imno;						/* image number of corresponding sky frame */
  float *p_sky;						/* pointer to corresponding sky frame */
  char sky_name[85];					/* name of sky frame */
} BOOKPAGE;						/* synonym for structure */



/* SUMMATION PROCESS */
/* structure for bookkeeping of images */
typedef struct {

  char ima_name[85];					/* image name */

  int ima_number;					/* image number */

  double alpha;						/* RA in degrees */
  double declination;				     	/* declination in degrees */

  float delta_x;					/* pixel_offset in x compared to first ima */
  float delta_y;					/* pixel_offset in y */

  float org_delta_x;					/* original pixel_offset in x */
  float org_delta_y;					/* original pixel_offset in y */


  float countlevel;					/* median countlevel of frame */

  float seeing;						/* approximate seeing in frame */

} SUM_PAGE;


/* structure for output frames handling */
typedef struct {

  char name[90];					/* name of frame */

  float *p_frame;					/* pointer to mapped data */

  int ima_number;					/* image number of frame */

} FRAME_INFO;


/* structure for the sum-image and related information */
typedef struct {

  /* bookkeeping */					
  int total_sum;					/* # of summed images */
  int total_opened;

  float sum_level;					/* sum of countlevels */

  SUM_PAGE reference;					/* info on reference frame */

  /* info on designated output frames */
  FRAME_INFO master;					/* master output frame */
  FRAME_INFO dirty;					/* output frame for real sum */
  FRAME_INFO difference;				/* output frame for difference */

  SUM_PAGE info[1000];					/* 1000 info pages */

  /* images */
  float master_sum[PIX_AXIS][PIX_AXIS];		     	/* master image for summation */

  float dirty_sum[PIX_AXIS][PIX_AXIS];			/* frame for real sum */
							/*(without cosmics removal) */
  float scratch_sum[PIX_AXIS][PIX_AXIS];		/* buffer array for intermediate storage */

  /* overflow test */
  int overflow;						/* to test overflow of arrays */

} DEEP_STRUCT;


/* structure for the cosmics removal process */
typedef struct {

  /* bookkeeping */
  int top_stack;					/* keep track of top level */

  SUM_PAGE info[N_AVERAGE_MAX];

  /* image stack: images stored in a series, instead of pixels in series */
  /* float cosmics_stack[N_AVERAGE_MAX][PIX_AXIS][PIX_AXIS]; */

} COSMICS_STACK;



/* structure for storing object information */
typedef struct {

  /* position */
  float x_pos;
  float y_pos;

  /* FWHM */
  float fwhm;

  /* central intensity */
  float intensity;

} OBJECT_INFO;
