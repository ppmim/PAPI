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

  float countlevel;					/* median countlevel of frame */

} SUM_PAGE;


/* structure for output frames handling */
typedef struct {

  char name[90];					/* name of frame */

  float *p_frame;					/* pointer to mapped data */

  int ima_number;					/* image number of frame */

} FRAME_INFO;


/* structure for the sum-image and related information */
typedef struct {
  /* images */
  float master_sum[PIX_AXIS][PIX_AXIS];		     	/* master image for summation */

  float dirty_sum[PIX_AXIS][PIX_AXIS];			/* frame for real sum */
							/*(without cosmics removal) */
  float scratch_sum[PIX_AXIS][PIX_AXIS];		/* buffer array for intermediate storage */

  /* info on designated output frames */
  FRAME_INFO master;					/* master output frame */
  FRAME_INFO dirty;					/* output frame for real sum */
  FRAME_INFO difference;				/* output frame for difference */

  /* bookkeeping */					
  SUM_PAGE info[100];					/* 100 info pages */

  SUM_PAGE reference;					/* info on reference frame */

  float sum_level;					/* sum of countlevels */

  int total_sum;					/* # of summed images */
  int total_opened;

} DEEP_STRUCT;


/* structure for the cosmics removal process */
typedef struct {
  /* image stack: images stored in a series, instead of pixels in series */
  float cosmics_stack[N_AVERAGE_MAX][PIX_AXIS][PIX_AXIS];

  /* bookkeeping */
  SUM_PAGE info[N_AVERAGE_MAX];

  int top_stack;					/* keep track of top level */

} COSMICS_STACK;
