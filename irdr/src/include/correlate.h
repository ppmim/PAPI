/* correlate.h -- header file for correlate.c */

#define MAXNCC 512     /* max dimension of cross-correlation image [pixels] */
/*#define MAXNCC 6000 */ /*jmiguel test */

extern float
correlate (int *xlist, int *ylist, float *pixlist, int nlist, float *img, 
   int nx, int ny, int ixoff, int iyoff, float *xoff, float *yoff, int maxerr);
