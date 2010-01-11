/* cube.h -- header file for cube.c */

#define MAXNPLANES 999         /* maximum number of images planes in cube */

extern float * cube_median(float *planes[MAXNPLANES], int np, int nx, int ny,
                           float *scale, int offset);

extern float * cube_mean_nw(float *planes[MAXNPLANES], int np, int nx, int ny,
                            float *scale, int offset);

extern float * cube_mean_sw(float *planes[MAXNPLANES], float *weights, int np, 
                            int nx, int ny, float *scale, int offset);

extern float * cube_mean(float *planes[MAXNPLANES], float *wplanes[MAXNPLANES],
         int np, int nx, int ny, float **wplanesum, float *scale, int offset);

extern float *cube_sigma(float *planes[MAXNPLANES], int np, int nx, int ny,
                         float *scale, int offset);
