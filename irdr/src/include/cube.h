/* cube.h -- header file for cube.c */

#define MAXNPLANES 999         /* maximum number of images planes in cube */

extern float * cube_median_cl(float *planes[MAXNPLANES], int np, int nx, int ny,
                           float *scale, int offset);

extern float * cube_median(float *planes[MAXNPLANES], int np, int nx, int ny,
                           float *scale, int offset);

extern float * cube_sum(float *planes[MAXNPLANES], int np, int nx, int ny);
                           
extern float * cube_mean_nw(float *planes[MAXNPLANES], int np, int nx, int ny,
                            float *scale, int offset);

extern float * cube_mean_sw(float *planes[MAXNPLANES], float *weights, int np, 
                            int nx, int ny, float *scale, int offset);

extern float * cube_mean(float *planes[MAXNPLANES], float *wplanes[MAXNPLANES],
         int np, int nx, int ny, float **wplanesum, float *scale, int offset);

extern float *cube_sigma(float *planes[MAXNPLANES], int np, int nx, int ny,
                         float *scale, int offset);

extern float * cube_mean_min(float *planes[MAXNPLANES], int np, int nx, int ny,
                           float *scale, int offset, int N);

extern float * cube_mean_min_w(float *planes[MAXNPLANES], float *wplanes[MAXNPLANES],
         int np, int nx, int ny, float **wplanesum, float *scale, int offset, int N);

extern float * cube_median_min(float *planes[MAXNPLANES], int np, int nx, int ny,
                           float *scale, int offset, int N);

extern float * cube_mean_max(float *planes[MAXNPLANES], int np, int nx, int ny,
                           float *scale, int offset, int N);
