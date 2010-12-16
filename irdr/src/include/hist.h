/* hist.h -- header file for hist.c */

extern float histmode(float *sigma);

extern float 
histcalc(unsigned int *img, int nx, int ny, int nxb, int nyb, float *sig);

extern float 
histcalcf(float *fimg, int nx, int ny, int nxb, int nyb, float *sig);

extern float histcalca(float *a, int n);
