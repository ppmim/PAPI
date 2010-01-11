/* mean.h -- header file for mean.c */

extern float mean_nw(float *arr, int n, float nsig);

extern float mean(float *arr, float *w, int n, float nsig, float *sumw);

extern float stdev(float *arr, int n);
