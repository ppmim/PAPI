/* median.h -- header file for median.c */

extern float median(const float *arr, int n);

extern float median_double(const double *arr, int n);

extern float sample_median(const float *img, float *sigma, int n);

extern float median_absdev(const float *arr, float m, int n);

extern double median_absdev_double(const double *arr, float m, int n);
