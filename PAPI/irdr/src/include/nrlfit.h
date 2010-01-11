/* nrlfit.h -- header file for nrlfit.c */

extern void lfit(double *x, double *y, double *sig, int ndat, double *a, 
	int *ia, int ma, double **covar, double *chisq, void (*funcs)(double,double *,int));
extern void covsrt(double **covar, int ma, int *ia, int mfit);
extern void gaussj(double **a, int n, double **b, int m);
extern void fpoly(double x, double *p, int np);
extern double *vector(int nl);
extern double **matrix(int nl, int nc);
