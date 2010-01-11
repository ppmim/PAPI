/*
 * lfit
 * from numerical recipies in C
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrlfit.h"


#define SQR(a) ((a)*(a))
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void lfit(x,y,sig,ndat,a,ia,ma,covar,chisq,funcs)
int ndat,ma,ia[];
double x[],y[],sig[],a[],**covar,*chisq;
void (*funcs)(double,double *,int); 
{
    int i,j,k,l,m,mfit=0;
    double ym,wt,sum,sig2i,**beta,*afunc;

    /* allocation of beta[ma][1] */
    beta=(double**) malloc(ma*sizeof(double*));
    if(!beta) {fprintf(stderr,"malloc failed\n"); exit(1);}
    for(i=0;i<ma;i++) {
        beta[i]=(double*) malloc(1*sizeof(double));
        if(!beta[i]) {fprintf(stderr,"malloc failed\n"); exit(1);}
    }

    /* allocation of afunc[ma] */
    afunc=(double*) malloc(ma*sizeof(double));
    if(!afunc) {fprintf(stderr,"malloc failed\n"); exit(1);}

    for(j=0;j<ma;j++)
        if (ia[j]) mfit++;
    if (mfit==0) {fprintf(stderr,"lfit: no parameters to be fitted\n"); exit(1);}
    /* Initialise the (symmetric) matrix */
    for(j=0;j<mfit;j++) {
        for(k=0;k<mfit;k++) covar[j][k]=0.0;
        beta[j][0]=0.0;
    }
    /* Loop over data to accumulate coefficients of the normal equations */
    for(i=0;i<ndat;i++) {
        (*funcs)(x[i],afunc,ma);
        ym=y[i];
        /* subtract off dependences on known pieces of the fitting function */
        if (mfit<ma) {
            for(j=0;j<ma;j++)
                if (!ia[j]) ym -= a[j]*afunc[j];
        }
        sig2i=1.0/SQR(sig[i]);
        for (j=-1,l=0;l<ma;l++) {
            if (ia[l]) {
                wt=afunc[l]*sig2i;
                for (j++,k=-1,m=0;m<=l;m++)
                    if (ia[m]) covar[j][++k] += wt*afunc[m];
                beta[j][0] += ym*wt;
            }
        }
    }
    /* Fill in above the diagonal from symmetry */
    for(j=1;j<mfit;j++)
        for(k=0;k<j;k++)
            covar[k][j]=covar[j][k];
    /* Matrix solution */
    gaussj(covar,mfit,beta,1);
    for(j=-1,l=0;l<ma;l++)
        if (ia[l]) a[l]=beta[++j][0];
    /* Evaluate chisq of the fit */
    *chisq=0.0;
    for(i=0;i<ndat;i++) {
        (*funcs)(x[i],afunc,ma);
        for(sum=0.0,j=0;j<ma;j++) sum += a[j]*afunc[j];
        *chisq += SQR((y[i]-sum)/sig[i]);
    }
    /* Sort covariance matrix to true order of fitting coefficients */    
    covsrt(covar,ma,ia,mfit);

    /* free memory */    
    free(afunc);
    for(i=0;i<ma;i++) free(beta[i]);
    free(beta);
}

void covsrt(covar,ma,ia,mfit)
double **covar;
int ma,ia[],mfit;
{
    int i,j,k;
    double swap;

    for(i=mfit;i<ma;i++)
        for(j=0;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
    k=mfit;
    for(j=ma-1;j>=0;j--) {
        if (ia[j]) {
            for(i=0;i<ma;i++) SWAP(covar[i][k],covar[i][j]);
            for(i=0;i<ma;i++) SWAP(covar[k][i],covar[j][i]);
            k--;
        }
    }
}

void gaussj(a,n,b,m)
double **a,**b;
int n,m;
{
    int *indxc,*indxr,*ipiv;
    int i,icol,irow,j,k,l,ll;
    double big,dum,pivinv,swap;

    indxc=(int*)malloc(n*sizeof(int));
    indxr=(int*)malloc(n*sizeof(int));
    ipiv=(int*)malloc(n*sizeof(int));
    if((!indxc)||(!indxr)||(!ipiv)) {
        fprintf(stderr,"malloc failed\n"); exit(1);
    }

    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if (ipiv[j] != 1)
                for (k=0;k<n;k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big=fabs(a[j][k]);
                            irow=j;
                            icol=k;
                        }
                    } 
                }
        ++(ipiv[icol]);
        if (irow != icol) {
                for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
                for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
        }
        indxr[i]=irow;
        indxc[i]=icol;
        if (a[icol][icol] == 0.0) {
            fprintf(stderr,"gaussj: Singular Matrix\n"); exit(1);
        }
        pivinv=1.0/a[icol][icol];
        a[icol][icol]=1.0;
        for (l=0;l<n;l++) a[icol][l] *= pivinv;
        for (l=0;l<m;l++) b[icol][l] *= pivinv;
        for (ll=0;ll<n;ll++)
                if (ll != icol) {
                        dum=a[ll][icol];
                        a[ll][icol]=0.0;
                        for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
                        for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
                }
    }
    for (l=n-1;l>=0;l--) {
            if (indxr[l] != indxc[l])
                    for (k=0;k<n;k++)
                            SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
    free(ipiv);
    free(indxr);
    free(indxc);
}

/* 
 * Fitting routine for a polynomial of degree np-1, 
 * with coefficients in the array p[0..np-1]
 */
void fpoly(double x, double p[], int np)
{
    int j;
    p[0]=1.0;
    for(j=1;j<np;j++) p[j]=p[j-1]*x;
}

double *vector(int nl)
{
    double *v;
    v=(double*)malloc(nl*sizeof(double));
    if (!v) { fprintf(stderr,"malloc failed\n"); exit(1); }
    return v;
}

double **matrix(int nl, int nc)
{
    double **m;
    int i;
    m=(double**)malloc(nl*sizeof(double));
    if (!m) { fprintf(stderr,"malloc failed\n"); exit(1); }
    for(i=0;i<nl;i++) {
        m[i]=(double*)malloc(nc*sizeof(double));
        if (!m[i]) { fprintf(stderr,"malloc failed\n"); exit(1); }
    }
    return m;
}

#undef SWAP
#undef SQR
