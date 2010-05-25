/* hist.c -- fast, iterative histogram analysis for u_short & int data */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include "irdr.h"

#define MAXNITER  5
#define MAXNBINS  400000         /* modified by jmiguel ~April-2009 to manage to works with O2000 data files */
#define NSIG      5.0
#define TOL       0.01           /* convergence criterion to stop iteration */
#define BLOCKSIZE 256           /* default image block size for processing */

static unsigned int *bins, *hist;      /* counts in each histogram bin */
static unsigned int *cbins;           /* cumulative counts up to each bin */
static float *wcbins;   /* cumulative counts weighted by bin number */
static int minbin, maxbin, nbins;

static float histclip(int lcut, int hcut, float *med, float *sig);
static void histaccum(void);
static void histinit(void);
static void histfill(unsigned int *a, int n);

/* histinit: initialize histogram */
static void histinit(void)
{
    minbin = MAXNBINS;
    maxbin = 0;
    memset(hist, 0, MAXNBINS * sizeof(int));

}

/* histfill: increment histogram */
static void histfill(unsigned int *a, int n)
{
    int i;

    if (n < 1)
        return;

    for (i = 0; i < n; i++) {                       /* fill histogram */
        if (a[i]<0 ) printf("\nDEBUG- Valor negativo !!!!  no es posible ????");
        else if (a[i]>=400000);
        else 
        {
            hist[a[i]]++;

            if (a[i] > maxbin)
                maxbin = a[i];

            if (a[i] < minbin)
                minbin = a[i];
        }
    }

    nbins = maxbin - minbin + 1;
    bins = hist + minbin;
    


}

/* histmode: robust histogram mode from iterative k-sigma clipping */
extern float histmode(float *sigma)
{
    int i, lcut = 0, hcut = nbins - 1;
    float avg = 0.0, med, sig, prevsig = 9e9;

    histaccum();

    for (i = 0; i < MAXNITER; i++) {                /* iterative clipping */
        avg = histclip(lcut, hcut, &med, &sig);
        
        if (sig < 0.1 || fabs(sig / prevsig - 1.0) < TOL)
            break;

        lcut = max((int)(med - NSIG * sig - 0.5), 0);
        hcut = min((int)(med + NSIG * sig + 0.5), nbins - 1);
        prevsig = sig;
    }

    *sigma = sig;
    
    /*if (med<=0) printf("\nDEBUG- AQUI ESTA EL PROBLEMA (med=0)=> mode<0   !!!!! avg=%f,  med=%f LCUT=%d, HCUT=%d, nbins=%d\n\n ", avg, med,lcut, hcut, nbins);*/
    
    /* jmim 13-04-2009 */
    /*if (med<=0) return 0;
    else return (avg > med) ? (3.0 * med - 2.0 * avg) : avg; */   /* mode approx. */
    /* end jmim */
    
    
    return (avg > med) ? (3.0 * med - 2.0 * avg) : avg;   /* mode approx. */
}

/* histclip: get statistics of the clipped histogram (bins lcut, hcut) */
static float histclip(int lcut, int hcut, float *med, float *sig)
{
    int dn, dn0;
    double q25, wdn, wdn0;
    int i;

    dn0 = (lcut > 0) ? cbins[lcut-1] : 0;     /* total counts below bin lcut */
    wdn0 = (lcut > 0) ? wcbins[lcut-1] : 0.0;

    dn = cbins[hcut] - dn0;               /* total counts, bins lcut to hcut */
    wdn = wcbins[hcut] - wdn0;

    q25 = bisearch(dn0 + dn/4, cbins, nbins);         /* binary search */
    *med = bisearch(dn0 + dn/2, cbins, nbins);
    /*DEBUG if (*med==0)
    {
        printf("\nDEBUG- dn0=%d, dn=%d, nbins=%d  cbins[%d]=%d\n", dn0, dn, nbins, hcut, cbins[hcut]);
        for (i=0; i<100; i++) printf("cbins[%d]=%d", i, cbins[i]);
    }
    */
    *sig = (*med - q25) / 0.6745;

    return (dn > 0) ? (wdn / dn) : lcut;              /* return mean */
}

/* histaccum: create cumulative histograms */
static void histaccum(void)
{
    int i;

    cbins[0] = bins[0];

    wcbins[0] = 0;

    for (i = 1; i < nbins; i++) {
        cbins[i] = cbins[i-1] + bins[i];
        /*printf("\nDEBUG-  cbins[%d]=%d", i, cbins[i]);*/
        wcbins[i] = (float)i * bins[i] + wcbins[i-1];
        /*printf("\nI=%d", i);*/
    }
}

/* histcalc: calculate image mode and sigma using histogram analysis */
extern float histcalc(unsigned int *img, int nx, int ny, int nxb, int nyb,
                      float *sig)
{
    int i, j, k, nb, n = 0;
    float mode, *modes, *sigmas;


    cbins = (unsigned int *) emalloc( MAXNBINS* sizeof(unsigned int));
    hist = (unsigned int *) emalloc( MAXNBINS* sizeof(unsigned int));
    wcbins = (float *) emalloc( MAXNBINS* sizeof(float));

    if (nxb < 1 || nyb < 1)                 /* use default block size */
        nxb = nyb = BLOCKSIZE;

    nxb = min(nxb, nx);
    nyb = min(nyb, ny);
    nb  = nx / nxb * ny / nyb;               /* number of image blocks */

    modes  = (float *) emalloc(nb * sizeof(float));
    sigmas = (float *) emalloc(nb * sizeof(float));
    
    for (i = 0; (i+nyb) <= ny; i += nyb) {
        for (j = 0; (j+nxb) <= nx; j += nxb) {       /* foreach image block */
            histinit();

            for (k = i; k < i + nyb; k++)
              histfill(img + k * nx + j, nxb);         /* fill histogram */

            modes[n] = minbin + histmode(&sigmas[n]);
            /*if (modes[n]<0) modes[n]=0; 
            printf("\nDEBUG- mode[%d]= %f", n, modes[n]);*/
            n++;
        }
    }
    
        
    mode = mean_nw(modes, n, 5.0);
    /*printf ("\nDEBUG- MODE=%f\n", mode);
    if (mode <0 ) exit(0); DEBUG */

    if (sig != NULL)
        *sig = mean_nw(sigmas, n, 5.0);

    free(modes);  free(sigmas); free(cbins); free(hist); free(wcbins);

    return mode;
}

/* histcalcf: round float image to u_short and call histcalc() */
extern float histcalcf(float *fimg, int nx, int ny, int nxb, int nyb,
                       float *sig)
{
    float mode;
    unsigned int *img = longint(fimg, nx, ny);
    int i=0;
    
    /*for ( i=0; i<10; i++)
        printf("\nIMG[%d]=%d", i, img[i]);

    */
    mode = histcalc(img, nx, ny, nxb, nyb, sig);

    free(img);

    return mode;
}

/* histcalca: histogram analysis for float array */
extern float histcalca(float *a, int n)
{
    int i;
    float bkg, sig;

    bkg = sample_median(a, &sig, n);

    minbin = max(0, (bkg - 10 * sig));
    maxbin = min((MAXNBINS - 1), (bkg + 10 * sig));
    nbins  = maxbin - minbin + 1;
    bins   = hist + minbin;

    memset(bins, 0, nbins * sizeof(int));       /* initialize histogram */

    for (i = 0; i < n; i++)                        /* fill histogram */
        hist[(unsigned int)(a[i] + 0.5)]++;        /* round */

    return minbin + histmode(&sig);
}
