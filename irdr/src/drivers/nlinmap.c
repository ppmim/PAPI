/* 
 * - Computes the non-linearity coefficients a1 and a2 :
 * Fcorr = Fm * ( 1 + a1*Fm + a2*fm^2 )
 * Don't take points over satlim into consideration
 * Use exptime_ref as flux reference : runs are supposed to be taken with an
 * exptime_ref exposure between other exposures for calibration
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "irdr.h"
#include "nrlfit.h"

static void usage(void);
static char *fn[MAXNPLANES]; 
#define ma 3   /* number of coefficients for the regression (a0+a1*F+a2*F^2)*/
#define mdat 50

/***********************/
/* PROGRAMME PRINCIPAL */
/***********************/

int main (int argc, char **argv)
{
   double **covar, chisq, *fm, *feOfm, *sig, *a;
   char *listfn, *outdatfn, *nlincoeffn;
   int nx,ny,ndat,ia[ma],i,k,m,nplanes, kdat;
   float satlim, refexptime, expT, expTmem;
   float fref;
   float Fref1, Fref2, Fmem, F, *img;
   float fmTab[mdat], feOfmTab[mdat]; 
   FILE *OUT;

   /** read arguments **/
   if (argc!=6) usage();
   listfn = argv[1];
   satlim = atof(argv[2]);
   refexptime = atof(argv[3]);
   outdatfn = argv[4];
   nlincoeffn = argv[5];
  
   /** image ratios inititalisations **/ 
   if ((nplanes = readlist(listfn,fn,NULL,NULL,NULL, MAXNPLANES))<1)
       eprintf("%s: no valid image planes",argv[0]);

   /** compute image ratios **/
   Fmem=0; Fref1=0; Fref2=0;
   kdat=0;
   for(k=0;k<nplanes;k++) {
       img = readfits(fn[k],&nx,&ny,NULL,NULL);
       F = mean_nw(img, nx*ny, 3);
       free(img);
       if (get_key_float(fn[k],"EXPTIME",&expT)<0) 
           eprintf("nlinmap: %s doesn't have EXPTIME fitskey",fn[k]);
         if(expT==refexptime) {
           Fref1 = Fref2;
           Fref2 = F;
           if((Fmem!=0)&&(Fref1*Fref2!=0)) {
              fref = (Fref1+Fref2)/2.;
              fmTab[kdat] = Fmem;
              feOfmTab[kdat] = fref/refexptime/Fmem*expTmem;
              kdat++;
              Fmem=0;
           }
        }
        else {
            Fmem = F;
            expTmem = expT;
        }
   }
   if (Fmem!=0) {  /* last one to treat */
           fref = Fref2;
           fmTab[kdat] = Fmem;
           feOfmTab[kdat] = fref/refexptime/Fmem*expTmem;
           kdat++;
   }

   if(kdat>mdat) 
       eprintf("nlinmap: too many data points. change mdat value");


   if ((OUT = fopen(outdatfn, "w")) == NULL)
       eprintf("nlinmap: can't open %s\n",outdatfn);
   for(k=0;k<kdat;k++) { 
           fprintf(OUT,"%g  %g\n",fmTab[k],feOfmTab[k]);
   }
   fclose(OUT);


   /** lfit initialisations **/
   a = vector(ma);
   ia[0]=0; /* forcing Fc=Fm for no flux */
   for (m=1;m<ma;m++) ia[m]=1;
   a[0]=1.0;

   fm = vector(kdat);
   feOfm = vector(kdat);
   sig = vector(kdat);
   covar = matrix(kdat,kdat);
   for(i=0;i<kdat;i++) sig[i]=1.0;

   /** compute linearity coefficients */
   ndat = 0;
   for(k=0;k<kdat;k++) {
      if (fmTab[k]<satlim) {
           fm[ndat] = (double) fmTab[k];
           feOfm[ndat] = (double) feOfmTab[k];
           ndat++; 
       }
   }
   lfit(fm, feOfm, sig, ndat, a, ia, ma, covar, &chisq, fpoly);
   printf("nlinmap: a1 = %g\n",a[1]);
   printf("nlinmap: a2 = %g\n",a[2]);

   if((OUT = fopen(nlincoeffn,"a")) == NULL)
       eprintf("nlinmap: can't open %s\n",nlincoeffn);
   fprintf(OUT,"%g %g\n",a[1],a[2]);
   fclose(OUT);
           
   return 0;
}

/* print out usage and exit */
static void usage(void)
{
    static char *usage = "\n"
    "nlinmap - compute non-linearity coefficients\n\n"
    "usage: nlinmap listfn satlim refexptime outdatfn nlincoeffn\n"
    "where\tlistfn - contains list of FITS images\n"
        "\tsatlim - saturation limit\n"
        "\trefexptime - exptime used as reference\n"
        "\toutdatfn - filename of out data file, contains (fm, feOfm)\n"
        "\tlincoeffn - filename where coeffs will be added\n\n"
    "example: nlinmap filelist 40000 4.0 feOfm.c1.dat nlincoeff.dat\n\n";

    printf("%s", usage);
    exit(0);
}



