/* my_mode.c -- given a FITS file compute the mode */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "irdr.h"

static void usage(void);
#define MAXNBINS 400000

int main(int argc, char *argv[])
{
    int i, nx, ny;
    float min=+99999999;
    float max=-99999999;
    float *fimg;
    char fitsfilename[256];
    unsigned int *hist;
    long maxfreq=-1;
    int mode=-1;
    unsigned int *img = NULL;
    
    
    if (argc!=2)
        usage();
    
    hist = (unsigned int *) emalloc( MAXNBINS* sizeof(unsigned int));
    memset (hist, 0, MAXNBINS * sizeof(int));

    
    fimg = readfits(argv[1], &nx, &ny, NULL, NULL);
    img = longint(fimg, nx, ny);
    
    /* Fill the histogram */
    for (i=0; i < nx*ny; i++)
    {
        if ( img[i]>400000 || img[i]<=0 ) ;
        else
        {   
            hist[img[i]]++;
                       
            if (img[i] > max)
                max = img[i];

            if (img[i] < min)
                min = img[i];           
        }
    }
    
    /* compute the max. frequency (mode)*/
    maxfreq=-1;
    mode=-1;
    for (i=0; i <MAXNBINS ; i++)
    {
        if ( hist[i]> maxfreq ) {maxfreq=hist[i];mode=i;}
    }
    
    /* Show results */
    printf("\n MODE = %d ", mode);
    printf("\n MIN  = %f ", min);
    printf("\n MAX  = %f ", max);
    
    free(hist);
    
    return 0;
}

/* print usage and exit */
static void usage(void)
{
    static char *usage = "\n"
    "my_mode - compute image statistics mode \n\n"
    "usage: my_mode fitsfile "
    "example: my_mode gain.fits \n\n";

    printf("%s", usage);
    exit(0);
}
