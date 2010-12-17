/* skysub.c -- do sky subtraction and optional destriping */

/* 
 * input mask.fits is merged object mask and bpm.  the output bkg level is 
 * unchanged, it is the sky structure that is taken out.  should interpolate 
 * sky frame where sky weight is 0, but right now just assuming sky level
 * at those pixels equal to the mode of the sky frame.
 */

#include <stdio.h>
#include <stdlib.h>
#include "eprintf.h"
#include "stripe.h"
#include "fitsIO.h"
#include "hist.h"

/* sky subtraction and image destriping using object masking */
extern float *skysub(float *img, int nx, int ny, float bkg, float *bpm, 
                     float *sky, float *skyw, float *mask, char *type)
{
    int i;
    float skybkg, *imgout = (float *) emalloc(nx * ny * sizeof(float));

    skybkg = histcalcf(sky, nx, ny, -1, -1, NULL);
    
    printf ("\nSKYSUB ---skybkg== %f\n", skybkg);

    /* subtract out sky structure where sky image is valid (weight > 0) */

    for (i = 0; i < nx * ny; i++){
        if (skyw[i] > 0)
             imgout[i] = img[i] + (skybkg - sky[i]);  /* subt. sky structure */
	    /*imgout[i] = img[i] - sky[i]; jmiguel-tests */
        else
            imgout[i] = img[i];
        /*if (img[i]>220000)   
            printf("\nDEBUG> i=%d, img[i]=%f, skybkg=%f , sky[i]=%f IMGOUT=%f\n\n\n", i, img[i], skybkg, sky[i], imgout[i]);
        */   
	}
    /* do image destriping, mask indicates bad pix and object pix */

    destripe(imgout, mask, nx, ny, bkg, type);

    for (i = 0; i < nx * ny; i++)              /* set bad pixels to bkg lvl */
        if (bpm[i] <= 0)
            imgout[i] = bkg;

    return imgout;
}

/* sky subtraction and image destriping, no object masking */
extern float *skysub_nomask(float *img, int nx, int ny, float bkg, float *bpm,
                            float *sky, char *type)
{
    int i;
    float *imgout, skybkg;

    imgout = (float *) emalloc(nx * ny * sizeof(float));

    skybkg = histcalcf(sky, nx, ny, -1, -1, NULL);

    for (i = 0; i < nx * ny; i++)
        imgout[i] = img[i] + (skybkg - sky[i]);
	/* imgout[i] = img[i] - sky[i]; jmiguel-test */

    destripe(imgout, bpm, nx, ny, bkg, type);

    for (i = 0; i < nx * ny; i++)              /* set bad pixels to bkg lvl */
        if (bpm[i] <= 0)
            imgout[i] = bkg;

    return imgout;
}
