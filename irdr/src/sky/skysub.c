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
    
    /*printf ("\nSKYSUB ---skybkg== %f\n", skybkg);*/

    /* subtract out sky structure where sky image is valid (weight > 0) */

    for (i = 0; i < nx * ny; i++)
    {
        if (bpm[i] <= 0)
            imgout[i] = bkg;  /* set bad pixels to bkg lvl */
        else
        {    
            if (skyw[i] > 0)
                imgout[i] = img[i] + (skybkg - sky[i]);  /* subt. sky structure and add constant (skybkg) to preserve original count level */
            else
                imgout[i] = img[i];
        }
	}
	
    /* do image destriping, mask indicates bad pix and object pix */

    destripe(imgout, mask, nx, ny, bkg, type);

    /*
    for (i = 0; i < nx * ny; i++)              
        if (bpm[i] <= 0)
            imgout[i] = bkg;
    */
    return imgout;
}

/* sky subtraction and image destriping, no object masking */
/* NOTE: in order to preserve the original count level, a constant (mode)
   representing the median sky level is added to all pixels.
*/
extern float *skysub_nomask(float *img, int nx, int ny, float bkg, float *bpm,
                            float *sky, char *type)
{
    int i;
    float *imgout, skybkg;

    imgout = (float *) emalloc(nx * ny * sizeof(float));

    skybkg = histcalcf(sky, nx, ny, -1, -1, NULL);

    for (i = 0; i < nx * ny; i++)
    {	
    	if (bpm[i] <= 0)
            imgout[i] = bkg;  /* set bad pixels to bkg lvl */
    	
    	else{
            imgout[i] = img[i] + (skybkg - sky[i]); /* add constant (skybkg) to preserve original count level */
                
    	    /* only for debug !
    	    if (skybkg-sky[i]<0)
    	       printf("\n IMG= %f  skbkg= %f SKY= %f (d=%f)", img[i],  skybkg, sky[i], skybkg-sky[i]);
    	    */
    	}
    }
    
    destripe(imgout, bpm, nx, ny, bkg, type);

    	
    return imgout;
}

