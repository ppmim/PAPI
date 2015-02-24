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
    float skybkg, l_sky, *imgout = (float *) emalloc(nx * ny * sizeof(float));

    skybkg = histcalcf(sky, nx, ny, -1, -1, NULL);
    
    
    /*printf ("\nSKYSUB ---skybkg== %f\n", skybkg);*/

    /* subtract out sky structure where sky image is valid (weight > 0) */

    for (i = 0; i < nx * ny; i++)
    {
        if (bpm[i] <= 0)
            imgout[i] = bkg;  /* set bad pixels to bkg lvl */
            /* JMIM: I thik bad pixels should be set to __local__ bkg level to take into account when the pixel is in a star !! */
            /* compute local bkg, but probably only valid for isolated badpixels. !!*/
            /* TBC --> next probably only works fine for isolated pixels !!
            if ( (i-1) > 0 && (i+1) < nx * ny){
                l_sky = skybkg - (sky[i-1] + sky[i+1]) / 2.0 ;
                imgout[i] = (img[i-1] + img[i+1]) / 2.0 + l_sky ;  
            }
            else 
                imgout[i] = bkg;
            */
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
                
    	}
    }
    
    destripe(imgout, bpm, nx, ny, bkg, type);

    	
    return imgout;
}

