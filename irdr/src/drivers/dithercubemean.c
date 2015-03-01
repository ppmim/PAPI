/* dithercubemean.c -- coadd a stack of dithered FITS images */

/* 
 * Procedure:
 * Read FITS data planes and dither x,y pixel offsets listed in filelist
 * Generate corresponding weight image per data image
 * Calculate scale factors (zero offsets) to normalize image bkg levels
 * Calculate border size in pixels using dither offsets
 * Shift images and weight maps to register
 * Calculate the robust mean plane using weight maps and zero offsets
 *
 * Warning:
 * Need enough memory for all image planes and weight planes  
 * Expects images and gain map of same dimensions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "irdr.h"

static char *fn [MAXNPLANES];          /* FITS image file per image plane */
static float *data [MAXNPLANES];       /* pointers to image planes */
static float *wdata [MAXNPLANES];      /* pointers to weight map planes */
static float scale [MAXNPLANES];       /* scale factor per image plane */
static float xshift [MAXNPLANES];      /* image shift relative to first image */
static float yshift [MAXNPLANES];

static void usage(void);

int main(int argc, char *argv[])
{
    int i, nx, ny, nplanes, border;
    /* int xbelow=10, xabove=10, ybelow=10, yabove=10; */
    float *plane, *meanplane, *wplane, *sumwplanes, *gainmap;
    float bkg, sig = 0.0, avgscale = 0.0;
    int sum_flag = 0; /* if =1, compute the simple arithmetic sum of the planes */
    
    if (argc < 5)
        usage();
    
    if (argc == 6)
        /*read gainmap, altought it is not required !*/
        sum_flag = 1;
     

    nplanes = readlist(argv[1], fn, NULL, xshift, yshift, MAXNPLANES);

    if (nplanes < 1)
        eprintf("%s: no valid image planes\n", argv[0]);

    gainmap = readfits(argv[2], &nx, &ny, NULL, NULL);

    /* per togliere i bordi unitili
    border = new_get_border(xshift, yshift, nplanes, &xbelow, &xabove, &ybelow, &yabove);
    printf("Adding border of %d %d %d %d pixels\n", 
		    xbelow, xabove, ybelow, yabove);
    */

    border = get_border(xshift, yshift, nplanes);

    printf("Adding border of %d pixels\n", border);

    for (i = 0; i < nplanes; i++) {                  /* for each image plane */
        plane = readfits(fn[i], &nx, &ny, &bkg, &sig);

        if (bkg <= 0 || sig <= 0)
            eprintf("%s: bkg/sig <= 0 for %s BKG=%f, SIG=%f\n", argv[0], fn[i], bkg, sig);

        avgscale += (scale[i] = bkg);              /* store scale factors */

        wplane = getwmap(fn[i], nx, ny, gainmap, sig);

        data[i] = shift_image(plane, wplane, nx, ny, border, -xshift[i], 
                              -yshift[i], &wdata[i]); 

        /* per togliere i bordi inutili
	    * data[i] = new_shift_image(plane, wplane, nx, ny, 
			 xbelow, xabove, ybelow, yabove,
			 -xshift[i], -yshift[i], &wdata[i]);
	    */

        free(wplane);  free(plane);

        printf("%s %f %f %f %f\n", fn[i], bkg, sig, xshift[i], yshift[i]);
    }

    avgscale /= nplanes;

    for (i = 0; i < nplanes; i++)          /* use zero offsets to normalize */
        scale[i] = avgscale - scale[i];

    nx = nx + (2 * border);         /* image size of coadded dither set */
    ny = ny + (2 * border); 

    /* per togliere i bordi inutili
    nx = nx + xabove - xbelow;      * image size of coadded dither set *
    ny = ny + yabove - ybelow;      * image size of coadded dither set *
    */

    if (!sum_flag)
        meanplane = cube_mean(data, wdata, nplanes, nx, ny, &sumwplanes, scale, 1); /*ORIGINAL */
    else
        /* simple arithmetic sum of planes, no gainmap is used. */
        meanplane = cube_sum(data, nplanes, nx, ny);
   
    printf("Lets write the sum ...\n");
    
    writefits(argv[3], fn[0], (char*)meanplane, -32, nx, ny);
    
    if (sum_flag)
    {
        float i_time = 1.0;
        if (get_key_float(fn[0], "ITIME", &i_time) < 0) {
            fprintf(stderr, "get_wcs: unable to read ITIME in: %s\n", fn);
        }
        put_key_float(argv[3], "EXPTIME", nplanes*i_time);
    }
    put_key_int(argv[3], "NCOMBINE", nplanes);          /* update FITS hdr */
    put_key_float(argv[3], "DATAMODE", avgscale);

    /* write weight map */
    if (!sum_flag) writefits(argv[4], argv[3], (char*)sumwplanes, -32, nx, ny);
    
    return 0;
}

/* print out usage and exit */
static void usage(void)
{
    static char *usage = "\n"
    "dithercubemean - register dither frames and calculate mean plane\n\n"
    "usage: dithercubemean listfn gainfn outfn outwfn\n\n"
    "where listfn - contains list of: FITS_filename dither_x_off dither_y_off\n"
    "      gainfn - filename of FITS gainmap (normalized flat field)\n"
    "      outfn  - filename for output coadded FITS image\n"
    "      outwfn - filename for output coadded FITS weight image\n\n"
    "      sum    - (optional) flag to indicate to perform simple arithmetic sum of planes"
    "example: dithercubemean filelist gain.fits coadd.fits weight.fits\n\n";

    printf("%s", usage);
    exit(0);
}
