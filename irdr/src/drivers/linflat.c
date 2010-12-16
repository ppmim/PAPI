/*
 * linflat.c -- apply :
 * - non-linearity correction
 * - flatfield correction
 *
 * image = image * ( 1 + a1*image + a2*image^2 )
 * image = image / flat
 *
 * filelist is an ASCII list containing inputfn outputfn 
 * 
 * May 2002
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

static void usage(void);

int main(int argc, char *argv[])
{
    int i, j, nx, ny;
    char line[256], infn[256], outfn[256];
    unsigned short *shortimg;
    float *flat, *img, a1, a2, rap;
    FILE *fp;

    if (argc != 5)
        usage();

    a1 = atof(argv[2]);
    a2 = atof(argv[3]);
    flat = readfits(argv[4], &nx, &ny, NULL, NULL);

    /* take inverse of flat, faster to multiply */
    for (i = 0; i < nx * ny; i++)       
        if (flat[i] > 0.0)
            flat[i] = 1.0 / flat[i];
        else
            flat[i] = 0.0;

    if ((fp = fopen(argv[1], "r")) == NULL)
        eprintf("%s: failed opening: %s\n", argv[0], argv[1]);

    /* for each image of the filelist */
    while (fgets(line, sizeof(line), fp) != NULL) { 
        if (sscanf(line, "%s %s", infn, outfn) != 2)
            eprintf("%s: check list format\n", argv[0]);

        img = readfits(infn, &nx, &ny, NULL, NULL);

        printf("%s %d %d\n", infn, nx, ny);

        for (j = 0; j < nx * ny; j++) {
            /* linearity correction factor */
            rap = 1.0 + a1*img[j] + a2*img[j]*img[j];
            if(rap<1.0) rap=1.0;
            /* apply all corrections */
            img[j] = flat[j] * (img[j]*rap);
        }

        shortimg = shortint(img, nx, ny);

        writefits(outfn, infn, (char*)shortimg, 16, nx, ny);

        free(shortimg);  free(img);
    }

    return 0;
}

/* print out usage and exit */
static void usage(void)
{
    static char *usage = "\n"
    "linflat - do linearity and flatfield correction\n\n"
    "usage: darkflat listfn a1 a2 darkfn gainfn\n\n"
    "where listfn    - list of FITS images to flatten and output names:\n"
    "                  input_FITS_filename output_FITS_filename\n"
    "      a1        - first non-linearity correction coefficient\n"
    "      a2        - second non-linearity correction coefficient\n"
    "      gainfn    - filename of FITS gainmap (normalized flat field)\n\n";
 
    printf("%s", usage);
    exit(0);
}
