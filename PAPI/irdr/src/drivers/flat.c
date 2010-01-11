/* flat.c -- apply flatfield correction */

/* filelist is an ASCII list containing inputfn outputfn */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

static void usage(void);

int main(int argc, char *argv[])
{
    int i, j, nx, ny;
    char line[256], infn[256], outfn[256];
    unsigned short *shortimg;
    float *flat, *img;
    FILE *fp;

    if (argc != 3)
        usage();

    flat = readfits(argv[2], &nx, &ny, NULL, NULL);

    for (i = 0; i < nx * ny; i++)       /* take inverse, multiply is faster */
        if (flat[i] > 0.0)
            flat[i] = 1.0 / flat[i];
        else
            flat[i] = 0.0;

    if ((fp = fopen(argv[1], "r")) == NULL)
        eprintf("%s: failed opening: %s\n", argv[0], argv[1]);

    while (fgets(line, sizeof(line), fp) != NULL) {        /* read filelist */
        if (sscanf(line, "%s %s", infn, outfn) != 2)
            eprintf("%s: check list format\n", argv[0]);

        img = readfits(infn, &nx, &ny, NULL, NULL);

        printf("%s %d %d\n", infn, nx, ny);

        for (j = 0; j < nx * ny; j++)
            img[j] = flat[j] * img[j];

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
    "flat - do flatfield correction\n\n"
    "usage: flat listfn gainfn\n\n"
    "where listfn - list of FITS images to flatten and output names:\n"
    "               input_FITS_filename output_FITS_filename\n"
    "      gainfn - filename of FITS gainmap (normalized flat field)\n\n"
    "example: flat filelist gain.fits\n\n";
 
    printf("%s", usage);
    exit(0);
}
