/* dark.c -- apply dark subtraction  */

/* image = (image - dark) */

/* filelist is an ASCII list containing inputfn outputfn */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

static void usage(void);

int main(int argc, char *argv[])
{
    int j, nx, ny;
    char line[256], infn[256], outfn[256];
    /*unsigned short *shortimg;*/
    float *dark, *img;
    FILE *fp;

    if (argc != 3)
        usage();

    dark = readfits(argv[2], &nx, &ny, NULL, NULL);


    if ((fp = fopen(argv[1], "r")) == NULL)
        eprintf("%s: failed opening: %s\n", argv[0], argv[1]);

    while (fgets(line, sizeof(line), fp) != NULL) {        /* read filelist */
        if (sscanf(line, "%s %s", infn, outfn) != 2)
            eprintf("%s: check list format\n", argv[0]);

        img = readfits(infn, &nx, &ny, NULL, NULL);

        printf("%s %d %d\n", infn, nx, ny);

        for (j = 0; j < nx * ny; j++)
            img[j] = (img[j] - dark[j]);

        /*shortimg = shortint(img, nx, ny);*/

        /* For PANIC, we need 32 bits images, so we write  -32 (float) FITS*/
        writefits(outfn, infn, (char*)img, -32, nx, ny);

        /*free(shortimg);*/  free(img);
    }

    return 0;
}

/* print out usage and exit */
static void usage(void)
{
    static char *usage = "\n"
    "dark - do dark subtraction \n\n"
    "usage: dark listfn darkfn \n\n"
    "where listfn - list of FITS images to dark subtract and output names:\n"
    "               input_FITS_filename output_FITS_filename\n"
    "      darkfn - filename of FITS dark image\n"
    "example: dark filelist dark.fits\n\n";
 
    printf("%s", usage);
    exit(0);
}
