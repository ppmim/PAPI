/* imcalc.c -- image1 op image2 */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int i, nx, ny;
    float *img1, *img2, *img3;

    if (argc != 5)
        eprintf("Usage: %s a.fits op b.fits out.fits\n", argv[0]);

    img1 = readfits(argv[1], &nx, &ny, NULL, NULL);
    img2 = readfits(argv[3], &nx, &ny, NULL, NULL);
    img3 = (float *) emalloc(nx * ny * sizeof(float));

    if (argv[2][0] == '+') {
        for (i = 0; i < nx * ny; i++)
            img3[i] = img1[i] + img2[i];

    } else if (argv[2][0] == '-') {
        for (i = 0; i < nx * ny; i++)
            img3[i] = img1[i] - img2[i];

    } else if (argv[2][0] == '*') {
        for (i = 0; i < nx * ny; i++)
            img3[i] = img1[i] * img2[i];

    } else if (argv[2][0] == '/') {
        for (i = 0; i < nx * ny; i++) {
            if (img2[i] != 0)
                img3[i] = img1[i] / img2[i];
            else
                img3[i] = 1;
        }          

    } else {
        eprintf("%s: op should one of +,-,/,*\n", argv[0]);
    }

    writefits(argv[4], argv[1], (char*)img3, -32, nx, ny);

    return 0;
}
