/* stats.c - use hist.c to calculate image statistics */

/* 
 * calculate and optionally write to FITS header:
 * DATAMIN, DATAMAX, DATAMEAN, DATAMODE, DATAMED, DATASIG
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int i, j, nx, ny, updatehdr = 0;
    float *img, med, minv, maxv;
    float mean, mode, sig;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s [updatehdr] *.fits\n", argv[0]);
        exit(1);
    }
    
    printf("\tFILENAME\tMIN\tMAX\tMEAN\tMODE\tMED\tSIG\n");

    if (strcmp(argv[1], "updatehdr") == 0)
        updatehdr = 1;

    for (i = 1 + updatehdr; i < argc; i++) {
        img = readfits(argv[i], &nx, &ny, &mode, &sig);      /* MODE, SIGMA */

        minv = maxv = img[0];

        for (j = 1; j < nx * ny; j++) {                   /* MIN, MAX values */
            if (img[j] < minv)
                minv = img[j];

            if (img[j] > maxv)
                maxv = img[j];
        }

        med = median(img, nx * ny);                        /* MEDIAN */

        mean = 0.0;

        for (j = 0; j < nx * ny; j++)                     /* MEAN */
            mean += (float)img[j];

        mean /= (nx * ny);
        
        printf("\n\n");
        printf("\tFILENAME\tMIN\tMAX\tMEAN\tMODE\tMED\tSIG\n");
        printf("%s\t%f\t%f\t%9.2f\t%9.2f\t%f\t%9.2f\n", 
            argv[i], (float)minv, (float)maxv, mean, mode, (float)med, sig);

        if (updatehdr) {
            put_key_float(argv[i], "DATAMIN", (float)minv);
            put_key_float(argv[i], "DATAMAX", (float)maxv);
            put_key_float(argv[i], "DATAMEAN", (float)mean);
            put_key_float(argv[i], "DATAMODE", (float)mode);
            put_key_float(argv[i], "DATAMED", (float)med);
            put_key_float(argv[i], "DATASIG", (float)sig);
        }

        free(img);  img = NULL;
    }

    return 0;
}
