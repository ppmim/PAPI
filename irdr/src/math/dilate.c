/* dilate.c -- grow an object mask by multiplicative scaling of object sizes */

/*
 * object mask is an OBJECTS output image from SExtractor.  input scale is a
 * multiplicative scale factor: object regions expands by (1+scale).
 */

#include <stdio.h>
#include <stdlib.h>
#include "eprintf.h"
#include "dilate.h"
#include "minmax.h"

static void dilate_row(float *row, int nx, float scale);

/* grow object mask regions by a multiplicative factor of scale */ 
extern void dilate(float *mask, int nx, int ny, float scale)
{
    int i, j;
    float *col;

    printf("dilate: scale factor = %f\n", scale);

    if (scale <= 0)
        return;

    for (i = 0; i < ny; i++)                         /* dilate each row */
        dilate_row(mask + i*nx, nx, scale);

    col = (float*) emalloc(ny * sizeof(float));

    for (i = 0; i < nx; i++) {                       /* foreach column */
        for (j = 0; j < ny; j++)                     /* fill column buffer */
            col[j] = mask[j*nx+i];

        dilate_row(col, ny, scale);                  /* dilate each column */

        for (j = 0; j < ny; j++)
            if (col[j] > 0.0)                        /* update mask */
                mask[j*nx+i] = 1.0;
    }

    free(col);
}

/* apply object mask scaling to a 1-D array */
static void dilate_row(float *row, int nx, float scale)
{
    int j;

    for (j = 0; j < nx; j++) {
        if (row[j] > 0.0) {                           /* if mask region */
            int size, jbeg = j, jend;

            while (row[j++] > 0.0 && j < nx) ;        /* calculate extent */

            jend = j - 1;                             /* scale mask region */
            size = jend - jbeg + 1;
            jbeg = jbeg - scale * size / 2.0 + 0.5;    /* + 0.5 to round */
            jend = jend + scale * size / 2.0 + 0.5;
            jbeg = max(0, jbeg);
            jend = min(nx - 1, jend);

            for (j = jbeg; j <= jend; j++)            /* update (grow) mask */
                row[j] = 1.0;
        }
    }
}
