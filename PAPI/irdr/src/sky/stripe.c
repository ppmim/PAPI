/* stripe.c -- do image destriping by subtracting row/column modes */

/* input bpm indicates bad pixels and object pixels (weight 0.0) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "irdr.h"

#define NQUAD 4
#define QUADSIZE 1024 /*snap =>512*/
/* #define QUADSIZE 128 */

static float col_mode[QUADSIZE], row_mode[QUADSIZE];

static void fixmodes(float *arr, float bkg);
static void dorows(float *img, float *bpm, int nx, int ny, float bkg);
static void docols(float *img, float *bpm, int nx, int ny, float bkg);

/* image destriping, either rows, columns, or both */
extern void destripe(float *img, float *bpm, int nx, int ny, float bkg, 
                     char *type)
{

    /*printf("\nDestriping [%s] ....\n", type);*/
    if (!strcmp(type, "rowcol")) {
        dorows(img, bpm, nx, ny, bkg);
        docols(img, bpm, nx, ny, bkg);
    } else if (!strcmp(type, "colrow")) {
        docols(img, bpm, nx, ny, bkg);
        dorows(img, bpm, nx, ny, bkg);
    } else if (!strcmp(type, "row")) {
        dorows(img, bpm, nx, ny, bkg);
    } else if (!strcmp(type, "col")) {
        docols(img, bpm, nx, ny, bkg);
    } else {
        /* fprintf(stderr, "destripe: warning: no destripe done\n"); */
    }
}

/* destripe columns of an image: subtract column modes, quadrant by quadrant */
extern void docols(float *img, float *bpm, int nx, int ny, float bkg)
{
    int i, j, k;

    if (nx != (2 * QUADSIZE) || ny != (2 * QUADSIZE))
        eprintf("docols: expected img size %d\n", 2 * QUADSIZE);

    for (k = 0; k < NQUAD; k++) {                       /* foreach quadrant */
        int xoffset = (k % 2) * QUADSIZE;
        int yoffset = (k / 2) * QUADSIZE;
        static float buf[QUADSIZE];

        for (j = xoffset; j < xoffset + QUADSIZE; j++) {     /* each column */
            int nbuf = 0;

            for (i = yoffset; i < yoffset + QUADSIZE; i++)      /* each row */
                if (bpm[i*nx+j] > 0.0)
                    buf[nbuf++] = img[i*nx+j];          /* non-masked pixels */
            
            if (nbuf > QUADSIZE / 2)                     /* get column mode */
                col_mode[j-xoffset] = histcalca(buf, nbuf);
            else
                col_mode[j-xoffset] = 0;
        }

        fixmodes(col_mode, bkg);

        for (i = yoffset; i < yoffset + QUADSIZE; i++) {    /* do subtract */
            float *row = img + i * nx;

            for (j = xoffset; j < xoffset + QUADSIZE; j++)
                row[j] += (bkg - col_mode[j-xoffset]);
        }
    }
}

/* destripe rows of an image: subtract row modes, quadrant by quadrant */
extern void dorows(float *img, float *bpm, int nx, int ny, float bkg)
{
    int i, j, k;

    if (nx != (2 * QUADSIZE) || ny != (2 * QUADSIZE))
        eprintf("dorows: expected image size %d\n", 2 * QUADSIZE);

    for (k = 0; k < NQUAD; k++) {                       /* foreach quadrant */
        int xoffset = (k % 2) * QUADSIZE;
        int yoffset = (k / 2) * QUADSIZE;
        static float buf[QUADSIZE];

        for (i = yoffset; i < yoffset + QUADSIZE; i++) {        /* each row */
            int nbuf = 0;

            for (j = xoffset; j < xoffset + QUADSIZE; j++)   /* each column */
                if (bpm[i*nx+j] > 0)
                    buf[nbuf++] = img[i*nx+j];          /* non-masked pixels */
            
            if (nbuf > QUADSIZE / 2)                         /* get row mode */
                row_mode[i-yoffset] = histcalca(buf, nbuf);
            else
                row_mode[i-yoffset] = 0;
        }

        fixmodes(row_mode, bkg);

        for (i = yoffset; i < yoffset + QUADSIZE; i++) {     /* do subtract */
            float mode = row_mode[i-yoffset];
            float *row = img + i * nx;

            for (j = xoffset; j < xoffset + QUADSIZE; j++)
                row[j] += (bkg - mode);
        }
    }
}

/* fixmodes: fill in missing points with nearest neighbor interpolation */
static void fixmodes(float *arr, float bkg)
{
    int i, n = QUADSIZE;
    float fixmode = bkg;

    for (i = 0; i < n; i++)
        if (arr[i] > 0) {
            fixmode = arr[i];
            break;
        }

    for (i = 0; i < n; i++)
        if (arr[i] == 0)
            arr[i] = fixmode;
        else
            fixmode = arr[i];
}
