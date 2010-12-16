/* getwmap.c -- produce a weight map from exptime * gainmap / variance */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "irdr.h"

static float *readobjm(char *fn, int *nxm, int *nym);

/* getwmap: produce a weight map from exptime * gainmap / variance */
extern float *getwmap(char *fn, int nx, int ny, float *gain, float sigma)
{
    int i;
    float *wmap, exptime, variance, scale, ncombine;

    if (get_key_float(fn, "INT_TIME", &exptime) < 0)
        if (get_key_float(fn, "EXPTIME", &exptime) < 0)
            eprintf("getwmap: ERR: keyword INT_TIME/EXPTIME not found\n");

    if (get_key_float(fn, "NCOMBINE", &ncombine) < 0)
        ncombine = 1.0;

    variance = (sigma > 0.0) ? (sigma * sigma) : 1.0;

    scale = ncombine * exptime / variance;

    wmap = (float *) emalloc(nx * ny * sizeof(float));

    for (i = 0; i < nx * ny; i++)
        wmap[i] = scale * gain[i];

    printf("  getwmap: %s %f %f %f\n", fn, ncombine, exptime, variance);

    return wmap;
}

/* getmask: mask object pixels in a weight map */
extern float *getmask(float *wmap, int nx, int ny, char *objmfn, 
                      float xoff, float yoff)
{
    int i, j, ixoff, iyoff;
    int nxobj, nyobj, border;
    float *mask, *objmask;         /* mask is merged weight map and objmask */

    mask = (float *) emalloc(nx * ny * sizeof(float));

    objmask = readobjm(objmfn, &nxobj, &nyobj);

    border = (nxobj - nx) / 2;       /* objmask is from coadded dither set */

    /* offset of current frame into dither set objmask */

    ixoff = (int)(border - xoff);
    iyoff = (int)(border - yoff);

    for (i = 0; i < ny; i++) {
        float *maskrow = mask + (i * nx);
        float *wmaprow = wmap + (i * nx);
        float *objrow = objmask + ((i + iyoff) * nxobj + ixoff);

        for (j = 0; j < nx; j++)
            if (objrow[j] > 0)                 /* if object pixel... */
                maskrow[j] = 0.0;           /* mask out weightmap pixel */
            else
                maskrow[j] = wmaprow[j];
    }

    printf("  getmask: %d %d %d\n", border, ixoff, iyoff);

    return mask;
}

/* readobjm: read in next object mask if different than previously read mask */
static float *readobjm(char *fn, int *nxm, int *nym)
{
    static char *objmfn = NULL;
    static int nxobjm, nyobjm;
    static float *objm;

    if (objmfn == NULL || strcmp(fn, objmfn) != 0) {   /* objm not read yet? */
        if (objmfn != NULL) {
            free(objm); free(objmfn);
        }

        objm = readfits(fn, &nxobjm, &nyobjm, NULL, NULL);

        objmfn = estrdup(fn);
    }

    *nxm = nxobjm;
    *nym = nyobjm;

    printf("  readobjm: %s\n", fn);

    return objm;
}
