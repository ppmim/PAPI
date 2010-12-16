/* convert.c -- convert data type of image */

#include <stdio.h>
#include <stdlib.h>
#include "convert.h"
#include "eprintf.h"


/* 
 * shortint: convert a float image to a u_short image 
 */

extern unsigned short *shortint(float *fimg, int nx, int ny)
{
    int i, n = nx * ny;
    unsigned short *img;


    img = (unsigned short *) emalloc(n * sizeof(unsigned short));
    
    for (i = 0; i < n; i++)
        img[i] = (unsigned short) (fimg[i] + 0.5);

    /* free(fimg); */

    return img;
}
/*
 * longint: convert a float image to a u_int image
 */

extern unsigned int *longint(float *fimg, int nx, int ny)
{
    int i, n = nx * ny;
    unsigned int *img;

    img = (unsigned int *) emalloc(n * sizeof(unsigned int));

    for (i = 0; i < n; i++){
        img[i] = (unsigned int) (fimg[i] + 0.5);
    }

    /* free(fimg); */

    return img;
}
