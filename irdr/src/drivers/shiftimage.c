/* shiftimage.c -- test shift.c */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

int main(int argc, char *argv[])
{
    int nx, ny, border;
    float xshift, yshift, *img, *wimg, *imgout, *wimgout;

    if (argc != 8)
      eprintf("Usage: %s img.fits wimg.fits shift.fits wshift.fits border "
              "xshift yshift", argv[0]);

    img = readfits(argv[1], &nx, &ny, NULL, NULL);
    wimg = readfits(argv[2], &nx, &ny, NULL, NULL);

    border = atoi(argv[5]);
    xshift = (float)atof(argv[6]);
    yshift = (float)atof(argv[7]);

    imgout = shift_image(img, wimg, nx, ny, border, xshift, yshift, &wimgout);

    nx = nx + 2 * border;
    ny = ny + 2 * border;

    writefits(argv[3], argv[1], (char*)imgout, -32, nx, ny);
    writefits(argv[4], argv[1], (char*)wimgout, -32, nx, ny);
    
    return 0;
}
