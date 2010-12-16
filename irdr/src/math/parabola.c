/* parabola.c -- find coordinate of maximum from parabolic fit */

#include <stdio.h>
#include <stdlib.h>
#include "parabola.h"

/* 
 * fit a parabola to the peak (max) of a 1-D array and return parabola center
 * see Bevington, page 210
 */

extern float fit_parabola(const float *arr, int n)
{
    int i, x0 = 0;
    float d, maxv = arr[0];

    for (i = 0; i < n; i++)                   /* find peak */
        if (arr[i] > maxv)
            maxv = arr[x0 = i];

    if (x0 < 1)                             /* abort if peak at edge */
        return 0.0;

    if (x0 > n - 2)
        return n - 1;

    d = arr[x0+1] - (2.0 * arr[x0]) + arr[x0-1];

    return (d == 0) ? (float)x0 : (x0 + 0.5 - (arr[x0+1] - arr[x0]) / d);
}
