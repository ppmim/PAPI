/* bisearch.c -- binary search */

#include <stdio.h>
#include <stdlib.h>
#include "bisearch.h"

/* binary search of array a with linearly interpolated output index */
extern float bisearch(int x, const int *a, int n)
{
    int middle, left = 0, right = n - 1;
    int i;
    
    if (x <= a[left])
    {   /*printf("\nDEBUG- BINGOOOO !!! x=%d, a[left]=%d", x, a[left]);*/
        /*for ( i=0; i<n; i++) printf("\nDEBUG-a[%d]=%i", i, a[i]);*/
        return 0;
    }
    if (x > a[right])
        return n;

    while (right - left > 1) {
        middle = (right + left) / 2;

        if (x <= a[middle])
            right = middle;
        else
            left = middle;
    }

    /* add 0.5 to treat array as histogram with bins centered on integers */

    return (float)left + ((float)x - a[left]) / (a[right] - a[left]) + 0.5;
}
