/* kselect.c -- find element at kth position of sorted array */

#include <stdio.h>
#include <stdlib.h>
#include "kselect.h"

/*
 * given an array a of n elements, return the element that would be at 
 * position k, (0 <= k < n), if the array were sorted.  from Algorithms 
 * and Data Structures in C++ by Leendert Ammeraal, pg. 82.  O(n).
 *
 * NB: partially reorders data in array
 */

extern float kselect(float *a, int n, int k)
{
    while (n > 1) {
        int i = 0, j = n - 1;
        float x = a[j/2], w;

        do {
            while (a[i] < x) i++;
            while (a[j] > x) j--;
            if (i < j) { 
                w = a[i]; a[i] = a[j]; a[j] = w; 
            } else {
                if (i == j) i++;
                break;
            }
        } while (++i <= --j);

        if (k < i) 
            n = i; 
        else { 
            a += i; n -= i; k -= i; 
        }
    }

    return a[0];
}


/* 
 * kselect() version for double precision data 
 */

extern double kselect_double(double *a, int n, int k)
{
    while (n > 1) {
        int i = 0, j = n - 1;
        double x = a[j/2], w;

        do {
            while (a[i] < x) i++;
            while (a[j] > x) j--;
            if (i < j) { 
                w = a[i]; a[i] = a[j]; a[j] = w; 
            } else {
                if (i == j) i++;
                break;
            }
        } while (++i <= --j);

        if (k < i) 
            n = i; 
        else { 
            a += i; n -= i; k -= i; 
        }
    }

    return a[0];
}
