/* avgkey.c -- average the value of a FITS key over a list of FITS files */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

#define NSIG 5.0                        /* clip points > NSIG from median */

int main(int argc, char *argv[])
{
    int i, navg, nval = argc - 2;
    double *vals, med, sig, lcut, hcut, avg;

    if (argc < 4)
        eprintf("Usage: %s KEY *.fits\n", argv[0]);

    vals = (double *) emalloc(nval * sizeof(double));

    for (i = 2; i < argc; i++)
        if (get_key_double(argv[i], argv[1], &vals[i-2]) < 0)
            eprintf("%s: keyword %s not found\n", argv[0], argv[1]);

    med = median_double(vals, nval);
    sig = median_absdev_double(vals, med, nval) / 0.6745; 

    lcut = med - NSIG * sig;
    hcut = med + NSIG * sig;

    avg = 0.0; 
    navg = 0;

    for (i = 0; i < nval; i++)
        if (vals[i] > lcut && vals[i] < hcut) {
            avg += vals[i];
            navg++;
        }

    if (navg < 1)
        eprintf("%s: navg = %d\n", argv[0], navg);

    avg /= navg;
 
    for (i = 2; i < argc; i++)
        put_key_double(argv[i], argv[1], avg);

    return 0;
}
