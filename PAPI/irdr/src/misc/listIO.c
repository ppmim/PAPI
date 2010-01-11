/* listIO.c -- functions to read ASCII lists */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "listIO.h"
#include "eprintf.h"

/* 
 * readlist: read a list of string, string, float, float 
 */

extern int readlist(char *list, char *fn1[], char *fn2[], float *arr1, 
                    float *arr2, int nmax)
{
    int nfile = 0;
    static char line[512], file1[256], file2[256];
    float val1, val2;
    FILE *fp;

    if ((fp = fopen(list, "r")) == NULL)
        eprintf("readlist: failed opening: %s\n", list);

    while (fgets(line, sizeof(line), fp) != NULL) {
        if (fn1 != NULL && fn2 == NULL && arr1 == NULL && arr2 == NULL) {
            if (sscanf(line, "%s", file1) != 1)
                eprintf("readlist: check list format\n");

        } else if (fn1 != NULL && fn2 != NULL && arr1 == NULL && arr2 == NULL) {
            if (sscanf(line, "%s %s", file1, file2) != 2)
                eprintf("readlist: check list format\n");

        } else if (fn1 != NULL && fn2 == NULL && arr1 != NULL && arr2 != NULL) {
            if (sscanf(line, "%s %f %f", file1, &val1, &val2) != 3)
                eprintf("readlist: check list format\n");

        } else if (fn1 != NULL && fn2 != NULL && arr1 != NULL && arr2 != NULL) {
            if (sscanf(line, "%s %s %f %f", file1, file2, &val1, &val2) != 4)
                eprintf("readlist: check list format\n");

        } else {
            eprintf("readlist: unexpected list format\n");
        }

        fn1[nfile] = estrdup(file1);

        if (fn2 != NULL)
            fn2[nfile] = estrdup(file2);

        if (arr1 != NULL)
            arr1[nfile] = val1;

        if (arr2 != NULL)
            arr2[nfile] = val2;

        if (++nfile >= nmax)
            eprintf("readlist: file list too long, increase NMAX\n");
    }

    fclose(fp);

    return nfile;
}
