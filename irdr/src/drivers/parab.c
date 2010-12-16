/* parab.c -- test parabola.c */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

#define N 3

int main(int argc, char *argv[])
{
    float arr1[N] = {20330., 21435., 21181.};

    printf("%f\n", fit_parabola(arr1, N));
    
    return 0;
}
