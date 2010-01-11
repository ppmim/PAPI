#include <stdio.h>
#include <iostream.h>

main ()
{
    float **matrix;
    int i=0,j=0;

    /* Alloc memory */
    matrix = new float *[2048];
    for (i=0;i<2048;i++)
        matrix[i]=new float[2048];
        
    for (i=0;i<2048;i++)
        for (j=0;j<2048;j++)
            if (matrix[i][j]>0)
                printf("VAAA\n");
            else matrix[i][j]=1; /*printf("VAL= %f\n", matrix[i][j]);*/
                
    /* Free memory*/
    for (i=0;i<2048;i++)
        delete [] matrix[i];
    delete [] matrix;
    
                
    return 0;
}

