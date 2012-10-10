User Guide
==========

This will contain instructions for end users of the application.


Example of dither data full reduction using PAPI
================================================

Dithered science observations can be reduced using the PAPI recipes using the 
main script 'papi.py'.

The procedure implements a two-pass reduction to improve the estimation of the 
background. In the second pass the objects detected in the combined image are
masked before computing the background.


The raw data are supposed to be located at 

        /data/example_1 
        
and the generated files will be stored at 

        /data/out_ex_1  

1) First, one might want to have a look the the input frames first:

    $papi.py 