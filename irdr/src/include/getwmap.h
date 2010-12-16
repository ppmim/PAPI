/* getwmap.h -- header file for getwmap.c */

extern float *getwmap(char *fn, int nx, int ny, float *gain, float sigma);

extern float * getmask(float *wmap, int nx, int ny, char *objmfn, 
                       float xoff, float yoff);
