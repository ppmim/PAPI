/* skysub.h -- header file for skysub.c */

extern float * 
skysub(float *img, int nx, int ny, float bkg, float *bpm, 
       float *sky, float *skyw, float *mask, char *type);

extern float *
skysub_nomask(float *img, int nx, int ny, float bkg, float *bpm, 
              float *sky, char *type);
