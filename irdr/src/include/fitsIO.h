/* fitsIO.h -- header file for fitsIO.c */

#include "fitsfile.h"              /* from WCSTools */
#include "fitshead.h"

extern float *readfits(char *fn, int *nx, int *ny, float *bkg, float *sig);

extern void writefits(char *fn, char *fnhdr, char *img, int bpp,int nx,int ny);

extern char *get_key_str(char *fn, char *key);

extern int get_key_int(char *fn, char *key, int *val);

extern int get_key_double(char *fn, char *key, double *val);

extern int get_key_float(char *fn, char *key, float *val);

extern void put_key_int(char *fn, char *key, int val);

extern void put_key_double(char *fn, char *key, double val);

extern void put_key_float(char *fn, char *key, float val);

extern void put_key_str(char *fn, char *key, char *val);

extern int get_wcs(char *fn, double *ra, double *dec, double *scale, double *posang);

extern int 
put_wcs(char *fn, double ra, double dec, double scale, int nx, int ny);

