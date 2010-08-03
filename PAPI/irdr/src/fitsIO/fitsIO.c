/* fitsIO.c -- high level wrapper for WCSTools FITS routines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "fitsIO.h"
#include "eprintf.h"
#include "hist.h"

#define ARCSECPERDEG 3600.0

/* writefits: write a FITS file, copy the header from file fnhdr */
extern void writefits(char *fn, char *fnhdr, char *img, int bpp, int nx, int ny)
{
    int lhead, nbhead;
    char *hdr;

    if ((hdr = fitsrhead(fnhdr, &lhead, &nbhead)) == NULL)
        eprintf("writefits: fitsrhead failed: %s\n", fnhdr);

    if (hputi4(hdr, "BITPIX", bpp) < 0)
        eprintf("writefits: write BITPIX fail\n");

    if (hputi4(hdr, "NAXIS1", nx) < 0)
        eprintf("writefits: write NAXIS1 fail\n");

    if (hputi4(hdr, "NAXIS2", ny) < 0)
        eprintf("writefits: write NAXIS2 fail\n");

    if (! fitswimage(fn, hdr, (char *)img))
        eprintf("writefits: fitswimage failed\n");

    if (bpp == 16) {                            /* unsigned short data */
        put_key_float(fn, "BZERO", 32768.0);
        put_key_float(fn, "BSCALE", 1.0);
    } else {
        put_key_float(fn, "BZERO", 0.0);
        put_key_float(fn, "BSCALE", 1.0);
    }

    free(hdr);
}

/* readfits: read a FITS image, convert to floating point */
extern float *readfits(char *fn, int *nx, int *ny, float *bkg, float *sig)
{
    int i, lhead, nbhead, bitpix, npix;
    float *img = NULL, bzero = 0.0, bscale = 1.0;
    char *hdr;

    if ((hdr = fitsrhead(fn, &lhead, &nbhead)) == NULL)
        eprintf("readfits: fitsrhead failed\n");

    if (! hgeti4(hdr, "NAXIS1", nx))
        eprintf("readfits: get NAXIS1 failed: %s\n", fn);

    if (! hgeti4(hdr, "NAXIS2", ny))
        eprintf("readfits: get NAXIS2 failed: %s\n", fn);

    if (! hgeti4(hdr, "BITPIX", &bitpix))
        eprintf("readfits: get BITPIX failed: %s\n", fn);

    if (! hgetr4(hdr, "BZERO", &bzero))
        bzero = 0.0;

    if (! hgetr4(hdr, "BSCALE", &bscale))
        bscale = 1.0;

    npix = (*nx) * (*ny);

    if (bitpix == -32) {                                /* read float image */
        if ((img = (float *)fitsrimage(fn, nbhead, hdr)) == NULL)
            eprintf("readfits: fitsrimage failed\n");

        if (bzero != 0.0 || bscale != 1.0)
            for (i = 0; i < npix; i++)
                img[i] = bscale * img[i] + bzero;

        if (bkg != NULL)
            *bkg = histcalcf(img, *nx, *ny, -1, -1, sig);

    } else if (bitpix == 16) {         /* (un)signed short, convert to float */
        img = (float *) emalloc(npix * sizeof(float));

        if (bzero == 32768.0 && bscale == 1.0) {         /* unsigned short */
            unsigned short *buf;

            if ((buf = (unsigned short *)fitsrimage(fn, nbhead, hdr)) == NULL)
                eprintf("readfits: fitsrimage failed\n");

            if (bkg != NULL)
                *bkg = histcalc(buf, *nx, *ny, -1, -1, sig);

            for (i = 0; i < npix; i++)
                img[i] = (float)buf[i];

            free(buf);

        } else {
            signed short *buf;

            if ((buf = (signed short *)fitsrimage(fn, nbhead, hdr)) == NULL)
                eprintf("readfits: fitsrimage failed\n");

            if (bzero != 0.0 || bscale != 1.0) {    
                for (i = 0; i < npix; i++)
                    img[i] = bscale * buf[i] + bzero;
            } else {
                for (i = 0; i < npix; i++)
                    img[i] = (float)buf[i];
            }

            if (bkg != NULL)
                *bkg = histcalcf(img, *nx, *ny, -1, -1, sig);

            free(buf);
        }
    } else if (bitpix == 32) {         /* signed integer point data, convert to single precision floating */
        img = (float *) emalloc(npix * sizeof(float));

        if (bzero == 32768.0 && bscale == 1.0) {         /* unsigned int */
            unsigned int *buf; /* unsigned int */

            if ((buf = (unsigned int *)fitsrimage(fn, nbhead, hdr)) == NULL)
                eprintf("readfits: fitsrimage failed\n");

            if (bkg != NULL)
                *bkg = histcalc(buf, *nx, *ny, -1, -1, sig);

            for (i = 0; i < npix; i++)
                img[i] = (float)buf[i];

            free(buf);

        } else {
            int *buf;

            if ((buf = (int *)fitsrimage(fn, nbhead, hdr)) == NULL)
                eprintf("readfits: fitsrimage failed\n");

            if (bzero != 0.0 || bscale != 1.0) {
                for (i = 0; i < npix; i++)
                    img[i] = bscale * buf[i] + bzero;
            } else {
                for (i = 0; i < npix; i++)
                    img[i] = (float)buf[i];
            }

            if (bkg != NULL)
                *bkg = histcalcf(img, *nx, *ny, -1, -1, sig);

            free(buf);
        }


    } else {                                  /* can't convert to float type */
        eprintf("readfits: expected bitpix 16, -32 or 32\n");
    }

    free(hdr);

    return img;
}

/* get_key_str: read FITS keyword value into a string */
extern char *get_key_str(char *fn, char *key)
{
    int lhead, nbhead;
    char *hdr, *str, *tmp;

    if ((hdr = fitsrhead(fn, &lhead, &nbhead)) == NULL)
        eprintf("get_key_str: fitsrhead failed\n");

    if ((tmp = hgetc(hdr, key)) == NULL)
        return NULL;

    str = estrdup(tmp);             /* make copy, hgetc uses static buffer */

    free(hdr);

    return str;
}

/* get_key_int: read an integer FITS keyword value */
extern int get_key_int(char *fn, char *key, int *val)
{
    char *str = get_key_str(fn, key);

    if (str == NULL)
        return -1;

    if (sscanf(str, "%d", val) != 1)
        return -1;

    return 0;
}

/* get_key_double: read a double FITS keyword value */
extern int get_key_double(char *fn, char *key, double *val)
{
    char *str = get_key_str(fn, key);

    if (str == NULL)
        return -1;

    if (sscanf(str, "%lf", val) != 1)
        return -1;

    return 0;
}

/* get_key_float: read a float FITS keyword value */
extern int get_key_float(char *fn, char *key, float *val)
{
    char *str;

    if ((str = get_key_str(fn, key)) == NULL)
        return -1;

    if (sscanf(str, "%f", val) != 1)
        return -1;

    return 0;
}

/* put_key_int: add/modify a integer keyword in a FITS file header */
extern void put_key_int(char *fn, char *key, int val)
{
    char str[80];
    sprintf(str, "%d", val);
    put_key_str(fn, key, str);
}

/* put_key_double: add/modify a double keyword in a FITS file header */
extern void put_key_double(char *fn, char *key, double val)
{
    char str[80];
    sprintf(str, "%.10f", val);
    put_key_str(fn, key, str);
}

/* put_key_float: add/modify a float keyword in a FITS file header */
extern void put_key_float(char *fn, char *key, float val)
{
    char str[80];
    sprintf(str, "%f", val);
    put_key_str(fn, key, str);
}

/* put_key_str: write a FITS keyword (does not add '' around val) */
extern void put_key_str(char *fn, char *key, char *val)
{
    int lhead, nbhead, growhdr;
    char *hdr, *img = NULL;
    FILE *fp = NULL;

    if ((hdr = fitsrhead(fn, &lhead, &nbhead)) == NULL)
        eprintf("put_key_str: fitsrhead failed\n");

    growhdr = (!strncmp("END", hdr+strlen(hdr)-80, 3)) ? 1 : 0;

    if (growhdr) {       /* need to enlarge header, rewrite entire FITS file */
        if ((img = fitsrimage(fn, nbhead, hdr)) == NULL)
            eprintf("put_key_str: fitsrimage failed\n");

        if (hputc(hdr, key, val) < 0)
            eprintf("put_key_str: hputc fail: %s %s\n", key, val);

        if (!fitswimage(fn, hdr, img))
            eprintf("put_key: fitswimage failed\n");

        free(img);
    } else {                    /* else just overwrite hdr with modified hdr */
        if ((fp = fopen(fn, "r+")) == NULL)
            eprintf("put_key_str: failed opening %s\n", fn);

        if (hputc(hdr, key, val) < 0)
            eprintf("put_key_str: hputc fail: %s %s\n", key, val);

        if (fwrite(hdr, 1, strlen(hdr), fp) != strlen(hdr))
            eprintf("put_key_str: re-write header failed\n");

        fclose(fp);
    }

    free(hdr);

    return;
}

/* get_wcs: read RA, DEC, SCALE and POSANG from the header of a FITS file */
extern int get_wcs(char *fn, double *ra, double *dec, double *scale, double *posang)
{
    int lhead, nbhead;
/*  float DUPscale = 0.196, INTscale = 0.457;    */
    float LFscale = 0.25, SFscale=0.13, ARNscale=1.0;
/*  char *hdr, *tel, *cam, *objectiv;   */
    char *hdr, *objectiv; 
    char *instrument="0000000000000000000";


    if ((hdr = fitsrhead(fn, &lhead, &nbhead)) == NULL)
        eprintf("get_wcs: fitsrhead failed\n");


    /* -------------- Chech instrument ------------------- */

    /* if (hgets(hdr,"INSTRUME",20,instrument)==1) {	*/

    instrument = hgetc(hdr,"INSTRUME");

    /* --------------------  NICS  ---------------------- */

      if (strncmp(instrument,"NICS",4)==0)  {

	    if (! hgetra(hdr, "RA", ra)) {
	        fprintf(stderr, "get_wcs: unable to read RA in: %s\n", fn);
	        return -1;
	    }
	
	    if (! hgetdec(hdr, "DEC", dec)) {
	        fprintf(stderr, "get_wcs: unable to read DEC in: %s\n", fn);
	        return -1;
	    }
	
	
	    if (hgetr8(hdr, "POSANG", posang)==0) {
	        free(hdr);
	        return 0;
	    }
	
		/* objectiv = hgetc(hdr,"OBJECTIV"); */

		objectiv = hgetc(hdr,"MOTOR6");

	    /* if (hgets(hdr, "OBJECTIV", 80, objectiv)==0) {
	        free(hdr);
	        return 0;
	    }*/
	
	    
	    if (strncmp(objectiv,"LF",2)==0) 
	          *scale = LFscale;
	      else 
	          *scale = SFscale;
      }	    


    /* --------------------  ARNICA --------------------- */

      else if (strncmp(instrument,"ARNI",4)==0)  { 

	    *scale = ARNscale;
	    *posang = -90.0;      /* because N is left and E is down */

	    if (! hgetra(hdr, "RA", ra)) {
	        fprintf(stderr, "get_wcs: unable to read RA in: %s\n", fn);
	        return -1;
	    }
	
	    if (! hgetdec(hdr, "DEC", dec)) {
	        fprintf(stderr, "get_wcs: unable to read DEC in: %s\n", fn);
	        return -1;
	    }

		*dec = -1.0 * *dec; /* because N is left and E is down */
	
	  } 
    /* --------------------  2MASS  --------------------- */

      else if (strncmp(instrument,"2MAS",4)==0)  { 

		  printf ("2MASS\n");

	    *scale = 1.0;
	    *posang = 0.0;       

	    if (! hgetra(hdr, "CRVAL1", ra)) {
	        fprintf(stderr, "get_wcs: unable to read RA in: %s\n", fn);
	        return -1;
	    }
	
	    if (! hgetdec(hdr, "CRVAL2", dec)) {
	        fprintf(stderr, "get_wcs: unable to read DEC in: %s\n", fn);
	        return -1;
	    }
		  printf ("RA=%g  DEC=%g  SC=%g  ANG=%g\n", *ra, *dec, *scale, *posang);
	  }
      
    /* ---------------- WIRCAM ----------------- */
      else if (strncmp(instrument,"WIRCam",6)==0 || strncmp(instrument,"HAWKI",5)==0 )  { 

        /*printf ("WIRCam\n");*/

        *scale = 0.3;
        *posang = 0.0;       

        if (! hgetra(hdr, "CRVAL1", ra)) {
            fprintf(stderr, "get_wcs: unable to read RA in: %s\n", fn);
            return -1;
        }
    
        if (! hgetdec(hdr, "CRVAL2", dec)) {
            fprintf(stderr, "get_wcs: unable to read DEC in: %s\n", fn);
            return -1;
        }
        /*printf ("RA=%g  DEC=%g  SC=%g  ANG=%g\n", *ra, *dec, *scale, *posang);*/
      }
      /* ---------------- OMEGA2000 ----------------- */
      else if (strncmp(instrument,"Omega2000",9)==0){
	
	    fprintf(stderr, "I don't know %s\n", instrument);
	      
	    if (! hgetra(hdr, "RA", ra)) {
	        fprintf(stderr, "get_wcs: unable to read RA in: %s\n", fn);
	        return -1; 
	    }
	
	    if (! hgetdec(hdr, "DEC", dec)) {
	        fprintf(stderr, "get_wcs: unable to read DEC in: %s\n", fn);
	        return -1; 
	    }
	
	    if (hgetr8(hdr, "PIXSCALE", scale)==0) {
	        fprintf(stderr, "get_wcs: unable to read SCALE in: %s\n", fn);
	        free(hdr);
	        return -1; 
	    }


	    if (hgetr8(hdr, "ROT-RTA" /*"POSANG"*/, posang)==0) {
	        fprintf(stderr, "get_wcs: unable to read POSANG in: %s\n", fn);
	        free(hdr);
	        return -1; 
	    }  

	  /*fprintf(stderr, "Read RA=%g  DEC=%g  SC=%g  ANG=%g\n", *ra, *dec, *scale, *posang);*/ 

	 return 0;
	}
    /* ---------------- OTHER INSTRUMENT ----------------- */
    else {
	 fprintf(stderr, "UNKNOWN INSTRUMENT in %s\n", fn);
	 return -1;
    } 
/*
    if (hgetr8(hdr, "CDELT", scale) || hgetr8(hdr, "CDELT1", scale)) {
        *scale = *scale * ARCSECPERDEG;
        free(hdr);
        return 0;
    }

    if ((cam = hgetc(hdr, "CAMERA")) == NULL)
        eprintf("get_wcs: failed reading image scale / CAMERA keywords");
   

    if (strcmp(cam, "cirsi") && strcmp(cam, "CIRSI"))         if not CIRSI    
        eprintf("get_wcs: failed reading image scale");

    if ((tel = hgetc(hdr, "TELESCOP")) == NULL)
        eprintf("get_wcs: failed reading image scale / TELESCOP keywords");

    if (!strcmp("INT", tel))                                   INT data   
        *scale = INTscale;
    else if (!strcmp("LCO", tel))                           du Pont data    
        *scale = DUPscale;
    else
        eprintf("Expected TELESCOP keyword starts with LCO or INT\n");
*/


    free(hdr);

    return 0;
}

/* put_wcs: write rough WCS to FITS header (for CIRSI) */
extern int put_wcs(char *fn, double ra, double dec, double scale, 
                   int nx, int ny)
{
    int lhead, nbhead;
    double degscale = scale / ARCSECPERDEG;
    char *hdr, *img, str[81], comment[81];

    if ((hdr = fitsrhead(fn, &lhead, &nbhead)) == NULL)
        eprintf("put_wcs: fitsrhead failed\n");

    if (hputr8(hdr, "RA", ra) < 0) {
        
	fprintf(stderr, "put_wcs: failed writing RA\n");
        exit(1);
    }

    ra2str(str, 30, ra, 3);

    sprintf(comment, "[degrees]  (%s)        ", str);

    if (hputcom(hdr, "RA", comment) < 0) {
        fprintf(stderr, "put_wcs: write comment fail\n");
        exit(1);
    }

    if (hputr8(hdr, "DEC", dec) < 0) {
        fprintf(stderr, "put_wcs: failed writing DEC\n");
        exit(1);
    }

    dec2str(str, 30, dec, 2);

    sprintf(comment, "[degrees]  (%s)        ", str);

    if (hputcom(hdr, "DEC", comment) < 0) {
        fprintf(stderr, "put_wcs: write comment fail\n");
        exit(1);
    }

    if (hputr8(hdr, "SCALE", scale) < 0 ||
        hputr8(hdr, "CRPIX1", nx/2.) < 0 ||
        hputr8(hdr, "CRPIX2", ny/2.) < 0 ||
        hputr8(hdr, "CRVAL1", ra) < 0 ||
        hputr8(hdr, "CRVAL2", dec) < 0 ||
        hputs(hdr, "CTYPE1", "RA---TAN") < 0 ||
        hputs(hdr, "CTYPE2", "DEC--TAN") < 0 ||
        /*hputs(hdr, "RADECSYS", "FK5") < 0 ||
        hputs(hdr, "WAT0_001", "system=image") < 0 ||
        hputs(hdr, "WAT1_001", "wtype=tan axtype=ra") < 0 ||
        hputs(hdr, "WAT2_001", "wtype=tan axtype=dec") < 0 ||
        hputi4(hdr, "WCSDIM", 2) < 0 || */
        hputr8(hdr, "CD1_1", -degscale) < 0 ||
        hputr8(hdr, "CD1_2", 0.0) < 0 ||
        hputr8(hdr, "CD2_1", 0.0) < 0 ||
        hputr8(hdr, "CD2_2", degscale) < 0 )
        eprintf("put_wcs: failed writing WCS keys\n");

    if ((img = fitsrimage(fn, nbhead, hdr)) == NULL)
        eprintf("put_wcs: fitsrimage failed\n");

    if (!fitswimage(fn, hdr, img))
        eprintf("put_wcs: fitswimage failed\n");

    free(hdr);

    return 0;
}


