/* skyfilter.c -- given a list of FITS files do running sky subtraction */

/*
 * Procedure:
 *   Call cube_mean() to create a sky frame from the 2*hwid nearest frames
 *   Call skysub() to subtract the sky frame from the current frame
 *   Write the sky-subtracted frames to disk
 *
 * If the mask option is selected then the input filelist should contain:
 * imagefn objfn xshift yshift, where objfn is the object mask from the 
 * coadded dither set, and xshift,yshift are the dither offsets per frame.
 * If the nomask option is selected then the filelist should contain: imagefn.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "irdr.h"

#define MAXHWID 20                     /* max hwidth in frames of sky filter */

static char *fn [MAXNPLANES];          /* FITS image file per image plane */
static char *mfn [MAXNPLANES];         /* FITS objmask file per image plane */
static float *data [MAXNPLANES];       /* pointers to image planes */
static float *wdata [MAXNPLANES];      /* pointers to mask/weight map planes */
static float xshift [MAXNPLANES];      /* dither offset in x direction */
static float yshift [MAXNPLANES];      /* dither offset in y direction */
static float bkgs [MAXNPLANES];        /* bkg level per image plane */
static float sigs [MAXNPLANES];        /* sigma per image plane */
static float *gainmap;

static float *dbuf [2*MAXHWID];        /* store image ptrs to make sky frame */
static float *wbuf [2*MAXHWID];        /* collect corresponding mask ptrs */
static float scale [2*MAXHWID];        /* normalizations to make sky frame */

static char *outfn(char *fn);
static void readdata(int i, int usemask);
static void freedata(int i, int usemask);

static void usage(void);

int main(int argc, char *argv[])
{
    int i, nx, ny, nplanes;
    int skymodel = 1; /* 1= median&mean, 2= minimun */
    int hwid, skybeg, usemask = 0;
    unsigned int *skysubimg = NULL;
    float *sky = NULL, *skyw = NULL, *fimg;
    char aux[256];
    int fix_type=1;
    
    fprintf(stderr, "\nStart of skyfilter...\b");
    
    if (argc < 6)
        usage();

    if (!strcmp(argv[4], "mask"))
        usemask = 1;
    else if (!strcmp(argv[4], "nomask"))
        usemask = 0;
    else
        usage();

    if (argc==7) /* skymodel is an optional argument */
        if (!strcmp(argv[6], "median") || !strcmp(argv[6], "mean"))
            skymodel = 1;
        else if (!strcmp(argv[6], "min"))
            skymodel = 2;
        else
            usage();
    else skymodel = 1;    

    /* fix type */
    if (argc==8)
        fix_type = atoi(argv[7]); /* if >0, some BadPix fix is done (replace with backlvl)*/
    else
        fix_type = 0; /* default: no BadPix fix is done (NaNs) */
    
    if (usemask)
        nplanes = readlist(argv[1], fn, mfn, xshift, yshift, MAXNPLANES);
    else
        nplanes = readlist(argv[1], fn, NULL, NULL, NULL, MAXNPLANES);


    if (nplanes < 1)
        eprintf("[%s] Error: no valid image planes", argv[0]);

    gainmap = readfits(argv[2], &nx, &ny, NULL, NULL);

    hwid = atoi(argv[3]);

    if (hwid > MAXHWID) {
      hwid = MAXHWID;
      printf("[skyfilter] HALFNSKY reduced to %d \n", hwid);
    }


    if (hwid > MAXHWID || (2 * hwid + 1) > nplanes){
      eprintf("[skyfilter] Error: found wrong number of sky frames --> hwid %d, MAXHWID %d, nplanes %d\n", hwid, MAXHWID, nplanes);
      return -1;
    }


    /* Read the initial sliding window */
    for (i = 0; i < (2 * hwid + 1); i++)  {
	    printf("Nplanes: %d  i: %d \n", nplanes, i);
	    if (i<nplanes) readdata(i, usemask);
	}


    for (i = 0, skybeg = 0; i < nplanes; i++) {
        int j, nsky = 0, skyend = skybeg + 2 * hwid;
        float avgscale = 0.0;

	    if (skyend>=nplanes) skyend=nplanes-1;
	    printf("Image: %d   Sky: ", i);

        for (j = skybeg; j <= skyend; j++) {  /* collect adjacent frame ptrs */
            if (j != i) {                             /* skip current frame */
		printf (" %d", j);
                dbuf[nsky] = data[j];

                if (usemask)
                    wbuf[nsky] = wdata[j];

                avgscale += (scale[nsky] = bkgs[j]);
                nsky++;
            }
        }
	printf (" \n");

        if (nsky==0) {
            eprintf("[skyfilter] Error: not enought number of sky frames");
            return -1;
        }
        
        avgscale /= (float) nsky;

        for (j = 0; j < nsky; j++){
            scale[j] = avgscale - scale[j];
            printf("\nSCALE[%d]=%f", j, scale[j]);
        }
        
        if (usemask) {
            if (skymodel==1)
                /* NOTE: (clipped) mean provide a higher S/N ratio than median */
                sky = cube_mean(dbuf, wbuf, nsky, nx, ny, &skyw, scale, 1);
            else
                sky = cube_mean_min_w(dbuf, wbuf, nsky, nx, ny, &skyw, scale, 1, nsky/2);

#ifdef SKY_DEBUG
            strcpy(aux, "/data2/tmp/sky_");
            strcat(aux, basename(fn[i]));
            writefits(aux, fn[i], (char*)sky, -32, nx, ny);
                
            strcpy(aux, "/data2/tmp/skyw_");
            strcat(aux, basename(fn[i]));
            writefits(aux, fn[i], (char*)skyw, -32, nx, ny);
                
#endif
            fimg = skysub(data[i], nx, ny, bkgs[i], gainmap, sky, skyw, 
                            wdata[i], argv[5], fix_type);
        } else {
            if (skymodel==1) 
                sky = cube_median(dbuf, nsky, nx, ny, scale, 1);
            else 
                /* find mean value of the nsky/2-smallest values */
                /* sky = cube_mean_min(dbuf, nsky, nx, ny, scale, 1, nsky/2);*/
                
                /* find median value of the nsky/2-smallest values */
                sky = cube_median_min(dbuf, nsky, nx, ny, scale, 1, nsky/2);
                 
#ifdef SKY_DEBUG
           strcpy(aux, "/data2/tmp/sky_");
           strcat(aux, basename(fn[i]));
           writefits(aux, fn[i], (char*)sky, -32, nx, ny); 
#endif
           
           fimg = skysub_nomask(data[i], nx, ny, bkgs[i], gainmap, sky, 
                                   argv[5], fix_type);
        }

        /*skysubimg = longint(fimg, nx, ny);*/
        /* For PANIC, we need 32 bits images, so we write  -32 (float) FITS*/
        writefits(outfn(fn[i]), fn[i], (char*)fimg, -32, nx, ny);
        
        /* Add HISTORY to header */
        if (usemask) put_key_str(outfn(fn[i]), "HISTORY", "Sky subtracted with object mask");
        else put_key_str(outfn(fn[i]), "HISTORY", "Sky subtracted with NO object mask");
        
        /* free memory */
        free(sky);  /*free(skysubimg);*/  free(fimg);

        if (skyw != NULL)
            free(skyw);

        if (i >= hwid && i < (nplanes - hwid - 1)) {  /* move sliding window */
            freedata(i - hwid, usemask);
            readdata(i + hwid + 1, usemask);
            skybeg++;
        }
    }

    return 0;
}

/* readdata: load in data for image plane i */
static void readdata(int i, int usemask)
{
    int nx, ny;
    char aux[128];
    
    data[i] = readfits(fn[i], &nx, &ny, &bkgs[i], &sigs[i]);  /* image plane */

    if (bkgs[i] <= 0 || sigs[i] <= 0)
        eprintf("[skyfilter::readdata] Error: file %s, bkg %f, sig %f\n", fn[i], bkgs[i], sigs[i]);

    if (usemask) {
        /* getwmap: produce a weight map from exptime * gainmap / variance */
        float *wmap = getwmap(fn[i], nx, ny, gainmap, sigs[i]);
        /* getmask: mask 'object pixels' & 'bad pixels' (gain) in a weight map 
         * taking into account offsets of the main 'objectMask';
           Thus, wdata includes 'object & bad' pixels.
         */
        wdata[i] = getmask(wmap, nx, ny, mfn[i], xshift[i], yshift[i]);

#ifdef SKY_DEBUG
        sprintf(aux, "/data2/tmp/wmask_%d.fits", i);
        writefits(aux, fn[i], (char*)wdata[i], -32, nx, ny);
#endif
        free(wmap);
    }

    fprintf(stderr, " Reading %s bkg=%f sigs=%f\n", fn[i], bkgs[i], sigs[i]);
}

/* freedata: free data for image plane i */
static void freedata(int i, int usemask)
{
    free(data[i]);  data[i] = NULL;

    if (usemask) {
        free(wdata[i]);  wdata[i] = NULL;
    }
}

/* outfn: create output fn for sky subtracted image */
static char *outfn(char *fn)
{
    static char *ext = ".skysub";
    static char buf[256];

    if (strlen(fn) + strlen(ext) >= 256)
        eprintf("outfn: increase OUTLEN");

    sprintf(buf, "%s%s", fn, ext);

    return buf;
}

/* print usage and exit */
static void usage(void)
{
    static char *usage = "\n"
    "skyfilter - do running sky frame subtraction\n\n"
    "usage: skyfilter listfn gainfn hwidth mask|nomask "
    "row|col|rowcol|colrow|none skymode fix_type\n\n"
    "where listfn - if object masking is used, then listfn should contain:\n"
    "               img_filename objmask_filename dither_x_off dither_y_off\n"
    "               where objmask is the master object mask per dither set\n"
    "               and the dither offsets are given in pixels.  eg:\n"
    "               img1set1.fits set1objmask.fits   0.0   0.0\n"
    "               img2set1.fits set1objmask.fits -24.5  29.4\n"
    "               img3set1.fits set1objmask.fits -34.0 -19.3\n"
    "               img1set2.fits set2objmask.fits   0.0   0.0\n"
    "               img2set2.fits set2objmask.fits  14.0 -29.2\n"
    "               for no masking, listfn should contain:\n"
    "               image_filename\n"
    "      gainfn - filename of FITS gainmap (normalized flat field)\n"
    "      hwidth - half width of sky filter in frames\n"
    "      mask   - either mask for object masking or nomask\n"
    "      row|co - destriping correction with the following possibilities:\n"
    "               row for row offsets,\n"
    "               col for column offsets,\n"
    "               rowcol for row offsets then column offsets,\n"
    "               colrow for column offsets then row offsets,\n"
    "               none for no correction\n\n"
    "     skymodel- sky background model to use (optional, defaul=median"
    "               'median'/'mean' for normal/sparse fields"
    "               'min' for crowded fields"
    "     fix_type- type of behaviour with bad pixel (1=replace with bcklvl, 0=nan)"
    "example: skyfilter filelist gain.fits 4 mask rowcol mean 1\n\n";

    printf("%s", usage);
    exit(0);
}
