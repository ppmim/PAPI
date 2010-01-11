/* skyfilteronoff.c -- given a list of FITS files do running sky subtraction 
 *                assuming a ON-OFF mosaic */

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

#define MAXHWID 10                     /* max hwidth in frames of sky filter */

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
    int hwid, skybeg, usemask = 0;
    unsigned short *skysubimg = NULL;
    float *sky = NULL, *skyw = NULL, *fimg;

    if (argc != 6)
        usage();

    if (!strcmp(argv[4], "mask"))
        usemask = 1;
    else if (!strcmp(argv[4], "nomask"))
        usemask = 0;
    else
        usage();

    if (usemask)
        nplanes = readlist(argv[1], fn, mfn, xshift, yshift, MAXNPLANES);
    else
        nplanes = readlist(argv[1], fn, NULL, NULL, NULL, MAXNPLANES);

    if (nplanes < 1)
        eprintf("%s: no valid image planes", argv[0]);

    gainmap = readfits(argv[2], &nx, &ny, NULL, NULL);

    hwid = atoi(argv[3]);

    if (hwid > MAXHWID) {
      hwid=MAXHWID;
      printf("HALFNSKY reduced to %d \n", hwid);
    }

    /* maxwid = (nplanes+2)/4;
    if (hwid > maxwid) {
      hwid=maxwid;
      printf("HALFNSKY reduced to %d \n", hwid);
    } */

    /* if (hwid > MAXHWID || (2 * hwid + 1) > nplanes)
      eprintf("hwid %d, MAXHWID %d, nplanes %d\n", hwid, MAXHWID, nplanes); */


    for (i = 0; i < (4 * hwid); i++)  {
	    if (i<nplanes) readdata(i, usemask);
	}

    for (i = 0, skybeg = 1; i < nplanes; i+=2) {
        int j, nsky = 0, skyend = skybeg + 4 * hwid - 1;
        float avgscale = 0.0;

	if (skyend>=nplanes) skyend=nplanes-1;
	printf("Image: %d   Sky: ", i);

        for (j = skybeg; j <= skyend; j+=2) {  /* collect adjacent frame ptrs */
            if (j != i) {                             /* skip current frame */
                dbuf[nsky] = data[j];

                if (usemask)
                    wbuf[nsky] = wdata[j];

                avgscale += (scale[nsky] = bkgs[j]);
		printf (" %d ", j);
                nsky++;
            }
        }
	printf (" \n");

        avgscale /= (float) nsky;

        for (j = 0; j < nsky; j++) {
            scale[j] = avgscale - scale[j];
	}

        if (usemask) {
            sky = cube_mean(dbuf, wbuf, nsky, nx, ny, &skyw, scale, 1);
            fimg = skysub(data[i], nx, ny, bkgs[i], gainmap, sky, skyw, 
                            wdata[i], argv[5]);
        } else {
            sky = cube_median(dbuf, nsky, nx, ny, scale, 1);
            fimg = skysub_nomask(data[i], nx, ny, bkgs[i], gainmap, sky, 
                                   argv[5]);
        }

        skysubimg = shortint(fimg, nx, ny);

        writefits(outfn(fn[i]), fn[i], (char*)skysubimg, 16, nx, ny);

        free(sky);  free(skysubimg);  free(fimg);

        if (skyw != NULL)
            free(skyw);

        if (i >= 2*hwid && i<(nplanes-2*hwid-1)) {  /* move sliding window */
            freedata(i - 2*hwid , usemask); 
            freedata(i - 2*hwid +1, usemask); 
            readdata(i + 2*hwid , usemask);
            readdata(i + 2*hwid + 1, usemask);
            skybeg+=2;
        }

        if (i == (nplanes-3))   {
            readdata(i+2, usemask);
	}
    }

    return 0;
}

/* readdata: load in data for image plane i */
static void readdata(int i, int usemask)
{
    int nx, ny;

    data[i] = readfits(fn[i], &nx, &ny, &bkgs[i], &sigs[i]);  /* image plane */

    if (bkgs[i] <= 0 || sigs[i] <= 0)
        eprintf("readdata: ERR %s, bkg %f, sig %f\n", fn[i], bkgs[i], sigs[i]);

    if (usemask) {
        float *wmap = getwmap(fn[i], nx, ny, gainmap, sigs[i]);
        wdata[i] = getmask(wmap, nx, ny, mfn[i], xshift[i], yshift[i]);
        free(wmap);
    }

    fprintf(stderr, " Reading %s %f %f\n", fn[i], bkgs[i], sigs[i]);
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
    "row|col|rowcol|colrow|none\n\n"
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
    "example: skyfilter filelist gain.fits 4 mask rowcol\n\n";

    printf("%s", usage);
    exit(0);
}
