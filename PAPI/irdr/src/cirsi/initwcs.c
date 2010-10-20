/* initwcs.c -- set approximate WCS from RA/DEC/SCALE/NAXIS keywords */

/*
 * RA,DEC keywords in raw data indicate telescope focal plane center,
 * not chip center.  use the known chip offsets from the telescope pointing
 * position to estimate the image RA,DEC center, and call put_wcs()
 * with the image RA,DEC center to write a rough WCS for the image.
 *
 * initwcs.c has built-in defaults for the chip offsets and rotations for
 * both du Pont and INT data, but in the case that these are not correct
 * the user can create a file called "chipoffsets" containing chip_number,
 * ra_offset, dec_offset, and rotation (degrees), eg for CAHA 2.2m:
 * 1 -0.128  0.128  0.0
 * 2 -0.128 -0.128  0.0
 * 3  0.128 -0.128  0.0
 * 4  0.128  0.128  0.0
 *
 * initwcs assumes that RA,DEC keywords refer to telescope pointing
 * position, so do not run initwcs multiple times on the same image.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "irdr.h"

#define DEGRAD 0.017453293
#define NCHIP 4

/* OFFSETS (from telescope field center to chip center in degrees) */

static double CA35_ra_chip_offset [NCHIP] = {0.0, -0.054, 0.054, 0.054};
static double CA35_dec_chip_offset [NCHIP] = {0.0, -0.054, -0.054, 0.054};
static double CA22_ra_chip_offset [NCHIP] = {0.128, 0.128, -0.128, -0.128};
static double CA22_dec_chip_offset [NCHIP] = {0.128, -0.128, -0.128, 0.128};

/* ROTATIONS (chip rotations in degrees) */

static double CA35_chip_rotation [NCHIP] = {0.0, -0.4, -0.2, -0.3};
static double CA22_chip_rotation [NCHIP] = {0.2, 0.35, 0.85, 0.75};

/* precession routine from WCSTools: */

extern void fk5prec(double epoch0, double epoch1, double *ra, double *dec);

static void read_offsets(void);

int main(int argc, char *argv[])
{
    int i, nx, ny, chip = -1;
    double ra, dec, scale, posang, epoch0 = 0.0;
    char *tel=NULL;
    char *inst=NULL;

    if (argc < 2)
        eprintf("Usage: %s *.fits\n", argv[0]);

    read_offsets();
 
    for (i = 1; i < argc; i++) {
        /*get_key_int(argv[i], "CHIP", &chip); ONLY FOR DEMO !!!!!!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!*/
        chip = 1;
        
        if (chip < 1 || chip > NCHIP)
            eprintf("%s: unexpected chip number: %d\n", chip);

        if (get_key_int(argv[i], "NAXIS1", &nx) < 0)
            eprintf("%s: failed reading NAXIS1\n", argv[0]);

        if (get_key_int(argv[i], "NAXIS2", &ny) < 0)
            eprintf("%s: failed reading NAXIS2\n", argv[0]);
        

        if (get_wcs(argv[i], &ra, &dec, &scale, &posang)==-1)
            eprintf("%s: get_wcs failed\n", argv[0]);
  
        if (get_key_double(argv[i], "EQUINOX", &epoch0) < 0 || epoch0 < 1990) {
            fprintf(stderr, "WARN: fail reading EQUINOX, assume 2000.0\n");
            epoch0 = 2000.0;
        }

        if ((tel = get_key_str(argv[i], "TELESCOP")) == NULL)
            fprintf(stderr, "WARN: %s: failed reading TELESCOP keyword\n", argv[0]);
        else
        {
            if (strncmp("CA 3.5m", tel, 7) == 0) {                            /* CAHA 3.5 */
                ra  += CA35_ra_chip_offset[chip-1] / cos(DEGRAD * fabs(dec));
                dec += CA35_dec_chip_offset[chip-1];
                put_key_double(argv[i], "CHIPROT", CA35_chip_rotation[chip-1]);
    
            } else if (strncmp("CA 2.2m", tel, 7) == 0) {                     /* CAHA 2.2 */
                ra  += CA22_ra_chip_offset[chip-1] / cos(DEGRAD * fabs(dec));
                dec += CA22_dec_chip_offset[chip-1];
                put_key_double(argv[i], "CHIPROT", CA22_chip_rotation[chip-1]);
    
            } 
        }
        /* Init WCS for the instrument O2000 */
        
        if ((inst=get_key_str(argv[i], "INSTRUME")) ==  NULL)
            eprintf("[ERROR] %s: failed reading INSTRUME keyword\n", argv[0]);
            
        if (strncmp(inst,"Omega2000",9)==0) 
        {
            fk5prec(epoch0, 2000.0, &ra, &dec);

	        fprintf(stderr,"\n-->new  RA=%f, newDEC=%f\n", ra, dec);

            put_wcs(argv[i], ra, dec, scale, nx, ny);

            put_key_double(argv[i], "EQUINOX", 2000.0);
        }
        else if (strncmp(inst,"HAWKI",5)==0)
        {
            fprintf(stderr, "With HAWKI instrument no aditional WCS initilization required ...\n");
        }
        else
        /* in any case, we try ...*/
        {
            fprintf(stderr, "WARN: initwcs found a non O2000 instrument ...\n");
            fk5prec(epoch0, 2000.0, &ra, &dec);

            fprintf(stderr,"\n-->new  RA=%f, newDEC=%f\n", ra, dec);

            put_wcs(argv[i], ra, dec, scale, nx, ny);

            put_key_double(argv[i], "EQUINOX", 2000.0);
        
        }
    }
    return 0;
}

/* read chip offsets and rotations from the file "chipoffsets" if it exists */
static void read_offsets(void)
{
    int chipno;
    float dra, ddec, rot;
    char line[256], *fn = "chipoffsets";
    FILE *fp;

    if ((fp = fopen(fn, "r")) == NULL)
        return;

    fprintf(stderr, "reading chip offsets...\n");

    while (fgets(line, sizeof(line), fp) != NULL) {
        if (sscanf(line, "%d %f %f %f", &chipno, &dra, &ddec, &rot) != 4)
            eprintf("read_offsets: check list format\n");

        if (chipno < 1 || chipno > NCHIP)
            eprintf("read_offsets: bad chip no.: %d\n", chipno);

        /* just overwrite both sets */
        CA35_ra_chip_offset[chipno-1] = CA22_ra_chip_offset[chipno-1] = dra;
        CA35_dec_chip_offset[chipno-1] = CA22_dec_chip_offset[chipno-1] = ddec;
        CA35_chip_rotation[chipno-1] = CA22_chip_rotation[chipno-1] = rot;
    }
}
