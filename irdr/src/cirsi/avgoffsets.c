/* avgoffsets.c -- average the dither offsets from the four chips */

/*
 * read text files called offsets.rXtoY.c[1-4] (X and Y are the beginning 
 * and ending run numbers of the dither set as input on the command line), 
 * take the average of the dither offset measurements from the four chips, 
 * and write the averages to offsets.rXtoY
 */

#include <stdio.h>
#include <stdlib.h>
#include "irdr.h"

#define MINW   0.4             /* min confidence to use offset measurement */
#define NCHIPS 4

static float xoff[NCHIPS], yoff[NCHIPS];       /* dither offset per chip */
static float woff[NCHIPS];                  /* offset confidence, 0 to 1 */

int main(int argc, char *argv[])
{
    int i, j, k, run, noff;
    char outfn[256], infn[256], fn[256], line[256];
    FILE *outfp, *infp;
    float x, y, w;

    if (argc != 3)
        eprintf("Usage: %s run_begin run_end\n", argv[0]);

    sprintf(outfn, "offsets.r%sto%s", argv[1], argv[2]);

    if ((outfp = fopen(outfn, "w")) == NULL)
        eprintf("%s: failed opening: %s\n", argv[0], outfn);

    for (i = atoi(argv[1]); i <= atoi(argv[2]); i++) {     /* foreach run */
        noff = 0;       /* collect offset measurement per chip for this run */

        for (j = 1; j <= NCHIPS; j++) {
            sprintf(infn, "offsets.r%sto%s.c%d", argv[1], argv[2], j);

            if ((infp = fopen(infn, "r")) == NULL) {
                fprintf(stderr, "%s: warn: open fail: %s\n", argv[0], infn);
                continue;
            }

            while (fgets(line, sizeof(line), infp) != NULL) {
                if (sscanf(line, "%s %f %f %f", fn, &x, &y, &w) != 4)
                    eprintf("%s: sscanf: %s\n", argv[0], line);

                if (sscanf(fn, "./irx_%d", &run) != 1)
                    eprintf("%s: sscanf run: %s\n", argv[0], fn);
 
                if (run == i) {          /* found offsets entry for this run */
                    xoff[noff] = x;
                    yoff[noff] = y;
                    woff[noff] = w;

                    if (++noff > NCHIPS)
                        eprintf("%s: noff > NCHIPS\n", argv[0]);
                }
            }

            fclose(infp);
        }

        if (noff > 0) {                  /* write avg offsets for this run */
            int nbuf = 0;
            float bufx[NCHIPS], bufy[NCHIPS];
            float sumx = 0.0, sumy = 0.0, sumw = 0.0;

            for (k = 0; k < noff; k++) {
                if (woff[k] > MINW) {
                    sumx += woff[k] * xoff[k];
                    sumy += woff[k] * yoff[k];
                    sumw += woff[k];
                    bufx[nbuf] = xoff[k];
                    bufy[nbuf++] = yoff[k];
                }
            }

            if (sumw > 0.0) {
                sumx /= sumw;
                sumy /= sumw;
                fprintf(outfp, "%d %f %f %f %f %f\n", i, sumx, sumy, sumw,
                        stdev(bufx, nbuf), stdev(bufy, nbuf));
            } else {
                fprintf(stderr, "%s: WARN offsets fail, run %d\n", argv[0], i);
            }
        }
    }

    return 0;
}
