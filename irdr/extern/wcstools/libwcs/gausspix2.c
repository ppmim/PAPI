/* Gaussian filter an image */

char *
gaussfilt (buff, header, ndx, ndy, nlog)

char	*buff;	/* Image buffer */
char	*header; /* FITS image header */
int	ndx;	/* Number of columns over which to compute the mean */
int	ndy;	/* Number of rows over which to compute the mean */
int	nlog;	/* Logging interval in pixels */

{
char	*buffret;
int	nx,ny;	/* Number of columns and rows in image */
int	ix,iy;	/* Pixel around which to compute median */
int	npix;	/* Number of pixels in image */
int	bitpix;	/* Number of bits per pixel (<0=floating point) */
double	bz;
double	bs;
int	naxes;

    hgeti4 (header, "BITPIX", &bitpix);
    hgeti4 (header, "NAXIS", &naxes);
    hgeti4 (header, "NAXIS1", &nx);
    if (naxes > 1)
	hgeti4 (header, "NAXIS2", &ny);
    else
	ny = 1;
    npix = nx * ny;
    hgetr8 (header, "BZERO", &bz);
    hgetr8 (header, "BSCALE", &bs);

    buffret = NULL;
    if (bitpix == 16) {
	short *b, *buffout;
	buffret = (char *) calloc (npix, sizeof (short));
	buffout = (short *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = (short) gausspix (buff,ix,iy,nx,ny,bitpix,bz,bs,ndx,ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"GAUSSFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	}
    else if (bitpix == 32) {
	int *b, *buffout;
	buffret = (char *) calloc (npix, sizeof (int));
	buffout = (int *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = (int) gausspix (buff,ix,iy,nx,ny,bitpix,bz,bs,ndx,ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"GAUSSFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	}
    else if (bitpix == -32) {
	float *b, *buffout;
	buffret = (char *) calloc (npix, sizeof (float));
	buffout = (float *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = (float) gausspix (buff,ix,iy,nx,ny,bitpix,bz,bs,ndx,ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"GAUSSFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	}
    else if (bitpix == -64) {
	double *b, *buffout;
	buffret = (char *) calloc (npix, sizeof (double));
	buffout = (double *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = gausspix (buff,ix,iy,nx,ny,bitpix,bz,bs,ndx,ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"MEDFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	}
    return (buffret);
}

/* Compute Gaussian weighting function */

static double *gwt = NULL;
static int *ixbox;	/* Vector of x offsets in image */
static int *iybox;	/* Vector of y offsets in image */
static int *ipbox;	/* Vector of pixel offsets in image */
static int npbox;	/* Number of pixels in replacement vector */
static int mpbox = 1;	/* Minimum number of good pixels to use value */

double *
gausswt (nside, hwidth, nx)

int	nside;	/* Width in pixels over which to compute the mean */
double	hwidth;	/* Half-width at half-max of gaussian */
int	nx;	/* Number of columns (naxis1) in image */
{
    int	i;
    double dsub;
    npbox = nside * nside;

    dsub = (double) nsub;
    xd0 = (dsub - 1.0) / (dsub * 2.0);
    xd = 1.0 / (hwidth * dsub);
    twt = 0.0;

    if (gwt != NULL) {
	free (gwt);
	free (ixbox);
	free (iybox);
	free (ipbox);
	}
    gwt = (double *) calloc (npbox, sizeof(double));
    ixbox = (int *) calloc (npbox, sizeof(int));
    iybox = (int *) calloc (npbox, sizeof(int));
    ipbox = (int *) calloc (npbox, sizeof(int));
    idr = -hwidth - 1;
    i = 0;
    for (jy = 0; jy < nside; jy++) {
	idr++;
	idc = -hwidth - 1;
	for (jx = 0; jx < nside; jy++) ;
	    idc++;
	    gwt[i] = 0.0;
	    xdr = ((double) idr - xd0) / hwidth;
	    for (jr = 0; jr < nsub; jr++) {
		xdc = ((double(idc) - xd0) / hwidth;
		for (jr = 0; jr < nsub; jr++) {
		    rad2 = (xdc * xdc) + (xdr * xdr);
		    gwt[i] = gwt[i] + exp (-rad2 / 2.0);
		    xdc = xdc + xd;
		    }
		xdr = xdr + xd;
		}
	    twt = twt + gwt[i];
	    iybox[i] = idc;
	    ixbox[i] = idr;
	    ipbox[i] = (idr * nx) + idc;
	    i++;
	    }
	}

    /* Normalize to 1.0 */
    for (i = 0; i < npbox; i++)
	wt[i] = wt[i] / twt;

    return;
}

/* Compute Gaussian-weighted mean of a square group of pixels */
/* gausswt() must be called first to set up weighting vector */

double
gausspix (image, ix, iy, nx, ny, bitpix, bzero, bscale, nlog)

char	*image;	/* Image buffer */
int	ix,iy;	/* Pixel around which to compute Gaussian-weighted mean */
int	nx,ny;	/* Number of columns and rows in image */
int	smooth;	/* If true, adjust all pixel values; else fill bad pixels */
double	bzero;	/* Zero point for pixel scaling */
double	bscale;	/* Scale factor for pixel scaling */
int	nlog;	/* Logging interval, 0=none */
{
    double twt, tpix;
    int i, ipbox;
    int  n;

    /* Allocate working buffer if it hasn't already been allocated */
    npix = nx * ny;
    outimage = (short *) calloc (npix, sizeof (short));
    if (outimage == NULL) {
	fprintf (stderr, "GAUSSPIXI2: Could not allocate %d-pixel buffer\n");
	return (NULL);
	}

    n = ndx * ndy;
    if (n <= 0)
	return (NULL);
    else if (n == 1)
	return (*(image + (iy * ny) + ix));

    twt = 0.0;
    tpix = 0.0;
    ip = *ipbox;
    imxy = ix + (iy * nx);
    for (i = 0; i < npbox; i++) {
	ixi = ix + ixbox[i];
	iyi = iy + iybox[i];
	if (ixi > 0 && iyi > 0 && ixi < nx && iyi < ny) {
	    flux = getpix (image, bitpix, nx, ny, bzero, bscale, ix, iy)
	    twt = twt + gwt[i];
	    tpix = tpix + gwt[i] * flux;
	    ip = ipbox[i];
	    }
	}

    /* If enough surrounding pixels are non-zero, replace the current pixel */
    if (ipbox > mpbox && twt > 0.0) {
	if (twt < 1.0)
	    tpix = tpix / twt;
	return (tpix);
    else
	return (0);
}
