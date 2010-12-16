int ncb;	/* number of columns in image buffer array
int ncol;	/* number of columns to read
int nrow;	/* number of rows to read
int nbox;	/* radius in pixels of box
double wbox;	/* half-width at half-max of gaussian
int allpix;	/* .true. to adjust all pixels, else only bad ones
int fixpix;	/* .true. to adjust only bad pixels
int badval;	/* 1 if replacing pixels with single bad value
double fluxmin,fluxmax; /* min and max allowable values , or
			min and max disallowed values for single bad value */
double errflux;	/* delta flux for single bad value (for floating point)
double wt(225);	/* radial gaussian weights
int icbox(225);	/* offsets of box in columns
int irbox(225);	/* offsets of box in rows
int ibbox(225);	/* offsets of box in buffer
int ibufx(50000);	/* ring buffer for replacement pixels
char string[80];	/* string for output
int lbox;	/* logging interval during fill
char null;
int mpbox,npbox,nbpix,nside,i,idr,ir,idc,ic,ibr,jc,jr;
int minrow,maxrow,mincol,maxcol,npix,ibp0,ibpr,nrep,nsub;
int ibrin,ibrout,ipix,irow,icol,ibp,ibpx,ibpx0,ipbox,icx,irx;
double rad2,xdr,xdc,xd0,xd;
double flux,flux0,twt,tpix,skipval;
int diag;

    null = (char) 0;
    call pgetl (parms,'diag',diag)
    if (diag)
	printf ("Subroutine imfilt\n");

    nbox = 1;
    call pgeti4 (parms,'nbox',nbox)
    nside = 1 + (2 * nbox);
    npbox = nside * nside;
    if (nside < 2) {
	fprintf (stderr,"imfilt: Gaussian box too small, side = %d < 2\n",nside);
	return;
	}
    else if (nside > 15) {
	fprintf (stderr, "imfilt: Gaussian box too big, side = %d > 15\n,nside);
	return;
	}

/* Halfwidth in pixels of gaussian weighting */
    wbox = 1.0d0
    call pgetr8 (parms,'wbox',wbox)
    if (wbox .le. 0.d0) wbox = 1.d0

/* Image header information */
    hgeti4 (header,'NAXIS1',ncb);
    hgeti4 (header,'BITPIX',nbpix);

    skipval = -9999.0;
    call pgetr8 (parms,'skipval',skipval)

/* Logging interval */
    lbox = 0
    call pgeti4 (parms,'lbox',lbox)

    allpix = 0
    call pgetl (parms,'allpix',allpix)
    if (allpix) {
	fprintf (stderr, "Image smoothed with %d x %d %6.2f-pixel gaussian box\n",
		 nside,nside,wbox,mpbox);
	fixpix = 0;
	}
    else {
	fprintf (stderr,'image filled using %d x %d %6.2f-pixel box, %d adj\n",
		 nside,nside,wbox,mpbox);
	fixpix = 1;
	}
    if (diag) call pstring (string)
    hputs (header,"HISTORY", string);

/* Compute Gaussian weighting and offset tables (row, column, buffer) */
    if (filter == GAUSSIAN) {
	i = 0;
	idr = -nbox - 1
	xd0 = dble (nsub - 1) / dble (nsub * 2)
	xd = 1.d0 / (wbox * dble (nsub))
	twt = 0.d0
	do ir = 1, nside
	idr = idr + 1
	idc = -nbox - 1
	do ic = 1, nside
    	i = i + 1
    	idc = idc + 1
    	wt(i) = 0.d0
    	xdr = (dble (idr) - xd0) / wbox
    	do jr = 1, nsub
    	    xdc = (dble (idc)  - xd0) / wbox
    	    do jc = 1, nsub
    		rad2 = xdc*xdc + xdr*xdr
    		wt(i) = wt(i) + dexp (-rad2 / 2.d0)
    		xdc = xdc + xd
    		enddo
    	    xdr = xdr + xd
    	    enddo
    	twt = twt + wt(i)
    	icbox(i) = idc
    	irbox(i) = idr
    	ibbox(i) = (idr * ncb) + idc
    	if (diag) write(*,*) ir,ic,rad2,wt(i),ibbox(i)
    	enddo
	enddo

    for (i = 0; i < npbox; i++) {
	wt(i) = wt(i) / twt
	}

/* Get region of image to modify */
    hgeti4 (header,'naxis1',ncol);
    hgeti4 (header,'naxis2',nrow);
    if (diag) {
       fprintf (stderr, "'imbox:  columns',mincol,' to',maxcol,', rows',
     1    	      minrow,' to',maxrow
       endif
    if (nrow <= 0 || ncol <= 0) return

/* initialize buffer pointers */
    npix = nrow * ncol
    ibp0 = (minrow - 2) * ncb
    ibpr = (minrow - 2) * ncb + mincol
    nrep = 0
    ir = -nbox - 1
    ibrin = 0
    ibrout = 0
    ipix = 0

/* loop through rows of image */
    for (irow = 0; irow < nrow; irow++) {
	ibp0 = ibp0 + ncb;
	ibrin = ibrin + 1;
	if (ibrin > nside)
	    ibrin = 1;
	ibpx0 = ncb * (ibrin - 1);
/*	write(*,*) irow,ir,ibrin,ibrout,ibp0,ibpx0 */

/* Loop through columns of image */
	for (icol = 0; icol < ncol; icol++) {
	    ibp = ibp0 + icol;
	    ipix = ipix + 1;
	    call imgetpix (ibuf,ibp,flux0,nbpix)
	    ibpx = ibpx0 + icol;

	    /* Compute replacement value */
    		if ((flux0.gt.fluxmin).and.(flux0.lt.fluxmax)) {
    		    call imsetpix (ibufx,ibpx,flux0,nbpix)
    		    go to 590
    		elseif (diag) {
    		    write(*,*) 'imbox: ',irow,icol,' dropped'
    		}
    	    }
    	    }

    	twt = 0.d0
    	tpix = 0.d0
    	ipbox = 0
    	do i = 1, npbox
    	    icx = icol + icbox(i)
    	    irx = irow + irbox(i)
    	    if ((icx .ge. mincol) .and. (icx .le. maxcol) .and.
     1    		(irx .ge. minrow) .and. (irx .le. maxrow)) {
    		call imgetpix (ibuf,ibp+ibbox(i),flux,nbpix)
    		if (badval) {
    		    if ((flux.lt.fluxmin).or.(flux.gt.fluxmax)) {
    			ipbox = ipbox + 1
    			twt = twt + wt(i)
    			tpix = tpix + wt(i) * flux
    			endif
    		elseif ((flux.ge.fluxmin).and.(flux.le.fluxmax)) {
    		    ipbox = ipbox + 1
    		    twt = twt + wt(i)
    		    tpix = tpix + wt(i) * flux
c    		    write(*,*) irow,icol,irx,icx,flux,wt(i)
    		endif
    		endif
    	    enddo

/* If enough of surrounding pixels are non-zero, replace the current pixel */
    	if ((ipbox .gt. mpbox) .and. (twt .gt. 0)) {
    	    if (nbpix .gt. 0) {
    		flux = idnint (tpix)
    	    else
    		flux = tpix
    	    endif
    	    call imsetpix (ibufx,ibpx,flux,nbpix)
    	    nrep = nrep + 1
    	else
    	    call imsetpix (ibufx,ibpx,flux0,nbpix)
    	endif

/* Log progress */
590    	if ((lbox .gt. 0) .and. (mod (ipix,lbox) .eq. 0)) {
    	    write(string,595) nrep,ipix,npix
595    	    format(i7,i9,i9,'$')
    	    call pline (string)
    	    endif

    	enddo

/* Write fixed row back into image if it won't be used any more */
	ir = ir + 1
	if (ir .gt. 0) {
    	ibrout = ibrout + 1
    	if (ibrout .gt. nside) ibrout = 1
    	ibpx = ncb * (ibrout - 1) + mincol
    	ibpr = ibpr + ncb 
    	call immove (ibufx,ibpx,ibuf,ibpr,ncol,nbpix)
    	endif

	enddo

    ibrout = ibrout + 1
    do ibr = ibrout, nside
	ir = ir + 1
	ibpr = ibpr + ncb
	ibpx = ncb * (ibr - 1) + mincol
	call immove (ibufx,ibpx,ibuf,ibpr,ncol,nbpix)
c	write(*,*) irow,ir,ibr
	enddo

    write(string,*) nrep,' /',npix,' pixels replaced$'
    if ((diag) .or. (lbox .gt. 0)) call pstring (string)
    call hputs (header,'history',string)
