 	program findobj

c	23.9.92   H. Hippelein
c       18.2.94   adapted to PHOTOMETRY version 
c       18.7.02   framesize increased  to 10x10k
c                 max number of objects to 200
c	19.7.02   adapted for stand-alone in OMEGA2000
c	12.12.02  changed to handle .fits files :RF
c	13.01.03  changed path of include files for CA :RF
	

c	logical		null
	byte            ichar
	character*1	type
	character*8     label
	character*64	table
	character*64	in_ima
	character*80    message, cnt, idt
	integer*4	knu, kun, dnu, dun, nelem, belem, stat
	integer*4	np(2), lin_pix(2)
	integer		ex_stat, iav, n_stars(3), m_field(5)
	real*4		frstat(3), sav_off(2), xymove(6), psf(4)
	real*4		x(200), y(200), f(200), wida(200), widb(200)
	real*4		peak(200), back(200), pa(200), xy(2), z(3), mag(200)
	real*8		start(2), step(2), radec(6)

	include 'st_def.inc'
	include 'st_dat.inc'
	
	common	/vmr	/madrid(1)

	PARAMETER	(MAX=1000)
	data             kun, knu, dun, dnu / 0, 0, 0, 0 /
	data             m_field / 0, 0, 0, 0, 0 /
	data             sav_off / 0., 0. /
	data             xymove / 0.,0.,0.,0.,0.,0. /
	data             radec / 0.,0.,0.,0.,0.,0. /

	call sts pro ('findobj')

c	message = 'Operation not permitted !'
c	message = 'Cannot read input parameters !'

c	get catalog

	call stkrdc ('in_ima', 1,1,64, iav, in_ima, kun, knu, stat)
        do l = 1, len(in_ima)
            if (in_ima(l:l) .ne. ' ') l_list = l
        end do
 

c       get desired number of objects (default=10)

	call stkrdi ('n_obj', 1,1, iav, nobj0, kun, knu, stat)

c	loop over input frames

	    call stiget (in_ima, d_r4_format, f_i_mode, f_ima_type,
	1	2, naxis, np, start, step, idt, cnt, ic, nc, stat)
	    np1 = np(1)
	    np2 = np(2)

c	statistics of input frame

	nelem = 0
        call stdfnd (nc, 'fr_stat', type, nelem, belem, stat)
	    if  (type .ne. 'R') then
		call statis (madrid(ic), np1, np2, frstat)
	        call stdwrr (nc, 'fr_stat', frstat, 1,3, dun, stat)
	    else 
c	        call stdrdr (nc, 'fr_stat', 1,3,3, frstat, dun, dnu, stat)
c	        if  (frstat(3) .eq. 0.) then
		call statis (madrid(ic), np1, np2, frstat)
	        call stdwrr (nc, 'fr_stat', frstat, 1,3, dun, stat)
c		end if
	    end if
	    type *, ' frstat =', frstat
		   
c	find objects above threshold

c	    nobj0 = 10
c	    type *, nobj0
	    thresh = 10. * frstat(3) * np1/1000
201	    nobj = nobj0
	    call detect (madrid(ic), np1, np2, frstat, thresh, nobj, x, y, f)

	    type *, nobj, ' objects found at threshold', thresh

c	analyse size and shape

	    call analyse (madrid(ic), np1, np2, nobj, x, y, f, seeing, 
	1		  frstat(3), peak, back, wida, widb, pa)
c	    if (nobj .lt. 10) then
c		n0 = 1.5 * n0
c		if (n0 .lt. 80) goto  201
c	    end if

	    psf(1) = 2.
	    psf(2) = seeing * 1.414
	    psf(3) = seeing * 1.414
	    psf(4) = 0.
	   call stdwrr (nc, 'psf_find', psf, 1,4, dun, stat)

c	create output table

	    table = in_ima

c	changed to make compatible with .fits files
	    do i = 1, 60
	       if (table(i:i+3) .eq. ".bdf") table(i:i+3) = ".tbl"
	       if (table(i:i+4) .eq. ".fits") table(i+5:i+8) = ".tbl"
	    end do

	    idt = in_ima
	    n_stars(1) = 1
	    n_stars(2) = nobj
	    n_stars(3) = 0
	    call xy_crea (table, nobj, idt, n_stars, m_field, sav_off,
	1	          start, step, radec, xymove, nt, stat)
	    if (stat .ne. 0)  goto 9999

	    type *, table
	    type 1110
 1110	    format ('  i   xpos    ypos    mag    fwhma  fwhmb    PA',
	1	    '    type')

c       transform to world coordinates and write in table

	    do iobj = 1, nobj
	        ichar = '1'
		mag(iobj) = 30. - 2.5 * log10(f(iobj))
	    	x(iobj) = (x(iobj)-1) * step(1) + start(1)
	    	y(iobj) = (y(iobj)-1) * step(2) + start(2)
		label = 's'
		size = sqrt (abs(wida(iobj) * widb(iobj)))
		if (size .le. 0.8 * seeing) label = 'n'
		if (size .gt. 1.15 * seeing) label = 'n'
		wida(iobj) = abs (2.35 * wida(iobj) * step(1))
		widb(iobj) = abs (2.35 * widb(iobj) * step(2))
	        type 1111, iobj, x(iobj), y(iobj), mag(iobj), 
	1	    wida(iobj), widb(iobj), pa(iobj), label
 1111	        format (i3, 2F8.2, F8.2, 2F7.2, F8.2, a11)
		xy(1) = x(iobj)
		xy(2) = y(iobj)
		z(1) = mag(iobj)
	        z(2) = sqrt(wida(iobj)*widb(iobj))
		z(3) = pa(iobj)
		zero = 0.
	    
		call xy_writ (nt, iobj, iobj, ichar, label, lin_pix, xy, 
	1		      zero, z, stat)
	   end do

	   see_fwhm = seeing * 2.35
	   see_arcs = abs (see_fwhm * 0.45)
	   type 915, see_fwhm, see_arcs
 915	   format (' seeing(fwhm) =', F5.2, ' pixels = ', F5.2, ' arcs')

	   call tbtclo (nt, stat)
	   call stfclo (nc, stat)


9999	ex_stat = stat
c	call epilog_px ('findobj', ex_stat, message)
	call stsepi
	call exit ()

	end



	subroutine statis (a, np1, np2, frstat)

	real*4		a(np1,np2), frstat(3), n(100000)
	real*8          sa, saa, amean

	i1 = max (1, np1/5)
	i2 = np1 - i1 + 1
	ii = max (8, np1/100)
	j1 = max (1, np2/5)
	j2 = np2 - j1 + 1
	jj = max (8, np2/100)

c       find extremes and total dispersion

	nsamp = 0
	amin = a(i1,j1)
	amax = a(i1,j1)
	sa = 0.
	saa = 0.
	do  j = j1, j2, jj
	    do  i = i1, i2, ii
	        aij = a(i,j)
                if (aij .ne. 0.0) then 
	            if (aij .lt. amin) amin = aij
	            if (aij .gt. amax) amax = aij
		    nsamp = nsamp + 1
		    sa = sa + aij
		    saa = saa + aij * aij
		end if
	    end do
	end do
	amean = sa / nsamp
	arg = saa - amean*sa
	frstat(2) = sqrt((saa - amean*sa) / nsamp)
	amin = amean - 2. * frstat(2)
	amax = amean + 2. * frstat(2)

c       make histogram

	nbin = nsamp
	do  ibin = 1, nbin
	    n(ibin) = 0
	end do
	binsize = (amax - amin) / nbin
	do  j = j1, j2, jj
	    do  i = i1, i2, ii
                if (a(i,j) .ne. 0.0) then
	            ibin = (a(i,j) - amin) / binsize
		    ibin = max (1, min (ibin, nbin))
		    n(ibin) = n(ibin) + 1
                end if
	    end do
	end do

c       find 5% bin and 95% bin

	n05 = 0.05 * nsamp
	n95 = 0.95 * nsamp
	nsum = 0
	do  ibin = 1, nbin
	    nsum = nsum + n(ibin)
	    if (nsum .le. n05) imin = ibin
	    if (nsum .le. n95) imax = ibin + 3
	end do
	bmin = amin + imin * binsize
	bmax = amin + imax * binsize

c       find mean and scatter within new borders

	na = 0
	sa = 0.
	saa = 0.
	do  j = j1, j2, jj
	    do  i = i1, i2, ii
	        aij = a(i,j)
	        if (aij .gt. bmin .and. aij .lt. bmax .and. aij .ne. 0.0) then
		    na = na + 1
		    sa = sa + aij
		    saa = saa + aij * aij
		end if
	    end do
	end do
	amean = sa / na
	frstat(1) = amean
	frstat(3) = sqrt((saa - amean*sa) / na)
c	frstat(3) = 1.2 * frstat(3)

	return
	end





	subroutine detect (b, n1, n2, frstat, thresh, nobj, x, y, f)

	real*4		b(n1,n2), h(10000,10000)
	real*4		frstat(3), x(1000), y(1000), f(1000), ia(10000), ib(10000)
	logical		edge

c	type *, 'run detect'
	nmin = 0.8 * nobj
	nmax = min (1.3 * nobj, 200)

10	continue
	type *, 'thresh =', thresh
	do  j = 1, n2
c	   type *, n2, j
	    do  i = 1, n1
		h(i,j) = b(i,j) - frstat(1)
	    end do
	end do

	iobj = 0
	do  j0 = 1, n2
c	    type *, j0, iobj
	    do  i0 = 1, n1
		edge = .false.
		if (h(i0, j0) .ge. thresh) then
		    i = i0
		    j = j0
		    id = 1
		    jd = 0
		    jmax = j
		    ia(j) = i
		    ib(j) = i
90		    i1 = i + id
		    j1 = j + jd
		    i2 = i1 + jd
		    j2 = j1 - id
c		    type *, iobj, i, j, i1, j1, i2, j2
		    if (h(i1,j1) .lt. thresh .or. i1 .lt. 1 .or.
	1		i1 .gt. n1 .or. j1 .lt. 1 .or. j1 .gt. n2) then
			ihelp = id
			id = -jd
			jd = ihelp
		    else if (h(i2,j2) .lt. thresh .or. i2 .lt. 1 .or. 
	1		i2 .gt. n1 .or. j2 .lt. 1 .or. j2 .gt. n2) then
			i = i1
			j = j1
		    else
			i = i2
			j = j2
			ihelp = id
			id = jd
			jd = -ihelp
		    end if
		    if (j .gt. jmax) then
			ia(j) = i
			ib(j) = i
		    else
			ia(j) = min(i, ia(j))
			ib(j) = max(i, ib(j))
		    end if
		    jmax = max(j, jmax)

		    if (i .eq. 1 .or. i .eq. n1 .or.
	1		j .eq. 1 .or. j .eq. n2) edge = .true.

c	check if isophote back at origin

		    if (i .ne. i0 .or. j .ne. j0 .or.
	1		id .ne. 0 .or. jd .ne. -1) goto 90

c	accumulate object pixels and set values of help frame to 0

		    npix = 0
		    sh = 0.
		    sih = 0
		    sjh = 0
		    do  j = j0, jmax
			do  i = ia(j), ib(j)
			    if (h(i,j) .gt. thresh) then
				npix = npix + 1
				sih = sih + i * h(i,j)
				sjh = sjh + j * h(i,j)
				sh = sh + h(i,j)
				h(i,j) = 0.
			    end if
			end do
		    end do

c 	center of mass, only when npix > 10

		    if (npix .le. 10 .or. edge .eq. .true.) goto 95
		    iobj = iobj + 1

c	if nobj > nmax try with higher threshold

		    if (iobj .gt. nmax) then 
	    		thresh = 1.5 * thresh
	    		goto 10
		    end if

		    x(iobj) = float(sih) / sh
		    y(iobj) = float(sjh) / sh
		    f(iobj) = sh
		end if
95	    end do
	end do

c	if nobj < nmin retry with lower threshold

	if (iobj .lt. nmin) then
	   type *, 'threshold, nobj:', thresh, iobj
	   thresh = 0.6 * thresh 
	   goto 10
	end if
	nobj = iobj

c	type *, ' nobj = ', iobj, ' threshold = ', thresh

	return 
	end



	subroutine analyse (a, n1, n2, nobj, x, y, f, seeing, rms, peak, 
	1		back, wida, widb, pa)

	integer*4       mobj(20)
	real*4		a(n1,n2), x(200), y(200), f(200), peak(200)
	real*4		back(200), wida(200), widb(200), pa(200)
	real*4          size(200), kappa

c	determine center of mass, size and PA by momentum analysis

	g = 0.01745

c       2 iterartions: first a seeing FWHM of 5 pixels is used to 
c       set apertures. Then the seeing is approximated from the 
c       median of the object sizes. In the second iteration the 
c       estimated seeing is used for the apertures. 

	seeing = 5 / 2.35
	do m = 1, 2

c	determine background in ring aperture

	   ri = 2.7 * seeing
	   ro = 4.0 * seeing
	   ri = 3.0 * seeing
	   ro = 4.0 * seeing
	   ri2 = ri**2
	   ro2 = ro**2
	   do iobj = 1, nobj
	      itera = 0
 20	      itera = itera + 1
	      i1 = max (1,  min (nint (x(iobj) - ro), n1))
	      i2 = min (n1, max (nint (x(iobj) + ro), 1))
	      j1 = max (1,  min (nint (y(iobj) - ro), n2))
	      j2 = min (n2, max (nint (y(iobj) + ro), 1))
c             type *, x(iobj), y(iobj), i1,i2,j1,j2
	      sigm = 1000000.
	      do  iter = 1, 1
		 n = 0
		 sa = 0.
		 saa = 0.
		 do j = j1, j2
		    do i = i1, i2
		       r2 = (i-x(iobj))**2 + (j-y(iobj))**2
		       if (r2 .lt. ri2 .or. r2 .gt. ro2) goto 10
		       aij = a(i,j)
		       if (abs(aij-aver) .gt. 2. * sigm) goto 10
		       n = n + 1
		       sa = sa + aij
		       saa= saa + aij**2
 10		    end do
		 end do
		 aver = sa / n
		 sigm = sqrt ((saa - aver*sa) / n)
	      end do
	      back(iobj) = aver

c	analyse object

	      rho = 2.5 * seeing
	      rho = 3.0 * seeing
	      rho2 = rho**2
	      i1 = max (1,  min (nint (x(iobj) - rho), n1))
	      i2 = min (n1, max (nint (x(iobj) + rho), 1))
	      j1 = max (1,  min (nint (y(iobj) - rho), n2))
	      j2 = min (n2, max (nint (y(iobj) + rho), 1))

	      na = 0
	      sa = 0.
	      sai = 0.
	      saj = 0.
	      do j = j2, j1, -1
		 dy = j - y(iobj)
		 do i = i1, i2
		    dx = i - x(iobj)
		    if ((dx**2 + dy**2) .gt. rho2) goto 110
		    na = na + 1
		    aij = a(i,j) - back(iobj)
		    sa = sa + aij
		    sai = sai + aij * i
		    saj = saj + aij * j
 110		 end do
	      end do
	    
	      f(iobj) = sa
	      cnti = sai / sa
	      cntj = saj / sa
	      d = sqrt ((cnti-x(iobj))**2 + (cntj-y(iobj))**2)
	      if (d .le. .7*seeing .or. itera .eq. 1 .and. d .le. seeing) then
		 x(iobj) = cnti
		 y(iobj) = cntj
	      else
		 if (itera .ge. 5) then
		    f(iobj) = 0.
		    goto 900
		 end if
		 x(iobj) = x(iobj) + .2 * (cnti-x(iobj)) / abs(cnti-x(iobj))
		 y(iobj) = y(iobj) + .2 * (cntj-y(iobj)) / abs(cntj-y(iobj))
		 goto 20
	      end if
c	      type *, iobj, x(iobj), y(iobj)

c	elongation, pa,  width

	      n = 0
	      sa = 0.
	      saii = 0.
	      sajj = 0.
	      sakk = 0.
	      do j = j1, j2
		 dy = (j - y(iobj))
		 dyy = dy * dy
		 do i = i1, i2
		    dx = (i - x(iobj))
		    dxx = dx * dx
		    if ((dxx + dyy) .gt. rho2) goto 120
		    aij = a(i,j) - back(iobj)
		    n = n + 1
		    sa = sa + aij
		    saii = saii + aij * dxx
		    sajj = sajj + aij * dyy
		    sakk = sakk + aij * dx * dy
 120		 end do
	      end do
	      saii = max(saii, 0.001)
	      sajj = max(sajj, 0.001)
	      saiijj = saii - sajj
	      if (abs(saiijj) .le. 0.001) saiijj = 0.00001

	      size(iobj) = sqrt ((saii + sajj) / 2. / sa) 
	      elo = sqrt ((saii-sajj)**2 + 4.*sakk**2) / (saii + sajj)
	      wida(iobj) = size(iobj) * (1. + elo)
	      widb(iobj) = size(iobj) * (1. - elo)
	      if (widb(iobj) .lt. 0.) widb(iobj) = -widb(iobj)

	      theta = atan2 (2. * sakk, saiijj) / 2. / g
	      pa(iobj) = theta - 90. 
	      if (pa(iobj). lt. 0.) pa(iobj) = pa(iobj) + 180.
c	      type *, iobj, wida(iobj), widb(iobj), size(iobj)

 900	   end do

c	   if (m .eq. 2) goto 950

c       determine seeing

c       find start value from 
c       histogram in objects per seeing bin of 0.2

	   do ibin = 1, 20
	      mobj(ibin) = 0
	   end do
	   do iobj = 1, nobj
	      ibin = size(iobj) / 0.2
	      mobj(ibin) = mobj(ibin) + 1
c	      ibin = (size(iobj)+0.2) / 0.2
c	      mobj(ibin) = mobj(ibin) + 1
	   end do
	   mmax = 0
	   do ibin = 1, 20
	      if (mobj(ibin) .gt. mmax) then 
		 mmax = mobj(ibin)
		 mbin = ibin
		 seeing = 0.2 * ibin - 0.1 
	      end if
	   end do
c	   type *, (size(i), i = 1, nobj)
c	   type *, mobj, mbin, seeing
	   sigs = 0.2
	   kappa = 1.6
	   do iter = 1, 10
 901	      n = 0
	      ss = 0.
	      sss = 0.
	      sigkap = kappa * sigs
	      do i = 1, nobj
		 if (abs(size(i) - seeing) .lt. kappa * sigs) then
		    n = n + 1
		    ss = ss + size(i)
		    sss = sss + size(i)**2
		 end if
 910	      end do
	      if (n .lt. 5) then
                sigs = 1.5 * sigs 
		 goto 901
	      end if
	      seeing = ss / n
c	      type *, iter, n, seeing, sigs
	      sigs = sqrt ((sss - seeing * ss) / n)
	   end do
	   see_fwhm = seeing * 2.35

 950	   continue
	end do

c	goto 960
	do iobj = 1, nobj
	   siglev = (f(iobj) / rms) * (min (widb(iobj) / seeing, 1.))**2
c          type 1111, iobj, x(iobj), y(iobj), f(iobj), 
c	1	wida(iobj), widb(iobj), pa(iobj), siglev
 1111	   format (i3, 2F8.2, F10.0, 2F7.2, 2F8.2)
	end do

c       reject cosmics and very faint objects

960	do iobj = 1, nobj
	   siglev = (f(iobj) / rms) * (min(wida(iobj) / seeing, 1.))**2
	   if (siglev .lt. 50.) then
c             type 1111, iobj, x(iobj), y(iobj), f(iobj), wida(iobj), 
c	              1widb(iobj), pa(iobj)
	      do jobj = iobj, nobj-1
		 x(jobj) = x(jobj+1)
		 y(jobj) = y(jobj+1)
		 f(jobj) = f(jobj+1)
		 back(jobj) = back(jobj+1)
		 size(jobj) = size(jobj+1)
		 wida(jobj) = wida(jobj+1)
		 widb(jobj) = widb(jobj+1)
		 pa(jobj) = pa(jobj+1)
	      end do
	      nobj = nobj - 1
	      goto 960
	   end if
	end do

	return
	end
