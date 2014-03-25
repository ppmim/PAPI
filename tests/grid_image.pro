;+
; NAME:
;       GRID_IMAGE
;
; PURPOSE:
;   	Creates the image of a regular grid of points or segments.
;		The grid can be flat or distorted (barrel or pincushion) by a cubic distortion.
;		Mainly useful to test distortion programs.
;
; CALLING SEQUENCE:
;
;		Result = GRID_IMAGE(sizex, sizey, nx, ny, border, DISTORTION=distortion, FWHM=fwhm, $
;						TRASL=trasl, ROT=rot, CENTRAL=central, SQUARES=squares, $
;						HIGH_PRECISION=HIGH_PRECISION, [xdef, ydef], [xs, ys])
;
; INPUTS:
;       sizex:     x size of the final image
;
;       sizey:     y size of the final image
;
;       nx, ny:    number of grid points in x and y (greater or equal 2)
;
;       border:    optional border to leave around the grid (pixels)
;
;
; KEYWORDS:
;
;       DISTORTION: distortion at the image corner:
;			      (if < 0 you have a barrel distortion, if > 0 a pincushion distortion,
;                 if = 0 there is no distortion).
;				  Distortion is computed respect to the image center, or optical axis (see AXIS keyword).
;				  NOTE: image center is defined as (lato/2.0-0.5): i.e. it is the center of the central pixel for odd images
;				  and the corner between the 4 central pixels for even images.
;
;		AXIS:     two element vector with the (x,y) coordinate of the optical axis (i.e. center of distortion).
;				  If not set, the image center is taken as center of distortion.

;       FWHM:     full width at half maximum of the gaussian PSF to convolve with the grid points
;				  NOTE: when FWHM is used and bordo=0, the centers of edge points are at the edges of
;				  the pixels areas, instead when FWHM=0 the single points are centered at the centers of
;				  the pixels at the border of the image. This means that in the two cases the centers xdef.ydef
;				  are NOT the same! In fact, in the FHWM=0 case border is set to 0.5 to have points at the centers
;				  of the edge pixels.
;				  NOTE: if FWHM=0 the points position is discretized to the pixels so that we suggest to always use
;				  a non zero FWHM (it can also be fractional if small points are needed).
;
;       TRASL:    2-element vector of the form [offx, offy] with the x and y translation of optical
;				  axis, i.e. the center of the distortion.
;				  NOTE: the order is first Rotation, then Translation.
;
;       ROT:	  angle of rotation around image center, expressed in degrees, counter-clockwise.
;				  NOTE: the order is first Rotation, then Translation.
;
;       CENTRAL:  if set, the central point of the grid is enhanced (value = 2.0 instead of 1.0)
;
;		SQUARES:  Set this keyword to produce an image of a reticle instead of single points.
;				  In this case the fwhm value is the lines width.
;				  Attention: the grid nodes will be connected by straight segments, i.e. the nodes
;				  positions follow the distortion, but the segments are not distorted.
;
;		HIGH_PRECISION: integer value between 2 and 10: if provided, the image of the grid
;					is computed with the given sub-pixel precision; only useful if FWHM is set.
;					The odd greater or equal to this number will be used.
;					NOTE: High value of this keyword imply that the computation can be very slow.
;
; OUTPUTS:
;
;       The function returns a float 2-D array of sizex times sizey with value 0 but in the points,
;		which value is 1.0.
;		If FHHM is not zero then the total energy of each point is 1.0
;
; OPTIONAL OUTPUTS:
;
;       xdef,ydef:	x and y coordinates of the final grid (pixels)
;
;       xs,ys:		x and y coordinates of the grid without distortion and translation (pixels)
;
; PROCEDURE:
;
;		A lens optical distortion is a radial deformation which can be written in power series with
;		odd exponents, starting with r^3. Here only this cubic term is considered for the distortion.
;
; CALLS:
;
;		The function uses FSHIFT and SEGMENT, plus Astrolib library.
;
; MODIFICATION HISTORY:
;       February 2004, G. Li Causi, Rome Astronomical Observatory
;       August 2004, G. Li Causi - Added ROT keyword
;		July 2006, G.Li Causi - correct center of image definition and use bordo=0.5 when FWHM=0.
;		August 2008, G.Li Causi - make first Rotation, then Translation; added AXIS keyword.
;
; 		Planned Improvements:
;		- implement the PSF also near the image edges.
;		- Restore DEFORM keyword to also distort the gaussian PSF.
;-


FUNCTION GRID_IMAGE, latox, latoy, nx, ny, bordo_in, DISTORTION=distortion, FWHM=fwhm, DEFORM=deform, $
		TRASL=trasl, ROT=rot, AXIS=axis, CENTRAL=central, SQUARES=squares, HIGH_PRECISION=HIGH_PRECISION, xdef,ydef,xs,ys


ON_ERROR, 2

;****************
;PARAMETERS CHECK
;****************
IF NOT KEYWORD_SET(FWHM) THEN bordo = bordo_in + .5 ELSE bordo = bordo_in	;+.5 if single points, in order to have them at center of edge pixels

IF NOT KEYWORD_SET(FWHM) AND (latox/2. NE latox/2 OR latoy/2. NE latoy/2) THEN print, 'WARNING: FHWM=0 used with even-sized images yields sampling errors on points position...'

IF latox LE 0 OR latoy LE 0 THEN Message, 'LATOX and LATOY must be greater than 0'

IF nx GT (latox-2*bordo) OR ny GT (latoy -2*bordo) OR nx LT 2 OR ny LT 2 THEN $
								Message, 'NX and NY must be > 2 and < (LATO-2*BORDO)'


IF keyword_set(HIGH_PRECISION) THEN HIGH_PRECISION=odd_ge(FIX(HIGH_PRECISION)) > 1

;**********
;GRID SETUP
;**********
lx2 = latox/2. - .5			;central pixel, the -.5 is to center both for odd and even images
ly2 = latoy/2. - .5

passo_x = double(latox - 2 * bordo) / (nx-1)		;grid gap
passo_y = double(latoy - 2 * bordo) / (ny-1)

x = dindgen(nx) * passo_x + bordo - .5		;the -.5 is to have points at the very edges of the image when bordo=0.
y = dindgen(ny) * passo_y + bordo - .5

nstars = nx * ny			;number of grid points or reticle nodes
xs = dblarr(nstars)
ys = dblarr(nstars)


;Reforming points coordinates:
nn = 0
FOR n1 = 0, ny-1 DO BEGIN
	FOR n2 = 0, nx-1 DO BEGIN
  		xs[nn] = x[n2]
 		ys[nn] = y[n1]
 		nn = nn + 1
  ENDFOR
ENDFOR


;*************
;GRID ROTATION
;*************
IF KEYWORD_SET(ROT) THEN BEGIN		;rotazione
	th = rot * !DTOR
	FOR n = 0, nn-1 DO BEGIN
  		dx = xs[n] - lx2
  		dy = ys[n] - ly2
  		ro = SQRT(dx^2 + dy^2)
  		angle = ATAN(dy, dx) + th
  		xs[n] = ro * cos(angle) + lx2
 		ys[n] = ro * sin(angle) + ly2
	ENDFOR
ENDIF


;****************
;GRID TRANSLATION
;****************
IF KEYWORD_SET(TRASL) THEN BEGIN		;traslazione
	sx = double(trasl[0])
	sy = double(trasl[1])
ENDIF ELSE BEGIN
	sx = 0d
	sy = 0d
ENDELSE

xs = xs + sx
ys = ys + sy



;******************
;IMAGE CONSTRUCTION
;******************
;img = fltarr(latox,latoy)
;FOR i=0, nstars-1 DO BEGIN			;set fluxes
;	IF xs[i] LE latox-1 AND  ys[i] LE latoy-1 AND xs[i] GE 0 AND ys[i] GE 0 THEN BEGIN
;		img[xs[i], ys[i]] = 1.0
;	ENDIF
;ENDFOR


img = readfits('/tmp/prueba1.fits')


nc = (ny/2) * nx + nx/2
IF KEYWORD_SET(CENTRAL) AND xs[nc] LE latox-1 AND  ys[nc] LE latoy-1 AND $
		xs[nc] GT 0 AND ys[nc] GT 0 THEN img[xs[nc],ys[nc]]=2.0		;set central flux

IF NOT KEYWORD_SET(DISTORTION) AND NOT KEYWORD_SET(DEFORM) AND NOT keyword_set(FWHM) AND NOT keyword_set(squares) THEN RETURN, img	;image of 1 pixel points


;***************
;GRID DISTORTION
;***************
IF KEYWORD_SET(DISTORTION) THEN BEGIN

	img = fltarr(latox,latoy)

	IF NOT KEYWORD_SET(AXIS) THEN BEGIN
		ax = lx2
		ay = ly2
	ENDIF ELSE BEGIN
		ax = axis[0]
		ay = axis[1]
	ENDELSE

	xs = xs - ax	;punti rispetto all'asse ottico
	ys = ys - ay

	r = sqrt(xs^2 + ys^2)
	fi = atan(ys/xs)

	wge = where(xs eq 0 and ys ge 0)
	wlt = where(xs eq 0 and ys lt 0)
	IF max(wge) GE 0 THEN fi(wge) = 90 * !dtor
	IF max(wlt) GE 0 THEN fi(wlt) = -90 * !dtor

	ng=(where(xs lt 0))
	IF max(ng) GE 0 THEN fi(ng) = fi(ng) + 180 * !dtor

	IF KEYWORD_SET(DEFORM) THEN distortion = -distortion

	;Cubic radial distortion:
	r_corner = sqrt(lx2^2 + ly2^2)
	aa = distortion / r_corner^2			;distortion = ((r_corner + aa * r_corner^3) - r_corner) / r_corner
	rdef = r + aa * r^3

	xdef = rdef * cos(fi) + ax
	ydef = rdef * sin(fi) + ay

ENDIF ELSE BEGIN

	xdef = xs
	ydef = ys

ENDELSE


;***************
;PSF CONVOLUTION
;***************
IF KEYWORD_SET(FWHM) THEN BEGIN


	;IF keyword_set(DEFORM) THEN BEGIN		;<-------ricontrollare deform
	;	;DEFORM:   if set, the PSF of the grid points are distorted with the grid.
	;	;
	;
	;	;*****************
	;	;GAUSSIAN + DEFORM
	;	;*****************
	;	; (Chiama la function REPIX con la keyword /FAST.)
	;
	;	scal = 1.
	;	frame = 1*fwhm*scal
	;	lxs = latox*scal
	;	lys = latoy*scal
	;	xss = xs*scal
	;	yss = ys*scal
	;	lx2s = lx2*scal
	;	ly2s = ly2*scal
	;
	;	llx = lxs + 2*frame
	;	lly = lys + 2*frame
	;	img_big = fltarr(llx, lly)						;bordo per la convoluzione
	;	img_big(frame:frame+lxs-1,frame:frame+lys-1) = rebin(temporary(img), lxs,lys)
	;	stella = psf(fwhm*scal, 2 * frame)
	;
	;   print, 'Convolving with gaussian...'
	;	img_big = convol(temporary(img_big), stella, /CENTER)
	;
 	;    print, 'Computing deformation...'
	;	xdef = xdef*scal + frame
	;	ydef = ydef*scal + frame
	;
	;	xx = MCS( xdef, xs+lx2, ys+ly2, BOUNDS=[0, 0, llx, lly], nx=(llx+1), ny=(lly+1) )
  	;   yy = MCS( ydef, xs+lx2, ys+ly2, BOUNDS=[0, 0, llx, lly], nx=(llx+1), ny=(lly+1) )
	;
    ;	print, 'Deforming grid...'
    ;	img_big = repix(temporary(img_big), xx,yy,FAST=0)
	;
	;	img_big = temporary(img_big(frame:frame+lxs-1,frame:frame+lys-1))			;ritaglio l'immagine finale
	;	img = rebin(temporary(img_big), latox, latoy)
	;
  	;ENDIF ELSE BEGIN

		scal = 1.
		IF KEYWORD_SET(HIGH_PRECISION) THEN scal = float(HIGH_PRECISION)		;fattore di resampling (scal=6 -> precisione nel centro = 0.0005 pixel in simulazioni)

		lxs = latox * scal
		lys = latoy * scal

		img_big = fltarr(lxs, lys)

		IF KEYWORD_SET(SQUARES) THEN BEGIN

			;*****************
			;SQUARES NO DEFORM
			;*****************

		  	FOR nyy=1, ny-1 DO BEGIN
				n1 = nyy * nx + indgen(nx)
		  		xd1 = (xdef(n1) * scal)
		  		yd1 = (ydef(n1) * scal)
				n0 = (nyy-1) * nx + indgen(nx)
		  		xd0 = (xdef(n0) * scal)
		  		yd0 = (ydef(n0) * scal)
			  	img_big = segment(img_big, xd0, yd0, xd1, yd1, 1.0)
			ENDFOR
		  	FOR nxx=1, nx-1 DO BEGIN
				n1 = indgen(ny) * nx + nxx
		  		xd1 = (xdef(n1) * scal)
		  		yd1 = (ydef(n1) * scal)
				n0 = indgen(ny) * nx + (nxx-1)
		  		xd0 = (xdef(n0) * scal)
		  		yd0 = (ydef(n0) * scal)
			  	img_big = segment(img_big, xd0, yd0, xd1, yd1, 1.0)
			ENDFOR

			fwhm_factor = 5.
			ls = odd_ge(fwhm_factor * fwhm * scal) > 3		;box per la PSF: odd size,  >3 to always have a 2D array
			pp = psf_Gaussian(NPIXEL=[ls, ls], FWHM=fwhm*scal, /NORMALIZE)
			img_big = convol(img_big, pp)


		ENDIF ELSE BEGIN

			;******************
			;GAUSSIAN NO DEFORM
			;******************

			fwhm_factor = 5.
			ls = odd_ge(fwhm_factor * fwhm * scal) > 3		;box per la PSF: odd size, >3 to always have a 2D array
	     	stella = psf_Gaussian(NPIXEL=[ls,ls], FWHM=fwhm*scal, /NORMALIZE)

		  	FOR n=0, nstars-1 DO BEGIN
		  		xds = xdef[n] * scal
		  		yds = ydef[n] * scal
		  		fxds = fix(xds)
		  		fyds = fix(yds)

		  		IF (fxds-ls/2+scal/2) GE 0 AND (fyds-ls/2+scal/2) GE 0 AND (fxds+ls/2+scal/2) LT (lxs-1) AND (fyds+ls/2+scal/2) LT (lys-1) THEN BEGIN
		  			;Aggiungo scal/2 per avere il centro della stella al centro del pixel originale di img, che e' indicizzato al suo angolo in basso a sinistra
			  		img_big(fxds-ls/2+scal/2:fxds+ls/2+scal/2, fyds-ls/2+scal/2:fyds+ls/2+scal/2) = $
			  			img_big(fxds-ls/2+scal/2:fxds+ls/2+scal/2, fyds-ls/2+scal/2:fyds+ls/2+scal/2) + fshift(stella, xds-fxds, yds-fyds)
			  		IF KEYWORD_SET(central) AND n EQ nc THEN $
				  		img_big(fxds-ls/2+scal/2:fxds+ls/2+scal/2, fyds-ls/2+scal/2:fyds+ls/2+scal/2) = $
				  			img_big(fxds-ls/2+scal/2:fxds+ls/2+scal/2, fyds-ls/2+scal/2:fyds+ls/2+scal/2) + fshift(stella, xds-fxds, yds-fyds)
				ENDIF


		  	ENDFOR

		ENDELSE

		img = rebin(temporary(img_big), latox, latoy)

  	;ENDELSE



ENDIF ELSE BEGIN		;no fwhm



	IF KEYWORD_SET(SQUARES) THEN BEGIN

		;*****************
		;SQUARES NO DEFORM
		;*****************

	  	FOR nyy=1, ny-1 DO BEGIN
			n1 = nyy * nx + indgen(nx)
	  		xd1 = (xdef(n1))
	  		yd1 = (ydef(n1))
			n0 = (nyy-1) * nx + indgen(nx)
	  		xd0 = (xdef(n0))
	  		yd0 = (ydef(n0))
		  	img = segment(img, xd0, yd0, xd1, yd1, 1.0)
		ENDFOR
	  	FOR nxx=1, nx-1 DO BEGIN
			n1 = indgen(ny) * nx + nxx
	  		xd1 = (xdef(n1))
	  		yd1 = (ydef(n1))
			n0 = indgen(ny) * nx + (nxx-1)
	  		xd0 = (xdef(n0))
	  		yd0 = (ydef(n0))
		  	img = segment(img, xd0, yd0, xd1, yd1, 1.0)
		ENDFOR

		;FOR n = 0, nstars-1 DO BEGIN
		;	IF xdef(n) LE latox-1 AND  ydef(n) LE latoy-1 AND  xdef(n) GE 0 AND $
		;			ydef(n) GE 0 THEN img(xdef(n), ydef(n)) = 1.0
		;ENDFOR

		  	img = segment(img, xd0, yd0, xd1, yd1, 1.0)
		  	img = segment(img, xd0, yd0, xd1, yd1, 1.0)

	ENDIF ELSE BEGIN

		;*************
		;SIMPLE POINTS
		;*************

		;IF keyword_set(deform) THEN BEGIN	;<------?????			;deform but no fwhm

		;	xx = MCS( xdef, xs+lx2+sx, ys+ly2+sy, BOUNDS=[0, 0, latox, latoy], nx=(latox+1), ny=(latoy+1) )
		;  	yy = MCS( ydef, xs+lx2+sx, ys+ly2+sy, BOUNDS=[0, 0, latox, latoy], nx=(latox+1), ny=(latoy+1) )

	  	;	print, 'Deforming grid...'
	  	;	img_big = repix(temporary(img), xx,yy,/FAST)

		;ENDIF ELSE BEGIN								;no deform and no fwhm

			FOR n = 0, nstars-1 DO BEGIN
				IF xdef(n) LE latox-1 AND  ydef(n) LE latoy-1 AND  xdef(n) GE 0 AND $
						ydef(n) GE 0 THEN img(xdef(n), ydef(n)) = 1.0
			ENDFOR

		;ENDELSE

	ENDELSE

ENDELSE


RETURN, img

END
