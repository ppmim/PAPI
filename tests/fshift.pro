;+
; NAME:
;       FSHIFT
;
; PURPOSE:
;       Shifts an image of fractional offsets.
;
; CATEGORY:
;       IMAGE PROCESSING.
;
; CALLING SEQUENCE:
;
;       Result = FSHIFT(array, offx, offy [, EDGE_WRAP=edge_wrap] [,OUT=out])
;
; INPUTS:
;       array:  A 2 dimensional image to shift.
;
;       offx:   x offset
;
;       offy:   y offset
;
; KEYWORDS:
;       EDGE_WRAP:	If set, then the outside border is taken from the opposite border.
;
;       OUT:	Value for the border if neither EDGE_WRAP or EDGE_TRUNCATE are set.
;
; OUTPUTS:
;       Rresult: The image shifted of offx,offy.
;
; PROCEDURE:
;      	A convolution with a 2x1 matrix is used for x and then for y shift.
;		This procedure does conserve the flux.
;
; MODIFICATION HISTORY:
;       Feb 2004 - 	Gianluca Li Causi, INAF - Rome Astronomical Observatory
;					licausi@mporzio.astro.it
;					http://www.mporzio.astro.it/~licausi/
;-

FUNCTION FSHIFT, im, offx, offy, EDGE_WRAP=edge_wrap, out=out

on_error, 2

im = reform(im)
s = size(im)

lx = s(1) - 1
ly = s(2) - 1

imx_s = float(im)

iox = fix(offx)					;parte intera
ioy = fix(offy)
fracx = offx - iox				;frazione ulteriore
fracy = offy - ioy

IF NOT KEYWORD_SET(OUT) THEN out = 0.

;Parte intera dello shift
IF abs(iox) LE lx AND abs(ioy) LE ly THEN imx_s = shift(imx_s, iox, ioy)

imbig = fltarr(s[1]+2, s[2]+2) + out
imbig[1:s[1], 1:s[2]] = imx_s

IF KEYWORD_SET(EDGE_WRAP) THEN BEGIN
	IF offx gt 0 THEN imbig(0, *) = imbig(lx+1, *)
	IF offx lt 0 THEN imbig(lx+2, *) = imbig(1, *)
	IF offy gt 0 THEN imbig(*,0) = imbig(*, ly+1)
	IF offy lt 0 THEN imbig(*,ly+2) = imbig(*, 1)
ENDIF

IF NOT KEYWORD_SET(EDGE_WRAP) THEN BEGIN
	if iox gt 0 then imbig(0:iox, *) = out
	if iox lt 0 then imbig(2+lx+iox:1+lx, *) = out
	if ioy gt 0 then imbig(*,0:ioy) = out
	if ioy lt 0 then imbig(*,2+ly+ioy:1+ly) = out
ENDIF


;Parte frazionaria dello shift
IF fracx NE 0 THEN BEGIN
	matx = [fracx, 1.-fracx, 0] * (fracx GT 0) + [0, 1.+fracx, -fracx] * (fracx LT 0)
	imbig = convol(imbig, matx)
ENDIF
IF fracy NE 0 THEN BEGIN
	maty = [fracy, 1.-fracy, 0] * (fracy GT 0) + [0, 1.+fracy, -fracy] * (fracy LT 0)
	maty = transpose(maty)
	imbig = convol(imbig, maty)
ENDIF

imx_s = imbig[1:s[1], 1:s[2]]

RETURN, imx_s

END