
;
; read SExtractor catalog file and FITS header keywords to return array of:
; fwhm, elongation, position angle, nobjs, datamode, datasig, airmass
;

function readinfo, fn, x, y, mag
    nparm = 7
    catfn = fn + '.cat'

    if ((findfile(fn))(0) eq '' or (findfile(catfn))(0) eq '') then $
        return, replicate(0.0, nparm)

;    readcol, catfn, x, y, mag, fwhm, elong, theta, format='D,D,D,D,D,D', /sil
    readcat, catfn, x, y, mag, fwhm, elong, theta

    hdr = headfits(fn)

    return, [sxpar(hdr, 'DATAMODE'), sxpar(hdr, 'DATASIG'), $
             sxpar(hdr, 'AIRMASS'), median(fwhm), median(elong), $
             median(theta), n_elements(x)]
end

;
; read a SExtractor catalog

pro readcat, fn, x, y, mag, fwhm, elong, theta

nmax = 3000
x = fltarr(nmax)
y = fltarr(nmax)
mag = fltarr(nmax)
fwhm = fltarr(nmax)
elong = fltarr(nmax)
theta = fltarr(nmax)

xtmp = 0.0
ytmp = 0.0
magtmp = 0.0
fwhmtmp = 0.0
elongtmp = 0.0
thetatmp = 0.0

openr, input, fn, /get_lun
n = 0
while (not eof(input)) do begin
    readf,input,xtmp,ytmp,magtmp,fwhmtmp,elongtmp,thetatmp
    x(n) = xtmp
    y(n) = ytmp
    mag(n) = magtmp
    fwhm(n) = fwhmtmp
    elong(n) = elongtmp
    theta(n) = thetatmp
    n = n + 1
endwhile

if (n lt 1) then n = 1

x = x(0:n-1)
y = y(0:n-1)
mag = mag(0:n-1)
fwhm = fwhm(0:n-1)
elong = elong(0:n-1)
theta = theta(0:n-1)

free_lun, input
return
end

;
; match x,y,mag lists and output median and stdev mag difference
;

function match,x1,y1,m1,x2,y2,m2
    n = 0
    minsep = 10.0                     ; min separation in pixels
    dmag = fltarr(n_elements(x1))

    for i = 0L, n_elements(x1) - 1 do begin
        dists = sqrt((x1(i)-x2)^2. + (y1(i)-y2)^2.)
        nearest = where(dists eq min(dists))
        j = nearest(0)

        if (dists(j) lt minsep) then begin
            dmag(n) = m1(i) - m2(j)
            n = n + 1
        endif
    endfor

    if (n lt 3) then return, [0.0, 0.0]

    dmag = dmag(0:n-1)

    return, [median(dmag), robust_sigma(dmag)]
end

;
; read SExtractor catalogs of each loop, combined loops, and coadded dither
; sets to make plots of FWHM, elongation, datamode, and datasig verse run no.
;

pro verify, chipno, parmno, legend=legend

color, 'white', 0   &  white  = 0 
color, 'black', 1   &  black  = 1
color, 'red', 2     &  red    = 2
color, 'green', 3   &  green  = 3
color, 'blue', 4    &  blue   = 4
color, 'orange', 5  &  orange = 5
color, 'yellow', 6  &  yellow = 6

readcol, 'dithers', runbegs, runends, format='I,I', /silent

nloop   = 3
maxnrun = 999
ndither = n_elements(runbegs)

parnames = ['Data MODE', 'Data STDEV', 'Airmass', 'FWHM [pix]', $
            'Elongation (a/b)', 'Position Angle', 'Nobjs', 'Mag Offsets', $
            'Mag Stdev', 'Dither X STDEV', 'Dither Y STDEV']

npar     = n_elements(parnames)
parloops = fltarr(npar, maxnrun, nloop)          ; single loops
parloopc = fltarr(npar, maxnrun)                 ; combined loops
parcdset = fltarr(npar, ndither)                 ; coadded dither set
xtnames  = replicate(' ', ndither)
xtvals   = intarr(ndither)

irun = 0
ninfo = 7                                  ; output array size of readinfo()
chipstr = string(chipno, format='(I1)')
fstr = '(A,X,9(F9.2,X))'

openw, output, 'verify.c'+chipstr, /get_lun
printf,output,'Filename',parnames,format='(10(A,X))'

for iset = 0, n_elements(runbegs) - 1 do begin
    runbeg = runbegs(iset)
    runend = runends(iset)
    rbstr = strcompress(string(runbeg),/remove)
    restr = strcompress(string(runend), /remove)
    fn = 'irx.r' + rbstr + 'to' + restr + '.c' + chipstr + '.fits'
    xtnames(iset) = rbstr
    xtvals(iset) = irun

    parcdset(0:ninfo-1, iset) = readinfo(fn,xd,yd,magd)
    printf,output,fn+'.cat',parcdset(0:ninfo-1,iset),format=fstr

    offsetsfn = 'offsets.r' + rbstr + 'to' + restr
    readcol,offsetsfn,offrun,offx,offy,matchfrac,offxsig,offysig,/silent
    border = fix(max([abs(offx),abs(offy)]) + 3.5)       ; see shift.c
    border = border + (4 - (border mod 4))

    for run = runbeg, runend do begin
        runstr = string(run, format='(I5.5)')
        fn = 'irx_' + runstr + '_c' + chipstr + '.fits.skysub'
        parloopc(0:ninfo-1, irun) = readinfo(fn,x,y,mag)
        printf,output,fn+'.cat',parloopc(0:ninfo-1,irun),format=fstr

        use = where(run eq offrun, count)

        if (count eq 1) then begin
            use = use(0)
            x = x + border - offx(use)
            y = y + border - offy(use)
            parloopc(ninfo:ninfo+1, irun) = match(x,y,mag,xd,yd,magd)
            parloopc(ninfo+2:ninfo+3, irun) = [offxsig(use),offysig(use)]
        endif else begin
            print, 'ERR: run not in offsets file: ', run
        endelse

        for loop = 1, nloop do begin
            loopstr = string(loop, format='(I3.3)')
            fn = 'irx_'+runstr+'_c' + chipstr + '_' + loopstr + '.fits.skysub'
            parloops(0:ninfo-1, irun, loop-1) = readinfo(fn,x,y,mag)
            printf,output,fn+'.cat',parloops(0:ninfo-1,irun,loop-1),format=fstr
        endfor

        irun = irun + 1
    endfor
endfor

parloops = parloops(*,0:irun-1,*)
parloopc = parloopc(*,0:irun-1)


if (parmno lt 5) then begin

  maxv1 = max(parloops(parmno,*,*))
  maxv2 = max(parloopc(parmno,*))
  maxv3 = max(parcdset(parmno,*))
  maxv  = max([maxv1,maxv2,maxv3])

  plot,parloops(parmno,*,0),xtickname=xtnames,xtickv=xtvals,xticks=1>(iset-1), $
    psym=2,title='Chip'+chipstr,xtitle='Run Number',ytitle=parnames(parmno), $
    yrange=[0,maxv]
  oplot,parloops(parmno,*,0),psym=2,color=red
  oplot,parloops(parmno,*,1),psym=1,color=orange
  oplot,parloops(parmno,*,2),psym=4,color=yellow
  oplot,parloopc(parmno,*),psym=5,color=green
  oplot,xtvals,parcdset(parmno,*),psym=6,color=blue
endif

if (parmno eq 5) then begin
  plot,parloops(parmno,*,0),xtickname=xtnames,xtickv=xtvals,xticks=1>(iset-1), $
    psym=2,title='Chip'+chipstr,xtitle='Run Number',ytitle=parnames(parmno), $
    yrange=[-90,90],ystyle=1
  oplot,parloops(parmno,*,0),psym=2,color=red
  oplot,parloops(parmno,*,1),psym=1,color=orange
  oplot,parloops(parmno,*,2),psym=4,color=yellow
  oplot,parloopc(parmno,*),psym=5,color=green
  oplot,xtvals,parcdset(parmno,*),psym=6,color=blue
endif

if (parmno eq 6) then begin
  plot,parloops(parmno,*,0),xtickname=xtnames,xtickv=xtvals,xticks=1>(iset-1), $
    psym=2,title='Chip'+chipstr,xtitle='Run Number',ytitle=parnames(parmno), $
    yrange=[0,max(parcdset(parmno,*)+5)],ystyle=1
  oplot,parloops(parmno,*,0),psym=2,color=red
  oplot,parloops(parmno,*,1),psym=1,color=orange
  oplot,parloops(parmno,*,2),psym=4,color=yellow
  oplot,parloopc(parmno,*),psym=5,color=green
  oplot,xtvals,parcdset(parmno,*),psym=6,color=blue
endif

if (parmno eq 7 or parmno eq 8) then begin
  ploterr,parloopc(7,*),parloopc(8,*),xtickname=xtnames, $
    xtickv=xtvals,xticks=1>(iset-1),psym=5,title='Chip'+chipstr, $
    xtitle='Run Number',ytitle=parnames(parmno),color=black
  oploterr,parloopc(7,*),parloopc(8,*),color=green,psym=5
endif

if (parmno eq 9 or parmno eq 10) then begin
  plot,parloopc(parmno,*,0),xtickname=xtnames,xtickv=xtvals,xticks=1>(iset-1), $
    psym=5,title='Chip'+chipstr,xtitle='Run Number',ytitle=parnames(parmno)
  oplot,parloopc(parmno,*),psym=5,color=green
endif

if (keyword_set(legend)) then begin
  if (parmno eq 9 or parmno eq 10) then begin
    legend, ['loopcombined'],psym=[5], colors=[green], $
          textcolors=[black],/top,/right
  endif else if (parmno ne 7) then begin
    legend, ['loop 1', 'loop 2', 'loop 3', 'loopcombined', 'ditherset'], $
          psym=[2,1,4,5,6], colors=[red,orange,yellow,green,blue], $
          textcolors=[black,black,black,black,black],/bottom,/right
  endif
endif

free_lun, output

end
