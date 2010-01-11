pro doverify
;
; run verify.pro on all four chips
;

parnames = ['bkg', 'sig', 'airmass', 'fwhm', 'elong', 'theta', 'nobjs', $
            'magoffset', 'magstdev', 'ditherxsig', 'ditherysig']

cd, current=cwd

dirs = strsplit(cwd, '/', /extract)
ndirs = n_elements(dirs)
basenm = dirs(ndirs-2) + '_' + dirs(ndirs-1)

!p.multi=[0,2,2]

openw,1,basenm + '.html'
printf,1,'<body bgcolor="#FFFFF">'

for parmno = 0, n_elements(parnames)-1 do begin

  if (parmno ne 8) then begin
    verify,1,parmno,/legend
    verify,2,parmno
    verify,3,parmno
    verify,4,parmno

    xyouts,0.5,0.5,cwd,align=0.5,/normal

    tvlct,tempr,tempg,tempb,/get

    gifname = 'gifs/' + basenm + '.' + parnames(parmno) + '.gif'

    write_gif, gifname, tvrd(), tempr, tempg, tempb

    printf,1,'<img src="' + gifname + '">'
  endif

endfor


close,1

end
