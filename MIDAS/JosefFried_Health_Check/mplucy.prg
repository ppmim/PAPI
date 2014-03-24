! MPLUCY.PRG
!
! Point-sources-plus-smooth-background restoration.
!
! Richard Hook, ST-ECF, April 1993.
! Revised for V0.4, October 1993
! Hidden parameter handling added - November 1993
!
! Procedure for control within MIDAS
!
define/par p1 ? ima "Input data image: "
define/local data/c/1/80 {p1}
!
define/par p2 ? ima "PSF image(s): "
define/local psfs/c/1/80 {p2}
!
define/par p3 ? n "Fractional strength of entropy term: "
define/local fracent/D/1/1 {p3}
!
define/par p4 ? n "Number of Iterations: "
define/local niter/i/1/1 {p4}
!
define/par p5 ? ima "Deconvolved output image: "
define/local decon/c/1/80 {p5}
!
define/par p6 ? ima "File containing list of star positions: "
define/local starlist/c/1/80 {p6}
!
define/par p7 ? ima "Output photometric text table: "
define/local outtab/c/1/80 {p7}
!
define/par p8 ? cha "Check hidden parameter values (y/n)? "
define/local checkpar/c/1/1 {p8}
!
! We have got all the normal parameters on the command line,
! we now check and display the "hidden" parameters
!
!
wk stat/i/1/1 0
!
if checkpar .eq. "y" then
   write/out "-Verifying hidden parameter values:"
   write/key stat/i/1/1 0
   compute/key stat = m$existk("verbose")
   if stat .ne. 1 write/key verbose/c/1/1 y
   inquire/key verbose "Give verbose diagnostics ({verbose})? "
   !
   compute/key stat = m$existk("accel")
   if stat .ne. 1 write/key accel/c/1/1 y
   inquire/key accel "Use accelerated method ({accel})? "
   !
   compute/key stat = m$existk("fprior")
   if stat .ne. 1 write/key fprior/c/1/1 y
   inquire/key fprior "Use floating prior ({fprior})? "
   !
   compute/key stat = m$existk("Xsubsam")
   if stat .ne. 1 write/key Xsubsam/i/1/1 1
   inquire/key Xsubsam "X subsampling ({Xsubsam})? "
   !
   compute/key stat = m$existk("Ysubsam")
   if stat .ne. 1 write/key Ysubsam/i/1/1 1
   inquire/key Ysubsam "Y subsampling ({Ysubsam})? "
   !
   compute/key stat = m$existk("modpsf")
   if stat .ne. 1 write/key modpsf/c/1/1 n
   inquire/key modpsf "Modify the PSF iteratively ({modpsf})? "
   !
   write/key s/i/1/1 1
   write/key t/i/1/1 1
   !
   if modpsf .eq. "y" then
      compute/key stat = m$existk("outpsf")
      if stat .ne. 1 write/key outpsf/c/1/80 newpsf
      inquire/key outpsf "Name for output modified PSF image ({outpsf})? "
      compute/key s = {aux_mode(7)}+1
      compute/key t = 80-{s}+1
      if s .gt. 1 write/key outpsf/c/{s}/{t} " "
   endif
   !
   compute/key stat = m$existk("back")
   if stat .ne. 1 write/key back/c/1/80 back  
   inquire/key back "Name for output background image ({back})? "
   compute/key s = {aux_mode(7)}+1
   compute/key t = 80-{s}+1
   if s .gt. 1 write/key back/c/{s}/{t} " "
   !
   compute/key stat = m$existk("magzero")
   if stat .ne. 1 write/key magzero/d/1/1 20.0
   inquire/key magzero "Zero point for magnitudes ({magzero})? "
   !
   compute/key stat = m$existk("skernel")
   if stat .ne. 1 write/key skernel/d/1/1 1.0
   inquire/key skernel "Smoothing kernel sigma ({skernel})? "

else
   write/out "-Hidden parameters are as follows:"
   !
   compute/key stat = m$existk("verbose")
   if stat .eq. 1 then
      write/out "verbose flag is set to: " {verbose}
   else
      write/key verbose/c/1/1 Y
      write/out "setting verbose flag to the default: y"
   endif
!
compute/key stat = m$existk("accel")
if stat .eq. 1 then
   write/out "accel flag is set to: " {accel}
else
   write/key accel/c/1/1 Y
   write/out "setting accel flag to the default: y"
endif
!
compute/key stat = m$existk("fprior")
if stat .eq. 1 then
   write/out "fprior flag is set to: " {fprior}
else
   write/key fprior/c/1/1 Y
   write/out "setting fprior flag to the default: y"
endif
!
compute/key stat = m$existk("Xsubsam")
if stat .eq. 1 then
   write/out "Xsubsam is set to: " {Xsubsam}
else
   write/key Xsubsam/i/1/1 1
   write/out "setting Xsubsam to the default: 1"
endif
!
compute/key stat = m$existk("Ysubsam")
if stat .eq. 1 then
   write/out "Ysubsam is set to: " {Ysubsam}
else
   write/key Ysubsam/i/1/1 1
   write/out "setting Ysubsam to the default: 1"
endif
!
compute/key stat = m$existk("modpsf")
if stat .eq. 1 then
   write/out "modpsf flag is set to: " {modpsf}
else
   write/key modpsf/c/1/1 N
   write/out "setting modpsf flag to the default: n"
endif
if modpsf .eq. "y" then
   compute/key stat = m$existk("outpsf")
   if stat .eq. 1 then
      write/out "outpsf image name is set to: " {outpsf}
   else
      write/key outpsf/c/1/80 "newpsf"
      write/out "setting outpsf image name to the default: outpsf"
   endif
endif
!
compute/key stat = m$existk("back")
if stat .eq. 1 then
   write/out "back output image name is set to: " {back}
else
   write/key back/c/1/80 " " all
   write/out "setting output image name to the default: <null>"
endif
!
compute/key stat = m$existk("magzero")
if stat .eq. 1 then
   write/out "magzero is set to: " {magzero}
else
   write/key magzero/d/1/1 20.0
   write/out "setting magzero to the default of 20.0"
endif
!
compute/key stat = m$existk("skernel")
if stat .eq. 1 then
   write/out "skernel is set to: " {skernel}
else
   write/key skernel/d/1/1 1.0
   write/out "setting skernel to the default of 1.0"
endif
endif
!
write/out "-Starting execution"
!
run FMP:mplucy.exe
!
copy/dd {p1} lhcuts {p5} lhcuts
copy/dd {p1} scale {p5} scale 
wd {p5} mplucy/c/1/80 " mplucy {p1} {p2} {p3} {p4} {p5} {p6} {p7}"
copy/kd skernel {p5} skernel/d/1/1
copy/dd {p1} ident {p5} ident
load {p5}
