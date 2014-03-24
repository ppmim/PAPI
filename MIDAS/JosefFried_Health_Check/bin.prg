!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! bin in out mode npix_x,npix_y parameter 
!
!	binning of frames
!
!	in,out        = in,out frames
!	npix_x,npix_y = size of result frame
!
!	mode = bicu   bicubic spline interpolation
!	       xy     monotonic interploation first x then y
!              yx     monotonic interploation first y then x
!              mean   monotonic interploation 0.5*(xy+yx)
!              integ  integer shift and rebinning
!              linear double linear in x and y with rotation and shear
!
!       parameter
!	       bicu  :
!	       xy,yx,mean :
!              integ  xstart,ystart,ibinx,ibiny
!              linear xstart,ystart,dx,dy,phi(rad)
!              
!-----------------------------------------------
!
defi/par p1 z ima "inputframe=?"
defi/par p2 z1 ima "outputframe=?"
defi/par p3 mean c "mode=?"
defi/par p4 1000,1000 n "npo_x,npo_y=?"
defi/par p5 1,1,1,1,0 n "params=?"
WRITE/KEYW IN_A/C/1/60 'P1'
WRITE/KEYW OUT_A/C/1/60 'P2'
WRITE/KEYW action/C/1/8 'P3'
WRITE/KEYW INPUTI/I/1/5 'P4'
WRITE/KEYW INPUTR/R/1/5 'P5'
RUN FMP:bin
cop/dd {p1} ident {p2} ident
