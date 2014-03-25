DEFINE/PAR P1 ? ? "catalog.dat"         !Catalog of focusfield, first column RA in deg,second DEC in deg
DEFINE/PAR P2 ? I/A "Image.bdf"                 !Image taken of desired focus field, pipelined 
DEFINE/PAR P3 " " C "Telescope position RA hh:mm:ss.s"
DEFINE/PAR P4 " " C "Telescope position DEC dd:mmm:ss.s"
DEFINE/PAR P5 test C "Name of focus field"

DEFINE/LOCAL size/R/1/1 7
DEFINE/LOCAL bin/R/1/1 1
DEFINE/LOCAL dummy_RA/C/1/16 {P3}
DEFINE/LOCAL dummy_DEC/C/1/16 {P4}
define/local alpha0/R/1/1 0
define/local delta0/R/1/1 0

alpha0 = ({dummy_RA(1:2)} + {dummy_RA(4:5)}/60 + {dummy_RA(7:10)}/3600 )*15
delta0 = {dummy_DEC(1:2)} + {dummy_DEC(4:5)}/60 + {dummy_DEC(7:10)}/3600 

@@ O2K_UTIL:/obs_macros/new/ad2xy {P1} {alpha0} {delta0}



alpha0 = 1024
delta0 = 1024

COM/tab cat  :x_pix = (:x_rad *206265/0.44942  + {alpha0} )
COM/tab cat  :y_pix = (:y_rad *206265/0.44942  + {delta0} )
COM/TAB cat :XSTART = :x_pix - {size}
COM/TAB cat :XEND = :x_pix + {size}
COM/TAB cat :YSTART = :y_pix - {size}
COM/TAB cat :YEND = :y_pix + {size}

CENTER/GAUSS {P2},cat cat

SELECT/TAB cat :XCEN.NE.NULL

CLEAR/CHANNEL OVERLAY

load/ima {P2} sc=-3,a ce=c

load/tab cat :XCEN :YCEN ? ? ? red
load/tab cat :x_pix :y_pix ? ? 5 blue

CREATE/TAB USER_{P5} 2 500
COPY/TT cat :XCEN USER_{P5} :X 
COPY/TT cat :YCEN USER_{P5} :Y 
!WRITE/DESC {P5}.tbl RA/C/1/16 {dummy_RA}
!WRITE/DESC {P5}.tbl DEC/C/1/16 {dummy_DEC}

WRITE/OUT "If you are not satisfied with the identification, please run prg again with different alpha, delta!" 

return
