DEFI/PAR P1 z C " table_in=?"
DEFI/PAR P2 z2 C "table_out=?"
!DEFI/PAR P3 575.,489.,700,20,.41,2.21,2.515E4 N "xc,yc,rad_max,nbin,pix/",pix/Mpc,gal/sqdeg=?"
!DEFI/PAR P4 2,3,11,1000,1000 N " x_col,y_col,weicol,npx,npy=?"
DEFI/PAR P3 ? N "xc,yc,rad_max,nbin,pix/",pix/Mpc,gal/sqdeg=?"
DEFI/PAR P4 ? N " x_col,y_col,npx,npy=?"



WRITE/KEYW INPUTR/R/1/20 'P3'
WRITE/KEYW INPUTI/I/1/6 'P4'
RUN FMP:ring
