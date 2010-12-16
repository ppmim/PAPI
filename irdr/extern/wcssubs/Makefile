CFLAGS= -g

OBJS =  wcsinit.o wcs.o wcscon.o slasubs.o fitsfile.o imhfile.o \
	hget.o hput.o iget.o imio.o worldpos.o tnxpos.o dsspos.o platepos.o \
        wcslib.o lin.o cel.o proj.o sph.o wcstrig.o dateutil.o

libwcs.a:	$(OBJS)
	ar rv $@ $?
	ranlib $@

cel.o:		wcslib.h
fitsfile.o:	fitsfile.h
hget.o:		fitshead.h
hput.o:		fitshead.h
iget.o:		fitshead.h
imhfile.o:	fitsfile.h
imio.o:		fitshead.h
lin.o:		wcslib.h
platepos.o:	wcs.h
proj.o:		wcslib.h
sph.o:		wcslib.h
tnxpos.o:	wcs.h
wcs.o:		fitshead.h wcs.h
wcsinit.o:	fitshead.h wcs.h
wcscon.o:	wcs.h
wcslib.o:	wcslib.h
wcstrig.o:	wcslib.h
worldpos.o:	wcs.h
dateutil.o:	fitsfile.h
cel.h:		wcslib.h
fitsfile.h:	fitshead.h imio.h
proj.h:		wcslib.h
wcs.h:		wcslib.h fitshead.h
