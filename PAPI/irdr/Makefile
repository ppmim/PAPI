#
# irdr global Makefile (for gcc)
#

CC        = gcc
OUTPUT    = -o $@
IFLAGS    = -I./src/include -I./extern/wcssubs
CFLAGS    = -O3 -Wall -ansi -pedantic
# CFLAGS    = -g 
RM        = rm -f
RANLIB    = true
AR        = ar
ARFLAGS   = rcv
LDFLAGS   = -L./lib -L./extern/wcssubs
LDLIBS    = -lirdr -lc -lm -lwcs
COMPFLAGS = $(CFLAGS) $(IFLAGS)
LINKFLAGS = $(LDFLAGS) $(LDLIBS)
COMPILE.c = $(CC) $(CFLAGS) $(IFLAGS) -c
LINK.c    = $(CC) $(CFLAGS) $(IFLAGS) $(LDFLAGS)

.c:
	$(LINK.c) $(OUTPUT) $< $(LDLIBS)
.c.o:
	$(COMPILE.c) $(OUTPUT) $<
.c.a:
	$(COMPILE.c) -o $% $<
	$(AR) $(ARFLAGS) $@ $

default: irdr

all: clean wcs irdr

irdr: cleanirdr irdrlib ./bin/statsmask ./bin/dithers \
	./bin/cubemean ./bin/gainmap ./bin/mask ./bin/floatimage \
	./bin/normalize ./bin/initwcs ./bin/unshiftimage \
	./bin/shiftimage ./bin/dithercubemean ./bin/darkflat \
	./bin/offsets ./bin/offsetsnowcs ./bin/flat ./bin/imcalc \
	./bin/mosaic \
	./bin/destripe ./bin/dilate ./bin/sky ./bin/parab \
	./bin/printwcs ./bin/avgkey ./bin/stats ./bin/avgoffsets \
	./bin/skyfilter ./bin/deshadow ./bin/imagegrid ./bin/defringe \
	./bin/skyfilteronoff ./bin/skyfilteroffon ./bin/skyfilter_general \
	./bin/linflat ./bin/lindarkflat ./bin/dark ./bin/skyfilter_single ./bin/my_mode ./bin/nlinmap

OBJ_FITSIO = src/fitsIO/fitsIO.o

OBJ_MATH = src/math/median.o \
		src/math/mean.o \
		src/math/dilate.o \
		src/math/bisearch.o \
		src/math/kselect.o \
		src/math/correlate.o \
		src/math/parabola.o \
		src/math/shift.o \
		src/math/hist.o \
		src/math/nrlfit.o

OBJ_SKY = src/sky/skysub.o \
		src/sky/getwmap.o \
		src/sky/stripe.o

OBJ_CIRSI = 

OBJ_MISC = src/misc/eprintf.o \
		src/misc/convert.o \
		src/misc/listIO.o

OBJ_CUBE = src/cube/cube.o

OBJ_IRDR = $(OBJ_FITSIO) $(OBJ_MATH) $(OBJ_CUBE) $(OBJ_MISC) \
		$(OBJ_OFFSETS) $(OBJ_REG) $(OBJ_CIRSI) $(OBJ_SKY)

clean: cleanirdr cleanwcs

cleanirdr:
	$(RM) src/*/*.o bin/* lib/*.a

cleanwcs:
	$(RM) extern/wcssubs/*.o
	$(RM) extern/wcssubs/*.a
	$(RM) extern/wcstools/*.o
	$(RM) extern/wcstools/*.a
	$(RM) extern/wcstools/libwcs/*.o
	$(RM) extern/wcstools/libwcs/*.a
	$(RM) extern/wcstools/bin/*

wcs: 
	@(cd extern/wcssubs ; make)
	@(cd extern/wcstools; make all)

match: 
	@(cd extern/match ; rm -f *.o match ; make)

#
# Build irdr library
#

irdrlib: ./lib/libirdr.a

./lib/libirdr.a: $(OBJ_IRDR)
	$(AR) $(ARFLAGS) ./lib/libirdr.a $(OBJ_IRDR)
	$(RANLIB) ./lib/libirdr.a

#
# Build irdr executables
#

./bin/stats: src/drivers/stats.c
	$(CC) $(COMPFLAGS) src/drivers/stats.c \
		-o bin/stats $(LINKFLAGS)

./bin/imcalc: src/drivers/imcalc.c
	$(CC) $(COMPFLAGS) src/drivers/imcalc.c \
		-o bin/imcalc $(LINKFLAGS)

./bin/floatimage: src/drivers/floatimage.c
	$(CC) $(COMPFLAGS) src/drivers/floatimage.c \
		-o bin/floatimage $(LINKFLAGS)

./bin/statsmask: src/drivers/statsmask.c
	$(CC) $(COMPFLAGS) src/drivers/statsmask.c \
		-o bin/statsmask $(LINKFLAGS)

./bin/normalize: src/drivers/normalize.c
	$(CC) $(COMPFLAGS) src/drivers/normalize.c \
		-o bin/normalize $(LINKFLAGS)

./bin/cubemean: src/drivers/cubemean.c
	$(CC) $(COMPFLAGS) src/drivers/cubemean.c \
		-o bin/cubemean $(LINKFLAGS)

./bin/dithercubemean: src/drivers/dithercubemean.c
	$(CC) $(COMPFLAGS) src/drivers/dithercubemean.c \
		-o bin/dithercubemean $(LINKFLAGS)

./bin/mask: src/drivers/mask.c
	$(CC) $(COMPFLAGS) src/drivers/mask.c \
		-o bin/mask $(LINKFLAGS)

./bin/destripe: src/drivers/destripe.c
	$(CC) $(COMPFLAGS) src/drivers/destripe.c \
		-o bin/destripe $(LINKFLAGS)

./bin/flat: src/drivers/flat.c
	$(CC) $(COMPFLAGS) src/drivers/flat.c \
		-o bin/flat $(LINKFLAGS)

./bin/darkflat: src/drivers/darkflat.c
	$(CC) $(COMPFLAGS) src/drivers/darkflat.c \
		-o bin/darkflat $(LINKFLAGS)

./bin/gainmap: src/drivers/gainmap.c
	$(CC) $(COMPFLAGS) src/drivers/gainmap.c \
		-o bin/gainmap $(LINKFLAGS)

./bin/shiftimage: src/drivers/shiftimage.c
	$(CC) $(COMPFLAGS) src/drivers/shiftimage.c \
		-o bin/shiftimage $(LINKFLAGS)

./bin/unshiftimage: src/drivers/unshiftimage.c
	$(CC) $(COMPFLAGS) src/drivers/unshiftimage.c \
		-o bin/unshiftimage $(LINKFLAGS)

./bin/dilate: src/drivers/dilate.c
	$(CC) $(COMPFLAGS) src/drivers/dilate.c \
		-o bin/dilate $(LINKFLAGS)

./bin/offsets: src/drivers/offsets.c
	$(CC) $(COMPFLAGS) src/drivers/offsets.c \
		-o bin/offsets $(LINKFLAGS)

./bin/mosaic: src/drivers/mosaic.c
	$(CC) $(COMPFLAGS) src/drivers/mosaic.c \
		-o bin/mosaic $(LINKFLAGS)

./bin/offsetsnowcs: src/drivers/offsetsnowcs.c
	$(CC) $(COMPFLAGS) src/drivers/offsetsnowcs.c \
		-o bin/offsetsnowcs $(LINKFLAGS)

./bin/imagegrid: src/drivers/imagegrid.c
	$(CC) $(COMPFLAGS) src/drivers/imagegrid.c \
		-o bin/imagegrid $(LINKFLAGS)

./bin/initwcs: src/cirsi/initwcs.c
	$(CC) $(COMPFLAGS) src/cirsi/initwcs.c \
		-o bin/initwcs $(LINKFLAGS)

./bin/printwcs: src/cirsi/printwcs.c
	$(CC) $(COMPFLAGS) src/cirsi/printwcs.c \
		-o bin/printwcs $(LINKFLAGS)

./bin/avgkey: src/drivers/avgkey.c
	$(CC) $(COMPFLAGS) src/drivers/avgkey.c \
		-o bin/avgkey $(LINKFLAGS)

./bin/deshadow: src/cirsi/deshadow.c
	$(CC) $(COMPFLAGS) src/cirsi/deshadow.c \
		-o bin/deshadow $(LINKFLAGS)

./bin/avgoffsets: src/cirsi/avgoffsets.c
	$(CC) $(COMPFLAGS) src/cirsi/avgoffsets.c \
		-o bin/avgoffsets $(LINKFLAGS)

./bin/dithers: src/cirsi/dithers.c
	$(CC) $(COMPFLAGS) src/cirsi/dithers.c \
		-o bin/dithers $(LINKFLAGS)

./bin/parab: src/drivers/parab.c
	$(CC) $(COMPFLAGS) src/drivers/parab.c \
		-o bin/parab $(LINKFLAGS)

./bin/skyfilter: src/drivers/skyfilter.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilter.c \
		-o bin/skyfilter $(LINKFLAGS)

./bin/skyfilteronoff: src/drivers/skyfilteronoff.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilteronoff.c \
		-o bin/skyfilteronoff $(LINKFLAGS)

./bin/skyfilteroffon: src/drivers/skyfilteroffon.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilteroffon.c \
		-o bin/skyfilteroffon $(LINKFLAGS)

./bin/sky: src/drivers/sky.c
	$(CC) $(COMPFLAGS) src/drivers/sky.c \
		-o bin/sky $(LINKFLAGS)

./bin/defringe: src/cirsi/defringe.c
	$(CC) $(COMPFLAGS) src/cirsi/defringe.c \
		-o bin/defringe $(LINKFLAGS)

./bin/linflat: src/drivers/linflat.c
	$(CC) $(COMPFLAGS) src/drivers/linflat.c \
		-o bin/linflat $(LINKFLAGS)

./bin/lindarkflat: src/drivers/lindarkflat.c
	$(CC) $(COMPFLAGS) src/drivers/lindarkflat.c \
		-o bin/lindarkflat $(LINKFLAGS)

./bin/dark: src/drivers/dark.c
	$(CC) $(COMPFLAGS) src/drivers/dark.c \
		-o bin/dark $(LINKFLAGS)

./bin/skyfilter_single: src/drivers/skyfilter_single.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilter_single.c \
		-o bin/skyfilter_single $(LINKFLAGS)

./bin/skyfilter_general: src/drivers/skyfilter_general.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilter_general.c \
		-o bin/skyfilter_general $(LINKFLAGS)
		
./bin/my_mode: src/drivers/my_mode.c
	$(CC) $(COMPFLAGS) src/drivers/my_mode.c \
		-o bin/my_mode $(LINKFLAGS)

./bin/my_mode: src/drivers/nlinmap.c
	$(CC) $(COMPFLAGS) src/drivers/nlinmap.c \
		-o bin/nlinmap $(LINKFLAGS)		
