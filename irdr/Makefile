#
# irdr global Makefile (for gcc)
#

CC        = gcc
OUTPUT    = -o $@
IFLAGS    = -I./src/include  -I./extern/wcstools -I./extern/wcstools/libwcs
CFLAGS    = -O3 -Wall -ansi -pedantic
# CFLAGS    = -g 
RM        = rm -f
MKDIR	  = mkdir -p
RANLIB    = true
AR        = ar
ARFLAGS   = rcv
LDFLAGS   = -L./lib   -L./extern/wcstools/libwcs
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


# It has been modified to only have the binaries used by PAPI
irdr: cleanirdr irdrlib  \
	./bin/dithercubemean ./bin/offsets ./bin/cubemean \
	./bin/skyfilter ./bin/skyfilteronoff ./bin/skyfilteroffon ./bin/skyfilter_general \
	./bin/skyfilter_single \
	./bin/my_mode

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
	$(MKDIR) lib bin

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

./bin/cubemean: src/drivers/cubemean.c
	$(CC) $(COMPFLAGS) src/drivers/cubemean.c \
		-o bin/cubemean $(LINKFLAGS)

./bin/dithercubemean: src/drivers/dithercubemean.c
	$(CC) $(COMPFLAGS) src/drivers/dithercubemean.c \
		-o bin/dithercubemean $(LINKFLAGS)

./bin/offsets: src/drivers/offsets.c
	$(CC) $(COMPFLAGS) src/drivers/offsets.c \
		-o bin/offsets $(LINKFLAGS)

./bin/skyfilter: src/drivers/skyfilter.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilter.c \
		-o bin/skyfilter $(LINKFLAGS)

./bin/skyfilteronoff: src/drivers/skyfilteronoff.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilteronoff.c \
		-o bin/skyfilteronoff $(LINKFLAGS)

./bin/skyfilteroffon: src/drivers/skyfilteroffon.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilteroffon.c \
		-o bin/skyfilteroffon $(LINKFLAGS)

./bin/skyfilter_single: src/drivers/skyfilter_single.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilter_single.c \
		-o bin/skyfilter_single $(LINKFLAGS)

./bin/skyfilter_general: src/drivers/skyfilter_general.c
	$(CC) $(COMPFLAGS) src/drivers/skyfilter_general.c \
		-o bin/skyfilter_general $(LINKFLAGS)
		
./bin/my_mode: src/drivers/my_mode.c
	$(CC) $(COMPFLAGS) src/drivers/my_mode.c \
		-o bin/my_mode $(LINKFLAGS)

