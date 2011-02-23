#	Makefile for OMEGA_pipeline	
#	Rene Fassbender	23.08.02
#	
#	All secundary modules are linked together with the main module
#
#	Vers: 2.0	22.04.04	changed for use outside MPIA
#______________________________________________________________________________________
#

MID_INCLUDE=/disk-a/caha/panic/DEVELOP/MIDAS/08FEBpl1.1/incl
MID_LIB=/disk-a/caha/panic/DEVELOP/MIDAS/08FEBpl1.1/lib
CC = gcc -I$(MID_INCLUDE)
CFLAGS = -Wall -ggdb -L $MID_LIB

BIBS =-L$(MID_LIB) -lmidas -lm  
#-lsocket -lnsl
#/photo/source/NUM_RECIP/NUM_RECIP.a /photo/source/SUB/phot_sub.a \
#/photo/source/LIB/phot_lib.a
#-ldisp -lgen -lgmidas -lidi -lsubmid /mpisys/lib/nag/libnag.a

O2K_UTIL=/disk-a/caha/panic/DEVELOP/PIPELINE/Ref/O2k/o2k_utilities/
# paths of the directories
SRCDIR = ${O2K_UTIL}/pipeline/C_code/
EXEDIR = ${O2K_UTIL}/pipeline/exe/omega_pipeline.exe
OBJDIR = ${O2K_UTIL}/pipeline/obj/



# all modules must be listed here
OBJLIST =  OMEGA_pipeline.o flat_correction.o subframe_statistics.o normalize_frame.o extract_name.o load_stack.o\
 get_sky.o save_sky.o subtract_sky.o copy_frame.o original_level.o minimum_sky.o sky_clipping.o frame_info.o \
sum_preparation.o make_cleansum.o dirty_sum.o put_tosum.o save_dirty.o find_object.o cut_frame.o \
badpixel_correction.o elaborate_sum.o


OBJ_PATH_LIST =  $(OBJDIR)OMEGA_pipeline.o $(OBJDIR)flat_correction.o $(OBJDIR)subframe_statistics.o \
$(OBJDIR)normalize_frame.o $(OBJDIR)extract_name.o $(OBJDIR)load_stack.o $(OBJDIR)get_sky.o $(OBJDIR)save_sky.o \
$(OBJDIR)subtract_sky.o $(OBJDIR)copy_frame.o $(OBJDIR)original_level.o $(OBJDIR)minimum_sky.o \
$(OBJDIR)sky_clipping.o $(OBJDIR)frame_info.o $(OBJDIR)sum_preparation.o $(OBJDIR)make_cleansum.o \
$(OBJDIR)dirty_sum.o $(OBJDIR)put_tosum.o $(OBJDIR)save_dirty.o $(OBJDIR)find_object.o $(OBJDIR)cut_frame.o \
$(OBJDIR)badpixel_correction.o $(OBJDIR)elaborate_sum.o




# all : OMEGA_pipeline test_function

test_function : test_function.c
	$(CC) $(CFLAGS) test_function.c -o ${O2K_UTIL}/pipeline/exe/test_function $(BIBS) 

OMEGA_pipeline : $(OBJLIST)
	$(CC) $(OBJ_PATH_LIST) -o  $(BIBS)


OMEGA_pipeline.o : $(SRCDIR)OMEGA_pipeline.c $(SRCDIR)pipe_para.h $(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h 
	$(CC) $(CFLAGS) -c $(SRCDIR)OMEGA_pipeline.c -o $(OBJDIR)OMEGA_pipeline.o



# create the secundary module object files with all dependancies

# Module 1
flat_correction.o : $(SRCDIR)flat_correction.c 		
	$(CC) $(CFLAGS) -c $(SRCDIR)flat_correction.c -o $(OBJDIR)flat_correction.o


# Module 2
subframe_statistics.o : $(SRCDIR)subframe_statistics.c $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h 
	$(CC) $(CFLAGS) -c $(SRCDIR)subframe_statistics.c -o $(OBJDIR)subframe_statistics.o


# Module 3
normalize_frame.o : $(SRCDIR)normalize_frame.c $(SRCDIR)pipe_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)normalize_frame.c -o $(OBJDIR)normalize_frame.o


# Module 4
extract_name.o : $(SRCDIR)extract_name.c $(SRCDIR)pipe_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)extract_name.c -o $(OBJDIR)extract_name.o


# Module 5
load_stack.o : $(SRCDIR)load_stack.c $(SRCDIR)pipe_para.h $(SRCDIR)OPRIME_para.h \
$(SRCDIR)O2000_para.h $(SRCDIR)struct_def.h
	$(CC) $(CFLAGS) -c $(SRCDIR)load_stack.c -o $(OBJDIR)load_stack.o


# Module 6
get_sky.o : $(SRCDIR)get_sky.c $(SRCDIR)pipe_para.h $(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h 
	$(CC) $(CFLAGS) -c $(SRCDIR)get_sky.c -o $(OBJDIR)get_sky.o


# Module 7
save_sky.o : $(SRCDIR)save_sky.c $(SRCDIR)pipe_para.h $(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h 
	$(CC) $(CFLAGS) -c $(SRCDIR)save_sky.c -o $(OBJDIR)save_sky.o


# Module 8
subtract_sky.o : $(SRCDIR)subtract_sky.c $(SRCDIR)pipe_para.h $(SRCDIR)OPRIME_para.h \
$(SRCDIR)O2000_para.h $(SRCDIR)struct_def.h
	$(CC) $(CFLAGS) -c $(SRCDIR)subtract_sky.c -o $(OBJDIR)subtract_sky.o


# Module 9
copy_frame.o : $(SRCDIR)copy_frame.c 		
	$(CC) $(CFLAGS) -c $(SRCDIR)copy_frame.c -o $(OBJDIR)copy_frame.o


# Module 10
original_level.o : $(SRCDIR)original_level.c 		
	$(CC) $(CFLAGS) -c $(SRCDIR)original_level.c -o $(OBJDIR)original_level.o


# Module 11
minimum_sky.o : $(SRCDIR)minimum_sky.c $(SRCDIR)pipe_para.h $(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)minimum_sky.c -o $(OBJDIR)minimum_sky.o


# Module 12
sky_clipping.o : $(SRCDIR)sky_clipping.c $(SRCDIR)pipe_para.h $(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)sky_clipping.c -o $(OBJDIR)sky_clipping.o


# Module 13
frame_info.o : $(SRCDIR)frame_info.c $(SRCDIR)struct_def.h
	$(CC) $(CFLAGS) -c $(SRCDIR)frame_info.c -o $(OBJDIR)frame_info.o


# Module 14
sum_preparation.o : $(SRCDIR)sum_preparation.c $(SRCDIR)struct_def.h  $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)sum_preparation.c -o $(OBJDIR)sum_preparation.o


# Module 15
make_cleansum.o : $(SRCDIR)make_cleansum.c $(SRCDIR)struct_def.h  $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)make_cleansum.c -o $(OBJDIR)make_cleansum.o


# Module 16
dirty_sum.o : $(SRCDIR)dirty_sum.c $(SRCDIR)struct_def.h  $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)dirty_sum.c -o $(OBJDIR)dirty_sum.o


# Module 17
put_tosum.o : $(SRCDIR)put_tosum.c $(SRCDIR)struct_def.h  $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)put_tosum.c -o $(OBJDIR)put_tosum.o


# Module 18
save_dirty.o : $(SRCDIR)save_dirty.c $(SRCDIR)struct_def.h  $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)save_dirty.c -o $(OBJDIR)save_dirty.o


# Module 19
find_object.o : $(SRCDIR)find_object.c $(SRCDIR)struct_def.h  $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)find_object.c -o $(OBJDIR)find_object.o


# Module 20
cut_frame.o : $(SRCDIR)cut_frame.c $(SRCDIR)struct_def.h  $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)cut_frame.c -o $(OBJDIR)cut_frame.o


# Module 21
badpixel_correction.o : $(SRCDIR)badpixel_correction.c $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)badpixel_correction.c -o $(OBJDIR)badpixel_correction.o


# Module 22
elaborate_sum.o : $(SRCDIR)elaborate_sum.c $(SRCDIR)struct_def.h  $(SRCDIR)pipe_para.h \
$(SRCDIR)OPRIME_para.h $(SRCDIR)O2000_para.h
	$(CC) $(CFLAGS) -c $(SRCDIR)elaborate_sum.c -o $(OBJDIR)elaborate_sum.o






# clean up command
clean :
	rm -f $(OBJDIR)*.o





