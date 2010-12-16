iraf=/disk-a/caha/iraf/irafbin/
arch=.suse
#-----------------------------------------------------------------------------
#
# SNAP 1.3 Makefile 
#
# Filippo Mannucci, Maj 2006
#
#-----------------------------------------------------------------------------
# parameters: edit the following lines

# 1- Full link to the FORTRAN compiler (IRAF might hide it).
#    Get it with, for example : which f77 
     FF = /usr/bin/gfortran
#    FF = /usr/bin/g77-33
     CC = /usr/bin/cc


# 2- IRAF STSDAS directory. 
#     Under IRAF get the stsdas base directory by the command: show stsdas
#     Examples:
#     STSDAS = /iraf/extern/stsdas#         	Linux, "Linux for Astronomy"
#     STSDAS = /scisoft/iraf/extern/stsdas/bin# Linux, Scisoft
      STSDAS = /scisoft/iraf/extern/stsdas#     Linux, Scisoft
#     STSDAS = /iraf/stsdas#                	SOLARIS
#     STSDAS = /usr/im/stsdas2.1.1/#			SOLARIS


# 3- If you need to install SExtractor, uncomment one of the following:
#     SEXMACHINE = aix		# (for IBMs RS6000 running AIX)
#     SEXMACHINE = alpha		# (for DEC-ALPHAs with Digital UNIX)
#     SEXMACHINE = hpux		# (for HP/UX systems)
#     SEXMACHINE = linuxpc		# (for PCs running LINUX, using gcc)
      SEXMACHINE = linuxp2		# (for Pentium2/3/4 PCs running LINUX)
#     SEXMACHINE = linuxk7		# (for Athlon PCs running LINUX)
#     SEXMACHINE = sgi		# (for SGI platforms)
#     SEXMACHINE = solaris		# (for SUN-Solaris machines)
#     SEXMACHINE = sunos		# (for SUN-OS platforms)
#     SEXMACHINE = ultrix		# (for DEC stations running ULTRI)



# END OF SECTION TO EDIT
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# other definitions

RM       = rm -f
SNAPDIR  = $(PWD)
IRAFDIR  = $(iraf)
ARCHIT   = $(arch)
IRAFBIN  = $(IRAFDIR)bin$(ARCHIT)
STSDASLIB = $(STSDAS)/bin$(ARCHIT)
  

include Def/Definitions$(ARCHIT)


#-------------------------------------------------------------------
# targets

default: all

show:
	@echo IRAF=$(iraf)
	@echo STSDAS=$(STSDAS)
	@echo FF=$(FF)
	@echo CC=$(CC)

#all: clean crt_nics newirdr snap nics arnica isaac sofi cl 
all: clean  crt_nics newirdr snap nics arnica isaac sofi cl

clean:
	$(RM) ./irdr/bin/*
	$(RM) ./snap
	$(RM) ./*.cl

crt_nics: 
	$(FF) -o crt_nics ./src/crt_nics$(ARCHIT).f 

newirdr:    
	cd irdr; make all CC=$(CC)

snap: 
	chmod u+x clearsnap

nics: generic/nics.generic
	$(AWK) '{if (gsub(/SNAPDIR/, "$(SNAPDIR)")>=0) print $0}' \
	generic/nics.generic | \
	$(AWK) '{if (gsub(/IRAFDIR/, "$(IRAFDIR)")>=0) print $0}' | \
	$(AWK) '{if (gsub(/AWK/, "$(AWK)")>=0) print $0}' | \
	$(AWK) '{if (gsub(/ARCHIT/, "$(ARCHIT)")>=0) print $0}'   | \
	$(AWK) '{if (gsub(/IRAFARCH ./, "IRAFARCH ")>=0) print $0}' \
	> nics
	chmod u+x nics
	cp nics snap

arnica: generic/arnica.generic
	$(AWK) '{if (gsub(/SNAPDIR/, "$(SNAPDIR)")>=0) print $0}' \
	generic/arnica.generic | \
	$(AWK) '{if (gsub(/IRAFDIR/, "$(IRAFDIR)")>=0) print $0}' | \
	$(AWK) '{if (gsub(/AWK/, "$(AWK)")>=0) print $0}' | \
	$(AWK) '{if (gsub(/ARCHIT/, "$(ARCHIT)")>=0) print $0}'   | \
	$(AWK) '{if (gsub(/IRAFARCH ./, "IRAFARCH ")>=0) print $0}' \
	> arnica
	chmod u+x arnica

isaac: generic/isaac.generic
	$(AWK) '{if (gsub(/SNAPDIR/, "$(SNAPDIR)")>=0) print $0}' \
	generic/isaac.generic | \
	$(AWK) '{if (gsub(/IRAFDIR/, "$(IRAFDIR)")>=0) print $0}' | \
	$(AWK) '{if (gsub(/AWK/, "$(AWK)")>=0) print $0}' | \
	$(AWK) '{if (gsub(/ARCHIT/, "$(ARCHIT)")>=0) print $0}'   | \
	$(AWK) '{if (gsub(/IRAFARCH ./, "IRAFARCH ")>=0) print $0}' \
	> isaac
	chmod u+x isaac

sofi: generic/sofi.generic
	$(AWK) '{if (gsub(/SNAPDIR/, "$(SNAPDIR)")>=0) print $0}' \
	generic/sofi.generic | \
	$(AWK) '{if (gsub(/IRAFDIR/, "$(IRAFDIR)")>=0) print $0}' | \
	$(AWK) '{if (gsub(/AWK/, "$(AWK)")>=0) print $0}' | \
	$(AWK) '{if (gsub(/ARCHIT/, "$(ARCHIT)")>=0) print $0}'   | \
	$(AWK) '{if (gsub(/IRAFARCH ./, "IRAFARCH ")>=0) print $0}' \
	> sofi
	chmod u+x sofi

cl: 
	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/objmask.cl.generic > objmask.cl
	chmod u+x objmask.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/head-nics.cl.generic > head-nics.cl
	chmod u+x head-nics.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/head-arnica.cl.generic > head-arnica.cl
	chmod u+x head-arnica.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/head-isaac.cl.generic > head-isaac.cl
	chmod u+x head-isaac.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/mkflat.cl.generic > mkflat.cl
	chmod u+x mkflat.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/distor.cl.generic | \
	$(AWK) '{if (gsub(/SNAPDIR/, "$(SNAPDIR)")>=0) print $0}' > distor.cl
	chmod u+x distor.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/imrep.cl.generic > imrep.cl
	chmod u+x imrep.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/dimtrim.cl.generic > dimtrim.cl
	chmod u+x dimtrim.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/imtrim.cl.generic > imtrim.cl
	chmod u+x imtrim.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/dgain.cl.generic > dgain.cl
	chmod u+x dgain.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/imcomb.cl.generic > imcomb.cl
	chmod u+x imcomb.cl

	$(AWK) '{if (gsub(/IRAFBIN/, "$(IRAFBIN)")>=0) print $0}' \
	generic/addmask.cl.generic | \
	$(AWK) '{if (gsub(/SNAPDIR/, "$(SNAPDIR)")>=0) print $0}' > addmask.cl
	chmod u+x addmask.cl


sex:
	tar xvf sex_2.2.2.tar
	cd sextractor2.2.2/source ; \
	make SEXMACHINE=$(SEXMACHINE) ; \
	mv sex ../.. ; \
	make clean 
	rm -r sextractor2.2.2; \



