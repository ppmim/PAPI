#	Script zum Erstellen eines ausfuehrbaren Moduls 	17-Mar-93  HJR
#		for Calar Alto	(RF)
#_____________________________________________________________________________
#
f77 -e -fast -f -o {$argv[1]}.exe {$argv[1]}.f \
/calar1/photo/source/SUB/phot_sub.a \
/calar1/photo/source/LIB/phot_lib.a \
/calar1/photo/source/NUM_RECIP/NUM_RECIP.a \
-L${MID_LIB} \
-R${MID_LIB} \
-ldisp \
-lgen \
-lgmidas \
-lidi \
-lsubmid \
-lmidas \
-lsocket -lnsl \
/mpisys/lib/nag/libnag.a
#
if ( $status == 0 )  then
	echo "            ..... $argv[1] created and moved to user area"
endif
