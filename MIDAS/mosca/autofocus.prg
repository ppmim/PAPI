!
! autofocus.prg
!  
! for taking automatic focus series
!
!  JWF oct 97 
!  modified Jul 98
!  modified Nov.98/U.T.


if p1(1:4) .eq. "help" then
	write/out @@ autofocus [1.focus value],[step],[#steps] [exptime]
	write/out
	write/out take automatically a focus series  
	write/out either MOSCA or the telescope can be focussed
	write/out outputfile = test0001
	write/out
	return
endif


defi/local device/c/1/20 "M                   "
inquire/key device "focus MOSCA [def] or TELESCOPE [M/T] ? "

device = m$upper(device)

defi/local text/c/1/10 "              "


define/par p2 ? n "exposure time ? "


! focus limits for MOSCA and telescope 
defi/local moscalim/r/1/2 6400.,9000.
defi/local tellim/r/1/2 20.,40.
defi/local foclim/r/1/1 0.

if device .eq. "M"  then
	write/out configure MOSCA ...
	$ /disk-a/staff/MOSCA/scripts/mosca_calib_mirror in wait
	$ /disk-a/staff/MOSCA/scripts/mosca_lamp_cont on 60% wait
	$ /disk-a/staff/MOSCA/scripts/mosca_mask 2 wait

	write/out configure CCD ...
        $ camera_dewar1_exptime {p2}
        $ camera_dewar1_exptype science
        $ camera_expmode test
        write/out select filter (normally filter+grism = free) 
	write/out select area on CCD (for SITE16a:x1=800 x2=1200 y1=1400 y2=2400 )
	inquire/key text " Hit <CR> when ready ! "
	write/out
	write/out NOTE: focus and step must be in micrometers (eg. 8300,30,7) !
	write/out

	define/par p1 ? n "expected focus, step, # steps = ? "
	defi/local foc/r/1/3 {p1}
	defi/local nfoc/i/1/1 {foc(3)}
	write/key focser/r/1/3 0.,0.,0.
	focser(1) = foc(1)-nfoc/2*foc(2)
	focser(2) = foc(2)
	focser(3) = nfoc
	
	
! check entered values 
	
	foclim = foc(1)-nfoc/2*foc(2)
	if foclim .lt. moscalim(1)  then
		write/out ERROR : focus to small - exiting ...
		$ audioplay -v 80 /disk-a/staff/MOSCA/mosca_prg/crash.au
		goto exit
	endif
	foclim = foc(1)+nfoc/2*foc(2)
	if foclim .gt. moscalim(2)  then
		write/out ERROR : focus to large - exiting ...
		$ audioplay -v 80 /disk-a/staff/MOSCA/mosca_prg/crash.au
		goto exit
	endif


 

endif



if device .eq. "T"  then
	write/out configure MOSCA ...
	$ /disk-a/staff/MOSCA/scripts/mosca_calib_mirror out wait
	$ /disk-a/staff/MOSCA/scripts/mosca_mask 5 wait
	write/out configure CCD ...
	defi/par p2 ? n "exposure time = ? "
        $ camera_dewar1_exptime {p2}
        $ camera_dewar1_exptype science
        $ camera_expmode test
        write/out Select filter and area of CCD ...
	inquire/key text " Hit <CR> when ready ! "
	write/out  Note: telescope focus must be in millimeters (eg. 30.4) 

	define/par p1 ? n "expected focus, step, # steps = ? "
	defi/local foc/r/1/3 {p1}
	defi/local nfoc/i/1/1 {foc(3)}
	write/key focser/r/1/3 0.,0.,0.
	focser(1) = foc(1)-nfoc/2*foc(2)
	focser(2) = foc(2)
	focser(3) = nfoc	
! check entered values 

	if foc(1) .gt. 40.0 then
		write/out ERROR: telescope focus must be in millimeters (eg. 30.4)... exiting
		$ audioplay -v 80 /disk-a/staff/MOSCA/mosca_prg/crash.au
		goto exit
	endif
	
	foclim = foc(1)-nfoc/2*foc(2)
	if foclim .lt. tellim(1)  then
		write/out ERROR : focus to small - exiting ...
		$ audioplay -v 80 /disk-a/staff/MOSCA/mosca_prg/crash.au
		goto exit
	endif
	foclim = foc(1)+nfoc/2*foc(2)
	if foclim .gt. tellim(2)  then
		write/out ERROR : focus to large - exiting ...
		$ audioplay -v 80 /disk-a/staff/MOSCA/mosca_prg/crash.au
		goto exit
	endif

	


endif





! set CCD

write/out Configuring CCD ...

$ camera_expmode test
$ camera_chip1_rowshifts 50
set/format f2.0
$ camera_dewar1_f_loops {foc(3)}

$ camera_dewar1_exptype focus
$ camera_objectname "focus series"



defi/local i/i/1/1 0
defi/local focnew/r/1/1 0



do i = 1 nfoc
	set/format f6.3
	focnew = focser(1)+(i-1)*focser(2)
	if device(1:1) .eq. "M" then
		if focnew .gt. moscalim(1) .and. focnew .lt. moscalim(2) then 
			write/out new MOSCA focus = {focnew}
			set/format f10.0	
			$/disk-a/staff/MOSCA/scripts/mosca_afocus {focnew}
			wait/secs 3.
		else 
			write/out ERROR: MOSCA focus out of range - exiting ...
		endif
	else
		if focnew .gt. tellim(1) .and. focnew .lt. tellim(2) then
			write/out new telescope focus = {focnew}
			set/format f6.3	
			$/disk-a/staff/TECS35/scripts/t3_afocus {focnew}
               		 wait/secs 3.
		else
			write/out ERROR: MOSCA focus out of range - exiting ...
		endif
		
        endif

	set/format i2
	write/out starting exposure {i} of {nfoc}
	$ camera_dewar1_start wait_fast
enddo

write/out 
write/out exposure done !
write/out
if device(1:1) .eq. "M" then 
	write/out a good observer always turns off the calibration lamps ...
	$ /disk-a/staff/MOSCA/scripts/mosca_lamp_cont off wait
endif

write/out
write/out Load image test0001 into display and 
write/out Use @@ focus to analyze this series!
write/out Do not forget to set CCD to full frame!!!
write/out

wait/secs 3.
! set camera to normal mode
$camera_expmode normal


exit:

