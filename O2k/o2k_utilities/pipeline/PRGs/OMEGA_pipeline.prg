!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!.IDENT		OMEGA_pipeline.prg
!
!.CALL		@@ OMEGA_pipeline [image_cat] [frames] [sky_mode] [sum] ...
!				  ...[kappa_sum] [cuts] [flags] [flatfiels]
!
!.PARAMETERS	1) image_cat: name (and path if not local) of image catalog to be reduced
!			      automatic = online reduction with auto-icat extraction
!		   label: ICATALOG
!		2) frames: # of frames used for sky determination before/after target
!		   label: FRAMES		
!		3) sky_mode: mode for the sky determination
!		   0,n = minimum mode: take average of n smallest values as sky
!		   1 = median mode: take real median as sky
!		   2,k = outlier mode: elimate outliers first (more than k sigma), then
! 			 take median as sky
!		   label: MODE		Ex: mode=min,4  or  mode=2,3.5
!		4) sum: [n_average, action_flag, sum_save_flag]
!			parameters for summation of single images
!		   n_average =# of images for which median is taken in sum_process   
!		   action_flag: 0 = do reduction + summation of images;
!				1 = summation only ; 2 = reduction only
!		   sum_save_flag: 0 = write out every n_th sum-frame;
!				  1 = save only final master_sum;
!				  2 = save @ end: master_sum + real_sum + difference
!		   label:SUM
!		5) kappa_sum: parameter for the cosmics removal. In the summation process
!		              all values > med_level+(kappa_sum*sigma) are clipped off
!		   label: kappa_sum
!		6) cuts: values for LHCUTS of output frames in units of sigma
!		   label: CUTS	   Ex: 1,3 means CUT_MIN=med -1*sigma; CUT_MAX=med+3sigma
!		7) flags: save , position , screen_output
!		  a) save: flag for sky saving
!		   0 = sky will not be saved; 1 = sky for each frame will be saved
!		  b) position: flag for position where flatfield correction takes place
!		   0 = no correction; 1 = correction at start; 2 = correction at end 
!		  c) screen_output: flag for the amount of pipeline feedback on screen
!		   0=no output (log only); 1=display output
!		   label:FLAGS
!			
!		8) flatfield, badpixelmask, dark_frame: names of flatfield,
!			bad-pixel-mask and dark_frame to be used for reduction
!		   OR path and name of ascii file which contains the calibration frames. 
!		      A filename is indicated by a preceding "&". 
!		      Example: flatfield=&/disk-a/o2k/calibration.cal
!		   label: FLATFIELD
!
!.AUTHOR	Rene Fassbender, MPIA - Heidelberg
!.ENVIRONMENT	MIDAS					
!.MODULE	OMEGA2000 Pipeline
!.KEYWORDS	IR pipeline
!.PURPOSE	do pipeline setup and call OMEGA_pipeline.exe 
!.COMMENTS		
!.VERSION	1.00		26.08.02
!		1.1		18.10.02	added final parameters, setup for FITS
!		1.2 		30.10.02	implemented different sky_modes
!		1.3		13.11.02	implemented cuts for LHCUTS
!		2.0		20.11.02	implemented summation parameters
!		2.1		04.03.03	included badpixelmask
!		2.2		20.03.03	added P5 = kappa_sum
!		3.0		18.03.03	add dark_current frame
!		4.0		07.04.03	implemented online reduction preparation
!		4.1		08.04.03	allow icat with path information
!		5.0		12.06.03	final modifications: use environment
!						O2K_UTIL, use fire35 filepath,
!						allow ASCII-file for calibration info
!----------------------------------------------------------------------------------------



!+++++++++++++ 1 PARAMETERS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	! define labels for crossreference
crossref icatalog frames mode sum kappa_sum cuts flags flatfield

	! define parameters and their default values 
define/parameter P1 automatic C		! icat
define/parameter P2 3 N			! frames
define/parameter P3 1 ?			! sky mode, can be mixed para
define/parameter P4 7,0,0 N		! sum
define/parameter P5 10 N 		! kappa_sum
define/parameter P6 3,10 N		! default for LHCUTS
define/parameter P7 0,1,2 N		! flags
define/parameter P8 noinput C		! calibration frames
	



!++++++++++++++ 2 HELP TEXT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if P1(1:4) .eq. "help" then

write/out
write/out "This is the HELP PAGE for OMEGA_PIPELINE:"
write/out
write/out "call: @@ OMEGA_pipeline icatalog frames position sum kappa_sum cuts flags flatfield"
write/out
write/out "The command line parameters are:"
write/out "1) image_cat: name (and path if not local) of image catalog to be reduced"
write/out "		 automatic = online reduction with auto-icat extraction"
write/out "label: ICATALOG"
write/out "2) frames: # of frames used for sky determination before/after target"
write/out "label: FRAMES"		
write/out "3) sky_mode: mode for the sky determination"
write/out "   0,n = minimum mode: take average of n smallest values as sky"
write/out "   1 = median mode: take real median as sky"
write/out "   2,k = outlier mode: elimate outliers first (more than k sigma), then"
write/out "	   take median as sky"
write/out "   label: MODE" 		Ex: mode=min,4  or  mode=2,3.5	
write/out "4) sum: [n_average, action_flag, sum_save_flag]"
write/out "        parameters for summation of single images"
write/out "   n_average =# of images for which median is taken in sum_process "  
write/out "   action_flag: 0 = do reduction + summation of images;"
write/out "                1 = summation only ; 2 = reduction only"
write/out "   sum_save_flag: 0 = write out every n_th sum-frame;"
write/out "                  1 = save only final master_sum;"
write/out "                  2 = save @ end: master_sum + real_sum + difference"
write/out "   label:SUM"
write/out "5) kappa_sum: parameter for the cosmics removal. In the summation process"
write/out "		 all values > med_level+(kappa_sum*sigma) are clipped off"
write/out "   label: kappa_sum"
write/out "6) cuts: values for LHCUTS of output frames in units of sigma"
write/out "   label: CUTS	Ex: 1,3 means CUT_MIN=med -1*sigma; CUT_MAX=med+3sigma"
write/out "7) flags: save , position , screen_output"
write/out "  a) save: flag for sky saving"
write/out "     0 = sky will not be saved; 1 = sky for each frame will be saved"
write/out "  b) position: flag for position where flatfield correction takes place:"
write/out "     0 = no correction; 1 = correction at start; 2 = correction at end" 
write/out "  c) screen_output: flag for the amount of pipeline feedback on screen"
write/out "     0=no output (log only); 1=display output"
write/out "     label:FLAGS"
write/out
write/out "8) flatfield, badpixelmask, dark_frame: names of flatfield,bad-pixel-mask"
write/out "                                  and dark_frame to be used for reduction"
write/out "   OR path and name of ascii file which contains the calibration frames." 
write/out "      A filename is indicated by a preceding & ." 
write/out "      Example: flatfield=&/disk-a/o2k/calibration.cal"
write/out "   label: FLATFIELD"
write/out

return			! return to MIDAS session
endif

	


!++++++++++++++ 3 PIPELINE SETUP +++++++++++++++++++++++++++++++++++++++++++++++++++++++


	! absolute path were online information files are stored
	! define keyword which contains the path where all relevant files are
define/local file_path/c/1/60 "/disk-a/o2k/tmp"



	! Write number parameters into I/R-keywords to make them visible to the pipeline
define/local frames_key/I/1/1 {P2} ? +lower_levels

	! Write several keywords to make different MODE (P3) inputs possible 
define/local mode_key_n/I/1/1 4 ? +lower_levels		!default for n: 4
define/local mode_key_kappa/R/1/1 3 ? +lower_levels	!default for kappa: 3
define/local mode_key_sky/I/1/1 1 ? +lower_levels	!default for sky_mode: 1
define/local mode_key_i/I/1/2 0 ? +lower_levels	!
define/local mode_key_r/R/1/2 0 ? +lower_levels	!

	! Write keywords for the summation process
define/local p4_key/I/1/3 {P4} ? +lower_levels
define/local n_average_key/I/1/1
define/local action_flag_key/I/1/1
define/local sum_save_key/I/1/1

n_average_key = p4_key(1)
action_flag_key = p4_key(2)
sum_save_key = p4_key(3)


define/local kappa_sum_key/R/1/1 {P5} ? +lower_levels
define/local cuts_key/R/1/2 {P6} ? +lower_levels		
define/local flag_key/I/1/3 {P7} ? +lower_levels


!..........................ADJUST DEFAULTS
	! Calibration frames with defaults --> will be used as default in case
	! not all frames were explicitly specified in Case 3 (see below)
define/local flatfield_key/C/1/80 "norm_flat.fits" ? +lower_levels
define/local badpixel_key/C/1/80 "badpixelmask.fits" ? +lower_levels
define/local darkframe_key/C/1/80 "dark_extract.fits" ? +lower_levels

	! Keyword for name of ascii calibration file
define/local cal_file/C/1/100 "O2K_UTIL:/pipeline/CAL/calibration.cal" 


!........................................



	! Define additional variables
define/local komma_index/I/1/1 1			! index variable
define/local komma_index_2/I/1/1 1			! index variable
define/local two_kommas/I/1/1 0				! flag for 2 kommas
define/local slash_index/I/1/1				! /-position in icat

	!flag for online reduction: 0=no ; 1=yes
define/local online_flag/i/1/1 0 ? +lower_levels	! flag for online reduction 

	! flag for additional path information in icat: 0=no path ; 1=with path
define/local path_flag/I/1/1 0 ? +lower_levels


	! keywords for icat extraction and online preparation
	! more keywords for icat, path, integration time
define/local last_icat/c/1/80 {P1} ? +lower_levels	! for icat with path
define/local last_icat_name/c/1/80 {P1} ? +lower_levels	! for icat without path
define/local int_time/r/1/1 0  ?  +lower_levels		! integration time
define/local time_aux/c/1/10 				! auxiliary buffer




	! Set Midas to FITS mode: all new frames are in .fits format and can be updated
set/midas newfile=fits f_update=yes


	! Set up screen output option
if flag_key(3) .eq. 0 then	!---------

write/out
write/out "Screen output is set to LOGONLY"
write/out

	! screen output is only logged
set/midas output=logonly
endif				!---------



	! Extract sky mode from command line parameter
	! change keyword P3 to lower case
P3 = m$lower(P3)

	! test wether sky_mode = minimum
if m$index(P3,"0") .eq. 1 then			!........
write/key mode_key_i/I/1/2 {P3}
mode_key_sky = 0
mode_key_n = mode_key_i(2)
endif						!........

	! test wether sky_mode = outlier
if m$index(P3,"2") .eq. 1 then			!,,,,,,,,
write/key mode_key_r/R/1/2 {P3}
mode_key_sky = 2
mode_key_kappa = mode_key_r(2)
endif						!,,,,,,,,

	! test for string "min"
if m$index(P3,"min") .ge. 1 then		!........
mode_key_sky = 0
komma_index = m$index(P3,",")			! look for the ','
komma_index = komma_index+1

	if komma_index .ge. 3 then
	write/key mode_key_n {P3({komma_index}:{komma_index})}
	endif

endif						!........


	! test for string "out"
if m$index(P3,"out") .ge. 1 then		!--------
mode_key_sky = 2
komma_index = m$index(P3,",")			! look for the ','
komma_index = komma_index+1

	if komma_index .ge. 3 then
	write/key mode_key_kappa {P3({komma_index}:)}
	endif

endif						!--------




	! check wether icat contains path-information
slash_index = m$indexb(P1,"/")		! last slash in icat


if slash_index .ge. 1 then		! if / was found

	! set path_flag to 1
  write/keyword path_flag 1
	
	! save whole path with name
  last_icat = P1

	! set slash_index to first character if icat name
  slash_index = slash_index+1

	! extract catalog name only
  P1 = P1({slash_index}:) 

endif




!+++++++++++++++ 3b ONLINE REDUCTION PREPERATION +++++++++++++++++++++++++++++++++++++++

	
if P1(1:4) .eq. "auto" .and. path_flag .eq. 0 then	! online if icat=auto

	! file IDs
  define/local file_id_1/i/1/2				! file id for open/file
  define/local file_id_2/i/1/2				! file id for open/file
  define/local file_id_3/i/1/2				! file id for open/file


	! set keyword online_flag to 1
  write/keyword online_flag 1



	! read icat, path, integration time from files
	! file contains path and name of active image catalog
  open/file {file_path}/MacroLstIcat read file_id_1
  read/file {file_id_1(1)} last_icat
  close/file {file_id_1(1)}


	! file which contains the integration time per image
  open/file {file_path}/MacroIntTime read file_id_2
  read/file {file_id_2(1)} time_aux
  close/file {file_id_2(1)}

!  write/keyword int_time {time_aux}


	! file contains only name of active image catalog
  open/file {file_path}/MacroLstIcatName read file_id_3
  read/file {file_id_3(1)} last_icat_name
  close/file {file_id_3(1)}


	! write icat name in P1
  write/keyword P1 {last_icat_name}



	! get automatic flatfield: open frame, check filter, write right name in key
	! implemented below

endif





!§§§§§§§§§§§§§§§§ 3c CALIBRATION FRAMES §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

	! define flag for proper calibration frame selection
define/local calsel_flag/i/1/1 0	! 0 = start value ; 1 = noinput given ; 2 = file
define/local line_buffer/c/1/100	! buffer to read lines from file
define/local calfile_id/i/1/2		! ID for calibration file
define/local filter_key/c/1/10		! keyword for filter name
define/local first_image/c/1/100	! keyword for name of first image
define/local loop_cal/i/1/1		! loop variable
define/local cal_check/i/1/1 0		! test variable
define/local line_test/i/1/1 0		! count read in file lines
define/local line_test_2/i/1/1 0	! count read in file lines



	! DEFAULT
	! Case 1: noinput is given --> use default calibration file
	
if P8(1:7) .eq. "noinput" then
	 
  calsel_flag = 1	! set flag to 1

endif


	! Case 2: a calibration ascii file is specified --> extract calibration frames
if P8(1:1) .eq. "&" then
	 
  calsel_flag = 2	! set flag to 2

	! write calibration file name in keyword, flush first
  write/keyword cal_file "                                                            "
  write/key cal_file {P8(2:)}	! start with second character (& is first)


endif


	!-----------------------------------------------------------------------------
	! read file for case 1 and 2
if calsel_flag .gt. 0 then		

	! get filter information from first entry in image catalog
  set/icat {last_icat}

    copy/dkey #1 filter filter_key/c/1/10

  clear/icat



	! open calibration file
  open/file {cal_file} read calfile_id

	! check return value
  if calfile_id(1) .lt. 0 then
    write/out "WARNING: Error occurred while opening the calibration file" 
  endif



	!###########################################
  do loop_cal = 1 100 1		! read max 100 lines

	! reset line buffer
    write/keyword line_buffer " " all

	! read in lines
    read/file {calfile_id(1)} line_buffer


    line_test = line_test+1				! count read in lines


    if calfile_id(2) .lt. 0 goto out_loop		! quit loop if end of file

    if line_buffer(1:1) .eq. "!" goto end_loop		! comment line


    line_test_2 = line_test_2+1				! test 


	! check lines for different calibration frames 
    if line_buffer(1:6) .eq. "BPM = " then

	write/out "BPM found..."
   	write/key badpixel_key/c/1/80 {line_buffer(7:)}
	cal_check = cal_check+1	! record found calibration name
	goto end_loop		! continue with next line


    elseif line_buffer(1:7) .eq. "DARK = " then
	
	write/out "DARK found..."
	write/key darkframe_key/c/1/80 {line_buffer(8:)}	
	cal_check = cal_check+1	! record found calibration name
	goto end_loop		! continue with next line


    elseif line_buffer(1:7) .eq. "FLAT = " then
	
	write/out "FLAT found..."
	write/key flatfield_key/c/1/80 {line_buffer(8:)}
	cal_check = cal_check+1	! record found calibration name
	goto out_loop		! quit loop after FLAT was found
	

    elseif line_buffer(1:10) .eq. "{filter_key(1:10)}" then  

	write/out "Flatfield for current FILTER found..."
	write/key flatfield_key/c/1/80 {line_buffer(14:)}	  
	cal_check = cal_check+1	! record found calibration name
	goto out_loop		! quit loop after filter flat was found

    endif



    end_loop:			! continue with next line
  enddo
	!###########################################


  out_loop:			! label to quit loop
  close/file {calfile_id(1)}	! close file


  if cal_check .ne. 3 then	! #calibration frames differs from 3

	write/out "WARNING: The number of calibration frames identified is not 3 ..."

  endif


set/format I1

write/out "{line_test} lines were read in from file,{line_test_2} were not comment lines" 
write/out 
write/out 


endif 	! end of file reading
	!----------------------------------------------------------------------------



	! Case 3: the calibration frames are given as parameter
	! Extract names of calibration frames
	! If not input was given use default values 
if P8(1:7) .ne. "noinput" .and. calsel_flag .eq. 0 then

	! check for komma
	komma_index = m$index(P8,",")		! first komma
	komma_index_2 = m$indexb(P8,",")	! last komma

	if komma_index .ne. komma_index_2 then
	  two_kommas = 1			! set two_komma flag if kommas different
	endif


	if komma_index .ge. 1 .and. komma_index .le. 81 then
		! komma found

	   if P8(1:1) .ne. "?" then			! not default
		komma_index = komma_index-1		! set to last character
							! of flat_name	
		flatfield_key = P8(1:{komma_index})	! new flatfield name	
	   endif


		! reset komma_index
	   komma_index = m$index(P8,",")
	   komma_index = komma_index+1		! set to first character of badpixel name


	   if P8({komma_index}:{komma_index}) .ne. "?" then 	! not default

		if two_kommas .eq. 0 then		! last calibration frame name
		  badpixel_key = P8({komma_index}:)	! new bad-pixel name
		else
		  komma_index_2 = komma_index_2-1	! set to last character
		  badpixel_key = P8({komma_index}:{komma_index_2})
		endif
		
	   endif


		! now dark_current frame
	  if two_kommas .eq. 1 then			! second komma was detected

		komma_index_2 = m$indexb(P8,",")	! last komma
		komma_index_2 = komma_index_2+1		! set to first character

		if P8({komma_index_2}:{komma_index_2}) .ne. "?" then 	! not default
		  darkframe_key = P8({komma_index_2}:)
		endif

	  endif


	else
	   flatfield_key = P8		! input is flatfield name
		
	endif

endif

!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§







!+++++++++++++++ 4 RUN PIPELINE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! as long as pipeline is running, the PRG has no control

!run /disk-v/fassbend/Programme/exe/omega_pipeline.exe
run O2K_UTIL:/pipeline/exe/omega_pipeline.exe



!+++++++++++++++ 5 SHUT DOWN PIPELINE ++++++++++++++++++++++++++++++++++++++++++++++++++


write/out
write/out "Here is the PRG again.... "
write/out


	! Set mode back to MIDAS default (bdf frames)
set/midas work_env=midas

	! Set back screen output
set/midas output=yes


