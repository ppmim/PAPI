!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! .COPYRIGHT	(C) MPIA
! .IDENT	domeflats.prg
!
! .CALL		o2k/domeflats level images object lamp_start
!				  
! .PARAMETERS
!		P1: total countlevel of all images (lamp on)
!		P2: flag for image saving --> 0=only one sum image is saved
!					      1=all single images are saved
!		P3:object = object name for identifier
!		P4:lamp_start = flatfield lamp to be tested first: 1-5 
!		                for lamp 5: 5,level = 0,..,9		
!
! .PURPOSE	Take dome flatfields with OMEGA2000
! .ENVIRONMENT	MIDAS
! .AUTHOR	Rene Fassbender     
! .KEYWORDS	flatfields
! .COMMENTS	Note: The flatfield lamp GUI "ffl" has to be inactive when 
!		      starting the PRG			
!
! .CONVENTIONS	!! --> command of line necessary
!		!* --> command not needed but could be useful for testing
!
! .VERSION	1.0	16.03.03
!		1.1	10.03.03	added sync before last in test-run
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!!!!!!!!!!!!!!!!!! SETUP AND HELP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! command line parameter setup
crossref level images object lamp_start


	! help text
if p1(1:4) .eq. "help" then

write/out
write/out "domeflats.prg"
write/out 
write/out "call: o2k/domeflats level images object lamp_start"
write/out
write/out "This program takes dome flatfields for OMEGA2000:"
write/out 
write/out "The command line parameters are:"
write/out "P1: total countlevel of all images with lamp on"
write/out "P2: flag for image saving --> 0=only one sum image is saved"
write/out "			         1=all single images are save"
write/out "P3:object = object name for identifier"
write/out "P4:lamp_start = flatfield lamp to be tested first: 1-5" 
write/out "		   for lamp 5: 5,level = 0,..,9"
write/out		
write/out "Note: The flatfield-lamp-GUI "ffl" has to be inactive when" 
write/out "      starting the PRG"		
write/out

return		! back to MIDAS session

endif	



	! define command line parameters
define/parameter P1 ? N "Enter the total count level:"
define/parameter P2 1 N "Enter the flag for image saving:"
define/parameter P3 ? C "Enter an object name:"
define/parameter P4 5,5 N "Enter the flatfield lamp that is tested first:"


	! write parameters in keywords
define/local countlevel/r/1/1 {P1}
define/local image_flag/i/1/1 {P2}
define/local lamp_aux/i/1/2 {P4}
define/local start_lamp/i/1/1 {lamp_aux(1)}
define/local start_level/i/1/1  {lamp_aux(2)}

if start_lamp .le. 0 .or. start_lamp .ge. 6 then	! only lamps 1-5 allowed
  start_lamp = 5
  start_level = 5
  write/out "The entered lamp was not allowed. Default 5,5 will be used..."
endif


set/format I1 F4.2

	! display input parameters
write/out
write/out "Flatfields will now be taken with the current filter..."	
write/out "The total countlevel  will be about: {countlevel}"
write/out "The flag for image saving is (0=one ; 1=all): {image_flag}"
write/out "The identifier will contain the object name:  {P3}"
write/out "The first flatfield lamp to be tested will be: {start_lamp}"
write/out

set/format

!
!!!!!!!!!!!!!!!!!! KEYWORDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! INTERNAL PARAMETERS
	! wanted countlevel for single exposure: min,max 
define/local ok_level/r/1/2 1000,20000			! ok level
define/local good_level/r/1/2 2000,10000		! good level,<=best_level 
define/local best_level/r/1/1 10000			! norm level



	! measured level
define/local actual_level/r/1/5			! store levels for each test
define/local break_flag/i/1/1 0			! if = 1 break test image taking
define/local lampfive_flag/i/1/1 0		! if = 1 lamp 5 is the right lamp
define/local reached_level/r/1/1		! measured level

	! name and path of frame
define/local pathname_ima/c/1/200

	! integration times
define/local single_time/r/1/1
define/local act_time/r/1/1
define/local repetitions/i/1/1
	
	! for return value of abort-file-check
define/local abort_check/i/1/1		

	! lamp and level for lamp 5
define/local lamp/i/1/1 
define/local lamp_level/i/1/1 5

	! loop keywords
define/local loop/i/1/1 
	


	! remove existing abort file
	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("{geirslstabort}")
if abort_check .eq. 1 then
  $rm {geirslstabort} 
endif





!!!!!!!!!!!!!!!!! FIND LAMP AND INTEGRATION TIME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write/out
write/out "Test images are now taken to determine the proper flatfield lamp..."
write/out

set/format I1

	! set camera parameters
$cmd_panic_new crep 1				! 1 image
$cmd_panic_new itime 1				! 1 second
$cmd_panic_new object test images
$cmd_panic_new sync

lamp = start_lamp					! initialize lamp
lamp_level = start_level				! initialize lamp level



!----------------------------------------------------------------------------------
	! take test images
do loop = 1 5					! 5 loops is max	

	! check lampfive flag
  if lampfive_flag .ge. 1 then			! lamp 5 is right
    lamp = 5
  endif

  write/out "Flatfield lamp {lamp} is tested..."

	! test lamps 5,4,3,2,1
  if lamp .eq. 5 then
    $ flats L{lamp} ON {lamp_level}
    write/out "...with lamp level {lamp_level}"
  else
    $ flats L{lamp} ON
  endif

  wait/secs 1					! wait 1 second for lamp


  $cmd_panic_new read

  $cmd_panic_new sync


	! abort check: 0=does not exist ; 1=abort file exists
  abort_check = m$exist("{geirslstabort}")
  if abort_check .eq. 1 then
    write/out "Program is aborted..."
    $rm {geirslstabort}
    return
  endif


  $cmd_panic_new save -c testframe.fits

  $cmd_panic_new sync


  $ flats ALLOFF				! turns all lamps off
 

	! get file name and path	
  $cmd_panic_new last	! writes last filename in file geirsLstFile
	! writes last filename in keyword pathname_ima
  write/keyword pathname_ima </disk-a/o2k/tmp/geirsLstFile 


	! do statistics in test image
  set/midas output=logonly
  statistics/image {pathname_ima} [600,600:800,800] ? ? RNN 	! reduced stat.
  actual_level({loop}) = outputr(3)	! mean value in subframe
  set/midas output=yes
  write/out "Current test level:{actual_level({loop})}"


	! ****************** start level analysis *************************************
	! analyze level
  if actual_level({loop}) .lt. ok_level(1) then 	! 1) <minimum �����������������

    if lamp .ge. 2 then				! 1a) if higher lamp is available
      lamp = lamp-1				! try next higher lamp
    else					! 1b) lamp 1 is already used
      break_flag = 1
      lamp = 1					! lamp 1 is right
    endif


  elseif actual_level({loop}) .gt. ok_level(2) then	! 2) >maximum ����������������

    if lamp .eq. 5 then				! 2a) already lowest lamp
      lamp_level = lamp_level-1			! decrease lamp level
      lampfive_flag = 1			! level 5 is right
    elseif lamp .eq. 4 then			! 2b) 4 is too high
      lamp_level = 5				! increase lamp level
      lampfive_flag = 1			! level 5    
    else					! 2c) 1,2,3 are too high
      lamp = lamp+1				! go to weaker lamp
 
      if loop .ge. 2 then			! lower lamp already tested
        loop = loop-1				! set level postion right
        break_flag = 1				! set break flag
      endif
 
    endif


  else						! 3) level in specified ok_area ������

    if lamp .eq. 5 then				! 3a) lamp five is right
      lampfive_flag = 1

      if actual_level({loop}) .gt. good_level(2) then	! > ideal
        lamp_level = lamp_level-1			! go one level down
      elseif actual_level({loop}) .lt. good_level(1) then	! < ideal
        lamp_level = lamp_level+1				! go one level up
      else 
        break_flag = 1					! break
      endif

    else						! 3b) other than lamp 5 is ok

      if actual_level({loop}) .lt. good_level(2) then 	! real good level
        break_flag = 1					! break
      else
        lamp = lamp+1					! take next weaker lamp
 
        if loop .ge. 2 then				! if lower lamp was tested 
          loop = loop-1					! set level postion right
          break_flag = 1
        endif

      endif


    endif

  endif
	!********************* end of level analysis *******************************


	! break if good level was found
  if break_flag .eq. 1 then
    goto loopbreak				! break
  endif


enddo
!----------------------------------------------------------------------------------


loopbreak:					! break label


	! delete test image
$rm {pathname_ima}

$cmd_panic_new itime -stdout | write/key act_time/r/1/1

write/out
write/out "Flatfield lamp determination done..."
write/out


	! calculate integration time for single image
single_time = best_level/actual_level({loop})*act_time

	! condition should only be true for lamp 5
if single_time .lt. 1 then			! int time too short
  lamp = 5
  lamp_level = 0
  single_time = 1				! minimum integration
endif

set/format I1 f6.1

write/out
write/out "The appropriate lamp for the current filter is: {lamp}"

if lamp .eq. 5 then
  write/out "...the lamp level is:{lamp_level}"
endif




!!!!!!!!!!!!!!!!!!!!! TAKE FLATFIELDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


write/out "The integration time for a single integration is set to:{single_time}"
write/out "The flatfield lamp is now turned on and the integration starts..."


	! turn proper lamp on	
if lamp .eq. 5 then				! lamp 5 was picked
  $ flats L{lamp} ON {lamp_level}
else
  $ flats L{lamp} ON 				! lamp 1-4 was picked
endif

wait/secs 2

	! calculate repititions
repetitions = countlevel/best_level


	! write initial image descriptors
$cmd_panic_new counter DITH_NO clear	! clear dither counter 
$cmd_panic_new counter POINT_NO clear	! clear pointing no 
$cmd_panic_new counter EXPO_NO clear	! clear exposure counter-->EXPO_NO=1



	! set single image parameters
set/format I1

if image_flag .eq. 0 then			! one integrated image
  $cmd_panic_new crep {repetitions}
else					! save single images
  $cmd_panic_new crep 1
endif


$cmd_panic_new itime {single_time}
$cmd_panic_new sync




!+++++++++++++++++++ single flatfield image light on  +++++++++++++++++++++++++++++
if image_flag .eq. 0 then			! one integrated image

	! write object
  $cmd_panic_new object {P3} on

  $cmd_panic_new read

  $cmd_panic_new sync read

  $cmd_panic_new save -i

  $cmd_panic_new sync


	! get file name and path	
  $cmd_panic_new last	! writes last filename in file geirsLstFile
	! writes last filename in keyword pathname_ima
  write/keyword pathname_ima </disk-a/o2k/tmp/geirsLstFile 


	! do statistics in test image
  set/midas output=logonly
  statistics/image {pathname_ima} [600,600:800,800] ? ? RNN 	! reduced stat.
  reached_level = outputr(3)				! mean value in subframe
  set/midas output=yes

  write/out
  write/out "One integrated flatfield was taken..."
  write/out "...with actual countlevel: {reached_level}"
  write/out



!++++++++++++++++++ several flatfield images +++++++++++++++++++++++++++++++++++++
else					! save single images


!-------------------------------------------------------------
	! do loop over number of images
  do loop = 1 {repetitions} 1

    set/format I1
    write/out "Taking flatfield {loop} of {repetitions}..."	

	! write object
    $cmd_panic_new object {P3}:{loop}/{repetitions} on

    $cmd_panic_new read

    $cmd_panic_new sync


	! abort check: 0=does not exist ; 1=abort file exists
    abort_check = m$exist("{geirslstabort}")
    if abort_check .eq. 1 then
      write/out "Program is aborted..."
      $rm {geirslstabort}
      return
    endif

 
   $cmd_panic_new save

  enddo
!-------------------------------------------------------------

 
  $cmd_panic_new sync


	! get file name and path	
  $cmd_panic_new last	! writes last filename in file geirsLstFile
	! writes last filename in keyword pathname_ima
  write/keyword pathname_ima </disk-a/o2k/tmp/geirsLstFile 


	! do statistics in test image
  set/midas output=logonly
  statistics/image {pathname_ima} [600,600:800,800] ? ? RNN 	! reduced stat.
  reached_level = outputr(3)				! mean value in subframe
  set/midas output=yes

  set/format I1 f6.1

  write/out
  write/out "{repetitions} flatfields were taken..."
  write/out "...with actual countlevel for each flatfield: {reached_level}"
  write/out
  write/out "Now taking same sequence with light off ..."

endif


$ flats ALLOFF				! turns all lamps off

wait/secs 5

!+++++++++++++++++++ single flatfield image light off  ++++++++++++++++++++++++++++
if image_flag .eq. 0 then			! one integrated image

	! write object
  $cmd_panic_new object {P3} off

  $cmd_panic_new read

  $cmd_panic_new sync read

  $cmd_panic_new save -i

  $cmd_panic_new sync


	! get file name and path	
  $cmd_panic_new last	! writes last filename in file geirsLstFile
	! writes last filename in keyword pathname_ima
  write/keyword pathname_ima </disk-a/o2k/tmp/geirsLstFile 


	! do statistics in test image
  set/midas output=logonly
  statistics/image {pathname_ima} [600,600:800,800] ? ? RNN 	! reduced stat.
  reached_level = outputr(3)				! mean value in subframe
  set/midas output=yes

  write/out
  write/out "One integrated flatfield was taken with lamp off ..."
  write/out "                          ...with actual countlevel: {reached_level}"
  write/out



!++++++++++++++++++ several flatfield images +++++++++++++++++++++++++++++++++++++
else					! save single images


!-------------------------------------------------------------
	! do loop over number of images
  do loop = 1 {repetitions} 1

    set/format I1
    write/out "Taking flatfield {loop} of {repetitions} with lamp off ..."	

	! write object
    $cmd_panic_new object {P3}:{loop}/{repetitions} off

    $cmd_panic_new read

    $cmd_panic_new sync


	! abort check: 0=does not exist ; 1=abort file exists
    abort_check = m$exist("{geirslstabort}")
    if abort_check .eq. 1 then
      write/out "Program is aborted..."
      $rm {geirslstabort}
      return
    endif

 
   $cmd_panic_new save

  enddo
!-------------------------------------------------------------

 
  $cmd_panic_new sync



!!!!!!!!!!!!!!!!! CLOSE DOWN THINGS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write/out
write/out "Flatfields are now finished..." 
write/out


