!	s k y f l a t s . p r g		28-Sep-04	HJR
!
!  Series of flatfield exposures on the sky including dawn
!

crossref ident level coadds number offsets check saturation time

define/par P1 ? ? "Identifier: "
if p1(1:4) .eq. "help" then
   write/out
   write/out "skyflats.prg"
   write/out "call: o2k/skyflats [ident] [max_counts] [n_exp] [n_coadd]" 
   write/out "	      [read-mode] [dRA,dDEC]  [saturation]  [time]"
   write/out 
   write/out "The command line parameters are:"
   write/out "  P1 = identifier (use quotes if blanks are included)"
   write/out "  P2 = maximum exposure level in single exposure"
   write/out "  P3 = number of frames to be added up in memory before image is stored"
   write/out "  P4 = number of images to be stored"
   write/out "  P5 = telescope offsets [arcsec] in RA,DEC between stored frames"
   write/out "  P6 = check exposure level after each image [check/no_check]"
   write/out "  P7 = saturation level, monitoring level"
   write/out "  P8 = time of day (dusk / dawn)"
   write/out "To abort a sequence press ABORT in the GEIRS gui! Never use ^C"
   write/out
   goto exit
endif

define/par P2 20000 N/A "Maximum exposure level per frame [counts] : "
define/par P3 2 N/A "Number of co-adds : "
define/par P4 5 N/A "Number of images : "
define/par P5 0,30 N/A "Offsets in RA,DEC [arcsec] : "
define/par P6 check C/A "Level check after each exposure [check/no_check] : "
define/par P7 30000,45000 N/A "Saturation and monitoring level : "
define/par P8 dusk C/A "Evening or morning twilight [dusk / dawn] : "

if mid_session .ne. 32  then
   write/out "Please use OBSERVING (blue) MIDAS window to start skyflats !"
   $ auplay /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
   goto exit
endif


define/local abort_check/i/1/1 0	! for return value of abort-file-check
	! remove existing abort file
	! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
if abort_check .eq. 1 then
  $rm /disk-a/o2k/tmp/geirsLstAbort 
endif

P6 = m$lower(P6)
P8 = m$lower(P8)

define/local i/i/1/1 0
define/local loop/i/1/1 0
define/local n_exp/i/1/1 {P4}
define/local max_lev/r/1/1 {P2}
define/local curr_lev/r/1/2 0,0
define/local act_lev/r/1/1 0
define/local sat_lev/r/1/2 {P7}
define/local exp_time/r/1/1 0
define/local min_time/r/1/1 0
define/local start_time/r/1/1 60
define/local test_time/r/1/1 0
define/local offset/r/1/2 {P5}
define/local coadds/i/1/1 {P3}
define/local answer/c/1/8 " "
define/local true_time/r/1/1 0.

! setup camera for test exposure

$cmd_o2000 crep 1
$cmd_o2000 object test {P1}
$cmd_o2000 itime 0
$cmd_o2000 itime -stdout | write/key min_time/r/1/1
$cmd_o2000 sync

if P8(1:4) .eq. "dusk"  then   ! ==================== d u s k ================

   testexp_1:

   ! abort check: 0=does not exist ; 1=abort file exists
   abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
   if abort_check .eq. 1 then
      write/out "Program is aborted..."
	$auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
	$rm /disk-a/o2k/tmp/geirsLstAbort 	! remove file again
  	goto exit
   endif

   $cmd_o2000 read
   $cmd_o2000 sync
   $cmd_o2000 median -stdout -raw | write/key curr_lev/r/1/2
   act_lev = curr_lev(2)-curr_lev(1)

   if curr_lev(2) .ge. sat_lev(1) then
      if curr_lev(2) .le. sat_lev(2) then
         write/out "         sky still too bright, continue monitoring ..."
         goto testexp_1	     ! continue monitoring
      else
         write/out
         write/out "         Sky much too bright. Detector saturated ! Abort ..."
         write/out
         $auplay /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
         goto exit
      endif
   else
      set/format I1 F5.0
      write/out "         Minimum exposure time ({min_time}) --> {act_lev} cts/pixel"
      exp_time = max_lev/act_lev*min_time
      set/format F5.1
      write/out "         Sequence requires {exp_time}sec frame integration time."
      if exp_time .gt. 60 then
         $auplay /disk-a/staff/GEIRS/SOUNDS/doorbell.au
         inquire/key answer "Proceed with exposure [y=def/n] ?"
         answer = m$lower(answer)
         if answer(1:1) .eq. "n"  then
            write/out "         Abort. No exposure taken!"
            $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
           goto exit
         else
            write/out "         Taking flats now..."
         endif
      endif
   endif


   if min_time .gt. exp_time  then
      write/out "Too bright, cannot expose shorter than {min_time}"
      goto testexp_1
   endif
 
! take flatfield exposures
   $auplay /disk-a/staff/GEIRS/SOUNDS/whistle.au
   set/format I1 F5.1
   

   do i = 1 {n_exp}
      true_time = m$nint(exp_time*10.)/10.
      $cmd_o2000 object {P1}
      $cmd_o2000 itime {true_time}
      $cmd_o2000 crep {coadds}
      $cmd_o2000 sync
      $cmd_o2000 read
      $cmd_o2000 sync
         ! abort check: 0=does not exist ; 1=abort file exists
      abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
      if abort_check .eq. 1 then
         write/out "Program is aborted..."
         $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
         $rm /disk-a/o2k/tmp/geirsLstAbort 	! remove file again
         goto exit
      endif
      if coadds .eq. 1 then
         $cmd_o2000 save
      else
         $cmd_o2000 save -i
      endif

      if i .lt. n_exp then
         $cmd_o2000 tele relative {offset(1)} {offset(2)}
         $cmd_o2000 sync tele
         if P6(1:2) .ne. "no" then
            $cmd_o2000 object test {P1}
            $cmd_o2000 itime 0
            $cmd_o2000 crep 1
            $cmd_o2000 sync
            $cmd_o2000 read
            $cmd_o2000 sync
            $cmd_o2000 median -stdout -raw | write/key curr_lev/r/1/2
            act_lev = curr_lev(2)-curr_lev(1)
            exp_time = max_lev/act_lev*min_time
            write/out "         Level = {act_lev} --> exposure time for frame {i}+1 = {exp_time}sec"
            if exp_time .gt. 60 then
               $auplay /disk-a/staff/GEIRS/SOUNDS/doorbell.au
               inquire/key answer "Proceed with exposure [y=def/n] ?"
               answer = m$lower(answer)
               if answer(1:1) .eq. "n"  then
                  write/out "         Abort. No exposure taken!"
                  $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
                  goto exit
               else
                  write/out "         Taking next flat ..."
               endif
            endif
         endif
      endif
   enddo
endif

if P8(1:4) .eq. "dawn"  then   ! ==================== d a w n ================

   testexp_2:

   ! abort check: 0=does not exist ; 1=abort file exists
   abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
   if abort_check .eq. 1 then
      write/out "Program is aborted..."
	$auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
	$rm /disk-a/o2k/tmp/geirsLstAbort 	! remove file again
  	goto exit
   endif

   $cmd_o2000 read
   $cmd_o2000 sync
   $cmd_o2000 median -stdout -raw | write/key curr_lev/r/1/2
   act_lev = curr_lev(2)-curr_lev(1)

   if curr_lev(2) .ge. sat_lev(1) then
      write/out
      write/out "         Sky already too bright. Detector saturated ! Abort ..."
      write/out
      $auplay /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
      goto exit
   else
      set/format F5.0
      write/out "         Minimum exposure time ({min_time})--> {act_lev} cts/pixel"
      exp_time = max_lev/act_lev*min_time
      set/format F5.1
      write/out "         Sequence requires {exp_time}sec frame integration time."
      if loop .gt. 0  then
         if exp_time .gt. start_time  then
            goto testexp_2
         else
            write/out "         Starting to take flats ..."
            $auplay /disk-a/staff/GEIRS/SOUNDS/whistle.au
            goto take_data
         endif
      else
         loop = 1
         if exp_time .gt. start_time then
            $auplay /disk-a/staff/GEIRS/SOUNDS/doorbell.au
            write/out
            try_again:
            write/out "         To stop loop type ABORT, else ..."
            inquire/key answer "         Exposure time (sec) to start flats [60]?"
            i = m$len(answer)
            answer = m$upper(answer)
            if i .gt. 0  then
               if answer(1:5) .eq. "ABORT"  then
                  write/out
                  write/out "         Flat sequence aborted by user !"
                  $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
                  goto exit
               else
                  i = m$tstno(answer)
                  if i .eq. 0 then
                     write/out "         Input error, must be number ..."
                     $auplay /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
                     goto try_again
                  endif
                  start_time = {answer}
                  if exp_time .le. start_time  then
                     write/out "         Taking flats now..."
                  else
                     write/out "         Waiting until it is bright enough!"
                     $auplay /disk-a/staff/GEIRS/SOUNDS/whistle.au
                     goto testexp_2
                  endif
               endif
            else
               if exp_time .le. 60  then
                  write/out "         Taking flats now..."
               else
                  start_time = 60
                  write/out "         Waiting until it is bright enough!"
                  $auplay /disk-a/staff/GEIRS/SOUNDS/whistle.au
                  goto testexp_2
               endif
            endif
         endif
      endif
endif

! take flatfield exposures

   take_data:

   $auplay /disk-a/staff/GEIRS/SOUNDS/whistle.au
   set/format i2 F5.1

   

   do i = 1 {n_exp}
      true_time = m$nint(exp_time*10.)/10.
      $cmd_o2000 object {P1}
      $cmd_o2000 itime {true_time}
      $cmd_o2000 crep {coadds}
      $cmd_o2000 sync
      $cmd_o2000 read
      $cmd_o2000 sync
         ! abort check: 0=does not exist ; 1=abort file exists
      abort_check = m$exist("/disk-a/o2k/tmp/geirsLstAbort")
      if abort_check .eq. 1 then
         write/out "Program is aborted..."
         $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
         $rm /disk-a/o2k/tmp/geirsLstAbort 	! remove file again
         goto exit
      endif
      if coadds .eq. 1 then
         $cmd_o2000 save
      else
         $cmd_o2000 save -i
      endif

      if i .lt. n_exp then
         $cmd_o2000 tele relative {offset(1)} {offset(2)}
         $cmd_o2000 sync tele
         if P6(1:2) .ne. "no" then
            $cmd_o2000 object test {P1}
            $cmd_o2000 itime 0
            $cmd_o2000 crep 1
            $cmd_o2000 sync
            $cmd_o2000 read
            $cmd_o2000 sync
            $cmd_o2000 median -stdout -raw | write/key curr_lev/r/1/2
            act_lev = curr_lev(2)-curr_lev(1)
            exp_time = max_lev/act_lev*min_time
            if curr_lev(2) .gt. sat_lev(1) then
               write/out
               write/out "         It is now too bright. Sequence aborted ..."
               $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
               goto exit
            endif
            if exp_time .lt. min_time then
               write/out
               write/out "         It is now too bright. Sequence aborted ..."
               $auplay /disk-a/staff/GEIRS/SOUNDS/crash.au
               goto exit
            else
               write/out "         Level = {act_lev} --> exposure time adjusted to {exp_time}sec for frame {i}+1."
            endif
         endif
      endif
   enddo
endif

write/out  "         All done ..."
$auplay /disk-a/staff/GEIRS/SOUNDS/gong.au

exit:
return
