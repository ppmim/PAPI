! OMEGA2000 observation plan           P.I.     date

! The following is an example of how a P.I. might specify 
! his observations with OMEGA2000 to be executed in service mode.
!
! Please do not change neither the prolog nor the epilog!
! Please enter a unique name for this OB in the line starting define/local OB_name ...
!
! For each task (between the ++++++++ lines) you specifiy the following:
!
!        Comment line(s), summarizing the task (starting with exclamation mark)
!        Label with task number (task_ followed by three digits),  e.g.  task_004:
!        Line to write entry for this task into the log-file ({log_entry})
!        What you want to do (either $cmd_o2000 commands or o2k/...,
!                    or MIDAS manipulations like loading an image into the display etc.
!        Exit the procedure (goto exit)
!
! It is assumed that you are familiar with the latest edition of the OMEGA2000 manual !
!

!--------------------------------------------------------------------------
define/par P1 ? ? "OB_name : "                         !
define/par P2 ? ? "Specify task number : "             !
                                                       !
define/local task/c/1/16 "task_"                       !  p r o l o g
define/local task_i/i/1/1 {P2}                         !
define/local last_ima/c/1/80 " "                       !
define/local answer/c/1/16 " "                         !     d o
define/local OB_name/c/1/80 {P1}                       !
define/local i/I/1/1 0                                 !
define/local fctrl/i/1/2 0,0                           !    n o t 
define/local isodate/c/1/32 " "                        !
isodate = m$isodate()                                  !
open/file {P1}.log WRITE fctrl                         !  c h a n g e 
write/file {fctrl(1)} {isodate} Log for OB {P1}        !
                                                       !
set/format I3                                          !    t h e s e
                                                       !
write/key task/c/6/3 {task_i}                          !
define/local log_entry/c/1/80 "write/file {fctrl(1)} - !    l i n e s !!
                  {isodate} {OB_name} {task} started"  !
                                                       !
goto {task}                                            !
                                                       !
!-------------------------------------------------------------------------

! Dark exposures

task_001:
{log_entry}

o2k/dark 200 5				! 5 dark exposures up to 200sec, log spacing

goto exit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Dome flats

task_002:
{log_entry}

$cmd_o2000  filter H			! select H filter
o2k/domeflats 100000 1 "dome flat" 5,0	! take individual flats, sum = 100000 with lamp 5

goto exit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! sky flats

task_003:
{log_entry}

$cmd_o2000 tele absolute 12 13 14 -0 3 15 2000	! move telescope to flat field position

o2k/skyflats  H   20000 3   2  "sky flat H"	! take 3 flats (each sum of 2 exposures)

goto exit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! take science exposures

task_004:
{log_entry}

o2k/aequi "12 13 14.0 -0 54 32 2000" H flat_H 1234,1456    ! acquire field
if Q1(1:4) .eq. "stop"  then
   write/out "Acquisition aborted"
   write/out "No data taken"
   goto exit
endif

o2k/dither  my_catalogue_1   1500,60,3  "my object"        ! take first series
o2k/dither  my_catalogue_1   1500,60,3  "my object" 26     ! continue series

goto exit
!---------------------------------------------------------------------------------
                                                       !
exit:                                                  !
                                                       !        e p i l o g
close/file {fctrl}   ! close log-file                  !
                                                       !  d o  n o t  c h a n g e
return                                                 !
!---------------------------------------------------------------------------------
