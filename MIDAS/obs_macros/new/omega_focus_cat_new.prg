!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! .COPYRIGHT    (C) MPIA
! .IDENT    omega_focus.prg
!
! .CALL     @@ omega_focus [root] [start_index] [focus,step] [ima_number] 
!       [itime_total,single] [action_flag] [object_number] [boxsize]
!
! .PARAMETERS   P1=root: root name for images, if omitted filename in GUI is used
!             (for auto-name-finding-->4 digits before .fits are expected)
!       P2=index: index of first image in series --> ima_name = root0001
!       P3=focus [focus,step]: estimated focus value, stepsize   
!                   both in microns
!       P4=number: number of images for focus series
!       P5=itime [total,single]: total integration time, single int_time
!       P6=action: action_flag,graph_output
!              action_flag --> 0=take images, anylize, focus telescope
!                      1=analize only   
!                      2=anyalize and focus telescope
!                      3=take images and analize
!              graph_output--> 0=only final graph is shown
!                      1=show all 4 graphs  
!       P7=catalog: Name of the catalogue used for focussing
!       P8=boxsize: total boxsize in arcsecs around found objects
!
! .PURPOSE  A procedure that determines the best focus for OMEGA2000
! .ENVIRONMENT  MIDAS
! .AUTHOR   Rene Fassbender     
! .KEYWORDS telescope focus
! .COMMENTS PRG has to be started in directory where images are stored
!       This version uses an automatic "find object" routine
!   I   in addition, the selection proccess is automated
!       Three selections will be used: 
!               1) object classification of "find_object"
!               2) check for saturated objects
!               3) galaxy check with an intensity plot
!       The final plot is saved under: masterframe_plot.ps.
!
! .CONVENTIONS  !! --> command of line necessary
!       !* --> command not needed but could be useful for testing
!
! .VERSION  1.0 06/2002     original routine: omega2002.prg
!       2.0 12.12.02    prepare telescope application: .fits imas,
!                   command line parameters, telescope commands
!       2.1 12.01.03    some corrections with root name
!       3.0 13.01.03    parameter adaption for OMEGA2000
!       3.1 15.01.03    changed focus cmds to millimeters units
!                   inroduced pixelscale for O2k
!       3.11    17.01.03    corrected seeing-->missing brackets
!                   and edge selectionfocus
!       3.2 18.01.03    save final plot
!       4.0 12.03.03    implemented check of telescope returns,
!                   graphics output option, abort check
!                   image descriptors
!       4.1     13.03.03    implement automatic filename detection
!       4.2 16.03.03    check before adjusting focus,remove tel_ret
!       4.3 18.03.03    changed default values and add focus check
!       4.31    11.04.03    corrected prompt option and graph default
!       4.4 12.04.03    add time and date in final graph
!       4.5 05.06.03    use environment variable as entry point
!       4.6 06.05.05    change of master possible
!       5.0 09.05.05    focus routine uses catalogues as refernce 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!!!!!!!!!!!!!!!!!!!!! SETUP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cross reference parameter labels
crossref root index focus number itime action caalog boxsize


    ! define parameters and default values
define/parameter P1 + ?         ! root
define/parameter P2 1 C/C       ! start_index
define/parameter P3 ? N/A       ! focus
define/parameter P4 9 N/A       ! imagenumber
define/parameter P5 20,2 N/A        ! integration time
define/parameter P6 0,1 N/C     ! action
define/parameter P7 ? C/A      ! catalog
define/parameter P8 18 N/C      ! boxsize



!!!!!!!!!!!!!!!!!!!!!! INTRODUCTION + HELP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if p1(1:4) .eq. "help" then

write/out
write/out "omega_focus.prg"
write/out "call: @@ omega_focus [root] [start_index] [focus,step] [ima_number]" 
write/out "       [itime_total,single] [action_flag] [object_number] [boxsize]"
write/out 
write/out "The command line parameters are:"
write/out " P1=root: root name for images, if omitted filename in GUI is used"
write/out "         (for auto-name-finding --> 4 digits before .fits are expected)"
write/out " P2=index,master: index of first image in series --> ima_name = root0001"
write/out "                  number of master frame, default is best focus" 
write/out 
write/out " P3=focus [focus,step]: estimated focus value, stepsize "  
write/out "                 both in microns"
write/out " P4=number: number of images for focus series"
write/out " P5=itime [total,single]: total integration time, single int_time"
write/out "                              both in seconds" 
write/out " P6=action: action_flag,graph_output"
write/out "        action_flag --> 0=take images, analyse, focus telescope"
write/out "                1=analyse only"  
write/out "                2=analyse and focus telescope"
write/out "                3=take images and analyse"
write/out "        graph_output--> 0=shows only final graph"
write/out "                1=shows all 4 graphs"
write/out " P7=catalog: name of the catalog used for focussing, e.g. 9h, 18h..."
write/out " P8=boxsize: total boxsize in arcsecs around found objects"
write/out
write/out "This program will find the best telescope focus for OMEGA2000:"
write/out 
write/out "The PRG has to be started in the directory where the images are stored."
write/out
write/out "The image names have to be of the following format: root????.fits,"
write/out "where ???? are consecutive integers, e.g. focus0097, focus0098,..."
write/out
write/out "To find the objects, the images are compared to a catalogue. 
write/out
write/out "Appropriate stellar objects will be selected."
write/out "CENTER/GAUSS will determine the FWHM in the x and y direction "
write/out 
write/out "The final plot is saved under: masterframe_plot.ps" 
write/out
write/out "@@ omega_focus    : take default values for all parameters"
write/out
write/out " -->if root is not specified, filenames from GUI are taken for new images,"
write/out "    or default_root = focus is used analysis with old images "
write/out
write/out

return      ! back to MIDAS session

endif




!!!!!!!!!!!!!!!!!!!!!! DETECTOR PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initialization in section MORE DEFINITIONS
    ! Define a saturation cut-off for the selection
define/local pixelsaturation/r/1/1      ! cut-off = 40,000 counts per ima

    ! Initialization in STELLAR OBJECT SELECTION 4)
    ! Define a FWHM cut-off to select galaxies
define/local galaxycutoff/r/2/1         ! low (1) and high (2) cut-off




!!!!!!!!!!!!!!!!!!!!!!!!  KEYWORDS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fundamental parameters have to be specified by user: # of images, 
! best estimated focus value, stepsize of telescope, name of first image

    ! Define local keywords for the above parameters
define/local imagename1/c/1/60          ! name of first image
define/local imagenumber/i/1/1          ! default is 7 images
define/local focusvalue/i/1/1           ! initial focus value in microns
define/local stepsize/i/1/1             ! stepsize in microns
define/local dum_focus/r/1/2            ! help key for parameter input
define/local int_time/r/1/2             ! total, single integration time
define/local action_flag/i/1/1          ! defines action 
define/local start_index/i/1/1          ! dummy for start index
define/local i_p/i/1/1                  ! test if change of master
define/local reference/i/1/1            ! new master index

    ! keywords for automatic name finding
define/local last_name/c/1/60           ! last image name
define/local slash_pos/i/1/1            ! position of last slash in name
define/local dot_pos/i/1/1              ! position of last period in name
define/local num_pos/i/1/1              ! position of first number in name
define/local last_index/i/1/1           ! key for name extraction
define/local pathname_ima/c/1/200       ! for path and name of first image
define/local autoname_flag/i/1/1 0      ! flag for automatic name finding
                                        ! 0=no autoname ; 1=automatic name

define/local makeodd/r/1/1              ! dummy to make imagenumber odd
define/local master_index/i/1/1         ! index of master frame
define/local loop/i/1/1                 ! loop parameter
define/local root_save/c/1/80           ! dummy for root name

define/local graph_aux/i/1/2 0,1        ! dummy for graphics output
define/local graph_output/i/1/1 1       ! keyword for graphics flag
define/local foc_time/c/1/30            ! for storing date and time
define/local tel_return/i/1/1 0           ! telescope status


    ! Define the master frame
    ! Define the number of objects for find_object
    ! Define the boxsize for the region of interest around the selected objects
define/local masterframe/c/1/60     
define/local catalog/c/1/60     ! default = 40 objects  
define/local boxsize/r/1/1      ! default value = 18 (world coordinates)
define/local halfbox/r/1/1      ! half the boxsize

    ! Define the pixelscale for omega2000
define/local pixelscale/r/1/1 0.45



!!!!!!!!!!!!!!!!!!!!! INITIALIZE KEYWORDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

set/format

write/key imagenumber {P4}

    ! Make sure the number of images is odd and >= 3; 
    ! otherwise set to 7 or made odd
if imagenumber .lt. 3 then
imagenumber = 7
write/out "Imagenumber was set to 7, since input was too small... "
endif

makeodd = (imagenumber-0.5)/2
imagenumber = 2*m$nint(makeodd)+1   ! transfer even # to next odd number
                    ! m$int() gives back nearest interger



    ! get master 
i_p = m$index(P2,",")
if i_p .eq. 0 then
   write/keyword start_index {P2}
    reference =  0
else
   i_p = i_p - 1
   write/keyword start_index {P2(1:{i_p})}
   i_p = i_p + 2
   reference =  {P4({i_p}:)}
endif

write/keyword dum_focus/r/1/2 {P3}
write/keyword focusvalue {dum_focus(1)}
write/keyword stepsize {dum_focus(2)}

write/keyword int_time {P5}
write/keyword action_flag {P6}
write/keyword graph_aux/i/1/2 {P6}
write/keyword catalog {P7}
write/keyword boxsize {P8}

!if m$exist("/disk-a/o2k/focus/{catalog}.tbl") .eq. 0 then
if m$exist("{catalog}.tbl") .eq. 0 then
    write/out "Catalog not found, please check location of {catalog}.tbl"
!    write/out "It should be localized in /disk-a/o2k/focus/"
!    write/out "We highly recommend to stick to the designated focus fields!"
!    write/out "However, there is the possibility to create your own focus field:"
!    write/out "@@ O2K_UTIL:/obs_macros/crea_foc_cat"
    $play -q /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
    goto exit
endif    

crea/tab focus_cat 2 500

!Copy/table /disk-a/o2k/focus/{catalog} focus_cat
Copy/table {catalog} focus_cat

    ! initialize graphics output flag
graph_output = graph_aux(2)

    ! compute index for master image
if reference .eq. 0 then
    master_index = m$nint(makeodd)+{start_index}
else
    master_index = {reference}
endif

    ! do the root name last
if P1(1:1) .ne. "+" then        ! real root name

write/keyword imagename1 {P1}{start_index}.fits
write/keyword masterframe {P1}{master_index}.fits

else                    ! if "go", use default value

write/keyword P1 focus
write/keyword imagename1 focus0001.fits
write/keyword masterframe focus{master_index}.fits
  
write/keyword autoname_flag 1       ! set autoname flag
  
endif


    ! initialize time and date
foc_time = m$time()








    ! check telescope parameters for valid input
if focusvalue .lt. 15000 .or. focusvalue .gt. 40000 then
  write/out
  write/out "The estimated telescope focus of {focusvalue} is outside the valid range"
  write/out "focus limits of [15000,40000]! The program is aborted..."
  write/out
  $play -q /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
  return
endif


if stepsize .lt. 10 .or. stepsize .gt. 1000 then
  write/out
  write/out "The focus steps of {stepsize} are outside the limits of [10,1000]"
  write/out "The stepsize is set to the default value of 200..."
  stepsize = 200
  $play -q /disk-a/staff/GEIRS/SOUNDS/sorrydave.au
endif





!!!!!!!!!!!!!!!!!!!!!!!! DISPLAY PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

set/format I1 F5.2

write/out
write/out "The following parameters are used:"
write/out "Automatic name finding enabled (0=no;1=yes):{autoname_flag}"
write/out "Estimated focus value in microns: {focusvalue}"
write/out "Stepsize in microns: {stepsize}"
write/out "Number of images: {imagenumber}"
write/out "Total integration time per image in seconds: {int_time(1)}"
write/out "Integration time for single frame in seconds: {int_time(2)}"
write/out "--> the following three items are not active if set to auto-name-finding "
write/out "Root name: {P1}"
write/out "Name of first image: {imagename1}"
write/out "Name of masterframe: {masterframe}"
write/out "Action_flag is set to: {action_flag}"
write/out "(0=take_ima+analyze+move_focus; 1=anal.only; 2=anal.+focus; 3=take+anal)"
write/out "Graph_output is set to: {graph_output}"
write/out "(0=only final graph is shown ; 1=all 4 graphs are shown)"
write/out "Used catalogue: {P7}"
write/out "Total boxsize in arcsecs around objects: {boxsize}"
write/out





!!!!!!!!!!!!!!!!!!!!!!!! TAKE IMAGES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

define/local repetition/i/1/1       ! # of repetitions for a single sum image
repetition = int_time(1)/int_time(2)+0.5

define/local abort_check/i/1/1 0    ! for return value of abort-file-check


!++++++++++++++++++++++++ START action_flag=0 +++++++++++++++++++++++++++++++++++++
    ! start data taking if action flag is set accordingly
if action_flag .eq. 0 .or. action_flag .eq. 3 then

!move telescope to exact position needed by catalogue!!

DEFINE/LOC ra/C/1/60 " "
DEFINE/LOC dec/C/1/60 " "
COPY/DKEY /disk-a/o2k/focus/{P7}.tbl RA ra
COPY/DKEY /disk-a/o2k/focus/{P7}.tbl DEC dec
!   $ {tecs_script}/caget -t mpia_t3_delta | write/keyword deltazero(1)Maybe this one will come in handy...
!    $ {tecs_script}/caget -t mpia_t3_alpha | write/keyword alphazero(1)

$ {tecs_script}/t_posit {ra} {dec} 2000.0
! müsste das nicht die aktuelle Zeit sein??
write/out "Telescope is moved to coordinates of the catalogue {P7}: {ra} , {dec} "

write/out "{imagenumber} images are now taken at the current observing position..."
write/out   
    

    ! remove existing abort file
    ! abort check: 0=does not exist ; 1=abort file exists
abort_check = m$exist("{geirslstabort}")
if abort_check .eq. 1 then
  $rm {geirslstabort} 
endif


    ! define keywords for the data taking
define/local first_focus/r/1/1      ! first absolute focus position 
first_focus = focusvalue - (stepsize * m$nint(makeodd))     ! in microns
first_focus = first_focus/1000      ! in millimeters as needed for t_afocus


define/local real_step/r/1/1 {stepsize} ! real keyword for stepsize
real_step = real_step/1000      ! step in millimeters

define/local current_name/c/1/60    ! key for current image name
current_name = imagename1

define/local current_index/i/1/1
write/keyword current_index {start_index}


    ! set telescope to first absolute focus position
set/format f6.3
write/out
write/out "First image taken at focus position: {first_focus} mm"


$cmd_panic_new sync         ! wait until telescope processes are finished

    ! first focus position
$ {tecs_script}/t_afocus {first_focus} -
        | awk '{if(NR==1){print $1}}' | write/keyword tel_return    ! pipe return
      if tel_return .lt. 0 then
         write/out "ERROR: Telescope return value for t_afocus signals an error..."
         write/out "...the program is aborted"
         $play -q /disk-a/staff/GEIRS/SOUNDS/crash.au
         goto exit
      endif 

wait/sec 10.    ' to ensure telescope focus is really set

    ! set the repetition parameter and integration time
set/format I1
$cmd_panic_new crep {repetition}
$cmd_panic_new itime {int_time(2)}
$cmd_panic_new sync

    ! initialize image descriptors
$cmd_panic_new counter DITH_NO clear        ! set dither counter to 0
$cmd_panic_new counter POINT_NO clear       ! set pointing no to 0
$cmd_panic_new counter EXPO_NO clear        ! reset exposure counter  



!-----------------------------------------------------------
    ! do loop over number of images
do loop = 1 {imagenumber} 1

    set/format f6.3

    ! telescope macro
    $cmd_panic_new object focus: {loop} of {imagenumber}
    
    $cmd_panic_new read
    $cmd_panic_new sync

        ! abort check: 0=does not exist ; 1=abort file exists
    abort_check = m$exist("{geirslstabort}")
    if abort_check .eq. 1 then
      write/out "Program is aborted..." 
      $rm {geirslstabort}     ! remove file again
      $play -q /disk-a/staff/GEIRS/SOUNDS/crash.au
      return
    endif

    if autoname_flag .eq. 0 then    
      $cmd_panic_new save -i {current_name} ! file with specified name is saved
    else    
      $cmd_panic_new save -i            ! file with name in GUI is saved
    endif
    
    
    set/format f5.3

    $ {tecs_script}/t_dfocus {real_step} -
        | awk '{if(NR==1){print $1}}' | write/keyword tel_return    ! pipe return
      if tel_return .lt. 0 then
         write/out "ERROR: Telescope return value for t_dfocus signals an error..."
         write/out "...the program is aborted"
         $play -q /disk-a/staff/GEIRS/SOUNDS/crash.au
         goto exit
      endif 

    
    if autoname_flag .eq. 0 then
      write/out "Image {current_name} was saved..."
    endif


    if loop .lt. {imagenumber} then
    write/out "         {loop} image done. Focus was increased by {real_step} mm for next image..."
    endif

    ! get name for next image
    set/format I4
    current_index = current_index+1     ! increment index
    write/keyword current_name {P1}{current_index}.fits

enddo
!------------------------------------------------------------


write/keyword autoname_flag 2           ! 2-->autoname was used 

$cmd_panic_new sync

write/out
write/out "Taking of focus frames done..."
write/out
$play -q /disk-a/staff/GEIRS/SOUNDS/whistle.au


endif   ! end of image taking 
!++++++++++++++++++++++++ END action_flag=0 +++++++++++++++++++++++++++++++++++++




!!!!!!!!!!!!!!!!!!!!!!!!! MORE DEFINITIONS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! define saturation cutoff = number of integrated images * saturation for single image
    ! 40,000 is saturation limit per single image
    ! repetition is the number of integrated images
pixelsaturation = repetition * 30000



define/local halfnumber/i/1/1       ! # of points on each side of the
halfnumber = (imagenumber-1)/2      ! best focus: 7 images --> 3
halfbox = boxsize/(2*{pixelscale})  ! for o2k in arcsecs
boxsize = boxsize/pixelscale        ! for o2k in arcsecs

    ! save the root name before it gets lost
write/keyword root_save/c/1/80 {P1}


! Create a table "storeresults" for saving the average fwhm and the stddev
! in the x and y direction of each image
! Format: #rows= #images ; #columns= 5
! Column names: :xavfwhm, :yavfwhm, :xstddev, ::ystddev
! :xscale holds the matching x-values for the final plot


create/table storeresults 5 {imagenumber} null
    ! there is a bug that sets keyword P1 back to the default value
    
!*write/out "P1 in 294 ist {P1}"

create/column storeresults :xavfwhm
create/column storeresults :yavfwhm
create/column storeresults :xstddev
create/column storeresults :ystddev
create/column storeresults :xscale
!*read/table storeresults           ! Test


    ! reinitialize the root name
write/key P1/c/1/80 {root_save}
    
!*write/out "P1 in 298 ist {P1}"



!!!!!!!!!!!!!!!!!!! EXTRACT AUTOMATIC NAME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! if true, automatic-name-finding was used
    ! the names for the analysis have to be extracted now
if autoname_flag .eq. 2 then

  $cmd_panic_new last   ! writes last filename in file geirsLstFile

    ! writes last filename in keyword pathname_ima
  write/keyword pathname_ima </disk-a/o2k/tmp/geirsLstFile  

    ! get position of last / in patname_ima
  slash_pos = m$indexb(pathname_ima,"/")
  slash_pos = slash_pos+1           ! position of first root character

    ! store last name in keyword
  last_name = pathname_ima({slash_pos}:)
  write/out "Extracted last name is: {last_name}"

  dot_pos = m$indexb(last_name,".") ! position of last period
  num_pos = dot_pos-4           ! position of first digit
  dot_pos = dot_pos-1           ! now position of last digit

  last_index = {last_name({num_pos}:{dot_pos})}  ! index of last image

  num_pos = num_pos-1           ! position of last root character

  P1 = last_name(1:{num_pos})       ! root

    ! build first image name and mastername
  start_index = last_index-imagenumber+1
  if reference .eq. 0 then 
    master_index = last_index-((imagenumber-1)/2)
  else 
    master_index = last_index - imagenumber + 1 + reference - start_index
  endif
  
  set/format I4

  imagename1 = "{P1}{start_index}.fits"
  masterframe = "{P1}{master_index}.fits"

  write/out "The extracted names for first and master image are:"
  write/out "{imagename1} and {masterframe}"

endif


!*return    ! test




!!!!!!!!!!!!!!!!!!!  COMPARISON WITH CATALOGUE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CREA/TAB gausstable 5 500

compute/table focus_cat :xstart = :X-{halfbox}
compute/table focus_cat :xend = :X+{halfbox}
compute/table focus_cat :ystart = :Y-{halfbox}
compute/table focus_cat :yend = :Y+{halfbox}

CENTER/GAUSS {masterframe},focus_cat gausstable
SELECT/TAB gausstable :XCEN.NE.NULL

set/midas output=logonly
Stat/tab gausstable :XCEN
set/midas output=yes 
 
if outputi(2) .lt. 3 then
    WRITE/OUT "There are not enough stars found to calculate focus, something is wrong"
    $play -q /disk-a/staff/GEIRS/SOUNDS/crash.au
    !goto exit 
endif

clear/channel overlay                   ! clears the cursor marks

    ! mark the found objects on the display
load/table gausstable :xcen :ycen ? 1 3 3 -1
    ! x_pos,y_pos, circles(1), size(3), red(3), do not connect(-1)
load/tab focus_cat :X :Y ? 1 5 blue -1    
    

!!!!!!!!!!!!!!!!!!!! HISTOGRAM OBJECTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! convert to arcsecs
compute/table gausstable :xfwhm = :xfwhm*{pixelscale}
compute/table gausstable :yfwhm = :yfwhm*{pixelscale}


    ! Plot XFWHM and YFWHM in a histogramm, if graphics output is wanted
if graph_output .eq. 1 then

  create/graphics           ! open a graphics window    

  set/graphics colour=1         ! black lines for xfwhm plot

  plot/histogram gausstable :xfwhm  ? 0.1,0,5 ? ? 1,1,45
        ! automatic scaling ; binsize 0.1 , min=0, max=5 arcsec
        ! style: staircase steps with hashing

  label/graphic "X-,Y-FWHM-Histogram" 0.6,0.9,no 0 1 0
  label/graphic "of all objects" 0.6,0.85,no 0 1 0
    ! Label: position in no=normalized coord, 0 angle, 1=size, centerd
        
  set/graphics colour=2         ! red lines for overplot
  overplot/histogram gausstable :yfwhm ? 0.1 ? ? 1,1,-60
                    ! bin_size 0.1
endif

    
    ! Make a new column :averagefwhm that holds the average of the two FWHM
compute/table gausstable :averagefwhm = sqrt(:xfwhm*:yfwhm)


!!!!!!!!!!!!!!!!!!  HISTOGRAM: SELECTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Plot a FWHM histogramm of the selected objects

if graph_output .eq. 1 then 

  create/graphics 2         ! new graphics window for selected objects
                    ! repeat plotting procedure
  set/graphics colour=1         ! black lines for xfwhm plot
  plot/histogram gausstable :averagefwhm  ? 0.1,0,5 ? ? 1,1,45
        ! automatic scaling ; binsize 0.1 , min=0, max=5 arcsec
        ! style: staircase steps with hashing
        
  label/graphic "Averaged FWHM" 0.6,0.9,no 0 1 0
  label/graphic "of selected objects" 0.6,0.85,no 0 1 0
    ! Label: position in no=normalized coord, 0 angle, 1=size, centerd

endif



!!!!!!!!!!!!!!!!!!  APPLY TO ALL FRAMES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The selected regions of the master frame will now be applied to all frames.
! Center/Gauss will calculate the fwhm in every frame.
! The average of the FWHM  in the x and y direction of every frame is 
! calculated

write/out
write/out "The selected regions of the master frame are now applied to"
write/out "all frames. The average of the fwhm in the x and y direction"
write/out "is calculated and then displayed..."
write/out
write/out "...the results for all frames starting with first are:"
write/out
write/out
    

set/format

    ! A DO loop will apply center/gauss to all images and store the average
    ! x,y-fwhm value and the stdev in the table "storeresults"

define/local count/i/1/1 1      ! define a loop variable
define/local name_length/i/1/1      ! holds # of characters of frame name
define/local imagename/c/1/60 {imagename1}      
    ! keyword that holds current image name, 
    !initialized with name of first image

name_length = m$strlen(imagename)
name_length = name_length-8 ! go to beginning of the number index

define/local runnumber/i/1/1 {start_index}  ! first index
    ! framenumber of focus0097.fits



!----------------------------------------------------------------------------------

do count = 1 {imagenumber} 1        ! loopvar=start, end, step

!*read/keywort imagename            ! check imagename

set/midas output=logonly

center/gauss {imagename},gausstable frameresults
    ! use columns x,ystart & x,yend of find_object_results 
    ! on the image focusxxx
    ! and store results in the temporary table frameresults

    ! convert to arcsecs
compute/table frameresults :xfwhm = :xfwhm*{pixelscale}
compute/table frameresults :yfwhm = :yfwhm*{pixelscale}


!*set/midas output=yes      ! if active, results are displayed on terminal
!*read/table frameresults :xfwhm,yfwhm  ! Check output


    ! find the average of each column
    ! first for :xfwhm

statistics/table frameresults :xfwhm 
                ! calculates: min, max, mean(3), stddev(4)
                ! of column xfwhm; stored in  real keyword
                ! OUTPUTR

    ! average and stddev is stored
copy/kt outputr/r/3/2 storeresults :xavfwhm :xstddev @{count}


    ! now for :yfwhm
statistics/table frameresults :yfwhm 
    
    ! average and stddev is stored
copy/kt outputr/r/3/2 storeresults :yavfwhm :ystddev @{count}

    ! Imagename of next image has to be written in keyword "imagename"
set/format I4               ! Integer has to be: 0001
runnumber = {runnumber}+1
!*read/keyword runnumber
write/keyword imagename/c/{name_length}/4 {runnumber}

enddo                   ! end of do loop

!----------------------------------------------------------------------------------


set/midas output=yes            ! turn terminal log back on
read/table storeresults         ! check table entries
set/format              ! back to default format




!!!!!!!!!!!!!!!!!!!!!  PLOT RESULTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The averaged x,y-fwhm values for each frame are plotted
! REGRESSION/POLY is used to fit a parabola to the values
! The minimum of the parabola determines the best focus position (x-axis)
! and the seeing (y-axis)
! The seeing is the arithmetic mean of the two minima ( x & y-fwhm)
! The best focus position is the geometric mean of the two minima

write/out
write/out
write/out "The results are now plotted. Two parabolas are fitted and from"
write/out "the minima the best focus value and the seeing are determined..."
write/out
write/out


    ! First determine appropriate x and y scale

    ! Define keywords for the scaling

define/local xmin/r/1/1 -{halfnumber}
define/local xmax/r/1/1 {halfnumber}
define/local ymin/r/1/1 0
define/local ymax/r/1/1 5
define/local xoffset/r/1/1 1            ! x-offset for graph
define/local yoffset/r/1/1 0.5          ! y-offset for graph


set/midas output=logonly

statistics/table storeresults :xavfwhm      ! get min & max y-values

ymin = outputr(1)
ymax = outputr(2)

statistics/table storeresults :yavfwhm      ! get min & max y-values

set/midas output=yes

if outputr(1) .lt. ymin ymin = outputr(1)   ! choose lowest min
if outputr(2) .gt. ymax ymax = outputr(2)   ! choose highest max

ymin = ymin-yoffset             ! include an offset
ymax = ymax+yoffset
xmin = xmin-xoffset
xmax = xmax+xoffset             ! final scale parameters 


    ! set up an appropriate coordinate system 
create/graphics 3               ! use graphics window 3
set/graphics colour=1               ! back to black lines

plot/axes {xmin},{xmax},100,100 {ymin},{ymax} ?-
 "Focus position (micron)" "FWHM (arcsec)" 
        ! 100=distance between big tickmarks; 100= distance between small
        ! --> large distance, so tickmarks do not show up
        ! Will be changed to different units later
        
        
    
    ! Plot the data points

    ! First for XFWHM
set/graphics ltype=0 stype=4            ! no line ; triangles(4) as symbols

compute/table storeresults :xscale = sequence-1-{halfnumber}    
    ! vector for the x-axis input, symmetric around 0       


!*read/table storeresults 

overplot/table storeresults :xscale :xavfwhm
        ! x-scale is in arbitrary units e.g. -3,..,0,...,3
        
    
    ! Now for YFWHM
set/graphics ltype=0 colour=2 stype=3       ! red squares ,no line
overplot/table storeresults :xscale :yavfwhm

!!!!!!!!!!!!!!!!!!!!! remove most discrepant points !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

compute/tab storeresults :diff_fwhm = :xavfwhm - :yavfwhm
set/midas output=logonly

stat/tab storeresults :diff_fwhm
define/local var_fwhm/r/1/1 0
var_fwhm = 2.5*outputr(4)
sel/tab storeresults abs(:diff_fwhm).le.{var_fwhm}
set/midas output=yes

!!!!!!!!!!!!!!!!!!!!!  FIT AND PLOT PARABOLAS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fit two parabolas to the data points

! Create a table "coefficients" with columns :xcoeff and :ycoeff where the
! coefficients of the fitted parabolas are stored
! The parabola has the form: y = a + b*x + c*x^2

create/table coefficients 2 3 null
create/column coefficients :xcoeff
create/column coefficients :ycoeff

set/midas output=logonly

    ! Parabola for :XFWHM
regression/polynominal storeresults :xavfwhm :xscale 2
                ! table, dependent var., independ. var, degree

    ! Save the coefficients
set/format
coefficients,:xcoeff,@1 = {outputd(1)}
coefficients,:xcoeff,@2 = {outputd(2)}
coefficients,:xcoeff,@3 = {outputd(3)}


    ! Parabola for :YFWHM
regression/polynominal storeresults :yavfwhm :xscale 2
                ! table, dependent var., independ. var, degree

    ! Save the coefficients
coefficients,:ycoeff,@1 = {outputd(1)}
coefficients,:ycoeff,@2 = {outputd(2)}
coefficients,:ycoeff,@3 = {outputd(3)}

set/midas output=yes

!*read/table coefficients


    ! Plot the fit-parabolas

    ! Create a table "fitpoints" with columns :absc, :xfit, :yfit
    ! that hold the x and y-axis points to be plotted
create/table fitpoints 3 200        ! 3 columns, 200 points each
create/column fitpoints :absc       ! absc = abscissa
create/column fitpoints :xfit
create/column fitpoints :yfit

    ! Get evenly spread points over the abscissa interval (xmin, xmax)
compute/table fitpoints :absc = {xmin}+(sequence*2*{xmax}/200)
    
    ! Compute the corresponding y-values for the XFWHM-parabola 
compute/table fitpoints :xfit = {coefficients,:xcoeff,1}+-
{coefficients,:xcoeff,2}*:absc+{coefficients,:xcoeff,3}*(:absc**2)
    
    ! And for the YFWHM-parabola
compute/table fitpoints :yfit = {coefficients,:ycoeff,1}+-
{coefficients,:ycoeff,2}*:absc+{coefficients,:ycoeff,3}*(:absc**2)  

!*read/table fitpoints


    ! Plotting

set/graphics ltype=1 colour=1 stype=0   ! solid black line, no symbols
overplot/table fitpoints :absc :xfit    ! XFWHM-parabola

set/graphics colour=2           ! solid red lines
overplot/table fitpoints :absc :yfit    ! YFWHM-parabola


    ! Relabel the x-axis with the real telescope 
    ! position in microns (instead of image number)
xmin = focusvalue-(m$abs(xmin)*stepsize)    ! m$abs() = absolute value
xmax = focusvalue+(xmax*stepsize)
!*read/keyword xmin,xmax

set/graphics colour=1           ! back to black
define/local smallticks/r/1/1       ! distance between small tickmarks
smallticks = stepsize/5         ! stepsize = dist between large tm
                        
overplot/axes {xmin},{xmax},{stepsize},{smallticks} {ymin},{ymax} ?-
 "Focus position (micron)" "FWHM (arcsec)" 





!!!!!!!!!!!!!!!!!  BEST FOCUS AND SEEING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Calculate the best focus position from the minimum of the parabola

    ! Keyword definitions
define/local seeingxy/r/2/1         ! y_min of the two parabolas
define/local seeing/r/1/1
define/local bestfocusxy/r/2/1          ! minimum of the two parabolas
define/local bestfocus/r/1/1            ! final result


    ! The minimum of the parabola y = a + b*x +c*x^2 is:
    ! x_min = -b / 2c
bestfocusxy(1) = -({coefficients,:xcoeff,2}/(2*{coefficients,:xcoeff,3}))
bestfocusxy(2) = -({coefficients,:ycoeff,2}/(2*{coefficients,:ycoeff,3}))


    ! Best focus position is arithmetic mean of the two minima
bestfocus = (bestfocusxy(1)+bestfocusxy(2))/2   ! still arbitrary units
bestfocus = (bestfocus*stepsize)+focusvalue     ! in microns 
bestfocusxy(1) = (bestfocusxy(1)*stepsize)+focusvalue 
bestfocusxy(2) = (bestfocusxy(2)*stepsize)+focusvalue 


    ! The y-value @ the minimum of a parabola is:
    ! y_min = a - b^2/4c
seeingxy(1) = {coefficients,:xcoeff,1}-(({coefficients,:xcoeff,2})**2/(4*{coefficients,:xcoeff,3}))
seeingxy(2) = {coefficients,:ycoeff,1}-(({coefficients,:ycoeff,2})**2/(4*{coefficients,:ycoeff,3}))

    ! The real seeing is the geometric mean of the two y_min values
seeing = m$sqrt(seeingxy(1)*seeingxy(2))



! Display the final results on the screen and in graphics window

set/format f6.3                 ! give out appropriate format

label/graphic "Seeing = {seeing} arcsec" 60,115,mm 0 1.5 1
    ! x, y_position in mm ; angle size pos_ind=starting @ position


    ! write time and date in graph
if action_flag .eq. 0 .or. action_flag .eq. 3 then  ! images were just taken

  label/graphic "Focus series was taken:" 60,102,mm 0 0.8 1
  label/graphic "{foc_time}" 115,102,mm 0 0.8 1

else

  label/graphic "Focus analysis was done:" 60,102,mm 0 0.8 1
  label/graphic "{foc_time}" 115,102,mm 0 0.8 1

endif

    
write/out
write/out
write/out "=================================================================="
write/out "The seeing is: {seeing} arcsecs (FWHM)"
set/format f6.1
write/out
write/out "The focus is at telescope position: {bestfocus} microns"
write/out "(in X = {bestfocusxy(1)}; in Y = {bestfocusxy(2)})"
write/out "=================================================================="
write/out

label/graphic "Focus = {bestfocus}" 60,108,mm 0 1.5 1

set/format                  ! back to default


    ! Label symbols

overplot/symbol 4 35,30,mm
label/graphic "x-fwhm" 40,30,mm 0 1 1
label/graphics "{masterframe}" 33.75,20,mm 0 0.8 1

set/graphics colour=2               ! red

overplot/symbol 3 35,25,mm
label/graphic "y-fwhm" 40,25,mm 0 1 1


    ! save the final plot
write/out
copy/graphics postscript
$mv postscript.ps {masterframe}_plot.ps




!!!!!!!!!!!!!!!!!!!! FOKUS TELESCOPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! adjust telescope focus position if action_flag is set
if action_flag .eq. 0 .or. action_flag .eq. 2 then

  set/format f6.1   

    ! do check before moving the focus
  define/local move_check/c/1/3 "no"

  write/out
  write/out

  $play -q /disk-a/staff/GEIRS/SOUNDS/doorbell.au
  inquire/keyword move_check -
        "Adjust telescope focus to {bestfocus} (y/n)? [{move_check}]:" flush
  write/out


  if move_check(1:1) .ne. "y" then  ! if not yes --> return
    write/out "Telescope focus was not adjusted..."
    write/out
    return
  endif


  write/out "The telescope focus is now adjusted to the best focus value: {bestfocus}"
  write/out

  set/format f6.3 

  bestfocus = bestfocus/1000    ! in millimeters

  $ {tecs_script}/t_afocus {bestfocus} -
        | awk '{if(NR==1){print $1}}' | write/keyword tel_return    ! pipe return
      if tel_return .lt. 0 then
         write/out "ERROR: Telescope return value for t_afocus signals an error..."
         write/out "... could not set focus to best value!"
         $play -q /disk-a/staff/GEIRS/SOUNDS/crash.au
         goto exit
      endif 

  set/format

  write/out "Focusing done..."
  write/out
exit:
endif
$play -q /disk-a/staff/GEIRS/SOUNDS/gong.au

return
