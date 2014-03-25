!
! OMEGA_tilt.prg
!
! A procedure that determines the detector tilt of OMEGA2000
! based on the focus-routine omega2002.prg	
!
! Rene Fassbender     12/2002
!
!
! This version uses an automatic "find object" routine
! In addition, the selection proccess is automated
! Three selections will be used: 1) object classification of "find_object"
!				 2) check for saturated objects
!				 3) galaxy check with an intensity plot



! Internal conventions: !! --> command of line necessary
!			!* --> command not needed but could be useful for testing

! 05.06.03	use environment variable O2K_UTIL for findobj path definition


!!!!!!!!!!!!!!!!!!!!!! INTRODUCTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write/out
write/out "This program will find the detector tilt for OMEGA2000:"
write/out 
write/out "The image names have to be of the following format: focus???.prg,"
write/out "where ??? are consecutive integers, e.g. focus097, focus098,..."
write/out
write/out "Objects on the master frame are detected with ´find_object´ "
write/out "and marked on the display"
write/out
write/out "Appropriate stellar objects will be selected."
write/out "CENTER/GAUSS will determine the FWHM in the x and y direction "
write/out 
write/out 
write/out "First, some important parameters have to be specified:"
write/out
write/out



!!!!!!!!!!!!!!!!!!!!!!!!  KEYWORDS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fundamental parameters have to be specified by user: # of images, 
! best estimated focus value, stepsize of telescope, name of first image

define/local count/i/1/1 1		! define a loop variable

	! Define local keywords for the above parameters
define/local imagename1/c/1/20 focus097	! name of first image
define/local imagenumber/i/1/1 7	! default is 7 images
define/local focusvalue/i/1/1 30000	! best initial focus value in microns
define/local stepsize/i/1/1 200		! stepsize in microns

define/local makeodd/r/1/1		! dummy to make imagenumber odd


	! Define the master frame
	! Define the number of objects for find_object
	! Define the boxsize for the region of interest around the selected objects
define/local masterframe/c/1/20 focus100 	
define/local objectnumber/i/1/1 200	! default = 40 objects  
define/local boxsize/r/1/1 32		! default value = 18 (pixel)
define/local halfbox/r/1/1 16		! half the boxsize

	! Define a saturation cut-off for the selection
define/local pixelsaturation/i/1/1 200000	! cut-off = 350000 counts / pix	
	! Define a FWHM cut-off to select galaxies
define/local galaxycutoff/r/2/1			! low (1) and high (2) cut-off


!!!!!!!!!!!!!!!!!!!!!!  INQUIRE KEYS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! A shortcut to skip inquiries
!*goto noinquire !**********


! Ask user for input of parameters

set/format I1 F3.0			! set format to appropriate output

inquire/keyword focusvalue "Enter your estimated focus value in microns-
 [{focusvalue}]:"

inquire/keyword imagenumber "Enter number of images you want to take-
 (5, 7, 9, ...)  [{imagenumber}]:"

inquire/keyword stepsize "Enter the stepsize in microns [{stepsize}]:"

inquire/keyword imagename1 "Enter the name of the first image [{imagename1}]:"

inquire/keyword masterframe "Enter the imagename of the master frame-
 [{masterframe}]:"


inquire/keyword objectnumber "Enter the number of objects you want to use-
 [{objectnumber}]:"
 
inquire/keyword boxsize "Enter the boxsize for the region of interest around-
 the objects (in arcsecs) [{boxsize}]:"	

set/format				! back to default format

!*noinquire: !**************


	! Make sure the number of images is odd and >= 3; 
	! otherwise set to 7 or made odd
if imagenumber .le. 3 imagenumber = 7
makeodd = (imagenumber-0.5)/2
imagenumber = 2*m$nint(makeodd)+1	! transfer even # to next odd number
					! m$int() gives back nearest interger



!!!!!!!!!!!!!!!!!!!!!!!! DEFINITIONS FOR TILT !!!!!!!!!!!!!!!!!!!!!!!!!!
	! Create image for results with name: tilt_frame 
create/image tilt_frame = {masterframe} 
	! initialized with 0

set/format I1	
	
	! Create a table for each frame to store the results
	! name: individual_+index
	! with column names: :x_pos , :y_pos, :x_fwhm, :y_fwhm
do count = 1 {imagenumber}
	create/table individual_{count} 4 null
	create/column individual_{count} :x_pos
	create/column individual_{count} :y_pos
	create/column individual_{count} :x_fwhm
	create/column individual_{count} :y_fwhm
enddo

set/format
	
	
	! Create a table for the final individual results	
	! name: ind_results
	! columns: x_pos, y_pos, x_focus, y_focus, ave_focus
create/table ind_results 5 null
create/column ind_results :x_pos
create/column ind_results :y_pos
create/column ind_results :x_focus
create/column ind_results :y_focus
create/column ind_results :ave_focus	


	! Create a table for the temporary storage of data belonging to one object
	! name: star_table
	! columns: x_fwhm, y_fwhm, xscale, yscale
create/table star_table 3 {imagenumber} null
create/column star_table :x_fwhm
create/column star_table :y_fwhm
create/column star_table :scale



	! Keyword for storage of # of selected objects
define/local sel_number/i/1/1 
define/local loop/i/1/1 		! second loop variable



!!!!!!!!!!!!!!!!!!!!!!!!! MORE DEFINITIONS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

define/local halfnumber/i/1/1		! # of points on each side of the
halfnumber = (imagenumber-1)/2		! best focus: 7 images --> 3
halfbox = boxsize/2


! Create a table "storeresults" for saving averaged results
! in the x and y direction of each image
! Format: #rows= #images ; #columns= 5
! Column names: :xavfwhm, :yavfwhm, :xstddev, :ystddev
! :xscale holds the matching x-values for the final plot

create/table storeresults 5 {imagenumber} null
create/column storeresults :xavfwhm
create/column storeresults :yavfwhm
create/column storeresults :xstddev
create/column storeresults :ystddev
create/column storeresults :xscale
!*read/table storeresults			! Test
 

write/out 
write/out "The found objects are marked with a circle (s=stellar, n=not stellar), "
write/out "the rectangle defines the region of interest around the objects"
write/out




!!!!!!!!!!!!!!!!!!!  "FIND OBJECT" ROUTINE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! findobj.prg takes the input image frame "focus100" and finds a given number 
! of objects in this frame. The output tabel "focus100.tbl contains the
! columns :x_mark, :y_mark with the position of the objects and a 
! classification of the object type in the column :label;
! s for stellar, n for not stellar

	
	! create a display window and load the master frame
create/display 0
assign/display d,0

load {masterframe}


	! Start automatic find_object routine
set/midas output=logonly
	
@@ O2K_UTIL:/obs_macros/Tilt/findobj {masterframe} {objectnumber}
!@@ /export/home/o2k/fassbender/omega2000_fokus/findobj {masterframe} {objectnumber}	
					! input image, # of objects to be found
					! output table = focus100.tbl

	! focus100.tbl holds the positions of the objects
	! To avoid confusion, a table "find_object_results" is created
	! The positions of the objects are copied to this table, as well as
	! the object type. All other later results will be appended 
	! "find_object_results", which is used to generate "gausstable" 
	! with center/gauss

set/midas output=yes	
write/out
write/out


create/table find_object_results 10 200 null

	! copy the three relevant columns to find_object_results and
	! name them: :xcenter, :ycenter, :type
copy/tt {masterframe} :x_mark find_object_results :xcenter
copy/tt {masterframe} :y_mark find_object_results :ycenter
copy/tt {masterframe} :label find_object_results :type


clear/channel overlay                 	! clears the cursor marks

	! mark the found objects on the display
load/table find_object_results :xcenter :ycenter :type 1 3 3 -1
	! x_pos,y_pos, object type, circles(1), size(3), red(3), do not connect(-1)


	! Define the boxes 
	! Add the columns :xstart, :xend, :ystart, :yend to the table
	! find:object_results
compute/table find_object_results :xstart = :xcenter-{halfbox}
compute/table find_object_results :xend = :xcenter+{halfbox}
compute/table find_object_results :ystart = :ycenter-{halfbox}
compute/table find_object_results :yend = :ycenter+{halfbox}

	!Display the defined rectangles on the screen
draw/rectangle find_object_results F N 4
	! table with input_specs, F=world coord, N= not filled, 3=green

	! Compute the FWHM for the histogram selection and store the results
	! in gausstable
set/midas output=logonly

center/gauss {masterframe},find_object_results gausstable

set/midas output=yes

	

!!!!!!!!!!!!!!!!!!!! HISTOGRAM OBJECTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Plot XFWHM and YFWHM in a histogramm

create/graphics 			! open a graphics window	

set/graphics colour=1			! black lines for xfwhm plot

plot/histogram gausstable :xfwhm  ? 0.1,0,5 ? ? 1,1,45
		! automatic scaling ; binsize 0.1 , min=0, max=5 arcsec
		! style: staircase steps with hashing

label/graphic "X-,Y-FWHM-Histogram" 0.6,0.9,no 0 1 0
label/graphic "of all objects" 0.6,0.85,no 0 1 0
	! Label: position in no=normalized coord, 0 angle, 1=size, centerd
		
set/graphics colour=2			! red lines for overplot
overplot/histogram gausstable :yfwhm ? 0.1 ? ? 1,1,-60
					! bin_size 0.1



!!!!!!!!!!!!!!!!!!!!!  INTENSITY PLOT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Plot FWHM versus central intensity to identify galaxis
	
	! Make a new column :averagefwhm that holds the average of the two FWHM
compute/table gausstable :averagefwhm = sqrt(:xfwhm*:yfwhm)
copy/tt gausstable :averagefwhm find_object_results :averagefwhm 
	! copy this column to find_object_resutls
	
create/graphics 1
set/graphics colour=1 ltype=0 stype=4	! black, no lines, triangles

plot/table gausstable :icent :averagefwhm

label/graphic "Intensity plot of all objects" 0.5,0.9,no 0 1 0
label/graphic "red lines = cuts" 0.5,0.85,no 0 1 0
	! Label: position in no=normalized coord, 0 angle, 1=size, centerd

	! Cutting is done in: SELECTION 4)



!!!!!!!!!!!!!!!!!!!!  STELLAR OBJECT SELECTION  !!!!!!!!!!!!!!!!!!!!!

	! 1) Select objects that are of type s=stellar----------------

select/table find_object_results :type.eq."s*"
	
	

	! 2) Throw out objects that are too close to the CCD edge------------

	! Define local keywords x_start, x_end, y_start, y_end that hold the 
	! image size information in world coordinates
define/local x_start/r/1/1	
define/local x_end/r/1/1
define/local y_start/r/1/1
define/local y_end/r/1/1

	! Use the image descriptors of focus100 to initialize the keywords
x_start = {{masterframe},start(1)}
x_end = {{masterframe},start(1)}+{{masterframe},step(1)}*{{masterframe},npix(1)}

y_start = {{masterframe},start(2)}
y_end = {{masterframe},start(2)}+{{masterframe},step(2)}*{{masterframe},npix(2)}

write/out
read/keyword x_start,x_end,y_start,y_end	! display image size in wc


	! Transfer the keywords to the cut-off values by adding and 
	! subtracting the boxsize
x_start = x_start+{boxsize}
x_end = x_end-{boxsize}
y_start = y_start+{boxsize}
y_end = y_end-{boxsize}


	! Throw out all objects that are within one boxsize of edge
select/table find_object_results sel.and.:xcenter.gt.{x_start}
select/table find_object_results sel.and.:xcenter.lt.{x_end}
select/table find_object_results sel.and.:ycenter.gt.{y_start}
select/table find_object_results sel.and.:ycenter.lt.{y_end}


	
	! 3) Throw out saturated objects-------------------------------

	! Get statistical properties on subframes and append results to
	! find_object_results. The maximum intensity is stored in 
	! the column :max
set/midas output=logonly

statistics/image {masterframe} find_object_results ? ? ? find_object_results,A

set/midas output=yes

 
	! The cut-off is defined in the keyword pixelsaturation = 350,000
write/out
write/out "The cut-off for saturated objects is set to: {pixelsaturation}"
write/out

	! Throw out all objects that have maximum intensities > cut-off
select/table find_object_results sel.and.:max.lt.{pixelsaturation}



	! 4) Make an intensity plot to identify galaxies; use cuts to 
	!    elminate them--------------------------------------------
	
statistics/table find_object_results :averagefwhm
	! mean=outputr(3) , stddev=outputr(4)


	!++++++++++ can be optimized +++++++++++++++++++++++++++++++++++
galaxycutoff(1) = outputr(3)-1.5*outputr(4)	! low cutoff = mean -stddev
galaxycutoff(2) = outputr(3)+1.0*outputr(4)	! high cutoff = mean+0.5 stddev


	
	! Plot the cut-offs in the graph
set/graphics colour=2 ltype=1 stype=0	! red line, no symbols

overplot/line 4 0,{galaxycutoff(1)},wo 450000,{galaxycutoff(1)},wo	! dashed line
overplot/line 4 0,{galaxycutoff(2)},wo 450000,{galaxycutoff(2)},wo	! dashed line


	! Now the selction: throw out all objects that are outside the cutoffs
select/table find_object_results sel.and.:averagefwhm.gt.{galaxycutoff(1)}
select/table find_object_results sel.and.:averagefwhm.lt.{galaxycutoff(2)}


set/format F4.3 I1

write/out
write/out
write/out "The cut-offs for the galaxy selection are..."
write/out "low: {galaxycutoff(1)} arcsec"
write/out "high: {galaxycutoff(2)} arcsec"
write/out


!*read/table find_object_results		! check the selected entries


	! mark the selected objects on the display

load/table find_object_results :xcenter :ycenter ? 8 4 5 -1
	! x_pos,y_pos, object type, crosses(8), size(4), blue(5), do not connect

	! store number of selected objects in sel_number
sel_number = {outputi(1)}


write/out
write/out
write/out "The selected table entries are marked with a blue cross on the screen."
write/out "{outputi(1)} stellar objects were selected for further processing..."
write/out

set/format



!!!!!!!!!!!!!!!!!!  HISTOGRAM: SELECTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Plot a FWHM histogramm of the selected objects
	
create/graphics 2			! new graphics window for selected objects
					! repeat plotting procedure
set/graphics colour=1			! black lines for xfwhm plot
plot/histogram find_object_results :averagefwhm  ? 0.1,0,5 ? ? 1,1,45
		! automatic scaling ; binsize 0.1 , min=0, max=5 arcsec
		! style: staircase steps with hashing
		
label/graphic "Averaged FWHM" 0.6,0.9,no 0 1 0
label/graphic "of selected objects" 0.6,0.85,no 0 1 0
	! Label: position in no=normalized coord, 0 angle, 1=size, centerd




!!!!!!!!!!!!!!!!!!  ALL FRAMES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
	

	! A DO loop will apply center/gauss to all images and store the average
	! x,y-fwhm value and the stdev in the table "storeresults"

define/local imagename/c/1/20 {imagename1}		
	! keyword that holds current image name, 
	!initialized with name of first image
define/local runnumber/i/1/1 {imagename1(6:8)}	! framenumber of focus097.bdf



!-----------------------------------------------------------------

do count = 1 {imagenumber} 1		! loopvar=start, end, step

set/format I1

!*read/keywort imagename			! check imagename

set/midas output=logonly

center/gauss {imagename},find_object_results frameresults
	! use columns x,ystart & x,yend of find_object_results 
	! on the image focusxxx
	! and store results in the temporary table frameresults


!*set/midas output=yes			! if active, results are displayed on terminal
!*read/table frameresults :xfwhm,yfwhm	! Check output


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


	! save the individual results in the corresponding table
copy/tt frameresults :xcen individual_{count} :x_pos
copy/tt frameresults :ycen individual_{count} :y_pos
copy/tt frameresults :xfwhm individual_{count} :x_fwhm
copy/tt frameresults :yfwhm individual_{count} :y_fwhm



	! Imagename of next image has to be written in keyword "imagename"
set/format I3				! Integer has to be: 001
runnumber = {runnumber}+1
!*read/keyword runnumber
write/keyword imagename/c/6/3 {runnumber}

enddo					! end of do loop

!-------------------------------------------------------------------

set/midas output=yes			! turn terminal log back on
read/table storeresults			! check table entries
set/format				! back to default format



!!!!!!!!!!!!!!!!!!!!!  PLOT RESULTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
define/local xoffset/r/1/1 1			! x-offset for graph
define/local yoffset/r/1/1 0.5			! y-offset for graph


set/midas output=logonly

statistics/table storeresults :xavfwhm		! get min & max y-values

ymin = outputr(1)
ymax = outputr(2)

statistics/table storeresults :yavfwhm		! get min & max y-values

set/midas output=yes

if outputr(1) .lt. ymin ymin = outputr(1)	! choose lowest min
if outputr(2) .gt. ymax ymax = outputr(2)	! choose highest max

ymin = ymin-yoffset				! include an offset
ymax = ymax+yoffset
xmin = xmin-xoffset
xmax = xmax+xoffset				! final scale parameters 


	! set up an appropriate coordinate system 
create/graphics 3				! use graphics window 3
set/graphics colour=1 				! back to black lines

plot/axes {xmin},{xmax},100,100 {ymin},{ymax} ?-
 "Telescope Position (microns)" "FWHM (arcsec)" 
 		! 100=distance between big tickmarks; 100= distance between small
 		! --> large distance, so tickmarks do not show up
 		! Will be changed to different units later
 		
 		
	
	! Plot the data points

	! First for XFWHM
set/graphics ltype=0 stype=4			! no line ; triangles(4) as symbols

compute/table storeresults :xscale = sequence-1-{halfnumber}	
	! vector for the x-axis input, symmetric around 0		


!*read/table storeresults 

overplot/table storeresults :xscale :xavfwhm
		! x-scale is in arbitrary units e.g. -3,..,0,...,3
		
	
	! Now for YFWHM
set/graphics ltype=0 colour=2 stype=3		! red squares ,no line
overplot/table storeresults :xscale :yavfwhm



!!!!!!!!!!!!!!!!!!!!!  FIT AND PLOT PARABOLAS  !!!!!!!!!!!!!!!!!!!!!!!!!!!

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
create/table fitpoints 3 200		! 3 columns, 200 points each
create/column fitpoints :absc		! absc = abscissa
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

set/graphics ltype=1 colour=1 stype=0	! solid black line, no symbols
overplot/table fitpoints :absc :xfit	! XFWHM-parabola

set/graphics colour=2			! solid red lines
overplot/table fitpoints :absc :yfit	! YFWHM-parabola


	! Relabel the x-axis with the real telescope 
	! position in microns (instead of image number)
xmin = focusvalue-(m$abs(xmin)*stepsize)	! m$abs() = absolute value
xmax = focusvalue+(xmax*stepsize)
!*read/keyword xmin,xmax

set/graphics colour=1			! back to black
define/local smallticks/r/1/1		! distance between small tickmarks
smallticks = stepsize/5			! stepsize = dist between large tm
						
overplot/axes {xmin},{xmax},{stepsize},{smallticks} {ymin},{ymax} ?-
 "Telescope Position (microns)" "FWHM (arcsec)" 



!!!!!!!!!!!!!!!!!  BEST FOCUS AND SEEING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Calculate the best focus position from the minimum of the parabola

	! Keyword definitions
define/local seeingxy/r/2/1			! y_min of the two parabolas
define/local seeing/r/1/1
define/local bestfocusxy/r/2/1			! minimum of the two parabolas
define/local bestfocus/r/1/1			! final result


	! The minimum of the parabola y = a + b*x +c*x^2 is:
	! x_min = -b / 2c
bestfocusxy(1) = -{coefficients,:xcoeff,2}/(2*{coefficients,:xcoeff,3})
bestfocusxy(2) = -{coefficients,:ycoeff,2}/(2*{coefficients,:ycoeff,3})


	! Best focus position is arithmetic mean of the two minima
bestfocus = (bestfocusxy(1)+bestfocusxy(2))/2	! still arbitrary units
bestfocus = (bestfocus*stepsize)+focusvalue  	! in microns 
bestfocusxy(1) = (bestfocusxy(1)*stepsize)+focusvalue 
bestfocusxy(2) = (bestfocusxy(2)*stepsize)+focusvalue 


	! The y-value @ the minimum of a parabola is:
	! y_min = a - b^2/4c
seeingxy(1) = {coefficients,:xcoeff,1}-{coefficients,:xcoeff,2}**2/(4*{coefficients,:xcoeff,3})
seeingxy(2) = {coefficients,:ycoeff,1}-{coefficients,:ycoeff,2}**2/(4*{coefficients,:ycoeff,3})

	! The real seeing is the geometric mean of the two y_min values
seeing = m$sqrt(seeingxy(1)*seeingxy(2))



! Display the final results on the screen and in graphics window

set/format f6.3					! give out appropriate format

label/graphic "Seeing = {seeing} arcsec" 60,115,mm 0 1.5 1
	! x, y_position in mm ; angle size pos_ind=starting @ position
	
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

set/format					! back to default


	! Label symbols

overplot/symbol 4 35,30,mm
label/graphic "x-fwhm" 40,30,mm 0 1 1

set/graphics colour=2				! red

overplot/symbol 3 35,25,mm
label/graphic "y-fwhm" 40,25,mm 0 1 1






!!!!!!!!!!!!!!!!!!!! NOW GET LOCAL FOCUS POSITIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! define local keywords
define/local starfocusxy/r/2/1
define/local starfocus/r/1/1
define/local xpix_start/i/1/1
define/local ypix_start/i/1/1
define/local xpix_end/i/1/1
define/local ypix_end/i/1/1
define/local store_value/r/1/1
define/local table_value/r/1/1


	! copy xscale (-3,-2,-1,0,1,2,3)column to table star_table
copy/tt storeresults :xscale star_table :scale


write/out 
write/out "The focus position for each star is now determined individually..."
write/out



	! DO loop over all selected objects
!----------------------------------------------------------------------------------
do loop = 1 {sel_number}

set/format I1

	! store position of next object in table ind_results
write/table ind_results :x_pos @{loop} {individual_1,:x_pos,{loop}}	
write/table ind_results :y_pos @{loop} {individual_1,:y_pos,{loop}}


	! get data for next object and store it in star_table
	!::::::::::::::::::::::::::::::::::::::::::::::::::::::
	do count = 1 {imagenumber}

	compute/keyword table_value = m$value(individual_{count},:x_fwhm,@{loop})
	write/table star_table :x_fwhm @{count} {table_value}

	compute/keyword table_value = m$value(individual_{count},:y_fwhm,@{loop})
	write/table star_table :y_fwhm @{count} {table_value}

	enddo
	!::::::::::::::::::::::::::::::::::::::::::::::::::::::


	! fit parabolas to the data
	! in x-direction
set/midas output=logonly
regression/polynominal star_table :x_fwhm :scale 2	
		      	! table, dependent var., independ. var, degree

	! Save the coefficients
set/format

coefficients,:xcoeff,@1 = {outputd(1)}
coefficients,:xcoeff,@2 = {outputd(2)}
coefficients,:xcoeff,@3 = {outputd(3)}


	! Parabola for :YFWHM
regression/polynominal star_table :y_fwhm :scale 2
		      	! table, dependent var., independ. var, degree

	! Save the coefficients
coefficients,:ycoeff,@1 = {outputd(1)}
coefficients,:ycoeff,@2 = {outputd(2)}
coefficients,:ycoeff,@3 = {outputd(3)}

set/midas output=yes


	! store minima of parabolas 
	! The minimum of the parabola y = a + b*x +c*x^2 is:
	! x_min = -b / 2c
starfocusxy(1) = -{coefficients,:xcoeff,2}/(2*{coefficients,:xcoeff,3})
starfocusxy(2) = -{coefficients,:ycoeff,2}/(2*{coefficients,:ycoeff,3})


	! Best focus position is arithmetic mean of the two minima
starfocus = (starfocusxy(1)+starfocusxy(2))/2	! still arbitrary units
starfocus = (starfocus*stepsize)+focusvalue  	! in microns 
starfocusxy(1) = (starfocusxy(1)*stepsize)+focusvalue 
starfocusxy(2) = (starfocusxy(2)*stepsize)+focusvalue 


	! take difference to average focus value in microns
starfocus = starfocus-bestfocus
starfocusxy(1) = starfocusxy(1)-bestfocusxy(1)
starfocusxy(2) = starfocusxy(2)-bestfocusxy(2)


	! Save results
write/table ind_results :x_focus @{loop} {starfocusxy(1)}
write/table ind_results :y_focus @{loop} {starfocusxy(2)}
write/table ind_results :ave_focus @{loop} {starfocus}


	! get pixel start and end in world coordinates
xpix_start = {individual_1,:x_pos,{loop}}-8
xpix_end = {individual_1,:x_pos,{loop}}+8
ypix_start = {individual_1,:y_pos,{loop}}-8
ypix_end = {individual_1,:y_pos,{loop}}+8


	! value to write into image: 9*9 pixel
store_value = starfocus +100000

	! Write results in output image
write/image tilt_frame [{xpix_start},{ypix_start}:{xpix_end},{ypix_end}] {store_value} all

		
enddo
!----------------------------------------------------------------------------------
	
	
	
	! display results
write/out
write/out "The differential focus results in microns are:"
write/out

read/table ind_results

set/format

	!load image
create/display 1
assign/display d,1
clear/channel overlay

load/image tilt_frame sc=-2 ce=c cuts=99800,100200
	
