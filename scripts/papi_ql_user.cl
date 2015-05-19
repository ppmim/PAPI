obsutil
if (access("home$focus_seq.txt"))
{
   starfocus(images="@home$focus_seq.txt", focus = "T_FOCUS", fstep = "", nexposures = 1, 
             coords = "mark1", wcs = "logical", size = "GFWHM", scale = 1, sbuffer = 10, 
             swidth=10, radius = 5, satura = INDEF, ignore_sat =no , 
             imagecur = "", display =yes, frame = 1, graphcur = "", 
             logfile = "starfocus.log")
             
    delete("home$focus_seq.txt", yes, verify=no)

    # Call next Python script to read the starfocus.log file and plot a more useful plot for the focus fitting.
    !/home/panic/bin/runStarfocus.py -o starfocus.log
}
else
{ 
   print ("\n")
   print("*******************************************************")
   print("ERROR:  Cannot find file iraf/focus_seq.txt ")
   print("*******************************************************")
   print ("\n")
   
}
;
