obsutil
if (access("/tmp/focus_seq.txt"))
{
   starfocus(images="@/tmp/focus_seq.txt", focus = "T_FOCUS", fstep = "", nexposures = 1, 
             coords = "mark1", wcs = "logical", size = "GFWHM", scale = 1, sbuffer = 10, 
             swidth= 10, radius = 5, satura = INDEF, ignore_sat =no , 
             imagecur = "", display =yes, frame = 1, graphcur = "", 
             logfile = "starfocus.log")
             
    delete("/tmp/focus_seq.txt", yes, verify=no)
}
else
{ 
   print ("\n")
   print("*******************************************************")
   print("ERROR:  Cannot find file /tmp/focus_seq1.txt ")
   print("*******************************************************")
   print ("\n")
   
}
;
