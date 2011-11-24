
!
!    km_infits.prg
!
!  08/Jan/99
!


def/par p1 ZYX ? " number of file >"
def/par p2 wfi0 ?  " output name >"
def/par p3 /data/ ?" directory of fits files >"
def/par p4 wfi ?  " input name >"
!
if P1 .eq. "ZYX" then
   writ/out " call:"
   writ/out "       INPUT/wfi  N_file  OUT_id  directory    IN_id"
   writ/out " e.g.  INPUT/wfi  7528    fram    /data/mine/  wfi"
   writ/out
   return
endif

def/loc N/I/1/1 {p1}

set/format I5

if p4 .ne. "wfi" then
	set/format I4
        intape/fits {N} {p2} {p3}{p4}{n}.fits 
else
        intape/fits {N} wfi {p3}{p4}{n}.fits
endif

!
return


