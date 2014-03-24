$! findobj.prg (for O2000 package)

! @@ findobj in_image

 defin/par p1 ? ? "Input image: "
 defin/par p2 10 num "Number of stars: "

 write/key in_ima/c/1/60 {p1}
 write/key n_obj/i/1/1 {p2}

 run O2K_UTIL:/obs_macros/Focus/findobj.exe

 del/key in_ima, n_obj

 return
