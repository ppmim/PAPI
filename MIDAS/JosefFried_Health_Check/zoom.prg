!
! zoom.prg
!
! jwf 4.12.1996
!
!
if p1(1:4) .eq. "help" then 
	write/out @@ zoom factor g/c 
	write/out zooms image by factor with object marked with cursor as center
	write/out c = center gauss  g = get cursor
	write/out
	return
endif




defi/par p1 ? n "zoom factor=?"
defi/par p2 g ima "g=def/c [cen/gauss or get/curs]"
if "{p2}" .eq. "c" then
   write/out center object
   cg
   load {idimemc} scale={p1} center={outputr(5)},{outputr(6)}
endif
if "{p2}" .eq. "g" then
   write/out mark position
   get/cursor
   load {idimemc} scale={p1} center={outputr(12)},{outputr(13)}
endif



exit:
