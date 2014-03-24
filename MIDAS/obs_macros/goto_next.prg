!  goto_next                HJR   19-Sep-03
!
!  Set telescope from center position to dither position n-1
!  in preparation for o2k/dither with a starting position other than 1

define/par P1 ? ? "Index of next exposure : "
define/par P2 0 N/A "Integer pixels [0] or non-integer [1] pixels : "


	! for non-integer pixel offsets
if {P2} .eq. 1 then

  define/local x_offset/i/1/20 0,100,-150,180,-230,-50,200,130,-380,230
  write/keyword x_offset/i/11/20 -210,150,-90,320,-50,-220,140,-240,40,260

  define/local y_offset/i/1/20 0,100,50,-220,-30,150,-200,180,-80,250,
  write/keyword y_offset/i/11/20 -70,-310,160,-100,290,-100,-50,-190,370,-400
  
  define/local X_back/I/1/1 -130
  define/local Y_back/i/1/1 200
	
	
	! define repetion offsets
  define/local x_repetition/i/1/20 0,-20,0,20,0,-30,0,30,0,-20,
  write/keyword x_repetition/i/11/20 -30,-40,-30,-20,0,20,30,40,30,20

  define/local y_repetition/i/1/20 20,0,-20,0,30,0,-30,0,40,30
  write/keyword y_repetition/i/11/20 20,0,-20,-30,-40,-30,-20,0,20,30

!..............................................................................

else	! integer pixel offsets	: multiples of 9 = 2pixels

  define/local x_offset/i/1/20 0,99,-153,180,-225,-45,198,126,-378,225
  write/keyword x_offset/i/11/20 -207,153,-90,315,-45,-216,144,-234,36,252

  define/local y_offset/i/1/20 0,99,54,-216,-27,153,-198,180,-90,252
  write/keyword y_offset/i/11/20 -72,-306,162,-99,288,-99,-54,-189,360,-396

  define/local X_back/I/1/1 -135
  define/local Y_back/i/1/1 198
	
	
	! define repetion offsets--> modified to yield sum=0
  define/local x_repetition/i/1/20 0,-18,0,18,0,-27,0,27,0,-18,
  write/keyword x_repetition/i/11/20 -27,-36,-27,-18,0,18,27,36,27,18

  define/local y_repetition/i/1/20 18,0,-18,0,27,0,-27,0,36,27
  write/keyword y_repetition/i/11/20 18,0,-18,-27,-36,-27,-18,0,18,27

  write/out "Integer pixel offsets were defined..."

endif

define/local repetitions/i/1/1 0
define/local n_dither/i/1/1 0
repetitions = ({P1}-1)/20
n_dither = {P1} - 1 - repetitions*20

define/local act_pos/i/1/2 0,0
define/local ix/i/1/1 0

act_pos(1) = x_repetition({repetitions})
act_pos(2) = y_repetition({repetitions})

do ix = 1 n_dither
	act_pos(1) = act_pos(1) + x_offset({ix})
	act_pos(2) = act_pos(2) + y_offset({ix})
enddo

set/format I1 F5.2
write/out
write/out "Repetitions {repetitions} / dither position {n_dither}"
write/out "Offset {act_pos(1)} / {act_pos(2)} 1/10 arcsec"
write/out

$ {tecs_script}/t_coord_system xy
$ {tecs_script}/t_offset {act_pos(1)} {act_pos(2)}

return
