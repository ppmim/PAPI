! Dither positions from file
!  Use P6 with file name to access own pattern of 20 positions
!  This will use the generic repetition pattern, if not wanted 
!         use second call with other dither pattern.
! To implement: Move fctrl for log to place after definition of dither pattern.
!               Logic for P6:
!               Try to open file P6, if successfull read data.
!               Otherwise it has to be 0 or 1 as before.


define/local fctrl/i/1/2 0,0

if m$exist("{P6}") .eq. 1 then
   define/local input_buffer/c/1/80 " "
   define/local X_offset/i/1/20 0 all
   define/local Y_offset/i/1/20 0 all
   define/local values/r/1/2 0,0
   define/local X_back/i/1/1 0
   define/local Y_back/i/1/1 0
   define/local i/i/1/1 0

   open/file {P6} read fctrl
   ! First line is a comment
   read/file {fctrl(1)} input_buffer

   set/format I1

   do i = 1 20
      read/file {fctrl(1)} input_buffer
      if m$index(input_buffer," ") .gt. 0  then
         write/out "         Data file for offsets must not contain blanks!"
         write/out "         Check line {i} and following..."
!         $play -q $GEIRS_DIR/SOUNDS/crash.au
         goto exit
      endif
      if fctrl(2) .eq. -1 then
         write/out "Unexpected EOF reached!   Abort ..."
!         $play -q $GEIRS_DIR/SOUNDS/crash.au
         goto exit
      else
         write/key values/r/1/2 {input_buffer}
         x_offset({i}) = {values(1)}*10
         y_offset({i}) = {values(2)}*10
         X_back = X_back + x_offset({i})
         Y_back = Y_back + y_offset({i})
      endif
   enddo
   exit:
   close/file {fctrl(1)}
else
!   as before
endif
