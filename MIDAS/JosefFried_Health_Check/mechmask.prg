!
! mechmask.prg
!
! computes CNC program file
!
! jwf apr 1998
!

defi/par p1 ? ima "table=?"

defi/par p2 ? n " output file [ 4 digit number only  neg.number -> exit]= ? "
if {p2} .lt. 0 goto exit


defi/local yps/r/1/1 0.
defi/local ype/r/1/1 0.
defi/local xp/r/1/1 0.
defi/local loop/i/1/1 1
defi/local numb/i/1/1 10
defi/local nslits/i/1/4 0,0
defi/local file/i/1/2 0,0

! used to correct length of slits for min. overlap
defi/local fras/r/1/1 0.25
fras = fras/4.

open/file {p2} write file

write/file {file} %{p2} ( -D11-000000- -06.05.98- Maske{p2})
write/file {file} ?
write/file {file} 0000
write/file {file} %{p2}*T

!define tools
write/file {file} T1 I/M:M R0.5 L0 FE300 FZ100 S2000 (Stichel 0.93)
write/file {file} T2 I/M:M R0.8 L11 FE300 FZ100 S1500 (Bohrer 1.7)
write/file {file} T3 I/M:M R0.2 L26.356 FE300 FZ100 S2000 (Bohrer 0.35)

write/file {file} ?
write/file {file} 0000
write/file {file} " "
write/file {file} %{p2}*% ( -D11-000000- -06.05.98- Maske{p2})


write/file {file} N10 G0 Z100 G17 T1 M3 M62

!  slits

sel/tab {p1}_mask :type.eq.0
cop/tab {p1}_mask z
nslits = {outputi(1)}
set/format i2 f6.3
write/out  {nslits}  slits

do loop = 1 {nslits}
	numb = {numb}+10
	xp = m$value(z,:x_m,@{loop}) 
	yps = m$value(z,:y_m_s,@{loop})-fras
	ype = m$value(z,:y_m_e,@{loop})+fras
	write/file {file} N{numb} G0 X{xp} Y{ype}
        numb = {numb}+10
	write/file {file} N{numb} G1 Z0
	numb = numb+10
	write/file {file} N{numb} G1 X{xp} Y{yps}
	numb = {numb}+10
	write/file {file} N{numb} G1 X{xp} Y{ype}
	numb = numb+10
	write/file {file} N{numb} G0 Z5
enddo

! Bohre Loecher an

sel/tab {p1}_mask :type.ne.0
cop/tab {p1}_mask z
nslits = {outputi(1)}
set/format i2 f6.3
write/out  {nslits}  holes

do loop = 1 {nslits}
	numb = {numb}+10
	xp = m$value(z,:x_m,@{loop}) 
	yps = m$value(z,:y_m,@{loop}) 
	write/file {file} N{numb} G0 X{xp} Y{yps}
        numb = {numb}+10
	write/file {file} N{numb} G1 Z0.3
	numb = numb+10
	write/file {file} N{numb} G0 Z5
enddo







! large holes for stars
! change tools
numb = numb+10
write/file {file} N{numb} G0 Z100
numb = numb+10
write/file {file} N{numb} M6
numb = numb+10
write/file {file} N{numb} G17 T2 M3 M62 [1.7]

sel/tab {p1}_mask :type.eq.1
defi/local nstars/i/1/1 {outputi(1)}
write/out large holes for {nstars} stars
cop/tab {p1}_mask z
do loop = 1 {nstars}
	xp = m$value(z,:x_m,@{loop}) 
	yps = m$value(z,:y_m,@{loop}) 
	numb = numb+10
	write/file {file} N{numb} G0 X{xp} Y{yps}
	numb = numb+10
	write/file {file} N{numb} G0 Z2
	numb = numb+10
	write/file {file} N{numb} G1 Z-2
	numb = numb+10
	write/file {file} N{numb} G0 Z10
enddo

! ref. holes
! change tools

numb = numb+10
write/file {file} N{numb} G0 Z100
numb = numb+10
write/file {file} N{numb} M6
numb = numb+10
write/file {file} N{numb} G17 T3 M3 M62 [0.35]


sel/tab {p1}_mask :type.eq.2
nstars = {outputi(1)}
write/out {nstars}  ref.holes...
cop/tab {p1}_mask z
defi/local xo/r/1/1 0.
defi/local yo/r/1/1 0.

do loop = 1 {nstars}
	xp = m$value(z,:x_m,@{loop}) 
	yps = m$value(z,:y_m,@{loop}) 
	numb = numb+10
	write/file {file} N{numb} G0 X{xp} Y{yps}
	numb = numb+10
	write/file {file} N{numb} G0 Z2
	numb = numb+10
	write/file {file} N{numb} G1 Z-0.3
	numb = numb+10
	write/file {file} N{numb} G0 Z10
enddo

numb = numb+10
write/file {file} N{numb} G0 Z100 M30
write/file {file} ?
write/file {file} A44C
close/file {file}
write/out output -> file {p2}



exit:
