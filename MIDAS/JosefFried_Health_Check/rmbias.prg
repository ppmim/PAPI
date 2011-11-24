
! rmbias.prg

defi/local chip1/c/1/1 "a"
define/loc chip2/c/1/1 "b"
define/loc chip3/c/1/1 "c"
define/loc chip4/c/1/1 "d"
define/loc chip5/c/1/1 "e"
define/loc chip6/c/1/1 "f"
define/loc chip7/c/1/1 "g"
define/loc chip8/c/1/1 "h"
define/loc n/i/1/1 1            !Laufvariable fuer 7 Einzelchips

!echo/on
def/par p1 ? C "enter filename to remove its bias from"

avera/col biasof{p1} = {p1}a @2099,@2137
filter/smooth biasof{p1} smoothbiasof{p1} 50
grow/ima imabiasof{p1} = biasof{p1} 0.,1.,2142
rotate/counter imabiasof{p1} rotbiasof{p1}
copy/dd {p1}a *,1 rotbiasof{p1}
comp/ima rb{p1}a = {p1}a -rotbiasof{p1}
$rm biasof{p1}.bdf
$rm smoothbiasof{p1}.bdf
$rm imabiasof{p1}.bdf
$rm rotbiasof{p1}.bdf




do n = 2 8
	
	set/format I1
	write/out chip {chip{n}}
	avera/col biasof{p1}{chip{n}} = {p1}{chip{n}} @2099,@2137
	filter/smooth biasof{p1}{chip{n}} smoothbiasof{p1}{chip{n}} 50
	grow/ima imabiasof{p1}{chip{n}} = biasof{p1}{chip{n}} 0.,1.,2142
	rotate/counter imabiasof{p1}{chip{n}} rotbiasof{p1}{chip{n}}
	copy/dd {p1}{chip{n}} *,1 rotbiasof{p1}{chip{n}}
	comp/ima rb{p1}{chip{n}} = {p1}{chip{n}} -rotbiasof{p1}{chip{n}}

$rm biasof{p1}{chip{n}}.bdf
$rm smoothbiasof{p1}{chip{n}}.bdf
$rm imabiasof{p1}{chip{n}}.bdf
$rm rotbiasof{p1}{chip{n}}.bdf
 
enddo

cop/dd {p1} *,3 rb{p1}a
