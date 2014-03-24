write/out "Mark object first, then sky background ..."

clear/chan over
get/curs CR_mon N ? 2,2
stat ? CR_mon ? ? FNNN CR_mon,A

define/local xstart/i/1/1 0
define/local ystart/i/1/1 0
define/local xend/i/1/1 0
define/local yend/i/1/1 0
define/local n_star/i/1/1 0
define/local n_back/i/1/1 0

define/local star/r/1/1 0
define/local back/r/1/1 0
define/local signal/r/1/1 0

xstart = m$value(CR_mon,:Xstartpix,@1)
ystart = m$value(CR_mon,:Ystartpix,@1)
xend = m$value(CR_mon,:Xendpix,@1)
yend = m$value(CR_mon,:Yendpix,@1)

n_star = (xstart-xend+1)*(ystart-yend+1)

xstart = m$value(CR_mon,:Xstartpix,@2)
ystart = m$value(CR_mon,:Ystartpix,@2)
xend = m$value(CR_mon,:Xendpix,@2)
yend = m$value(CR_mon,:Yendpix,@2)

n_back = (xstart-xend+1)*(ystart-yend+1)

star = m$value(CR_mon,:tot_intens,@1)
back = m$value(CR_mon,:tot_intens,@2)

signal = star - back*n_star/n_back

write/out "Counts in object = {signal}"

return

