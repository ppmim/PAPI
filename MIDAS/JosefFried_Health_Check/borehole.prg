!
!
!

write/out select area lower left,upper right:

get/cur cursor

defi/local ll/i/1/2 1,1
defi/local ur/i/1/2 1,1

ll(1) = m$value(cursor,#4,@1)
ll(2) = m$value(cursor,#5,@1)
ur(1) = m$value(cursor,#4,@2)
ur(2) = m$value(cursor,#5,@2)


stat/ima CUBE [{ll(1)},{ll(2)},<:{ur(1)},{ur(2)},>] ? ? ? ? P

