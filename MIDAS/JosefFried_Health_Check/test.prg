defi/par p1 test ima "root=?"
defi/par p2 1,21 n "files from,to=?"
defi/par p3 mean ima "mean outputfile=?"
defi/par p4 sigma ima "stdev outputfile=?"




defi/loc root/c/1/20 {p1}
defi/loc n/i/1/2 {p2}

defi/loc mean/c/1/20 {p3}
defi/loc sigma/c/1/20 {p4}

defi/loc nf/i/1/1 0
nf = n(2)-n(1)+1


comp mean = {root}{n(1)}*0.
comp sigma = {root}{n(1)}*0.




cop/dk {root}{n(1)} npix npix


defi/loc m/i/1/1 0

crea/ima z 2,{npix(1)},{nf}
comp z = z*0.


inputc = m$time()
write/out start: {inputc}

do m = 1 npix(2)  ! loop over rows in input images

    comp z = 0.
    do n = 1 nf
	extra/ima ext = {root}{n} [<,@{m}:>,@{m}]
	insert/ima ext z @1,@{n} 
    enddo
    
    stat/ima z COLUMN ? ? outtab=outm,image,Mean >NULL
    insert/image outm {mean} @1,@{m}
 
    stat/ima z COLUMN ? ? outtab=outs,image,Stddev >NULL
    insert/image outs {sigma} @1,@{m}

     
enddo

write/out output -> files {mean}, {sigma}
inputc = m$time()
write/out end : {inputc}



