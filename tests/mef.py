import pyfits


myflat=pyfits.open("/scr/HAWK-I/hawki-demo-0.3/raw/JITTER/HAWKI.2008-12-03T04:38:51.789.fits")
hdr = myflat[0].header.copy()
hdr.update('OBJECT','MASTER_GAINMAP')
hdr.add_history('Gain map based on myflat')
prihdr=pyfits.PrimaryHDU(None,hdr)
fo = pyfits.HDUList()

# Add primary header to output file...
fo.append(prihdr)


# Add each extension
for chip in range(0, 4):
    hdu = pyfits.ImageHDU(data=myflat[chip+1].data, header=myflat[chip+1].header)
    hdu.scale('float32') # importat to set first data type ??
    hdu.header.update('EXTVER',1)
    fo.append(hdu)

fo.writeto("/tmp/prue.fits")
fo.close()
del fo, hdu
