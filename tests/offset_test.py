import numpy
import astropy.io.fits as fits
from astropy import wcs
import math

def getWCSPointingOffsets(images_in, p_offsets_file="/tmp/offsets.pap"):

      # Init variables
      i = 0 
      offsets_mat = None
      pix_scale = 0.45
      ra0 = -1
      dec0 = -1
      x_pix = 0.0
      y_pix = 0.0
      offsets = numpy.zeros([len(images_in),2] , dtype=numpy.float32)
      
      # Reference image
      ref_image = images_in[0]
      print "REF_IMAGE=",ref_image
      try:
            h0 = fits.getheader(ref_image)
            w0 = wcs.WCS(h0)
            ra0 = w0.wcs_pix2world(x_pix, y_pix, 1)[0]
            dec0 = w0.wcs_pix2world(x_pix, y_pix, 1)[1]
            print "RA0= %s DEC0= %s"%(ra0,dec0)
      except Exception,e:
          raise e
        
      
      offset_txt_file = open(p_offsets_file, "w")
      for my_image in images_in:
        try:
              h = fits.getheader(my_image)
              w = wcs.WCS(h)
              ra = w.wcs_pix2world(x_pix, y_pix, 1)[0]
              dec = w.wcs_pix2world(x_pix, y_pix, 1)[1]
              print "RA[%d]= %s DEC[%d]= %s"%(i, ra, i, dec)
              # Assummed that North is up and East is left
              offsets[i][0] = ((ra - ra0)*3600*math.cos(dec0/57.296)) / float(pix_scale)
              offsets[i][1] = ((dec0 - dec)*3600) / float(pix_scale) 
              print "offset_ra  = %s"%offsets[i][0]
              print "offset_dec = %s"%offsets[i][1]
              offset_txt_file.write(my_image + "   " + "%.6f   %0.6f"%(offsets[i][0], offsets[i][1]))
              i+=1
        except Exception,e:
          raise e
        
      offset_txt_file.close()
      
      # Write out offsets to file
      # numpy.savetxt(p_offsets_file, offsets, fmt='%.6f')
      print offsets
      
      return offsets

images_in = ["/data2/out/Q04/NGC752_focus0001.Q04.skysub.ast.fits", 
             "/data2/out/Q04/NGC752_focus0002.Q04.skysub.ast.fits", 
             "/data2/out/Q04/NGC752_focus0003.Q04.skysub.ast.fits",
             "/data2/out/Q04/NGC752_focus0004.Q04.skysub.ast.fits",
             "/data2/out/Q04/NGC752_focus0005.Q04.skysub.ast.fits"]

getWCSPointingOffsets(images_in,
                              p_offsets_file="/tmp/offsets.pap")
