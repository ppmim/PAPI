from __future__ import (absolute_import, division,
                        unicode_literals)
import numpy as np
from astropy import log
import astropy.io.fits as fits
import sys

def _clean_masked_pixels(data, mask, size=5, exclude_mask=None):
    """
    Clean masked pixels in an image.  Each masked pixel is replaced by
    the median of unmasked pixels in a 2D window of ``size`` centered on
    it.  If all pixels in the window are masked, then the window is
    increased in size until unmasked pixels are found.

    Pixels in ``exclude_mask`` are not cleaned, but they are excluded
    when calculating the local median.
    """

    assert size % 2 == 1, 'size must be an odd integer'
    assert data.shape == mask.shape, \
        'mask must have the same shape as image'
    ny, nx = data.shape
    
    #mask_coords = np.argwhere(mask)
    mask_coords = np.argwhere(mask==0)
    
    if exclude_mask is not None:
        assert data.shape == exclude_mask.shape, \
            'exclude_mask must have the same shape as data'
        maskall = np.logical_or(mask, exclude_mask)
    else:
        maskall = mask
        
    #mask_idx = maskall.nonzero()
    mask_idx = (maskall==0)
    data_nanmask = data.copy()
    data_nanmask[mask_idx] = np.nan
    print "Number of BPs = %d" % np.isnan(data_nanmask).sum()

    nexpanded = 0
    for coord in mask_coords:
        y, x = coord
        median_val, expanded = _local_median(data_nanmask, x, y, nx, ny,
                                             size=size)
        data[y, x] = median_val
        if expanded:
            nexpanded += 1
    if nexpanded > 0:
        log.info('    Found {0} {1}x{1} masked regions while '
                 'cleaning.'.format(nexpanded, size))
    return data


def _local_median(data_nanmask, x, y, nx, ny, size=5, expanded=False):
    """Compute the local median in a 2D window, excluding NaN."""
    
    hy, hx = size // 2, size // 2
    x0, x1 = np.array([x - hx, x + hx + 1]).clip(0, nx)
    y0, y1 = np.array([y - hy, y + hy + 1]).clip(0, ny)
    region = data_nanmask[y0:y1, x0:x1].ravel()
    goodpixels = region[np.isfinite(region)]
    if len(goodpixels) > 0:
        median_val = np.median(goodpixels)
    else:
        newsize = size + 2     # keep size odd
        median_val, expanded = _local_median(data_nanmask, x, y, nx, ny,
                                             size=newsize, expanded=True)
    return median_val, expanded

def cleanBadPixels( input_image, bpm ):
    
    input_data, input_header = fits.getdata(input_image, header=True)
    input_mask = fits.getdata(bpm, header=False)
    
    new_data = _clean_masked_pixels(input_data, input_mask)
    
    # Write cleaned data
    fits.writeto(input_image.replace(".fits", "_c.fits"), 
                 new_data, header=input_header, clobber=True)
    
    return

if __name__ == "__main__":
    
    cleanBadPixels(sys.argv[1], sys.argv[2])
    
    sys.exit(0)
    