from astropy.io import fits
from astropy import wcs
import os
import yaml
import numpy as np


def main(surface_values, ximg_arcsec, yimg_arcsec, outdir):
    x = ximg_arcsec
    y = yimg_arcsec
    
    sv_ct_no0 = surface_values
    sv_ct_no0 = sv_ct_no0*(sv_ct_no0>0)

    w = wcs.WCS(naxis=2)
    
    # what is the center pixel of the XY grid.
    w.wcs.crpix = [x.shape[0]/2, y.shape[1]/2]
    
    # what is the galactic coordinate of that pixel.
    w.wcs.crval = [x.min()+np.abs(x.max() - x.min())/2, y.min()+np.abs(y.max() - y.min())/2]

    # what is the pixel scale in lon, lat.
    w.wcs.cdelt = np.array([0.2, 0.2])

    # write the HDU object WITH THE HEADER
    header = w.to_header()
    hdul = fits.PrimaryHDU(data=sv_ct_no0, header=header)
    hdul.writeto(os.path.join(outdir, 'fits/surface_image.fits'), overwrite=True)

    