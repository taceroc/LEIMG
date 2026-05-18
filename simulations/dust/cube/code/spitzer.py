import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate
from astropy.io import fits
from astropy import units as u
import astropy.cosmology.units as cu
from astropy.wcs import WCS
from astropy.coordinates import Angle

# from astropy import units as u
from astropy.nddata import CCDData

import ccdproc as ccdp
from photutils import detect_sources

from astropy.stats import sigma_clipped_stats
from photutils.datasets import load_star_image

from photutils.detection import DAOStarFinder

from astropy.visualization import SqrtStretch

from astropy.visualization.mpl_normalize import ImageNormalize

from photutils.aperture import CircularAperture



import sys
# sys.path.append('/content/drive/MyDrive/LE2023/dust/code')
path = "..\\Data\\tgt-1-selected_Post_BCDs\\r12914432\\ch3\\pbcd\\SPITZER_I3_12914432_0000_6_E8357568_maic.fits"
datas = fits.open(path)
header = datas[0].header

data = datas[0].data
cmax = np.nanpercentile(data.flatten(), 99, axis=0)
cmin = np.nanpercentile(data.flatten(), 1, axis=0)

data_nonan = data[500:5000, 500:-1000].copy()
plt.imshow(data_nonan, clim=(cmin,cmax), cmap="Greys_r", origin='lower')

### Indicate $\eta$ carinae location, as in [Light echoes reveal an unexpectedly cool η Carinae during its 19th-century Great Eruption, Rest et all 2011](https://arxiv.org/pdf/1112.2210.pdf)
rahms, dechms = ("10h44m12.127s", "-60d16m01.69s")
print(Angle(rahms).degree)
print(Angle(dechms).degree)

#convert to pixel location
ra = Angle(rahms).degree
dec = Angle(dechms).degree
# f = fits.open('f242_e0_c98.fits')
mywcs = WCS(header)
x, y = mywcs.all_world2pix(ra, dec, 1)
# print(x,y)

plt.imshow(data, clim=(cmin,cmax), cmap="Greys_r", origin="lower")
plt.scatter(x,y, s=80, facecolors='none', edgecolors='r')
dust_center_x = x - 300
dust_center_y = y + 600
plt.scatter(dust_center_x, dust_center_y, s = 80, facecolors='none', edgecolors='blue')
# plt.scatter(mywcs.all_world2pix(ra-0.5, dec-0.5, 1)[0],mywcs.all_world2pix(ra+0.5, dec, 1)[1],
            #  s=80, facecolors='none', edgecolors='r')
plt.xlim(5600, 7000)
plt.ylim(3500, 4500)
# Convert dust center to coordinates
ra_dust, dec_dust =  mywcs.all_pix2world(x - 300, y + 600, 0)
Angle(dec_dust, u.degree).to_string(unit=u.degree, sep=('deg', 'm', 's'))

## NEW FILE WITH THE CROP
newf = fits.PrimaryHDU()
newf.data = data[3400:4500, 5400:6700].copy()
newf.header = header.copy()
newf.header.update(mywcs[3400:4500, 5400:6700].to_header())
cmax_crop = np.nanpercentile(newf.data.flatten(), 99, axis=0)
cmin_crop = np.nanpercentile(newf.data.flatten(), 1, axis=0)
plt.imshow(newf.data, cmap="Greys_r", clim=(cmin,cmax), origin='lower')
mywcs_crop = WCS(newf.header)
x, y = mywcs_crop.all_world2pix(ra, dec, 1)

work_data = data[3400:4400, 5600:6500].copy()
mean, median, std = sigma_clipped_stats(work_data, sigma=3.0)  
print((mean, median, std))  
# find stars
daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)  
sources = daofind(work_data - median)  
for col in sources.colnames:  
    if col not in ('id', 'npix'):
        sources[col].info.format = '%.2f'  # for consistent table output

sources.pprint(max_width=76)

#mark stars
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=2.0)
norm = ImageNormalize(stretch=SqrtStretch())

plt.imshow(work_data, clim=(cmin,cmax), cmap="Greys_r", origin='lower')#cmap='Greys', origin='lower', norm=norm, interpolation='nearest')

apertures.plot(color='blue', lw=1.5, alpha=0.5)
# plt.show()

mask_data = work_data.copy()
for aperture in apertures:
    mask0 = aperture.to_mask().to_image(work_data.shape)
    masks = np.logical_not(mask0).astype(int)
    mask_data = mask_data * masks


#select a batch to interpolate the empty sopaces

data_chunk = mask_data.copy()[660:960,360:660]

def interpolate_missing_pixels(
        image: np.ndarray,
        method: str = 'cubic',
        fill_value: int = 0):
    """
    :param image: a 2D image
    :param method: interpolation method, one of
        'nearest', 'linear', 'cubic'.
    :param fill_value: which value to use for filling up data outside the
        convex hull of known pixel values.
        Default is 0, Has no effect for 'nearest'.
    :return: the image with missing values interpolated
    """
    from scipy import interpolate

    h, w = image.shape[:2]
    # vv = img[np.argwhere(img!=0)[:,0], np.argwhere(img!=0)[:,1]]
    vv = np.argwhere(image!=0)
    xx, yy = np.meshgrid(np.arange(w), np.arange(h))
    
    known_x = np.arange(w)[vv[:, 0]]
    known_y = np.arange(h)[vv[:, 1]]
    known_v = image[np.argwhere(image!=0)[:,0], np.argwhere(image!=0)[:,1]]
    missing_x = xx[np.argwhere(image==0)[:,0], np.argwhere(image==0)[:,1]]
    missing_y = yy[np.argwhere(image==0)[:,0], np.argwhere(image==0)[:,1]]

    interp_values = interpolate.griddata(
        (known_x, known_y), known_v, (missing_x, missing_y),
        method=method, fill_value=fill_value)

    interp_image = image.copy()
    interp_image[missing_y, missing_x] = interp_values

    return interp_image

inter_data_chunk = interpolate_missing_pixels(image = data_chunk, method="cubic")
data_chunk_norm = (inter_data_chunk - np.min(inter_data_chunk)) / (np.max(inter_data_chunk) - np.min(inter_data_chunk)) * inter_data_chunk.shape[0]
# np.save("..\\Data\\spitzer_batch1_2d.npy", data_chunk_norm)

data_chunk_small = mask_data.copy()[760:760+44,460:460+44]
data_chunk_norm = (data_chunk_small - np.min(data_chunk_small)) / (np.max(data_chunk_small) - np.min(data_chunk_small)) * data_chunk_small.shape[0]
# np.save("..\\Data\\spitzer_batch_2d.npy", data_chunk_norm)
