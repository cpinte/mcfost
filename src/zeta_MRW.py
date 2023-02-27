import numpy as  np
import matplotlib.pyplot as plt
from copy import deepcopy
from decimal import Decimal
from astropy.io import fits
from scipy import interpolate

def calculate_zeta(y, n_its):
    n = np.arange(1, n_its)
    zeta = np.sum(2*(-1)**(n+1)*y**(n**2))
    # print(np.sum(zeta)
    return np.min([zeta, 1])


y_values = np.linspace(0, 1, 1000, endpoint=True)
zeta_values = np.zeros_like(y_values)

for idx, y in enumerate(y_values):
    zeta_values[idx] = calculate_zeta(y, 1000)

data = np.stack([y_values, zeta_values], axis=1)

hdu = fits.PrimaryHDU(data)
hdu.header['NAXIS'] = 2

hdu.writeto('zeta_MRW.fits', overwrite=True)

