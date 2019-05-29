#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 12:41:20 2019

@author: jesper
"""

from astropy.io import fits
from scipy.optimize import curve_fit
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import time

plt.close('all')

hdu = fits.open("hd201585_harps-n_lines.fits")
data = hdu[1].data

lines = 1-data['LINES']
vel_vector = data['V'][0]
bjd = data['BJD']

def gauss(x, a, b, c, d):
    return a * np.exp(-(x-b)**2/(2 * c**2)) + d

print("Fitting gaussians...")
time.sleep(0.5)
vels = []
v_unc = []
for i in tqdm(range(len(lines))):
    fit, pov = curve_fit(gauss, vel_vector, lines[i])
    if i == 6:
        plt.figure()
        plt.plot(vel_vector, lines[i], 'r-', label='Data')
        plt.plot(np.linspace(-200, 200, 1000), gauss(np.linspace(-200, 200, 1000), *fit), 'b--', label='Gaussian fit')
        plt.xlabel('Velocity [km/s]')
        plt.ylabel('CCF')
        plt.legend()
        plt.tight_layout()
        plt.show()
    vels.append(fit[1])
    v_unc.append(np.sqrt(np.diag(pov))[1])
vels = np.array(vels)

plt.figure()
plt.errorbar(bjd, vels, yerr=v_unc, fmt='r.')
plt.xlabel('BJD')
plt.ylabel('Velocity [km/s]')
plt.tight_layout()
plt.show()