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

#Lines er CCF, Så det vi plotter er styrken af CCF'en som funktion af RV'en, til forskellige tidspunkter. 
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

#Vi splitter data op i de to observerede transitter
bjd1 = bjd[np.where(bjd < 2457950)]
bjd2 = bjd[np.where(bjd > 2457950)]
lines1 = lines[np.where(bjd < 2457950)]
lines2 = lines[np.where(bjd > 2457950)]

#Vi tager et gns. af alle 'lines' der vurderes til at ligge uden for transitten
OT_lines1 = np.vstack((lines1[0:7,:],lines1[-7:,:]))
OT_gns1 = np.mean(OT_lines1,0)

#Vi trækker dette gns. fra alle 'lines'
norm_lines1 = lines1-OT_gns1

plt.figure()
plt.plot(vel_vector,norm_lines1[0])
plt.plot(vel_vector,lines1[0])
plt.plot(vel_vector,OT_gns1)
plt.plot(vel_vector,norm_lines1[20])
plt.tight_layout()
plt.show()

locx = np.linspace(0,250,6)
tick_labx = np.round(np.linspace(vel_vector[200],vel_vector[450],6),2)
locy = np.array([0,10,20,30,40])
tick_laby = np.round(bjd1[locy]-2457939,2)

plt.figure(figsize=(7,7))
color_map = plt.imshow(norm_lines1[:,200:451], aspect = 'auto')
color_map.set_cmap('gray')
plt.xticks(locx,tick_labx)
plt.yticks(locy,tick_laby)
plt.xlabel('RV')
plt.ylabel('Time [BJD]-2457939')
plt.show()

