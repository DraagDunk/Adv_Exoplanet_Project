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
from rm_function import rm_function_fixed

plt.close('all')

hdu = fits.open("hd201585_harps-n_lines.fits")
data = hdu[1].data

#Lines er CCF, Så det vi plotter er styrken af CCF'en som funktion af RV'en, til forskellige tidspunkter. 
lines = 1-data['LINES']
vel_vector = data['V'][0]
bjd = data['BJD']

def gauss(x, a, b, c, d):
    return a * np.exp(-(x-b)**2/(2 * c**2)) + d

def inv_gauss(x, a, b, c, d):
    return -(a * np.exp(-(x-b)**2/(2 * c**2))) + d

print("Fitting gaussians...")
time.sleep(0.5)
#%%
vels = []
v_unc = []
for i in tqdm(range(len(lines))):
    fit, pov = curve_fit(gauss, vel_vector, lines[i])
    if i == 10:
        plt.figure()
        plt.plot(vel_vector, lines[i], 'r-', label='Data')
        plt.plot(np.linspace(vel_vector[0],vel_vector[-1], 1000), 
                gauss(np.linspace(vel_vector[0],vel_vector[-1], 1000), *fit), 'b--', label='Gaussian fit')
        plt.xlabel('Velocity [km/s]')
        plt.ylabel('CCF')
        plt.xlim([vel_vector[0],vel_vector[-1]])
        plt.title('CCF during transit')
        plt.legend()
        plt.tight_layout()
        #plt.savefig('figures/CCF_it_fit.png')
        plt.show()
    vels.append(fit[1])
    v_unc.append(np.sqrt(np.diag(pov))[1])
vels = np.array(vels)
v_unc = np.array(v_unc)

#%%

#Vi splitter data op i de to observerede transitter
bjd1 = bjd[np.where(bjd < 2457950)]
bjd2 = bjd[np.where(bjd > 2457950)]
lines1 = lines[np.where(bjd < 2457950)]
lines2 = lines[np.where(bjd > 2457950)]
vels1 = vels[np.where(bjd < 2457950)]
vels2 = vels[np.where(bjd > 2457950)]
yerr1 = v_unc[np.where(bjd < 2457950)]
yerr2 = v_unc[np.where(bjd > 2457950)]

plt.figure(figsize = (7,4))
plt.subplot(1,2,1)
plt.errorbar(bjd1, vels1, yerr=yerr1, fmt='r.')
plt.xlabel('BJD')
plt.ylabel('Velocity [km/s]')
plt.subplot(1,2,2)
plt.errorbar(bjd2, vels2, yerr=yerr2, fmt='r.')
plt.xlabel('BJD')
plt.suptitle('CCF centroids')
plt.savefig('figures/ccf_centroids.png')
plt.show()

#%%
#Vi tager et gns. af alle 'lines' der vurderes til at ligge uden for transitten
OT_lines1 = np.vstack((lines1[0:7,:],lines1[-7:,:]))
OT_gns1 = np.mean(OT_lines1,0)

#Vi trækker dette gns. fra alle 'lines'
norm_lines1 = lines1-OT_gns1

# Vi fitter omvendte gaussfunktioner til alle de normerede linjer.
RM_index = 20

RM_vel1 = []
RM_err1 = []
for i in range(len(norm_lines1)):
    line_max = np.where(0-norm_lines1[i] == np.max(0-norm_lines1[i][np.where((vel_vector > -100) & (vel_vector < 100))]))
    norm_fit, norm_pcov = curve_fit(gauss, vel_vector, 0-norm_lines1[i],
                                    bounds=([0.003,-100,1,-0.1],[0.02,100,20,0.1]),
                                    p0=[0.006, vel_vector[line_max], 10, 0])
    if i == RM_index:
        index_fit = norm_fit
    RM_vel1.append(norm_fit[1])
    RM_err1.append(np.sqrt(np.diag(norm_pcov))[1])
RM_vel1 = np.array(RM_vel1)
RM_err1 = np.array(RM_err1)

#Fit Gaussian to out-of-transit mean CCF

OT_fit, OT_pcov = curve_fit(gauss, vel_vector, OT_gns1)
sys_vel = OT_fit[1]
RM_vel1_cor = RM_vel1[9:-8]-sys_vel
bjd1_cor = bjd1[9:-8]
#%%
plt.figure(figsize=(8,5))
plt.subplot(2,1,1)
plt.plot(vel_vector,lines1[0],'-',label='CCF OT')
plt.plot(vel_vector,OT_gns1,'--',label='Mean CCF OT')
plt.plot(vel_vector, gauss(vel_vector, *OT_fit),'--', label = 'Fit to mean')
plt.legend()
plt.xlim([vel_vector[0],vel_vector[-1]])
plt.xticks([],[])
plt.title('Raw and Normalized CCFs')
plt.ylabel('CCF')
plt.subplot(2,1,2)
plt.plot(vel_vector,norm_lines1[20],'-',label='Norm. line IT')
plt.plot(vel_vector,norm_lines1[0],':',label='Norm. line OOT')
plt.xlim([vel_vector[0],vel_vector[-1]])
plt.legend()
plt.xlabel('RV [km/s]')
plt.ylabel('Norm. CCF')
plt.subplots_adjust(hspace=0)
#plt.savefig('figures/ccfs_norm.png')
plt.show()

#%%

locx = np.linspace(0,250,6)
tick_labx = np.round(np.linspace(vel_vector[200],vel_vector[450],6),2)
locy = np.array([0,10,20,30,40])
tick_laby = np.round(bjd1[locy]-2457939,2)

plt.figure(figsize=(8,5))
plt.plot(vel_vector, 0-norm_lines1[20],'-',label='Norm. CCF')
plt.plot(vel_vector, gauss(vel_vector, *index_fit), '-',label='Gaussian fit')
plt.xlabel('RV [km/s]')
plt.ylabel('Norm. CCF')
plt.legend()
plt.title('Fit to inv. norm. CCFs')
plt.tight_layout()
plt.savefig('figures/fitted_ccf.png')
plt.show()
#%%
plt.figure()
plt.errorbar(bjd1_cor, RM_vel1_cor, yerr=RM_err1[9:-8], fmt='.k', label='RM curve, corr.')
#plt.errorbar(bjd1, RM_vel1, yerr=RM_err1, fmt='.', label='RM curve')
plt.plot([bjd1[0]-0.02,bjd1[-1]+0.02],[0,0],'--k',linewidth=1)
plt.xlim([bjd1[0]-0.02,bjd1[-1]+0.02])
plt.ylim([min(RM_vel1)-5,max(RM_vel1)+5])
plt.xlabel('Time [BJD]')
plt.ylabel('Velocities [km/s]')
plt.title('RM Curve, Shifted IT CCFs')
#plt.savefig('figures/RM_shift.png')
plt.show()

mid_t = bjd1_cor[0] + 0.5*(bjd1_cor[-1] - bjd1_cor[0])
bjd1_norm = bjd1_cor-mid_t
#%%
rm_fit, rm_pcov = curve_fit(rm_function_fixed, bjd1_norm, RM_vel1_cor*1000,
                            bounds = ([-70,75,90],[-30,90,130]),
                            p0 = [-45,87,109])
#%%

#%%
plt.figure()
plt.errorbar(bjd1_cor, RM_vel1_cor*1000, yerr=RM_err1[9:-8], fmt='.k', label='RM curve, corr.')
#plt.errorbar(bjd1, RM_vel1*, yerr=RM_err1, fmt='.', label='RM curve')
plt.plot(bjd1_cor,rm_function_fixed(bjd1_norm,*rm_fit),'-', label='Fit')
plt.xlabel('Time [BJD]-2457939')
plt.ylabel('Velocities [km/s]')
plt.tight_layout()
plt.show()

#%%
#plt.close('all')

plt.figure()
for i in range(len(lines1)):
    plt.plot(vel_vector,lines1[i]+0.01*i,'-k')
plt.xlabel('Radial Velocity, [km/s]')
plt.ylabel('CCF')
plt.title('All CCFs')
plt.yticks([],[])
plt.tight_layout()
plt.xlim([vel_vector[0],vel_vector[-1]])
plt.savefig('figures/All_CCFs.png')
plt.show()

plt.figure(figsize=(7,7))
color_map = plt.imshow(np.flipud(norm_lines1[:,200:451]), aspect = 'auto')
color_map.set_cmap('gray')
plt.xticks(locx,tick_labx-round(sys_vel,2))
plt.yticks(locy,tick_laby)
plt.xlabel('RV')
plt.ylabel('Time [BJD]-2457939')
plt.savefig('figures/Colormap.png')
plt.show()




