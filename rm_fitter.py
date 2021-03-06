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
                            bounds = ([-70,75,90],[-20,90,130]),
                            p0 = [-45,87,109])
#%%

#%%
t_linspace = np.linspace(bjd1_cor[0]-0.03, bjd1_cor[-1]+0.03,200)
tzero_linspace = np.linspace(bjd1_norm[0]-0.03,bjd1_norm[-1]+0.03,200)

plt.figure()
plt.errorbar(bjd1_cor, RM_vel1_cor, yerr=RM_err1[9:-8], fmt='.k', label='RM curve, corr.')
#plt.errorbar(bjd1, RM_vel1*, yerr=RM_err1, fmt='.', label='RM curve')
plt.plot(t_linspace,rm_function_fixed(tzero_linspace,*rm_fit)/1000,'-', label='Fit')
plt.title('Data fitted to three parameters')
plt.xlabel('Time [BJD]')
plt.ylabel('Velocities [km/s]')
plt.xlim([bjd1_cor[0]-0.02,bjd1_cor[-1]+0.02])
plt.legend()
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


#%%
per = 2.148780  #period in days
t0_teo = 2457097.278    #T0 from archive
n = int((mid_t-t0_teo)/per)
t0 = t0_teo + n*per
dur = 0.169 #Transit duration in days
start_t = t0-dur/2
end_t = t0+dur/2
mid_t_index = np.where(abs(bjd1-t0)==min(abs(bjd1-t0)))[0][0]
start_index = np.where(abs(bjd1-start_t)==min(abs(bjd1-start_t)))[0][0]
end_index = np.where(abs(bjd1-end_t)==min(abs(bjd1-end_t)))[0][0]

dur_full = (per/np.pi)*np.arcsin((0.00977/0.043)*np.sqrt((1+(0.15/2.1)**2)+0.23**2)/np.sin(np.deg2rad(87)))
start_full = t0-dur_full/2
end_full = t0+dur_full/2
start_full_index = np.where(abs(bjd1-start_full)==min(abs(bjd1-start_full)))[0][0]
end_full_index = np.where(abs(bjd1-end_full)==min(abs(bjd1-end_full)))[0][0]

vel_cor = vel_vector[100:551]-sys_vel
rm_max_index = np.where(abs(vel_cor-109)==min(abs(vel_cor-109)))[0][0]
rm_min_index = np.where(abs(vel_cor+109)==min(abs(vel_cor+109)))[0][0]
rm_mid_index = np.where(abs(vel_cor)==min(abs(vel_cor)))[0][0]

locx = np.linspace(0,450,6)
tick_labx = np.round(np.linspace(vel_cor[0],vel_vector[-1],6),2)
locy = np.array([0,10,20,30,40])
tick_laby = np.round(bjd1[locy]-2457939,2)

plt.figure(figsize=(7,7))
color_map = plt.imshow(np.flipud(norm_lines1[:,100:551]), aspect = 'auto')
color_map.set_cmap('gray')
plt.plot([0,len(norm_lines1[0,100:550])],[mid_t_index,mid_t_index],'--r')
plt.plot([0,len(norm_lines1[0,100:550])],[start_index,start_index],'-k')
plt.plot([0,len(norm_lines1[0,100:550])],[end_index,end_index],'-k')
plt.plot([0,len(norm_lines1[0,100:550])],[start_full_index,start_full_index],'--k')
plt.plot([0,len(norm_lines1[0,100:550])],[end_full_index,end_full_index],'--k')
plt.plot([rm_max_index,rm_max_index],[0,len(norm_lines1)-1],'-r')
plt.plot([rm_min_index,rm_min_index],[0,len(norm_lines1)-1],'-r')
plt.plot([rm_mid_index,rm_mid_index],[0,len(norm_lines1)-1],'--r')
plt.xticks(locx,tick_labx-round(sys_vel,2))
plt.yticks(locy,tick_laby)
plt.xlabel('RV')
plt.ylabel('Time [BJD]-2457939')
plt.savefig('figures/Colormap.png')
plt.show()




