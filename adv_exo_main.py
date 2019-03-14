#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:25:32 2019

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random as rnd
from scipy.optimize import curve_fit

plt.close('all')

# Give stjerneparametre
R = 1
#r, lum = simple_limb(R)

# Definer arbejdsområde
lim = 401
x = np.linspace(-1, 1, lim)
y = np.linspace(-1, 1, lim)

# Lav stjerne
X = []
Y = []
for i in range(len(x)):
    for j in range(len(y)):
        if x[i]**2 + y[j]**2 <= R**2:
            X.append(x[i])
            Y.append(y[j])
            
X = np.array(X)
Y = np.array(Y)

# Limb darkening function:
mu = np.sqrt(1 - X**2 - Y**2)
mu[np.where(np.isnan(mu)==True)] = 0
u1 = 0.6
u2 = 0
L = 1 - u1 * (1 - mu) - u2 * (1 - mu)**2

# Definér planetparametre 
p_rad = rnd.randint(2, 50)*R/100
p_b = rnd.randint(0,99)/100
p_o = rnd.randint(-180, 180)
p_a = rnd.randint(20,1000)*R/10
p_e = 0
#rnd.randint(0,90)/100
p_i = np.pi/2
p_omega = np.pi/2
nu_max = np.arcsin(R*(1.1+p_rad/R)/p_a)
p_nu = np.linspace(-nu_max, nu_max, 101)
r_orb = p_a*(1-p_e**2)/(1+p_e*np.cos(p_nu))

# Byg vektorer af x- og -ykoordinater af planet
init_p_x = -r_orb*np.cos(p_omega + p_nu)
init_p_y = -r_orb*np.sin(p_omega + p_nu)*np.cos(p_i) + p_b

p_x = init_p_x * np.cos(p_o*(np.pi/180)) + init_p_y * np.sin(p_o*(np.pi/180))
p_y = - init_p_x * np.sin(p_o*(np.pi/180)) + init_p_y * np.cos(p_o*(np.pi/180))

# Plot stjerne
fig1 = plt.figure()
ax1 = plt.gca()
ax1.set_facecolor(('black'))
ax2 = fig1.add_subplot(111)
#line1 = ax2.scatter(X, Y, c=L, cmap=cm.gray)
line1, = ax2.plot(X, Y, 'r,')
plt.xlabel(r'x [$R_\odot$]')
plt.ylabel(r'y [$R_\odot$]')
plt.axis('equal')
plt.tight_layout()
plt.show()

## Plot stjerne med colormap
#figx = plt.figure()
#axx = plt.gca()
#axx.set_facecolor(('black'))
#plt.scatter(X, Y, marker=',', c=L, cmap=cm.gray)
##line1, = ax2.plot(X, Y, 'r,')
#plt.xlabel(r'x [$R_\odot$]')
#plt.ylabel(r'y [$R_\odot$]')
#plt.axis('equal')
#plt.tight_layout()
#plt.show()

# Byg X og Y værdier af stjernen for hvert bevægelsesskridt af planeten
X_p = []
Y_p = []
L_p = []
for i in range(len(p_x)):
    X_planet = X[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
    X_p.append(X_planet)
    Y_planet = Y[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
    Y_p.append(Y_planet)
    L_planet = L[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
    L_p.append(L_planet)
X_p = np.array(X_p)
Y_p = np.array(Y_p)
L_p = np.array(L_p)

# Kollaps lysstyrke fra stjernen ned på x-aksen i histogram.
L_sum_p = []
for i in range(len(X_p)):
    L_sum = []
    for j in x:
        L_sum.append(sum(L_p[i][np.where(X_p[i]==j)]))
    L_sum = np.array(L_sum)
    L_sum_p.append(L_sum)
L_sum_p = np.array(L_sum_p)

# Giv hver x-værdi en RV-hastighed.
RV_amp = 2 # km/s
RV_x = RV_amp * x/R

# For hver indgang i luminositetsobjektet, fit til en gausskurve
def gauss_func(x, a, b, c):
    return a * np.exp(-((x - b)**2)/(2 * c**2))
gaussians = []
centroids = []
for i in range(len(p_x)):
    f, fs = curve_fit(gauss_func, RV_x, L_sum_p[i])
    gaussians.append(f)
    centroids.append(f[1])

# Kollaps y-akse
fig2 = plt.figure()
ax3 = fig2.add_subplot(111)
line2, = ax3.plot(RV_x, L_sum_p[0], 'k-')
line3, = ax3.plot(RV_x, gauss_func(RV_x, *gaussians[0]), 'r-')
plt.xlabel('RV [km/s]')
plt.ylabel('Luminosity [unit]')
plt.ylim([0,1.1*max(L_sum_p[0])])
plt.xlim([-RV_amp, RV_amp])
plt.tight_layout()
plt.show()

# Opdatér plots
for i in range(len(X_p)):
    line1.set_xdata(X_p[i])
    line1.set_ydata(Y_p[i])
    fig1.canvas.draw()
    fig1.canvas.flush_events()
    line2.set_ydata(L_sum_p[i])
    fig2.canvas.draw()
    line3.set_ydata(gauss_func(RV_x, *gaussians[i]))
    fig2.canvas.flush_events()
    
# Plot Rossiter-McLaughlin kurve
plt.figure()
plt.plot(range(len(p_x)), centroids, 'k-')
plt.title('Rossiter-McLaughlin curve')
plt.xlabel('Step')
plt.ylabel('RV [km/s]')
plt.tight_layout()
plt.show()

#%% TODO liste

# Implement color/rotation
# Create RM diagram
# Implement as function of obliquity, rotation, planet radius & impact parameter
# Allersidst: Animér RM diagram