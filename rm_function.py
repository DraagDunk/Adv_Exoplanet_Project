#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:21:10 2019

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Inputs:
#   time: days
#   semimajor axis: stellar radii
#   mass: solar mass
#   obliquity: degrees
#   inclination: degrees
#   argument of periapsis: degrees
#   RV amplitude: km/s
def rm_function(t, t_p, a, e, m1, m2, rad_rat, obl, incl, omega, RV_amp, lim=401, return_grid=False):
    
    # Define values
    G = 6.647 * 10**(-11) # m³ kg⁻¹ s⁻²
    R = 695510000 # m
    # Recalculate parameters
    m1 = m1 * 2 * 10**30 # kg/m_sun
    m2 = m2 * 2 * 10**30 # kg/m_sun
    a = a * R # R_star
    t = t * 24 * 60**2 # s/days
    t_p = t_p * 24 * 60**2 # s/days
    obl = obl * (np.pi/180) # rad/deg
    omega = omega * (np.pi/180) # rad/deg
    incl = incl * (np.pi/180) # rad/deg
    RV_amp = RV_amp * 1000 # m/km
    
    # Calculate mean anomaly
    G_mu = G * m1 + G * m2
    Mt = np.sqrt((G_mu)/(a**3)) * (t - t_p)
    
    # Calculate eccentric anomaly
    E_0 = np.copy(Mt)
    E_prev = np.copy(E_0)
    E = np.zeros(len(E_prev))
    for i in range(len(E)):
        E[i] = E_prev[i] - ((E_prev[i] - e * np.sin(E_prev[i]) - Mt[i]) / (1 - e * np.sin(E_prev[i])))
        while abs(E[i] - E_prev[i]) > 10**(-8):
            E_prev[i] = np.copy(E[i])
            E[i] = E_prev[i] - ((E_prev[i] - e * np.sin(E_prev[i]) - Mt[i]) / (1 - e * np.sin(E_prev[i])))
    
    # Calculate true anomaly
    nu = 2 * np.arctan( ((1+e)/(1-e))**(1/2) * np.tan(E/2) )
    
    # Calculate separation
    r = a*(1-e**2)/(1+e*np.cos(nu))
    
    # Build vectors of x-, y- and z-coordinates of the planet
    init_p_x = -r*np.cos(omega + nu)
    init_p_y = -r*np.sin(omega + nu)*np.cos(incl)
    p_z = r*np.sin(omega + nu)*np.sin(incl)

    p_x = init_p_x * np.cos(obl) + init_p_y * np.sin(obl)
    p_y = - init_p_x * np.sin(obl) + init_p_y * np.cos(obl)
    
    # BUILD STAR

    # Define work area
    lim = lim
    x = np.linspace(-R, R, lim)
    y = np.linspace(-R, R, lim)

    # Build star
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
    
    # Calculate planet radius
    p_rad = R * rad_rat
    
    # Build X and Y of the star for each time-step of the planet
    X_p = []
    Y_p = []
    L_p = []
    for i in range(len(p_x)):
        if p_z[i] >= 0:
            X_planet = X[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
            X_p.append(X_planet)
            Y_planet = Y[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
            Y_p.append(Y_planet)
            L_planet = L[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
            L_p.append(L_planet)
        elif p_z[i] < 0:
            X_p.append(X)
            Y_p.append(Y)
            L_p.append(L)
    X_p = np.array(X_p)
    Y_p = np.array(Y_p)
    L_p = np.array(L_p)

    # Collapse luminosity along x-axis into a histogram
    L_sum_p = []
    for i in range(len(X_p)):
        L_sum = []
        for j in x:
            L_sum.append(sum(L_p[i][np.where(X_p[i]==j)]))
        L_sum = np.array(L_sum)
        L_sum_p.append(L_sum)
    L_sum_p = np.array(L_sum_p)

    # Give every x-value an RV-value
    RV_x = RV_amp * x/R

    # Fit every list in luminosity object to a gaussian
    def gauss_func(x, a, b, c):
        return a * np.exp(-((x - b)**2)/(2 * c**2))
    gaussians = []
    centroids = []
    for i in range(len(p_x)):
        f, fs = curve_fit(gauss_func, RV_x, L_sum_p[i])
        gaussians.append(f)
        centroids.append(f[1])
    
    if return_grid == False:
        return centroids
    elif return_grid == True:
        return centroids, gaussians, X, Y, RV_x, L_sum_p, X_p, Y_p

#t = np.linspace(-1,1, 200)
## t, t_p, a, e, m1, m2, rad_rat, obl, i, omega, RV_amp
#RV = rm_function(t, 0, 2, 0.1, 1, 0.001, 0.5, 90, 0, 90, 0, 40)
#
#plt.figure()
#plt.plot(t, RV, 'k-')
#plt.xlabel('time [days]')
#plt.ylabel('RV [km/s]')
#plt.tight_layout()
#plt.show()