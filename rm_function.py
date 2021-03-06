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
def rm_function(t, t_p, a, e, m1, m2, rad_rat, obl, incl, omega, RV_amp, R, lim=401, return_option=2):
    
    # Define values
    G = 6.647 * 10**(-11) # m³ kg⁻¹ s⁻²
    # Recalculate parameters(enheder er på omregningsfaktorer)
    m1 = m1 * 2 * 10**30 # kg/m_sun
    m2 = m2 * 2 * 10**30 # kg/m_sun
    a = a * R # R_star
    t = t * 24 * 60**2 # s/days
    t_p = t_p * 24 * 60**2 # s/days
    obl = obl * (np.pi/180) # rad/deg
    omega = omega * (np.pi/180) # rad/deg
    incl = incl * (np.pi/180) # rad/deg
    RV_amp = RV_amp * 1000 # m/km
    R = R
    
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
    
    X_dark = []
    L_dark = []
    for i in range(len(p_x)):
        if p_z[i] >= 0:
            X_planet = X[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
            X_p.append(X_planet)
            Y_planet = Y[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
            Y_p.append(Y_planet)
            L_planet = L[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
            L_p.append(L_planet)
            # Build luminosity behind planet
            X_behind = X[np.where((X-p_x[i])**2+(Y-p_y[i])**2<p_rad**2)]
            X_dark.append(X_behind)
            L_behind = L[np.where((X-p_x[i])**2+(Y-p_y[i])**2<p_rad**2)]
            L_dark.append(L_behind)
        elif p_z[i] < 0:
            X_p.append(X)
            Y_p.append(Y)
            L_p.append(L)
            X_dark.append(X)
            L_dark.append(np.zeros(X.shape))
    X_p = np.array(X_p)
    Y_p = np.array(Y_p)
    L_p = np.array(L_p)
    
    X_dark = np.array(X_dark)
    L_dark = np.array(L_dark)

    # Give every x-value an RV-value
    tail_x = np.linspace(-(1201/401) * R, (1201/401) * R, 1201)
    RV_x = (RV_amp * tail_x/R)
    
    # Create gaussian from microturbulence and instrumental error
    def alt_gauss(x, FWHM):
        return FWHM**(-1) * ((4 * np.log(2))/np.pi)**(1/2) * np.exp(-(4 * np.log(2) * x**2)/(FWHM**2))

    sigma_micro = 2000
    sigma_inst = 2500
    FWHM = np.sqrt(sigma_micro**2 + sigma_inst**2)
    conv_func = alt_gauss(RV_x, FWHM)
    
    # Collapse luminosity along x-axis into a histogram
    L_sum_p = []
    L_sum_dark = []
    for i in range(len(X_p)):
        L_sum = []
        dark_sum = []
        for j in x:
            L_sum.append(sum(L_p[i][np.where(X_p[i]==j)]))
            dark_sum.append(sum(L_dark[i][np.where(X_dark[i]==j)]))
        L_sum = np.array(L_sum)
        dark_sum = np.array(dark_sum)
        # Add zeros to both ends of L_sum
        L_sum = np.hstack((np.zeros(400), L_sum, np.zeros(400)))
        dark_sum = np.hstack((np.zeros(400), dark_sum, np.zeros(400)))
        L_fin = np.convolve(L_sum, conv_func, mode='same')
        dark_fin = np.convolve(dark_sum, conv_func, mode='same')
        L_sum_p.append(L_fin)
        L_sum_dark.append(dark_fin)
    L_sum_p = np.array(L_sum_p)
    L_sum_dark = np.array(L_sum_dark)

    # Fit every list in luminosity object to a gaussian
    def gauss_func(x, a, b, c):
        return a * np.exp(-((x - b)**2)/(2 * c**2))
    gaussians = []
    centroids = []
    centroids_avg = []
    for i in range(len(p_x)):
        f, fs = curve_fit(gauss_func, RV_x, L_sum_p[i], p0=[12, 0, 50])
        gaussians.append(f)
        centroids.append(f[1])
        avg_centroid = np.sum(RV_amp*X_p[i]/R)/len(X_p[i])
        centroids_avg.append(avg_centroid)
    centroids = np.array(centroids)
    centroids_avg = np.array(centroids_avg)
    
    d_gaussians = []
    d_centroids = []
    d_centroids_avg = []
    for i in range(len(p_x)):
        df, dfs = curve_fit(gauss_func, RV_x, L_sum_dark[i])
        d_gaussians.append(df)
        d_centroids.append(df[1])
        d_avg_centroid = np.sum(RV_amp*X_dark[i]/R)/len(X_dark[i])
        d_centroids_avg.append(d_avg_centroid)
    d_centroids = np.array(d_centroids)
    d_centroids_avg = np.array(d_centroids_avg)
        
    if return_option == 0:
        return centroids, centroids_avg
    elif return_option == 1:
        return centroids, centroids_avg, gaussians, X, Y, RV_x, L_sum_p, X_p, Y_p
    elif return_option == 2:
        return d_centroids, d_centroids_avg
    elif return_option == 3:
        return d_centroids, d_centroids_avg, d_gaussians, X, Y, RV_x, L_sum_dark, X_p, Y_p

#%% Function with fixed parameters
# Inputs:
#   time: days
#   semimajor axis: stellar radii
#   mass: solar mass
#   obliquity: degrees
#   inclination: degrees
#   argument of periapsis: degrees
#   RV amplitude: km/s

# Return_option 4 returns array that can be fitted to data

def rm_function_fixed(t, obl, incl, RV_amp, lim=401, return_option=4):
    
    # Define fixed parameters:
    t_p = 0
    a = 9.25
    e = 0
    m1 = 1.72
    m2 = 0.003532
    rad_rat = 0.0735
    omega = 90

    
    
    # Define values
    G = 6.647 * 10**(-11) # m³ kg⁻¹ s⁻²
    R =  695510000 # m
    # Recalculate parameters(enheder er på omregningsfaktorer)
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
    
    X_dark = []
    L_dark = []
    for i in range(len(p_x)):
        if p_z[i] >= 0:
            X_planet = X[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
            X_p.append(X_planet)
            Y_planet = Y[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
            Y_p.append(Y_planet)
            L_planet = L[np.where((X-p_x[i])**2+(Y-p_y[i])**2>p_rad**2)]
            L_p.append(L_planet)
            # Build luminosity behind planet
            X_behind = X[np.where((X-p_x[i])**2+(Y-p_y[i])**2<p_rad**2)]
            X_dark.append(X_behind)
            L_behind = L[np.where((X-p_x[i])**2+(Y-p_y[i])**2<p_rad**2)]
            L_dark.append(L_behind)
        elif p_z[i] < 0:
            X_p.append(X)
            Y_p.append(Y)
            L_p.append(L)
            X_dark.append(X)
            L_dark.append(np.zeros(X.shape))
    X_p = np.array(X_p)
    Y_p = np.array(Y_p)
    L_p = np.array(L_p)
    
    X_dark = np.array(X_dark)
    L_dark = np.array(L_dark)

    # Give every x-value an RV-value
    tail_x = np.linspace(-(1201/401) * R, (1201/401) * R, 1201)
    RV_x = (RV_amp * tail_x/R)
    
    # Create gaussian from microturbulence and instrumental error
    def alt_gauss(x, FWHM):
        return FWHM**(-1) * ((4 * np.log(2))/np.pi)**(1/2) * np.exp(-(4 * np.log(2) * x**2)/(FWHM**2))

    sigma_micro = 2000
    sigma_inst = 2500
    FWHM = np.sqrt(sigma_micro**2 + sigma_inst**2)
    conv_func = alt_gauss(RV_x, FWHM)
    
    # Collapse luminosity along x-axis into a histogram
    L_sum_p = []
    L_sum_dark = []
    for i in range(len(X_p)):
        L_sum = []
        dark_sum = []
        for j in x:
            L_sum.append(sum(L_p[i][np.where(X_p[i]==j)]))
            dark_sum.append(sum(L_dark[i][np.where(X_dark[i]==j)]))
        L_sum = np.array(L_sum)
        dark_sum = np.array(dark_sum)
        # Add zeros to both ends of L_sum
        L_sum = np.hstack((np.zeros(400), L_sum, np.zeros(400)))
        dark_sum = np.hstack((np.zeros(400), dark_sum, np.zeros(400)))
        L_fin = np.convolve(L_sum, conv_func, mode='same')
        dark_fin = np.convolve(dark_sum, conv_func, mode='same')
        L_sum_p.append(L_fin)
        L_sum_dark.append(dark_fin)
    L_sum_p = np.array(L_sum_p)
    L_sum_dark = np.array(L_sum_dark)

    # Fit every list in luminosity object to a gaussian
    def gauss_func(x, a, b, c):
        return a * np.exp(-((x - b)**2)/(2 * c**2))
    gaussians = []
    centroids = []
    centroids_avg = []
    for i in range(len(p_x)):
        f, fs = curve_fit(gauss_func, RV_x, L_sum_p[i], p0=[12, 0, 50])
        gaussians.append(f)
        centroids.append(f[1])
        avg_centroid = np.sum(RV_amp*X_p[i]/R)/len(X_p[i])
        centroids_avg.append(avg_centroid)
    centroids = np.array(centroids)
    centroids_avg = np.array(centroids_avg)
    
    d_gaussians = []
    d_centroids = []
    d_centroids_avg = []
    for i in range(len(p_x)):
        df, dfs = curve_fit(gauss_func, RV_x, L_sum_dark[i])
        d_gaussians.append(df)
        d_centroids.append(df[1])
        d_avg_centroid = np.sum(RV_amp*X_dark[i]/R)/len(X_dark[i])
        d_centroids_avg.append(d_avg_centroid)
    d_centroids = np.array(d_centroids)
    d_centroids_avg = np.array(d_centroids_avg)
    
    d_centroids_avg[np.where(np.isnan(d_centroids_avg))]=0
    
    if return_option == 0:
        return centroids, centroids_avg
    elif return_option == 1:
        return centroids, centroids_avg, gaussians, X, Y, RV_x, L_sum_p, X_p, Y_p
    elif return_option == 2:
        return d_centroids, d_centroids_avg
    elif return_option == 3:
        return d_centroids, d_centroids_avg, d_gaussians, X, Y, RV_x, L_sum_dark, X_p, Y_p
    elif return_option == 4:
        return d_centroids_avg
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