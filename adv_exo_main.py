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

plt.close('all')

# Give stjerneparametre
R = 1
#r, lum = simple_limb(R)

# Definer arbejdsområde
x = np.linspace(-1, 1, 501)
y = np.linspace(-1, 1, 501)

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
L = R**2-(X**2 + Y**2)

p_rad = rnd.randint(2, 50)/100
p_b = rnd.randint(0,99)/100
p_o = rnd.randint(-180, 180)
svar = [p_rad, p_o, p_b]
alpha = np.tan((p_o*np.pi)/180)
beta = p_b/np.cos((p_o*np.pi)/180)*(1+p_rad)
p_x = np.linspace( -np.sign(beta) * (1 + p_rad), np.sign(beta) * (1 + p_rad), 100)
move_vector = []

p_y = np.array([])
for i in p_x:
    p_y = np.append(p_y, alpha*i+beta)
            
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

L_sum_p = []
for i in range(len(X_p)):
    L_sum = []
    for j in x:
        L_sum.append(sum(L_p[i][np.where(X_p[i]==j)]))
    L_sum = np.array(L_sum)
    L_sum_p.append(L_sum)
    
L_sum_p = np.array(L_sum_p)

# Kollaps y-akse
fig2 = plt.figure()
ax3 = fig2.add_subplot(111)
line2, = ax3.plot(x, L_sum_p[0], 'r-')
plt.xlabel(r'x [$R_\odot$]')
plt.ylabel('Luminosity [unit]')
plt.ylim([0,350])
plt.xlim([-1, 1])
plt.tight_layout()
plt.show()

for i in range(len(X_p)):
    line1.set_xdata(X_p[i])
    line1.set_ydata(Y_p[i])
    fig1.canvas.draw()
    fig1.canvas.flush_events()
    line2.set_ydata(L_sum_p[i])
    fig2.canvas.draw()
    fig2.canvas.flush_events()

#%% TODO liste

# Implement color/rotation
# Create RM diagram
# Implement as function of obliquity, rotation, planet radius & impact parameter
# Allersidst: Animér RM diagram