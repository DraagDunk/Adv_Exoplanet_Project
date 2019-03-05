#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:25:32 2019

@author: jesper
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.close('all')

# Give stjerneparametre
R = 1
r, lum = simple_limb(R)

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
b = R**2-(X**2 + Y**2)
            
# Plot stjerne
plt.figure()
plt.scatter(X, Y, c=b, cmap=cm.gray)
plt.xlabel(r'x [$R_\odot$]')
plt.ylabel(r'y [$R_\odot$]')
plt.axis('equal')
plt.tight_layout()
plt.show()

# Kollaps y-akse
plt.figure()
plt.hist(X, bins=100)
plt.xlabel(r'x [$R_\odot$]')
plt.tight_layout()
plt.show()

#%% TODO liste

# Lav model for transit
# Kollaps stjerne for farver
# Ellers andet?