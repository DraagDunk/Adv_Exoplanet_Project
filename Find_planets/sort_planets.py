#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 12:00:42 2019

@author: dina
"""
import numpy as np
import matplotlib.pyplot as plt

planets = open('planets.txt','r')
lines = planets.readlines()
for i in range(len(lines)):
    lines[i] = lines[i].split(',')

sun_rad = 695.51*10**6 #[m]
jup_rad = 69.911*10**6 #[m]

RA_max = 21.0*15
RA_min = 11.0*15
DEC_min = -10.0
DEC_max = 90.0
mag_thres = 14.0
dep_thres = 0.001
vsini_thres = 0.0

tran_planets = []
for i in range(len(lines)):
    if lines[i][0] == '1':
        tran_planets.append(lines[i])

nice_planets = []
for i in range(len(tran_planets)):
    RA = float(tran_planets[i][8])
    DEC = float(tran_planets[i][9])
    mag_str = tran_planets[i][10]
    srad_str = tran_planets[i][13]
    prad_str = tran_planets[i][23].strip('\n')
    T0_str = tran_planets[i][19]
    per_str = tran_planets[i][4]
    tdur_str = tran_planets[i][18]
    vsini_str = tran_planets[i][21]
    if RA < RA_max and RA > RA_min and DEC < DEC_max and DEC > DEC_min and mag_str != '' and srad_str != '' and prad_str != '' and T0_str != '' and per_str != '' and tdur_str != '' and vsini_str != '':
        mag = float(mag_str)
        srad = float(srad_str)
        prad = float(prad_str)
        tran_dep = (prad*jup_rad/(srad*sun_rad))**2
        vsini = float(vsini_str)
        if mag < mag_thres and tran_dep >= dep_thres and vsini > vsini_thres:
            nice_planets.append(tran_planets[i])
    
T0s = np.zeros(len(nice_planets))
periods = np.zeros(len(nice_planets))
stars_file = open('file_input_STARS','w+')
bjd_file = open('file_input_BJD','w+')
for i in range(len(nice_planets)):
    name = nice_planets[i][2].replace(' ','')
    tdur = float(nice_planets[i][18])*24
    RA_h = float(nice_planets[i][8])/15
    DEC_deg = nice_planets[i][9]
    mag_st = nice_planets[i][10]
    fil = '1 ' + name + ' ' + str(tdur) + ' ' + str(RA_h) + ' ' + DEC_deg + ' ' + mag_st + ' ' + name + '.BJD'
    
    dat_file = open(name + '.dat','w+')
    dat_file.write(fil)
    
    T0s[i] = float(nice_planets[i][19])
    periods[i] = float(nice_planets[i][4])
    
    stars_file.write(name + '.dat\n')
    bjd_file.write(name + '.BJD\n')

T0_fil = open('t0s.txt','w+')
per_fil = open('periods.txt','w+')
for i in range(len(nice_planets)):
    T0_fil.write(str(T0s[i]) + '\n')
    per_fil.write(str(periods[i]) + '\n')
