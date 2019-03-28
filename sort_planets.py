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

depth_thres = 1
# =============================================================================
# RA_thres = 
# DEC_thres = 
# mag_thres = 
# =============================================================================

nice_planets = []

for i in range(len(lines)):
    if lines[i][0] == '1':
        nice_planets.append(lines[i])

