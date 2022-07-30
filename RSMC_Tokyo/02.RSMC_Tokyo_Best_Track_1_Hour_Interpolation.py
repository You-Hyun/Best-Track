#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 09:03:02 2022

@author: yhbaek
"""

import os 
import numpy as np
from scipy import interpolate 

# %% 

Main_DIR=('/home/yhbaek/DATA/Best-Track/RSMC_Tokyo/')
Rearrangement_DIR=(Main_DIR  + 'Data/02.Rearrangement/')
Interpolation_DIR=(Main_DIR  + 'Data/03.Interpolation/')

if not (os.path.exists(Interpolation_DIR)):
    os.makedirs(Interpolation_DIR)

# %%
for Year in np.arange(1950, 2023):
    for Num in np.arange(1, 100):
        TC_Num_file=(Rearrangement_DIR + str(Year) + '/' + str(Year) + '_' + str(Num).zfill(2) + '.npy')
        
        if (os.path.exists(TC_Num_file)): 
            Track=np.load(TC_Num_file)
            
# %% 1-h interpolation
            
        
    del(Num)
del(Year)        
        





# %% 1-h interpolation

x=np.where(np.isnan(Track[N, :, L])==False)[0]

if len(x)>=2:
    y=Track[N, x, L]
    f=interpolate.interp1d(x, y, kind='linear')
    
    xnew=np.arange(x[0], x[len(x)-1]+1)
    ynew=f(xnew)
    
    Track[N, x[0]:x[len(x)-1]+1, L]=ynew
    
    del(y, f, xnew, ynew)
del(x)
