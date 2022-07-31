#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 09:03:02 2022

@author: yhbaek
"""

import os 
import numpy as np
import pandas as pd
from scipy import interpolate 

# %% 

Main_DIR=('/Users/yhbaek/STUDY/2022/Best_Track/RSMC_Tokyo/')
Rearrangement_DIR=(Main_DIR  + 'Data/02.Rearrangement/')
Interpolation_DIR=(Main_DIR  + 'Data/03.Interpolation/')

if not (os.path.exists(Interpolation_DIR)):
    os.makedirs(Interpolation_DIR)

# %%
N=1
for Year in np.arange(1950, 2023):
    if not (os.path.exists(Interpolation_DIR + str(Year))):
        os.makedirs(Interpolation_DIR + str(Year))
            
    for Num in np.arange(1, 100):
        TC_Num_file=(Rearrangement_DIR + str(Year) + '/' + str(Year) + '_' + str(Num).zfill(2) + '.npy')
        
        if (os.path.exists(TC_Num_file)): 
            Track=np.load(TC_Num_file)

# %% Convert Gregorian date to Julian date
        
            str_date=pd.Timestamp(year=int(Track[0, 0]), 
                                  month = int(Track[0, 1]), 
                                  day = int(Track[0, 2]), 
                                  hour = int(Track[0, 3])).to_julian_date()
            end_date=pd.Timestamp(year=int(Track[-1, 0]), 
                                  month = int(Track[-1, 1]), 
                                  day = int(Track[-1, 2]), 
                                  hour = int(Track[-1, 3])).to_julian_date()   
                            
            Track_interp=np.zeros((len(np.arange(str_date, end_date+(0.041666667), 0.041666667)), 14))   # 0.041666667 = 1hour
            
            count=0
            for L in np.arange(str_date, end_date+(0.041666667), 0.041666667):                          # 0.041666667 = 1hour
                
                Track_interp[count, 0]=pd.to_datetime(L, unit = 'D', origin = 'julian').year
                Track_interp[count, 1]=pd.to_datetime(L, unit = 'D', origin = 'julian').month
                Track_interp[count, 2]=pd.to_datetime(L, unit = 'D', origin = 'julian').day
                Track_interp[count, 3]=pd.to_datetime(L, unit = 'D', origin = 'julian').hour
                
                Same_Time=np.where((Track[:, 0]==Track_interp[count, 0]) & 
                                   (Track[:, 1]==Track_interp[count, 1]) &
                                   (Track[:, 2]==Track_interp[count, 2]) &
                                   (Track[:, 3]==Track_interp[count, 3]))[0]
                if len(Same_Time)>0:
                    Track_interp[count, 4:]=Track[Same_Time, 4:]
                else:
                    Track_interp[count, 4:]=np.NaN
                
                count=count+1
                del(Same_Time)
            del(L)
            del(count)
            del(str_date, end_date)
            
            Track_interp[np.where(Track_interp[:, 4:]==0)[0], np.where(Track_interp[:, 4:]==0)[1]+4]=np.NaN
            
    # %% 1-h interpolation
            
            for Variable in np.arange(4, Track.shape[1], 1):
                x=np.where(np.isnan(Track_interp[:, Variable])==False)[0]
            
                if len(x)>=2:
                    y=Track_interp[x, Variable]
                    f=interpolate.interp1d(x, y, kind='linear')
                    
                    xnew=np.arange(x[0], x[len(x)-1]+1)
                    ynew=f(xnew)
                    
                    Track_interp[np.arange(x[0], x[len(x)-1]+1), Variable]=ynew
                
                    del(y, f, xnew, ynew)
                del(x)
            del(Variable)   
        
            np.save(Interpolation_DIR + str(Year) + '/' + str(Year) + '_' + str(Num).zfill(2) + '.npy', Track_interp)
            print('1-Hour Interp. - N: ', N, ' Year: ', str(Year), ' TC ID: ', str(Num).zfill(2))
            N=N+1
            del(Track, Track_interp)     
        del(TC_Num_file)
    del(Num)
del(Year)    
del(N)    
        

