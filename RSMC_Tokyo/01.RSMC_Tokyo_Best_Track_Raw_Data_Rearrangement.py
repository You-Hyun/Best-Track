#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 17:08:11 2022

@author: yhbaek
"""


import os 
import numpy as np

# %% 

Main_DIR=('/home/yhbaek/DATA/Best-Track/RSMC_Tokyo/')
Raw_DiR=(Main_DIR + 'Data/01.Raw/')
Rearrangement_DIR=(Main_DIR  + 'Data/02.Rearrangement/')

if not (os.path.exists(Rearrangement_DIR)):
    os.makedirs(Rearrangement_DIR)

# %% Read Header

track_file=open(Raw_DiR + 'bst_all_20220721.txt', 'r')


for N in range(2000):
    
    header=track_file.readline()
    
    if len(header)>0:
        Indicator=int(header[:5])
        
        if int(header[6:8])>=51:
            Year=int('19' + header[6:8])
        else:
            Year=int('20' + header[6:8])
            
        NumberID=int(header[8:10])
        Dataline=int(header[12:15])        
        print('N: ', N, ' Year: ', str(Year), ' TC ID: ', str(NumberID).zfill(2))
# %%
    
        if not (os.path.exists(Rearrangement_DIR + str(Year))):
            os.makedirs(Rearrangement_DIR + str(Year))
    
        Track=np.zeros((Dataline, 14))
            
        for L in range(Dataline):
                   
            line=track_file.readline()
            
            if int(line[:2])>=51:
                Track[L, 0]=int('19' + line[:2])                        #Year
            else:
                Track[L, 0]=int('20' + line[:2])                        #Year
                
            Track[L, 1]=int(line[2:4])                                  #Month
            Track[L, 2]=int(line[4:6])                                  #Day
            Track[L, 3]=int(line[6:8])                                  #Hour
            Track[L, 4]=int(line[15:18])/10                             #Lat
            Track[L, 5]=int(line[19:23])/10                             #Lon
            Track[L, 6]=int(line[24:28])                                #Press
            
            if (line[33:36]=='   ')==True:
                Track[L, 7]=np.NaN
            else:
                Track[L, 7]=int(line[33:36])                            #Int
           
            if (line[41:42]==' ')==True:
                Track[L, 8]=np.NaN
            else:
                Track[L, 8]=int(line[41:42])                            #R50_dir
                
            if (line[42:46]=='    ')==True:
                Track[L, 9]=np.NaN
            else:
                Track[L, 9]=int(line[42:46])                            #R50_long
                
            if (line[47:51]=='    ')==True:
                Track[L, 10]=np.NaN
            else:
                Track[L, 10]=int(line[47:51])                          #R50_Short               
                
            if (line[52:53]==' ')==True:
                Track[L, 11]=np.NaN
            else:
                Track[L, 11]=int(line[52:53])                           #R30_dir   
             
            if (line[53:57]=='    ')==True:
                Track[L, 12]=np.NaN
            else:
                Track[L, 12]=int(line[53:57])                           #R30_long
                
            if (line[58:62]=='    ')==True:
                Track[L, 13]=np.NaN
            else:
                Track[L, 13]=int(line[58:62])                           #R30_Short    
                
                
            del(line)
        del(L)   
        
        
        np.save(Rearrangement_DIR + str(Year) + '/' + str(Year) + '_' + str(NumberID).zfill(2) + '.npy', Track)
        del(Track)
        del(Indicator, NumberID, Dataline, Year)
        
del(N)    
del(header)

track_file.close()
del(track_file)
# %%
