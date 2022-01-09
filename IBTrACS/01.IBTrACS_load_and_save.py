#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 14:49:51 2021

@author: yhbaek
"""

import os 
import numpy as np
from scipy import interpolate 
import julian
from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import wget
# %% Data load
DIR=('/home/yhbaek/STUDY/DATA/Best-Track/IBTrACS/')

IBTrACS_path='https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/netcdf/IBTrACS.ALL.v04r00.nc'

if os.path.isfile(DIR + 'IBTrACS.ALL.v04r00.nc'):
    os.remove(DIR + 'IBTrACS.ALL.v04r00.nc')

wget.download(IBTrACS_path, DIR)

track_file = Dataset(DIR + 'IBTrACS.ALL.v04r00.nc',"r")
#print(track_file.variables.keys())

src='usa'
source = np.array(track_file.variables['source_' + src][:])
basin = np.array(track_file.variables['basin'][:])
time = np.array(track_file.variables['time'][:])
sid = np.array(track_file.variables['sid'][:])
dist2land = np.array(track_file.variables['dist2land'][:], dtype=np.float32)

r34 = np.array(track_file.variables[src + '_r34'][:], dtype=np.float32)
r50 = np.array(track_file.variables[src + '_r50'][:], dtype=np.float32)
r64 = np.array(track_file.variables[src + '_r64'][:], dtype=np.float32)
roci = np.array(track_file.variables[src + '_roci'][:], dtype=np.float32)
poci = np.array(track_file.variables[src + '_poci'][:], dtype=np.float32)
rmw = np.array(track_file.variables[src + '_rmw'][:], dtype=np.float32)
mws = np.array(track_file.variables[src + '_wind'][:], dtype=np.float32)
pres = np.array(track_file.variables[src + '_pres'][:], dtype=np.float32)
eye = np.array(track_file.variables[src + '_eye'][:], dtype=np.float32)

lon = np.array(track_file.variables[src + '_lon'][:], dtype=np.float32)
lat = np.array(track_file.variables[src + '_lat'][:], dtype=np.float32)

del(src)

time[np.where(time==-9.999e+06)]=np.NaN
dist2land[np.where(dist2land==-9999)]=np.NaN

r34[np.where(r34==-9999)]=np.NaN
r50[np.where(r50==-9999)]=np.NaN
r64[np.where(r64==-9999)]=np.NaN
rmw[np.where(rmw==-9999)]=np.NaN
mws[np.where(mws==-9999)]=np.NaN
roci[np.where(roci==-9999)]=np.NaN
poci[np.where(poci==-9999)]=np.NaN
pres[np.where(pres==-9999)]=np.NaN
eye[np.where(eye==-9999)]=np.NaN

lon[np.where(lon==-9999)]=np.NaN
lat[np.where(lat==-9999)]=np.NaN

del(track_file)


# %% Find WP region
find_basin=np.zeros(len(time[:, 0]))        
for N in range(len(time[:, 0])):
    if len(np.where( (basin[N, :, 0]==b'W') & (basin[N, :, 1]==b'P'))[0])>0: # WP
        find_basin[N]=1
del(N)        
    
time=time[np.where(find_basin==1)[0], :]
dist2land=dist2land[np.where(find_basin==1)[0], :]

lon=lon[np.where(find_basin==1)[0], :]
lat=lat[np.where(find_basin==1)[0], :]
mws=mws[np.where(find_basin==1)[0], :]
pres=pres[np.where(find_basin==1)[0], :]
rmw=rmw[np.where(find_basin==1)[0], :]
roci=roci[np.where(find_basin==1)[0], :]
poci=poci[np.where(find_basin==1)[0], :]
eye=eye[np.where(find_basin==1)[0], :]

r34=r34[np.where(find_basin==1)[0], :, :]
r50=r50[np.where(find_basin==1)[0], :, :]
r64=r64[np.where(find_basin==1)[0], :, :]


source=source[np.where(find_basin==1)[0], :]
basin=basin[np.where(find_basin==1)[0], :]
sid=sid[np.where(find_basin==1)[0], :]

del(find_basin)
# %% Find source
find_source=np.zeros(len(time[:, 0]))        
for N in range(len(time[:, 0])):
    if len(np.where( (source[N, 0]==b'b') & (source[N, 1]==b'w') & (source[N, 2]==b'p'))[0])>0: # bwp
        find_source[N]=1
    elif len(np.where( (source[N, 0]==b'B') & (source[N, 1]==b'O') & (source[N, 2]==b'T'))[0])>0: # tcvitals
        find_source[N]=2       
del(N)        


#agency=1     # 1=JTWC, 2=tcvitals

time=time[np.where( (find_source==1) | (find_source==2) )[0], :]
dist2land=dist2land[np.where( (find_source==1) | (find_source==2) )[0], :]

lon=lon[np.where( (find_source==1) | (find_source==2) )[0], :]
lat=lat[np.where( (find_source==1) | (find_source==2) )[0], :]
mws=mws[np.where( (find_source==1) | (find_source==2) )[0], :]
pres=pres[np.where( (find_source==1) | (find_source==2) )[0], :]
rmw=rmw[np.where( (find_source==1) | (find_source==2) )[0], :]
roci=roci[np.where( (find_source==1) | (find_source==2) )[0], :]
poci=poci[np.where( (find_source==1) | (find_source==2) )[0], :]
eye=eye[np.where( (find_source==1) | (find_source==2) )[0], :]

r34=r34[np.where( (find_source==1) | (find_source==2) )[0], :, :]
r50=r50[np.where( (find_source==1) | (find_source==2) )[0], :, :]
r64=r64[np.where( (find_source==1) | (find_source==2) )[0], :, :]

source=source[np.where( (find_source==1) | (find_source==2) )[0], :]
basin=basin[np.where( (find_source==1) | (find_source==2) )[0], :]
sid=sid[np.where( (find_source==1) | (find_source==2) )[0], :]

del(find_source)



# %% Remove 0 or 1 length data
find_length=np.zeros(len(time[:, 0]))      
for N in range(len(time[:, 0])):
    find_length[N]=len(np.where(np.isnan(time[N, :])==False)[0])
del(N)    

time=time[np.where(find_length>1)[0], :]
dist2land=dist2land[np.where(find_length>1)[0], :]

lon=lon[np.where(find_length>1)[0], :]
lat=lat[np.where(find_length>1)[0], :]
mws=mws[np.where(find_length>1)[0], :]
pres=pres[np.where(find_length>1)[0], :]
rmw=rmw[np.where(find_length>1)[0], :]
roci=roci[np.where(find_length>1)[0], :]
poci=poci[np.where(find_length>1)[0], :]
eye=eye[np.where(find_length>1)[0], :]

r34=r34[np.where(find_length>1)[0], :, :]
r50=r50[np.where(find_length>1)[0], :, :]
r64=r64[np.where(find_length>1)[0], :, :]

source=source[np.where(find_length>1)[0], :]
basin=basin[np.where(find_length>1)[0], :]
sid=sid[np.where(find_length>1)[0], :]

del(find_length)
# %% Julian date to Gregorian date
Year=np.zeros(time.shape, dtype=float)
Month=np.zeros(time.shape, dtype=float)
Day=np.zeros(time.shape, dtype=float)
Hour=np.zeros(time.shape, dtype=float)

for N in range(len(time[:, 0])):
    print('01.Julian date to Gregorian date -', N, '/', len(time[:, 0])-1)
    for L in range(len(time[0, :])):
        if np.isnan(time[N, L])==False:
            Year[N, L]=julian.from_jd(time[N, L], fmt='mjd').year
            Month[N, L]=julian.from_jd(time[N, L], fmt='mjd').month
            Day[N, L]=julian.from_jd(time[N, L], fmt='mjd').day
            Hour[N, L]=julian.from_jd(time[N, L], fmt='mjd').hour
        elif np.isnan(time[N, L])==True:
            Year[N, L]=np.NaN
            Month[N, L]=np.NaN
            Day[N, L]=np.NaN
            Hour[N, L]=np.NaN
    del(L)
del(N)    
        

# %% Make Track file
Variable=29
Track=np.zeros((len(time[:, 0]), 1000, Variable), dtype=np.float64)
Track[np.where(Track==0)]=np.NaN

time_factor=0.04166666667
for N in range(len(time[:, 0])):
    print('02.Make Track file -', N, '/', len(time[:, 0])-1)
    
    ID=int(sid[N, 0]+sid[N, 1]+sid[N, 2]+sid[N, 3]+sid[N, 4]+sid[N, 5]+sid[N, 6]+sid[N, 8]+sid[N, 9]+sid[N, 10]+
    sid[N, 11]+sid[N, 12])
    
    Count=0
    isnan_time=time[N, np.where(np.isnan(time[N, :])==False)]
    for L in np.arange(isnan_time[0, 0], isnan_time[0, len(isnan_time[0, :])-1], time_factor):
        Track[N, Count, 0]=ID                                      # SID    [0]
        Track[N, Count, 1]=L                                       # Julday [1]
        Track[N, Count, 2]=julian.from_jd(L, fmt='mjd').year       # Year   [2]
        Track[N, Count, 3]=julian.from_jd(L, fmt='mjd').month      # Month  [3]
        Track[N, Count, 4]=julian.from_jd(L, fmt='mjd').day        # Day    [4]
        Track[N, Count, 5]=julian.from_jd(L, fmt='mjd').hour       # Hour   [5]

        Count=Count+1
    del(L)        
    del(Count, isnan_time, ID)
del(N)
del(time_factor)

print(int(Track[0, 0, 0]))

# %% Fill in the variables
for N in range(len(time[:, 0])):
    print('03.Fill in the Variables -', N, '/', len(time[:, 0])-1)
    for L in range(len(time[0, :])):
        find_same_time=np.where((Year[N, L]==Track[N, :, 2]) & (Month[N, L]==Track[N, :, 3]) & (Day[N, L]==Track[N, :, 4]) & (Hour[N, L]==Track[N, :, 5]))[0]
        Track[N, find_same_time, 6]=lon[N, L]                      # Lon         [6]
        Track[N, find_same_time, 7]=lat[N, L]                      # Lat         [7]
        Track[N, find_same_time, 8]=mws[N, L]                      # MWS         [8]
        Track[N, find_same_time, 9]=dist2land[N, L]                # DIST2LAND   [9]
        Track[N, find_same_time, 10]=pres[N, L]                    # Pres        [10]
        Track[N, find_same_time, 11]=poci[N, L]                    # POCI        [11]
        Track[N, find_same_time, 12]=eye[N, L]                     # EYE         [12]
        
        Track[N, find_same_time, 13]=rmw[N, L]                     # RMW         [13]
        Track[N, find_same_time, 14]=roci[N, L]                    # ROCI        [14]
        
        Track[N, find_same_time, 15:19]=r34[N, L, :]               # R34      [15-18]
        Track[N, find_same_time, 19:23]=r50[N, L, :]               # R50      [19-22]
        Track[N, find_same_time, 23:27]=r64[N, L, :]               # R64      [23-26]
    del(L)
    del(find_same_time)
del(N)    
        
del(Year, Month, Day, Hour, time)
del(sid, basin, lon, lat, dist2land, mws, pres, eye, rmw, roci, poci, r34, r50, r64)
# %% 1-h interpolation
for N in range(len(Track[:, 0, 1])):
    print('04.1-h interpolation -', N, '/', len(Track[:, 0, 1])-1)
    for L in np.arange(6, Variable, 1):
        
        x=np.where(np.isnan(Track[N, :, L])==False)[0]
        
        if len(x)>=2:
            y=Track[N, x, L]
            f=interpolate.interp1d(x, y, kind='linear')
            
            xnew=np.arange(x[0], x[len(x)-1]+1)
            ynew=f(xnew)
            
            Track[N, x[0]:x[len(x)-1]+1, L]=ynew
            
            del(y, f, xnew, ynew)
        del(x)
    del(L)   
del(N)        

del(Variable)
# %% Heading angle

for N in range(len(Track[:, 0, 0])):
    print('05.Heading angle -', N, '/', len(Track[:, 0, 0])-1)
    for L in np.arange(1, len(Track[0, :, 0])):
                
        dist_x=Track[N, L, 6]-Track[N, L-1, 6]
        dist_y=Track[N, L, 7]-Track[N, L-1, 7]

        if (dist_x > 0) & (dist_y >= 0):             # Q1
            if (dist_y!=0):
                Track[N, L, 27]=90-(math.degrees(math.atan(abs(dist_y)/abs(dist_x))))   # Heading angle [27]
            elif (dist_y==0):   
                Track[N, L, 27]=90
                
        elif (dist_x <= 0) & (dist_y > 0):           # Q2
            if (dist_x!=0):
                Track[N, L, 27]=270+(math.degrees(math.atan(abs(dist_y)/abs(dist_x))))  # Heading angle [27]
            elif (dist_x==0):
                Track[N, L, 27]=0

        elif (dist_x < 0) & (dist_y <= 0):           # Q3
            if (dist_y!=0):
                Track[N, L, 27]=270-(math.degrees(math.atan(abs(dist_y)/abs(dist_x))))  # Heading angle [27]
            elif (dist_y==0):
                Track[N, L, 27]=270
                
        elif (dist_x >= 0) & (dist_y < 0):           # Q4
            if (dist_x!=0):
                Track[N, L, 27]=90+(math.degrees(math.atan(abs(dist_y)/abs(dist_x))))   # Heading angle [27]
            elif (dist_x==0):
                Track[N, L, 27]=180
                
        del(dist_x, dist_y)
    del(L)  
del(N)


# %% Source information

source_info=np.zeros(len(Track[:, 0, 0]))  
      
for N in range(len(Track[:, 0, 0])):
    if len(np.where( (source[N, 0]==b'b') & (source[N, 1]==b'w') & (source[N, 2]==b'p'))[0])>0: # bwp
        source_info[N]=1
    elif len(np.where( (source[N, 0]==b'B') & (source[N, 1]==b'O') & (source[N, 2]==b'T'))[0])>0: # tcvitals
        source_info[N]=2       
del(N)  

for N in range(len(Track[:, 0, 0])):
    Track[N, np.where( np.isnan(Track[N, :, 0])==0 )[0], 28]=source_info[N]             # Source information [28] 1=bwp, 2=tcvitals
del(N)

del(source, source_info)
# %% Save data
from datetime import datetime
now = datetime.now().strftime("%Y%m%d")

np.save(DIR + 'IBTrACS_' + now, Track)

Track_JTWC=Track[ np.where( Track[:, 0, 28]==1 )[0], :, : ]
Track_TCvitals=Track[ np.where( Track[:, 0, 28]==2 )[0], :, : ]

np.save(DIR + 'IBTrACS_JTWC_' + now, Track_JTWC)
np.save(DIR + 'IBTrACS_TCvitals_' + now, Track_TCvitals)
 
del(Track_JTWC, Track_TCvitals)

# %%
select_date=np.zeros((len(Track[:, 0, 0]), 7))
for L in range(len(Track[:, 0, 0])):
    find_data=Track[L, np.where(np.isnan(Track[L, :, 0])==0), 2:5].squeeze()
    select_date[L, :3]=find_data[0, :]
    select_date[L, 3:6]=find_data[len(find_data[:, 0])-1, :]
    select_date[L, 6]=Track[L, 0, 28]
    del(find_data)
del(L)    
    
np.save(DIR + 'IBTrACS_date_' + now, select_date)
    
