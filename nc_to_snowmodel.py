from netCDF4 import Dataset
import get_archived as ga
import numpy as np
import pandas as pd
import datetime as dt
import geopandas as gpd
from shapely.geometry import Point

##-----------------------------------------------
#User Inputs
inputFilename = '/storage/dylan/HRRR/processed/TUM_HRRR_March_2017.nc'
outputFilename = '/storage/dylan/HRRR/processed/HRRR_March_2017.dat'
projection_Lon = -97.5
projection_Lat = 38.5
##-----------------------------------------------


## Function to take netCDF file and write data to met forcing file for SnowModel


hrrrNCDF = Dataset(inputFilename, 'r', format='NETCDF4_CLASSIC')

x = hrrrNCDF.variables['x'][:]
y = hrrrNCDF.variables['y'][:]
elev = hrrrNCDF.variables['elev'][:,:]

time_units = hrrrNCDF.variables['time'].units
time_str = str(time_units).split()[2]
[dates, hrs] = time_str.split('T')
[yr, mon, day] = dates.split('-')
hr = hrs.split(':')[0]

num_hrs = len(hrrrNCDF.dimensions['time'])


lat = hrrrNCDF.variables['lat'][:,:]
lon = hrrrNCDF.variables['lon'][:,:]

row_to_shapely = lambda row, col : Point(lon[row,col],lat[row,col])
points = [row_to_shapely(i,j) for i in range(0,lat.shape[0]) for j in range(0,lon.shape[1])]
met_df = pd.DataFrame(points,columns=['geometry'])
met_gdf = gpd.GeoDataFrame(met_df,geometry='geometry')

#set latlon_gdf to projection used in grib file so it can be converted
lat_lon_string = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84"
met_gdf.crs = lat_lon_string
proj_gdf = met_gdf.to_crs({'init':'epsg:26911'})#,str(epsg_code))})

Y_UTM  = np.reshape(proj_gdf['geometry'].y.values,(y.size,x.size))
X_UTM  = np.reshape(proj_gdf['geometry'].x.values,(y.size,x.size))

Tair = np.zeros((num_hrs,lat.shape[0],lat.shape[1]))

Tair[:,:,:] = hrrrNCDF.variables['Temperature_height_above_ground'][:,0,:,:]
#Tair = Tair-273.15

Press = hrrrNCDF.variables['PSFC'][:,0,:,:]
Precip = hrrrNCDF.variables['RAINCV'][:,0,:,:]


u10 = hrrrNCDF.variables['U-component_of_wind_height_above_ground'][:,0,:,:]
v10 = hrrrNCDF.variables['V-component_of_wind_height_above_ground'][:,0,:,:]

WindSpd = np.sqrt((u10**2)+(v10**2))
WindDir = np.zeros(u10.shape)
#Conversion of wind direction from Projected coordinate system in HRRR to true-north
#Following fortran code on HRRR website
for i in range(0,x.size):
	for j in range(0,y.size):
		angle2 = np.sin(np.deg2rad(projection_Lat))*(lon[j,i]-projection_Lon)*0.017453
		sinx2 = np.sin(angle2)
		cosx2 = np.cos(angle2)
		u_corr = (cosx2*u10[:,j,i])+(sinx2*v10[:,j,i])
		v_corr = (sinx2*u10[:,j,i])+(cosx2*v10[:,j,i])
		dir_corr = np.rad2deg(np.arctan2(v_corr,u_corr))
		WindDir[:,j,i] = np.where(dir_corr > 0,dir_corr,dir_corr+360)[:]

SpecHum = hrrrNCDF.variables['Q2'][:,0,:,:]
#RH from q -- taken from Justin Pflug's Fortran script. Temp must be in C for this
    
A = 611.21
B = 17.502
C = 240.97
    
es = A * np.exp((B * (Tair))/(C + (Tair)))
e = SpecHum * Press / (0.622 + 0.378 * SpecHum)
    
RH = 100*(e/es)
RH = np.where(RH <= 100,RH,100)
RH = np.where(RH > 0,RH,0) 

x = hrrrNCDF.variables['x'][:]
y = hrrrNCDF.variables['y'][:]
elev = hrrrNCDF.variables['elev'][:,:]

time_units = hrrrNCDF.variables['time'].units
time_str = str(time_units).split()[2]
[dates, hrs] = time_str.split('T')
[yr, mon, day] = dates.split('-')
hr = hrs.split(':')[0]

num_hrs = len(hrrrNCDF.dimensions['time'])

time_series = pd.date_range(start=dt.datetime(int(yr),int(mon),int(day),int(hr)),periods=num_hrs,freq='1H')
date_array = ga.date_range_2_array(time_series)

num_stations = (x.size*y.size)

stationIDs = np.arange(1,num_stations+1)
stationIDs = np.reshape(stationIDs,(x.size,y.size))

f = open(outputFilename,'w+')
for t in range (0,num_hrs):
	f.write('\t%d\n'%num_stations)
	for j in range(0,x.size):
		for i in range(0,y.size):
			dataLine = ('%4d %2d %2d %.2f %d %.1f %.1f %.1f %.2f %.2f %.2f %.2f %.2f\n'%\
				(int(date_array.iloc[t,0]),int(date_array.iloc[t,1]),int(date_array.iloc[t,2]),\
				int(date_array.iloc[t,3]),stationIDs[j,i],X_UTM[i,j],Y_UTM[i,j],elev[i,j],\
				Tair[t,i,j],\
				RH[t,i,j],\
				WindSpd[t,i,j],\
				WindDir[t,i,j],\
				Precip[t,i,j]))
			f.write(dataLine)		

f.close()
