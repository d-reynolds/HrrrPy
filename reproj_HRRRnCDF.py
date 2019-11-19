# Script to download archived HRRR data and clip to supplied bounds
#
# RELEASE NOTES
#   Version 1.0 Written by Dylan Reynolds (reyno18@uw.edu), Feb 2019)


import get_archived as ga
import pygrib
from netCDF4 import Dataset
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon
import numpy as np
import gdal
import os
import xarray as xr
import cfgrib
##-------------------------------------------------------------------------------------------------------
##USER DEFINED INPUTS
##-------------------------------------------------------------------------------------------------------
storage_path = '/storage/dylan/HRRR/processed/'
filename = 'Stitched_SM_WY2017.nc'
##-------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------

inputFile = xr.open_dataset(storage_path + filename)
domain_lat = inputFile.latitude.values
domain_lon = inputFile.longitude.values
NUM_HRS = inputFile.time.shape[0]

print('Data clipped, reprojecting to Lambert Conformal coordinates...')
#use geopandas to reporoject lat-lon in local UTM projection

row_to_shapely = lambda row, col : Point(domain_lon[row,col],domain_lat[row,col])
points = [row_to_shapely(i,j) for i in range(0,domain_lon.shape[0]) for j in range(0,np.transpose(domain_lon).shape[0])]
met_df = pd.DataFrame(points,columns=['geometry'])
met_gdf = gpd.GeoDataFrame(met_df,geometry='geometry')

#set latlon_gdf to projection used in grib file so it can be converted
lcc_proj_string = "+proj=lcc +a=6371229 +b=6371229 +lon_0=262.5 +lat_0=38.5 +lat_1=38.5 +lat_2=38.5"
lat_lon_string = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84"
met_gdf.crs = lat_lon_string
proj_gdf = met_gdf.to_crs(lcc_proj_string)#,str(epsg_code))})

print('Creating netCDF file structure...')

#Create netCDF structure and initilaize with dimensions given by coordinate lengths and Time length of data
hrrrNCDF = Dataset(("/storage/dylan/HRRR/processed/reproj_"+filename), "w", format="NETCDF4_CLASSIC")
hrrrNCDF.createDimension('x',np.transpose(domain_lon).shape[0])
hrrrNCDF.createDimension('y',domain_lat.shape[0])
hrrrNCDF.createDimension('height_above_ground1',1)
hrrrNCDF.createDimension('height_above_ground2',1)
hrrrNCDF.createDimension('time',NUM_HRS)

#Fill in netCDF structure with necesarry variables and var-attributes

Temperature_height_above_ground = \
hrrrNCDF.createVariable('Temperature_height_above_ground','f4',("time",'height_above_ground2','y','x',))
Temperature_height_above_ground.units = "K"
Temperature_height_above_ground.missing_value = -9999
Temperature_height_above_ground.long_name = "Temperature @ height_above_ground"
Temperature_height_above_ground.grid_mapping = "LCC_projection"
Temperature_height_above_ground.coordinates = "time_run time height_above_ground2 y x "

time = hrrrNCDF.createVariable('time','f8',("time",))
time.long_name = "Forecast time for ForecastModelRunCollection"
time.standard_name = "time"
time.units = ("hours since " + str(inputFile.time.values[0]))
time._CoordinateAxisType = "Time"

time_run = hrrrNCDF.createVariable('time_run','f8',("time",))
time_run.long_name = "run times for coordinate = time"
time_run.standard_name = "forecast_reference_time"
time_run.units = ("hours since " + str(inputFile.time.values[0]))
time_run._CoordinateAxisType = "RunTime"

height_above_ground2 = \
hrrrNCDF.createVariable('height_above_ground2','f8',('height_above_ground2',))
height_above_ground2.units = "m"
height_above_ground2.long_name = "Specified height level above ground"
height_above_ground2.positive = "up"
height_above_ground2._CoordinateAxisType = "Height"
height_above_ground2._CoordinateZisPositive = "up"

y = hrrrNCDF.createVariable('y','f8',('y',))
y.untis = "m"
y.standard_name = "projection_y_coordinate"
y._CoordianteAxisType = "GeoY"

x = hrrrNCDF.createVariable('x','f8',('x',))
x.units = "m"
x.standard_name = "projection_x_coordinate"
x._CoordianteAxisType = "GeoX"

lat = hrrrNCDF.createVariable('lat','f8',('y','x',))
lat.standard_name = "latitude"

lon = hrrrNCDF.createVariable('lon','f8',('y','x',))
lon.standard_name = "longitude"

elev = hrrrNCDF.createVariable('elev','f8',('y','x',))
elev.units = "m"
elev.long_name = "Height of terrain at model coordinate"

LCC_projection = hrrrNCDF.createVariable('LCC_projection','S1')
LCC_projection.grid_mapping_name = "lambert_conformal_conic"
LCC_projection.longitude_of_central_meridian = 262.5
LCC_projection.latitude_of_projection_origin = 38.5
LCC_projection.standard_parallel = [38.5, 38.5]
LCC_projection.semi_minor_axis = 6371229
LCC_projection.semi_major_axis = 6371229
LCC_projection._CoordinateTransformType = "Projection"
LCC_projection._CoordinateAxisTypes = "GeoX GeoY"

V_component_of_wind_height_above_ground = \
hrrrNCDF.createVariable('V-component_of_wind_height_above_ground','f4',("time",'height_above_ground1','y','x',))
V_component_of_wind_height_above_ground.missing_value = np.nan
V_component_of_wind_height_above_ground.units = "m s-1"
V_component_of_wind_height_above_ground.grid_mapping = "LCC_projection"
V_component_of_wind_height_above_ground.coordinates = "time_run time height_above_ground1 y x "

U_component_of_wind_height_above_ground = \
hrrrNCDF.createVariable('U-component_of_wind_height_above_ground','f4',("time",'height_above_ground1','y','x',))
U_component_of_wind_height_above_ground.missing_value = np.nan
U_component_of_wind_height_above_ground.units = "m s-1"
U_component_of_wind_height_above_ground.grid_mapping = "LCC_projection"
U_component_of_wind_height_above_ground.coordinates = "time_run time height_above_ground1 y x "

height_above_ground1 = \
hrrrNCDF.createVariable('height_above_ground1','f8',('height_above_ground1',))
height_above_ground1.units = "m"
height_above_ground1.long_name = "Specified height level above ground"
height_above_ground1.positive = "up"
height_above_ground1.GRIB_level_type = "103"
height_above_ground1._CoordianteAxisType = "Height"
height_above_ground1._CoordianteZisPositive = "up"


Total_cloud_cover = \
hrrrNCDF.createVariable('Total_cloud_cover','f4',("time",'y','x',))
Total_cloud_cover.long_name = "Total cloud cover @ entire atmosphere"
Total_cloud_cover.missing_value = -9999
Total_cloud_cover.units = "percent"
Total_cloud_cover.grid_mapping = "LCC_projection"
Total_cloud_cover.coordinates = "time_run time y x "

PSFC = \
hrrrNCDF.createVariable('PSFC','f8',("time",'height_above_ground2','y','x',))
PSFC.long_name = "PSFC"
PSFC.units = "Pa"
PSFC.grid_mapping = "LCC_projection"
PSFC.coordinates = "time height_above_ground2 y x "

SW_in= \
hrrrNCDF.createVariable('SW_in','f8',("time",'height_above_ground2','y','x',))
SW_in.long_name = "Incoming shortwave radiation"
SW_in.units = "W/m^2"
SW_in.grid_mapping = "LCC_projection"
SW_in.coordinates = "time height_above_ground2 y x "

LW_in= \
hrrrNCDF.createVariable('LW_in','f8',("time",'height_above_ground2','y','x',))
LW_in.long_name = "Incoming longwave radiation"
LW_in.units = "W/m^2"
LW_in.grid_mapping = "LCC_projection"
LW_in.coordinates = "time height_above_ground2 y x "

Q2= \
hrrrNCDF.createVariable('Q2','f8',("time",'height_above_ground2','y','x',))
Q2.long_name = "Specific Humidity at height above ground"
Q2.units = "kg/kg"
Q2.grid_mapping = "LCC_projection"
Q2.coordinates = "time height_above_ground2 y x "

RAINCV = \
hrrrNCDF.createVariable('RAINCV','f8',("time",'height_above_ground2','y','x',))
RAINCV.long_name = "Hourly precipitation rate"
RAINCV.units = "mm/hr"
RAINCV.grid_mapping = "LCC_projection"
RAINCV.coordinates = "time height_above_ground2 y x "

#Add global attributes

hrrrNCDF.Conventions = "CF-1.6"
hrrrNCDF.Generating_Model = "HRRR"
hrrrNCDF.source = "Archived HRRR"
hrrrNCDF.institution = "University of Washington"
hrrrNCDF.featureType = "GRID"
hrrrNCDF.cdm_data_type = "GRID"

print('Writing data to netCDF file...')
#Add data to variables

if np.any(np.isnan(inputFile.t2m)): print('temp has nans')
if np.any(np.isnan(inputFile.tcc)): print('tcc has nans')
if np.any(np.isnan(inputFile.sp)): print('sp has nans')
if np.any(np.isnan(inputFile.q)): print('q has nans')
if np.any(np.isnan(inputFile.dswrf)): print('sw has nans')
if np.any(np.isnan(inputFile.dlwrf)): print('lw has nans')
if np.any(np.isnan(inputFile.prate)): print('precip has nans')
if np.any(np.isnan(inputFile.u10)): print('u has nans')
if np.any(np.isnan(inputFile.v10)): print('v has nans')

LCC_projection[:] = 'd'
time[:]  = range(0,NUM_HRS)
time_run[:] = 1
height_above_ground2[:]  = 2
height_above_ground1[:]  = 10
#reshape x and y's for easy assignment to variables
tempYs  = np.reshape(proj_gdf['geometry'].y.values,domain_lat.shape)
tempXs  = np.reshape(proj_gdf['geometry'].x.values,domain_lon.shape)
x[:] = tempXs[0,:]
y[:] = tempYs[:,0]
lat[:,:] = domain_lat
lon[:,:] = domain_lon
elev[:,:] = inputFile.orog.isel(time=0)
Temperature_height_above_ground[:,0,:,:] = inputFile.t2m
U_component_of_wind_height_above_ground[:,0,:,:]  = inputFile.u10
V_component_of_wind_height_above_ground[:,0,:,:]  = inputFile.v10
Total_cloud_cover[:,:,:]  = inputFile.tcc
PSFC[:,0,:,:]  = inputFile.sp
Q2[:,0,:,:]  = inputFile.q
SW_in[:,0,:,:]  = inputFile.dswrf
LW_in[:,0,:,:]  = inputFile.dlwrf
RAINCV[:,0,:,:]  = inputFile.prate
print('Done.')


	
