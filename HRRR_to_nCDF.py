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

##-------------------------------------------------------------------------------------------------------
##USER DEFINED INPUTS
##-------------------------------------------------------------------------------------------------------

LAT_MIN=37.7
LAT_MAX=38.3
LON_MIN=-119.85
LON_MAX=-119.15

#Desired times of HRRR data
start_yr = 2017
end_yr = 2017
start_mon = 3
end_mon = 3
start_day = 19
end_day = 19
start_hr = 13
end_hr = 13

#path to your local storage directory. Script will create a temp workspace and output netCDF workspace
STORAGE_PATH = '/storage/dylan'
out_filename = 'TUM_HRRR_QuickTest.nc'
saveHRRR = False
##-------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------


#Indexes of grib messages in HRRR .grib2 output file
U_10_INDEX = 59
V_10_INDEX = 60
T_2_INDEX = 54
CLOUDCOVER_INDEX = 86
PRESSSURF_INDEX = 45
SPECHUM_INDEX = 56
PRECIP_INDEX = 63
OROGRAPHY_INDEX = 46

#Setup file directories for script
raw_dir = (STORAGE_PATH+"/HRRR/raw")
processed_dir = (STORAGE_PATH+"/HRRR/processed")
if not os.path.isdir(raw_dir):
	os.mkdir(raw_dir)
if not os.path.isdir(processed_dir):
	os.mkdir(processed_dir)

#Get list of dates in desired interval for download
dates = ga.format_dates(start_yr,start_mon,start_day,start_hr,end_yr,end_mon,end_day,end_hr)
NUM_HRS = dates.shape[0]


print('Beginning HRRR downloads...')
for t in range(0,NUM_HRS):
	print('Downloading %d of %d...' % ((t+1),NUM_HRS))
	#Download archived HRRR data
	filename = ga.get_archived(STORAGE_PATH,dates.iloc[t,0],dates.iloc[t,1],dates.iloc[t,2],dates.iloc[t,3])

	grbs = pygrib.open(filename)

	U_10m_temp, lats, lons = grbs.message(U_10_INDEX).data()
	V_10m_temp, lats, lons = grbs.message(V_10_INDEX).data()
	T_2m_temp, lats, lons = grbs.message(T_2_INDEX).data()
	CloudCover_temp, lats, lons = grbs.message(CLOUDCOVER_INDEX).data()
	Press_Surf_temp, lats, lons = grbs.message(PRESSSURF_INDEX).data()
	SpecHum_temp, lats, lons = grbs.message(SPECHUM_INDEX).data()
	PrecipRate_temp, lats, lons = grbs.message(PRECIP_INDEX).data()
	
	lat_min_cord = 0
	lon_min_cord = 0
	lat_max_cord = 0
	lon_max_cord = 0
	
	#If it is the first download, create variable structures and find bounds
	if t == 0:
		print('Projection parameters are: ',grbs.message(63).projparams)
		print('Finding domain bounds...')
	
		LLC = [0,0]
		URC = [0,0]
		LRC = [0,0]
		ULC = [0,0]
		
		#functions to find distance of given coordinate to lower-left-corner & upper-right-corner of domain
		LLC_dist = lambda lat, lon : np.sqrt((abs(LAT_MIN-lat)**2)+(abs(LON_MIN-lon)**2))
		URC_dist = lambda lat, lon : np.sqrt((abs(LAT_MAX-lat)**2)+(abs(LON_MAX-lon)**2))		
		LRC_dist = lambda lat, lon : np.sqrt((abs(LAT_MIN-lat)**2)+(abs(LON_MAX-lon)**2))
		ULC_dist = lambda lat, lon : np.sqrt((abs(LAT_MAX-lat)**2)+(abs(LON_MIN-lon)**2))

		LLclosest = 10000
		URclosest = 10000
		LRclosest = 10000
		ULclosest = 10000

		for j in range(0,lats.shape[0]):
			for i in range(0,lons.shape[1]):
				if (LLC_dist(lats[j,i],lons[j,i]) < LLclosest):
					LLclosest = LLC_dist(lats[j,i],lons[j,i])
					LLC = [j,i]
				if (URC_dist(lats[j,i],lons[j,i]) < URclosest):
					URclosest = URC_dist(lats[j,i],lons[j,i])
					URC = [j,i]
				if (LRC_dist(lats[j,i],lons[j,i]) < LRclosest):
					LRclosest = LRC_dist(lats[j,i],lons[j,i])
					LRC = [j,i]
				if (ULC_dist(lats[j,i],lons[j,i]) < ULclosest):
					ULclosest = ULC_dist(lats[j,i],lons[j,i])
					ULC = [j,i]

		lat_max_cord = max(URC[0],ULC[0])
		lat_min_cord = min(LRC[0],LLC[0])
		lon_max_cord = max(LRC[1],URC[1])
		lon_min_cord = min(LLC[1],ULC[1])

		num_rows = lat_max_cord-lat_min_cord
		num_cols = lon_max_cord-lon_min_cord

		U_10m = np.zeros(shape=(NUM_HRS,num_rows,num_cols))
		V_10m = np.zeros(shape=(NUM_HRS,num_rows,num_cols))
		T_2m = np.zeros(shape=(NUM_HRS,num_rows,num_cols))
		CloudCover = np.zeros(shape=(NUM_HRS,num_rows,num_cols))
		Press_Surf = np.zeros(shape=(NUM_HRS,num_rows,num_cols))
		SpecHum = np.zeros(shape=(NUM_HRS,num_rows,num_cols))
		PrecipRate = np.zeros(shape=(NUM_HRS,num_rows,num_cols))

		orography, lats, lons = grbs.message(OROGRAPHY_INDEX).data()

		domain_lat = lats[lat_min_cord:lat_max_cord,lon_min_cord:lon_max_cord]
		domain_lon = lons[lat_min_cord:lat_max_cord,lon_min_cord:lon_max_cord]
		domain_elev = orography[lat_min_cord:lat_max_cord,lon_min_cord:lon_max_cord]

	for i in range(lat_min_cord,lat_max_cord):
		for j  in range(lon_min_cord,lon_max_cord):
			U_10m[t,(i-lat_min_cord),(j-lon_min_cord)]=U_10m_temp[i,j]
			V_10m[t,(i-lat_min_cord),(j-lon_min_cord)]=V_10m_temp[i,j]
			T_2m[t,(i-lat_min_cord),(j-lon_min_cord)]=T_2m_temp[i,j]
			CloudCover[t,(i-lat_min_cord),(j-lon_min_cord)]=CloudCover_temp[i,j]
			Press_Surf[t,(i-lat_min_cord),(j-lon_min_cord)]=Press_Surf_temp[i,j]
			SpecHum[t,(i-lat_min_cord),(j-lon_min_cord)]=SpecHum_temp[i,j]
			PrecipRate[t,(i-lat_min_cord),(j-lon_min_cord)]=PrecipRate_temp[i,j]
	grbs.close()
	if not saveHRRR : os.remove(filename)

print('Data clipped, reprojecting to Lambert Conformal coordinates...')
#use geopandas to reporoject lat-lon in local UTM projection
row_to_shapely = lambda row, col : Point(domain_lon[row,col],domain_lat[row,col])
points = [row_to_shapely(i,j) for i in range(0,domain_lon.shape[0]) for j in range(0,domain_lon.shape[1])]
met_df = pd.DataFrame(points,columns=['geometry'])
met_gdf = gpd.GeoDataFrame(met_df,geometry='geometry')

#set latlon_gdf to projection used in grib file so it can be converted
lcc_proj_string = "+proj=lcc +a=6371229 +b=6371229 +lon_0=262.5 +lat_0=38.5 +lat_1=38.5 +lat_2=38.5"
lat_lon_string = "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84"
met_gdf.crs = lat_lon_string
proj_gdf = met_gdf.to_crs(lcc_proj_string)#,str(epsg_code))})

print('Creating netCDF file structure...')

#Create netCDF structure and initilaize with dimensions given by coordinate lengths and Time length of data
hrrrNCDF = Dataset(("/storage/dylan/HRRR/processed/"+out_filename), "w", format="NETCDF4_CLASSIC")

hrrrNCDF.createDimension('x',domain_lon.shape[1])
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
time.units = ("hours since %d-%02d-%02dT%02d:00:00Z"%\
        (int(dates.iloc[0,0]),int(dates.iloc[0,1]),int(dates.iloc[0,2]),int(dates.iloc[0,3])))
time._CoordinateAxisType = "Time"

time_run = hrrrNCDF.createVariable('time_run','f8',("time",))
time_run.long_name = "run times for coordinate = time"
time_run.standard_name = "forecast_reference_time"
time_run.units = ("hours since %d-%02d-%02dT%02d:00:00Z"%\
	(int(dates.iloc[0,0]),int(dates.iloc[0,1]),int(dates.iloc[0,2]),int(dates.iloc[0,3])))
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
print(np.min(T_2m))
print(np.max(T_2m))
#Add data to variables
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
elev[:,:] = domain_elev
Temperature_height_above_ground[:,0,:,:] = T_2m
U_component_of_wind_height_above_ground[:,0,:,:]  = U_10m
V_component_of_wind_height_above_ground[:,0,:,:]  = V_10m
Total_cloud_cover[:,:,:]  = CloudCover
PSFC[:,0,:,:]  = Press_Surf
Q2[:,0,:,:]  = SpecHum
#Convert from HRRR units of kg/m^2/s to mm/hr
RAINCV[:,0,:,:]  = (PrecipRate*3600)

print('Done.')


	
