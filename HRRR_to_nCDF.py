# Script to download archived HRRR data and clip to supplied bounds
#
# RELEASE NOTES
#   Version 1.0 Written by Dylan Reynolds (reyno18@uw.edu), Feb 2019)


import get_archived as ga
import numpy as np
import os
import xarray as xr
import cfgrib
import sys

##-------------------------------------------------------------------------------------------------------
##USER DEFINED INPUTS
##-------------------------------------------------------------------------------------------------------

#TUM cords
LAT_MIN=37.7
LAT_MAX=38.3
LON_MIN=-119.85
LON_MAX=-119.15


#Desired times of HRRR data
start_yr = 2018
end_yr = 2018
start_mon = 6
end_mon = 10
start_day = 1
end_day = 24
start_hr = 0
end_hr = 0



#path to your local storage directory. Script will create a temp workspace and output netCDF workspace
STORAGE_PATH = '/storage/dylan'
out_filename = 'Killarney_data.nc'

#Flag that determines if raw HRRR grib file should be kept after processing
saveHRRR = False

#Flag to set if we are restarting a failed HRRR download. If so, the user should changing the starting
#date and set restart to True. This will begin processing at the new start date and append to the old record
restart = False
##-------------------------------------------------------------------------------------------------------
##End User Defined Inputs
##-------------------------------------------------------------------------------------------------------

##-------------------------------------------------------------------------------------------------------
##Command-line arguments section
##-------------------------------------------------------------------------------------------------------
#This script can be run with arguments from the command-line. Check to see if we have any, and process if so
if (len(sys.argv) > 1):
	#There are user-supplied arguments -- let's get to it!
	if not(len(sys.argv) == 15):
		print('Not enough input variables have been supplied to the script.')
		print('Continuing with user-defined inputs hard-coded in script.')

	else:
		LAT_MAX = float(sys.argv[1])
		LAT_MIN = float(sys.argv[2])
		LON_MAX = float(sys.argv[3])
		LON_MIN = float(sys.argv[4])

		start_yr = int(sys.argv[5])
		start_mon = int(sys.argv[6])
		start_day = int(sys.argv[7])
		start_hr = int(sys.argv[8])

		end_yr = int(sys.argv[9])
		end_mon = int(sys.argv[10])
		end_day = int(sys.argv[11])
		end_hr = int(sys.argv[12])
		
		STORAGE_PATH = str(sys.argv[13])
		out_filename = str(sys.argv[14])+".nc"

##-------------------------------------------------------------------------------------------------------
##End command-line arguments section
##-------------------------------------------------------------------------------------------------------


#Setup file directories for script
raw_dir = (STORAGE_PATH+"HRRR/raw/")
processed_dir = (STORAGE_PATH+"HRRR/processed/")
if not os.path.isdir(raw_dir):
	os.makedirs(raw_dir)
if not os.path.isdir(processed_dir):
	os.makedirs(processed_dir)

print('created dirs')

#Get list of dates in desired interval for download
dates = ga.format_dates(start_yr,start_mon,start_day,start_hr,end_yr,end_mon,end_day,end_hr)
NUM_HRS = dates.shape[0]
print('made dates')

print('Beginning HRRR downloads...')
for t in range(0,NUM_HRS):
	badRecordFlag = False
	print('Downloading %d of %d...' % ((t+1),NUM_HRS))
	#Download archived HRRR data
	filename = ga.get_archived(raw_dir,dates.iloc[t,0],dates.iloc[t,1],dates.iloc[t,2],dates.iloc[t,3])
	if (filename == ""):
		badRecordFlag = True	

	print('starting to open grib files')
	try:
		HRRR_surf = xr.open_dataset(filename, engine='cfgrib',\
		backend_kwargs=(dict(filter_by_keys={'typeOfLevel': 'surface', 'stepType': 'instant'},indexpath='')))

		HRRR_2m = xr.open_dataset(filename, engine='cfgrib',\
		backend_kwargs=(dict(filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 2},indexpath='')))
	
		HRRR_10m = xr.open_dataset(filename, engine='cfgrib',\
		backend_kwargs=(dict(filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 10},indexpath='')))
	
		HRRR_atm = xr.open_dataset(filename, engine='cfgrib',\
		backend_kwargs=(dict(filter_by_keys={'typeOfLevel': 'atmosphere', 'level': 0},indexpath='')))

		HRRR_10m = HRRR_10m.drop(['heightAboveGround'])
		HRRR_2m = HRRR_2m.drop(['heightAboveGround'])
		#      Convert precip from mm/s to mm/hr
		HRRR_surf.prate.values = HRRR_surf.prate.values*3600.0
		
		merged_vars = xr.merge([HRRR_10m.u10, \
			HRRR_10m.v10,\
		 	HRRR_2m.t2m,\
		 	HRRR_atm.tcc,\
		 	HRRR_surf.sp,\
		 	HRRR_2m.q,\
			HRRR_surf.prate,\
		 	HRRR_surf.dswrf,\
		 	HRRR_surf.dlwrf,\
		 	HRRR_surf.orog])
		merged_vars = merged_vars.drop(['atmosphere','surface','step','valid_time'])

		#Longitude comes out of HRRR in 0-360, convert to -180 - 180
		merged_vars.longitude.values = merged_vars.longitude.values-360
		if t==0:
			subset_vars = merged_vars.where((merged_vars.latitude > LAT_MIN)&\
			 (merged_vars.latitude < LAT_MAX)&\
			 (merged_vars.longitude > LON_MIN)&\
			 (merged_vars.longitude < LON_MAX), drop=True)

			subset_lat = subset_vars.latitude
			subset_lon = subset_vars.longitude

		record_vars = merged_vars.where(merged_vars.latitude.isin(subset_lat), drop=True)
	except:
                badRecordFlag = True


	#If there is no cahced HRRR file, or we couldn't read the HRRR file, fill in the record with -9999 values
	if badRecordFlag:
		print('Error reading grib file, filling record with NaNs')
		with xr.open_dataset(processed_dir+out_filename,engine='netcdf4') as processed_HRRR:
			#Get last time slice from processed_HRRR file
			record_vars = processed_HRRR.isel(time=[-1])
			record_vars['u10'].values.fill('-9999.0')
			record_vars['v10'].values.fill('-9999.0')
			record_vars['t2m'].values.fill('-9999.0')
			record_vars['tcc'].values.fill('-9999.0')
			record_vars['sp'].values.fill('-9999.0')
			record_vars['q'].values.fill('-9999.0')
			record_vars['prate'].values.fill('-9999.0')
			record_vars['dswrf'].values.fill('-9999.0')
			record_vars['dlwrf'].values.fill('-9999.0')
			record_vars['time'].values = record_vars['time'].values + np.timedelta64(1, 'h') 
	print('Datasubset, saving...')

	if t==0 and restart != True:
		#Set up encoding parameters to use compression when writing netcdf file
		comp = dict(zlib=True, complevel=5, _FillValue=-9999.0)
		encoding = {var: comp for var in record_vars.data_vars}
		
		record_vars.time.encoding = {'zlib':True,'complevel':5,'units':('hours since '+str(record_vars.time.values))}
		record_vars.to_netcdf(processed_dir+out_filename,encoding=encoding)
	else:
		with xr.open_dataset(processed_dir+out_filename,engine='netcdf4') as processed_HRRR:
			combined = xr.concat([processed_HRRR, record_vars], dim='time')
		combined.to_netcdf(processed_dir+out_filename,'w')

	if not(saveHRRR) and filename: 
		os.remove(filename)

print('Done!')
	
