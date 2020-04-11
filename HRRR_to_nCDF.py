# Script to download archived HRRR data and clip to supplied bounds
# 
# Thank you to Brian Blaylock for the use of the "download_HRRR_variable_from_pando" function
# HRRR archive doi: 10.7278/S5JQ0Z5B
#
# RELEASE NOTES
#   Version 1.0 Written by Dylan Reynolds (reyno18@uw.edu), Feb 2019)


import get_archived as ga
import numpy as np
import os
import xarray as xr
import cfgrib
import sys
import tracemalloc
import download_HRRR_variable_from_pando as dHRRR
import glob
from datetime import date

##-------------------------------------------------------------------------------------------------------
##HELPER FUNCTIONS
##-------------------------------------------------------------------------------------------------------

#Used to clear storage directory at end of each loop
def clearDir(inDir):
	files = glob.glob(inDir+'*')
	for f in files:
		os.remove(f)

#Some of the HRRR levels share heightAboveGround and have conflicting values. Remove from record before merging
def preprocessOpen(ds):
	if 'heightAboveGround' in ds.coords:
		ds = ds.drop(['heightAboveGround'])
	if 'valid_time' in ds.coords:
		ds = ds.drop(['valid_time'])
	if 'step' in ds.coords:
		ds = ds.drop(['step'])
	#If we are handling the precip field, shift time forward one hour since it is a forecast
	#This is necesarry to merge all variables along the time dimension
	if 'prate' in ds.data_vars:
		ds.time.values = ds.time.values + np.timedelta64(1,'h')
	return ds

#Function to make sure netCDF file is complete before saving
def checkComplete(ds):
	completeList = ['u10','v10','t2m','tcc','sp','q','prate','dswrf','dlwrf']
	returnList = []
	for var in completeList:
		if not(var in ds.data_vars):
			returnList = np.append(returnList, var)
	return returnList
##-------------------------------------------------------------------------------------------------------
##USER DEFINED INPUTS
##-------------------------------------------------------------------------------------------------------

#TUM cords
LAT_MIN=37.6
LAT_MAX=37.7
LON_MIN=-119.25
LON_MAX=-119.15


#Desired times of HRRR data
start_yr = 2016
end_yr = 2016
start_mon = 11
end_mon = 11
start_day = 12
end_day = 24
start_hr = 12
end_hr = 0



#path to your local storage directory. Script will create a temp workspace and output netCDF workspace
STORAGE_PATH = './'
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
#Precip dates are set 1 hour back since we have to read a forecast file 1 hour in the future
prec_dates = ga.format_dates(start_yr,start_mon,start_day,start_hr,end_yr,end_mon,end_day,end_hr,offset=-1)

NUM_HRS = dates.shape[0]

print('made dates')


#clear raw directory
clearDir(raw_dir)

#GRIB keys in HRRR .grib file for the downloader to read
grib_vars = ['TMP:2 m','UGRD:10 m','VGRD:10 m','SPFH:2 m','PRATE:surface',\
	'DSWRF:surface','DLWRF:surface','HGT:surface','PRES:surface','TCDC:entire']

print('Beginning HRRR downloads...')
for t in range(0,NUM_HRS):
	badRecordFlag=False
	print('Downloading %d of %d...' % ((t+1),NUM_HRS))
	#Download archived HRRR data
	cur_date = date(int(dates.iloc[t,0]),int(dates.iloc[t,1]),int(dates.iloc[t,2]))
	prec_cur_date = date(int(prec_dates.iloc[t,0]),int(prec_dates.iloc[t,1]),int(prec_dates.iloc[t,2]))
	print('starting HRRR downloader')
	try:
		#Download all desired variables from HRRR archive
		for var in grib_vars:
			if ('PRATE' in var):
				dHRRR.download_HRRR_variable_from_pando(prec_cur_date,var,hours=[int(prec_dates.iloc[t,3])],fxx=[1],outdir=raw_dir)
			else:
				dHRRR.download_HRRR_variable_from_pando(cur_date,var,hours=[int(dates.iloc[t,3])],fxx=[0],outdir=raw_dir) 
	except:
		badRecordFlag = True

	print('ending HRRR downloader')
	
	record_vars = xr.Dataset()
	print('starting to open grib files')
	try:
		#Open all downloaded grib files together
		merged_vars = xr.open_mfdataset(raw_dir+'*.grib2',engine='cfgrib',combine='by_coords',preprocess=preprocessOpen)
		print('merged all grib files')	
		#Convert longitude from 0-360 to -180 - 180
		merged_vars.longitude.values = merged_vars.longitude.values-360
		
		#Convert precip from mm/sec to mm/hr
		merged_vars.prate.values = merged_vars.prate.values*3600.0
		
		if t==0:
			#Handle case of a point
			if LAT_MIN == LAT_MAX and LON_MIN == LON_MAX:
				dist = lambda lat1, lon1, lat2, lon2 : np.sqrt(lat1**2-lat2**2+lon1**2-lon2**2)
				closest = 10000
				for lat, lon in zip(merged_vars.latitude.values.flat, merged_vars.longitude.values.flat):
					new = dist(LAT_MAX,LON_MAX,lat,lon)
					if new < closest:
						closest = new
						subset_lat = lat
			else: 
				subset_vars = merged_vars.where((merged_vars.latitude > LAT_MIN)&\
                         		(merged_vars.latitude < LAT_MAX)&\
                         		(merged_vars.longitude > LON_MIN)&\
                         		(merged_vars.longitude < LON_MAX), drop=True)
				subset_lat = subset_vars.latitude
		record_vars = merged_vars.where(merged_vars.latitude.isin(subset_lat), drop=True)
		#Handle case that model grid changes over the course of the archive, giving no points when using initial indices
		if record_vars.x.values.size == 0:
			subset_vars = merged_vars.where((merged_vars.latitude > LAT_MIN)&\
                                        (merged_vars.latitude < LAT_MAX)&\
                                        (merged_vars.longitude > LON_MIN)&\
                                        (merged_vars.longitude < LON_MAX), drop=True)
			subset_lat = subset_vars.latitude
			record_vars = merged_vars.where(merged_vars.latitude.isin(subset_lat), drop=True)
			#Try again -- if obtaining a new point didnt work, it's just a bad record
			if record_vars.x.values.size == 0:
				badRecordFlag = True
	except:
		badRecordFlag = True
	fill_vars = checkComplete(record_vars)
	
	#If there is no cahced HRRR file, or we couldn't read the HRRR file, fill in the record with -9999 values
	if (len(fill_vars) or badRecordFlag):
		print('Error reading grib file, filling record with NaNs')
		with xr.open_dataset(processed_dir+out_filename,engine='netcdf4') as processed_HRRR:
			#Get last time slice from processed_HRRR file
			fill_ds = processed_HRRR.isel(time=[-1])
			for var in fill_vars:
				fill_ds[var].values.fill('-9999.0')
				record_vars[var] = fill_ds[var]
			record_vars['time'] = fill_ds['time']
			record_vars['time'].values = fill_ds['time'].values + np.timedelta64(1, 'h')
	print('Datasubset, saving...')
	print(record_vars)
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
		del combined
	if not(saveHRRR): 
		clearDir(raw_dir)
        
        
print('Done!')
	
