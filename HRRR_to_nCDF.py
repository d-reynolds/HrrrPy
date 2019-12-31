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
		return ds.drop(['heightAboveGround'])
	else:
		return ds

##-------------------------------------------------------------------------------------------------------
##USER DEFINED INPUTS
##-------------------------------------------------------------------------------------------------------

#TUM cords
LAT_MIN=37.7
LAT_MAX=37.7
LON_MIN=-119.15
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
NUM_HRS = dates.shape[0]

print('made dates')


#clear raw directory
clearDir(raw_dir)

#GRIB keys in HRRR .grib file for the downloader to read
grib_vars = ['TMP:2 m','UGRD:10 m','VGRD:10 m','SPFH:2 m','PRATE:surface',\
	'DSWRF:surface','DLWRF:surface','HGT:surface','PRES:surface','TCDC:entire']

print('Beginning HRRR downloads...')
for t in range(0,NUM_HRS):
	badRecordFlag = False
	print('Downloading %d of %d...' % ((t+1),NUM_HRS))
	#Download archived HRRR data
	cur_date = date(int(dates.iloc[t,0]),int(dates.iloc[t,1]),int(dates.iloc[t,2]))
	print('starting HRRR downloader')
	try:
		#Download all desired variables from HRRR archive
		for var in grib_vars:
			dHRRR.download_HRRR_variable_from_pando(cur_date,var,hours=[int(dates.iloc[t,3])],outdir=raw_dir) 
	except:
		badRecordFlag = True

	print('ending HRRR downloader')
	

	print('starting to open grib files')
	try:
		#Open all downloaded grib files together
		merged_vars = xr.open_mfdataset(raw_dir+'*.grib2',engine='cfgrib',combine='by_coords',preprocess=preprocessOpen)
		
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
		del combined
	if not(saveHRRR): 
		clearDir(raw_dir)
        
        
print('Done!')
	
