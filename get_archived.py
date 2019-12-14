
# Helper functions to download archived HRRR data
#
# RELEASE NOTES
#   Version 1.0 Written by Dylan Reynolds (reyno18@uw.edu), Feb 2019)

import urllib.request
import pandas as pd
import datetime as dt
import numpy as np


def get_archived(storage_path,year, month, day, hour):

	#Ingests the year, month, day, and hour of a desired archived HRRR run.
	#storage_path -- path to the base directory where the HRRR files should be stored
	#year -- year of HRRR run
	#month -- month of desired HRRR run
	#day -- day of desired HRRR run
	#hour -- hour of desired HRRR run

	year = str(year)
	month = str(month)
	day = str(day)
	hour = str(hour)

	if len(year) < 4:
		year = '20'+year
	if len(month) < 2:
		month = '0'+month
	if len(day) < 2:
		day = '0'+day
	if len(hour) < 2:
		hour = '0'+hour

	dwnload_str = 'https://pando-rgw02.chpc.utah.edu/hrrr/sfc/'+year+month+day+'/hrrr.t'+hour+'z.wrfsfcf01.grib2';
	filename = storage_path+'temp_'+year+month+day+hour+'.grib2'
	
	try:
		urllib.request.urlretrieve(dwnload_str,filename)
	except:
		print("No HRRR file archived at: "+dwnload_str)
		return ""
	return filename


def format_dates(st_yr,st_mon,st_day,st_hr,end_yr,end_mon,end_day,end_hr):
	#Helper function to main script. Ingests start and end date, returning a numeric array of dates by
	#year, month, day, hour
	
	pd_range = pd.date_range(start=dt.datetime(st_yr,st_mon,st_day,st_hr),\
				end=dt.datetime(end_yr,end_mon,end_day,end_hr),\
				freq='1H')
	return date_range_2_array(pd_range)

def date_range_2_array(dRange):
        
	#Helper function which takes in a pandas date range and splits it into numeric dates by year, month
	#day, hour

	num_hrs = dRange.shape[0]
	out_dates = pd.DataFrame(np.zeros(shape=(num_hrs,4)),columns=['year','month','day','hour'])
	
	for i in range(0,num_hrs):
                split_date_str = str(dRange[i]).split()
                dates = split_date_str[0].split('-')
                hr = split_date_str[1].split(':')
                out_dates.iloc[i,0] = dates[0]
                out_dates.iloc[i,1] = dates[1]
                out_dates.iloc[i,2] = dates[2]
                out_dates.iloc[i,3] = hr[0]
	return out_dates


