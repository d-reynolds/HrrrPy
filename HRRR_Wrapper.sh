#!/bin/bash

set -e

if [[ $OSTYPE == *'linux'* ]]
then
	isLinux=true
else
	isLinux=false
fi

echo
echo
echo Welcome to the HRRR-Downloader!
echo This tool was developed by
echo Dylan Reynolds at the University of Washington
echo \(Contact: reyno18@uw.edu\)
echo 
echo This tool is made possible by the
echo University of Utah\'s HRRR archive
echo managed by Brian Blaylock
echo \(http://hrrr.chpc.utah.edu/\)
echo
echo ----------------------------------------------------
echo Answer these questions about your region of interest
echo ----------------------------------------------------

read -p "MAX Latitude:" max_lat

if (( $(bc <<< "$max_lat > 47.84") || $(bc <<< "$max_lat < 21.14") ))
then
	echo Latitudes must be between 47.84 and 21.14 degrees
	exit 1
fi

read -p "MIN Latitude:" min_lat

if (( $(bc <<< "$min_lat > 47.84") || $(bc <<< "$min_lat < 21.14") ))
then 
        echo Latitudes must be between 47.84 and 21.14 degrees
	exit 1
fi

if (( $(bc <<< "$max_lat < $min_lat") ))
then 
	echo The maximum latitude must be greater than the minimum latitude
        exit 1
fi

read -p "MAX Longitude:" max_lon

if (( $(bc <<< "$max_lon > -60.9") || $(bc <<< "$max_lon < -122.7") ))
then
        echo Longitudes must be between -60.9 and -122.7 degrees
        exit 1
fi

read -p "MIN Longitude:" min_lon


if (( $(bc <<< "$min_lon > -60.9") || $(bc <<< "$min_lon < -122.7") ))
then
        echo Longitudes must be between -60.9 and -122.7 degrees
        exit 1
fi

if (( $(bc <<< "$max_lon < $min_lon") ))
then
        echo The maximum longitude must be greater than the minimum longitude
        exit 1
fi

echo
echo --------------------------------------------
echo Answer these questions about your date range
echo --------------------------------------------

currentDate=$(date +%s)


read -p "START year(YYYY):" st_yr

read -p "START month(mm):" st_mon

read -p "START day(dd):" st_day

read -p "START hour(HH):" st_hr

if $isLinux
then
	startDate=$(date -d"${st_yr}/${st_mon}/${st_day} ${st_hr}:00:00" +%s)
	archiveDate=$(date -d"2016/07/15 00:00:00" +%s)
else
	startDate=$(date -j -f '%Y%m%d%H%M%S' ${st_yr}${st_mon}${st_day}${st_hr}0000 +%s)
	archiveDate=$(date -j -f '%Y%m%d%H%M' 201607150000 +%s)
fi


if [ $startDate -lt $archiveDate ]
then
	echo HRRR archive starts on July 15th, 2016
	exit 1
fi

read -p "END year(YYYY):" ed_yr

read -p "END month(mm):" ed_mon

read -p "END day(dd):" ed_day

read -p "END hour(HH):" ed_hr

if $isLinux
then
	endDate=$(date -d"${ed_yr}/${ed_mon}/${ed_day} ${ed_hr}:00:00" +%s)
else
	endDate=$(date -j -f '%Y%m%d%H%M%S' ${ed_yr}${ed_mon}${ed_day}${ed_hr}0000 +%s)
fi


if [ $endDate -gt $currentDate ]
then
        echo End date cannot be greater than current date
        exit 1
fi

if [ $endDate -lt $startDate ]
then
	echo End date must be greater than start date
	exit 1
fi

echo ---------------------------------
echo
echo Couple more questions... 
echo
echo What directory should we save to? \(include trailing \"/\"\)
read save_dir

mkdir -p $save_dir

echo
echo What file name should we use? \(w/o extension\)
read fname

echo
echo Should we run in the background \([y]/[n]\)?
read backgroundQ

echo -----------------------------------
echo
echo Press "Enter" to start
read dummy



if [[ $backgroundQ == *'y'* ]]
then
	nohup python HRRR_to_nCDF.py $max_lat $min_lat $max_lon $min_lon $st_yr $st_mon $st_day $st_hr $ed_yr $ed_mon $ed_day $ed_hr $save_dir $fname > ${save_dir}+${fname}+'.out' &

elif [[ $backgroundQ == *'n'* ]]
then
	python HRRR_to_nCDF.py $max_lat $min_lat $max_lon $min_lon $st_yr $st_mon $st_day $st_hr $ed_yr $ed_mon $ed_day $ed_hr $save_dir $fname
fi
