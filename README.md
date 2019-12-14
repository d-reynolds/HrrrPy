# HrrrPy
Python tools for downloading archived HRRR data from the University of Utah HRRR archive (http://hrrr.chpc.utah.edu/)

Opening and subsetting HRRR data to desired variables is currently slow due to the CFgrib Xarray driver's current options for opening grib files with heterogeneous messages. This driver still allows for a more robust indexing of variables by name instead of message number, and informative errors. If you plan to download large spatial/temporal datasets, allocate a few days for the script to run.

To dowload archived HRRR data:
* Download this repository 
* Run "HRRR_Wrapper.sh" from the directory
* Follow the prompts 

To reproject HRRR data for input to WindNinja/MicroMet:
* Edit the "User Defined Variables" Section of reproj_HRRRnCDF.py to reflect the name of the HRRR nCDF file
* Run "python reproj_HRRRnCDF.py"

To process reprojected HRRR data into SnowModel .dat file:
* Edit the "User Defined Variables" Section of nCDF_to_snowmodel.py 
* Run "python nCDF_to_snowmodel.py"
