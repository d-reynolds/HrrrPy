# HrrrPy
Python tools for downloading archived HRRR data from the University of Utah HRRR archive (http://hrrr.chpc.utah.edu/)

To dowload archived HRRR data:
*Edit the "User Defined Variables" Section of HRRR_to_nCDF.py to reflect the desired lat-lon and date range
*Run "python HRRR_to_nCDF.py"

To reproject HRRR data for input to WindNinja/MicroMet:
*Edit the "User Defined Variables" Section of reproj_HRRRnCDF.py to reflect the name of the HRRR nCDF file
*Run "python reproj_HRRRnCDF.py"

To process reprojected HRRR data into SnowModel .dat file:
*Edit the "User Defined Variables" Section of nCDF_to_snowmodel.py 
*Run "python nCDF_to_snowmodel.py"
