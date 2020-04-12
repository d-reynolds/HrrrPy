# HrrrPy
Python tools for downloading archived HRRR data from the University of Utah HRRR archive (http://hrrr.chpc.utah.edu/)

The main python script currently pulls: Temperature, Specific Humidity, LW, SW, 10-m winds, Precipitation, Surface Pressure, and Total Cloud Cover. If you plan to download large spatial/temporal datasets, allocate a few days for the script to run.

Commands to setup the project for use:
```bash
git clone https://github.com/d-reynolds/HrrrPy.git
cd HrrrPy
conda env create -f HrrrPy.yml
conda activate HrrrPyEnv
chmod +x HRRR_Wrapper.sh
```
 
To dowload archived HRRR data:
* Run "./HRRR_Wrapper.sh" from the directory
* Follow the prompts 

To reproject HRRR data for input to WindNinja/MicroMet:
* Edit the "User Defined Variables" Section of reproj_HRRRnCDF.py to reflect the name of the HRRR nCDF file
* Run "python reproj_HRRRnCDF.py"

To process reprojected HRRR data into SnowModel .dat file:
* Edit the "User Defined Variables" Section of nCDF_to_snowmodel.py 
* Run "python nCDF_to_snowmodel.py"
