#!/bin/csh

cd /app/work_area/

mkdir temp_stuff

# Random name generation 
set random_name = "`openssl rand -hex 4`"

# To set local environment variables and what not.
setenv PFILES "$HOME/multipipe_$random_name/paramfiles"
setenv GLOG_log_dir "$HOME/multipipe_$random_name/log"
mkdir -p $PFILES $GLOG_log_dir

cp UVIT_DriverModule.par $AS1/uvit/paramfiles/
cp -R $AS1/uvit/paramfiles/* $PFILES

# fake db. 
touch /tmp/USNOA2_VIS_GALEX_NUV_FUV_catalogue.db

# Run the driver module
UVIT_DriverModule n LEVL1AS1UVT*tar_V*

# To remove the temporary log and paramfiles directory.
rm -rf $HOME/multipipe_$random_name

# To remove unwanted files
rm -rf output_RAPC*
rm -rf AUX4*
rm -rf *_level*
rm -f *_EDITED_L2
rm -f *PCONLY*
rm -f Differences.txt
rm -f Driver_reference.txt
rm -f Driver_Total_NUV_SCIENCEDATAFILE.txt
rm -f xytheta.txt
rm -f ZeroCentroid.txt
rm -f core.*
rm -rf app

# To remove the temporary things.
rm -rf /tmp/*
rm -rf temp_stuff


