#!/bin/bash

##############################################################################
#
# Author: Jon Badalamenti, Bond lab, University of Minnesota
#
# Description:
#     Script for moving PacBio secondary analysis runs generated in 
# SMRT Portal to a timestamped directory outside of the ~/smrtanalysis 
# installation tree
#
# Usage: pbsave.sh
#    
# Execute this script BEFORE uninstalling SMRT Analysis from your 
# home directory; otherwise all secondary analysis data in 
# ~/smrtanalysis/userdata/jobs will be deleted!
#
############################################################################## 

# unset and set timestamped foldername variable
# this destination directory for saved PacBio data has the format 
# YYYY-MM-DD_HHMM in 24-hour time

unset foldername
foldername=$(date "+%Y-%m-%d_%H%M")

# create timestamped directory in ~/pb_data_dump if it does not already exist 
# and move contents of ~/smrtanalysis/userdata/jobs/ to this timestamped 
# directory
if [ -d ~/pb_data_dump ];
then
	mkdir ~/pb_data_dump/"$foldername"
	mv ~/smrtanalysis/userdata/jobs/* ~/pb_data_dump/"$foldername"
else
	mkdir -p  ~/pb_data_dump/"$foldername"
	mv ~/smrtanalysis/userdata/jobs/* ~/pb_data_dump/"$foldername"  
fi

# clear foldername variable
unset foldername

exit 0
