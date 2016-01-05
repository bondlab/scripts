#!/bin/bash

##############################################################################
#
# Author: Jon Badalamenti, Bond lab, University of Minnesota
#   Slightly modified by Kevin Silverstein, MSI, University of Minnesota
#
# Description:
#     A script for linking the name of parent SMRT cell data directories to 
# their corresponding sample based on parsing the metadata.xml file
#
# Usage: get_smrtcell_info.sh
#
# Execute this script within the parent directory containing SMRT cell 
# directories. These directories should have names beginning with A through H,
# followed by two digits, then an underscore, and then one digit. For example:
#  A01_1
#
# This script creates a text file named smrtcell_info and prints the contents 
# of this file to standard out
#
############################################################################## 

# make sure the file smrtcell_info doesn't exist.  If it does remove it
if [ -f smrtcell_info ];
then
        echo removing old smrtcell_info file
        /bin/rm smrtcell_info
fi

# unset and define SMRTCELLS variable as list of contents of current directory
unset SMRTCELLS
SMRTCELLS=$(ls)

# descend into each SMRT cell directory and extract readable information from 
# its metadata.xml file
for s in $SMRTCELLS; do 
	echo parsing metadata for $s;
	cd $s;
	# collect the run date and write to temporary file
	xml_grep WhenStarted --text_only *.metadata.xml | head -n 1 | sed 's/^\(.*\)T.*/\1/' >> parsed;
	# collect the SMRT cell directory name and append temporary file
	xml_grep CollectionPathUri --text_only *.metadata.xml | sed 's/.*\([A-Z][0-9][0-9]_[0-9]\)\/$/\1/' >> parsed; 
	# collect the corresponding sample name and append temporary file
	xml_grep Name --text_only *.metadata.xml | head -n 2 | sed -n '2{p;q}' >>parsed; 
	# collapse parsed data onto one line
	tr '\n' ' ' < parsed > one_line; 
	cd ..; 
done
echo && echo

# concatenate SMRT cell information from all SMRT cell directories into one 
# file containing one line per SMRT cell
for s in $SMRTCELLS; do 
	(cat "$s/one_line"; echo) >> smrtcell_info; 
done

# remove temporary files
for s in $SMRTCELLS; do 
	rm $s/one_line $s/parsed; 
done
	
# print contents of smrtcell_info
INFO=$(<smrtcell_info)
printf "$INFO\n"

exit 0
