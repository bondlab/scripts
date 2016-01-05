#!/bin/bash

##############################################################################
# Author: Jon Badalamenti, Bond Lab, University of Minnesota
#
# Description: 
#     Script for extracting unique read IDs from an unmappedSubreads.fasta
#	  file, so as to generate a whitelist for re-running analyses on a subset of reads
#
# See https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/ \
#     HGAP-Whitelisting-Tutorial
#
# Usage: make_whitelist.sh
#     
# must be run within the /data directory corresponding to any secondary
# analysis job that creates an unmappedSubreads.fasta file (e.g. HGAP,
# resequencing, base mod detection). Example: 
# 
# ~/smrtanalysis/userdata/jobs/01N/01NXXX/data/unmappedSubreads.fasta
#
# this directory must contain an unmappedSubreads.fasta file
# requires an active MSI isub session: https://www.msi.umn.edu/isub
#
# Especially useful for re-running HGAP on a subset of reads 
#
############################################################################## 
set -e

# evaluate if unmappedSubreads.fasta file exists
if [ -f ./unmappedSubreads.fasta ];
then
	# extract headers from unmappedSubreads.fasta
	# remove leading '>'
	# remove trailing coordinates
	# sort and identify unique read IDs
	# write output to whitelist.list
	echo extracting read IDs from unmapped subreads...
	grep '>' unmappedSubreads.fasta | sed 's/>//' | sed 's/\//\t/2' | cut -f1 | sort | uniq > whitelist.list
	# collect total number of reads
	TOTAL=$(wc -l whitelist.list | sed 's/ whitelist.list//')
	echo
	echo done. $TOTAL unique reads extracted.
	echo whitelist of unmapped subread IDs written to whitelist.list
	
else
        echo
        echo ERROR: unmapped subreads FASTA file does not exist in this directory
        echo navigate to directory containing unmappedSubreads.fasta and re-run script
        echo exiting
        echo 
fi

exit 0