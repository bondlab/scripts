from sys import argv
import sys
import os
from os import rename
import re
import glob
import fileinput
from Bio.Seq import UnknownSeq
from Bio.SeqUtils import molecular_weight
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO

script, filename = argv # not sure what this does but I need it

txt = open(filename)

print "Opening your file %r:" % filename # check I opened the right file


all_hemes = []
hemestotal = 0
records = list(SeqIO.parse(filename, "fasta")) # counts number of > records

for seq_record in SeqIO.parse(filename, "fasta"):
	hemes = "blank"
	count = 0 	# resets the count and the fields
	while hemes != "heme":
		hemes = seq_record.description.split()[count] # each description is numbered, spaces separate
		count= count + 1	
	all_hemes.append(seq_record.description.split()[count - 2]) # keeps list of hemes, counting is hard
	hemestotal = hemestotal + float(seq_record.description.split()[count - 2]) # counts hemes as integer
print "List of hemes", all_hemes

average = hemestotal / (len(records))

print "average hemes", average

# This is going to get ugly, but keep going using the same pattern to find kDa
# It gets twitchy when the values are decimals

all_kda = []
kdatotal = 0
records = list(SeqIO.parse(filename, "fasta")) # counts number of > records

for seq_record in SeqIO.parse(filename, "fasta"):
	all_kda = "blank"
	count = 0 	# resets the count and the fields
	while all_kda != "kDa,":
		all_kda = seq_record.description.split()[count] # changes the variable to next word
		count= count + 1	
		# punted on figuring on how to make an array of kDa
	kdatotal = kdatotal + float(seq_record.description.split()[count - 2]) # counts kDa as integer
	
print "total kDa", kdatotal

averagekda = kdatotal / (len(records))
print "average size", averagekda, "kDa"

finalratio = averagekda / average

print "kDa:heme=", finalratio # the ratio of kDa to hemes

# write output string to output file as I figure this out

OutFileName = filename + "." + str(average) + ".processed.txt"
OutFile = open(OutFileName, 'w')
OutFile.write(filename)
OutFile.write(" ")
OutFile.write(str(average))
OutFile.write(" ")
OutFile.write(str(averagekda))
OutFile.write(" ")
OutFile.write(str(finalratio))
OutFile.write("\n")

# Something is wrong with the average heme calc - it seems to take the lowest or round down?

#for f in $(ls /folder/folder/*.faa); do python fasta.py $f;done

OutFile.close()


print "Done"