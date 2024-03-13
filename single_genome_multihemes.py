#!/usr/bin/env python

import sys
import os														# 'rename' will be accessed via os.rename
import re
import glob
import fileinput
from Bio import SeqIO, BiopythonParserWarning, BiopythonWarning	 # Importing from Bio in one line
from Bio.Seq import UnknownSeq
from Bio.SeqUtils import molecular_weight
import argparse													# 'RawTextHelpFormatter' will be accessed via argparse.RawTextHelpFormatter
import warnings

from Bio import BiopythonWarning

# Suppress the BiopythonWarning about partial codons - may need to modify script to make all CDS be multiples of 3 in the future
warnings.simplefilter('ignore', BiopythonWarning)


# script help and usage

parser = argparse.ArgumentParser(
	description='This script parses a Genbank file and writes tab-delimited statistics of putative multiheme c-type '
				'cytochromes (3 or more CXXCH motifs; if this qualification is met, it also counts CXXXCH motifs).\n\n'
				'NOTE: Organisms with multiple Genbank records (e.g., those with multiple chromosomes or plasmids) '
				'should be concatenated into a single .gbk file before executing this script. For example:\n'
				'% cat NC_000001.gbk NC_000002.gbk [...NC_00000n.gbk] > concatenated.gbk\n',
	epilog='Modified and maintained by Bond Lab, University of Minnesota (https://bondlab.umn.edu)\n',
	formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('genbank_file', help='Genbank file, ideally containing translated amino acid sequences, but '
										 'translation will be attempted if not present. Requires SOURCE annotation '
										 'to be present in the header. For example:\n\n'
										 'LOCUS		  NC_002939			   3814128 bp	 DNA	 circular CON 16-MAY-2014\n'
										 'DEFINITION  Geobacter sulfurreducens PCA chromosome, complete genome.\n'
										 'ACCESSION	  NC_002939\nVERSION	 NC_002939.5  GI:400756305\n'
										 'DBLINK	  BioProject: PRJNA57743\nKEYWORDS	  RefSeq.\n'
										 'SOURCE	  Geobacter sulfurreducens PCA <--- *MUST BE PRESENT*\n	 '
										 'ORGANISM	Geobacter sulfurreducens PCA <--- *MUST BE PRESENT*')
args = parser.parse_args()


# open the file, and read the input file name for use in the output filename
genbankFile = open(sys.argv[1], 'r')
InputFileName = str(sys.argv[1])


# create empty lists of multiheme cytochromes and feature types
multihemeCytochromes = []
featureList = []
clusters = []
OutputString = ""


# parse input genbank file with SeqIO
for sequenceRecord in SeqIO.parse(genbankFile, "genbank"): 

	# define output filename as the organism (SOURCE from genbank file) and accession number
	# strip trailing parentheticals such "(complete genome)" as well as strain designations
	# replace spaces with underscore
	
	OrganismName = str(sequenceRecord.annotations["source"]).split(", ")[0].split(" (")[0].split(" DSM ")[0].split(" ATCC ")[0].replace(" ", "_").replace(".", "")
	AccessionNumber = str(sequenceRecord.id).strip("['], ")
	RecordName = str(sequenceRecord.id).strip("['], ")
	
	# append filename with .cytochromes.txt
	OutFileName = "%s.%s.cytochromes.txt" % (OrganismName, AccessionNumber)
	
	# define output file handle
	OutFile = open(OutFileName, 'w')
	
	
	# read each feature in the genbank file
	for feature in sequenceRecord.features:
	
		# capture all feature types and append to a long featureList
		featureList.append(str(feature.type))
		
		#takes elements one at a time, features are things like 'gene', then 'CDS', then 'misc_feature'
		
		# find unique feature types
		featureTypes = set(featureList)
		
		if feature.type == "CDS": 
			
		# Translate protein sequence on the fly if not already present
			if 'translation' not in feature.qualifiers:
				DNAseq = feature.extract(sequenceRecord.seq)
				
				# Use transl_table default of table 1
				# transl_table = feature.qualifiers.get("transl_table", "1").strip("[']")
				feature.qualifiers["translation"] = DNAseq.translate(table=1, to_stop=True)
			
			# identify all cannonical CXXCH motifs in protein coding sequences
			CXXCHmotifs = re.findall('C..CH', str(feature.qualifiers["translation"]))
			
			if len(CXXCHmotifs) >= 3:	# in this case, multiheme is defined as 3 or more hemes. Add it to the list and scan for extra features.
			
				# append the list of multiheme cytochromes each time another is identified
				multihemeCytochromes.append(feature)
				
				# identify all CXXXCH motifs in protein coding sequences if they already have at least 3 CXXCH motifs
				CXXXCHmotifs = re.findall('C...CH', str(feature.qualifiers["translation"]))
				hemeBindingMotifs = CXXCHmotifs + CXXXCHmotifs
				
				# identify all CXCH condensed motifs in protein coding sequences if they already have at least 3 CXXCH motifs
				CXCHmotifs = re.findall('C.CH', str(feature.qualifiers["translation"]))
				
				# identify all omcZ-like CX14CH motifs in protein coding sequences if they already have at least 3 CXXCH motifs
				CXXXXXXCHmotifs = re.findall('C...........[!=C]..CH', str(feature.qualifiers["translation"]))
				hemeBindingMotifs = CXCHmotifs + CXXCHmotifs + CXXXXXXCHmotifs + CXXXCHmotifs
				
				# calculate data about each cytochrome for final output
				
				# capture locus tag as string with [' locus tag '] characters removed
				locusTag = str(feature.qualifiers["locus_tag"]).strip("[']")
				
				# count number of motifs found in each protein coding sequence
				motifCount = len(hemeBindingMotifs)
				
				# count number of non-standard C--H motifs found in each protein coding sequence
				oddmotifCount = len(CXCHmotifs + CXXXXXXCHmotifs + CXXXCHmotifs)
				
				# capture translated amino acid sequence as string with [' AASEQ '] characters removed
				AAseq = str(feature.qualifiers["translation"]).strip("[']")
				
				# determine length of amino acid sequence **[' AASEQ '] characters must be removed for an accurate count!
				AAlength = len(str(feature.qualifiers["translation"]).strip("[']"))
				
				# if no ambiguous AAs present, calculate each cytochrome's molecular weight
				AmbiguousAA = re.findall('[BXZJUO]', str(feature.qualifiers["translation"]))
				
				if not AmbiguousAA:
				
					MolecularWeight = molecular_weight(AAseq, "protein")
				
				# calculate heme density as number of hemes per kDa
				HemeDensity = (float(motifCount) / MolecularWeight) * 1000
		
				if multihemeCytochromes:
					
					# add gene name to definition line if present in CDS feature qualifiers
					if 'gene' in feature.qualifiers:
					
						# define GeneName variable
						GeneName = str(feature.qualifiers["gene"]).strip("[']")
						
						# define output string in tab-delimited format if cytochromes were predicted
						OutputString = "%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%1.2f\t%1.3f\t%s\t%s\t%s" % (locusTag, str(feature.qualifiers["product"]).strip("[']"), GeneName, motifCount, len(CXXCHmotifs), len(CXXXCHmotifs), oddmotifCount, AAlength, float(MolecularWeight / 1000), HemeDensity, OrganismName.replace("_", " ").replace("sp ", "sp. "), RecordName, AAseq)

						# write output string to output file
						OutFile.write(OutputString)
						OutFile.write("\n")
				
					if not 'gene' in feature.qualifiers:  # in case it doesn't have a gene name, leave that field blank
							
						# define output string if cytochromes were predicted
						OutputString = "%s\t%s\t \t%i\t%i\t%i\t%i\t%i\t%1.2f\t%1.3f\t%s\t%s\t%s" % (locusTag, str(feature.qualifiers["product"]).strip("[']"), motifCount, len(CXXCHmotifs), len(CXXXCHmotifs), oddmotifCount, AAlength, float(MolecularWeight / 1000), HemeDensity, OrganismName.replace("_", " ").replace("sp ", "sp. "), RecordName, AAseq)
						
						# write output string to output file
						OutFile.write(OutputString)
						OutFile.write("\n")
						
	# errors

	# raise error message and exit script if no annotations present
	if 'CDS' not in featureTypes:
		print("%s\tERROR: Genbank file contains no annotations" % OrganismName.replace("_", " ").replace("sp ", "sp. "))
		os.remove(OutFileName)
		sys.exit(1)	 # Use sys.exit(1) to indicate error condition
	
	# raise error message and exit script if Genbank file contains no DNA sequence to translate from on the fly
	if type(sequenceRecord.seq) == UnknownSeq:
		print("%s\tERROR: Genbank file contains unknown DNA sequence" % OrganismName.replace("_", " ").replace("sp ", "sp. "))
		os.remove(OutFileName)
		sys.exit(1)	 # Use sys.exit(1) to indicate error condition



#
#
# done finding cytochromes
#
# report statistics and clean up
#
# if the genome was in multiple contigs, take all the files and concat into one, then add headers
#

if multihemeCytochromes:

	# simplify final output filename - strip off old extension
	RenamedOutFile = InputFileName[:-4] + ".cytochromes.txt"
	
	# Remove output file if it already exists from a previous execution of the script
	if os.path.isfile(RenamedOutFile):
		os.remove(RenamedOutFile)
	OutputFileList = glob.glob(""+OrganismName+"*.cytochromes.txt")
	
	# if multiple .txt files were generated, i.e. if genome exists in >1 contig, concatenate output files to clean up
	if len(OutputFileList) > 1:
		ConcatenatedOutFile = "%s.cytochromes.txt" % InputFileName[:-4]
						
		#check if ConcatenatedOutFile already exists and remove it
		if os.path.isfile(ConcatenatedOutFile):
			os.remove(ConcatenatedOutFile)
		
		# loop over OutputFileList if more than one .txt file generated, and concatenate .txt files, but ignore any files that are empty
		for txtFile in OutputFileList:
			if os.stat(txtFile).st_size > 0:
				
				# write a single concatenated .txt file 
				os.system("cat "+txtFile+" >> "+ConcatenatedOutFile+"")

			# remove individual txt files produced from each contig
			os.remove(txtFile)
	
		
		# display number of cytochromes identified for those watching at home and output filename
		
		print("%s\t%i cytochromes" % (OrganismName.replace("_", " ").replace("sp ", "sp. "), len(multihemeCytochromes)))
		
		
		# add column headers to concatenated output file 
		
		headers = 'LOCUS_TAG PRODUCT GENE HEMES CXXCH CXXXCH ODDMOTIF AA kDa HEMES/kDa ORGANISM ACCESSION SEQUENCE'.split()
		for line in fileinput.input([ConcatenatedOutFile], inplace=True):
			if fileinput.isfirstline():
				print('\t'.join(headers))
			print(line.rstrip()),				# removes extra linebreak caused by reprinting line
		
		# close input and output files
		genbankFile.close()
		OutFile.close()
	#
	#
	#	
	#  same final file handling, but executed only if genome exists as single contig
	#
	#
	
	if len(OutputFileList) == 1:

		# rename output file
		# rename output file
		os.rename(OutFileName, RenamedOutFile)


		# display number of cytochromes identified and output filename
		print("%s\t%i cytochromes" % (OrganismName.replace("_", " ").replace("sp ", "sp. "), len(multihemeCytochromes)))
		# print "amino acid sequences of multiheme cytochromes written to %s" % RenamedOutFile
		
		# add column headers to concatenated output file
		headers = 'LOCUS_TAG PRODUCT GENE HEMES CXXCH CXXXCH ODDMOTIF AA kDa HEMES/kDa ORGANISM ACCESSION SEQUENCE'.split()
		for line in fileinput.input([RenamedOutFile], inplace=True):
			if fileinput.isfirstline():
				print('\t'.join(headers))
			print(line.rstrip()),			# removes extra linebreak caused by reprinting line
			
		# close input and output files
		genbankFile.close()
		OutFile.close()
		
		
	
# report error if NO multiheme cytochromes were identified

if not multihemeCytochromes:
	print("%s\t%i cytochromes" % (OrganismName.replace("_", " ").replace("sp ", "sp. "),len(multihemeCytochromes)))
	
	#remove empty output file(s)
	EmptyFiles = glob.glob(""+OrganismName+"*.cytochromes.txt")
	for txtFile in EmptyFiles:
		if os.stat(txtFile).st_size == 0:
			os.remove(txtFile)
	
	genbankFile.close()

