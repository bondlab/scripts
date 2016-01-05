#!/usr/bin/env python

import sys
import os
from os import rename
import re
import glob
import fileinput
from Bio import SeqIO
from Bio.Seq import UnknownSeq
from Bio.SeqUtils import molecular_weight
import argparse
from argparse import RawTextHelpFormatter

# silence Biopython warnings of improper indentation of Genbank files and CDS that are not multiples of 3
import warnings
from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonWarning)

# script help and usage
parser=argparse.ArgumentParser(
    description='This script parses a Genbank file and writes tab-delimited statistics of putative multiheme c-type\ncytochromes (3 or more CXXCH motifs; if this qualification met, counts also CXXXCH motifs)\n\nNOTE: Organisms with multiple Genbank records (e.g. those with multiple chromosomes or plasmids)\nshould be concatenated into a single .gbk file before executing this script. For example:\n% cat NC_000001.gbk NC_000002.gbk [...NC_00000n.gbk] > concatenated.gbk\n\nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)', 
    epilog='Author: Jon Badalamenti, Bond Lab, University of Minnesota (http://www.thebondlab.org)\nJune 2015\n \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('[GENBANK FILE]', help='Genbank file containing translated amino acid sequences. Requires SOURCE annotation\nto be present in the header. For example:\n \nLOCUS       NC_002939            3814128 bp    DNA     circular CON 16-MAY-2014\nDEFINITION  Geobacter sulfurreducens PCA chromosome, complete genome.\nACCESSION   NC_002939\nVERSION     NC_002939.5  GI:400756305\nDBLINK      BioProject: PRJNA57743\nKEYWORDS    RefSeq.\nSOURCE      Geobacter sulfurreducens PCA <--- *MUST BE PRESENT* \n  ORGANISM  Geobacter sulfurreducens PCA <--- *MUST BE PRESENT*')
args=parser.parse_args()

# define input file handle
genbankFile = open(sys.argv[1], 'r')

# create empty lists of multiheme cytochromes and feature types
multihemeCytochromes = []
featureList = []

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
	
	# loop over features in genbank file
	for feature in sequenceRecord.features:
	
		# capture all feature types and append to featureList
		featureList.append(str(feature.type))
		
		# find unique feature types
		featureTypes = set(featureList)
		
		if feature.type == "CDS":
			
			# translate protein sequence on the fly if not already present
			if 'translation' not in feature.qualifiers:
				DNAseq = feature.extract(sequenceRecord.seq)
				feature.qualifiers["translation"] = DNAseq.translate(table=str(feature.qualifiers["transl_table"]).strip("[']"), to_stop=True)
			
			# identify all CXXCH motifs in protein coding sequences
			CXXCHmotifs = re.findall('C..CH', str(feature.qualifiers["translation"]))
			
			if len(CXXCHmotifs) >= 3:
				# append the list of multiheme cytochromes each time another is identified
				multihemeCytochromes.append(feature)
				
				# identify all CXXXCH motifs in protein coding sequences if they already have at least 3 CXXCH motifs
				CXXXCHmotifs = re.findall('C...CH', str(feature.qualifiers["translation"]))
				hemeBindingMotifs = CXXCHmotifs + CXXXCHmotifs
				
				# define variables for printing output
				# cast locus tag as string with [' locus tag '] characters removed
				locusTag = str(feature.qualifiers["locus_tag"]).strip("[']")
				
				# count number of CXXCH motifs found in each protein coding sequence
				motifCount = len(hemeBindingMotifs)
				
				# cast translated amino acid sequence as string with [' AASEQ '] characters removed
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
					
					# add gene name to FASTA definition line if present in CDS feature qualifiers
					if 'gene' in feature.qualifiers:
					
						# define GeneName variable
						GeneName = str(feature.qualifiers["gene"]).strip("[']")
						
						# define output string in FASTA format if cytochromes were predicted
						OutputString = "%s\t%s\t%s\t%i\t%i\t%i\t%i\t%1.2f\t%1.3f\t%s\t%s\t%s" % (locusTag, str(feature.qualifiers["product"]).strip("[']"), GeneName, motifCount, len(CXXCHmotifs), len(CXXXCHmotifs), AAlength, float(MolecularWeight / 1000), HemeDensity, OrganismName.replace("_", " ").replace("sp ", "sp. "), RecordName, AAseq)
				
						# write output string to output file
						OutFile.write(OutputString + "\n")
				
					if not 'gene' in feature.qualifiers:
							
						# define output string in FASTA format if cytochromes were predicted
						OutputString = "%s\t%s\t \t%i\t%i\t%i\t%i\t%1.2f\t%1.3f\t%s\t%s\t%s" % (locusTag, str(feature.qualifiers["product"]).strip("[']"), motifCount, len(CXXCHmotifs), len(CXXXCHmotifs), AAlength, float(MolecularWeight / 1000), HemeDensity, OrganismName.replace("_", " ").replace("sp ", "sp. "), RecordName, AAseq)
				
						# write output string to output file
						OutFile.write(OutputString + "\n")

	# raise error message and exit script if no annotations present
	if 'CDS' not in featureTypes:
		print "%s\tERROR: Genbank file contains no annotation" % OrganismName.replace("_", " ").replace("sp ", "sp. ")
		os.remove(OutFileName)
		exit()
	
	# raise error message and exit script if Genbank file contains no DNA sequence to translate from on the fly
	if type(sequenceRecord.seq) == UnknownSeq:
		print "%s\tERROR: Genbank file contains unknown DNA sequence" % OrganismName.replace("_", " ").replace("sp ", "sp. ")
		os.remove(OutFileName)
		exit()

# report statistics if multiheme cytochromes were identified
if multihemeCytochromes:

	# simplify final output filename
	RenamedOutFile = OrganismName + ".cytochromes.txt"
	
	# Remove output file if it already exists from a previous execution of the script
	if os.path.isfile(RenamedOutFile):
		os.remove(RenamedOutFile)
	OutputFileList = glob.glob(""+OrganismName+"*.cytochromes.txt")
	
	# if multiple .txt files were generated, i.e. if genome exists in >1 contig, concatenate output files
	if len(OutputFileList) > 1:
		ConcatenatedOutFile = "%s.cytochromes.concatenated.txt" % OrganismName
						
		#check if ConcatenatedOutFile already exists and remove it
		if os.path.isfile(ConcatenatedOutFile):
			os.remove(ConcatenatedOutFile)
		
		# loop over OutputFileList if more than one .txt file generated, and concatenate .txt files only if not they are empty
		for txtFile in OutputFileList:
			if os.stat(txtFile).st_size > 0:
				
				# write concatenated .txt file	
				os.system("cat "+txtFile+" >> "+ConcatenatedOutFile+"")

			# remove individual txt files produced from each contig
			os.remove(txtFile)
		
		# display number of cytochromes identified and output filename
		print "%s\t%i cytochromes" % (OrganismName.replace("_", " ").replace("sp ", "sp. "), len(multihemeCytochromes))
		# print "amino acid sequences of multiheme cytochromes written to %s" % ConcatenatedOutFile
		
		# add column headers to concatenated output file
		headers = 'LOCUS_TAG PRODUCT GENE HEMES CXXCH CXXXCH AA kDa HEMES/kDa ORGANISM ACCESSION SEQUENCE'.split()
		for line in fileinput.input([ConcatenatedOutFile], inplace=True):
			if fileinput.isfirstline():
				print '\t'.join(headers)
			print line,
		
		# close input and output files
		genbankFile.close()
		OutFile.close()
		
	# execute only if genome exists as single contig
	if len(OutputFileList) == 1:

		# rename output file
		rename(OutFileName, RenamedOutFile)

		# display number of cytochromes identified and output filename
		print "%s\t%i cytochromes" % (OrganismName.replace("_", " ").replace("sp ", "sp. "), len(multihemeCytochromes))
		# print "amino acid sequences of multiheme cytochromes written to %s" % RenamedOutFile
		
		# add column headers to concatenated output file
		headers = 'LOCUS_TAG PRODUCT GENE HEMES CXXCH CXXXCH AA kDa HEMES/kDa ORGANISM ACCESSION SEQUENCE'.split()
		for line in fileinput.input([RenamedOutFile], inplace=True):
			if fileinput.isfirstline():
				print '\t'.join(headers)
			print line,
			
		# close input and output files
		genbankFile.close()
		OutFile.close()
	
# report error if NO multiheme cytochromes were identified
if not multihemeCytochromes:
	print "%s\t%i cytochromes" % (OrganismName.replace("_", " ").replace("sp ", "sp. "),len(multihemeCytochromes))
	
	#remove empty output file(s)
	EmptyFiles = glob.glob(""+OrganismName+"*.cytochromes.txt")
	for txtFile in EmptyFiles:
		if os.stat(txtFile).st_size == 0:
			os.remove(txtFile)
	
	genbankFile.close()