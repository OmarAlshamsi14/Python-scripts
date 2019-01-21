#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Alphabet import IUPAC
import sys
import re
import argparse
from argparse import RawTextHelpFormatter

# This script is run in command line like so: > python pilin_stats.py 'fastafile.fa' , name of output file
        

# function that generates the sequence of the mature pilin protein
def mature_pilin_sequence(seq):
    m = re.search(r'[GAS][ACFGILMNPQSTVWY][ACFGILMNPQSTVWY][ACFGILMNPQSTVWY][ACFGILMNPQSTVWY][DE]', seq)
    if m:
        cleavage_motif = m.group()
        cleavage_pos = m.start()
        print('Found cleavage site at position %i with sequence %s' % (cleavage_pos, cleavage_motif))
        mature_seq = seq[cleavage_pos + 1:]
        return mature_seq

#function that counts the aromatic amino acids in the mautre pilin sequence
def count_aromatics(protein):
    tyrosines = protein.upper().count('Y')
    phenylalanines = protein.upper().count('F')
    tryptophans = protein.upper().count('W')
    histidines = protein.upper().count('H')
    aromatics = tyrosines + phenylalanines + tryptophans + histidines
    return aromatics, tyrosines, phenylalanines, histidines, tryptophans


#function that calculates the largest aromatic free gap in the amino acid sequence
def aromatic_gap(protein):
    gaps = re.split(r'[YFWH]', protein.upper())
    lengths = []
    for gap in gaps:
        gap_length = len(gap)
        lengths.append(gap_length)
    lengths.sort()
    largest_gap = lengths[-1]
    return largest_gap

    

# script help and usage
parser=argparse.ArgumentParser(
    description='This script takes an amino acid fasta file of putative type IV pilin proteins and writes a tab-delimited file containing the mature pilin sequence, the mature length, and the number and percentage of aromatic amino acids within the mature length of the pilin. \nRequires BioPython v. 1.65 or later (http://biopython.org/wiki/Download)', formatter_class=RawTextHelpFormatter)
parser.add_argument('inputfile.fa', help='Amino acid fasta file containing putative type IV pilin proteins')
parser.add_argument('outfile',  help='Output file containing your pilin stats in a tab-delimited table')
args=parser.parse_args()


#define the input and output files
Inputfile = sys.argv[1]
Output_name = sys.argv[2]

#open outfile to write
output = open(Output_name + '.txt', 'w')
headers = 'SequenceID\tAromatic_count\tmature_pilin length\tpercent_aromatics\tlargest gap\ttyrosines\tphenylalanines\thistidines\ttryptophans'
output.write(headers + '\n')
# open our fasta input file

epili = []
with open(Inputfile, 'rU') as f:
        # iterate over fasta file and parse out sequences
	for record in SeqIO.parse(f, "fasta", IUPAC.protein ):
            print('Analyzing sequence: ' + record.id)
            
            # Generate the full sequences
            full_sequence = (record.seq).upper().rstrip('\n')
            
            #generate the mature sequences
            mature_sequence = mature_pilin_sequence(str(full_sequence))
            
            print('Mature seuqencs is: ' + mature_sequence)
            
            #count the aromatics of the mature sequence
            aromatics = count_aromatics(mature_sequence)[0]
            
            # count number of tyrosines
            tyrosines = count_aromatics(mature_sequence)[1]
            
            #count the number of phenylalanines
            phenylalanines = count_aromatics(mature_sequence)[2]
            
            #count the number of histidines
            histidines = count_aromatics(mature_sequence)[3]
            
            #count the number of tryptophans
            tryptophans = count_aromatics(mature_sequence)[4]
            
            #determine the length of the full sequence
            mature_length = len(mature_sequence)
            
            #calculate the percent of aromatic amino acids in the mature sequence
            aromatic_percent = (aromatics / mature_length) * 100
            
            #calculates largest aromatic free gap
            free_gap = aromatic_gap(mature_sequence)
            
            print('Aromatic percentage is: ' + str(round(aromatic_percent, 2)))
            
            # writes out put in tab delimited format, make sure order of vaules here matches order of headers variable deifined above
            output.write('%s\t%i\t%i\t%.2f\t%i\t%i\t%i\t%i\t%i\t%s' % (record.id, aromatics, mature_length, aromatic_percent, free_gap, tyrosines, phenylalanines, histidines, tryptophans, record.seq) + '\n')
            if aromatic_percent > 9.80:
                epili.append(record)
SeqIO.write(epili, Output_name + "_epili.fasta", "fasta")
print(epili)
#close output file
output.close()

