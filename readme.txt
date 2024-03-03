Welcome to the GeneSweeper!
Written by: Abdulaziz Alajlan
dump all CSV files in the same folder as the program.py file in order for code to run properly
Sometimes entering product name gives incorrect output. if entering the full name of the product doesn't work enter a uniqe snippet of the product name
The link file is used to get the meta gene details from the JGI website so that a fasta file can be created
For the Links output file the result might be different depending on your product. Use the Link_Gen.py file to find out which part of the line variable is the gene length so you can get the correct link files

Libraries used:
tkinter
from tkinter, simpledialog
from tkinter, messagebox
os
csv

Functions:

renamer() [No parameters]:
Renaming files from 1 to the number of files available.

filter_and_append_to_csv(Input_file, output_file, product):
Creates an empty CSV file, gives appropriate headers, runs through all numbered files and appends them to csv files on the condition that the product name is the same as the inputted product name given by the user.

linkgen()[No parameters]:
Creates an empty CSV file, gives appropriate headers, adds links to JGI pages by putting gene and scaffold IDs in the link address.

delete_csv_files()[No parameters]:
Deletes all CSV files besides output and link files.

renamer() [Not used]:
Renames output files to have the product name be incorporated to file name.

fileautomation():
Runs all functions while creating the file address for the filter_and_append_to_csv() parameters.

Process:

Dump all CSV files into the program folder.
Run the GeneSweeper to get the output file and a links file.
Open each link in a links file with gene sequences length of longer than 1000, copy the amino acid sequence from the JGI website, and paste it into a fasta file.
