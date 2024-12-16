#!/usr/bin/python3
# -*- coding : utf-8 -*-

#licence auteur et version

__authors__ = "Loïk Galtier"
__Version__ = "1.0.0"
__Licence__ = ("This program is free software: you can redistribute it and/or modify"
               " it under the terms of the GNU General Public License as published by)"
               " the Free Software Foundation, either version 3 of the License, or "
               "(at your option) any later version. "
               "This program is distributed in the hope that it will be useful,"
               " but WITHOUT ANY WARRANTY; without even the implied warranty of"
               " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the "
               " GNU General Public License for more details. "
               "You should have received a copy of the GNU General Public License"
               " along with this program. If not, see <https://www.gnu.org/licenses/>.")

import os
import re
import csv
import sys

import pandas as pd

####### Help #######
def show_help():
	print("To run this script:")
	print("SamMapScript.py \t will ask some questions to launch the analysis")
	print("SamMapScript.py -h \t will display help")
	print("SamMapScript.py [file] \t will start the analysis of the file")
	print("SamMapScript.py [file] [option...] \t will start the analysis of the file with the following options:\n\n")
	
	print("--[Options]-- \n")
	print("quality[int] \t sets the minimum read quality to analyze to 'int'\t(reads with quality < quality_min are ignored as if they never existed)")
	print("step_qual[int] \t defines the step size for quality intervals as 'int'")
	print("specific_sequence[seq] \t restricts the analysis to the reference sequence 'seq'; to specify multiple sequences: specific_sequence[seq1,seq2] (be careful with spaces)\t(reads with other sequence are ignored as if they never existed) \n")
	
	print("To run specific steps only:")
	print("all\t Every possible action.")
	print("statistic_map\t Provides the quantities and percentages of mapped reads, semi-mapped reads (mapped but with a CIGAR containing other letters than 'M'), and unmapped reads, according to the reference sequences")
	print("analyse_flag\t Provides a human readable analysis of the flags.")
	print("count_flag\t\t Displays the counts of flag bits across all sequences.")
	print("analyse_quality\t Shows the distribution of quality scores in intervals with a default step size of 10 (or adjustable using step_qual[]).")
	print("analyse_cigar\t Presents the counts of CIGAR letters across all sequences combined. \n")
	
	print("exemple : /home/user/Downloads/SamMapScript.py /home/user/Downloads/mapping.sam quality[5] step_qual[15] specific_sequence[reference1,reference2] statistic_map count_flag")
	sys.exit()


######## Setting the file ########
# check the selected file exist and have the .sam
def find_file():
	if len(sys.argv) == 1:                                                  # If only the SamMapScript.py argument
		no_file_found = True                                                # Set a boolean for the while after
		while no_file_found:                                                #
			path = input("What is the path to the SAM file : ?")            # Ask the path
			if os.path.exists(path) and re.search(r".sam", path):    # If the path exist AND the .sam extension is here (frm the path in the input)
				no_file_found = False                                       # Quit the while
				print("A SAM file has been found")                          # A little message
				return path                                                 # Save the path to the file
			else:
				print("No SAM file found, please try again \n")             # Ask again for a file.
				
	else:                                                                   # Else if more argument
		if os.path.exists(sys.argv[1]) and re.search(r".sam", sys.argv[1]):  # if the path exist AND the .sam extension is here (from the path pu in argument)
			print("A SAM file has been found")                              # A little message
			return sys.argv[1]                                              # Save the path to the file
		else:
			print("Error: No SAM file has been found \n")                   # No sam
			sys.exit()                                                      # Quit the program

#### Define quality ####
def ask_quality():
	if len(sys.argv) == 1:                                                  # If only the SamMapScript.py argument
		no_answer = True                                                    # Set a boolean for the while after
		while no_answer:                                                    #
			need_quality_min = input("Do you want to apply a minimum quality filter for the rest of the processing: (Y/N)").upper()  # Ask if needed a quality minimal
			if need_quality_min == "Y" or "y" or "yes":                     # If user say yea
				quality_min_input = input("What minimum value did you want to apply : [0-255] \n")  # Ask the quality filter
				if quality_min_input.isdigit():                             # If is a number
					print("All operations will be performed on reads with a quality greater than: " + quality_min_input + "\n")
					no_answer = False                                       # exit while (needed ?)
					return int(quality_min_input)                           # save the quality minimal
				else:                                                       # if not a number
					print("This is not a correct value \n")                 # A lttle message and ask again
			elif need_quality_min == "N" or "n" or "no":                    # If user say nah
				no_answer = False                                           # exit while
				print("\n")                                                 # message for beauty
				return 0                                                    # quality min equal 0 so no limit (except negative)
			else:                                                           # If user say "don't know" or "hello there"
				print("This is not a correct value \n")                     # not correct value, try again
	else:
		return 0                                                            # If more argument (but not quality_min[]) set no limit (0)

######## Extract part and use header ########
def sam_reading(sam_file_path):
	with open(sam_file_path, "r") as sam_file:                              # Open the sam file wich Reading text
		
		sequence_names = []                                                 # Create a list for putting name of reference sequence dictionary
		flags = []                                                          # Create a list for putting the flags
		rnames = []                                                         # Create a list for putting reference name
		mapqs = []                                                          # Create a list for putting quality
		cigars = []                                                         # Create a list for putting Cigar
		
		number_of_read_total = 0                                            # Create a variable for counting the total of read
		
		hd_table = {}                                                       # Create a dictionary for the header part of the file
		sq_table = {}                                                       # Create a dictionary for the sequence part of the file
		pg_table = {}                                                       # Create a dictionary for the program part of the file
		
		
		for line in sam_file:                                               # for every line
			if line.startswith("@"):                                        # If is a header line (starting with @)
				if line.startswith("@HD"):                                  # If it is the header line
					vN = re.search(r"VN:([^\t|\n]+)", line)          # Find the motif : "VN:([^\t|\n]+)" for the version of format
					sO = re.search(r"SO:([^\t|\n]+)", line)          # Find the motif : "SO:([^\t|\n]+)" for the order of alignment (unknown,unsorted,queryname,coordinate)
					if vN:                                                  #
						hd_table[vN.group(1)] = ""                          #
						if sO:                                              #
							hd_table[vN.group(1)]= [sO.group(1)]            #
				
				elif line.startswith("@PG"):                                # If it is the program line
					idnamepg = re.search(r"ID:([^\t|\n]+)", line)    # Find the motif : "ID:([^\t|\n]+)" for the id of the program
					versionpg = re.search(r"VN:([^\t|\n]+)", line)   # Find the motif : "VN:([^\t|\n]+)" for the version of the program
					if idnamepg:
						pg_table[idnamepg.group(1)] = ""
						if versionpg:
							pg_table[idnamepg.group(1)] = versionpg.group(1)
						
				elif line.startswith("@SQ"):                                # If it is the sequence line
					sN = re.search(r"SN:([^\t|\n]+)", line)          # Find the motif : "SN:([^\t|\n]+)" for the name of the sequence
					lN = re.search(r"LN:([^\t|\n]+)", line)         # Find the motif : "LN:([^\t|\n]+)" for the lenght of the sequence
					if sN and lN:
						sequence_names.append(sN.group(1))
						sq_table[sN.group(1)] = lN.group(1)
						
				elif line.startswith("@RG"):                                # If it is the Read group line
					idgroup = re.search(r"ID:([^\t|\n]+)", line)
					if idgroup:
						print("A reading group has been formed, its unique ID is: " + idgroup.group(1) + "\n")
						save_file.write("A reading group has been formed, its unique ID is: " + idgroup.group(1) + "\n")
		sam_file.seek(0)                                                    # For stopping the loop
		
		print("_________________ \n Header :\n")
		save_file.write("_________________ \n Header :\n")
		
		if hd_table:
			print("The data in this file is :")
			save_file.write("The data in this file is :\n")
			data = [(f"{cle} |", f"{val} |") for cle, val in hd_table.items()]
			t_data = pd.DataFrame(data, columns=["Format Version |", "Order |"])
			print(t_data, "\n")
			save_file.write(f"{t_data.to_string(index=False)} \n\n")
			
		if sq_table:
			print("Reference sequences were found :")
			save_file.write("Reference sequences were found : \n")
			data = [(f"{cle} |", f"{val} |") for cle, val in sq_table.items()]
			t_data = pd.DataFrame(data, columns=["sequence name |", "Length |"])
			print(t_data, "\n")
			save_file.write(f"{t_data.to_string(index=False)} \n\n")
			
		
		if pg_table:
			print("Programs were used in it :")
			save_file.write("Programs were used in it :\n")
			data = [(f"{cle} |", f"{val} |") for cle, val in pg_table.items()]
			t_data = pd.DataFrame(data, columns=["Programme id |", "Version |"])
			print(t_data, "\n")
			save_file.write(f"{t_data.to_string(index=False)} \n\n")
			
		sam_reader = csv.reader(sam_file, delimiter='\t')           # Reads the SAM file, splitting each line by tabs
		for row in sam_reader:                                      # Iterates over each row in the SAM file
			if row[0].startswith("@"):                              # Skips header lines that start with "@"
				continue
			if quality_min <= int(row[4]) < 255:                    # Filters alignments based on a minimum quality and excludes more than 254
				if (len(specific_sequence) > 0 and row[2] in specific_sequence) or len(specific_sequence) == 0:  # Checks if the aligned sequence is a specific sequence or if no specific is provided
					if len(row) >= 5:                               # Check if all the part exist
						flag = int(row[1])                          # Extracts the FLAG value from the alignment
						flags.append(flag)
						
						rname = row[2]                              # Extracts the reference sequence name
						rnames.append(rname)
						
						mapq = row[4]                               # Extracts the mapping quality
						mapqs.append(mapq)
						
						cigar = row[5]                              # Extracts the CIGAR
						cigars.append(cigar)
						
						number_of_read_total += 1                   # Increments the counter for the total number of alignments
		if len(specific_sequence) > 0:                              # If a specific sequences:
			return specific_sequence, flags, rnames, mapqs, cigars, number_of_read_total
		else:
			return sequence_names, flags, rnames, mapqs, cigars, number_of_read_total


######## Convert to binary ########
def flags_to_binary(flags):
	for i in range(len(flags)):                                     # Loop that will go from 0 to the size of the flags
		flags[i] = bin(flags[i])                                    # Convert to binary
		flags[i] = flags[i][2:].zfill(12)                           # Remove first 2 characters of flag (0 and b) & fill with 0 on left up to 12 # if more than 12, will pass with
	return flags                                                    # Returning flags in binary

######## Statistic of alignment ########
def statistic_of_mapped_reads(binary_flag):
	nb_mapped_table = {}
	nb_semi_mapped_table = {}
	nb_unmapped = 0

	# mapped and semi-mapped part
	for sqname in specific_sequence:
		for i in range(len(binary_flag)):                           # Loop that iterates from 0 to the amount of reads
			if rname[i] == sqname:                                  # check if read is mapped to reference sequence
				flag = binary_flag[i]                               # Put a binary flag in the flag variable
				if flag[-2] == "1" and flag[-3] == "0":             # Bit 2 code for "aligned" info. If "1" = aligned. Bit 3 code for "unmapped"
					if not re.fullmatch(r"(\d*M+)+", cigars[i]):  # If there are not only M, then there is an alignment problem.
						if sqname in nb_semi_mapped_table:          # If sqname exists in the semi-mapped table
							nb_semi_mapped_table[sqname] += 1       # Add 1 to count the number of reads with sqname key in nb_semi_mapped_table
						else:                                       # Otherwise initialize it
							nb_semi_mapped_table[sqname] = 1
					else:                                           # If there is no problem (only M)
						if sqname in nb_mapped_table:               # If sqname exists in the mapped table
							nb_mapped_table[sqname] += 1            # Add 1 to count the number of reads with sqname key in nb_mapped_table
						else:                                       # Otherwise initialize it
							nb_mapped_table[sqname] = 1
						

	
	# unmapped part
	for i in range(len(binary_flag)):                               # Loop that iterates from 0 to the amount of reads
		flag = binary_flag[i]                                       # Put a binary flag in the flag variable
		if flag[-3] == "1" or rname[i] == "*" or rname[i] not in specific_sequence:  # bit 3 codes for "read unmapped" info. If "1" = unmapped. If rname[i] is "*" or not in reference, the read is not mapped
			nb_unmapped += 1                                        # Add 1 for the number of unmapped
	
	print(f"_________________ \n Alignment Statistic: With minimum quality of : {quality_min}\n")
	save_file.write(f"_________________ \n Alignment Statistic: With minimum quality of : {quality_min}\n\n")
	print(f"Amongs {totalNumberOfRead} reads :")
	save_file.write(f"Amongs {totalNumberOfRead} reads : \n")
	
	if totalNumberOfRead != 0:                                      # If the total is not 0 for the division ( and the analyse is not empty)
		data1 = [(f"{cle} |", f"{val} |", f" {round(val / totalNumberOfRead, 4) * 100} |", f" mapped |") for cle, val in nb_mapped_table.items()]  #
		t_data1 = pd.DataFrame(data1, columns=[" Sequence reference |", "Quantity |", " % |", " mapped ?|"])
		t_data = pd.concat([t_data1], axis=0)
		
		data2 = [(f"{cle} |", f"{val} |", f" {round(val / totalNumberOfRead, 4) * 100} |", f" Semi-mapped |") for cle, val in nb_semi_mapped_table.items()]  #
		t_data2 = pd.DataFrame(data2, columns=[" Sequence reference |", "Quantity |", " % |", f" mapped ?|"])
		t_data = pd.concat([t_data1, t_data2], axis=0)
		
		if nb_unmapped > 0:
			data3 = [(" Non aligné |", f"{nb_unmapped} |", f" {round(nb_unmapped / totalNumberOfRead, 4) * 100} |", f" not mapped |")]  #
			t_data3 = pd.DataFrame(data3, columns=[" Sequence reference |", "Quantity |", " % |", f" mapped ?|"])
			if nb_semi_mapped_table:
				t_data = pd.concat([t_data1, t_data2, t_data3], axis=0)
			elif nb_mapped_table:
				t_data = pd.concat([t_data1, t_data3], axis=0)
		print(t_data.to_string(index=False), "\n")
		save_file.write(f"{t_data.to_string(index=False)} \n\n")

def dico_count_flag():
	
	flags_count = {}
	
	for i in range(len(binary_flags)):                              # Loop that iterates from 0 to the amount of reads
		flag = binary_flags[i]                                      # Put a binary flag in the flag variable
		
		for j in range(len(flag)):                                  # Loop that iterate from 0 to the length of the flag (12 bits min here)
			if flag[-j - 1] == "1":
				if (str(2 ** j) + " bits") not in flags_count:
					flags_count[str(2 ** j) + " bits"] = 1
				else:
					flags_count[str(2 ** j) + " bits"] += 1
			else:
				if (str(2 ** j) + " bits") not in flags_count:
					flags_count[str(2 ** j) + " bits"] = 0
	return flags_count
	
def analysis_flag():
	print(f"_________________ \nFlag analysis: with minimum quality {quality_min} and on reference sequences: {specific_sequence}\n")
	save_file.write(f"_________________ \nFlag analysis: with minimum quality {quality_min} and on reference sequences: {specific_sequence}\n\n")
	for i in range(len(flags_counts)):
		if i == 0:
			print(f"There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) of the reads are paired and so {100 - round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}% are independent. \n")
			save_file.write(f"There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) of the reads are paired and so {100 - round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}% are independent. \n\n")
		if i == 2:
			if flags_counts[str(2 ** i) + " bits"] == flags_counts[str(2 ** (i+1)) + " bits"]:
				print(f"{flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) reads and their associated reads are not sequence aligned. \n")
				save_file.write(f"{flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) reads and their associated reads are not sequence aligned. \n\n")
			else:
				print(f"Some reads are aligned but not their complement. \n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) unaligned reads and {flags_counts[str(2 ** (1+i)) + ' bits']} ({round(flags_counts[str(2 ** (1+i)) + ' bits'] / totalNumberOfRead * 100, 2)}%) associated unaligned reads.\n")
				save_file.write(f"Some reads are aligned but not their complement.\n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) unaligned reads and {flags_counts[str(2 ** (1+i)) + ' bits']} ({round(flags_counts[str(2 ** (1+i)) + ' bits'] / totalNumberOfRead * 100, 2)}%) associated unaligned reads.\n\n")
		
		if i == 4:
			if flags_counts[str(2 ** i) + " bits"] == flags_counts[str(2 ** (i+1)) + " bits"]:
				print(f"All reads and their associates are in the same orientation.\n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) in the complementary reverse direction and as many in the non-reverse direction.\n")
				save_file.write(f"All reads and their associates are in the same orientation.\n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) in the complementary reverse direction and as many in the non-reverse direction.\n\n")
			else:
				print(f"There are more reads in one direction than the other!\n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) in the reverse complementary direction but {flags_counts[str(2 ** (1+i)) + ' bits']} ({round(flags_counts[str(2 ** (i+1)) + ' bits'] / totalNumberOfRead * 100, 2)}%) in the non-reverse direction.\n")
				save_file.write(f"There are more reads in one direction than the other!\n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 2)}%) in the reverse complementary direction but {flags_counts[str(2 ** (1+i)) + ' bits']} ({round(flags_counts[str(2 ** (i+1)) + ' bits'] / totalNumberOfRead * 100, 2)}%) in the non-reverse direction.\n\n")
				
		if i == 6:
			if flags_counts[str(2 ** i) + " bits"] == flags_counts[str(2 ** (i + 1)) + " bits"]:
				print(f"All reads are in pairs.\n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 1)}%) in the first segment and as many as in the second.\n")
				save_file.write(f"All reads are in pairs.\n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 1)}%) in the first segment and as many as in the second.\n\n")
			
			else:
				print(f"There is a problem in the file with the parameters entered.\n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 1)}%) in first segment and {flags_counts[str(2 ** (1+i)) + ' bits']} ({round(flags_counts[str(2 ** (i+1)) + ' bits'] /totalNumberOfRead * 100, 2)}%) in seconds.\n")
				save_file.write(f"There is a problem in the file with the parameters entered.\n There are {flags_counts[str(2 ** i) + ' bits']} ({round(flags_counts[str(2 ** i) + ' bits'] / totalNumberOfRead * 100, 1)}%) in first segment and {flags_counts[str(2 ** (1+i)) + ' bits']} ({round(flags_counts[str(2 ** (i+1)) + ' bits'] /totalNumberOfRead * 100, 2)}%) in seconds.\n\n")

####### Bonus : It happens that if parameter "detail flag" is called, some value is hard-coded here, but this def is not included in the use of the initial project ########
def count_flag_number():

	print(f"_________________ \nQuantity of flags: with minimum quality of: {quality_min}\n")
	save_file.write(f"_________________ \nQuantity of flags: with minimum quality of: {quality_min}\n\n")
	pd.set_option('display.max_colwidth', 40)
	
	if totalNumberOfRead != 0:
		data1 = [(f"{cle} |",f"{val} |", f" {round(val / totalNumberOfRead * 100, 2)} |") for cle, val in flags_counts.items()]  #
		t_data1 = pd.DataFrame(data1, columns=["Bits |","Quantity |", " % |"])
		
		t_data = pd.concat([t_data1], axis=1)
		
		print(t_data.to_string(index=False))
		save_file.write(f"{t_data.to_string(index=False)} \n\n")
		
#### ####
def quality_analysis():
	print(f"\n_________________ \nQuality distribution for sequence: {specific_sequence}, with a minimum quality of {quality_min}, and in steps of {quality_step}\n")
	save_file.write(f"\n_________________ \nQuality distribution for sequence: {specific_sequence}, with a minimum quality of {quality_min}, and in steps of {quality_step}\n\n")
	qual_dico_sequence={}
	
	for sqname in specific_sequence:
		qual_dico = {
			'0 <= quality':0,
		}
		
		for i in range(len(binary_flags)):          # Loop that iterates from 0 to the amount of reads
			if rname[i] == sqname:                  # check if read is mapped to reference sequence
				if mapqs[i].isdigit() and mapqs[i] != '0':  # check if mapQ is a number
					j = 0
					while j < 156:
						tmp = str(j) + " < quality <= " + str(j + quality_step)
						if j < int(mapqs[i]) <= j + quality_step:
							if tmp in qual_dico:
								qual_dico[tmp] += 1  # Add 1 to count the number of reads
							else:
								qual_dico[tmp] = 1
							j = 156                 # To stop the loop
						else:
							if tmp not in qual_dico:
								qual_dico[tmp] = 0  # allows to initialize all values (keys), up to the largest
						j += quality_step
				else:
					qual_dico['0 <= quality'] += 1  # Add 1 to count the number of reads
		
		qual_dico_sequence[sqname] = qual_dico      # Stock
		

		data1 = [(f" |", f" |", f" |", f"  |")]  # remplacer par %
		t_data1 = pd.DataFrame(data1, columns=[" Sequence reference |", " Quality |", " Quantity |", " % |"])
		t_data = pd.concat([t_data1], axis=1)
		
		for upKey, upValue in qual_dico_sequence.items():
			if totalNumberOfRead != 0:
				data1 = [(f"{upKey} |", f"{key} |", f"{val} |", f" {round(val / totalNumberOfRead * 100, 4)} |") for key, val in upValue.items()] #remplacer par %
				t_data1 = pd.DataFrame(data1, columns=[" Sequence reference |"," Quality |", " Quantity |", " % |"])
				t_data = pd.concat([t_data1], axis=0)
		print(t_data)
		save_file.write(f"{t_data.to_string(index=False)} \n\n")
		
#### Bonus #####
def cigar_analysis():

	print("\n_________________ \n Quantité des cigars :\n")
	save_file.write("\n_________________ \n Quantité des cigars :\n")
	longueur_theorique = 0
	longueur_constate = 0
	
	
	cigar_dico = {}
	
	for sqname in specific_sequence:
		for i in range(len(binary_flags)):      # Loop that iterates from 0 to the amount of reads
			if rname[i] == sqname:              # check if read is mapped to reference sequence
				if re.search(r"(\*)", cigars[i]):   # If cigar is "*" ignore
					continue
				else:
					tmps = re.findall(fr"(\d*)([A-Za-z])", cigars[i])   # If cigar have a digit (or not) and one character
					for tmp in tmps:
						if tmp[0]:
							number = int(tmp[0])
						else:
							number= 1  # default, 1 if no digit
						letter = tmp[1]
						
						if letter in cigar_dico:
							cigar_dico[letter] += number  # Add number found earlier to count the number of read
						else:
							cigar_dico[letter] = number
						longueur_theorique += number
						
	if longueur_theorique != 0:
		data1 = [(f"{cle} |", f"{val} |", f"{round(val / longueur_theorique, 5) * 100} |") for cle, val in cigar_dico.items()]  #
		t_data1 = pd.DataFrame(data1, columns=[" Code |", " Valeur |", " % |"])
		t_data = pd.concat([t_data1], axis=1)
		print(t_data)
		save_file.write(t_data.to_string(index=False))
		save_file.write("\n")
	
####### Start ########

only_step = False                                               # Boolean for launching every function (False) or only selected step (True)
quality_set = False                                             # Boolean for asking the minimal quality (False) or not changing the set quality (True)

if len(sys.argv) == 1:                                          # If no argument (except the SamMapScript.py) :
	print("You did not enter an argument, you can answer the following question or use -h for more details or automate the process. \n") # A little message

for arg in enumerate(sys.argv[1:]):                             # Before all, for every argument except the SamMapScript.py
	if arg[1] in ("-h","--help"):                               # If argument is -h or --help
		show_help()                                             # show help message

#### Open file ####
sam_file_path = find_file()                                     # find_file to find the .sam file

quality_min = 0                                                 # Minimum quality for all other actions (reads with quality < quality_min are ignored as if they never existed)
quality_step = 10                                               # Step for quality slices in distribution
specific_sequence = []                                          # Sequence taken for all other actions (reads with other sequence are ignored as if they never existed)
size_flag = 12                                                  # Minimal size for binary flags

flags_counts = {}                                               # Count of every bits in flag found on the selection

nb_of_save = 0

while True:
	name_save_file = f"Result_n{str(nb_of_save)}_SamMapScript_V1_0_0.txt"
	if not os.path.exists(name_save_file):
		break
	nb_of_save += 1

save_file = open(name_save_file, 'a')                           # Open to append
save_file.write(str(sys.argv))
save_file.write(f"\n\n")

for arg in enumerate(sys.argv[1:]):                             # Check every argument except the first (SamMapScript.py)
	if arg[1].startswith("quality["):                           # if the argument start by "quality[" :
		qual = re.search(r"quality\[(\d+)\]",arg[1])     # Search with the pattern quality\[(\d+)\], to have a group 1 of one or more digit
		quality_min = int(qual.group(1))                        # Put this digit in quality_min
		print(f"The minimum quality will be {quality_min}")     # A little message
		save_file.write(f"The minimum quality will be {quality_min} \n")
		quality_set = True                                      # set the quality_set bool to true to not change the minimum quality after
		
	
	
	if arg[1].startswith("step_qual["):                         # if the argument start by "step_qual[" :
		qualstep = re.search(r"step_qual\[(\d+)\]",arg[1]) # Search with the pattern step_qual\[(\d+)\], to have a group 1 of one or more digit
		quality_step = int(qualstep.group(1))                   # Put this digit in quality_step
		print(f"The analysis of the distribution of qualities will be done in steps of {quality_step}") # A little message
		save_file.write(f"The analysis of the distribution of qualities will be done in steps of {quality_step}\n")
	
	if arg[1].startswith("specific_sequence["):                 # If the argument start by "specific_sequence[" :
		seq = re.search(r"specific_sequence\[(.+)\]",arg[1])    # Search with the pattern specific_sequence\[(.+)\], to have a group 1 of one or more caractére
		if seq:                                                 # If seq is not null
			specific_sequence = seq.group(1).split(',')         # Create a list with every word separate by a comma
			print(f"Only reads aligned to the sequence: {specific_sequence} will be used.") #A little message
			save_file.write(f"Only reads aligned to the sequence: {specific_sequence} will be used.\n")

if not quality_set:                                             # If quality set is equal at False
	ask_quality()                                               # Ask if we want a quality limiter (if only one argument)
	
specific_sequence, flags, rname, mapqs, cigars, totalNumberOfRead = sam_reading(sam_file_path)      # sam_reading function takes the path as parameter and returns the flags, qualities, sequence names, CIGARs and the total quantity of reads
binary_flags = flags_to_binary(flags)                                                               # Convert flags to 12-bit binary

for arg in enumerate(sys.argv[1:]):                             # Check every argument except the first (SamMapScript.py) # Second time because we need the precedent variable
	if arg[1] == "all":
		only_step = True
		flags_counts = dico_count_flag()                            # Launch the function "dico_count_flag" who set the dictionary "flags_counts" for the next action
		statistic_of_mapped_reads(binary_flags)                     # Launch the function "statistic_of_mapped_reads" who take the flags in binary to show the distribution of mapped (or not) reads amongs the sequence
		analysis_flag()                                             # Launch the function "analysis_flag" who present the bits on a human language
		count_flag_number()                                         # Launch the function "count_flag_number" who show in détail the number of every bit
		quality_analysis()                                          # Launch the function "quality_analysis" who show the distribution of quality
		cigar_analysis()                                            # Launch the function "cigar_analysis" who present every number of letter in CIGAR
	
	if arg[1] == "statistic_map":                               # If the argument start by "statistic_map" :
		statistic_of_mapped_reads(binary_flags)                 # Launch the function "statistic_of_mapped_reads" who take the flags in binary to show the distribution of mapped (or not) reads amongs the sequence
		only_step = True                                        # Set the only_step bool to true to not launch other actions that specified
		
	if arg[1] == "count_flag" or arg[1] == "analyse_flag":      # If the argument start by "count_flag" or "analyse_flag" :
		flags_counts = dico_count_flag()                        # Launch the function "dico_count_flag" who set the dictionary "flags_counts" for the next action
		only_step = True                                        # Set the only_step bool to true to not launch other actions that specified
		
	if arg[1] == "analyse_flag":                                # If the argument start by "analyse_flag" :
		analysis_flag()                                         # Launch the function "analysis_flag" who present the bits on a human language
		
	if arg[1] == "count_flag":                                  # If the argument start by "count_flag" :
		count_flag_number()                                     # Launch the function "count_flag_number" who show in détail the number of every bit
		
	if arg[1] == "analyse_quality":                             # If the argument start by "analyse_quality" :
		quality_analysis()                                      # Launch the function "quality_analysis" who show the distribution of quality
		only_step = True                                        # Set the only_step bool to true to not launch other actions that specified
		
	if arg[1] == "analyse_cigar":                               # If the argument start by "analyse_cigar" :
		cigar_analysis()                                        # Launch the function "cigar_analysis" who present every number of letter in CIGAR
		only_step = True                                        # Set the only_step bool to true to not launch other actions that specified
		
#### Number of reads ####
if not only_step:                                               # If the bool only_step is False:
	flags_counts = dico_count_flag()                            # Launch the function "dico_count_flag" who set the dictionary "flags_counts" for the next action
	statistic_of_mapped_reads(binary_flags)                     # Launch the function "statistic_of_mapped_reads" who take the flags in binary to show the distribution of mapped (or not) reads amongs the sequence
	analysis_flag()                                             # Launch the function "analysis_flag" who present the bits on a human language
	#count_flag_number()                                        # Launch the function "count_flag_number" who show in détail the number of every bit
	quality_analysis()                                          # Launch the function "quality_analysis" who show the distribution of quality
	#cigar_analysis()                                           # Launch the function "cigar_analysis" who present every number of letter in CIGAR


###### That's all folks ########
