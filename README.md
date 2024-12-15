# Sam Reading in Bioinformatics
The following scripts were created for the Bioinformatics project in the Systems UE (HAI724I) of the Bioinformatics master's degree at the University of Montpellier.
They answer the following questions in a reference .sam file:
How many reads are mapped? 
How are the reads (and pairs of reads) mapped?  
Where are the reads mapped? Is the alignment uniform along the reference sequence? 
What is the quality of the reads mapped?

# Installation
To run the following scripts, you need the following libraries: CSV, Sys, Os and re (included by default in the python library).
You need also panda (install with PyPI :**pip install pandas**)

# Utilisation

Use SamMapScript.py -h for more help

you can launch the script without parameter, to go on a step by step tutorial
Or you can put parameter to automate :
	SamMapScript.py [file] [option...] will start the analysis of the file with the following options:
    --[Options]--
    quality[int] 	 sets the minimum read quality to analyze to 'int'	(reads with quality < quality_min are ignored as if they never existed)
    step_qual[int] 	 defines the step size for quality intervals as 'int'
    specific_sequence[seq] 	 restricts the analysis to the reference sequence 'seq'; to specify multiple sequences: specific_sequence[seq1,seq2] (be careful with spaces)	(reads with other sequence are ignored as if they never existed)
        To run specific steps only:
    statistic_map	 Provides the quantities and percentages of mapped reads, semi-mapped reads (mapped but with a CIGAR containing other letters than 'M'), and unmapped reads, according to the reference sequences
    analyse_flag	 Provides a human readable analysis of the flags.
    count_flag		 Displays the counts of flag bits across all sequences.
    analyse_quality	 Shows the distribution of quality scores in intervals with a default step size of 10 (or adjustable using step_qual[]).
    analyse_cigar	 Presents the counts of CIGAR letters across all sequences combined. 
    
exemple : /home/user/Téléchargements/SamMapScript.py /home/user/Téléchargements/mapping.sam quality[5] step_qual[15] specific_sequence[reference1,reference2] statistic_map count_flag


