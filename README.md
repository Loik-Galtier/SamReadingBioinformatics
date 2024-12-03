# Sam Reading in Bioinformatics
The following scripts were created for the Bioinformatics project in the Systems UE (HAI724I) of the Bioinformatics master's degree at the University of Montpellier.
They answer the following questions in a reference .sam file:
How many reads are mapped? 
How are the reads (and pairs of reads) mapped?  
Where are the reads mapped? Is the alignment uniform along the reference sequence? 
What is the quality of the reads mapped?

# Installation
To run the following scripts, you need the following libraries: CSV, Sys, Os and re (included by default in the python library).


# Utilisation

Use SamMapScript.py -h for more help

you can launch the script without parameter, to go on a step by step tutorial
Or you can put parameter to automate :
SamMapScript.py [file...] [Quality Min] [Alignement need] [Specific operation only]

