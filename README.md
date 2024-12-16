# Sam Reading in Bioinformatics
*By Loik Galtier*
---
## Intro
The following script were created for the Bioinformatics project of the EU Systems (HAI724I) of the Bioinformatics master of the University of Montpellier.
They produce data, to answer the following questions in a reference .sam file but can be use for more:
How many reads are mapped?
How are the reads (and pairs of reads) mapped?
Where are the reads mapped? Is the alignment uniform along the reference sequence?
What is the quality of the mapped reads?
It will present the results in table form !

## Installation
To run the following scripts, you need python 3.12. Maybe you need PIP or conda to install library.
You can use the following command to get the project :
```git clone https://github.com/Loik-Galtier/SamReadingBioinformatics.git```
You need the following libraries: CSV, Sys, Os and re (included by default in the python library).
You also need panda, install it with PyPI: 
```pip install pandas```

## Usage

```Use SamMapScript.py -h``` for more help

you can run the script without parameters, to follow a step-by-step tutorial
Or you can give a parameter to the automation:
```SamMapScript.py [file] [option...]```
will start the analysis of the file with the following options:

--[Options]--
    **quality[int]** : sets the minimum read quality to be analyzed to 'int' (reads with quality < quality_min are ignored as if they never existed)
    **step_qual[int]** : sets the step size for quality intervals to 'int'
    **specific_sequence[seq]** : limits the analysis to the reference sequence 'seq'; to specify multiple sequences: specific_sequence[seq1,seq2] (watch out for spaces) (reads with another sequence are ignored as if they never existed) had never existed)

To run only specific steps:
    **all** : all the next action.
    **statistic_map** : Provides the amounts and percentages of mapped reads, semi-mapped reads (mapped but with a CIGAR containing letters other than "M"), and unmapped reads, based on the reference sequences
    **analyze_flag** : Provides a human-readable analysis of the flags.
    **count_flag** : Displays the number of flag bits across all sequences.
    **analyse_quality** : Displays the distribution of quality scores in intervals with a default step size of 10 (or adjustable using step_qual[]).
    **analyze_cigar** : Displays the number of CIGAR letters across all sequences combined.

example: ```SamReadingBioinformatics/SamMapScript.py /home/user/Downloads/mapping.sam quality[5] step_qual[15] specific_sequence[reference1,reference2] statistic_map count_flag```


## Validation
You can run with the file test.sam in the github and compare your result with the **Result_testsam_SamMapScript_V1_0_0.txt**

---
## Licence
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by) the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. "
You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
