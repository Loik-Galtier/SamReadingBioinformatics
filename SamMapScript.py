#!/usr/bin/python3
#-*- coding : utf-8 -*-
import os
import re
import sys


#check the selected file exist and have the .sam
def findFile ():
    noFileFound = True
    while noFileFound :
        name = input("Path and Name of the file : ?") #Ask the path
        if os.path.exists(name) and re.search(r".sam", name): # if the path exist AND the .sam exist
            noFileFound = False #quit the while
            print("There is a file") # continue -> Ask what question to respond (1, 2, 3, 4)
            return name #return the file open
        else:
            print("Don't have a file") # Ask again for a file.

def chooseParameter():
    print("Choose a option to add for all next operation:")
    print("0 : no special parameter")
    print("1 : Only mapped")
    print("2 : no unmapped")
    print("3 : Max quality of X")
    parameter = input("?")
    return parameter

def flagBits(flag) : #Used same def of the template because is good and if we have to keep the table of 2¹⁶ for analyzing, we gonna do the same part (of convert in bits and set to the 12 line of the table) #look after for change ?

    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB)
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.

    return flagB

def CountReadFlag():
    Occurence = {} #dico
    with open(name, "rt") as SamFile:
        for line in SamFile:
            if not re.search(r"^@", line): #if is not a line commentary
                Field_line = line.split("\t")  # separate for every collone
                if (Field_line[1] in Occurence):
                    Occurence[Field_line[1]] += 1
                else:
                    Occurence[Field_line[1]] = 1
        print(Occurence)

def WhoUnMapped():
    Count = 0
    with open(name, "rt") as SamFile:
        for line in SamFile:
            if not re.search(r"^@", line): #if is not a line commentary
                Field_line = line.split("\t")  # separate for every collone
                flag = flagBits(Field_line[1])  # take the 2nd (0 1) field and pu it in binnary

                if int(flag[-3]) == 1:  # If the 3rd bits (from the last) is = to 1, so the flag is unmapped
                    Count += 1  # add one Unmapped count

        return Count

### Start ###
name = findFile()
parameter = chooseParameter()
print(name)
CountReadFlag()
print(WhoUnMapped())
