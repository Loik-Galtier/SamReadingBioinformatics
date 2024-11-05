#!/usr/bin/python3
#-*- coding : utf-8 -*-

import os
import re
import sys


#check the selected file exist and have the .sam
def FindFile ():
    noFileFound = True
    while noFileFound :
        name = input("Path and Name of the file : ?")
        if os.path.exists(name) and re.search(r".sam", name):
            noFileFound = False
            print("There is a file") # continue -> Ask what question to respond (1, 2, 3, 4)
            with open(name, "a+") as SamFile:
                return SamFile

        else:
            print("Don't have a file") # Ask again for a file.





SamFile = FindFile()

