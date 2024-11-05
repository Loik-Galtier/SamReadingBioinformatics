#!/usr/bin/python3
#-*- coding : utf-8 -*-

import os
import re

#aSearching a library wich for browsing folder, mb Tkinter bu lighter if possible

#check the selected file exist and have the .sam
def FindFile ():
    name = input("Path and Name of the file : ?")
    if os.path.exists(name) and re.search(r".sam", name):
        print("There is a file") # continue -> Ask what question to respond (1, 2, 3, 4)
    else:
        print("Don't have a file") # Ask again for a file.
        FindFile()




FindFile()


