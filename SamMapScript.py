
import os
import re


#add Tkinter in the futur for all graphics and add drag and drop
name = input("Name of the file : ?") #file name (can work with alone)
path = input("File directory path : ?") #direcory WIHOUT the file name (cant work without name)


#check the selected path + file exist and have the .sam  #wich tkinter put, the name can be change without changing the path. #or not.
if os.path.exists(path + name) and re.search(r".sam", name):
    print("There is a file")
else:
    print("Don't have a file")
#need to create def for next time


#say if you prefer Aicha that i do on english or french