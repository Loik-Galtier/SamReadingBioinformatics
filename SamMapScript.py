
import os
import re

#aSearching a library wich for browsing folder, mb Tkinter bu lighter if possible

name = input("Path and Name of the file : ?") #Simplify


#check the selected file exist and have the .sam
if os.path.exists(name) and re.search(r".sam", name):
    print("There is a file")
else:
    print("Don't have a file")
#need to create def for next time



#say if you prefer Aicha that i do on english or french