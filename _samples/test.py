#!/usr/bin/env python

# loading a Tkinter listbox with data and selecting a data line

from Tkinter import *

import math


def get_list(event):
    """
    function to read the listbox selection
    and put the result in a label
    """
    # get line index tuple
    sel = listbox1.curselection()
    # get the text
    seltext = listbox1.get(int(sel[0]))
    label1.config(text=seltext)

# create the sample data file
str1 = """ethyl alcohol
ethanol
ethyl hydroxide
hydroxyethane
methyl hydroxymethane
ethoxy hydride
gin
bourbon
rum
schnaps
"""
fout = open("chem_data.txt", "w")
fout.write(str1)
fout.close()

# read the data file into a list
fin = open("chem_data.txt", "r")
chem_list = fin.readlines()
fin.close()
# strip the trailing newline char
chem_list = [chem.rstrip() for chem in chem_list]

root = Tk()
# create the listbox (note that size is in characters)
listbox1 = Listbox(root, width=50, height=6)
listbox1.grid(row=0, column=0)

# create a vertical scrollbar to the right of the listbox
yscroll = Scrollbar(command=listbox1.yview, orient=VERTICAL)
yscroll.grid(row=0, column=1, sticky=N+S)
listbox1.configure(yscrollcommand=yscroll.set)

# label to display selection
label1 = Label(root, text='Click on an item in the list')
label1.grid(row=1, column=0)

# load the listbox with data
for item in chem_list:
    listbox1.insert(END, item)

# left mouse click on a list item to display selection
listbox1.bind('<ButtonRelease-1>', get_list)

a=sqrt(2.1)*pow(2,9393)
print a

root.mainloop()

# Here's where the program starts
if __name__ == "__main__": 

  # Check if the minimum required external packages are installed
  checkSetup()

  
  # Run it!
  run()

  # That's it - BYE!



