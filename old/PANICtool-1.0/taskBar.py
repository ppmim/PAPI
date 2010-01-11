################################################################################
#
# FIEStool
#
# taskBar.py
#
# Last update 23/05/2005
#
################################################################################

"""
   Provide task- and status-bar for pipeFrame or calibFrame.
"""


# Import necessary modules

import popUp

from Tkinter import *

################################################################################



################################################################################

class taskBar(Frame):

  "Create a task bar with status frame for a given 'task'"
  
  def __init__(self, parent, task, checkbuttonvar=None, bgcolor=None):
  
    "Setup the initial bar"

    # Save task as global variable for this class
    self.task = task


    # Inialize as if this were an instance of the Frame class
    Frame.__init__(self, parent, bg=bgcolor)

    # Save the default background color for future reference
    self.defaultbgcolor = self.cget("bg")

    # And make 'self' frame flexible in width
    self.columnconfigure(0, weight=0, minsize=200)
    self.columnconfigure(1, weight=1)

    # And construct the bar
    self.makeBar(checkbuttonvar)


  def makeBar(self, checkbuttonvar=None):

    """
      Create the bar from scratch. Setting a 'checkbuttonvar' Tkinter variable
      will provide an on/off button and hold its status.
    """

    # If a Tkinter integer variable is provided, create a checkbox and attach
    # the variable to the checkbox status. The attached text is not clickable.
    if (checkbuttonvar) :
      self.button = Checkbutton(self,
			highlightthickness = 0,
			text     = self.task.buttonText,
			variable = checkbuttonvar,
			anchor   = W,
			bg       = self.defaultbgcolor,
			command  = self.highlight_off)
    else :
      # If no variable is given, put the buttontext of the task in a button
      self.button = Button(self,
			highlightthickness = 0,
			bg       = self.defaultbgcolor,
		        text     = self.task.buttonText)


    # Put the button in the frame and attach a popup help text to it
    self.button.grid(row=0, column=0, sticky=E+W)
    popUp.popUp(self, self.button, self.task.buttonText, self.task.__doc__)

    # Create a default status label
    self.status = Label(self, anchor=W, text='Not done...', bg=self.defaultbgcolor)
    self.status.grid(row=0, column=1, columnspan=2, sticky=E+W)
    # And attach a help popup
    popUp.popUp(self, self.status, title="Status",
               text="This field shows the latest status message related to this task")

    # Save the button and the label in self.elements for future reference
    self.elements = [self.button, self.status]

    # There is no '.grid()' call here. The bar is to be gridded by the caller!


  def highlight(self, color='#339933'):    # '#339933' is a nice green colour

    "Highlight the task bar by changing the background"

    self.config(bg = color)
    for element in self.elements:
      element.config(bg = color)


  def highlight_off(self):

    "Reset the background color to the default value"

    self.config(bg = self.defaultbgcolor)
    for element in self.elements:
      element.config(bg = self.defaultbgcolor)


  def set_label(self, text=''):

    "Change the value of the label to 'text'"
    
    self.status.config(text = text)


  def disable(self):

    "Disable the task bar"

    self.button.config(state=DISABLED)
  

  def enable(self):

    "Enable the task bar"

    self.button.config(state=NORMAL)
