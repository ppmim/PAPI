################################################################################
#
# FIEStool
#
# popUp.py
#
# Last update 23/05/2005
#
################################################################################

"""
   Provide simple popup windows holding help texts and attach it to
   existing Tkinter objects on screen
"""

# Import necessary modules

import messageLog
import string

from Tkinter import *

################################################################################

# Set local variables :
# Default frame width, height and font

framewidth=40
frameheight=5
framefont='courier 10'

################################################################################

class popUp:

  def __init__(self, parent, object, title='No title', text='No text defined'):

    "Setup the popup window with given 'title' and 'text', and bind it to 'object'"

    # The popup itself is a toplevel window
    self.top = Toplevel(parent, borderwidth=3, bg='blue')
    # and should not have any window manager decorations
    self.top.overrideredirect(1)

    # Create a text container...
    self.popwindow = Text(self.top, font=framefont, width=framewidth, height=frameheight, wrap='word')
    self.popwindow.config(state=NORMAL)
    # ...and fill it with a title line
    self.popwindow.insert(END, '%s\n\n' % title.center(framewidth))

    # Proprocess the text line by removing all excessive whitespace of 'text'
    if type(text) != type(""): text = ""
    newtext = string.join(string.split(text))

    # Insert text in content widget
    self.popwindow.insert(END, newtext)
    self.popwindow.config(state=DISABLED)
    self.popwindow.grid(sticky=E+W)

    # Bind popup methods to the object
    object.bind('<Button-3>', self.show)
    object.bind('<ButtonRelease-3>', self.hide)

    # Hide the popup from view (initial state)
    self.top.withdraw()


  def show(self, event):

    "Display the popup on screen, at the current mouse position"

    self.top.geometry('+%i+%i' % (event.x_root, event.y_root))
    self.top.deiconify()

    # Automatically adjust vertical size if not all text is visible
    while self.popwindow.yview() != (0.0, 1.0):

      currentheight = int(self.popwindow.cget('height'))
      self.popwindow.config(height=currentheight+1)
      # Redisplay the window in order to update the result of yview()
      self.top.update_idletasks()

      # There is an upper limit to the window height. Mainly to catch a runaway.      
      if currentheight > 100:
        messageLog.put('Popup window size exceeds acceptable range!')
	break


  def hide(self, event):

    "Hide the popup from view"

    self.top.withdraw()
