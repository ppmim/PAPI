################################################################################
#
# FIEStool
#
# configFrame.py
#
# Last update 05/10/2005
#
################################################################################

"""
   Provide a window to interactively edit (a given list of) settings defined by
   the 'config' module. This window can either be a top-level window, or
   embedded in another window. 
"""

# Import external modules

import config
import Dialog
import myDialog
import messageLog
import copy
import os
import popUp

from Tkinter import *

################################################################################

class AutoScrollbar(Scrollbar):
    # a scrollbar that hides itself if it's not needed.  only
    # works if you use the grid geometry manager.
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            # grid_remove is currently missing from Tkinter!
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        Scrollbar.set(self, lo, hi)
    def pack(self, **kw):
        raise TclError, "cannot use pack with this widget"
    def place(self, **kw):
        raise TclError, "cannot use place with this widget"

################################################################################

class window(Frame):

  """
     Provide a window to edit settings defined in the 'config' module
  """

  def __init__(self, parent, title='', optionlist=[], embed=False, width=200):

    # 'parent' is the Tkinkter object under which this window is
    # hierarchially placed. 'Title' defines the frame header and 'optionlist'
    # the list of option names that are to be listed in this frame. Setting
    # 'embed' to 1 will make this frame as an embeddable object (an instance
    # of Tkinter.Frame that should be gridded into a Tkinter container object).

    # embed and title are global variables in this class
  
    self.title = title
    self.embed = embed
    self.minwidth = width

    # Check if frame is to be embedded

    if not self.embed :

      # If not, create the top level window (hidden from view)
      self.top = Toplevel(parent)
      self.top.withdraw()
      
      # Give the top-level frame a title
      self.top.title(self.title)
      
      # ...and reset the title parameter, so that it will no be displayed
      # once again...
      self.title=''
      
      # Make this top-level frame flexible in width
      self.top.columnconfigure(0, weight=1)
      self.top.rowconfigure(0, weight=1)

      # Clicking 'close window' on this frame calls the self.hide method
      self.top.protocol("WM_DELETE_WINDOW", self.hide)

      # Create a scrollbar that will disappear if not necessary
      vscrollbar = AutoScrollbar(self.top)
      vscrollbar.grid(row=0, column=1, sticky=N+S)

      # Create a Canvas object that can contain the whole frame
      self.canvas = Canvas(self.top, yscrollcommand=vscrollbar.set)
      self.canvas.columnconfigure(0, weight=1)
      self.canvas.rowconfigure(0, weight=1)

      # Attach the method to position the canvas to the scrollbar
      vscrollbar.config(command=self.canvas.yview)


    else:

      # This frame is to be embedded in 'parent'
      self.canvas = parent

    # Inialize as if this were an instance of the Frame class
    Frame.__init__(self, self.canvas)

    # And make also the 'self' frame flexible in width
    self.columnconfigure(0, weight=1)

    if not self.embed:
      self.rowconfigure(0, weight=1)

    # Optionlist contains the names of the options to be displayed in this
    # frame. If it is not defined, obtain _all_ options from 'config'. Note
    # that in that case the list will be unsorted.
    if optionlist: self.optionlist = optionlist
    else: self.optionlist = config.options.keys()

    # Make a (deep) copy of the optionlist (for when 'reset' is called)
    self.oldoptions = {}
    for key in self.optionlist:
      self.oldoptions[key] = copy.deepcopy(config.options[key])

    # Build the actual frame
    self.makeFrame()

    # And update the frame
    self.update()

    if not self.embed:
      # If the frame is not embedded, display it
      self.canvas.create_window(0, 0, anchor=NW, window=self)

    else:
      # Otherwise, leave the displaying (gridding) for later, and
      # start an automatic periodic update (because settings may change
      # while the frame is visible)
      self.autoupdate()


  def makeFrame(self):

    "Construct the frame from scratch"

    # Put all data in a root frame that is embedded in the 'self' Frame object.
    # In this way, it is possible to destroy and rebuid the frame object,
    # without destroying the 'self' frame object.
    self.root = Frame(self)

    # All columns of the root frame should be flexible in width
    self.root.columnconfigure(0, weight=1, minsize=200)
    self.root.columnconfigure(1, weight=1, minsize=self.minwidth)
    self.root.columnconfigure(2, weight=1, minsize=self.minwidth)

    # Initialize (reset) the list of button objects (for future reference)
    self.optionbuttonlist = {}

    # Row counter
    i = 1

    # If a title is given, display a label spanning the entire width of
    # the frame
    if self.title:
      label = Label(self.root, text=self.title)
      label.grid(columnspan=3, row=i, sticky=E+W)
      # Next row
      i = i + 1


    # Loop through the list of options
    for key in self.optionlist:

      # Get the value of the named option
      try:
        option = config.options[key]
      except KeyError:
        # Diaplay an error in the log if the option does not exist
        messageLog.put_error("Cannot find option '%s' in current configuration!" % key)
	# Skip the creation of an optionButton and continue with the next option
	continue

      # Create an optionButton (button and value field) for this option
      button = optionButton(self.root, option, width=self.minwidth)
      button.grid(column=0, columnspan=3, row=i, sticky=E+W)

      # Save this Tkinter object in 'optionbuttonlist' for future reference
      self.optionbuttonlist[key] = button

      # Newt row
      i = i + 1

    # Only if this is a top-level window, add 'Save', 'Load' and 'Reset'
    # buttons, and an 'Close window' button to hide the window.
    if not self.embed:

      button = Button(self.root, text='Save settings')
      # Attach the 'saveConfig' method to this button
      button.config(command=self.saveConfig)
      button.grid(column=0, row=i, sticky=E+W)
      # Add pop-up help to this button
      popUp.popUp(self.root, button, title='Save settings',
                  text="""Save the settings above (and only these) to file.
		       """)

      button = Button(self.root, text='Load settings')
      # Attach the 'loadConfig' method to this button
      button.config(command = self.loadConfig)
      button.grid(column=1, row=i, sticky=E+W)
      # Add pop-up help to this button
      popUp.popUp(self.root, button, title='Load settings',
                  text="""Load settings above (and only these) from file.
		       """)

      button = Button(self.root, text='Reset')
      # Attach the 'resetConfig' method to this button
      button.config(command = self.resetConfig)
      button.grid(column=2, row=i, sticky=E+W)
      # Add pop-up help to this button
      popUp.popUp(self.root, button, title='Reset',
                  text="""Revert to the most recently saved or loaded settings.
		       """)

      # Next row
      i = i + 1

      # Add the 'Close window' button, and attach the hide method
      button = Button(self.root, text='Close window', command=self.hide)
      button.grid(row=i, columnspan=3, sticky=E+W)



    # Display the contents of the roor window
    self.root.grid(sticky=N+S+E+W)


  def show(self):

    """
       If embedded, display this window. If top-level, then deiconify and set
       focus and grab to this window
    """

    if self.embed:
      self.grid(sticky=E+W)
    else:
      self.top.deiconify()
      self.top.focus_set()
      self.top.grab_set()
      self.top.wait_visibility()
    
      self.canvas.grid(row=0, column=0, sticky=N+S+E+W)

      # Determine width and heigth of the contained frame
      canvaschild  = self.canvas.children.values()[0]
      neededwidth  = canvaschild.winfo_width()
      neededheight = canvaschild.winfo_height()

      maxheight = self.canvas.winfo_screenheight()
      maxwidth  = self.canvas.winfo_screenwidth()

      neededheight = min(maxheight, neededheight)
      neededwidth  = min(maxwidth,  neededwidth)

      # Set the canvas size and enable scrolling
      self.canvas.config(width=neededwidth, height=neededheight)
      self.canvas.config(scrollregion=self.canvas.bbox("all"))

    # After displaying, update the contents
    self.update()


  def hide(self):
  
    "Hide the window from view"

    if self.embed:
      # Probably never called (requires makeFrame to reconstruct the contents)
      self.grid_forget()
    else:
      self.top.grab_release()
      self.top.withdraw()


  def update(self):
  
    "Update the contents in the window"

    # Loop over the available optionButtons, and call their update method
    for name in self.optionbuttonlist.iterkeys(): 
      self.optionbuttonlist[name].option = config.options[name]
      self.optionbuttonlist[name].update()


  def autoupdate(self):

    "Periodically update the window, because settings may change while on display"

    # Update
    self.update()

    # And reschedule itself in one second from now
    self.after(1000, self.autoupdate)


  def saveConfig(self, outfile=None):

    "Save the displayed configuration data to disk"

    if outfile is None:

      # Ask for filename if not yet given
      outfile = myDialog.SaveFile(parent=self, initialdir=config.options['config_dir'].value,
                                  filetypes=[('Config files', '*.cfg'), ('All files', '*')])


    if outfile:

      # Call config.save method for the listed options (only!)
      # jmiguel-19-03-2008: Modified to vave all options, not only the listed in main frame
      config.save(outfile) #, optionlist=self.optionlist)

      # Store current options in 'oldoptions', in case restore is issued later
      for key in self.optionlist:
	self.oldoptions[key] = copy.deepcopy(config.options[key])


  def loadConfig(self, infile=None):

    "Load the displayed configuration options from disk"

    if infile is None:

      # Ask for filename if not yet given
      infile = myDialog.LoadFile(parent=self, initialdir=config.options['config_dir'].value,
                                 filetypes=[('Config files', '*.cfg'), ('All files', '*')])


    if infile:

      # Read listed options (only!) from file into config object
      infile = os.path.abspath(infile)
      config.load(infile, optionlist=self.optionlist)

      # Store current options in 'oldoptions', in case restore is issued later
      for key in self.optionlist:
	self.oldoptions[key] = copy.deepcopy(config.options[key])

      # Update frame with new values
      self.update()


  def resetConfig(self):

    "Revert option values to that of last load or save"

    # Check with user through small pop-up window 
    answer = Dialog.Dialog(self, title='', bitmap='',
           strings=('Yes', 'No'), text='Revert to last load/save?', default=1)

    if answer.num == 0:

      # Restore the config options to the old values in oldoptions
      for key in self.optionlist:
        config.options[key] = copy.deepcopy(self.oldoptions[key])

      messageLog.put('Restored configuration from last load/save', 3)

      # Update frame with 'new' values
      self.update()

################################################################################

class optionButton(Frame):

  """
     Provides a Tkinter.Frame instance containing a Button and a Label showing
     the current value of one single option. Clicking the Button will edit the
     value.
  """

  # Constant for this class :

  # A dictionary with the myDialog methods to call when a
  # parameter of certain type is to be edited
  functionlist = {'directory'		:	myDialog.GetDirectory,
    		  'file'		:	myDialog.LoadFile,
    		  'fitsfile'		:	myDialog.LoadFITSFile,
    		  'outfile'		:	myDialog.SaveFile,
    		  'outfitsfile'		:	myDialog.SaveFITSFile,
		  'filelist'		:	myDialog.LoadFiles,
		  'fitsfilelist'	:	myDialog.LoadFITSFiles,
		  'integer'		:	myDialog.GetInteger,
		  'other'		:	myDialog.GetOther,
		 }
		 # NB: Selection with 'save file' will not ask for permission
		 #     when overwriting a file... Something to improve?


  def __init__(self, parent, option, width=200):

    # Option and parent are global variables in this class
    self.option   = option
    self.minwidth = width
    self.parent   = parent

    # Inialize as if this were an instance of the Frame class.
    Frame.__init__(self, parent)

    # And make also the 'self' frame flexible in width
    self.columnconfigure(0, weight=0, minsize=200)
    self.columnconfigure(1, weight=1, minsize=self.minwidth)

    # Now construct te button
    self.makeButton()

    # Attach the edit method corresponding to this option's type to self.editfunc
    self.editfunc = self.functionlist[self.option.type]

    # And update this button
    self.update()


  def makeButton(self):

    "Construct an optionButton from scratch"

    # Add a button in the leftmost column, and attach the edit method to it
    self.button = Button(self, text=self.option.name, command=self.edit)
    # Display the button
    self.button.grid(row=0, column=0, sticky=E+W)

    # Add pop-up help to this button (text comes from option.desc)
    buttonPopup = popUp.popUp(self, self.button, title=self.option.name, text=self.option.desc)

    # Make a (disabled) Entry field
    self.label = Entry(self, state=DISABLED, width=0)

    # Try to make the text color of this item black.
    # For compatibility purposes with Python versions < 2.3.
    try:
      self.label.config(disabledforeground='black')
    except TclError: pass

    # Display the label in the right hand column
    self.label.grid(row=0, column=1, sticky=E+W)


  def edit(self):

    """
       Prompt the user for the new value of this optionButton and update the
       value in memory
    """

    # If the option is an input directory, its initial value should be the
    # current value. Other options do not require this
    if self.option.indir:
      if self.option.type in ('fitsfile', 'fitsfilelist', 'outfitsfile'):
        pattern=config.options['filename_filter'].value
      else:
        pattern=None
      newval = self.editfunc(parent=self, initialdir=self.option.indir.value, default=self.option.value, pattern=pattern)
    else:
      newval = self.editfunc(parent=self, initialvalue=self.option.value)

    # If a new value was obtained, adjust the corresponding option
    if newval is not None:
      self.option.value = newval

    # And update the value of this optionButton
    self.update()


  def update(self):

    "Update the value displayed in the label of this optionButton"

    # Change the label temporarily to 'normal' state, to allow
    # this routine to modify its contents
    self.label.config(state=NORMAL)

    # Remove the current text
    self.label.delete(0, END)

    # Write the current value to the label
    self.label.insert(END, self.option.nicevalue())

    # And restore the state to 'disabled', to prohibit the user from changing
    # this field directly
    self.label.config(state=DISABLED)
