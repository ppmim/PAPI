################################################################################
#
# PANICtool
#
# autoQueue.py
#
# Last update 19/05/2008
#
################################################################################

"""
   Manages the queues of frames that are waiting to be processed or have been
   processed. Automatic updating of the queues when data appears or disappears
   from disk is also implemented.
"""

# Import external modules

import os
import copy
import config
import configFrame
import fnmatch

import myDialog
import messageLog
import frameLog
import popUp
import fits
import fileList
import misc.dataset

from Tkinter import *


###############################################################
# Data base ( we use the same instance misc.dataset.filesDB)
import misc.dataset

#######################################################
  
################################################################################

# Define local constants

autocheckwaitingtime = 5000	# Waiting time between automatic queue updates

################################################################################

class window(Frame):

  """
     Provides a window to manage the unprocessed and processed queues, as
     well as methods to automatically update these queues. Class inherits
     from Tkinter's 'Frame' class.
  """

  def __init__(self, parent):

    # Inialize as if this were an instance of the Frame class.
    Frame.__init__(self, parent)

    # Make this frame flexible in width
    self.columnconfigure(0, weight=1)

    # Define a Tkinter string variable containing the standard filter pattern
    # for automatic addition of files to the 'waiting list'
    self.filterpattern = StringVar(master=self)
    self.filterpattern.set('*.fits')

    # Define Tkinter integer variables containing the number of unprocessed
    # files, and the number of reduced files.
    self.queuesize = IntVar(master=self)
    self.queuesize.set(0)

    self.reducedsize = IntVar(master=self)
    self.queuesize.set(0)

    # Define a Tkinter integer variable that defines if aumatic queue
    # updating should occur or not.
    self.autocheck = IntVar(master=self)
    self.autocheck.set(0)

    # Define the two lists containing the filenames of unprocessed and reduced
    # files.
    self.newfiles     = []
    self.reducedfiles = []

    # A local variable that will contain the ID of the waiting loop.
    self.checkID = None

    # Now construct the frame
    self.makeFrame()


  def makeFrame(self):

    "Construct the frame from scratch"

    # Put all data in a root frame that is embedded in the 'self' Frame object.
    # In this way, it is possible to destroy and rebuid the frame object,
    # without destroying the 'self' frame object.
    
    self.root = Frame(self)
    
    # Give the root frame a minimum size, and make it flexible in the
    # middle column
    self.root.columnconfigure(0, weight=0, minsize=200)
    self.root.columnconfigure(1, weight=1)
    self.root.columnconfigure(2, weight=1)


    # Row counter
    i = 1

    # Label for first row
    Label(self.root, anchor=W, text='Number of files in queue', width=25).grid(row=i,
                                                     column=0, sticky=E+W)

    # Text entry for first row, containing the number of unprocessed files
    nqueue = Entry(self.root, textvariable=self.queuesize, state=DISABLED)

    # Try to make the text color of this item black.
    # For compatibility purposes with Python versions < 2.3.
    try:
      nqueue.config(disabledforeground='black')
    except TclError: pass

    # And display this label   
    nqueue.grid(row=i, column=1, sticky=E+W)

    # Button to modify the selection of unprocessed files
    queueButton = Button(self.root, text='Edit queue', width=25)
    queueButton.config(command = self.editNewFiles)
    queueButton.grid(row=i, column=2, sticky=E+W)

    # And attach popup-help to this button.
    queuePopup = popUp.popUp(self.root, queueButton, title='Edit files in queue',
                          text="""Select or deselect files that are pending to be
			          processed. Selecting files that are already processed
				  will cause them to be reprocessed.
			       """)


    # Next row
    i = i + 1

    # Label for second row
    Label(self.root, anchor=W, text='Number of processed files', width=25).grid(row=i,
                                                     column=0, sticky=E+W)

    # Text entry for second row, containing the number of processed files
    nproc = Entry(self.root, textvariable=self.reducedsize, state=DISABLED)

    # Try to make the text color of this item black.
    # For compatibility purposes with Python versions < 2.3.
    try:
      nproc.config(disabledforeground='black')
    except TclError: pass

    # And display this label   
    nproc.grid(row=i, column=1, sticky=E+W)

    # Button to modify the selection of reduced files
    listButton = Button(self.root, text='Edit list of processed files', width=25)
    listButton.config(command = self.editReducedFiles)
    listButton.grid(row=i, column=2, sticky=E+W)

    # And attach popup-help to this button.
    listPopup = popUp.popUp(self.root, listButton, title='Edit list of processed files',
                          text="""View and change the list of processed files. Deselecting
			          processed files will NOT cause them to be reprocessed.
			       """)


    # Next row
    i = i + 1


    filterEntry = configFrame.window(self.root, optionlist=['filename_filter'], embed=True, width=100)
    filterEntry.grid(row=i, column=0, columnspan=2, sticky=E+W)

    # Label for third row
    #Label(self.root, width=25, anchor=W, text='Filename filter').grid(row=i,
    #                                                 column=0, sticky=E+W)

    # Text entry for third row, containing the filter pattern used when
    # automatically selecting items for the unprocessed list.
    #filterEntry = Entry(self.root, textvariable=self.filterpattern)
    #filterEntry.grid(row=i, column=1, sticky=E+W)

    # Attach popup-help to this entry.
    filterPopup = popUp.popUp(self.root, filterEntry, title='Filename filter',
                          text="""When AutoCheck is running, files that match this filter
			          will automatically be appended to the queue of files
				  waiting to be processed. 
			       """)


    # Checkbutton to enable/disable automatic checking for the arrival of
    # unprocessed files that match the pattern in the entry above
    autoButton = Checkbutton(self.root, width=25, text='Autocheck for new files',
                               variable=self.autocheck)
    autoButton.grid(row=i, column=2, sticky=E+W)

    # Associate with this checkbox the method that starts automatich checking
    autoButton.config(command = self.autoCheckFiles)

    # Attach popup-help to this entry.
    autoPopup = popUp.popUp(self.root, autoButton, title='Filename filter',
                          text="""Turning this option on will start periodic checking
			          of the arrival of new files in the input directory. If
				  a newly created file matches the filename filter, it will
				  be appended to the queue of files waiting to be processed.
				  Checking this box will not automatically start the data
				  processing. Use the 'START processing' button for that.
			       """)

    # Finally, display the entire frame (that is, the 'root' frame)
    self.root.grid(sticky=E+W)



  def updateQueueSize(self):

    "Update the lengths of the unprocessed and processed queues on screen"

    messageLog.put('Updating queue sizes', 9)
    self.queuesize.set(len(self.newfiles))
    self.reducedsize.set(len(self.reducedfiles))



  def autoCheckFiles(self):

    """
       Invoked when the autocheck checkbox is toggled. Depending upon whether
       it was selected or deselected, will start or stop the automatic
       checking for new files for the unprocessed queue.
    """

    # Is the autocheck button (still) selected?
    if self.autocheck.get():
    
      # Yes, it was...

      # Make sure that inpath and outpath are different. If they are the
      # identical, a newly reduced file in the 'output' directory would appear
      # to be a new, unreduced file in the 'input' directory, triggering a bad
      # runaway.

      inpath       =  os.path.normpath(config.options['inpath'].value)
      outpath      =  os.path.normpath(config.options['outpath'].value)

      if inpath == outpath:
	messageLog.put_error('Putting "Output directory" equal to "Input directory" will cause')
	messageLog.put_error('runaway problems in combination with "Autocheck for new files".')
	messageLog.put_error('Autochecking stopped!')
        # In this case, force an untoggle of the button and return
        self.autocheck.set(0)
	return

      # If this is the first time this routine is called, prepare the
      # directory list.
      if (self.checkID is None):
      
        # Store the current contents of the input directory for future reference
        self.dirlist =  [os.path.join(inpath, file) for file in os.listdir(inpath)]
        messageLog.put('Autochecking started', 5)

      # Check for new files in the input directory
      self.findNewFiles()
      
      # Call this routine again later
      self.checkID = self.after(autocheckwaitingtime, self.autoCheckFiles)


    else:

      # Well, if the autocheck button is not selected (anymore)

      # Try to cancel the next call to this routine (if there is one scheduled)
      try: self.after_cancel(self.checkID)
      except: pass

      # Reset the ID of the waiting loop
      self.checkID = None

      messageLog.put('Autochecking stopped', 5)



  def findNewFiles(self):

    """
       Find files that are not yet listed in the queue of unprocessed
       files, and match these with the filename filter pattern
    """

    
    # Retrieve the current value of the filename pattern from the widget
    pattern = config.options['filename_filter'].value

    # Retrieve the current value of the data input directory
    inpath  = config.options['inpath'].value

    # Read the directory contents
    contents = [os.path.join(inpath, file) for file in os.listdir(inpath)]

    # Check the obtained list of files agains the existing directory list
    # Remove files that already existed in the directory list

    # ORDER of if-statements is important!

    iterlist = copy.copy(self.dirlist)
    for file in iterlist:

      if file not in contents:
        # Hmm... a strange situation. Apparently a file listed in self.dirlist
	# DISappeared from the directory. Adjust the lists accordingly
        messageLog.put_warning('File %s disappeared from directory - updating lists' % file, 5)
        self.dirlist.remove(file)

        ## new---
        #fc=fits.FitsFile(file)
        #fc.recognize()
        #fits.detected_files.remove( file )
        # Show into the frameLog
        #frameLog.put_warning (os.path.basename(file) + " --- " + "Removed !" )
        ## wen---

        # Do NOT swap the following two statements!

        # Is this file already in the list of processed files?
	if file in self.reducedfiles: self.reducedfiles.remove(file)

        # Is this file already in the list of unprocessed files?
	if file in self.newfiles:     self.newfiles.remove(file)


      if file     in contents:
        # Normal situation, all files in self.dirlist are also in the current
	# directory listing. Remove these files one-by-one from the list, so that
	# the remaining files are those that are the new files in the input
	# directory (this time).
        contents.remove(file)

    # Now loop over the remaining files
    for file in contents:

      # And append these to 'self.dirlist' for future reference
      self.dirlist.append(file)

      # Only look at the filename (disregard from directory path)
      basename = os.path.basename(file)

      # Make sure that (1) the 'file' is not a directory entry, (2) that it
      # matches the filename filter pattern, and (3) that it not yet in
      # the unprocessed or processed list.
      if   ( (not os.path.isdir(file))
         and (fnmatch.fnmatch(basename, pattern))
         and (file not in self.reducedfiles)
         and (file not in self.newfiles) ):

        messageLog.put('Appending file to the queue', 9)

        ## new---

        fc=fits.FitsFile(file)
        #fc.recognize()
        #fits.detected_files.append( file )
        # Append to the frameLog
        fileList.put (fc.filename , fc.type )
        misc.dataset.filesDB.insert(file)
        misc.dataset.filesDB.ListDataSet()
        
        ## wen---
        
        # Only then, append the file to the list of files to process.
        self.newfiles.append(file)


    # Finally, update queue sizes on screen
    self.updateQueueSize()



  def editNewFiles(self):

    "Open dialog to manually (de)select files in the unprocessed list"

    newlist  = myDialog.LoadFITSFiles(parent=self,
                    initialdir=config.options['inpath'].value,
		    headlist=config.options['fitsheaders'].value,
		    pattern=config.options['filename_filter'].value,
                    default=self.newfiles)

    if newlist is not None:

      self.newfiles = newlist

      # If the file previously appeared in the other list, remove it from
      # that list
      for file in newlist:
        if file in self.reducedfiles: self.reducedfiles.remove(file)

      # Update queue sizes on screen
      self.updateQueueSize()



  def editReducedFiles(self):

    "Open dialog to manually (de)select files in the processed list"

    newlist  = myDialog.LoadFITSFiles(parent=self,
                    initialdir=config.options['inpath'].value,
		    headlist=config.options['fitsheaders'].value,
		    pattern=config.options['filename_filter'].value,
                    default=self.reducedfiles)

    if newlist is not None:

      self.reducedfiles = newlist

      # If the file previously appeared in the other list, remove it from
      # that list
      for file in newlist:
        if file in self.newfiles: self.newfiles.remove(file)

      # Update queue sizes on screen
      self.updateQueueSize()
