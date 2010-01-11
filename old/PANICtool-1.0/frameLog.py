################################################################################
#
# PANICtool
#
# frameLog.py
#
# Last update 07/04/2008
#
################################################################################

"""
   Provides facility for detected frames , showing the frame name and type, as
   well as the option to display a selected frame.
"""

# import external modules

import config
import myDialog
import Queue
import popUp

from Dialog import Dialog
from time import strftime
from Tkinter import *

################################################################################

# Define local constants

# 'log' is an instance of Queue, in which messages temporarily can be stored
# until they are retrieved and displayed by the logging window.

# The major advantage with this implementation is the fact that 'log' is an 
# instance, and not a variable. Therefore, in other modules that import
# frameLog, frameLog.log will refer to exactly the same instance (and not
# a local variable)

log      = Queue.Queue()

autocheckwaitingtime =  1000	# Waiting time between automatic updates
autologsavewait      =  5000	# Waiting time between log saves
autologsavefile      = "lastsession.log"

# Three helper routines to put messages, errors or warnings in the log queue

def put(message, loglevel=0, errorflag=0):
  log.put((message, loglevel, errorflag))

def put_error(message, loglevel=0, errorflag=1):
  log.put((message, loglevel, errorflag))

def put_warning(message, loglevel=0, errorflag=2):
  log.put((message, loglevel, errorflag))


################################################################################

class window(Frame):

  "Provide a window to display a (filtered) list of messages"

  def __init__(self, parent):

    # Inialize as if this were an instance of the Frame class
    Frame.__init__(self, parent)
    
    # Define the initial loglevel (used for filtering) to be 5. The value
    # is stored in a Tkinter integer variable, so it can easily be edited or
    # implemented by other widgets

    self.loglevel = IntVar(self)
    self.loglevel.set(5)

    # Make this frame flexible in width and height
    # (The messageLog eidget is the only widget flexible in height)
    self.columnconfigure(0, weight=1)
    self.rowconfigure(0, weight=1)

    # Build the frame
    self.makeFrame()

    # Local variables that contain the ID of the waiting loop
    self.updateID = None
    self.autologsaveID = None

    # And start automatic updating
    self.autoupdate()
    self.autologsave()


  def makeFrame(self):

    "Build the log window from scratch"

    # Put all data in a root frame that is embedded in the 'self' Frame object.
    # In this way, it is possible to destroy and rebuid the frame object,
    # without destroying the 'self' frame object.

    self.root = Frame(self)

    # All columns of the root frame should be flexible in width
    self.root.columnconfigure(0, weight=0, minsize=120)
    self.root.columnconfigure(1, weight=1)
    self.root.columnconfigure(2, weight=1)

    # The root frame should also be flexible in height
    self.root.rowconfigure(1, weight=1)

    # Create a text box that will contain the log data
    self.box = Text(self.root)
    
    # Attach scrollbars to the x- and y-axes
    self.xscrollbar = Scrollbar(self.root, orient=HORIZONTAL)
    self.yscrollbar = Scrollbar(self.root)
    self.box.config(height=20,
                    wrap=NONE,
                    yscrollcommand=self.yscrollbar.set,
		    xscrollcommand=self.xscrollbar.set)
    self.xscrollbar.config(command=self.box.xview)
    self.yscrollbar.config(command=self.box.yview)

    # Do not allow the user to edit the contents of the box
    self.box.config(state=DISABLED)

    # Provide pop-up help
    boxPopup = popUp.popUp(self.root, self.box, title='FrameLog window',
                           text="""In this window all new frames detected are displayed. The type
                                   of the detected frame will also be showed.
				""")

    # Display the box and the scrollbars in the window
    self.box.grid(row=1, column=0, columnspan=2, sticky=N+E+S+W)
    self.yscrollbar.grid(row=1, column=2, sticky=N+S) 
    self.xscrollbar.grid(row=2, column=0, columnspan=3, sticky=E+W)

    # Test a Listbox
    self.list = Listbox(self.root, selectmode=EXTENDED)
    self.list.insert(END,"una casa")
    self.box.insert(END,"dos cosas")
    self.list.insert(END,"tres cosas")
    #self.list.pack(fill=BOTH, expand=1)
    #self.current = None
    #self.poll() # start polling the list

    # Display the window
    self.root.grid(sticky=N+E+S+W)


  def update(self):

    "Read any pending messages from the queue, and display these in the window"

    while not log.empty():
      # Extract the message, and optional level/error arguments from queue
      (message, level, errorflag) = log.get()

      # Write the text in the text box
      self.write(message, level=level, errorflag=errorflag)


  def autoupdate(self):

    "Periodically update the window, because settings may change while on display"

    # Try to cancel the next call to this routine (if there is one scheduled)
    # In this way, there can never be more than one waiting loop active
    try: self.after_cancel(self.updateID)
    except: pass

    # Update
    self.update()

    # And reschedule itself in one second from now
    self.updateID = self.after(autocheckwaitingtime, self.autoupdate)


  def autologsave(self, tag_begin="1.0"):

    "Periodically save the log to a file"

    # Try to cancel the next call to this routine (if there is one scheduled)
    # In this way, there can never be more than one waiting loop active
    try: self.after_cancel(self.autologsaveID)
    except:
      if self.autologsaveID == None:
	# This was apparently the first time the routine was called,
	# so reinitialize the autosave logfile
	filehandle = open(autologsavefile, mode='w')
	filehandle.flush()
	filehandle.close()

    # Save log to file
    self.save(tag_begin, filename=autologsavefile, mode='a')
    tag_begin = self.box.index('end - 1 line')


    # And reschedule itself in one second from now
    self.autologsaveID = self.after(autologsavewait, self.autologsave, (tag_begin))


  def write(self, message, level=0, errorflag=0):

    "Write the provided message into the message log"


    # Only insert this message is it's 'level' is less than or equal to loglevel,
    if (level <= self.loglevel.get()):

       # Depending on  and whether this is a message, error or warning, set
       # the text color and separation strings. Different separation strings
       # make it possible to distinguish between messages and errors (warnings)
       # when colors are not available (as may be the case in saved log files).
      if errorflag == 1:
        textcolor='red'
	sepstring='***'
      elif errorflag > 1:
        textcolor='#339933' # which is a dark green color
	sepstring='-*-'
      else:
        textcolor='black'
	sepstring='---'

      # Create a text string to be displayed in the box
      # strftime() gives a nicely formatted time stamp
      textstring = "%s %s (%i) %s\n" % (strftime('%H:%M:%S'), sepstring, level, message)

      # Prepare to change the text in the box
      self.box.config(state=NORMAL)
      
      # Make the end of the log file visible
      self.box.yview(MOVETO, 1.0)

      # Find where the current text ends (= beginning of the new message)
      oldend = self.box.index('end - 1 line')

      # Insert the string      
      self.box.insert(END, textstring)
      
      # Find the end of the message in the box
      newend = self.box.index('end - 1 line')
      
      # Tag the message in the box and give it the appropriate color
      self.box.tag_add(textstring, oldend, newend)
      self.box.tag_config(textstring, foreground=textcolor)

      # Disable changes to the box
      self.box.config(state=DISABLED)


  def clear(self):

    "Clear the message log window. Gives option to save to file first."

    # Make sure the user really wants this
    answer = Dialog(self, title='Clear log file?', bitmap='', default=2,
           strings=('Clear log', 'Save to File', 'Abort'),
	   text='Do your really want to clear the log?\nData will be lost!')

    if answer.num == 0:
      # OK, delete all
      self.box.config(state=NORMAL)
      self.box.delete("1.0", END)
      self.box.config(state=DISABLED)
    if answer.num == 1:
      # Save to file first
      self.save()


  def save(self, tag_begin="1.0", tag_end=END, filename=None, mode='w'):

    "Save contents of the message log window to file"

    # Ask the user for an appropriate file to save to
    if not filename: filename = myDialog.SaveFile(self, initialdir=config.options['outpath'].value)
    
    # Still no filename means user cancelled -> return!
    if not filename: return

    if tag_end == END: tag_end = self.box.index('end - 1 line')

    filehandle = open(filename, mode=mode)
    filehandle.write(self.box.get(tag_begin, tag_end))
    filehandle.close()
