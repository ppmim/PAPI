################################################################################
#
# PANIClook
#
# pipeFrame.py
#
# Last update 27/02/2008
#
################################################################################

"""
   Provide a frame containing a list of tasks to perform. Tasks can be
   (de)selected manually. This module also contains routines to control the
   scheduling and execution of frame reductions using threads. Parts of this
   module will run as a child thread, in order to relieve load from the main
   thread, of which the main task is to keep the GUI updated.
"""

# Import external modules

import config
import tasks
import messageLog
import taskManager
import taskBar
import fileUtils
import os
import Queue
import threading

from Tkinter import *
from time import sleep


# Define constants for this module

# Available data reduction modes
#availableModes = ('QuickLook', 'Advanced', 'DoubleSpec')

# For each mode, define which tasks (from tasks.py) are to be included in
# the pipeline frame.

#tasknames = {}
#tasknames['QuickLook'] =	['autoload',  'preproc', 'subtractbias', 'flatfield',
#                                'plotcross', 'extspec',
#                                 'blazecorr', 'addwave']
#tasknames['Advanced']  =	['autoload', 'preproc', 'subtractbias',
#			         'scattering',  'flatfield', 'plotcross',
#			         'extspec',  'blazecorr',
#		         'addwave', 'ordermerge']
#tasknames['DoubleSpec']  =	['autoload', 'preproc', 'subtractbias',
#			         'scattering', 'flatfield', 'plotcross', 'getspecshift',
#		         'extspec',  'blazecorr',
#			         'addwave', 'adjustwave', #'ordermerge'
#			 ]


################################################################################

class Busy(Exception):

  """
     Custom error class, raised when the pipeline is busy and cannot accept
     new jobs. Generates no output
  """

  def __init__(self, *args):
    self.args = args


################################################################################

class window(Frame):

  """
     Provide a frame displaying the tasks that are to be executed for the
     current reduction mode defined in 'self.mode'.
  """

  def __init__(self, parent, GUIupdater, bgcolor=None):

    # 'parent' is the Tkinkter object under which this window is
    # hierarchially placed. 'GUIupdater' contains a Queue.Queue object into
    # which calls for updating the GUI may be placed (useful when running
    # as a thread). These variables are global in this class

    self.parent = parent
    self.GUIupdater = GUIupdater

    self.defaultbgcolor = bgcolor

    # Inialize as if this were an instance of the Frame class
    Frame.__init__(self, parent)
    # And make 'self' frame flexible in width
    self.columnconfigure(0, weight=1)

    # Define a dictionary that will contain the tasks which the user has
    # selected
    self.selectedtasks = {}
    
    # Container for taskBar objects
    self.taskbars = {}
    
    # Current list of tasks to display
    self.tasklist = []

    # Define a StringVar that will contain the current reduction mode
    self.mode = StringVar(self)

    # Take initial reduction mode from current configuration
    self.mode.set(config.options['currentmode'].value)

    # Build the actual frame
    self.makeFrame()

    # Create an Queue object that will contain the result of child threads
    self.outQueue = Queue.Queue()


  def makeFrame(self):

    "Construct the frame from scratch"

    # Put all data in a root frame that is embedded in the 'self' Frame object.
    # In this way, it is possible to destroy and rebuid the frame object,
    # without destroying the 'self' frame object.
    self.root = Frame(self, bg=self.defaultbgcolor)

    # The root frame should be flexible in width
    self.root.columnconfigure(0, weight=1)

    # Get the value of the current reduction mode and save in the config option
    self.currentmode = self.mode.get()
    config.options['currentmode'].set(self.currentmode)


    # Put the info on the tasks and modes into a variable that will be
    # 'unpacked' below
    modesandtasks = config.options['reductionmodes'].value

    self.tasklist = []
    task_is_selected = {}

    try:
      modesandtasks[self.currentmode]
    except KeyError:
      messageLog.put_error('No reduction modes are defined!')
      messageLog.put_error('You will not be able to do any reductions')
      messageLog.put_error('Update or reinstall your default configuration file')
      return

    # Obtain direct references to the class objects that are named in
    # 'tasknames' for the current mode, and if they are selected or not
    for (taskname, selected) in modesandtasks[self.currentmode]:
      self.tasklist.append( getattr(tasks, taskname) )
      task_is_selected[taskname] = selected

    # Row counter
    i = 1

    # Loop over tasks to create the taskBars
    for task in self.tasklist:

      # Create an integer that will be bound to the checkbutton in the taskBar
      self.selectedtasks[task.name] = IntVar(master=self.root)

      # Set the task to 'selected' or not
      if task_is_selected[task.name]:
        self.selectedtasks[task.name].set(1)

      # Create a 'taskBar' object, which consists of a checkbutton and a
      # status line. Store the taskBar in 'self.taskbars' for future reference
      self.taskbars[task.name] = taskBar.taskBar(self.root, task, self.selectedtasks[task.name], bgcolor=self.defaultbgcolor)
      self.taskbars[task.name].grid(sticky=E+W)

      # Next row
      i = i + 1

    # And display the contents of the root frame
    self.root.grid(sticky=E+W)


  def disable(self):

    "Make all displayed taskbars inactive"

    for task in self.tasklist:
      self.taskbars[task.name].disable()


  def enable(self):

    "Make all displayed taskbars active"

    for task in self.tasklist:
      self.taskbars[task.name].enable()


  def clear(self):

    "Clear the text of all displayed taskbars"

    for task in self.tasklist:
      self.taskbars[task.name].set_label()


  def reset(self):

    "Reset the contents of the 'self' frame"

    # Empty container objects
    self.selectedtasks = {}
    self.taskbars = {}
    self.tasklist = []

    # Destroy the root frame and rebuild it
    self.root.destroy()
    self.makeFrame()


  def process_frame(self, frame):

    """
       Prepare a child thread to do the reduction of a single frame. (Will
       run in main thread)
    """

    try:
      if self.pipethread.isAlive():
        # If busy reducing a frame, say so!
        raise Busy
    except AttributeError:
      # Will be raised if pipethread does not exist. This is OK, because
      # is means no reduction is going on
      pass


    # Make sure the input frames exists
    if not os.path.exists(frame):
      messageLog.put_error('Error processing %s' % frame)
      messageLog.put_error('Frame does not exist (anymore?). Skipped!')
      return

    # Put all the selected tasks (by their checkboxes) in the object 'todotasks'
    todotasks = []
    for task in self.tasklist:
      if self.selectedtasks[task.name].get():
	todotasks.append(task.name)

    # Start a child thread that will do the data 'crunching'
    messageLog.put('Starting pipeline thread', 9)
    self.pipethread = threading.Thread(None, self.run_pipe,
			 args=(frame, self.outQueue, todotasks, self.GUIupdater, self.currentmode))
    self.pipethread.start()


  def run_pipe(self, inframe, outQueue, todotasks, GUIupdater, mode):

    "Do the reduction a single frame. (Will run in child thread)"

    # Prepare local variables
    doneTasks = []
    outframe = ''
    exit = 0

    # Display an initial message in the logging window
    nchars = len('Processing %s' % (inframe))
    messageLog.put('='*nchars)
    messageLog.put('Processing %s' % (inframe))
    messageLog.put('='*nchars)

    # Remember original name of frame
    firstinframe = inframe

    # Split filename into path, base name and extension
    (path, base) = os.path.split(inframe)
    (outfilebase, extn)  = os.path.splitext(base)

    # Remove existing output files with same basename.
    fileUtils.removefiles(os.path.join(config.options['outpath'].value,
                                        '%s_step*.fits' % outfilebase))

    # Disable the contents of the 'self' frame
    GUIupdater.put((self.disable, (), {}))

    # Clear the old messages in the taskBars
    GUIupdater.put((self.clear, (), {}))

    # Can this safely be commented out? ONLY if the mode change is correctly
    # adjusted upon rebuilding the frame after a change of mode.
    #
    # Define a (new) configuration option that will contain the name
    # of the current reduction mode. This paramter will be referenced by the
    # IRAF-based tasks, and will search for the calibration files in the
    # corresponding directory "tasks/[mode]".
    #config.options['currentmode'].set(mode)

    # Keep track of number of reduction steps
    step = 0

    # Loop over the tasks
    for task in self.tasklist:

      # Increase task counter
      step = step + 1

      # Highlight the status bar corresponding to this task
      GUIupdater.put((self.taskbars[task.name].highlight, (), {}))

      # Skip the tasks if it wasn't in 'todotasks' (i.e. it was not selected
      # by the user).
      if not task.name in todotasks:
	GUIupdater.put((self.taskbars[task.name].set_label, (), {'text':'Skipped'}))
	GUIupdater.put((self.taskbars[task.name].highlight_off, (), {}))
	continue

      # Give the GUIupdater time to adjust the screen
      if task.inthread:
        while not GUIupdater.empty(): sleep(0.1)

      # Make sure any previous reduction steps required by the current task
      # have been performed before
      for neededTask in task.prereq:
        if not neededTask in doneTasks:
	  # If the reduction cannot continue, make it clear to the user
	  messageLog.put_error('Cannot execute "%s" without "%s"'
	            % (task.buttonText, getattr(tasks, neededTask).buttonText))
          GUIupdater.put((self.taskbars[task.name].set_label, (), {'text':'Stopped'}))
	  GUIupdater.put((self.taskbars[task.name].highlight, (), {'color':'red'}))
          exit = 1 # And make sure we exit the reduction loop (only 'break'
	           # would not be enough here)

      # Exit if needed
      if exit: break

      # Create name of output frame from step counter and step suffix
      if task.suffix:
	outframe = os.path.join(config.options['outpath'].value,
	                        '%s_step%03i_%s.fits' % (outfilebase, step, task.suffix))
      else:
	outframe = None

      # OK to start to run the task; reflect this in the taskBar label and log
      newlabel = 'Running "%s" with %s' % (task.buttonText, os.path.basename(inframe))
      GUIupdater.put((self.taskbars[task.name].set_label,
                     (), {'text':newlabel} ))
      messageLog.put(newlabel, 7)

      # Call the taskManager to do the task executing and retrieve possible
      # errors encountered. Any errors occuring while executing the task will
      # never end the reduction interface. Instead, a (descriptive?) error
      # will be shown in the log

      ######################################################################
      errorstring = taskManager.runTask(task, inframe, outframe)
      ######################################################################

      # Errors encountered? If yes, show them!
      if errorstring != None:
        if errorstring == '':
          GUIupdater.put((self.taskbars[task.name].set_label, (), {'text':'No further processing of this frame'}))
	  GUIupdater.put((self.taskbars[task.name].highlight, (), {'color':'orange'}))
	else:
	  GUIupdater.put((self.taskbars[task.name].set_label, (), {'text':errorstring}))
	  GUIupdater.put((self.taskbars[task.name].highlight, (), {'color':'red'}))
 	break
      else:
        # Everything went smoothly
        GUIupdater.put((self.taskbars[task.name].set_label, (), {'text':'Finished'}))
        GUIupdater.put((self.taskbars[task.name].highlight_off, (), {}))

      # Append task name to the list of successfully finished tasks
      doneTasks.append(task.name)

      # Only if the task we just executed returned an output file name, the
      # next task should use this as its input file name.
      if outframe and os.path.exists(outframe) and task.output:
	inframe = outframe
      else: outframe = inframe


    # After the loop over all tasks terminates, fill the outQueue, so that the
    # calling routine knows the reduction of the current frame is finished.
    self.outQueue.put((firstinframe, outframe))

    # Re-enable the 'self' frame
    GUIupdater.put((self.enable, (), {}))

    # And bring the good news to the user
    messageLog.put('Finished processing %s' % os.path.basename(firstinframe))



