################################################################################
#
# PANICtool
#
# calibFrame.py
#
# Last update 16/04/2008
#
################################################################################

"""
   Interactively calculate calibration frames (master biases, flats etc.) and
   definition files such as order locations and wavelength solutions.
"""

# Import external modules

import config
import tasks
import messageLog
import taskManager
import configFrame
import threading
import taskBar

from time import sleep
from Tkinter import *

################################################################################

# Define local constants

# There is currently only one available mode when calculating calibration
# frames and other related settings. Maybe more mode might be implemented in
# the future, but at the moment this parameter is only here to provide a
# parallel naming to the naming in 'pipeFrame'
availableModes = ('CalibCalc',)

# Include the following tasks from the 'tasks' module in the widget
# The 'run' method of these tasks should not expect any parameters
#tasknames = ['sumbias',    'sumflat',  'findord', 'findinterlacedord',
#             'plotorders', 'normflat', 'wavecal', 'interlacedwavecal', 'sumdarks' , 'subtractSky1','subtractSkyN' ]

tasknames = ['createMasterDark', 'createMasterDomeFlat', 'createMasterSkyFlat', 'createPixelMask','subtractSky1','subtractSkyN',"test","test2" ] 

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
     Provide a window to execute the calibration tasks, and set essential
     parameters, such as the location of input and output calibration frames.
     Class inherits from Tkinter's 'Frame' class.
  """

  def __init__(self, parent, GUIupdater, title=''):

    # Parent window and GUIupdater queue object are global variables
    # in this class. Parent is a Tkinter.Frame object, and GUIupdater a
    # Queue.Queue object. 
    self.parent = parent
    self.GUIupdater = GUIupdater

    # Set the current reduction mode (similar structure as in pipeFrame)
    self.currentmode = availableModes[0]

    # self.top is a frame under 'parent'
    self.top = Toplevel(parent)

    # Start this frame as hidden and without title
    self.top.withdraw()
    self.top.title(title)

    # Make this frame flexible in width
    self.top.columnconfigure(0, weight=1)
    self.top.rowconfigure(0, weight=1)

    # Clicking 'close window' on this frame calls the self.hide method
    self.top.protocol("WM_DELETE_WINDOW", self.hide)

    # Add a scrollbar
    vscrollbar = AutoScrollbar(self.top)
    vscrollbar.grid(row=0, column=1, sticky=N+S)

    # Create a Canvas object that can contain the whole frame
    self.canvas = Canvas(self.top, yscrollcommand=vscrollbar.set)
    self.canvas.columnconfigure(0, weight=1)
    self.canvas.rowconfigure(0, weight=1)

    # Attach the method to position the canvas to the scrollbar
    vscrollbar.config(command=self.canvas.yview)

    # Inialize as if this were an instance of the Frame class.
    Frame.__init__(self, self.canvas)
    # And make also the 'self' frame flexible in width
    self.columnconfigure(0, weight=1)

    # Initialize self.taskbars, a dictionary containing the task-bars
    # corresponding to each task. (Useful for changing its contents)
    self.taskbars = {}

    # Now construct the frame
    self.makeFrame()
    # Update with latest information 

    self.update()

    # Store it in a canvas object (to allow scrolling)
    self.canvas.create_window(0,0, anchor=NW, window=self)

    # Show frame (but it remains hidden)
    #self.grid(sticky=E+W)


  def makeFrame(self):
  
    "Construct the frame from scratch"

    # Put all data in a root frame that is embedded in the 'self' Frame object.
    # In this way, it is possible to destroy and rebuid the frame object,
    # without destroying the 'self' frame object.

    self.root = Frame(self)
    
    # Make also the root frame flexible in width
    self.root.columnconfigure(0, weight=1)

    # Row counter
    i = 0

    # Display a 'window' instance of the 'configFrame' class for a number
    # of input options.
    optionlist = ('inpath', 'filename_filter', 'darklist', 'flatlist', 'sciencelist')
    configFrame.window(self, optionlist=optionlist, title='Input options', embed=1).grid(row=i, stick=E+W)

    # Next row
    i = i + 1

    # Do the same for the output options
    optionlist = ('outpath', 'masterflat', 'masternormflat', 'mastersky')
    configFrame.window(self, optionlist=optionlist, title='Output options', embed=1).grid(row=i, stick=E+W)

    # Next row
    i = i + 1

    # Display a header label
    label = Label(self, text='Calculate calibration frames :')
    label.grid(row=i, columnspan=2, sticky=E+W)

    # Next row
    i = i + 1

    # Create a list containing each task mentioned in 'tasknames'
    self.tasklist = [ getattr(tasks, taskname) for taskname in tasknames ]
    
    # For each task, create a 'taskBar' object, which consists of (1) a
    # button to execute the task and (2) a status line.
    for task in self.tasklist:

      # Store in self.taskbars the taskBar instance that correspond to this
      # task and attach the 'runTask' method to the button.
      self.taskbars[task.name] = taskBar.taskBar(self.root, task)
      self.taskbars[task.name].button.config(command = lambda t = task: self.runTask(t))
      self.taskbars[task.name].grid(sticky=E+W)

      # Next row
      i = i + 1


    # At the bottom, display an 'OK' button that will close (hide) the window
    self.okbutton = Button(self.root, text='OK', command=self.hide)
    self.okbutton.grid(row=i, columnspan=2, sticky=E+W)

    # Display the contents of the root frame
    self.root.grid(sticky=E+W)


  def show(self):

    "Deiconify and set focus and grab to this window"

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

    self.top.grab_release()
    self.top.withdraw()


  def update(self):

    "Update the contents in the window"

    # Nothing needed for the moment. (Frame is always updated). Method may be
    # called by external (higher level) routines, so don't remove this method.
    pass


  def disable(self):

    "Disable the widget elements of this window"

    # The only things that makes sense to disable in this window are the
    # taskbars and the 'OK' button
    for taskbar in self.taskbars.values():
      taskbar.disable()

    self.okbutton.config(state=DISABLED)


  def enable(self):

    "Enable the widget elements of this window"

    # The only things that makes sense to enable in this window are the
    # taskbars and the 'OK' button
    for taskbar in self.taskbars.values():
      taskbar.enable()

    self.okbutton.config(state=NORMAL)
  

  def runTask(self, task):

    "Wrapper to run one given task"

    # Wrap into a thread if the tasks allows that
    if task.inthread:
      taskthread = threading.Thread(None, self._doTask,
                                args=(task, self.GUIupdater, self.currentmode))
      taskthread.start()
    else:
      # No thread, so run in foreground instead
      self._doTask(task, self.GUIupdater, self.currentmode)


  def _doTask(self, task, GUIupdater, mode):

    """
       Do the actual task execution. This method is called by runTask, and
       should not be called directly.
    """

    # Highlight the taskbar (or better: ask the GUIupdater to do that for me)
    GUIupdater.put((self.taskbars[task.name].highlight, (), {}))
    # Adjust the taskbar label
    GUIupdater.put((self.taskbars[task.name].set_label, (), {'text':'Executing "%s"' % task.buttonText}))

    # Disable the elements of this window (prohibit the user from messing around)
    GUIupdater.put((self.disable, (), {}))

    # Just before calling the task, save the current mode to the 'currentmode'
    # contifuration option.
    oldmode = config.options['currentmode'].value
    config.options['currentmode'].set(mode)

    # Give the GUIupdater time to adjust the screen
    if task.inthread:
      while not GUIupdater.empty(): sleep(0.1)

    # Finally... ask the taskManager to execute this task
    ##########################################################
    errorstring = taskManager.runTask(task)
    ##########################################################

    if errorstring:
      # If taskManager reported an error, display it in the taskbar status
      GUIupdater.put((self.taskbars[task.name].set_label, (), {'text':errorstring}))
      GUIupdater.put((self.taskbars[task.name].highlight, (), {'color':'red'}))
    else:
      # Normal termination of the task
      GUIupdater.put((self.taskbars[task.name].set_label, (), {'text':'Finished'}))
      GUIupdater.put((self.taskbars[task.name].highlight_off, (), {}))

    # Restore the 'currentmode' parameter to its old value
    config.options['currentmode'].set(oldmode)

    # Re-enable the elements of this window (restore user interaction)
    GUIupdater.put((self.enable, (), {}))

