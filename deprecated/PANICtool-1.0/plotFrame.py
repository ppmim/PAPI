################################################################################
#
# FIEStool
#
# plotFrame.py
#
# Last update 04/10/2005
#
################################################################################

"""
   Provide an interface to automatically or manually plot a reduced frame. 
"""

# Import external modules

import os
import config
import configFrame
import tasks
import taskManager
import myDialog
import messageLog
import popUp

from Tkinter import *
import tkFileDialog

################################################################################

class window(Frame):

  """
     Provide a window to manage the automatic or manual plotting of spectra. The
     user is able to provide axis boundaries, and can choose to use the
     quicklook plotting program 'Biggles', or the more extended IRAF interface.
  """

  def __init__(self, parent):

    # Inialize as if this were an instance of the Frame class.
    Frame.__init__(self, parent)

    # Make this frame flexible in width
    self.columnconfigure(0, weight=1)

    # Define a variable that will contain the name of the last plotted file
    self.lastplotname = None

    # Define a Tkinter integer variable that will be bound to a radiobutton
    self.usebiggles = IntVar(master=self)
    self.usebiggles.set(1)

    # Now construct the frame
    self.makeFrame()


  def makeFrame(self):

    "Construct the frame from scratch"

    # Put all data in a root frame that is embedded in the 'self' Frame object.
    # In this way, it is possible to destroy and rebuid the frame object,
    # without destroying the 'self' frame object.

    self.root = Frame(self)

    # Give the root frame a minimum size, and make it flexible in the
    # middle column (like other frames in this interface)
    self.root.columnconfigure(0, weight=1, minsize=200)
    self.root.columnconfigure(1, weight=1)

    # Row counter
    i = 1



    # Create a button that will invoke the 'plotSpec' method and  plot the
    # spectrum with Biggles. Provide an on-line help text for this button.
    replotbutton = Button(self.root, text='Replot last spectrum', command=self.plotSpec)
    replotbutton.grid(row=i, column=0, columnspan=2, sticky=E+W)
    replotbuttonPopup = popUp.popUp(self.root, replotbutton, title='Replot last spectrum',
                          text="""Replot the last reduced spectrum. You will be prompted
			          for a file if you have not reduced any spectra yet.
			       """)


    # Next row
    i = i + 1


    # Create a radiobutton for the option to plot spectra with the simple
    # 'Biggles' plotting package. Provide an on-line help text for this option.
    biggButton = Radiobutton(self.root, width=25, text='Plot with Biggles',
             variable=self.usebiggles, value=1, indicatoron=1, pady=4)
    biggButton.grid(row=i, column=0, sticky=E+W)
    biggButtonPopup = popUp.popUp(self.root, biggButton, title='Biggles',
                          text="""Biggles is a simple plotting package,
                                  producing non-interactive plots. Biggles'
                                  plots appear in popup windows that can be
                                  closed by clicking on them.
                               """)

    # Create another radiobutton for the option to plot spectra with the
    # 'splot'-function of IRAF. The call to this function will go through
    # PyRAF. Provide an on-line help text for this option.
    irafButton = Radiobutton(self.root, width=25, text='Plot with IRAF',
             variable=self.usebiggles, value=0, indicatoron=1, pady=4)
    irafButton.grid(row=i, column=1, sticky=E+W)
    irafButtonPopup = popUp.popUp(self.root, irafButton, title='IRAF plots',
                          text="""IRAF plots are produced with the 'splot' routine from
                                  the echelle reduction package. The plots appear in a popup window
                                  that can be closed by pressing 'q' with the mouse in
                                  the window. There are a large number of simple spectrum
                                  analysis tools available (press '?' for help).
                               """)

    # Next row
    i = i + 1


    # Create a configuration button for the starting wavelength
    wavebutton1 = configFrame.window(self.root, optionlist=['plot_startwave'],
                                     embed=True, width=20)
    wavebutton1.grid(row=i, column=0, sticky=E+W)


    # Idem for the ending wavelength
    wavebutton2 = configFrame.window(self.root, optionlist=['plot_endwave'],
                                     embed=True, width=20)
    wavebutton2.grid(row=i, column=1, sticky=E+W)

    # Next row
    i = i + 1


    # And a configuration button for plotting individual orders (not very common)
    orderbutton = configFrame.window(self.root, optionlist=['plot_defaultorder'],
                                     embed=True, width=20)
    orderbutton.grid(row=i, column=0, sticky=E+W)
    Label(self.root, anchor=W, text='(if no wavelengths available)').grid(row=i, column=1, stick=E+W)

    # Next row
    i = i + 1


    # Create a button that will invoke the 'plotSpec' method and plot the
    # spectrum with Biggles. Provide an on-line help text for this button.
    plotbutton = Button(self.root, text='Plot other spectrum', command=lambda:self.plotSpec(plotother=True) )
    plotbutton.grid(row=i, column=0, columnspan=2, sticky=E+W)
    plotbuttonPopup = popUp.popUp(self.root, plotbutton, title='Plot other spectrum',
                          text="""Plot a reduced spectrum using Biggles or IRAF. You will
			          be prompted for a file if 'Use default filename' is not
				  selected (or no valid filename is listed there).
			       """)


    # Next row
    i = i + 1


#    NOT YET IMPLEMENTED. DO WE REALLY WANT TO PROVIDE HARDCOPIES,
#    OR IS THAT NOT THE TASK OF QUICKLOOK REDUCTION SOFTWARE?
#
#    hardcopyButton = Button(self.root, text='Send last plot to printer',
#                           command=lambda: self.printSpec)
#    hardcopyButton.grid(row=i, column=1, sticky=E+W)
#    harccopyButtonPopup = popUp.popUp(self.root, hardcopyButton, title='Biggles',
#			  text="""Send the latest Biggles plot to the printer
#			       """)
#
#    i = i + 1




    # Finally, display the entire frame (that is, the 'root' frame)
    self.root.grid(sticky=E+W)


  def plotSpec(self, filename=None, plotother=False):

    """
       Determine the current state of the widget entries and plot the
       spectrum accordingly.
    """

    # Should I use the default name?
    if plotother:
      # No - put filename to None
      filename = None
    else:
      # Yes - read it from the entry
      filename = self.lastplotname


    # Now, if no filename is defined yet, open a dialog to ask for a file to be
    # plotted
    if not filename:
      filename = tkFileDialog.askopenfilename(parent=self, initialdir=config.options["outpath"].value)
      if filename: self.lastplotname = filename
      else: return


#   Not yet implemented :
#   LOOP OVER RANGES - WHEN THERE ARE SEVERAL RANGES DEFINED - DO WE WANT THIS?


    # Check that the file to plot (still) exists and then ask the task manager
    # to execute the 'plotspec' task.
    if os.path.exists(filename):
      taskManager.runTask(tasks.plotspec, filename, usebiggles=self.usebiggles.get())
    else:
      messageLog.put_warning('Could not open file for plotting - Skipped (%s)' % filename)



# NOT YET PROPERLY IMPLEMENTED...
  def printSpec(self):

    try:
      self.lastplot.write_eps('tmp.eps')
      messageLog.put('Sending latest plot to the printer')
    except:
      raise
