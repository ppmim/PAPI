################################################################################
#
# FIEStool
#
# myDialog.py
#
# Last update 29/09/2005
#
################################################################################

"""
   A not-so-nice module that modifies the file dialog methods in FileDialog
   (part of the Tkinter library distributed with Python) to contain fields for
   selected FITS headers. Also a new dialog for selecting several files at once
   is implemented. FileDialog itself is a badly comment module, so my comments
   are limited to the least necessary.
   
   Module PyFITS provides the interface to read FITS headers and is needed to
   use this module.
"""

# Import external modules

from Tkinter import *

import tkSimpleDialog
import tkFileDialog
import pyfits
import os

from FileDialog import *


################################################################################


class LoadFileDialog(FileDialog):

    """File selection dialog which checks that the file exists."""

    # Redefinition of LoadFileDialog in FileDialog.
    #
    # 'set_selection' now also updates the highlighted selection on screen
    # See also comments in LoadFilesDialog.go()

    title = "Select file"

    def ok_command(self):
        file = self.get_selection()
        if not os.path.isfile(file):
            self.master.bell()
        else:
            self.quit(file)

    def set_selection(self, file):
        self.selection.delete(0, END)
        self.selection.insert(END, os.path.join(self.directory, file))
        try:
          dirlist = list(self.files.get(0, END))
	  self.files.select_set(dirlist.index(file))
	except:
	  pass


################################################################################


class SaveFileDialog(FileDialog):

    """File selection dialog which checks that the file may be created."""

    # See comments in LoadFileDialog

    title = "Select File"

    def ok_command(self):
        file = self.get_selection()
#        if os.path.exists(file):
#            if os.path.isdir(file):
#                self.master.bell()
#                return
#            d = Dialog(self.top,
#                       title="Overwrite Existing File Question",
#                       text="Overwrite existing file %s?" % `file`,
#                       bitmap='questhead',
#                       default=1,
#                       strings=("Yes", "Cancel"))
#            if d.num != 0:
#                return
#        else:
#            head, tail = os.path.split(file)
#            if not os.path.isdir(head):
#                self.master.bell()
#                return
        self.quit(file)

    def set_selection(self, file):
        self.selection.delete(0, END)
        self.selection.insert(END, os.path.join(self.directory, file))
        try:
          dirlist = list(self.files.get(0, END))
	  self.files.select_set(dirlist.index(file))
	except:
	  pass


################################################################################


class LoadFilesDialog(FileDialog):

    """Multiple file selection dialog"""

    # Similar to LoadFileDialog, but will allow several files to be selected
    # at the same time. In this class, 'self.selection' is used in a different
    # way as in FileDialog. Here it can contain several files, and therefore
    # this widget element is not displayed on screen.

    title = "Select files"

    def __init__(self, master, title="Select files"):
        FileDialog.__init__(self, master, title=title)
	# Remove the selection bar from the screen.
	self.selection.pack_forget()
	self.files.config(selectmode = MULTIPLE)

    def go(self, dir_or_file=os.curdir, pattern="*", default="", key=None):
        # Remove the path from the default file names when asking for files
        if key and dialogstates.has_key(key):
            self.directory, pattern = dialogstates[key]
        else:
            dir_or_file = os.path.expanduser(dir_or_file)
            if os.path.isdir(dir_or_file):
                self.directory = dir_or_file
            else:
                self.directory, default = os.path.split(dir_or_file)
        # Set default selection to basenames of files
        default = [os.path.basename(file) for file in default]
        self.set_filter(self.directory, pattern)
        # The following two lines have been interchanged, because there is
	# no reason to set default before filtering
        self.filter_command()
        self.set_selection(default)
        self.selection.focus_set()
        self.top.grab_set()
        self.how = None
        self.master.mainloop()          # Exited by self.quit(how)
        if key:
            directory, pattern = self.get_filter()
#           if self.how:
#               directory = os.path.dirname(self.how)
            dialogstates[key] = directory, pattern
        self.top.destroy()
        return self.how

    def ok_command(self):
        file = [self.files.get(x) for x in self.files.curselection()]
        directory, pattern = self.get_filter()
	file = [os.path.join(directory, x) for x in file]
	self.quit(file)

    def set_selection(self, files):
        try:
          dirlist = list(self.files.get(0, END))
	  for file in files:
	    self.files.select_set(dirlist.index(file))
	except:
	  pass


###############################################################################


class FITSFileDialog(FileDialog):

    """File selection dialog with additional field containing FITS header"""

    # Very similar to FileDialog, except for the additional fields and the
    # difference in treating 'self.selection'

    def __init__(self, master, title=None, headlist=None):
    
        if title is None: title = self.title
        self.master = master

	# Type-check the headlist parameter and convert its contents
	# to capitalized single words
	if type(headlist) is type(''):
	  headlist = headlist.upper()
	  headlist = headlist.split()
	  headlist.reverse()
	self.headlist = tuple(headlist)

        self.directory = None

        self.top = Toplevel(master)
        self.top.title(title)
        self.top.iconname(title)

        self.botframe = Frame(self.top)
        self.botframe.pack(side=BOTTOM, fill=X)

        self.selection = Entry(self.top)
        self.selection.pack(side=BOTTOM, fill=X)
        self.selection.bind('<Return>', self.ok_event)

        self.filter = Entry(self.top)
        self.filter.pack(side=TOP, fill=X)
        self.filter.bind('<Return>', self.filter_command)

        self.headframe = Frame(self.top)
	self.headframe.pack(expand=NO, fill=BOTH)

        self.midframe = Frame(self.top)
        self.midframe.pack(expand=YES, fill=BOTH)

        self.filesbar = Scrollbar(self.midframe)
        self.filesbar.pack(side=RIGHT, fill=Y)

        # Similar to the way the directory and file lists are managed,
	# create a container for the FITS header related elements.
        self.fitshead = {}
	
	# Add the header fields and bind events to related methods
	for header in self.headlist:

	  self.fitshead[header] = Listbox(self.midframe, exportselection=0,
                               yscrollcommand=self.filesbar.set, height=20)
          self.fitshead[header].pack(side=RIGHT, expand=YES, fill=BOTH)
#         btags = self.fitshead[header].bindtags()
#         self.fitshead[header].bindtags(btags[1:] + btags[:1])
          self.fitshead[header].bind('<ButtonRelease-1>', self.fitshead_select_event)
#         self.fitshead[header].bind('<Double-ButtonRelease-1>', self.fitshead_double_event)
          self.fitshead[header].bind('<Button-4>', lambda e, s=self: s.adjust_scroll(SCROLL, -1, UNITS))
          self.fitshead[header].bind('<Button-5>', lambda e, s=self: s.adjust_scroll(SCROLL,  1, UNITS))

	  Label(self.headframe, text=str(header)).pack(side=RIGHT, expand=YES, fill=BOTH)


        self.files = Listbox(self.midframe, exportselection=0,
                             yscrollcommand=self.filesbar.set)
        self.files.pack(side=RIGHT, expand=YES, fill=BOTH)
#       btags = self.files.bindtags()
#       self.files.bindtags(btags[1:] + btags[:1])
        self.files.bind('<ButtonRelease-1>', self.files_select_event)
#       self.files.bind('<Double-ButtonRelease-1>', self.files_double_event)
        self.files.bind('<Button-4>', lambda e, s=self: s.adjust_scroll(SCROLL, -1, UNITS))
        self.files.bind('<Button-5>', lambda e, s=self: s.adjust_scroll(SCROLL,  1, UNITS))

	Label(self.headframe, text='Filename').pack(side=RIGHT, expand=YES, fill=BOTH)

        self.filesbar.config(command=self.adjust_scroll)

        self.dirsbar = Scrollbar(self.midframe)
        self.dirsbar.pack(side=LEFT, fill=Y)
        self.dirs = Listbox(self.midframe, exportselection=0,
                            yscrollcommand=(self.dirsbar, 'set'))
        self.dirs.pack(side=LEFT, expand=YES, fill=BOTH)
        self.dirsbar.config(command=(self.dirs, 'yview'))
        btags = self.dirs.bindtags()
        self.dirs.bindtags(btags[1:] + btags[:1])
        self.dirs.bind('<ButtonRelease-1>', self.dirs_select_event)
        self.dirs.bind('<Double-ButtonRelease-1>', self.dirs_double_event)

	Label(self.headframe, text='Directory').pack(side=LEFT, expand=YES, fill=BOTH)


        self.ok_button = Button(self.botframe,
                                 text="OK",
                                 command=self.ok_command)
        self.ok_button.pack(side=LEFT)
        self.filter_button = Button(self.botframe,
                                    text="Filter",
                                    command=self.filter_command)
        self.filter_button.pack(side=LEFT, expand=YES)
        self.cancel_button = Button(self.botframe,
                                    text="Cancel",
                                    command=self.cancel_command)
        self.cancel_button.pack(side=RIGHT)

        self.top.protocol('WM_DELETE_WINDOW', self.cancel_command)
        # XXX Are the following okay for a general audience?
        self.top.bind('<Alt-w>', self.cancel_command)
        self.top.bind('<Alt-W>', self.cancel_command)


    def adjust_scroll(self, *args):
        apply(self.files.yview, args)
	for header in self.headlist:
          apply(self.fitshead[header].yview, args)
        return "break"

    # New!
    def getFITShead(self, file, cardname):
        cardval = ""
	try:
	  img = pyfits.open(file)
	  cardval = img[0].header[cardname]
	  img.close()
	except: pass
	return cardval

    def filter_command(self, event=None):
        dir, pat = self.get_filter()
        try:
            names = os.listdir(dir)
        except os.error:
            self.master.bell()
            return
        self.directory = dir
        self.set_filter(dir, pat)
        names.sort()
        subdirs = [os.pardir]
        matchingfiles = []
	fullmatchingfiles = []
        for name in names:
            fullname = os.path.join(dir, name)
            if os.path.isdir(fullname):
                subdirs.append(name)
            elif fnmatch.fnmatch(name, pat):
                matchingfiles.append(name)
                fullmatchingfiles.append(fullname)
        self.dirs.delete(0, END)
        for name in subdirs:
            self.dirs.insert(END, name)
        self.files.delete(0, END)
	for header in self.headlist:
            self.fitshead[header].delete(0, END)
        for name in matchingfiles:
            self.files.insert(END, name)
	if self.headlist is not None:
	  for header in self.headlist:
            for name in fullmatchingfiles:
		self.fitshead[header].insert(END, self.getFITShead(name, header))
        head, tail = os.path.split(self.get_selection())
        if tail == os.curdir: tail = ''
        self.set_selection(tail)

#   # New! Similar to the existing dirs_double_event or files_double_event
#   def fitshead_double_event(self, event):
#       self.ok_command()

    # New! Similar to the existing dirs_select_event or files_select_event
    def fitshead_select_event(self, hdr):
        fidx = [ int(x) for x in self.files.curselection() ]
	for header in self.headlist:
          idx = [ int(x) for x in self.fitshead[header].curselection() ]
	  if idx != fidx: break
	self.files.select_clear(0, END)
	for i in idx: self.files.select_set(i)
        for header in self.headlist:
	  self.fitshead[header].select_clear(0, END)
	  for i in idx: self.fitshead[header].select_set(i)
        try:
          file = self.files.get(idx[0])
          self.set_selection(file)
	except IndexError: pass

    # Make sure that also the correct FITS headers are highlighted
    def files_select_event(self, event):
        idx = [ int(x) for x in self.files.curselection() ]
        for header in self.headlist:
	  self.fitshead[header].select_clear(0, END)
	  for i in idx: self.fitshead[header].select_set(i)
        try:
          file = self.files.get(idx[0])
          self.set_selection(file)
	except IndexError: pass


################################################################################


class LoadFITSFileDialog(FITSFileDialog):

    """File selection dialog which checks that the file exists."""

    # Similar to LoadFileDialog, but also display FITS headers

    title = "Select file"

    def ok_command(self):
        file = self.get_selection()
        if not os.path.isfile(file):
            self.master.bell()
        else:
            self.quit(file)

    def set_selection(self, file):
        self.selection.delete(0, END)
        self.selection.insert(END, os.path.join(self.directory, file))
        try:
          dirlist = list(self.files.get(0, END))
	  self.files.select_set(dirlist.index(file))
	  for header in self.headlist: self.fitshead[header].select_set(dirlist.index(file))
	except:
	  pass
  

class SaveFITSFileDialog(FITSFileDialog):

    """File selection dialog which checks that the file may be created."""

    # Similar to SaveFileDialog, but also display FITS files
    #
    # This dialog will be used to select a file to save to, which, of course
    # may be an non-existing file. The file will only be overwritten in a later
    # stage (by the reduction routines), so double-checking if the file exists
    # is eliminated

    title = "Select File"

    def ok_command(self):
        file = self.get_selection()
#
# (see comment above)
#
#        if os.path.exists(file):
#            if os.path.isdir(file):
#                self.master.bell()
#                return
#            d = Dialog(self.top,
#                       title="Overwrite Existing File Question",
#                       text="Overwrite existing file %s?" % `file`,
#                       bitmap='questhead',
#                       default=1,
#                       strings=("Yes", "Cancel"))
#            if d.num != 0:
#                return
#        else:
#            head, tail = os.path.split(file)
#            if not os.path.isdir(head):
#                self.master.bell()
#                return
        self.quit(file)

    def set_selection(self, file):
        self.selection.delete(0, END)
        self.selection.insert(END, os.path.join(self.directory, file))
        try:
          dirlist = list(self.files.get(0, END))
	  self.files.select_set(dirlist.index(file))
	  for header in self.headlist: self.fitshead[header].select_set(dirlist.index(file))
	except:
	  pass


################################################################################


class LoadFITSFilesDialog(FITSFileDialog):

    """Multiple file selection dialog"""

    # Similar to LoadFileDialog, but also display FITS headers and allow
    # multiple files to be selected at the same time.

    title = "Select files"

    def __init__(self, master, title=None, headlist=None):
        FITSFileDialog.__init__(self, master, title=title, headlist=headlist)
	# Remove the selection bar from the screen.
	self.selection.pack_forget()
	self.files.config(selectmode = MULTIPLE)
	for header in self.headlist:
	  self.fitshead[header].config(selectmode = MULTIPLE)

        self.selectall_button = Button(self.botframe,
                                    text="Select All",
                                    command=self.select_all)
        self.selectall_button.pack(side=RIGHT, expand=YES)
        self.selectnone_button = Button(self.botframe,
                                    text="Clear selection",
                                    command=self.select_none)
        self.selectnone_button.pack(side=RIGHT, expand=YES)


    def go(self, dir_or_file=os.curdir, pattern="*", default="", key=None):
        # Remove the path from the default file names when asking for files
        if key and dialogstates.has_key(key):
            self.directory, pattern = dialogstates[key]
        else:
            dir_or_file = os.path.expanduser(dir_or_file)
            if os.path.isdir(dir_or_file):
                self.directory = dir_or_file
            else:
                self.directory, default = os.path.split(dir_or_file)
        # Set default selection to basenames of files
        default = [os.path.basename(file) for file in default]
        self.set_filter(self.directory, pattern)
        # The following two lines have been interchanged, because there is
	# no reason to set default before filtering
        self.filter_command()
        self.set_selection(default)
        self.selection.focus_set()
        self.top.grab_set()
        self.how = None
        self.master.mainloop()          # Exited by self.quit(how)
        if key:
            directory, pattern = self.get_filter()
#           if self.how:
#               directory = os.path.dirname(self.how)
            dialogstates[key] = directory, pattern
        self.top.destroy()
        return self.how

    def ok_command(self):
        # Return the selected file(s) with full path name
        file = [self.files.get(x) for x in self.files.curselection()]
        directory, pattern = self.get_filter()
	file = [os.path.join(directory, x) for x in file]
	self.quit(file)

    # Make sure also the FITS header fields get selected
    def set_selection(self, files):
        try:
          dirlist = list(self.files.get(0, END))
	  for file in files:
	    self.files.select_set(dirlist.index(file))
	    for header in self.headlist: self.fitshead[header].select_set(dirlist.index(file))
	except:
	  pass

    # Make sure also the FITS header fields get selected
    def select_all(self):
        try:
          self.files.select_set(0, END)
	  for header in self.headlist: self.fitshead[header].select_set(0, END)
	except IndexError: pass

    # Make sure also the FITS header fields get selected
    def select_none(self):
        try:
          self.files.select_clear(0, END)
	  for header in self.headlist: self.fitshead[header].select_clear(0, END)
	except IndexError: pass



###############################################################################


# Some quick macros that simplify interaction with the above routines
# (Similar to definitions in FileDialog)


def LoadFile(parent=None, title=None, initialdir=os.curdir, **args):
    result = tkFileDialog.askopenfilename(parent=parent, title=title, initialdir=initialdir, **args)
    if not result: result = None
    else: result=str(result)
    return result

#    fd = LoadFileDialog(parent, title=title)
#    loadfile = fd.go(dir_or_file=initialdir, default=default, pattern=pattern, key=key)
#    return loadfile

def LoadFiles(parent=None, title=None, initialdir=os.curdir, **args):
    result = tkFileDialog.askopenfilenames(parent=parent, title=title, initialdir=initialdir, **args)
    if not result: result = None
    else: result=str(result)
    return result

#    fd = LoadFilesDialog(parent, title=title)
#    loadfiles = fd.go(dir_or_file=initialdir, default=default, pattern=pattern, key=key)
#    return loadfiles

def SaveFile(parent=None, title=None, initialdir=os.curdir, **args):
    result = tkFileDialog.asksaveasfilename(parent=parent, title=title, initialdir=initialdir, **args)
    if not result: result = None
    else: result=str(result)
    return result

#    fd = SaveFileDialog(parent, title=title)
#    savefile = fd.go(dir_or_file=initialdir, default=default, pattern=pattern, key=key)
#    return savefile


def LoadFITSFile(parent=None, title=None, initialdir=os.curdir, default="", pattern="*.fit*", headlist="OBJECT", key=None):
    fd = LoadFITSFileDialog(parent, title=title, headlist=headlist)
    loadfile = fd.go(dir_or_file=initialdir, default=default, pattern=pattern, key=key)
    return loadfile

def LoadFITSFiles(parent=None, title=None, initialdir=os.curdir, default="", pattern="*.fit*", headlist="OBJECT", key=None):
    fd = LoadFITSFilesDialog(parent, title=title, headlist=headlist)
    loadfile = fd.go(dir_or_file=initialdir, default=default, pattern=pattern, key=key)
    return loadfile

def SaveFITSFile(parent=None, title=None, initialdir=os.curdir, default="", pattern="*.fit*", headlist="OBJECT", key=None):
    fd = SaveFITSFileDialog(parent, title=title, headlist=headlist)
    savefile = fd.go(dir_or_file=initialdir, default=default, pattern=pattern, key=key)
    return savefile



# Added some methods from the tkFileDialog module. This makes this module more
# complete and gives a similar interface to the routines from both modules.
#
# The 'result = None' statements avoid that an empty string is returned. An
# empty string could be interpreted as a valid answer.

def GetDirectory(parent=None, initialdir=os.curdir, **opts):
  result = tkFileDialog.askdirectory(parent=parent, initialdir=initialdir)
  if not result: result = None
  return result

def GetOther(parent=None, prompt="", **opts):
  result = tkSimpleDialog.askstring('Give new value', prompt, parent=parent, **opts)
  if not result: result = None
  return result

def GetInteger(parent=None, prompt="", **opts):
  result = tkSimpleDialog.askinteger('Give new value (integer)', prompt, parent=parent, **opts)
  if not result: result = None
  return result


################################################################################

# Test method similar to the one in FileDialog

def test():
    """Simple test program."""
    root = Tk()
    root.withdraw()
    fd = LoadFiles(root)
    loadfiles = fd.go()
    print loadfiles


if __name__ == '__main__':
    test()
