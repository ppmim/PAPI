################################################################################
#
# FIEStool
#
# autoLoader.py
#
# Last update 29/09/2005
#
################################################################################

"""
   Performs automatic loading of configuration files based on the contents
   of FITS headers, and provide a GUI to define and edit selection rules..
"""

# Import external modules

import config
import copy
import os
import tkMessageBox
import pyfits
import myDialog
import popUp

from Tkinter import *

################################################################################

# Define local constants


# List of keys that define a rule, and the header values that will be displayed
# when editing these rules

keylist = ('rulenumber', 'fitshead', 'expression', 'iftrue', 'iffalse', 'filename')

header = {}
header['rulenumber']	=	'Rule no.'
header['fitshead']	=	'FITS header'
header['expression']	=	'Expression'
header['iftrue']	=	'If true'
header['iffalse']	=	'If false'
header['filename']	=	'Configuration file'


# Popup help texts for the autoLoader fields and the control buttons

helptext = {}
helptext['rulenumber']	= """Evaluation of header values and expressions will start from the
			     rule with the lowest number."""
helptext['fitshead']	= """When evaluating an expressiong, the value of this FITS header will
			     be determined. This field must not be empty."""
helptext['expression']	= """This expression will be appended to the value of the corresponding FITS
			     header. The expression must follow Python syntax, for example "== 2", or
			     "> 7". The action that is performed after evaluating the expression
			     depends on the value of the fields 'If true' and 'If false'."""
helptext['iftrue']	= """Determines the action taken if the rule's expression evaluates as TRUE.
			     If this field contains a rule number, further logical processing will
			     branch to this rule. If this field contains 'load', this rule's selected
			     configuration file will be loaded. If there is no value, no further
			     rules will be considered."""
helptext['iffalse']	= """Determines the action taken if the rule's expression evaluates as FALSE.
			     If this field contains a rule number, further logical processing will
			     branch to this rule. If this field contains 'load', this rule's selected
			     configuration file will be loaded. If there is no value, no further
			     rules will be considered."""
helptext['filename']	= """The name of a pipeline configuration file that may be loaded when
			     the field 'If true' or 'If false' contains the command 'load'."""
helptext['addrule']	= """Add another rule to the current set of rules."""
helptext['testrules']	= """Test the current set of rules using an existing FITS file. The
                             interface will report the resulting action, as well as any
			     encountered errors. Note that rules that are not evaluated are
			     also not tested."""
helptext['close']	= """Close this window. The current ruleset will be used if the 'AutoLoad'
			     task is called."""
helptext['refresh']	= """Redraw the current window."""


################################################################################

class RuleError(Exception):

  """
     Custom exception class in case an error occurs when processing rules.
  """

  def __init__(self, *args):
    self.args = args

################################################################################


class UserDict:

  """
     Custom class mimicing the behavior of a dictionary object.
     This class in necessary because inheriting methods from the dictionary
     object does not give all the needed functionality, and no further object
     can be attached.
  """

  def __init__(self): self.data = {}
  def __repr__(self): return repr(self.data)
  def __cmp__(self, dict):
          if type(dict) == type(self.data):
                  return cmp(self.data, dict)
          else:
                  return cmp(self.data, dict.data)
  def __len__(self): return len(self.data)
  def __getitem__(self, key): return self.data[key]
  def __setitem__(self, key, item): self.data[key] = item
  def __delitem__(self, key): del self.data[key]
  def keys(self): return self.data.keys()
  def items(self): return self.data.items()
  def values(self): return self.data.values()
  def has_key(self, key): return self.data.has_key(key)


################################################################################

class Rule:

  """
     Represents one single autoLoading rule. A valid rule consists of the name
     of a FITS keyword, an expression that can be parsed by Python (e.g. '==
     1'), the actions to be undertaken if the expression evaluates as true or
     as false, and (optionally) the name of the configuration file that should
     be loaded.
  """

  def __init__(self, init_val=("HEADER", "== 1", "", "", "")):

    (self.fitshead, self.expression, self.iftrue,
     self.iffalse, self.filename) = init_val 


  def __repr__(self):

    return repr( self.dump() )


  def dump(self):

    return (self.fitshead, self.expression, self.iftrue,
            self.iffalse, self.filename)

################################################################################

class ruleset:

  """
     Provides access and modification methods for a set of rules.
  """

  def __init__(self, initruleset=None):

    "Create a new set of rules"

    self.nrules = 0
    self.ruleset = dict()

    # Copy contents from old ruleset if provided
    if initruleset is not None: self.make(initruleset)


  def __repr__(self):
    return repr(self.dump())


  def __len__(self):
    return self.nrules


  def __getitem__(self, key):
    try:
      return self.ruleset[key]
    except KeyError: raise IndexError


  def values(self):
    return self.ruleset.values()

  def items(self):
    return self.ruleset.items()

  def keys(self):
    return self.ruleset.keys()


  def addrule(self, number, init_array=("", "", "", "", "")):

    "Add a new rule to the ruleset"

    self.ruleset[number] = Rule(init_array)
    self.nrules = self.nrules + 1


  def delrule(self, number):

    "Remove a given rule from the ruleset"

    del(self.ruleset[number])
    self.nrules = self.nrules - 1


  def clear(self):

    "Clear the contents of the ruleset"

    self.__init__()


# Deprecated method :
#
# def test(self):
#   "Test the validity of the ruleset"
#   for rule in self.ruleset.items():
#     try: result = rule.test()
#     except: raise


  def dump(self):

    "Return the current set of rules as a 2D array"

    outdata = []
    for (rulenumber, rule) in self.ruleset.items():
      outdata.append([rulenumber, rule.dump()])
    return outdata


  def make(self, indata):

    "Construct a new set of rules from an existing 2D array"

    self.clear()

    for (rulenumber, rule) in indata:
      self.addrule(rulenumber, rule)


  def run(self, fitsfile, number=0, count=0):

    """
       Process one single rule of the ruleset using the headers of a given FITS
       file, and either (1) execute this same function for another rule or (2)
       return the action to be undertaken. The action can be the name of a
       configuration file to load or 'None' if no further action is required.
    """

    # Keep track of how deep recursion goes
    count = count + 1

    # Leave if recursion seems to run too deep
    if (count > 100): raise RuleError, "Runaway error! Did more than 100 steps when parsing rules. Aborting."

    # Executing rule 0 means that the rule with the lowest number should be
    # executed. This will catch the situation when the first rule is not
    # number 0.
    if (number == 0): 
      try: number = min(self.ruleset.keys())
      except ValueError: raise RuleError, "There are no rules to consider."

    # Number may have been given as a string (normal during recursion)
    number = int(number)

    # If the rule does not exist, raise error
    try: thisrule = self.ruleset[number]
    except KeyError: raise RuleError, "Rule %i does not exist. Are you branching correctly?" % number

    # Open the FITS file for which the ruleset is executed and read the headers
    try: 
         infile = pyfits.open(fitsfile)
         headerlist = infile[0].header
    except : raise RuleError, "Cannot open or read the headers of FITS file %s" % fitsfile
    infile.close()

    # Test that the FITS header field is defined
    if not thisrule.fitshead: raise RuleError, "In rule %i: FITS header field empty" % number

    # Check that the value of the FITS header field makes sense
    try: headerval = headerlist[thisrule.fitshead.capitalize()]
    except KeyError: raise RuleError, "In rule %i: FITS header %s not found in %s" % (number, thisrule.fitshead.capitalize(), fitsfile)

    # Check that the FITS header has a value (to test against)
    if headerval == None: raise RuleError, "In rule %i: FITS header '%s' has no value in %s" % (number, thisrule.fitshead.capitalize(), fitsfile)

    # Construct the expression to be evaluated
    if type(headerval) == type(""):
      expression = '"%s" %s' % (headerval, thisrule.expression)
    else:
      expression = '%s %s' % (headerval, thisrule.expression)

    # Evaluate the expression
    try: result = eval(expression)
    except SyntaxError:
      raise RuleError, "In rule %i: Cannot parse expression (use Python logic) '%s'" % (number, expression)
    except Exception, errStr:
      raise RuleError, "In rule %i: (When testing '%s'): %s" % (number, expression, errStr)
    
    # If evaluation results in false, then take the defined action
    if result:
      if   thisrule.iftrue.isdigit():  return self.run(fitsfile, thisrule.iftrue, count)
      elif thisrule.iftrue == 'load':  return thisrule.filename
      elif thisrule.iftrue == 'abort':  return 'abort'
      elif thisrule.iftrue != '':
        raise RuleError, "In rule %i: Cannot interpret action '%s' in 'if true'. " % (number, thisrule.iftrue)

    # If evaluation results in true, then take the defined action
    elif not result:
      if   thisrule.iffalse.isdigit(): return self.run(fitsfile, thisrule.iffalse, count)
      elif thisrule.iffalse == 'load': return thisrule.filename
      elif thisrule.iffalse == 'abort': return 'abort'
      elif thisrule.iffalse != '':
        raise RuleError, "In rule %i: Cannot interpret action '%s' in 'if false'." % (number, thisrule.iffalse)

    # Nothing more to do... return without further actions.
    return None

################################################################################

class window:

  """
     Construct a window to modify and test a set of autoLoading rules
  """

  def __init__(self, parent, title='Configuration for autoloading configuration files'):

    """
       Define the initial ruleset (if given in configuration) and create the
       rule editing frame
    """

    self.parent = parent

    currentrules = config.options['autoloadrules'].value

    # self.ruleset   will contain the actual ruleset object.
    # self.tkentries will contain the Tkinter objects that contain
    #                the values on screen

    self.ruleset = ruleset(currentrules)
    self.tkentries = UserDict()

    # Create a new top-level frame, but keep it hidden
    self.root = Toplevel(self.parent)
    self.root.withdraw()
    self.root.title(title)

    # Make the frame
    self.makeFrame()

    # And fill it with the current values
    self.writeentries()


  def makeFrame(self):

    "Construct the actual frame"
    
    self.topframe = Frame(self.root)

    i = 0 # Current row
    j = 0 # Current column

    # Display the table headers
    for key in keylist:
      thislabel = Label(self.topframe, width=20, text=header[key])
      thislabel.grid(row=i, column=j, sticky=E+W)
      popUp.popUp(self.topframe, thislabel, title=header[key], text=helptext[key])
      j = j + 1
 
    i = i + 1

    # Check that some rules are defined
    if len(self.ruleset) == 0:

      # No rules defined (yet), so display a message about this
      Label(self.topframe, bg='green',
        text="There are no active rules. Press 'Add rule' to create a set of rules").grid(row=i,
	columnspan=7, sticky=E+W)
      i = i + 1

    else:


      # Display each rule    
      for (rulenumber, rule) in self.ruleset.items():

	self.tkentries[rulenumber] = Rule()

	j = 0

        # Make a numbered label for this rule
	ruleLabel = Label(self.topframe, width=5, text=str(rulenumber))
	ruleLabel.grid(row=i, column=j)
	# Attach a popup window with helptext
	popUp.popUp(self.topframe, ruleLabel,
                    title=header['rulenumber'], text=helptext['rulenumber'])

	j = j + 1

        # Create entry fields for each element of the rule
	for key in keylist[1:5]:
	  thisentry = Entry(self.topframe, width=15)
	  thisentry.grid(row=i, column=j)

	  # Attach each entry to the self.tkentries object (for future ref)
          setattr(self.tkentries[rulenumber], key, thisentry)
          popUp.popUp(self.topframe, thisentry, title=header[key], text=helptext[key])
          j = j + 1


        # Treat the config filename slightly differently (custom button)
	configfile = rule.filename
	if not configfile: configfile = "<None>"
	else : configfile = os.path.basename(configfile)

	configfilebutton = Button(self.topframe, width=25, text=configfile)
	configfilebutton.config(command=lambda n=rulenumber: self.pickFile(n))
	configfilebutton.grid(row=i, column=j)
	popUp.popUp(self.topframe, configfilebutton,
                    title=header['filename'], text=helptext['filename'])

        # Also here, attach the button to self.tkentries
	self.tkentries[rulenumber].filename = configfilebutton
	j = j + 1


        # Define an additional button to remove this rule
	delrulebutton = Button(self.topframe, width=15, text='Delete',
                               command=lambda n=rulenumber: self.delRule(n))
	delrulebutton.grid(row=i, column=j)
	popUp.popUp(self.topframe, delrulebutton, title="Delete rule",
                    text="Deletes this rule from the set of rules.")

	i = i + 1



    # At the bottom of the screen, add a set of useful buttons...

    addbutton  = Button(self.topframe, text='Add new rule')
    addbutton.grid(row=i, column=0, columnspan=2, sticky=E+W)
    addbutton.config(command = self.addRule)
    popUp.popUp(self.topframe, addbutton, title='Add new rule', text=helptext['addrule'])

    testbutton = Button(self.topframe, text='Test ruleset')
    testbutton.grid(row=i, column=2, columnspan=2, sticky=E+W)
    testbutton.config(command = self.testRuleSet)
    popUp.popUp(self.topframe, testbutton, title='Test ruleset', text=helptext['testrules'])

    closebutton = Button(self.topframe, text='Close window')
    closebutton.grid(row=i, column=4, columnspan=2, sticky=E+W)
    closebutton.config(command = self.hide)
    popUp.popUp(self.topframe, closebutton, title='Close window', text=helptext['close'])

    refrbutton = Button(self.topframe, text='Refresh')
    refrbutton.grid(row=i, column=6, columnspan=1, sticky=E+W)
    refrbutton.config(command = self.reset)
    popUp.popUp(self.topframe, refrbutton, title='Refresh', text=helptext['refresh'])


    # Display the actual frame (but may still be hidden from view)
    self.topframe.grid()


  def show(self):

    "Make window visible and give focus"

    self.root.deiconify()
    self.root.focus_set()
    self.root.grab_set()
    self.root.wait_visibility()
    

  def hide(self):

    "Hide window from view. Release grab and return to parent."

    self.root.grab_release()
    self.root.withdraw()

    # Also now, read the values of the entries on screen and store these
    self.readentries()


  def readentries(self):

    """
       Read the values that are currently in the widget, and store these values
       in the ruleset. Copy the existing ruleset to the system configuration
       object.
    """

    for (rulenumber, entries) in self.tkentries.items():

      for key in keylist[1:5]:

        # Get the value from the screen
        value = getattr(entries, key).get()
	# Store it in the ruleset
	setattr(self.ruleset[rulenumber], key, value)

    # Save the current ruleset in the configuration object
    config.options['autoloadrules'].value = self.ruleset.dump()


  def writeentries(self):

    """
       Take the values of the ruleset object and write these to the widget
       elements on the screen
    """

    for (rulenumber, ruleitems) in self.ruleset.items():

      for key in keylist[1:5]:

        # Get the value of the item in the ruleset
	value = getattr(ruleitems, key)

	# Write this value to the entry on screen
	thisentry = getattr(self.tkentries[rulenumber], key)
	thisentry.delete(0, END)
	thisentry.insert(END, value)


  def reset(self):
  
    "Refresh the window"

    self.readentries()
    self.topframe.destroy()
    self.makeFrame()
    self.writeentries()



  def testRuleSet(self):
  
    """
       Let the user select a FITS file to which the ruleset will be applied.
       Results are shown in a popup window.
    """

    # Select a FITS file
    fitsfile = myDialog.LoadFITSFile(parent=self.root, title='Select a FITS file for testing',
                    initialdir=config.options['inpath'].value,
                    headlist=config.options['fitsheaders'].value)

    # Did it succeed?
    if not fitsfile: return

    # Get the current values from the screen
    self.readentries()

    # Try to execute the ruleset. If an error occurs, warn the user about this
    try:
      loadfile = self.ruleset.run(fitsfile)
    except RuleError, errString:
      tkMessageBox.showwarning(message=errString)
      return

    # Ruleset execution terminated normally. Display the action that would
    # be taken
    if loadfile == 'abort':
      tkMessageBox.showinfo(message="I would have aborted processing this file")
    elif loadfile:
      tkMessageBox.showinfo(message="I would have loaded %s" % loadfile)
    else:
      tkMessageBox.showinfo(message="No more rules to consider. No new file would have been loaded")



  def delRule(self, rulenumber):
  
    "Remove one rule from the ruleset and update the screen"
  
    self.readentries()
    self.ruleset.delrule(rulenumber)
    del(self.tkentries[rulenumber])
    self.reset()



  def addRule(self):

    "Add a new rule to the ruleset and update the screen"

    self.readentries()

    # Give the new rule an unique number
    try: newrulenumber = max(self.ruleset.keys()) + 1
    except ValueError: newrulenumber = 1

    self.ruleset.addrule(newrulenumber)
    self.reset()


  def pickFile(self, rulenumber):

    "Select a configuration file"

    newfile = myDialog.LoadFile(parent=self.root, title='Select configuration file',
                    initialdir=config.options['config_dir'].value,
                    filetypes=[('Config files', '*.cfg'), ('All files', '*')])

    if newfile:
      self.tkentries[rulenumber].filename.config(text = os.path.basename(newfile))
      self.ruleset[rulenumber].filename = newfile



  def saveConfig(self, outfile=None):

    "Save the autoloader configuration data to disk"

    if outfile is None:

      # Ask for filename if not yet given
      outfile = myDialog.SaveFile(parent=self.root, title='Select autoloader configuration file',
                     initialdir=config.options['config_dir'].value,
                     filetypes=[('Autoloader config files', '*.acfg'), ('All files', '*')])


    if outfile:

      # Call config.save method for the listed options (only!)
      config.save_raw(outfile, optionlist=['autoloadrules'])


  def loadConfig(self, infile=None):

    "Load the autoloader configuration options from disk"

    if infile is None:

      # Ask for filename if not yet given
      infile = myDialog.LoadFile(parent=self.root, title='Select autoloader configuration file',
                     initialdir=config.options['config_dir'].value,
                     filetypes=[('Autoloader config files', '*.acfg'), ('All files', '*')])


    if infile:
      # Read listed options (only!) from file into config object
      infile = os.path.abspath(infile)
      config.load_raw(infile, optionlist=self.optionlist)

      # Update frame with new values
      self.reset()



################################################################################


