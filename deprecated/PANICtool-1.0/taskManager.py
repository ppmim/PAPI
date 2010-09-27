################################################################################
#
# FIEStool
#
# taskManager.py
#
# Last update 24/05/2005
#
################################################################################

"""
   Routine that manages and wraps around the different reduction tasks. Any
   errors that occur when executing tasks are caught and relayed to the message
   logging facility.
"""

# import external modules

import messageLog
import tasks
import traceback
import sys


################################################################################

def runTask(task, frame=None, *args, **kwds):

  """
     Execute 'task' with given '*args' and '**kwds'. Show task execution
     progress and all errors in messageLog. If an error occurs, will return
     a string stating so.
  """

  # Display status message that task is to be executed

  messageLog.put('Task manager invoked for %s with following keywords :' % task.name, 7)
  messageLog.put('%s %s' % (args, kwds), 7)

  try:
    # Execute the task
    if frame!=None:
      task().run(frame,*args, **kwds)
    else:
      task().run(*args, **kwds)

  except tasks.taskAbort, errString:
    return ''

  except tasks.taskError, errString:
    # Catch the 'standard' taskError, which indicates an 'expected' error
    messageLog.put('TaskManager: Encountered a taskError', 7)
    messageLog.put_error('An error occured in task "%s" :' % task.name)
    messageLog.put_error(errString)
    return 'Terminated by error'

  except Exception, inst:
    # Also catch all other errors and display the complete traceback log.
    # These errors are more severe and should not occur, although they will.
 
    messageLog.put('TaskManager: Encountered an unknown error', 7)
    messageLog.put_error('An unknown error occured in "%s" :' % task.name)
 
    # Get the error message in standard traceback format
    (exc_type, exc_value, exc_traceback) = sys.exc_info()
    errormessage = traceback.format_exception(exc_type, exc_value, exc_traceback)
    for message in errormessage:
      # Show all lines of traceback log
      for line in message.splitlines():
	messageLog.put_error(line)
    return 'Terminated by error'


  # If here, task has finished without errors
  messageLog.put('TaskManager: Task %s finished' % (task.name), 7)

  # No need to return error codes
  return None
