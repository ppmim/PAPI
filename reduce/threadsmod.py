#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2008-2015 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
#
# PAPI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Import requered modules

import threading
import time


import misc.display as display

#######################################
lock = threading.Lock()

class ExecTaskThread(threading.Thread):
    """ 
    Thread to execute a task and then signal with a event to a waiting 
    thread about the result of the task.
    """
    def __init__(self, task, task_info_list, *args):
    
        threading.Thread.__init__(self)
        self._task = task
        self._args = args
        #self._event=event # not used
        self._task_info = TaskInfo()
        self._task_info_list = task_info_list
    
    def run(self):
      
        lock.acquire()
        try:
            self._task_info._curr_status = "INITIATED"
            self._task_info._return      = self._task(*self._args)  # Execute the task
            self._task_info._exit_status = 0             # EXIT_SUCCESS, all was OK
            self._task_info._exc = None
        except Exception, e:
            self._task_info._curr_status = "FINISHED"
            self._task_info._return      = None
            self._task_info._exit_status = 1             # EXIT_FAILURE, some error happened
            self._task_info._exc = e  # Exception
            #raise e # If it is raised, I do not know where/how to catch it ??
        finally:
            self._task_info_list.append(self._task_info)
            lock.release()    
            #self._event.set() # signal for the waiting thread (consumer)

            
class WaitTaskThread(threading.Thread):
    """ 
    NOT USED until now
    NOT USED until now
    NOT USED until now
    
    Thread waiting for a task until a signal event is received
    
    NOT USED until now
    NOT USED until now
    NOT USED until now
    
    """
    def __init__(self, task_info, event):
    
        threading.Thread.__init__(self)
        self._event = event
        self._task_info = task_info
        self._run = True
                 
    def run(self):
      
        while self._run:
            print "WaitTaskThread waiting ..."
            self._event.wait()
            print "WaitTaskThread running ..."
            if self._task_info._exit_status == 0: # EXIT_SUCCESS, all was OK
                if self._task_info._return!=None:
                    display.showFrame(self._task_info._return)
            else:
                #QMessageBox.critical(None, "Error", "Error while running task ")
                print "[WaitTaskThread] Error while running task !"
                pass # nothing to do
            
            self._event.clear()   # restore the event condition

class TaskInfo():
    """ 
    Class where ExecTaskThread and WaitTaskThread share info and the result of execution
    """
    def __init__(self):
    
        self._name = None              # Name of the task executed
        self._curr_status = None       # Current status (INITIATED, FINISHED)
        self._exit_status = None       # 0 if all was OK or 1 if some error happened
        self._return = None            # Returned value from the task
        self._exc = None               # Exception launched when error happen
        
    def clear(self):
                  
        self._name = None
        self._exit_status = None
        self._return = None
        self._exc = None
