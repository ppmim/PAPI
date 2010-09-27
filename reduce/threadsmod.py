################################################################################
#
# PANICtool
#
# ReduceThread.py
#
# Author : jmiguel@iaa.e
#
# Last update 12/05/2008
#
################################################################################

# Import requered modules

import Pyro.naming, Pyro.core
from Pyro.errors import NamingError
import threading
import time

from qt import QMessageBox

import misc.display as display

#######################################

class ReduceThread(threading.Thread):  

    """
    Main threaded Client class to execute the call to the remote object for reduction  
    """
 
    def __init__(self, resource_id, source_frame, dark_frame, flat_frame, out_frame, output_dir='/tmp/' ):
        threading.Thread.__init__(self)
        
        self.resource_id  = resource_id 
        self.source_frame = source_frame
        self.dark_frame   = dark_frame
        self.flat_frame   = flat_frame
        self.out_frame    = out_frame
        self.output_dir   = output_dir
        self.error        = 0
        
    def run(self):  
        print "[ReduceThread]:Soy el hilo cliente ", self.resource_id
        
        # locate the NS
        locator = Pyro.naming.NameServerLocator()
        #print 'Searching Name Server...',
        ns = locator.getNS()
        
        # resolve the Pyro object
        #print 'finding object'
        try:
            name='sreduce_%d' %self.resource_id
            URI=ns.resolve(name)
            print 'URI:',URI
        except NamingError,x:
            print 'Couldn\'t find object, nameserver says:',x
            self.error=1
            raise
	
        # create a proxy for the Pyro object, and return that
        obj = Pyro.core.getProxyForURI(URI)
        #obj._setOneway('run')-->it would enable asynchronous call
        # Synchonous call to 'run'
        self.error = obj.run (self.source_frame, self.dark_frame, self.flat_frame, self.out_frame, self.output_dir, False)

        print "Thread::Run return", self.out_frame
        
class ExecTaskThread(threading.Thread):
    """ 
    Thread to execute a task and then signal with a event to a waiting thread about the result of the task
    """
    def __init__(self, task, task_info_list):
    
        threading.Thread.__init__(self)
        self._task=task
        #self._event=event # not used
        self._task_info = TaskInfo()
        self._task_info_list = task_info_list
    
    def run(self):
        
        try:
            try:
                self._task_info._curr_status = "INITIATED"
                self._task_info._return      = self._task()  # Execute the task
                self._task_info._exit_status = 0             # EXIT_SUCCESS, all was OK
            except Exception, e:
                self._task_info._curr_status = "FINISHED"
                self._task_info._return      = None
                self._task_info._exit_status = 1             # EXIT_FAILURE, some error happened
                self._task_info._exc = e
                raise e
        finally:
            self._task_info_list.append(self._task_info)
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
                QMessageBox.critical(None, "Error", "Error while running task ")
                pass # nothing to do
            
            self._event.clear()   # restore the event condition

class TaskInfo ():
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
                          
