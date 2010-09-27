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


#######################################

class ReduceThread(threading.Thread):  

    """
    Main threaded class to execute the call to the remote object for reduction  
    """
 
    def __init__(self, resource_id, source_frame, dark_frame, flat_frame, out_frame ):  
        threading.Thread.__init__(self)
        
        self.resource_id  = resource_id 
        self.source_frame = source_frame
        self.dark_frame   = dark_frame
        self.flat_frame   = flat_frame
        self.out_frame    = out_frame
        self.error        = 0
        
    def run(self):  
        #print "Soy el hilo", self.num
        
        # locate the NS
        locator = Pyro.naming.NameServerLocator()
        #print 'Searching Name Server...',
        ns = locator.getNS()
        
        # resolve the Pyro object
        #print 'finding object'
        try:
            name='sreduce_%d' %self.resource_id
            URI=ns.resolve(name)
            #print 'URI:',URI
        except NamingError,x:
            print 'Couldn\'t find object, nameserver says:',x
            self.error=1
            raise
	
        # create a proxy for the Pyro object, and return that
        obj = Pyro.core.getProxyForURI(URI)
        #obj._setOneway('run')-->it would enable asynchronous call
        # Synchonous call to 'run'
        self.error = obj.run (self.source_frame, self.dark_frame, self.flat_frame, self.out_frame, True)

        
        
