#########################################
##  Data Reduction Server Manager       #
##  -----------------------------       #
#########################################

import os
import Pyro.naming
import Pyro.core
from Pyro.errors import PyroError,NamingError

import reduce.reducemod

import logging
import misc.paLog

###### simple_reduce Pyro object

class simple_reduce(Pyro.core.ObjBase, reduce.reducemod.simple_reduce):
        pass

###### advanced reduce Pyro object
class advanced_reduce(Pyro.core.ObjBase, reduce.reducemod.simple_reduce):
	pass


## Number of pipeline server to launch; Normally, one for each frame detector
N_SERVERS=4

log=logging.getLogger("panic.server")

###### main server program
def main():

	misc.paLog.initLog("/tmp/panic.server.log", logging.DEBUG)

	log.info("Starting Multiple server engine with %d childs" %N_SERVERS)
	
	for i in range(N_SERVERS):
		pid = os.fork()
		if pid ==0:
			#Child
			server_child(i+1)

	#Wait until child end	
	os.wait()[0]
		
####### Child process
def server_child( i ):


	
	Pyro.config.PYRO_MULTITHREADED = 1
	
	Pyro.core.initServer()
        daemon = Pyro.core.Daemon()
        # locate the NS
        locator = Pyro.naming.NameServerLocator()
        print 'searching for Name Server...'
        ns = locator.getNS()
        daemon.useNameServer(ns)

        # connect a new object implementation (first unregister previous one)
        try:
                # 'sreduce' is the name by which our object will be known to the outside world
                ns.unregister('sreduce_%d' %i)
        except NamingError:
                pass

        # connect new object implementation
        daemon.connect(simple_reduce(),'sreduce_%d' %i)

        # enter the server loop.
        log.info( 'Server object "sreduce_%d" ready .' %i)
        daemon.requestLoop()
	

if __name__=="__main__":
        main()

        
