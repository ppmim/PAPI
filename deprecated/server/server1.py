import Pyro.naming
import Pyro.core
from Pyro.errors import PyroError,NamingError

import reduce

###### SimpleReduce Pyro object

class SimpleReduce(Pyro.core.ObjBase, reduce.SimpleReduce):
        pass


######## TLS INIT

def initTLS(tls):
	print "INIT TLS! TLS=",tls
	print "Setting counter for this TLS to 0."
	tls.counter=0
###### main server program

def main():

        Pyro.config.PYRO_MULTITHREADED = 1
	
        Pyro.core.initServer()
        daemon = Pyro.core.Daemon()
        # locate the NS
        locator = Pyro.naming.NameServerLocator()
        print 'searching for Name Server...'
        ns = locator.getNS()
        daemon.useNameServer(ns)

        # set TLS init func
        #daemon.setInitTLS(initTLS)
	
        # connect a new object implementation (first unregister previous one)
        try:
                # 'sreduce' is the name by which our object will be known to the outside world
                ns.unregister('sreduce_1')
        except NamingError:
                pass

        # connect new object implementation
        daemon.connect(SimpleReduce(),'sreduce_1')

        # enter the server loop.
        print 'Server object "sreduce_1" ready.'
        daemon.requestLoop()

if __name__=="__main__":
        main()

        
