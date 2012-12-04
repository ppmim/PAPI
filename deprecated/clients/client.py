
import Pyro.naming, Pyro.core
from Pyro.errors import NamingError
import threading
import time

N=4
#######################################

class MiThread(threading.Thread):  

	def __init__(self, num ):  
		threading.Thread.__init__(self)  
		self.num = num
		
	def run(self):  
		print "Soy el hilo", self.num

		# locate the NS
		locator = Pyro.naming.NameServerLocator()
		print 'Searching Name Server...',
		ns = locator.getNS()
		
		# resolve the Pyro object
		print 'finding object'
		try:
			name='sreduce_%d' %self.num
			URI=ns.resolve(name)
			print 'URI:',URI
		except NamingError,x:
			print 'Couldn\'t find object, nameserver says:',x
			raise SystemExit
	
		# create a proxy for the Pyro object, and return that
		obj = Pyro.core.getProxyForURI(URI)

		source_frame = '/disk-a/caha/panic/DATA/data_mat/orion002%d.fits' %self.num
		dark_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_dark.fits'
		flat_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_normflat.fits'
		out_frame    = '/disk-a/caha/panic/DATA/data_mat/out/prueba1.fits'
		appMask      = False
		
		obj.run (source_frame, dark_frame, flat_frame, out_frame, True)
		

###### main client program ###########

def main():

	start = time.time()
	
	threads = []
	for i in range(1,N+1):
		threads.append( MiThread(i) )
		threads[i-1].start()

	for t in threads:
		t.join()

	print "Finalizaron todas las HEBRAS en %f secs !!!" %(time.time()-start)

	#print obj.run(source_frame, dark_frame, flat_frame, out_frame)

############################################
	
if __name__=="__main__":
        main()




