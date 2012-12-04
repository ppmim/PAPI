
import Pyro.naming, Pyro.core
from Pyro.errors import NamingError

###### main client program ###########

def main():

	# locate the NS
	locator = Pyro.naming.NameServerLocator()
	print 'Searching Name Server...',
	ns = locator.getNS()
	
	# resolve the Pyro object
	print 'finding object'
	try:
		URI=ns.resolve('sreduce_3')
		print 'URI:',URI
	except NamingError,x:
		print 'Couldn\'t find object, nameserver says:',x
		raise SystemExit
	
	# create a proxy for the Pyro object, and return that
	obj = Pyro.core.getProxyForURI(URI)

	source_frame = '/disk-a/caha/panic/DATA/data_mat/QL1/orion0021_x4_3.fits'
	dark_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_dark.fits'
	flat_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_normflat.fits'
	out_frame    = '/disk-a/caha/panic/DATA/data_mat/out/prueba1.fits'
	appMask      = False
	print obj.run(source_frame, dark_frame, flat_frame, out_frame, True)

############################################
	
if __name__=="__main__":
        main()




