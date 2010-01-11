
import Pyro.naming, Pyro.core
from Pyro.errors import NamingError

import logging
import misc.paLog

log = logging.getLogger('panic.client')

###### main client program ###########

def main():

	misc.paLog.initLog('/tmp/client1.log',logging.DEBUG)
	
	log.info("Test client 1 start ....")

	# locate the NS
	locator = Pyro.naming.NameServerLocator()
	print 'Searching Name Server...',
	ns = locator.getNS()
	
	# resolve the Pyro object
	print 'finding object'
	try:
		URI=ns.resolve('sreduce_1')
		print 'URI:',URI
	except NamingError,x:
		print 'Couldn\'t find object, nameserver says:',x
		raise SystemExit
	
	# create a proxy for the Pyro object, and return that
	obj = Pyro.core.getProxyForURI(URI)

	source_frame = '/disk-a/caha/panic/DATA/data_mat/QL1/orion0021_x4_1.fits'
	dark_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_dark.fits'
	flat_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_normflat.fits'
	out_frame    = '/disk-a/caha/panic/DATA/data_mat/out/prueba1.fits'
	appMask      = False
	#obj._setOneway('run')
	print obj.run(source_frame, dark_frame, flat_frame, out_frame, '/tmp', False)
    #frame_in, master_dark, master_flat, result_frame="/tmp/result.fits", out_tmp_dir='/tmp/', appPixMask=False):
    

############################################
	
if __name__=="__main__":
        main()




